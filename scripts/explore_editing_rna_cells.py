"""Connect accepted RNA allele observations to barcode-level RNA exploration.

No original BAM is opened. Original normalized alleles and NTC counts are reused;
per-family states must reproduce the accepted own-guide counts exactly. Company
barcodes remain merged records, not independently verified physical cells.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import socket
import time

import numpy as np
import pandas as pd

from analyze_guided_rna_alleles import (
    Read, covers_reference, event_key, has_observed_base, raw_events,
)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def atomic(path, value):
    path = Path(path)
    temp = path.with_suffix(path.suffix + '.writing')
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding='utf-8')
    temp.replace(path)


def reverse_complement(value):
    return value.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def add_same_indel_support(cohort):
    """Distinguish >=2 molecules carrying the same indel from different indels."""
    maximum = defaultdict(int)
    for variant in cohort['variants']:
        if variant['kind'] not in ['INS', 'DEL']:
            continue
        for cell in variant['cell_support']:
            key = (variant['guide'], cell['cell_barcode'])
            maximum[key] = max(maximum[key], cell['alt'])
    by_guide = defaultdict(int)
    for cell in cohort['cells']:
        value = maximum[(cell['guide'], cell['cell_barcode'])]
        if value > cell['indel_families']:
            raise ValueError('Same-allele molecule count exceeds all indel families')
        cell['max_same_indel_families'] = value
        by_guide[cell['guide']] += int(value >= 2)
    for guide in cohort['guides']:
        guide['same_indel_ge2_records'] = by_guide[guide['guide']]
    return cohort


def cut_position(spacer, reference, start, fetch):
    """Canonical SpCas9 NGG PAM, checked against the actual reference sequence."""
    if reference == spacer:
        pam = fetch(start + 20, start + 23)
        return dict(strand='+', PAM=pam, cut0=start + 17 if pam[1:] == 'GG' else None)
    if reference == reverse_complement(spacer):
        pam = reverse_complement(fetch(start - 3, start))
        return dict(strand='-', PAM=pam, cut0=start + 3 if pam[1:] == 'GG' else None)
    raise ValueError('Spacer does not match frozen reference at its unique placement')


def family_state(observations, variant):
    chrom, pos, ref, alt = variant
    states = set()
    for read, events in observations:
        if variant in events:
            states.add('alt')
        else:
            value = covers_reference(read, pos, pos + len(ref))
            if value is not None:
                states.add('ref' if value == ref else 'other')
    return next(iter(states)) if len(states) == 1 else ('ambiguous' if states else None)


def vcf_by_id(path):
    result = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith('#'):
            continue
        a = line.split('\t')
        if a[2] in result:
            raise ValueError('Duplicate original variant ID')
        result[a[2]] = (a[0], int(a[1]) - 1, a[3], a[4])
    return result


class DigestReader(io.RawIOBase):
    def __init__(self, source):
        self.source = source
        self.digest = hashlib.sha256()

    def readable(self):
        return True

    def read(self, size=-1):
        value = self.source.read(size)
        self.digest.update(value)
        return value


def load_cache(path, receipt, selected_samples):
    """Hash compressed bytes during this necessary analysis read, not a second pass."""
    reads = {}
    with Path(path).open('rb') as source:
        hashed = DigestReader(source)
        with gzip.GzipFile(fileobj=hashed, mode='rb') as compressed:
            with io.TextIOWrapper(compressed) as f:
                head = json.loads(next(f))
                if head['region'] != receipt['region']:
                    raise ValueError('Region/cache contract mismatch')
                for line in f:
                    item = json.loads(line)
                    if item['kind'] != 'original_SAM_record':
                        raise ValueError('Unexpected cache member')
                    if item['sample'] not in selected_samples:
                        continue
                    read = Read(item['sam'], item['sample'])
                    if read.family is None:
                        continue
                    key = read.signature
                    if key not in reads or read.quality > reads[key].quality:
                        reads[key] = read
        # Gzip has consumed the compressed stream including its footer.
        if hashed.digest.hexdigest() != receipt['sha256']:
            raise ValueError('Consumed cache content differs from accepted receipt')
    return list(reads.values())


def extract_region(root, receipt, catalogue, wanted):
    import pysam
    region = receipt['region']
    selected = [g for g in region['guides'] if g['guide'] in wanted]
    path = root / 'alleles' / region['region_id']
    old = json.loads((path / 'result.json').read_text())
    if old['source_cache_sha256'] != receipt['sha256']:
        raise ValueError('Existing allele result bound to different cache')
    raw = vcf_by_id(path / 'observed_events.vcf')
    normalized = vcf_by_id(path / 'normalized_events.vcf')
    if set(raw) != set(normalized):
        raise ValueError('Saved normalization identities differ')
    lookup = {raw[k]: normalized[k] for k in raw}
    samples = {g['sample_id'] for g in selected}
    reads = load_cache(root / receipt['file'], receipt, samples)
    # Only the locus-sized reference fragments are read; no genome-wide rehash.
    with pysam.FastaFile(extract_region.reference) as fasta:
        from functools import lru_cache
        fetch = lru_cache(maxsize=8192)(fasta.fetch)
        families = defaultdict(list)
        for read in reads:
            keys = {event_key(e) for e in raw_events(read, region, fetch)}
            if not keys <= set(lookup):
                raise ValueError('New read event absent from accepted normalization')
            families[read.family].append((read, {lookup[k] for k in keys}))
        guide_rows, cell_rows, variant_rows, display = [], [], [], []
        for guide in selected:
            name = guide['guide']
            meta = catalogue[name]
            start, end, chrom = guide['start0'], guide['end0'], region['chrom']
            ref = fetch(chrom, start, end)
            spacer = meta['sequence']
            if extract_region.feature_sequence_orientation == 'reverse_complement_of_spacer_5prime_capture':
                spacer = reverse_complement(spacer)
            cut = cut_position(spacer, ref, start,
                               lambda a, b: fetch(chrom, a, b))
            own = {k: v for k, v in families.items() if k[0] == guide['sample_id']}
            cells = defaultdict(lambda: dict(regional_families=set(), cut_ref_families=set(),
                                             indel_families=set(), SNV_families=set(),
                                             ambiguous_families=set(), variants=set()))
            rows = [r for r in old['variants'] if r['guide'] == name]
            for family, observations in own.items():
                c = cells[family[1]]
                if any(has_observed_base(r, start, end) for r, _ in observations):
                    c['regional_families'].add(family)
                if cut['cut0'] is not None:
                    a, b = cut['cut0'] - 3, cut['cut0'] + 3
                    observed = {covers_reference(r, a, b) for r, _ in observations}
                    observed.discard(None)
                    if observed == {fetch(chrom, a, b)}:
                        c['cut_ref_families'].add(family)
                    elif len(observed) > 1:
                        c['ambiguous_families'].add(family)
            for variant in rows:
                key = event_key(variant)
                counts = dict(alt=0, ref=0, other=0, ambiguous=0)
                per_cell = defaultdict(lambda: dict(alt=0, ref=0, other=0, ambiguous=0))
                for family, observations in own.items():
                    state = family_state(observations, key)
                    if state is None:
                        continue
                    counts[state] += 1
                    per_cell[family[1]][state] += 1
                    c = cells[family[1]]
                    if state == 'alt':
                        c['indel_families' if variant['kind'] in ['INS', 'DEL'] else 'SNV_families'].add(family)
                        c['variants'].add(key)
                    elif state == 'ambiguous':
                        c['ambiguous_families'].add(family)
                for state in counts:
                    if counts[state] != variant['guide_' + state]:
                        raise ValueError(f'Own-family count not reproduced: {name} {key} {state}')
                if counts['alt']:
                    row = dict(variant, target=meta['target_gene_name'], target_gene_id=meta['target_gene_id'],
                               cut0=cut['cut0'], cut_offset=None if cut['cut0'] is None else key[1] + 1 - cut['cut0'],
                               length_divisible_by_3=bool(variant['length_change'] % 3 == 0),
                               primary_197=name in extract_region.primary,
                               target_ORF_verified=False,
                               cell_support=[dict(cell_barcode=cb, **v) for cb, v in per_cell.items()])
                    variant_rows.append(row)
            for cb, c in cells.items():
                cell_rows.append(dict(guide=name, cell_barcode=cb, target=meta['target_gene_name'],
                                      target_gene_id=meta['target_gene_id'],
                                      **{k: len(v) for k, v in c.items()}))
            guide_rows.append(dict(guide=name, target=meta['target_gene_name'], target_gene_id=meta['target_gene_id'],
                                   region_id=region['region_id'], chrom=chrom, start0=start, end0=end, **cut,
                                   primary_197=name in extract_region.primary, scan_capped=receipt['scan_capped'],
                                   assigned_QC_records=int(meta['single_guide_after_stage1_qc']),
                                   observed_records=sum(bool(c['regional_families'] or c['indel_families'] or c['SNV_families']) for c in cells.values()),
                                   cut_ref_records=sum(bool(c['cut_ref_families']) for c in cells.values()),
                                   indel_records=sum(bool(c['indel_families']) for c in cells.values()),
                                   SNV_records=sum(bool(c['SNV_families']) for c in cells.values()),
                                   indel_ge2_records=sum(len(c['indel_families']) >= 2 for c in cells.values()),
                                   variant_rows=sum(r['guide'] == name for r in variant_rows)))
            # Original read sequences, grouped by corrected RNA molecule, for honest local alignment plots.
            keep = []
            for family, observations in own.items():
                if not any(has_observed_base(r, start, end) or events for r, events in observations):
                    continue
                best = max((r for r, _ in observations), key=lambda r: r.quality)
                keep.append(dict(cell_barcode=family[1], UMI=family[2], GX=family[3],
                                 SAM=best.sam, PCR_observations=len(observations)))
            display.append(dict(guide=name, chrom=chrom, start0=start - 15, end0=end + 15,
                                reference=fetch(chrom, start - 15, end + 15), cut=cut,
                                molecules=keep, note='One original representative per RNA family; conflicts retained in quantitative states'))
    return dict(guides=guide_rows, cells=cell_rows, variants=variant_rows, display=display)


def sequence(cfg, progress):
    extract_region.feature_sequence_orientation = cfg['feature_sequence_orientation']
    if extract_region.feature_sequence_orientation not in ['spacer_forward', 'reverse_complement_of_spacer_5prime_capture']:
        raise ValueError('The assay-specific feature sequence orientation must be explicit')
    selected = json.loads(Path(cfg['selected_guides']).read_text())
    extract_region.primary = {r['guide'] for r in selected}
    if len(extract_region.primary) != 197:
        raise ValueError('Expected exact primary cohort of 197 sgRNAs')
    targets = {r['target_gene_id'] for r in selected}
    if len(targets) != 104:
        raise ValueError('Expected exact 104 target identities')
    catalogue = pd.read_csv(cfg['catalogue'], keep_default_na=False).set_index('guide').to_dict('index')
    wanted = {g for g, r in catalogue.items() if r['target_gene_id'] in targets}
    output = Path(cfg['output']) / 'sequence'
    output.mkdir(parents=True, exist_ok=True)
    for dataset in cfg['datasets']:
        root = Path(dataset['sequence_root'])
        checkpoint = json.loads((root / 'cache_checkpoint.json').read_text())
        extract_region.reference = cfg['reference']
        cohorts = []
        chosen = [r for r in checkpoint['regions'] if any(g['guide'] in wanted for g in r['region']['guides'])]
        for i, receipt in enumerate(chosen):
            destination = output / dataset['id'] / (receipt['region']['region_id'] + '.json')
            destination.parent.mkdir(exist_ok=True)
            if destination.exists():
                value = json.loads(destination.read_text())
            else:
                value = extract_region(root, receipt, catalogue, wanted)
                atomic(destination, value)
            cohorts.append(value)
            progress('sequence', dataset=dataset['id'], completed=i + 1, total=len(chosen))
        combined = {key: [r for c in cohorts for r in c[key]] for key in ['guides', 'cells', 'variants', 'display']}
        atomic(output / (dataset['id'] + '.json'), combined)
        if dataset['id'] == 'company' and {r['guide'] for r in combined['guides'] if r['primary_197']} != extract_region.primary:
            raise ValueError('Primary cohort not completely represented')


def covariance_pca(reference, query, components=10):
    """All records enter ordinary centered, unscaled PCA; no response labels enter."""
    from scipy.linalg import eigh
    n = len(reference) + len(query)
    summed = reference.sum(axis=0) + query.sum(axis=0)
    covariance = reference.T @ reference + query.T @ query - np.outer(summed, summed) / n
    covariance /= n - 1
    p = covariance.shape[0]
    values, vectors = eigh(covariance, subset_by_index=[max(0, p - components), p - 1], check_finite=False)
    order = np.argsort(values)[::-1]
    values, vectors = values[order], vectors[:, order]
    for j in range(vectors.shape[1]):
        anchor = np.argmax(np.abs(vectors[:, j]))
        vectors[:, j] *= 1 if vectors[anchor, j] >= 0 else -1
    mean = summed / n
    return (reference - mean) @ vectors, (query - mean) @ vectors, values / np.trace(covariance)


def rna(cfg, progress):
    import anndata as ad
    from scipy.linalg import eigh
    from threadpoolctl import threadpool_limits
    targets = {r['target_gene_id'] for r in json.loads(Path(cfg['selected_guides']).read_text())}
    catalogue = pd.read_csv(cfg['catalogue'], keep_default_na=False).set_index('guide')
    out = Path(cfg['output']) / 'RNA'
    out.mkdir(exist_ok=True)
    for dataset in cfg['datasets']:
        a = ad.read_h5ad(dataset['RNA_checkpoint'])
        source = pd.read_csv(dataset['RNA_cells'], dtype={'cell_barcode': str, 'guide_id': str}).set_index('cell_barcode', verify_integrity=True)
        meta = source.loc[a.obs_names.astype(str)].copy()
        for key in ['guide_id', 'role']:
            if not np.array_equal(meta[key].astype(str).to_numpy(), a.obs[key].astype(str).to_numpy()):
                raise ValueError('RNA/guide/record identities do not match')
        if not np.array_equal(meta.total_counts.to_numpy(dtype=float), a.obs.total_counts.to_numpy(dtype=float)):
            raise ValueError('Authoritative full-RNA totals differ')
        if a.shape[1] != 2000:
            raise ValueError('Frozen 2000-RNA panel changed')
        meta['target_gene_id'] = meta.guide_id.map(catalogue.target_gene_id)
        if meta.target_gene_id.isna().any():
            raise ValueError('Unmapped sgRNA identity')
        ntc = meta.role.eq('non_targeting_control').to_numpy()
        x = a.X.toarray().astype(np.float64)
        background = x[ntc]
        base_sum = background.sum(axis=0)
        with threadpool_limits(limits=2):
            cross = background.T @ background
        summaries = []
        dest = out / dataset['id']
        dest.mkdir(exist_ok=True)
        for i, target in enumerate(sorted(targets)):
            path = dest / (target + '.npz')
            mask = meta.target_gene_id.eq(target).to_numpy() & ~ntc
            query = x[mask]
            row = dict(target_gene_id=target, target=catalogue.loc[catalogue.target_gene_id.eq(target), 'target_gene_name'].iloc[0],
                       dataset=dataset['id'], records=int(mask.sum()), guides=int(meta.loc[mask, 'guide_id'].nunique()),
                       NTC_records=int(ntc.sum()), original_RNA_QC_records=len(source), available_RNA_records=len(meta),
                       target_RNA_removed_from_PCA=bool(target in a.var_names))
            if not path.exists() and len(query):
                use = np.asarray(a.var_names != target)
                n = len(background) + len(query)
                summed = base_sum + query.sum(axis=0)
                with threadpool_limits(limits=2):
                    covariance = cross + query.T @ query - np.outer(summed, summed) / n
                    covariance = covariance[np.ix_(use, use)] / (n - 1)
                    p = covariance.shape[0]
                    values, vectors = eigh(covariance, subset_by_index=[p - 10, p - 1], check_finite=False)
                    order = np.argsort(values)[::-1]
                    values, vectors = values[order], vectors[:, order]
                    for j in range(vectors.shape[1]):
                        anchor = np.argmax(np.abs(vectors[:, j])); vectors[:, j] *= 1 if vectors[anchor, j] >= 0 else -1
                    mean = summed[use] / n
                    bp = (background[:, use] - mean) @ vectors
                    qp = (query[:, use] - mean) @ vectors
                with path.with_suffix('.npz.writing').open('wb') as f:
                    np.savez_compressed(f, NTC_PCA=bp.astype(np.float32), target_PCA=qp.astype(np.float32),
                                        variance_ratio=values / np.trace(covariance),
                                        cell_barcode=meta.index.to_numpy(dtype=str)[mask],
                                        guide=meta.loc[mask, 'guide_id'].to_numpy(dtype=str),
                                        total_RNA=meta.loc[mask, 'total_counts'].to_numpy(dtype=np.float64))
                path.with_suffix('.npz.writing').replace(path)
            summaries.append(row)
            progress('RNA_PCA', dataset=dataset['id'], completed=i + 1, total=len(targets), target=row['target'], records=row['records'])
        pd.DataFrame(summaries).to_csv(out / (dataset['id'] + '_targets.csv'), index=False)


def main():
    p = argparse.ArgumentParser(); p.add_argument('--config', type=Path, required=True)
    p.add_argument('--phase', choices=['sequence', 'RNA', 'all'], default='all'); args = p.parse_args()
    cfg = json.loads(args.config.read_text())
    if os.getuid() != 2601 or socket.gethostname() != 'rhesusbase':
        raise ValueError('Wrong analysis host/owner')
    if hasattr(os, 'sched_setaffinity'):
        os.sched_setaffinity(0, sorted(os.sched_getaffinity(0))[:2])
    root = Path(cfg['output'])
    if not str(root).startswith('/home/user/data/wusy/260918_perturbseq/lab/analyses/method_comparison/20261103_'):
        raise ValueError('Wrong new data namespace')
    root.mkdir(parents=True, exist_ok=True)
    contract = root / 'analysis_contract.json'
    if contract.exists() and json.loads(contract.read_text()) != cfg:
        raise ValueError('Previous output belongs to a different fixed analysis contract')
    atomic(contract, cfg)
    for path, expected in cfg.get('input_sha256', {}).items():
        if sha(path) != expected:
            raise ValueError('Input content differs: ' + path)
    lock = root / 'active.lock'
    fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    stat = Path('/proc/self/stat').read_text(); startticks = int(stat[stat.rfind(')') + 2:].split()[19])
    os.write(fd, json.dumps(dict(pid=os.getpid(), uid=os.getuid(), startticks=startticks,
                                command=os.sys.argv, commit=cfg['commit'])).encode()); os.close(fd)
    def progress(stage, **values):
        data = dict(stage=stage, at=time.strftime('%Y-%m-%dT%H:%M:%S%z'), **values)
        atomic(root / 'status.json', data); print(json.dumps(data), flush=True)
    try:
        if args.phase in ['sequence', 'all']:
            sequence(cfg, progress)
        if args.phase in ['RNA', 'all']:
            rna(cfg, progress)
        progress('complete_pending_figures_and_interpretation')
    finally:
        lock.unlink()


if __name__ == '__main__':
    main()
