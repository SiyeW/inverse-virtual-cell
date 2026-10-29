"""Locate all design spacers with STAR and inventory indexed RNA coverage.

This is an observability assessment, not a variant caller or KO-efficiency
estimator. Mapping is exact genomic placement without an assumed PAM/cut offset.
Spliced RNA, barcode ambiguity and missing coverage are never called unedited.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1024*1024), b''):
            h.update(b)
    return h.hexdigest()


def atomic(path, value):
    path = Path(path)
    tmp = path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value, indent=2), encoding='utf8')
    tmp.replace(path)


def identity(path):
    s = Path(path).stat()
    return dict(bytes=s.st_size, mtime_ns=s.st_mtime_ns)


def mapping_command(star, genome, fasta, prefix):
    return [star, '--genomeDir', str(genome), '--readFilesIn', str(fasta),
            '--runThreadN', '2', '--genomeLoad', 'NoSharedMemory',
            '--outFileNamePrefix', str(prefix), '--outSAMtype', 'SAM',
            '--outSAMunmapped', 'Within', '--outSAMattributes', 'NH', 'HI', 'AS', 'nM',
            '--alignEndsType', 'EndToEnd', '--alignIntronMax', '1',
            '--alignSJoverhangMin', '999', '--alignSJDBoverhangMin', '999',
            '--outFilterMismatchNmax', '0', '--outFilterMatchNmin', '20',
            # STAR 2.7.2a reports AS=19 for an exact 20M alignment. Exactness
            # is defined by 20M/nM=0 plus direct FASTA comparison, not AS=20.
            '--outFilterScoreMin', '0', '--outFilterMultimapNmax', '1000',
            '--outSAMmultNmax', '-1', '--limitOutSAMoneReadBytes', '1000000',
            '--outTmpKeep', 'All']


def parse_mappings(path, rows, reference_fetch=None):
    hits = {f'g{i:04d}': set() for i in range(len(rows))}
    reported_nh = {name: 0 for name in hits}
    with Path(path).open() as f:
        for line in f:
            if line.startswith('@'):
                continue
            a = line.rstrip('\n').split('\t')
            if a[0] not in hits:
                raise ValueError('Unexpected mapped query identity')
            if int(a[1]) & 4:
                continue
            tags = {x.split(':', 2)[0]: x.split(':', 2)[2] for x in a[11:]}
            reported_nh[a[0]] = max(reported_nh[a[0]], int(tags.get('NH', 0)))
            if a[5] != '20M' or int(tags.get('nM', '-1')) != 0:
                continue
            start = int(a[3])-1
            if reference_fetch is not None:
                observed=reference_fetch(a[2],start,start+20).upper()
                query=rows[int(a[0][1:])]['sequence'].upper()
                if int(a[1]) & 16:
                    query=query.translate(str.maketrans('ACGT','TGCA'))[::-1]
                if observed!=query:
                    continue
            hits[a[0]].add((a[2], start, start+20, '-' if int(a[1]) & 16 else '+'))
    output = []
    for i, row in enumerate(rows):
        key = f'g{i:04d}'
        found = sorted(hits[key])
        result = dict(guide=row['guide'], target_gene_name=row['target_gene_name'],
                      source_control_type=row['source_control_type'], sequence=row['sequence'],
                      exact_genomic_placements=len(found), reported_alignment_count=reported_nh[key],
                      mapping_status='not_exactly_located', chromosome=None, start0=None, end0=None,
                      sequence_alignment_strand=None)
        if len(found) == 1 and reported_nh[key] == 1:
            result.update(mapping_status='unique_exact_genomic_placement',
                          chromosome=found[0][0], start0=found[0][1], end0=found[0][2],
                          sequence_alignment_strand=found[0][3])
        elif found or reported_nh[key] > 1:
            result['mapping_status'] = 'ambiguous_exact_placement'
        output.append(result)
    return output


def read_observation(read, start, end):
    """Count sequence opportunity; an N skip is never an observed deletion."""
    cursor = read.reference_start
    blocks, skipped = [], False
    for op, length in read.cigartuples or []:
        if op in [0, 7, 8]:
            blocks.append((cursor, cursor+length))
            cursor += length
        elif op in [2, 3]:
            if op == 3 and cursor < end and cursor+length > start:
                skipped = True
            cursor += length
    def overlap(left, right):
        return sum(max(0, min(right,b)-max(left,a)) for a,b in blocks)
    any_base = overlap(start,end) > 0
    spans = (not skipped and overlap(start-3,start)==3 and overlap(end,end+3)==3)
    return any_base, spans


def locus_coverage(bam, row, assigned, ntc_guides, max_reads):
    start, end = row['start0'], row['end0']
    result = dict(row, records_fetched=0, primary_MAPQ30_reads=0,
                  primary_reads_overlapping_spacer=0, primary_reads_spanning_spacer=0,
                  spanning_reads_without_CB=0, assigned_spanning_reads_without_UB=0,
                  scan_capped=False)
    families = {x:set() for x in ['own_guide','NTC','other_assigned']}
    barcodes = {x:set() for x in families}
    for read in bam.fetch(row['chromosome'],max(0,start-3),end+3):
        result['records_fetched'] += 1
        if result['records_fetched'] > max_reads:
            result['scan_capped'] = True
            break
        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.mapping_quality < 30:
            continue
        result['primary_MAPQ30_reads'] += 1
        any_base, spans = read_observation(read,start,end)
        result['primary_reads_overlapping_spacer'] += int(any_base)
        result['primary_reads_spanning_spacer'] += int(spans)
        if not spans:
            continue
        if not read.has_tag('CB'):
            result['spanning_reads_without_CB']+=1
            continue
        cb = read.get_tag('CB')
        guide = assigned.get(cb)
        if guide is None:
            continue
        category = 'own_guide' if guide == row['guide'] else ('NTC' if guide in ntc_guides else 'other_assigned')
        barcodes[category].add(cb)
        if read.has_tag('UB'):
            families[category].add((cb,read.get_tag('UB')))
        else:
            result['assigned_spanning_reads_without_UB']+=1
    for category in families:
        result[category+'_assigned_barcodes_spanning'] = len(barcodes[category])
        result[category+'_distinct_CB_UB_families_spanning'] = len(families[category])
    result['coverage_status'] = 'RNA_spanning_evidence_observed' if families['own_guide'] else 'no_spanning_RNA_family_in_assigned_records'
    if barcodes['own_guide'] and not families['own_guide']:
        result['coverage_status']='own_assigned_reads_span_but_UMI_identity_missing'
    if result['scan_capped']:
        result['coverage_status'] = 'read_cap_reached_lower_bound_only'
    return result


def run(config, progress):
    import pysam
    output = Path(config['output']);output.mkdir(parents=True,exist_ok=True)
    with Path(config['guides']).open() as f:
        rows = list(csv.DictReader(f))
    if digest(config['guides']) != config['guides_sha256']:
        raise ValueError('Guide source content differs')
    if len(rows)!=1196 or len({r['guide'] for r in rows})!=1196:
        raise ValueError('Unexpected design identity/coverage')
    if any(len(r['sequence'])!=20 or set(r['sequence'].upper())-set('ACGT') for r in rows):
        raise ValueError('Invalid 20nt spacer')
    for name in ['bam','bam_index','genome_fasta','genome_fasta_index']:
        if identity(config[name]) != config[name+'_identity']:
            raise ValueError(name+' identity differs')
    state = dict(guides_sha256=config['guides_sha256'],
                 bam_identity=config['bam_identity'],bam_index_identity=config['bam_index_identity'],
                 fasta_identity=config['genome_fasta_identity'],
                 max_reads_per_site=config.get('max_reads_per_site',100000),
                 method='exact_20nt_STAR_placement_then_indexed_RNA_observability_v1')
    mapping = output/'mapping';mapping.mkdir(exist_ok=True)
    placement_receipt = mapping/'placement_receipt.json'
    sam = mapping/'Aligned.out.sam'
    if placement_receipt.exists():
        p=json.loads(placement_receipt.read_text())
        if p['identity']!=state or p['sam_sha256']!=digest(sam):
            raise ValueError('Saved placement identity differs')
    else:
        if (mapping/'Log.out').exists():
            raise ValueError('Unaccepted previous alignment exists; preserve and diagnose')
        fasta=mapping/'guides.fa'
        fasta.write_text(''.join(f'>g{i:04d}\n{r["sequence"].upper()}\n' for i,r in enumerate(rows)))
        args=mapping_command(config['star'],config['star_index'],fasta,str(mapping)+'/')
        progress('locate_all_1196_design_spacers',reference_index_bytes=config['reference_index_bytes'])
        with (output/'STAR.stdout').open('w') as out,(output/'STAR.stderr').open('w') as err:
            subprocess.run(args,stdout=out,stderr=err,check=True,timeout=600)
        atomic(placement_receipt,dict(identity=state,sam_sha256=digest(sam),
                                     fasta_sha256=digest(fasta),command=args))
    with pysam.FastaFile(config['genome_fasta']) as reference:
        placements=parse_mappings(sam,rows,reference.fetch)
        reference_lengths=dict(zip(reference.references,reference.lengths))
    import gzip
    opener=gzip.open if str(config['assigned_cells']).endswith('.gz') else open
    if digest(config['assigned_cells'])!=config['assigned_cells_sha256']:
        raise ValueError('Cell assignment source differs')
    with opener(config['assigned_cells'],'rt') as f:
        cells=list(csv.DictReader(f))
    assigned={r['cell_barcode']:r['guide_id'] for r in cells}
    if len(assigned)!=len(cells):raise ValueError('Duplicated assigned barcode records')
    from collections import Counter
    assigned_counts=Counter(assigned.values())
    ntc={r['guide'] for r in rows if r['source_control_type']=='non_targeting_control'}
    state['assigned_cells_sha256']=config['assigned_cells_sha256']
    ledger=output/'coverage_checkpoint.json'
    saved=json.loads(ledger.read_text()) if ledger.exists() else dict(identity=state,rows=[])
    if saved['identity']!=state:raise ValueError('Coverage checkpoint source identity differs')
    done={r['guide']:r for r in saved['rows']}
    with pysam.AlignmentFile(config['bam'],'rb',index_filename=config['bam_index'],threads=1) as bam:
        if dict(zip(bam.references,bam.lengths))!=reference_lengths:
            raise ValueError('RNA BAM and placement reference contig identities differ')
        for i,row in enumerate(placements):
            if row['guide'] in done:continue
            row=dict(row,n_assigned_QC_records=assigned_counts[row['guide']])
            if not row['n_assigned_QC_records']:
                result=dict(row,coverage_status='not_assessed_no_assigned_QC_records')
            elif row['mapping_status']=='unique_exact_genomic_placement':
                result=locus_coverage(bam,row,assigned,ntc,config.get('max_reads_per_site',100000))
            else:
                result=dict(row,coverage_status='not_assessed_placement_unresolved')
            done[row['guide']]=result
            if (i+1)%20==0 or i==len(placements)-1:
                atomic(ledger,dict(identity=state,rows=list(done.values())))
                progress('indexed_RNA_coverage',completed_guides=len(done),total_guides=1196)
    for name in ['bam','bam_index','genome_fasta','genome_fasta_index']:
        if identity(config[name]) != config[name+'_identity']:raise ValueError('Input changed during coverage')
    results=[done[r['guide']] for r in rows]
    fields=list(dict.fromkeys(k for row in results for k in row))
    with (output/'all_design_guides_sequence_observability.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fields);w.writeheader();w.writerows(results)
    summary=dict(design_guides=1196,scope='company_merged_records_diagnostic',
                 uniquely_placed=sum(r['mapping_status']=='unique_exact_genomic_placement' for r in results),
                 targeting_unique=sum(r['source_control_type'] in ['targeting','positive_control'] and r['mapping_status']=='unique_exact_genomic_placement' for r in results),
                 targeting_with_own_assigned_spanning_CB_UB_families=sum(r['source_control_type'] in ['targeting','positive_control'] and r.get('own_guide_distinct_CB_UB_families_spanning',0)>0 for r in results),
                 capped_sites=sum(r.get('scan_capped',False) for r in results),
                 exact_cut_position_or_PAM_assumed=False,DNA_editing_efficiency_assessed=False,
                 variants_called=False,functional_KO_efficiency_assessed=False,
                 limitation='RNA observability over the full 20nt reference spacer, not KO efficiency; '
                 'no coverage is not an unedited call; no invented capture labels; CB/UB family counts are conditional on merged barcode records.')
    atomic(output/'summary.json',summary)
    atomic(output/'manifest.json',{p.name:dict(bytes=p.stat().st_size,sha256=digest(p)) for p in [output/'summary.json',output/'all_design_guides_sequence_observability.csv',ledger]})
    progress('complete_pending_interpretation',**summary)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--config',required=True);args=parser.parse_args()
    config=json.loads(Path(args.config).read_text());out=Path(config['output']);out.mkdir(parents=True,exist_ok=True)
    lock=out/'active.lock';fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY)
    os.write(fd,json.dumps(dict(pid=os.getpid(),uid=os.getuid(),config=args.config)).encode());os.close(fd)
    def progress(stage,**kw):
        value=dict(stage=stage,at=time.strftime('%Y-%m-%dT%H:%M:%S%z'),pid=os.getpid(),**kw)
        atomic(out/'status.json',value);print(json.dumps(value),flush=True)
    try:run(config,progress)
    except BaseException as e:progress('failed_pending_recovery',error=repr(e));raise
    finally:lock.unlink()


if __name__=='__main__':main()
