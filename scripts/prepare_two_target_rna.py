"""Prepare whole-gene counts and frozen assignments for the two-target RNA pilot."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import socket
import time

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from compare_established_perturbseq import load_counts, sha_file


def write_sparse_counts(counts, output):
    if counts.nnz >= np.iinfo(np.int32).max:
        raise ValueError('This prepared input exceeds R Matrix integer indexing')
    counts.sort_indices()
    arrays = [('counts.data.i32.gz', counts.data), ('counts.indices.i32.gz', counts.indices),
              ('counts.indptr.i32.gz', counts.indptr)]
    for name, values in arrays:
        partial = output / (name + '.partial')
        with gzip.open(partial, 'wb', compresslevel=3) as stream:
            stream.write(memoryview(values.astype('<i4', copy=False)))
        partial.replace(output / name)
    return [name for name, _ in arrays]


def prepare(config):
    if os.getuid() != 2601 or socket.gethostname() != 'rhesusbase':
        raise ValueError('The laboratory count preparation belongs on rhesusbase')
    os.sched_setaffinity(0, sorted(os.sched_getaffinity(0))[:2])
    output = Path(config['output'])
    data_root = Path('/home/user/data/wusy/260918_perturbseq').resolve()
    if not output.resolve().is_relative_to(data_root):
        raise ValueError('Prepared laboratory data must stay on the adopted data filesystem')
    output.mkdir(parents=True, exist_ok=True)
    proof_path = output / 'input_proof.json'
    contract_sha = hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()
    if proof_path.is_file():
        saved = json.loads(proof_path.read_text())
        if saved['contract_sha256'] != contract_sha:
            raise ValueError('Completed input belongs to a different contract')
        for name, item in saved['files'].items():
            if sha_file(output / name) != item['sha256']:
                raise ValueError('Completed prepared input differs: ' + name)
        print(json.dumps({'stage': 'reuse_verified_input', 'records': saved['records']}), flush=True)
        return saved
    if any(output.iterdir()):
        raise ValueError('Partial input directory is retained for diagnosis')

    def progress(stage, **values):
        print(json.dumps({'stage': stage, 'time': time.strftime('%Y-%m-%dT%H:%M:%S%z'), **values}), flush=True)
        if config.get('phase_file'):
            Path(config['phase_file']).write_text(stage + '\n')

    for name, expected in config['metadata_sha256'].items():
        if sha_file(name) != expected:
            raise ValueError('Accepted metadata differs: ' + name)
    qc = pd.read_csv(config['qc_table'], keep_default_na=False)
    assigned = pd.read_csv(config['assigned_cells'], keep_default_na=False)
    catalogue = pd.read_csv(config['catalogue'], keep_default_na=False)
    if not qc.cell_barcode.is_unique or not assigned.cell_barcode.is_unique or not catalogue.guide.is_unique:
        raise ValueError('Ambiguous count, assignment or design identity')
    indexed = qc.set_index('cell_barcode')
    frozen = indexed.reindex(assigned.cell_barcode)
    if not frozen.guide_id.tolist() == assigned.guide_id.tolist():
        raise ValueError('Frozen assignments differ from the accepted QC table')
    if not frozen.qc_pass_stage1.all() or not frozen.guide_num_features.eq(1).all():
        raise ValueError('The assigned input includes an unqualified record')
    if not np.array_equal(frozen.total_counts.to_numpy(), assigned.total_counts.to_numpy()):
        raise ValueError('Frozen full-RNA totals differ')
    assigned = assigned[['cell_barcode', 'guide_id', 'target', 'role', 'total_counts']].copy()
    for column in ('n_genes_by_counts', 'pct_counts_mt', 'raw_guide_umis', 'guide_umis'):
        assigned[column] = frozen[column].to_numpy()
    no_guide = qc.loc[qc.qc_pass_stage1 & qc.guide_status.eq('no_confident_call') & qc.raw_guide_umis.eq(0)].copy()
    no_guide['guide_id'] = ''
    no_guide['target'] = 'NO_GUIDE_DETECTED'
    no_guide['role'] = 'no_guide_detected'
    cells = pd.concat([assigned, no_guide[assigned.columns]], ignore_index=True)
    if len(assigned) != 72107 or len(no_guide) != 9094 or not cells.cell_barcode.is_unique:
        raise ValueError('The declared assigned and guide-free cohorts differ')
    cells_path = output / 'cells_for_count_read.csv.gz'
    cells.to_csv(cells_path, index=False)
    count_config = dict(config['count_input'], cells=str(cells_path), catalogue=config['catalogue'],
                        expected_cells=len(cells), allow_no_guide_records=True)
    with threadpool_limits(limits=2):
        counts, cells, genes, source_proof = load_counts(count_config, progress)
    annotation = catalogue.set_index('guide').reindex(cells.guide_id)
    cells['target_gene_id'] = annotation.target_gene_id.fillna('').to_numpy()
    cells['target_gene_name'] = annotation.target_gene_name.fillna('').to_numpy()
    cells['analysis_group'] = cells.target_gene_id
    cells.loc[cells.role.eq('non_targeting_control'), 'analysis_group'] = 'NTC'
    cells.loc[cells.role.eq('intergenic_control'), 'analysis_group'] = 'INTERGENIC'
    cells.loc[cells.role.eq('no_guide_detected'), 'analysis_group'] = 'NO_GUIDE'
    if cells.analysis_group.eq('').any():
        raise ValueError('Missing analysis-group identity')
    target_mask = cells.target_gene_name.isin(['TYMSOS', 'UBA1'])
    cells['pilot_record'] = target_mask | cells.role.isin([
        'non_targeting_control', 'intergenic_control', 'no_guide_detected'])
    expected = {'non_targeting_control': 19342, 'intergenic_control': 195, 'no_guide_detected': 9094}
    for role, number in expected.items():
        if int(cells.role.eq(role).sum()) != number:
            raise ValueError('Unexpected cohort size for ' + role)
    if int(target_mask.sum()) != 1219 or int(cells.pilot_record.sum()) != 29850:
        raise ValueError('The two-target pilot population differs')
    genes['seurat_key'] = [f'RNA{index:06d}' for index in range(len(genes))]
    genes.to_csv(output / 'genes.csv.gz', index=False)
    cells.to_csv(output / 'cells.csv.gz', index=False)
    catalogue[['guide', 'target_gene_id', 'target_gene_name', 'source_control_type']].to_csv(
        output / 'guide_design.csv', index=False)

    progress('write_RNA_sparse_arrays', records=len(cells), RNAs=len(genes), nnz=int(counts.nnz))
    # A cells-by-genes CSR is the same sparse storage as its genes-by-cells CSC
    # transpose. Gzipped binary vectors let both installed R versions read the
    # exact integer counts without a second HDF5 dependency or text conversion.
    pd.DataFrame([{'genes': len(genes), 'cells': len(cells), 'nnz': int(counts.nnz)}]).to_csv(
        output / 'matrix_dimensions.tsv', sep='\t', index=False)
    names = write_sparse_counts(counts, output) + [
        'cells.csv.gz', 'genes.csv.gz', 'guide_design.csv', 'matrix_dimensions.tsv']
    proof = dict(contract_sha256=contract_sha, records=len(cells), RNAs=len(genes),
                 assigned_records=len(assigned), guide_free_records=len(no_guide), pilot_records=29850,
                 matrix_shape_genes_by_cells=[len(genes), len(cells)], nnz=int(counts.nnz),
                 sparse_storage='R dgCMatrix slots: little-endian signed int32 counts, indices, pointers',
                 input_scope='company_merged_barcode_records_not_capture_resolved',
                 source=source_proof, metadata_sha256=config['metadata_sha256'],
                 files={name: {'bytes': (output/name).stat().st_size, 'sha256': sha_file(output/name)} for name in names})
    proof_path.write_text(json.dumps(proof, indent=2) + '\n')
    progress('whole_gene_input_complete', records=len(cells), RNAs=len(genes))
    return proof


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    arguments = parser.parse_args()
    prepare(json.loads(arguments.config.read_text()))
