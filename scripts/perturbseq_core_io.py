"""Small immutable-input contract shared by the adopted laboratory workflow."""
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from perturbseq_core_effects import TYPES, qc_flow


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        while block := stream.read(4 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf8')
    temporary.replace(path)


def identity(path):
    stat = Path(path).stat()
    return dict(bytes=stat.st_size, mtime_ns=stat.st_mtime_ns, device=stat.st_dev, inode=stat.st_ino)


def revision(expected, smoke=False):
    repo = Path(__file__).resolve().parents[1]
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip()
    dirty = bool(subprocess.check_output(['git', 'status', '--porcelain'], cwd=repo, text=True).strip())
    if commit != expected or (dirty and not smoke):
        raise ValueError('Use the configured fixed commit and a clean checkout for a formal run')
    return dict(commit=commit, working_tree_dirty=dirty, engineering_smoke=bool(smoke))


def load_config(path):
    path = Path(path).resolve()
    config = json.loads(path.read_text(encoding='utf8'))
    if config.get('schema_version') != 1:
        raise ValueError('Expected core-analysis configuration schema_version=1')
    for key in ['run_id', 'chip_id']:
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,95}', config[key]):
            raise ValueError('Unsafe run/chip identifier: ' + key)
    for key in ['source', 'catalogue']:
        spec = config[key]
        if not re.fullmatch('[0-9a-f]{64}', spec['sha256']):
            raise ValueError('A verified SHA256 is required: ' + key)
        spec['path'] = str((path.parent / spec['path']).resolve())
    for key in ['output', 'allowed_output_parent']:
        config[key] = str((path.parent / config[key]).resolve())
    if config['features']['mode'] == 'frozen_ids':
        spec = config['features']
        spec['path'] = str((path.parent / spec['path']).resolve())
        if not re.fullmatch('[0-9a-f]{64}', spec['sha256']):
            raise ValueError('Frozen feature list needs its own SHA256')
    elif config['features']['mode'] != 'pooled_raw_mean':
        raise ValueError('features.mode must be frozen_ids or pooled_raw_mean')
    for name in ['min_guide_cells', 'matched_cells', 'matched_repeats', 'threads', 'smoke_prefix_rows']:
        value = config[name]
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError('Positive integer required: ' + name)
    if config['min_guide_cells'] < config['matched_cells']:
        raise ValueError('min_guide_cells must be at least matched_cells')
    if not isinstance(config['seed'], int) or config['seed'] < 0:
        raise ValueError('Nonnegative integer seed required')
    # The accepted pseudocount is not silently exposed as a tuning parameter.
    if config.get('pseudocount', 1) != 1:
        raise ValueError('This workflow fixes the RNA pseudocount at 1')
    config['pseudocount'] = 1.
    cards = config.get('max_target_cards')
    if cards is not None and (isinstance(cards, bool) or not isinstance(cards, int) or cards < 0):
        raise ValueError('max_target_cards must be null or a nonnegative integer')
    config['max_target_cards'] = cards
    scope = config.get('input_scope', 'independent_capture')
    if scope not in ['independent_capture', 'company_merged_records_diagnostic']:
        raise ValueError('Unknown input_scope')
    notice = config.get('report_notice', '')
    if scope == 'company_merged_records_diagnostic' and (not isinstance(notice, str) or not notice.strip()):
        raise ValueError('Merged diagnostic input requires an explicit report_notice')
    config['input_scope'] = scope
    return config


def validate_metadata(obs, var):
    required = {'guide_id', 'perturbed_gene', 'is_control', 'guide_status', 'qc_pass_stage1',
                'main_analysis_stage1', 'total_counts', 'n_genes_by_counts', 'pct_counts_mt',
                'raw_guide_umis', 'guide_umis'}
    if not required.issubset(obs) or 'gene_name' not in var:
        raise ValueError('Use the documented all-cell RNA/guide QC H5AD contract')
    if obs.index.has_duplicates or var.index.has_duplicates or obs.index.isna().any() or var.index.isna().any():
        raise ValueError('Unique nonmissing barcode and RNA feature IDs required')
    if 'feature_type' in var and not var.feature_type.astype(str).eq('Gene Expression').all():
        raise ValueError('RNA-only source required; guide features cannot enter RNA effects')
    for column in ['is_control', 'qc_pass_stage1', 'main_analysis_stage1']:
        if obs[column].isna().any() or not pd.api.types.is_bool_dtype(obs[column]):
            raise ValueError('Actual boolean metadata required, not strings: ' + column)
    allowed = {'single_guide', 'multi_guide', 'no_confident_call'}
    if not set(obs.guide_status.astype(str)).issubset(allowed):
        raise ValueError('Unknown guide assignment status')
    for column in ['total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'raw_guide_umis', 'guide_umis']:
        values = obs[column].to_numpy(dtype=float)
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError('Invalid cell metric: ' + column)
    qc_flow(obs)


def read_source(spec, prefix_rows=None):
    """Read bounded CSR spans, preserve raw integer counts and the full denominator."""
    path = Path(spec['path'])
    before = identity(path)
    if before['bytes'] != spec['bytes']:
        raise ValueError('Source byte size differs')
    view = ad.read_h5ad(path, backed='r')
    try:
        all_obs, var = view.obs.copy(), view.var.copy()
    finally:
        view.file.close()
    if [len(all_obs), len(var)] != spec['shape']:
        raise ValueError('Source dimensions differ')
    validate_metadata(all_obs, var)
    # Reject cheap metadata/shape mistakes before a potentially expensive full read.
    if prefix_rows is None and sha256(path) != spec['sha256']:
        raise ValueError('Source SHA256 differs; no automatic repair or rerun')
    if prefix_rows is not None:
        all_obs = all_obs.iloc[:prefix_rows].copy()
    selected = all_obs.main_analysis_stage1.to_numpy()
    rows = np.flatnonzero(selected)
    if not len(rows):
        raise ValueError('No cells satisfy the frozen RNA-QC/single-guide policy')
    parts = []
    with h5py.File(path, 'r') as handle:
        x = handle['X']
        if x.attrs.get('encoding-type') != 'csr_matrix' or list(x.attrs['shape']) != spec['shape']:
            raise ValueError('Sparse raw RNA counts in X, CSR encoding required')
        for start in range(0, len(rows), 128):
            picked = rows[start:start + 128]
            first, stop = int(picked[0]), int(picked[-1]) + 1
            ptr = x['indptr'][first:stop + 1]
            lo, hi = int(ptr[0]), int(ptr[-1])
            data = x['data'][lo:hi]
            if not np.isfinite(data).all() or (data < 0).any() or (data != np.floor(data)).any():
                raise ValueError('Nonnegative finite original integer counts required')
            block = sparse.csr_matrix((data, x['indices'][lo:hi], ptr - lo), shape=(stop-first, len(var)))
            if not block.has_canonical_format:
                raise ValueError('Canonical CSR with unique feature indices required')
            if not np.array_equal(np.asarray(block.sum(1)).ravel(), all_obs.total_counts.to_numpy()[first:stop]):
                raise ValueError('Raw RNA totals disagree with source metadata')
            parts.append(block[picked-first])
    counts = sparse.vstack(parts, format='csr')
    obs = all_obs.iloc[rows].copy()
    obs['source_row'] = rows
    if (obs.total_counts <= 0).any() or identity(path) != before:
        raise ValueError('Nonpositive selected RNA totals or changed source identity')
    return counts, obs, all_obs, var, dict(before, sha256=spec['sha256'],
        full_hash_verified=prefix_rows is None, source_prefix_rows=len(all_obs), source_modified=False)


def read_catalogue(spec, all_obs, control_target):
    path = Path(spec['path'])
    if sha256(path) != spec['sha256']:
        raise ValueError('Catalogue SHA256 differs')
    library = pd.read_csv(path, dtype=str, keep_default_na=False)
    required = {'guide', 'target_gene_name', 'source_control_type'}
    if not required.issubset(library) or library.guide.duplicated().any() or (library[list(required)] == '').any().any():
        raise ValueError('Complete unique authoritative guide catalogue required')
    if not set(library.source_control_type).issubset(TYPES):
        raise ValueError('Unrecognized source control role')
    negative = library.source_control_type.isin(TYPES[:2])
    if not np.array_equal(negative, library.target_gene_name.eq(control_target)):
        raise ValueError('Negative catalogue roles and configured coarse control target conflict')
    if library.groupby('target_gene_name').source_control_type.nunique().gt(1).drop(control_target, errors='ignore').any():
        raise ValueError('A biological target has conflicting source roles')
    cells = all_obs.loc[all_obs.guide_status.astype(str).eq('single_guide')]
    lookup = library.set_index('guide')
    if not set(cells.guide_id.astype(str)).issubset(lookup.index):
        raise ValueError('Observed guide absent from the catalogue; no name guessing')
    expected_control = cells.guide_id.astype(str).map(lookup.source_control_type).isin(TYPES[:2])
    expected_targets = cells.guide_id.astype(str).map(lookup.target_gene_name).where(~expected_control, 'ctrl')
    if not np.array_equal(cells.is_control, expected_control) or not np.array_equal(cells.perturbed_gene.astype(str), expected_targets):
        raise ValueError('Catalogue/source cell assignment identity conflict')
    for key, mask in [('single_guide_before_qc', all_obs.guide_status.astype(str).eq('single_guide')),
                      ('single_guide_after_stage1_qc', all_obs.main_analysis_stage1)]:
        counts = all_obs.loc[mask, 'guide_id'].astype(str).value_counts()
        library[key] = library.guide.map(counts).fillna(0).astype(int)
    library['source_role_validated_in_this_assay'] = False
    return library.sort_values('guide').reset_index(drop=True)


def select_features(counts, var, spec):
    if spec['mode'] == 'frozen_ids':
        if sha256(spec['path']) != spec['sha256']:
            raise ValueError('Frozen RNA feature list SHA256 differs')
        table = pd.read_csv(spec['path'], dtype=str, keep_default_na=False)
        if 'gene_id' not in table or table.gene_id.duplicated().any() or table.gene_id.eq('').any():
            raise ValueError('Frozen feature CSV needs unique gene_id values')
        indices = var.index.astype(str).get_indexer(table.gene_id)
        if (indices < 0).any() or (len(indices) > 1 and (np.diff(indices) <= 0).any()):
            raise ValueError('Frozen RNA feature IDs must match source order exactly')
    else:
        cutoff = float(spec['raw_mean_gt'])
        if not np.isfinite(cutoff) or cutoff < 0:
            raise ValueError('Nonnegative pooled raw-count mean threshold required')
        indices = np.flatnonzero(np.asarray(counts.mean(0)).ravel() > cutoff)
    if len(indices) < 2:
        raise ValueError('Fewer than two descriptive RNA features; do not tune for positive results')
    return indices
