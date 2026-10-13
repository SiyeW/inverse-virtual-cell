"""Bounded real-data plumbing test; outputs are NOT scientific conclusions."""
import argparse
import json
from pathlib import Path
import subprocess
import time

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from analyze_expression_structure import blocks, sha256, atomic_json
from perturbation_statistics import energy_permutation, differential_expression_block, bh_adjust


def run(source, diagnostics, output):
    if output.exists():
        raise FileExistsError('Do not overwrite a previous smoke')
    output.mkdir(parents=True)
    started = time.monotonic()
    state_path = output / 'status.json'
    atomic_json(state_path, {'stage': 'verify_input', 'scientific_result': False})
    try:
        cp = json.loads((diagnostics / 'pca_checkpoint.json').read_text(encoding='utf8'))
        if sha256(source) != cp['input_sha256']:
            raise ValueError('Counts input identity mismatch')
        dstate = json.loads((diagnostics / 'expression_diagnostics_status.json').read_text(encoding='utf8'))
        pc_file = diagnostics / 'pca_coordinates.csv'
        if sha256(pc_file) != dstate['output_sha256'][pc_file.name]:
            raise ValueError('PCA identity mismatch')
        a = ad.read_h5ad(source, backed='r')
        obs, var = a.obs.copy(), a.var.copy()
        a.file.close()
        candidates = obs[obs.main_analysis_stage1 & ~obs.is_control]
        sizes = candidates.perturbed_gene.astype(str).value_counts()
        targets = sorted(sizes[sizes >= 20].index[:3])
        if len(targets) != 3:
            raise ValueError('Smoke requires 3 adequately represented targets')
        rng = np.random.default_rng(0)
        groups = {'ctrl': np.sort(rng.choice(np.flatnonzero(obs.main_analysis_stage1 & obs.is_control),128,replace=False))}
        for target in targets:
            pool = np.flatnonzero(obs.main_analysis_stage1 & obs.perturbed_gene.astype(str).eq(target))
            groups[target] = np.sort(rng.choice(pool,min(30,len(pool)),replace=False))
        rows = np.sort(np.concatenate(list(groups.values())))
        columns = np.asarray(cp['feature_indices'][:64],dtype=int)
        pieces = []
        with h5py.File(source,'r') as handle:
            for _, x in blocks(handle,rows,len(var),block_size=16):
                pieces.append(x[:,columns])
        counts = sparse.vstack(pieces,format='csr')
        values = counts.toarray().astype(float) * (10000/obs.total_counts.to_numpy()[rows])[:,None]
        index = {int(row):i for i,row in enumerate(rows)}
        pc = pd.read_csv(pc_file,index_col=0).iloc[:,:20]
        control = values[[index[int(row)] for row in groups['ctrl']]]
        energy_rows, gene_rows = [], []
        for target in targets:
            x = values[[index[int(row)] for row in groups[target]]]
            result = differential_expression_block(x,control)
            table = pd.DataFrame(result,index=var.index[columns])
            table['bh_q_among_64_smoke_genes_only'] = bh_adjust(table.pvalue.to_numpy())
            table['target'] = target
            table['scientific_result'] = False
            gene_rows.append(table)
            effect,pvalue = energy_permutation(pc.loc[obs.index[groups[target]]].to_numpy(),
                                               pc.loc[obs.index[groups['ctrl']]].to_numpy(),
                                               permutations=99)
            energy_rows.append({'target':target,'energy_v_statistic':effect,
                                'permutation_p_99_smoke_only':pvalue,'scientific_result':False})
        pd.concat(gene_rows).to_csv(output/'engineering_de_rows.csv')
        pd.DataFrame(energy_rows).to_csv(output/'engineering_energy_rows.csv',index=False)
        if not all(np.isfinite(pd.DataFrame(r).select_dtypes('number').to_numpy()).all()
                   for r in gene_rows):
            raise ValueError('Nonfinite smoke output')
        atomic_json(state_path, {'stage':'complete_engineering_smoke_only', 'scientific_result':False,
                    'elapsed_seconds':time.monotonic()-started,'source_sha256':cp['input_sha256'],
                    'code_commit':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
                    'controls':128,'targets':len(targets),'target_cells_max':30,'genes':64,'permutations':99,
                    'counts_representation':'source integer counts; normalized CP10k only in smoke array',
                    'energy_space':'diagnostic PCA, NOT the paper energy representation',
                    'main_analysis_set_finalized':False,'source_modified':False,
                    'output_sha256':{p.name:sha256(p) for p in output.glob('*.csv')}})
    except Exception as exc:
        atomic_json(state_path, {'stage':'failed','error':repr(exc),'scientific_result':False})
        raise


if __name__ == '__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--diagnostics',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    run(a.source,a.diagnostics,a.output)
