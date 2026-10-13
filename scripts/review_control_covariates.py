"""NTC-only sensitivity diagnostics; never residualize treatment expression.

Association with QC covariates is not evidence of causality. This script
does not select core controls, remove cells, or finalize the analysis set.
"""
import argparse
import json
from pathlib import Path
import subprocess
import time

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control
from threadpoolctl import threadpool_limits

from analyze_expression_structure import atomic_json, sha256
from diagnose_control_guides import centroid_null_test


def residualize(y, covariates):
    y, covariates = np.asarray(y,dtype=float),np.asarray(covariates,dtype=float)
    if y.ndim != 2 or covariates.ndim != 2 or len(y) != len(covariates):
        raise ValueError('Expected aligned 2D arrays')
    if not np.isfinite(y).all() or not np.isfinite(covariates).all():
        raise ValueError('Nonfinite covariates or coordinates')
    z = (covariates-covariates.mean(0))/np.maximum(covariates.std(0),1e-8)
    design = np.column_stack([np.ones(len(y)),z])
    coefficients,_,rank,_ = np.linalg.lstsq(design,y,rcond=None)
    return y-design@coefficients,int(rank)


def screen(values, guides, draws=1000):
    records=[]
    for i,guide in enumerate(sorted(set(guides))):
        mask=guides==guide
        row={'guide':guide,'cells':int(mask.sum())}
        if mask.sum()>=20:
            row['centroid_distance'],row['label_permutation_p']=centroid_null_test(values[mask],values[~mask],draws=draws,seed=i)
        records.append(row)
    result=pd.DataFrame(records).set_index('guide')
    tested=result.label_permutation_p.notna()
    result.loc[tested,'bh_q']=false_discovery_control(result.loc[tested,'label_permutation_p'].to_numpy())
    return result


def run(root, output, draws=1000):
    if draws<1000:
        raise ValueError('Control review requires at least 1000 permutations')
    if output.exists():
        raise FileExistsError('Preserve existing review')
    output.mkdir(parents=True)
    started=time.monotonic()
    state_path=output/'status.json'
    atomic_json(state_path,{'stage':'verify_input'})
    try:
        expr=root/'expression_diagnostics_v2'
        doublet=root/'doublet_diagnostics_full'
        original=root/'ntc_consistency'
        es=json.loads((expr/'expression_diagnostics_status.json').read_text())
        ds=json.loads((doublet/'status.json').read_text())
        ns=json.loads((original/'summary.json').read_text())
        checks=[(expr/'pca_coordinates.csv',es['output_sha256']['pca_coordinates.csv']),
                (doublet/'cell_doublet_diagnostics.csv.gz',ds['output_sha256']['cell_doublet_diagnostics.csv.gz']),
                (original/'ntc_guide_consistency.csv',ns['output_sha256'])]
        for path,digest in checks:
            if sha256(path)!=digest:
                raise ValueError(f'Input mismatch: {path.name}')
        obs=pd.read_csv(doublet/'cell_doublet_diagnostics.csv.gz',index_col=0)
        obs=obs[obs.is_control].copy()
        pc=pd.read_csv(expr/'pca_coordinates.csv',index_col=0).loc[obs.index].to_numpy()
        pc=(pc-pc.mean(0))/np.maximum(pc.std(0),1e-8)
        guides=obs.guide_id.to_numpy()
        covariates=np.column_stack([np.log1p(obs.total_counts),obs.pct_counts_mt,obs.doublet_score])
        atomic_json(state_path,{'stage':'quality_association_sensitivity','control_cells':len(obs)})
        with threadpool_limits(limits=1):
            residual,rank=residualize(pc,covariates)
            adjusted=screen(residual,guides,draws).add_prefix('qc_adjusted_')
            threshold=float(obs.doublet_score.quantile(.99))
            keep=obs.doublet_score.to_numpy()<=threshold
            trimmed=screen(pc[keep],guides[keep],draws).add_prefix('highscore_tail_excluded_')
            retested=screen(pc,guides,draws).add_prefix('baseline_retested_')
        table=pd.read_csv(original/'ntc_guide_consistency.csv').set_index('guide')
        table=table.join(adjusted).join(trimmed).join(retested)
        grouped=obs.groupby('guide_id',observed=True)
        table['median_doublet_score']=grouped.doublet_score.median()
        table['tail_sensitivity_cells_removed']=table.cells-table.highscore_tail_excluded_cells
        table.to_csv(output/'ntc_covariate_sensitivity.csv')
        flagged=table.bh_q<.05
        summary={'stage':'complete_diagnostics_only','elapsed_seconds':time.monotonic()-started,
                 'code_commit':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
                 'control_cells':len(obs),'permutations':draws,
                 'original_flagged_guides':int(flagged.sum()),
                 'baseline_retested_flagged_guides':int((table.baseline_retested_bh_q<.05).sum()),
                 'qc_adjusted_flagged_guides':int((table.qc_adjusted_bh_q<.05).sum()),
                 'tail_excluded_flagged_guides':int((table.highscore_tail_excluded_bh_q<.05).sum()),
                 'original_flags_persist_after_both':int((flagged & (table.qc_adjusted_bh_q<.05) & (table.highscore_tail_excluded_bh_q<.05)).sum()),
                 'ntc_tail_score_cutoff':threshold,'ntc_tail_cells':int((~keep).sum()),
                 'covariate_design_rank':rank,'filtering_applied_to_source':False,
                 'treatment_expression_residualized':False,
                 'interpretation':'NTC-only quality association sensitivity, not causal attribution or final core-control selection',
                 'output_sha256':sha256(output/'ntc_covariate_sensitivity.csv')}
        atomic_json(state_path,summary)
    except Exception as exc:
        atomic_json(state_path,{'stage':'failed','error':repr(exc)})
        raise


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--permutations',type=int,default=1000)
    a=p.parse_args()
    run(a.root,a.output,a.permutations)
