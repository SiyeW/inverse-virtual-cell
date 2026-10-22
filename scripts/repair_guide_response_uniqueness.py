"""Restore only incorrectly masked scalar responses in a distinct output.

The original sparse counts, full negative-control tests, and prior outputs
are untouched. Input reference names are checked across ALL stored RNA rows.
"""
from pathlib import Path
from collections import Counter
from datetime import datetime
import argparse,json,os,subprocess,time
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from threadpoolctl import threadpool_limits
from analyze_expression_structure import sha256
from copy_verified_delivery import atomic_json
from run_chip01_qc_supplement import write_response_agreement


def run(args):
    config=json.loads(args.config.read_text());parent=Path(config['output'])
    output=parent/'guide_agreement_global_ids_v2'
    if output.exists():raise FileExistsError('Distinct correction output required; do not overwrite')
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=Path(__file__).resolve().parents[1],text=True).strip()
    if commit!=args.expected_commit:raise ValueError('Correction must use the pinned origin-pushed commit')
    state=json.loads((parent/'status.json').read_text())
    if state['stage']!='complete' or (parent/'owner.lock').exists():raise ValueError('Parent computation not closed')
    selected=['tables/guide_response_profiles.npz','tables/guide_response_profiles_metadata.csv','tables/guide_coverage.csv']
    original={rel:sha256(parent/rel) for rel in selected}
    genes=pd.read_csv(config['gene_statistics'],index_col=0)
    names=genes.gene_name.astype('string').fillna('').to_numpy(dtype=str)
    global_counts=Counter(names)
    prior=np.load(parent/'tables/guide_response_profiles.npz')
    profiles=prior['values'].copy();before_profiles=profiles.copy()
    features=prior['feature_indices'];feature_ids=prior['gene_ids']
    if not np.array_equal(genes.index.to_numpy(dtype=str)[features],feature_ids):raise ValueError('Frozen feature IDs differ')
    records=pd.read_csv(parent/'tables/guide_response_profiles_metadata.csv')
    if not np.array_equal(records.guide.to_numpy(dtype=str),prior['guides']):raise ValueError('Guide order differs')
    bad=records.direct_target_feature_masked.to_numpy(dtype=bool)&np.array([global_counts[t]!=1 for t in records.target])
    if not bad.any():raise ValueError('No demonstrated global-ambiguity correction required')
    source=Path(config['source']);source_before=source.stat();started=time.monotonic()
    a=ad.read_h5ad(source,backed='r');all_obs=a.obs.copy();var=a.var.copy()
    if not np.array_equal(var.index.astype(str).to_numpy(),genes.index.astype(str).to_numpy()):raise ValueError('Source RNA order differs')
    obs=all_obs.loc[all_obs.main_analysis_stage1].copy()
    needed=all_obs.main_analysis_stage1 & all_obs.guide_id.astype(str).isin(records.loc[bad,'guide'])
    rows=np.flatnonzero(needed.to_numpy());row_meta=all_obs.iloc[rows]
    positions={target:int(np.flatnonzero(names[features]==target)[0]) for target in records.loc[bad,'target']}
    columns=np.array(sorted({int(features[p]) for p in positions.values()}))
    values=a[rows,columns].X
    values=values.toarray() if sparse.issparse(values) else np.asarray(values)
    a.file.close()
    reference_total=float(json.loads((parent/'guide_agreement_summary.json').read_text())['normalization_target'])
    normalized=values*(reference_total/row_meta.total_counts.to_numpy(dtype=float))[:,None]
    changes=[]
    for i in np.flatnonzero(bad):
        row=records.iloc[i];mask=row_meta.guide_id.astype(str).eq(row.guide).to_numpy()
        if int(mask.sum())!=int(row.cells):raise ValueError('Source guide cell identity differs')
        local=positions[row.target];global_index=int(features[local]);column=int(np.flatnonzero(columns==global_index)[0])
        mean=float(normalized[mask,column].mean())
        value=(mean-float(genes.iloc[global_index].core_control_normalized_mean))/float(genes.iloc[global_index].core_control_normalized_sd)
        changes.append(dict(guide=row.guide,target=row.target,feature_id=str(genes.index[global_index]),
                            global_name_matches=global_counts[row.target],prior_value=float(profiles[i,local]),restored_value=float(value)))
        profiles[i,local]=value;records.loc[i,'direct_target_feature_masked']=False
    changed=np.argwhere(profiles!=before_profiles)
    if len(changed)>int(bad.sum()):raise ValueError('Unexpected wider profile changes')
    output.mkdir();(output/'tables').mkdir();(output/'figures').mkdir()
    lock=output/'owner.lock';lock.write_text(json.dumps(dict(pid=os.getpid(),code_commit=commit)))
    try:
        import matplotlib
        matplotlib.use('Agg')
        coverage=pd.read_csv(parent/'tables/guide_coverage.csv')
        with threadpool_limits(config['threads']):
            summary=write_response_agreement(records,profiles,obs,coverage,features,reference_total,config,output)
        after=source.stat()
        if (source_before.st_size,source_before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns):raise ValueError('Source identity changed')
        for rel,expected in original.items():
            if sha256(parent/rel)!=expected:raise ValueError('Prior artifact changed')
        result=dict(code_commit=commit,parent_commit=state['code_commit'],prior_artifact_sha256=original,
            corrected_profiles=int(bad.sum()),changed_scalar_responses=len(changed),changes=changes,
            source_modified=False,fresh_source_sha_recomputed=False,parent_verified_source_sha256=config['source_sha256'],
            negative_control_results_recomputed=False,mask_uniqueness_scope='all stored RNA reference rows',
            summary=summary,elapsed_seconds=time.monotonic()-started)
        atomic_json(output/'correction_receipt.json',result)
        manifest={p.relative_to(output).as_posix():dict(bytes=p.stat().st_size,sha256=sha256(p))
                  for p in output.rglob('*') if p.is_file() and p.name!='owner.lock'}
        atomic_json(output/'output_manifest.json',manifest)
        atomic_json(output/'status.json',dict(stage='complete',time=datetime.now().astimezone().isoformat(),
            code_commit=commit,source_modified=False,corrected_profiles=int(bad.sum()),pid=os.getpid()))
        print(json.dumps(result),flush=True)
    finally:lock.unlink()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);p.add_argument('--expected-commit',required=True)
    run(p.parse_args())
