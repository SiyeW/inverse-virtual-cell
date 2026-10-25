"""Describe existing expression-cluster coverage; clusters are not capture or cell-type labels."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def analyze(cells, records, clusters):
    if cells.index.has_duplicates or clusters.Barcode.duplicated().any() or records.guide.duplicated().any():
        raise ValueError('Nonunique cell/guide/cluster identity')
    if clusters[['Barcode','Cluster']].isna().any().any() or not cells.index.isin(clusters.Barcode).all():
        raise ValueError('Missing cluster labels for analysis cells')
    # Reuse design identity/roles only. Old catalogue population counts are not inputs.
    lookup=records.loc[:,['guide','source_control_type']].set_index('guide')
    if not cells.guide_id.isin(lookup.index).all():raise ValueError('Unknown guide identity')
    labels=clusters.set_index('Barcode').Cluster.reindex(cells.index)
    roles=cells.guide_id.map(lookup.source_control_type)
    joined=pd.DataFrame(dict(guide=cells.guide_id,role=roles,cluster=labels),index=cells.index)
    table=pd.crosstab(joined.cluster,joined.role)
    table=table.reindex(sorted(table.index,key=lambda x:(0,int(x)) if str(x).isdigit() else (1,str(x))))
    reference=table['non_targeting_control']/table['non_targeting_control'].sum()
    guide_counts=pd.crosstab(joined.guide,joined.cluster).reindex(columns=table.index,fill_value=0)
    ntc=roles.eq('non_targeting_control')
    ntc_counts=table['non_targeting_control'].to_numpy(float)
    comparisons=[]
    for guide,row in guide_counts.iterrows():
        role=lookup.loc[guide,'source_control_type'];n=row.sum()
        ref=(ntc_counts-row.to_numpy())/(ntc_counts.sum()-n) if role=='non_targeting_control' else reference.to_numpy()
        fraction=row.to_numpy()/n
        comparisons.append(dict(guide=guide,role=role,cells=int(n),
            expression_cluster_total_variation=float(.5*np.abs(fraction-ref).sum()),
            reference_overlap_cells=0))
    coverage=table.copy();coverage['total_cells']=coverage.sum(axis=1)
    coverage['ntc_fraction_within_cluster']=coverage['non_targeting_control']/coverage.total_cells
    coverage['ntc_reference_fraction']=reference
    role_tv={}
    for role in table:
        p=table[role]/table[role].sum()
        role_tv[role]=float(.5*np.abs(p-reference).sum())
    summary=dict(analysis_cells=len(cells),source_cluster_rows=len(clusters),matched_analysis_cells=int(labels.notna().sum()),
        expression_clusters=len(table),ntc_cells=int(ntc.sum()),clusters_with_ntc=int((ntc_counts>0).sum()),
        cluster_ntc_fraction_range=[float(coverage.ntc_fraction_within_cluster.min()),float(coverage.ntc_fraction_within_cluster.max())],
        ntc_reference_cluster_counts=table['non_targeting_control'].astype(int).to_dict(),
        role_composition_tv_against_ntc=role_tv,clusters_are_validated_cell_types=False,
        clusters_are_capture_labels=False,coarse_coverage_not_conditional_exchangeability=True,
        new_clustering_or_model_fit=False,causal_effects_or_p_values=False,
        guide_counts_recomputed_from_analysis_cells=True,catalogue_roles_only=True,
        observed_guide_groups=len(guide_counts))
    return coverage,pd.DataFrame(comparisons),summary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['cells','guide-metadata','clusters','config','output']:p.add_argument('--'+name,required=True)
    a=p.parse_args();out=Path(a.output)
    if out.exists():raise ValueError('Unused output directory required')
    specs=json.loads(Path(a.config).read_text());paths=dict(cells=Path(a.cells),guide_metadata=Path(a.guide_metadata),clusters=Path(a.clusters))
    for name,path in paths.items():
        if hashlib.sha256(path.read_bytes()).hexdigest()!=specs[name]['sha256']:raise ValueError('Input content differs: '+name)
    cells=pd.read_csv(paths['cells'],index_col=0)
    records=pd.read_csv(paths['guide_metadata']);clusters=pd.read_csv(paths['clusters'],dtype={'Barcode':str,'Cluster':str})
    coverage,guides,summary=analyze(cells,records,clusters)
    summary.update(code_commit=specs['code_commit'],helper_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),sources=specs)
    out.mkdir();coverage.to_csv(out/'cluster_coverage.csv',index_label='expression_cluster');guides.to_csv(out/'guide_cluster_composition.csv',index=False)
    (out/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf8')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
