"""Inspect fixed-size energy sensitivity and existing control covariates."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
from compare_energy_fixed_n import fixed_rows
from guide_rna_energy import digest


def render(folder,points,clusters,output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    folder=Path(folder);output=Path(output);output.mkdir(exist_ok=True)
    a=pd.read_csv(folder/'same_cohort_all_cells.csv');b=pd.read_csv(folder/'fixed50.csv')
    if len(a)!=479 or a.guide.tolist()!=b.guide.tolist() or not b.n_each.eq(50).all():raise ValueError('Comparison cohort differs')
    for f in [a,b]:
        _,q,_,_=multipletests(f.pvalue,alpha=.05,method='fdr_bh')
        if not np.allclose(q,f.pvalue_adj,atol=1e-12):raise ValueError('Inconsistent within-cohort BH')
    if hashlib.sha256(Path(points).read_bytes()).hexdigest()!='74a9d0e4bd0854fe7a8b0e3b3ff0e507025bda37cca03a7797c4c3f1fd30395a':raise ValueError('Wrong point axes')
    if hashlib.sha256(Path(clusters).read_bytes()).hexdigest()!='74a7737f34f01d1ff3bce872c3a575afcc7f09f87242f09f11b4bd2e585de990':raise ValueError('Wrong coarse cluster source')
    with np.load(points,allow_pickle=False) as z:bc=z['cell_barcode'];labels=z['guide'];roles=z['role'];totals=z['total_RNA'];xy=z['UMAP']
    c=pd.read_csv(clusters).set_index('Barcode',verify_integrity=True)['Cluster'].reindex(bc)
    assert c.notna().all();c=c.to_numpy(dtype=int);categories=np.unique(c)
    ntc_pool=np.flatnonzero(roles=='non_targeting_control')
    ranked=np.array(sorted(ntc_pool,key=lambda i:hashlib.sha256(f'20261030|NTC_reference|{bc[i]}'.encode()).digest()))
    annotations=[]
    for r in b.itertuples():
        query=fixed_rows(np.flatnonzero(labels==r.guide),bc,r.guide,50,20261103);ref=ranked[labels[ranked]!=r.guide][:50]
        assert digest(dict(query=bc[query].tolist(),reference=bc[ref].tolist()))==r.sample_SHA256
        pq=np.array([(c[query]==k).mean() for k in categories]);pr=np.array([(c[ref]==k).mean() for k in categories])
        annotations.append(dict(guide=r.guide,role=r.role,query_median_UMI=float(np.median(totals[query])),reference_median_UMI=float(np.median(totals[ref])),
            absolute_log2_median_UMI_ratio=float(abs(np.log2(np.median(totals[query])/np.median(totals[ref])))),
            coarse_composition_total_variation=float(.5*np.abs(pq-pr).sum()),query_low_RNA_island_records=int((xy[query,0]>8).sum()),
            reference_low_RNA_island_records=int((xy[ref,0]>8).sum())))
    cov=pd.DataFrame(annotations).merge(b[['guide','energy_unbiased','pvalue','pvalue_adj']],on='guide',validate='one_to_one');cov.to_csv(output/'fixed50_covariate_diagnostics.csv',index=False)
    nc=cov.loc[cov.role.eq('non_targeting_control')];sig=nc.pvalue.lt(.05)
    diag=dict(groups=479,sample_identity_matches=True,coarse_groups_not_cell_types=True,covariates_not_regressed_or_filtered=True,
        NTC_energy_vs_abs_depth_difference_spearman=float(spearmanr(nc.energy_unbiased,nc.absolute_log2_median_UMI_ratio).statistic),
        NTC_energy_vs_coarse_composition_spearman=float(spearmanr(nc.energy_unbiased,nc.coarse_composition_total_variation).statistic),
        NTC_raw_signal_count=int(sig.sum()),NTC_groups=len(nc),
        NTC_raw_signal_median_depth_discrepancy=float(nc.loc[sig,'absolute_log2_median_UMI_ratio'].median()),
        NTC_other_median_depth_discrepancy=float(nc.loc[~sig,'absolute_log2_median_UMI_ratio'].median()),
        NTC_raw_signal_median_coarse_composition_TV=float(nc.loc[sig,'coarse_composition_total_variation'].median()),
        NTC_other_median_coarse_composition_TV=float(nc.loc[~sig,'coarse_composition_total_variation'].median()),
        causal_source_not_identified=True,all_original_results_preserved=True)
    (output/'covariate_summary.json').write_text(json.dumps(diag,indent=2))
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','DejaVu Sans'],'axes.unicode_minus':False,'font.size':11,'svg.fonttype':'none'})
    fig,axes=plt.subplots(1,2,figsize=(9.2,5.8),sharey=True);fig.subplots_adjust(left=.10,right=.97,bottom=.23,top=.82,wspace=.25)
    for ax,f,title in zip(axes,[a,b],['每条sgRNA使用全部合格细胞','每条sgRNA与参照各取50个细胞']):
        for i,(mask,color,label) in enumerate([(f.role.eq('non_targeting_control'),'#3d799d','NTC（99条）'),(f.role.isin(['targeting','positive_control']),'#c58d59','靶向sgRNA（379条）')]):
            g=f.loc[mask];y=g.energy_unbiased.to_numpy();jitter=np.array([int.from_bytes(hashlib.sha256(x.encode()).digest()[:4],'little')/2**32-.5 for x in g.guide])*.48
            ax.scatter(i+jitter,y,s=13,color=color,alpha=.40,linewidths=0,zorder=2)
            bp=ax.boxplot([y],positions=[i],widths=.52,showfliers=False,patch_artist=True,manage_ticks=False,
                medianprops=dict(color='#222222',linewidth=1.5),boxprops=dict(facecolor='none',edgecolor=color,linewidth=1.3),whiskerprops=dict(color=color),capprops=dict(color=color))
        ax.axhline(0,color='#8b9399',lw=.7);ax.set_xticks([0,1],['NTC（99条）','靶向sgRNA（379条）']);ax.set_title(title,fontsize=11,pad=12);ax.spines[['top','right']].set_visible(False)
    axes[0].set_ylabel('整体RNA差异（E-distance无偏估计）')
    fig.suptitle('同一批sgRNA：细胞数统一前后的RNA差异',fontsize=16,y=.96)
    fig.text(.5,.09,'每点一条sgRNA，黑线为中位数；两侧比较的是完全相同的sgRNA。\n少于50个细胞的sgRNA仍保留原结果，不进入本次对照；负值不表示RNA下调。',ha='center',va='center',fontsize=9.3,linespacing=1.8)
    for ext in ['png','svg']:fig.savefig(output/f'同细胞数RNA差异.{ext}',dpi=300,facecolor='white')
    plt.close(fig)
    print(json.dumps(diag,ensure_ascii=False))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--folder',type=Path,required=True);p.add_argument('--points',type=Path,required=True);p.add_argument('--clusters',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();render(a.folder,a.points,a.clusters,a.output)
