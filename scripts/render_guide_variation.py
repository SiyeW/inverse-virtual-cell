"""Two requested variation figures from fixed, verified result tables."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

DATASETS=['chip01','company','kolf']
LABELS={'chip01':'我们的首芯片','company':'公司合并矩阵','kolf':'论文 KOLF Strong'}
COLORS={'chip01':'#397ca9','company':'#c68a45','kolf':'#439477'}
KINDS=['within_guide','same_target','different_target','NTC']
CAPTIONS=['同一条sgRNA的细胞\n随机分成两组','靶向同一基因的\n不同sgRNA细胞组','靶向不同基因的\nsgRNA细胞组','检测到不同NTC sgRNA的\n细胞组']


def short_guide(value):
    value=str(value)
    return value if len(value)<=33 else value[:17]+'…'+value[-13:]


def validate_points(root,n):
    result={}
    for name in DATASETS:
        p=pd.read_csv(root/name/'variation_points.csv');p=p.loc[p.cells_per_group.eq(n)]
        if set(p.comparison)!=set(KINDS) or p.mean_absolute_difference.isna().any() or (p.mean_absolute_difference<0).any():raise ValueError('Invalid variation groups')
        sets=[set(p.loc[p.comparison.eq(k),'pair_id']) for k in KINDS]
        if any(x!=sets[0] for x in sets) or p[['comparison','pair_id']].duplicated().any():raise ValueError('Variation categories must use identical paired cohort')
        result[name]=p
    return result


def paper_batch_diagnostic(root,metadata,out,n):
    meta=pd.read_csv(metadata);meta=meta.loc[meta.source_control_type.eq('targeting')].set_index('guide')
    if not meta.index.is_unique:raise ValueError('A public targeting guide has multiple batch labels')
    plan=pd.read_csv(root/'kolf/comparison_plan.csv')
    for col in ['guide_a','guide_b','other_a','other_b']:
        plan[col+'_batch']=plan[col].map(meta.batch)
        if plan[col+'_batch'].isna().any():raise ValueError('Missing known public batch')
    plan['query_same_batch']=plan.guide_a_batch.eq(plan.guide_b_batch)
    plan['different_same_batch']=plan.other_a_batch.eq(plan.other_b_batch)
    plan['exact_batch_pattern_match']=[sorted([a,b])==sorted([c,d]) for a,b,c,d in zip(plan.guide_a_batch,plan.guide_b_batch,plan.other_a_batch,plan.other_b_batch)]
    values=pd.read_csv(root/'kolf/variation_points.csv');values=values.loc[values.cells_per_group.eq(n)].pivot(index='pair_id',columns='comparison',values='mean_absolute_difference')
    if set(values.index)!=set(plan.index):raise ValueError('Batch diagnostic changed source pairing')
    match=plan.exact_batch_pattern_match.to_numpy();selected=values.loc[match]
    summary=dict(all_pairs=len(plan),query_cross_batch_fraction=float((~plan.query_same_batch).mean()),
        different_cross_batch_fraction=float((~plan.different_same_batch).mean()),exact_batch_matched_pairs=int(match.sum()),
        matched_same_median=float(selected.same_target.median()),matched_different_median=float(selected.different_target.median()),
        matched_within_median=float(selected.within_guide.median()),matched_difference_minus_same=float((selected.different_target-selected.same_target).median()),
        matched_fraction_same_closer=float((selected.same_target<selected.different_target).mean()),
        RNA_values_or_references_changed=False,metadata_sha256=hashlib.sha256(Path(metadata).read_bytes()).hexdigest())
    plan.to_csv(out/'paper_pair_batch_annotations.csv',index=False)
    (out/'paper_batch_diagnostic.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary


def render(root,out,paper_metadata):
    root=Path(root);out=Path(out);out.mkdir(parents=True,exist_ok=False)
    font=Path('C:/Windows/Fonts/msyh.ttc')
    if not font.exists():raise ValueError('Readable Chinese font unavailable')
    font_manager.fontManager.addfont(str(font))
    plt.rcParams.update({'font.family':font_manager.FontProperties(fname=font).get_name(),'font.size':11,'axes.unicode_minus':False,'axes.spines.top':False,'axes.spines.right':False})
    summary=json.loads((root/'summary.json').read_text());n=int(summary['primary_n']);points=validate_points(root,n)
    batch=paper_batch_diagnostic(root,Path(paper_metadata),out,n)
    fig,ax=plt.subplots(figsize=(13.3,6.6),layout='constrained');rows=[];ymax=0
    fig.get_layout_engine().set(rect=(0,.09,1,.91))
    for j,name in enumerate(DATASETS):
        for i,kind in enumerate(KINDS):
            values=points[name].loc[points[name].comparison.eq(kind),'mean_absolute_difference'].to_numpy()
            pos=i+(j-1)*.24;box=ax.boxplot([values],positions=[pos],widths=.19,patch_artist=True,showfliers=False,manage_ticks=False,
                medianprops=dict(color=COLORS[name],lw=1.8),whiskerprops=dict(color=COLORS[name]),capprops=dict(color=COLORS[name]))
            box['boxes'][0].set(facecolor=COLORS[name],alpha=.2,edgecolor=COLORS[name])
            jitter=np.random.default_rng(90210+j*10+i).uniform(-.065,.065,len(values))
            ax.scatter(pos+jitter,values,s=9,alpha=.25,color=COLORS[name],edgecolors='none',zorder=2)
            ymax=max(ymax,float(values.max()));rows.append(dict(dataset=name,comparison=kind,pairs=len(values),median=float(np.median(values)),p10=float(np.quantile(values,.1)),p90=float(np.quantile(values,.9))))
        count=int(points[name].comparison.eq('same_target').sum())
        ax.plot([],[],marker='s',ls='',color=COLORS[name],label=f'{LABELS[name]}（{count}组配对）')
    ax.set(xticks=np.arange(4),xticklabels=CAPTIONS,xlim=(-.55,3.55),ylim=(0,ymax*1.1),
           ylabel='两组RNA变化的平均绝对差\n（log2单位；越低越接近）',title=f'每组{n}个细胞｜共同2,272项RNA作为基准')
    ax.legend(frameon=False,loc='upper right',fontsize=10)
    ax.grid(axis='y',alpha=.15,zorder=0)
    fig.text(.5,.025,'每点为一组配对8次抽样的中位数；组内项为两条guide分别随机拆组的差异再平均。\n'
             f'论文异靶配对{100*batch["different_cross_batch_fraction"]:.1f}%跨批（同靶{100*batch["query_cross_batch_fraction"]:.1f}%）；绝对值不直接代表实验质量。',
             ha='center',va='bottom',fontsize=9,color='#555555')
    fig.savefig(out/'01_variation_boxplots.png',dpi=180);fig.savefig(out/'01_variation_boxplots.svg');plt.close(fig)
    pd.DataFrame(rows).to_csv(out/'variation_summary.csv',index=False)
    case_data={};bound=0;case_summaries=[]
    for name in DATASETS:
        table=pd.read_csv(root/name/'RNA_case_points.csv.gz')
        if not table.cells_per_group.eq(n).all() or not np.isfinite(table[['x','y']]).all().all():raise ValueError('Invalid RNA case values')
        case_data[name]=table;bound=max(bound,float(np.abs(table[['x','y']].to_numpy()).max()))
    bound=max(.1,bound*1.07)
    fig,axes=plt.subplots(3,2,figsize=(13.1,13.6),sharex=True,sharey=True,layout='constrained')
    for row,name in enumerate(DATASETS):
        metadata=json.loads((root/name/'case_selection.json').read_text());records=pd.read_csv(root/name/'eligible_guides.csv').set_index('guide')
        for col,kind in enumerate(['same_target','different_target']):
            table=case_data[name].loc[case_data[name].comparison.eq(kind)];a,b=table.guide_a.iloc[0],table.guide_b.iloc[0]
            if not table.gene_id.is_unique:raise ValueError('Duplicated RNA in case scatter')
            ax=axes[row,col];ax.scatter(table.x,table.y,s=8,alpha=.3,color=COLORS[name],linewidths=0)
            ax.plot([-bound,bound],[-bound,bound],ls='--',color='#333333',lw=.8);ax.axhline(0,lw=.4,color='#aaaaaa');ax.axvline(0,lw=.4,color='#aaaaaa')
            mae=float(np.mean(np.abs(table.x-table.y)))
            ax.set(xlim=(-bound,bound),ylim=(-bound,bound),aspect='equal',
                   xlabel='第一组细胞的RNA变化（log2FC）',ylabel='第二组细胞的RNA变化（log2FC）',
                   title=LABELS[name]+'｜'+('同靶点' if kind=='same_target' else '不同靶点'))
            ax.text(.02,.97,f'第一组：{short_guide(a)}\n第二组：{short_guide(b)}\n每组{n}细胞；平均绝对差 {mae:.3f}',transform=ax.transAxes,
                    va='top',fontsize=8.5,bbox=dict(facecolor='white',edgecolor='none',alpha=.8))
            # Labels describe large differences only; they are not a DE test.
            top=table.assign(delta=np.abs(table.x-table.y)).nlargest(3,'delta')
            for k,r in enumerate(top.itertuples()):
                label=r.gene_name if pd.notna(r.gene_name) and str(r.gene_name)!='nan' else r.gene_id
                ax.annotate(str(label),(r.x,r.y),xytext=(5,5+8*k),textcoords='offset points',fontsize=7,alpha=.85)
            case_summaries.append(dict(dataset=name,comparison=kind,guide_a=a,guide_b=b,target_a=str(records.loc[a,'target']),target_b=str(records.loc[b,'target']),
                original_cells_a=int(records.loc[a,'cells']),original_cells_b=int(records.loc[b,'cells']),sampled_cells=n,repeat=1,RNA=len(table),MAE=mae))
    fig.suptitle('具体RNA变化：每个点是一个RNA；虚线为两组变化相同的 y=x',fontsize=14)
    fig.savefig(out/'02_RNA_pair_scatter.png',dpi=175);fig.savefig(out/'02_RNA_pair_scatter.svg');plt.close(fig)
    pd.DataFrame(case_summaries).to_csv(out/'case_summary.csv',index=False)
    n5=validate_points(root,n//2);sensitivity=[]
    for name in DATASETS:
        for kind in KINDS:
            a=points[name].loc[points[name].comparison.eq(kind)].set_index('pair_id');b=n5[name].loc[n5[name].comparison.eq(kind)].set_index('pair_id')
            if set(a.index)!=set(b.index):raise ValueError('Sample-size sensitivity changed guide-pair cohort')
            sensitivity.append(dict(dataset=name,comparison=kind,same_pairs=len(a),n_small=n//2,median_small=float(b.mean_absolute_difference.median()),n_primary=n,median_primary=float(a.mean_absolute_difference.median())))
    pd.DataFrame(sensitivity).to_csv(out/'same_cohort_sample_size_check.csv',index=False)
    provenance=dict(source_summary=summary,axes_shared=True,main_n=n,same_four_arm_pair_cohort=True,
                    case_repeat=1,figure1_each_point='median of8 draws; within is average of each query guide split distance',
                    direct_target_RNAs_not_used_for_distance=True,example_selection='coverage only; no common assessable target across all3 datasets',
                    plot_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),paper_batch_diagnostic=batch)
    (out/'plot_provenance.json').write_text(json.dumps(provenance,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    return rows


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--paper-metadata',type=Path,required=True);a=p.parse_args()
    print(json.dumps(render(a.source,a.output,a.paper_metadata),ensure_ascii=False))
