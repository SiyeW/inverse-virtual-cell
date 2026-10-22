"""Static, fixed-scope figures for a control-aware single-chip QC report."""
from pathlib import Path
from textwrap import fill
import numpy as np
import pandas as pd
from perturbseq_core_effects import TYPES, pair_effects

NAMES=dict(zip(TYPES,['NTC','Intergenic','Annotated positive','Other targeting']))
COLORS=dict(zip(TYPES,['#3c7c9f','#898292','#6e9d70','#cc8940']))


def save(fig,path):
    import matplotlib.pyplot as plt
    fig.savefig(path,dpi=160,facecolor='white',bbox_inches='tight')
    plt.close(fig)


def technical_figures(all_obs,overlay,flow,output):
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    x=np.arange(len(flow));axes[0].bar(x,flow.rna_qc_pass,label='RNA QC pass',color='#3c7c9f')
    axes[0].bar(x,flow.rna_qc_not_pass,bottom=flow.rna_qc_pass,label='RNA QC not pass',color='#c9ced0')
    axes[0].set(xticks=x,xticklabels=[s.replace('_','\n') for s in flow.guide_status],ylabel='Source cell records',title='Guide assignment and RNA QC')
    axes[0].legend(frameon=False)
    total=len(all_obs);rna=int(all_obs.qc_pass_stage1.astype(bool).sum());single=int(all_obs.guide_status.astype(str).eq('single_guide').sum());both=int(all_obs.main_analysis_stage1.astype(bool).sum())
    values=[total,rna,single,both]
    axes[1].bar(np.arange(4),values,color=['#9ba6ac','#7e9eaf','#70998d','#3c7c9f'])
    axes[1].set(xticks=np.arange(4),xticklabels=['All\nrecords','RNA QC\npass','Single guide\nbefore RNA QC','Both\nrequirements'],ylabel='Cell records',title='Intersecting requirements, not additive stages')
    for i,v in enumerate(values):axes[1].text(i,v,f'{v:,}',ha='center',va='bottom')
    save(fig,output/'figures/01_technical_qc.png')

    fig,axes=plt.subplots(1,3,figsize=(14,4.3),layout='constrained')
    for kind in TYPES:
        part=overlay.loc[overlay.source_control_type.eq(kind)];v=part.single_guide_after_stage1_qc.to_numpy()
        if not len(v):continue
        order=np.sort(v);axes[0].step(order,np.arange(1,len(v)+1)/len(v),where='post',color=COLORS[kind],label=f'{NAMES[kind]} ({len(v)})')
    axes[0].set(xscale='symlog',xlabel='RNA-QC single-guide cells per guide (zeros kept)',ylabel='Fraction of library guides',title='Coverage by source annotation')
    axes[0].legend(frameon=False,fontsize=8)
    counts=[overlay.loc[overlay.source_control_type.eq(k),'single_guide_after_stage1_qc'].sum() for k in TYPES]
    zeros=[100*overlay.loc[overlay.source_control_type.eq(k),'single_guide_after_stage1_qc'].eq(0).mean() for k in TYPES]
    for ax,values,title,ylabel in [(axes[1],counts,'Captured single-guide cells','RNA-QC cell records'),(axes[2],zeros,'Zero-coverage designs','Within-type guides (%)')]:
        ax.bar(np.arange(4),values,color=[COLORS[k] for k in TYPES]);ax.set(xticks=np.arange(4),xticklabels=['NTC','Intergenic','Annotated\npositive','Other\ntargeting'],title=title,ylabel=ylabel)
        for i,v in enumerate(values):
            label=(f'{v:,.0f}' if ax is axes[1] else f'{v:.1f}%') if np.isfinite(v) else 'N/A'
            ax.text(i,v if np.isfinite(v) else 0,label,ha='center',va='bottom',fontsize=9)
    save(fig,output/'figures/02_four_type_coverage.png')

    cell=all_obs.loc[all_obs.guide_status.astype(str).eq('single_guide')].copy()
    cell['source_control_type']=cell.guide_id.astype(str).map(overlay.set_index('guide').source_control_type)
    post=cell.loc[cell.qc_pass_stage1.astype(bool)]
    fig,axes=plt.subplots(2,3,figsize=(14,7.5),layout='constrained')
    for ax,metric,label in zip(axes.flat,['total_counts','n_genes_by_counts','pct_counts_mt','raw_guide_umis','guide_umis'],['RNA UMIs','Detected RNA features','Mitochondrial RNA (%)','Captured guide UMIs','Assigned guide UMIs']):
        values=[post.loc[post.source_control_type.eq(k),metric].to_numpy(float) for k in TYPES]
        boxes=ax.boxplot(values,tick_labels=['NTC','Intergenic','Annotated\npositive','Other\ntargeting'],patch_artist=True,showfliers=False)
        for patch,k in zip(boxes['boxes'],TYPES):patch.set_facecolor(COLORS[k])
        ax.set(title=label,ylabel=label)
        if metric!='pct_counts_mt':ax.set_yscale('symlog')
    values=[100*cell.loc[cell.source_control_type.eq(k),'qc_pass_stage1'].astype(bool).mean() for k in TYPES]
    axes[1,2].bar(np.arange(4),values,color=[COLORS[k] for k in TYPES]);axes[1,2].set(xticks=np.arange(4),xticklabels=['NTC','Intergenic','Annotated\npositive','Other\ntargeting'],ylim=(0,100),ylabel='RNA-QC retained (%)',title='Within single-guide calls')
    for i,v in enumerate(values):axes[1,2].text(i,v,f'{v:.1f}%',ha='center',va='bottom')
    fig.suptitle('Post-QC metrics; boxplot fliers omitted for display only; no new filtering')
    save(fig,output/'figures/03_four_type_cell_quality.png')


def effect_background_figures(records,effects,matched,metrics,examples,output):
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(13,4.5),layout='constrained')
    finite=matched[np.isfinite(matched)];bound=max(.05,float(np.max(np.abs(finite))))
    edges=np.linspace(-bound,bound,181)
    summaries=[]
    medians=metrics.groupby('guide',sort=True).rms_log2fc.median()
    for kind in TYPES:
        use=records.source_control_type.eq(kind).to_numpy();values=matched[use]
        if not len(values):continue
        weights=np.where(np.isfinite(values),1/np.maximum(np.isfinite(values).sum(1),1)[:,None],0)
        hist,_=np.histogram(values[np.isfinite(values)],bins=edges,weights=weights[np.isfinite(values)])
        axes[0].stairs(hist/(hist.sum()*np.diff(edges)),edges,color=COLORS[kind],label=f'{NAMES[kind]}: {int(use.sum())} guides')
        rms=np.sort(medians.loc[records.loc[use,'guide']].to_numpy())
        if len(rms)==1:axes[1].axvline(rms[0],color=COLORS[kind],ls=':',label=NAMES[kind]+' (one guide)')
        else:axes[1].step(rms,np.arange(1,len(rms)+1)/len(rms),where='post',color=COLORS[kind],label=NAMES[kind])
        flattened=values[np.isfinite(values)]
        summaries.append(dict(source_control_type=kind,guides=int(use.sum()),matched_cells=int(metrics.cells.iloc[0]),
             median_signed_log2fc=float(np.median(flattened)),mad_signed_log2fc=float(np.median(np.abs(flattened-np.median(flattened)))),
             median_guide_rms_log2fc=float(np.median(rms)),maximum_guide_rms_log2fc=float(np.max(rms))))
    axes[0].set(xlabel='RNA log2FC to nonoverlapping pure NTCs (+1 pseudocount)',ylabel='Guide-balanced RNA-feature density',title='First fixed equal-size subsample; full tails kept',yscale='symlog')
    axes[0].axvline(0,color='black',ls='--',lw=.8);axes[0].legend(frameon=False,fontsize=8)
    axes[1].set(xlabel='Per-guide RMS RNA log2FC (median of fixed subsamples)',ylabel='Fraction of eligible guides',title='Descriptive magnitudes, not p-values',ylim=(0,1.02))
    axes[1].legend(frameon=False,fontsize=8)
    fig.suptitle('Equal-size RNA responses; genes, guides and subsamples are not independent biological replicates')
    save(fig,output/'figures/05_matched_rna_effects.png')
    pd.DataFrame(summaries).to_csv(output/'tables/effect_background_summary.csv',index=False)
    rows=max(1,(len(examples)+1)//2)
    fig,axes=plt.subplots(rows,2,figsize=(12,3.7*rows),layout='constrained',squeeze=False)
    for ax,example in zip(axes.flat,examples):
        x=np.log1p(np.asarray(example['reference_mean']));y=np.log1p(np.asarray(example['mean']))
        ax.scatter(x,y,s=3,alpha=.35,color=example['color']);bound=max(float(x.max()),float(y.max()))
        ax.plot([0,bound],[0,bound],color='black',ls='--',lw=.8)
        ax.set(title=fill(example['label'],55),xlabel='Reference RNA mean (log1p)',ylabel='Sampled RNA mean (log1p)')
    for ax in axes.flat[len(examples):]:ax.set_visible(False)
    fig.suptitle('Same fixed RNA features and cell count; RNA expression, not Day0/Day10 guide abundance')
    save(fig,output/'figures/06_rna_mean_scatter.png')
    return pd.DataFrame(summaries)


def target_atlas(overlay,records,effects,support,pairs,output,max_cards=None,min_cells=10):
    import matplotlib.pyplot as plt
    lookup={g:i for i,g in enumerate(records.guide)}
    cards=[]
    lines=['# 逐靶点多sgRNA支持图册','',
      '本图册按靶点名称排序，展示所有具备至少两条可描述guide的靶点，并单列候选正控。其余靶点仍在完整表中，不把未评价当零效应。',
      '',f'同一芯片的guide不是独立生物学重复。RNA log2FC采用固定深度归一化后的均值比（伪计数1）；各guide至少{min_cells}个细胞才比较响应谱，不使用AD/MW显著性挑选靶点或RNA。',
      '','散点使用细胞最多的两条合格guide，以及全部既定RNA特征；不按响应强弱选guide。唯一直接靶RNA从下游响应比较中屏蔽，完整名称保留在表中。','',
      '[全部靶点及未评价原因](tables/target_support.csv)｜[全部guide覆盖与角色](tables/control_annotations.csv)｜[逐guide响应元数据](tables/guide_effect_metadata.csv)','']
    for number,row in enumerate(support.itertuples()):
        if row.assessment!='descriptive_assessable' and row.source_control_type!='positive_control':continue
        if max_cards is not None and len(cards)>=max_cards:break
        library=overlay.loc[overlay.target_gene_name.eq(row.target)].sort_values('guide')
        eligible=records.loc[records.target.eq(row.target)].sort_values(['cells','guide'],ascending=[False,True])
        fig,axes=plt.subplots(1,3,figsize=(13,4),layout='constrained')
        axes[0].bar(np.arange(len(library)),library.single_guide_after_stage1_qc,color=COLORS[row.source_control_type])
        axes[0].axhline(min_cells,color='black',ls='--',lw=.8)
        axes[0].set(xticks=np.arange(len(library)),xticklabels=[f'#{i+1}' for i in range(len(library))],ylabel='RNA-QC single-guide cells',title='All library guides, including zero coverage')
        if len(eligible)>1:
            values=effects[[lookup[g] for g in eligible.guide]]
            correlation=np.full((len(values),len(values)),np.nan)
            for i in range(len(values)):
                for j in range(len(values)):correlation[i,j]=pair_effects(values[i],values[j])['pearson']
            image=axes[1].imshow(correlation,vmin=-1,vmax=1,cmap='coolwarm')
            labels=[f'#{list(library.guide).index(g)+1}' for g in eligible.guide]
            axes[1].set(xticks=np.arange(len(labels)),xticklabels=labels,yticks=np.arange(len(labels)),yticklabels=labels,title='RNA log2FC Pearson r')
            for i in range(len(values)):
                for j in range(len(values)):axes[1].text(j,i,f'{correlation[i,j]:.2f}',ha='center',va='center',fontsize=8)
            fig.colorbar(image,ax=axes[1],fraction=.046)
            a,b=values[:2];keep=np.isfinite(a)&np.isfinite(b)
            axes[2].scatter(a[keep],b[keep],s=3,alpha=.3,color='#3c7c9f');bound=max(.05,float(np.max(np.abs(np.r_[a[keep],b[keep]]))))
            axes[2].plot([-bound,bound],[-bound,bound],color='black',ls='--',lw=.8);axes[2].axhline(0,lw=.5,color='grey');axes[2].axvline(0,lw=.5,color='grey')
            axes[2].set(xlabel=f'{labels[0]} RNA log2FC (n={eligible.cells.iloc[0]})',ylabel=f'{labels[1]} RNA log2FC (n={eligible.cells.iloc[1]})',title='Two highest-coverage guides; all fixed RNA features')
        else:
            for ax in axes[1:]:ax.axis('off');ax.text(.5,.5,'Fewer than two eligible guides\nNo pairwise agreement assessment',ha='center',va='center',transform=ax.transAxes)
        fig.suptitle(fill(row.target,100)+'\nDescriptive only; no hit or biological-replication claim',fontsize=11)
        rel=f'figures/targets/target_{number:03d}.png';save(fig,output/rel)
        cards.append(dict(target=row.target,figure=rel,source_control_type=row.source_control_type,assessment=row.assessment))
        lines += [f'## {number+1}. {row.target}','',f'![{row.target}的多guide支持]({rel})','',
          f'RNA质控后共{row.cells}个细胞，{row.eligible_guides}/{row.library_guides}条guide达到描述入口。最多细胞guide占该靶点细胞的{100*row.largest_guide_fraction:.1f}%；去掉该guide后剩{row.cells_without_largest}个细胞。' if row.cells else '当前没有通过RNA质控的单guide细胞，无法评价响应。','',
          '中间矩阵对角线是guide与自身相关，不是重复验证；右图只使用细胞最多的两条合格guide。', '',
          '| 图中编号 | sgRNA | 细胞数 |','|---|---|---:|']
        lines += [f'| #{i+1} | {r.guide} | {r.single_guide_after_stage1_qc} |' for i,r in enumerate(library.itertuples())]
        lines += ['', '覆盖不足或相关性低不等于靶点无功能；相同方向的小幅变化也不等于生物学重复验证。','']
    (output/'逐靶点支持图册.md').write_text('\n'.join(lines)+'\n',encoding='utf8')
    pd.DataFrame(cards,columns=['target','figure','source_control_type','assessment']).to_csv(output/'tables/target_figure_index.csv',index=False)
    return cards
