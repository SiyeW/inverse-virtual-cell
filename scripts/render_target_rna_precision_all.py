"""Show every originally measurable guide, not analyst-selected examples."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from render_target_rna_precision import validate_points

DATASETS=['chip01','company','kolf']
LABELS=['我们的首芯片','公司合并矩阵','论文 KOLF Strong']


def whole_rows(table):
    q=table.loc[table.cells.ge(10)&table.NTC_mean_per_10000.gt(0)&table.NTC_detected_fraction.ge(.1)].copy()
    q=q.sort_values(['target_gene_id','guide']).reset_index(drop=True)
    q.insert(0,'plot_row',np.arange(1,len(q)+1))
    return q


def render(original,precision,overview,out,font=None):
    if font:font_manager.fontManager.addfont(font);plt.rcParams['font.family']=font_manager.FontProperties(fname=font).get_name()
    else:plt.rcParams['font.family']='Microsoft YaHei'
    plt.rcParams['axes.unicode_minus']=False;out.mkdir(parents=True,exist_ok=True)
    fig=plt.figure(figsize=(17,16));ax=fig.add_axes([.015,.56,.97,.42]);ax.imshow(plt.imread(overview));ax.axis('off')
    fig.text(.5,.542,'下排：全部符合原蓝点条件的sgRNA区间',ha='center',fontsize=15)
    fig.text(.5,.524,'每一横行是一条sgRNA，按靶RNA编号和guide名称排序；三列的行号各自独立。',ha='center',fontsize=10)
    axes=[fig.add_axes([.065,.13,.26,.35]),fig.add_axes([.39,.13,.26,.35]),fig.add_axes([.715,.13,.26,.35])]
    summaries=[];orders=[];provenance={};cap=400.
    for d,label,ax in zip(DATASETS,LABELS,axes):
        old=pd.read_csv(original/d/'guide_target_RNA.csv');source=precision/d/'guide_target_RNA_precision.csv';table=pd.read_csv(source)
        validate_points(old,table);q=whole_rows(table);n=len(q);y=q.plot_row.to_numpy()
        assert n and not q.guide.duplicated().any()
        ci=q.guide_ci_status.eq('estimated').to_numpy();zero=q.guide_ci_status.eq('all_zero_no_RNA_mean_interval').to_numpy()
        if not np.isfinite(q[['NTC_equal_n_low','NTC_equal_n_high']].to_numpy()).all():raise ValueError('Unbounded null range in measurable plot cohort')
        ax.hlines(y,np.minimum(100*q.NTC_equal_n_low,cap),np.minimum(100*q.NTC_equal_n_high,cap),color='#b4b4b4',lw=.85,alpha=.75,zorder=1)
        ax.hlines(y[ci],np.minimum(100*q.loc[ci,'guide_ci_low'],cap),np.minimum(100*q.loc[ci,'guide_ci_high'],cap),color='#1a6fa5',lw=.75,alpha=.9,zorder=2)
        values=np.minimum(100*q.RNA_remaining_ratio.to_numpy(),cap)
        ax.scatter(values[~zero],y[~zero],s=4 if n>500 else 8,color='#14232e',linewidths=0,zorder=3)
        ax.scatter(values[zero],y[zero],s=5 if n>500 else 12,facecolors='white',edgecolors='#14232e',linewidths=.45,zorder=4)
        beyond=q.NTC_equal_n_high.gt(cap/100)|q.guide_ci_high.gt(cap/100)|q.RNA_remaining_ratio.gt(cap/100)
        ax.scatter(np.full(int(beyond.sum()),cap),y[beyond],s=7,c='#333333',marker='>',zorder=5)
        ax.axvline(100,ls='--',lw=.9,c='#333333');ax.set_xlim(-5,cap+8);ax.set_ylim(n+.5,.5)
        ticks=np.unique(np.rint(np.linspace(1,n,6)).astype(int));ax.set_yticks(ticks);ax.set_yticklabels(ticks)
        ax.set_ylabel('sgRNA排列行号',fontsize=10);ax.set_xlabel('靶RNA剩余量（%，100%=当前NTC均值）',fontsize=10)
        ax.set_title(label+f'：全部{n}条',fontsize=12,pad=12);ax.grid(axis='x',alpha=.12);ax.spines[['top','right']].set_visible(False)
        inside=q.RNA_remaining_ratio.ge(q.NTC_equal_n_low)&q.RNA_remaining_ratio.le(q.NTC_equal_n_high)
        covered=ci & q.guide_ci_low.le(1).to_numpy() & q.guide_ci_high.ge(1).to_numpy()
        summary=dict(dataset=d,measured_guides=n,estimated_CI=int(ci.sum()),all_zero_no_CI=int(zero.sum()),
                     point_inside_NTC95=int(inside.sum()),point_below_NTC95=int(q.RNA_remaining_ratio.lt(q.NTC_equal_n_low).sum()),
                     point_above_NTC95=int(q.RNA_remaining_ratio.gt(q.NTC_equal_n_high).sum()),CI_contains100=int(covered.sum()))
        summaries.append(summary)
        ax.text(.5,-.125,f'原点在NTC范围内：{inside.sum()}/{n}\n蓝线包含100%：{covered.sum()}/{ci.sum()}；全未检出无蓝线：{zero.sum()}条',transform=ax.transAxes,ha='center',va='top',fontsize=9)
        q.insert(0,'dataset_plot',d);orders.append(q)
        provenance[d]=dict(source_table_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),all_original_blue_points_included=True,
                           cohort='cells>=10; NTC_mean>0; NTC_detected_fraction>=0.10, identical to original blue-point rule',
                           sort='target_gene_id then guide, not RNA effect',plotted_guides=q.guide.tolist())
    handles=[Line2D([0],[0],color='#b4b4b4',lw=4,label='灰线：同n NTC组95%波动范围'),
             Line2D([0],[0],color='#1a6fa5',lw=1.4,label='蓝线：guide均值95%重抽样区间'),
             Line2D([0],[0],color='none',marker='o',markerfacecolor='#14232e',markeredgecolor='#14232e',label='点：原RNA剩余量')]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.049),ncol=3,frameon=False,fontsize=10)
    fig.text(.5,.014,'下排覆盖全部符合原蓝点条件的sgRNA：至少10个细胞、NTC检出率至少10%。灰点仍在上排，完整3216条结果保留在表。\n全未检出不画0±0；超400%用右箭头。区间是当前细胞集合的抽样精度，不是编辑率或独立实验重复；不据范围位置判有效/无效。',ha='center',fontsize=9)
    for ext in ['png','svg']:fig.savefig(out/('靶RNA总览与全体区间.'+ext),dpi=180)
    plt.close(fig)
    pd.DataFrame(summaries).to_csv(out.parent/'全体区间汇总.csv',index=False)
    pd.concat(orders,ignore_index=True).to_csv(out.parent/'下排全部sgRNA行号与区间.csv',index=False)
    proof=dict(original_overview_sha256=hashlib.sha256(overview.read_bytes()).hexdigest(),datasets=provenance,
               analyst_selected_cases=False,all_original_point_values_and_intervals_unchanged=True,scientific_recalculation=False)
    (out.parent/'whole_cohort_provenance.json').write_text(json.dumps(proof,ensure_ascii=False,indent=2),encoding='utf8')
    return summaries


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--original',type=Path,required=True);p.add_argument('--precision',type=Path,required=True);p.add_argument('--overview',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--font');a=p.parse_args()
    print(json.dumps(render(a.original,a.precision,a.overview,a.output,a.font),indent=2))
