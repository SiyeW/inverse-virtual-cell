"""An accepted overview and linked sampling interval details in one figure."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import font_manager

LAB_CASES=['EIF2AK3','GLI3','HES1','SMAD3','SOX2','UBA1']
PAPER_CASES=['RNF7','POU5F1']  # Cases fixed before the original reference analysis.
DATASETS=['chip01','company','kolf']
LABELS=['我们的首芯片','公司合并矩阵','论文 KOLF Strong']


def validate_points(old,new):
    if not old.guide.is_unique or not new.guide.is_unique or set(old.guide)!=set(new.guide):raise ValueError('Guide cohort changed')
    a=old.set_index('guide').sort_index();b=new.set_index('guide').sort_index()
    for col in ['cells','mean_per_10000','NTC_mean_per_10000','RNA_remaining_ratio','NTC_detected_fraction']:
        if not np.allclose(a[col],b[col],equal_nan=True,rtol=1e-12,atol=1e-12):raise ValueError('Original points changed')
    zero=b.guide_ci_status.eq('all_zero_no_RNA_mean_interval')
    if b.loc[zero,['guide_ci_low','guide_ci_high']].notna().any().any():raise ValueError('False precise interval for all-zero observations')
    eligible=b.cells.ge(10)&b.NTC_mean_per_10000.gt(0)
    if b.loc[eligible,'NTC_equal_n_low'].isna().any():raise ValueError('Missing NTC ranges')


def detail_rows(table,dataset):
    targets=LAB_CASES if dataset!='kolf' else PAPER_CASES
    q=table.loc[table.target.isin(targets)&table.cells.ge(10)&table.NTC_mean_per_10000.gt(0)].copy()
    q['order']=q.target.map({v:i for i,v in enumerate(targets)})
    return q.sort_values(['order','guide']).reset_index(drop=True)


def draw_details(axes,tables,title_prefix=True):
    selections={};cap=400.
    for d,label,ax in zip(DATASETS,LABELS,axes):
        rows=detail_rows(tables[d],d);selections[d]=rows.guide.tolist()
        n=len(rows)
        if n==0:
            ax.text(.5,.5,'预先指定的例子不满足覆盖条件',ha='center',transform=ax.transAxes);continue
        for k,row in enumerate(rows.itertuples()):
            lo=100*row.NTC_equal_n_low;hi=100*row.NTC_equal_n_high
            ax.plot([min(lo,cap),min(hi,cap)],[k,k],c='#b9b9b9',lw=6,alpha=.65,solid_capstyle='butt',zorder=1)
            if np.isfinite(row.guide_ci_low) and np.isfinite(row.guide_ci_high):
                ax.plot([min(100*row.guide_ci_low,cap),min(100*row.guide_ci_high,cap)],[k,k],c='#216a9b',lw=1.7,zorder=3)
                ax.plot([min(100*row.guide_ci_low,cap),min(100*row.guide_ci_high,cap)],[k,k],ls='none',marker='|',c='#216a9b',ms=6,zorder=3)
            value=min(100*row.RNA_remaining_ratio,cap)
            fill='white' if row.guide_ci_status=='all_zero_no_RNA_mean_interval' else '#182731'
            ax.scatter([value],[k],s=29,facecolor=fill,edgecolor='#182731',lw=.8,zorder=4)
            if row.source_control_type=='positive_control':ax.scatter([value],[k],s=66,marker='D',facecolors='none',edgecolors='#d17e24',lw=.9,zorder=5)
            if hi>cap or (np.isfinite(row.guide_ci_high) and 100*row.guide_ci_high>cap) or 100*row.RNA_remaining_ratio>cap:
                ax.plot(cap,k,marker='>',c='#555555',ms=5,zorder=6)
        ax.axvline(100,c='#333333',lw=.9,ls='--');ax.set_xlim(-5,cap+8)
        ax.set_ylim(n-.4,-.7);ax.set_yticks(range(n));ax.set_yticklabels([f'{r.guide}  n={r.cells}' for r in rows.itertuples()],fontsize=9)
        ax.tick_params(axis='y',length=0);ax.set_xlabel('靶RNA剩余量（%，100%=当前NTC均值）',fontsize=10)
        ax.set_title(label+('：项目预标注的6个靶点' if d!='kolf' else '：预先指定的RNF7、POU5F1'),fontsize=10,pad=12)
        ax.grid(axis='x',alpha=.15);ax.spines[['top','right','left']].set_visible(False)
        # Every displayed point remains traceable to the full table, no outcome selection.
    return selections


def render(oldroot,newroot,overview,out,font=None):
    if font:font_manager.fontManager.addfont(font);plt.rcParams['font.family']=font_manager.FontProperties(fname=font).get_name()
    else:plt.rcParams['font.family']='Microsoft YaHei'
    plt.rcParams['axes.unicode_minus']=False
    out.mkdir(parents=True,exist_ok=True);tables={};summaries=[]
    for d in DATASETS:
        old=pd.read_csv(oldroot/d/'guide_target_RNA.csv');new=pd.read_csv(newroot/d/'guide_target_RNA_precision.csv')
        validate_points(old,new);tables[d]=new
        p=new.loc[new.cells.ge(10)&new.NTC_mean_per_10000.gt(0)]
        measured=p.loc[p.NTC_detected_fraction.ge(.1)]
        summaries.append(dict(dataset=d,plotted_points=len(p),NTC_detected_ge10=len(measured),
                              estimated_intervals=int(measured.guide_ci_status.eq('estimated').sum()),
                              all_zero_no_RNA_CI=int(measured.guide_ci_status.eq('all_zero_no_RNA_mean_interval').sum()),
                              median_RNA_remaining_ratio=float(measured.RNA_remaining_ratio.median())))
    # The overview is an unmodified accepted image, not a regenerated analysis.
    fig=plt.figure(figsize=(17,15.5));top=fig.add_axes([.015,.55,.97,.43])
    top.imshow(plt.imread(overview));top.axis('off')
    fig.text(.5,.535,'下方：同一图中的具体sgRNA区间细节（完整逐guide结果见表）',ha='center',fontsize=14)
    axes=[fig.add_axes([.10,.12,.215,.355]),fig.add_axes([.43,.12,.215,.355]),fig.add_axes([.76,.12,.215,.355])]
    selections=draw_details(axes,tables)
    handles=[Line2D([0],[0],c='#b9b9b9',lw=6,label='灰粗线：相同细胞数的NTC组95%波动范围'),
             Line2D([0],[0],c='#216a9b',lw=1.8,label='蓝细线：guide均值的95%重抽样区间'),
             Line2D([0],[0],marker='o',c='none',markerfacecolor='#182731',markeredgecolor='#182731',label='点：原来的靶RNA剩余量')]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.07),ncol=1,frameon=False,fontsize=10)
    fig.text(.5,.018,'每条guide使用全部合格细胞；灰线按该guide的实际细胞数从NTC取组，参照排除这组NTC本身。\n蓝线固定原NTC均值，以SciPy BCa方法重抽guide细胞2000次；它是当前细胞集合的抽样精度，不是生物重复或编辑率。\n全未检出的guide不画0±0区间；下方例子在看结果前已指定，不能当全库代表。>400%用右箭头，完整值见表。',ha='center',fontsize=9)
    for ext in ['png','svg']:fig.savefig(out/('靶RNA总览与抽样区间.'+ext),dpi=180)
    plt.close(fig)
    pd.DataFrame(summaries).to_csv(out.parent/'figure_interval_summary.csv',index=False)
    pd.concat([detail_rows(tables[d],d) for d in DATASETS],ignore_index=True).to_csv(out.parent/'展示例子的完整区间.csv',index=False)
    provenance=dict(original_overview_sha256=hashlib.sha256(overview.read_bytes()).hexdigest(),original_points_unchanged=True,
                    detail_selection_lab=LAB_CASES,detail_selection_paper=PAPER_CASES,selected_guides=selections,
                    detail_selection_based_on_observed_response=False,NTC_pool_and_RNA_normalization_unchanged=True)
    (out.parent/'combination_provenance.json').write_text(json.dumps(provenance,indent=2),encoding='utf8')
    return summaries


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--original',type=Path,required=True);p.add_argument('--precision',type=Path,required=True);p.add_argument('--overview',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--font');a=p.parse_args()
    print(json.dumps(render(a.original,a.precision,a.overview,a.output,a.font),indent=2))
