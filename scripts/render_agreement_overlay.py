"""Render fixed-table cross-dataset comparisons; never select favorable results."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

DATASETS=['chip01','company','kolf']
LABELS={'chip01':'我们的首芯片','company':'公司合并矩阵','kolf':'论文 KOLF Strong'}
COLORS={'chip01':'#397ca9','company':'#c68a45','kolf':'#439477'}
METHODS={'old':'原方法：先求平均，再取log','new':'统一改法：先逐细胞取log，再求平均'}
COMPARISONS={'same_target':'同一基因的两条sgRNA','different_target':'不同基因的两条sgRNA','NTC':'两条NTC sgRNA'}


def load_tables(root,suffix=''):
    tables={}
    for name in DATASETS:
        for method in METHODS:
            p=root/f'{name}_{method}{suffix}_points.csv'
            table=pd.read_csv(p)
            if not table.target.is_unique or not np.isfinite(table[['largest_guide_fraction','median_pair_r']]).all().all():raise ValueError('Invalid target points')
            if ((table.largest_guide_fraction<=0)|(table.largest_guide_fraction>=1)).any():raise ValueError('Unassessable singleton target in scatter')
            tables[name,method]=table
        a=tables[name,'old'].set_index('target');b=tables[name,'new'].set_index('target')
        if set(a.index)!=set(b.index) or not np.array_equal(a.largest_guide_fraction,b.reindex(a.index).largest_guide_fraction):raise ValueError('Old/new points differ in target or cell coverage')
    return tables


def scatter(tables,path,title):
    fig,axes=plt.subplots(1,2,figsize=(12.7,5.1),sharex=True,sharey=True,layout='constrained')
    for ax,(method,heading) in zip(axes,METHODS.items()):
        for name in DATASETS[::-1]:
            t=tables[name,method];median=t.median_pair_r.median()
            ax.scatter(t.largest_guide_fraction,t.median_pair_r,s=20 if name!='kolf' else 14,
                alpha=.7 if name!='kolf' else .45,color=COLORS[name],linewidths=.3,edgecolors='white',
                label=f'{LABELS[name]}：{len(t)}靶点，中位数 {median:.3f}',zorder=3 if name!='kolf' else 2)
        ax.axhline(0,color='#777777',lw=.8,zorder=1)
        ax.set(title=heading,xlim=(0,1),ylim=(-1,1),xlabel='细胞最多的sgRNA占该靶点细胞数的比例')
        ax.legend(loc='upper right',frameon=False,fontsize=9)
        ax.spines[['top','right']].set_visible(False)
    axes[0].set_ylabel('同一靶点不同sgRNA的RNA变化相关\n（每个靶点内取配对相关系数的中位数）')
    fig.suptitle(title,fontsize=13)
    fig.savefig(path,dpi=180);fig.savefig(path.with_suffix('.svg'));plt.close(fig)


def backgrounds(root,path,suffix=''):
    fig,axes=plt.subplots(2,3,figsize=(13.7,7.5),sharex=True,sharey=True,layout='constrained')
    palette={'same_target':'#c68a45','different_target':'#7a8490','NTC':'#8a6aa0'}
    rows=[]
    for col,name in enumerate(DATASETS):
        for row,method in enumerate(METHODS):
            table=pd.read_csv(root/f'{name}_{method}{suffix}_pairs.csv');ax=axes[row,col]
            for kind in COMPARISONS:
                x=np.sort(table.loc[table.comparison.eq(kind),'pearson'].to_numpy())
                if not len(x) or not np.isfinite(x).all():raise ValueError('Missing fixed background')
                ax.step(np.r_[x[0],x],np.r_[0,np.arange(1,len(x)+1)/len(x)],where='post',color=palette[kind],lw=1.5,
                        label=f'{COMPARISONS[kind]}（中位数{float(np.median(x)):.3f}）')
                rows.append(dict(dataset=name,method=method,comparison=kind,pairs=len(x),median=float(np.median(x))))
            ax.axvline(0,color='#777777',lw=.6)
            ax.set(xlim=(-1,1),ylim=(0,1.02),title=LABELS[name]+'｜'+('原方法' if method=='old' else '先log后平均'))
            ax.legend(frameon=False,fontsize=8,loc='lower right');ax.spines[['top','right']].set_visible(False)
            if row==1:ax.set_xlabel('两条sgRNA的RNA变化相关系数')
            if col==0:ax.set_ylabel('不超过横轴相关系数的配对比例')
    fig.suptitle('原方法与统一改法：同基因、不同基因和NTC配对',fontsize=13)
    fig.savefig(path,dpi=180);fig.savefig(path.with_suffix('.svg'));plt.close(fig)
    return rows


def render(root,output):
    root=Path(root);output=Path(output);output.mkdir(exist_ok=False,parents=True)
    font=Path('C:/Windows/Fonts/msyh.ttc')
    if not font.exists():raise ValueError('Readable Chinese font unavailable')
    font_manager.fontManager.addfont(str(font));plt.rcParams.update({'font.family':font_manager.FontProperties(fname=font).get_name(),'font.size':10,'axes.unicode_minus':False})
    summary=json.loads((root/'summary.json').read_text())
    if summary['reference_both_methods']!='pooled NTC, whole own NTC guide excluded' or summary['batch_correction_in_main_comparison'] is not False:raise ValueError('Mixed methods cannot be overlaid')
    original=load_tables(root);common=load_tables(root,'_commonRNA')
    # Primary scope fixed before seeing results: intersect the original RNA sets.
    common_n=int(summary['common_RNA_count'])
    scatter(common,output/'01_shared_axes_common_RNA.png',f'我们的数据与论文同图比较｜共同{common_n:,}项RNA')
    scatter(original,output/'03_original_RNA_sets_sensitivity.png','保留各自原RNA集合的复核｜细胞和靶点与主图相同')
    comparisons=backgrounds(root,output/'02_same_different_NTC.png','_commonRNA')
    rows=[]
    for scope,tables in [('common_RNA',common),('original_RNA_sets',original)]:
        for name in DATASETS:
            for method in METHODS:
                t=tables[name,method]
                rows.append(dict(RNA_scope=scope,dataset=name,method=method,targets=len(t),median_target_r=float(t.median_pair_r.median())))
    pd.DataFrame(rows).to_csv(output/'target_summary.csv',index=False)
    pd.DataFrame(comparisons).to_csv(output/'background_summary.csv',index=False)
    provenance=dict(source_summary=summary,main_RNA_scope='intersection of three pre-existing RNA sets; selected before reading effects',
        same_eligible_targets_and_x_between_methods=True,main_axes=dict(x=[0,1],y=[-1,1]),
        publication_scope='Paper preselected Strong; not experiment quality ranking',
        source_files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in root.iterdir() if p.suffix=='.csv'})
    (output/'plot_provenance.json').write_text(json.dumps(provenance,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    return rows


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--source',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    print(json.dumps(render(args.source,args.output),ensure_ascii=False))
