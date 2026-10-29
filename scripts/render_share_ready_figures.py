"""Restrained, share-ready versions of the two accepted result figures.

This module only redraws accepted tables. It never reads count matrices,
changes cohorts, estimates effects, or performs a new statistical test.
"""
from pathlib import Path
import argparse,json,hashlib
from types import SimpleNamespace
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager,ticker
from matplotlib.lines import Line2D
from render_guide_variation import validate_points
from render_target_rna import select_points

DATASETS=['chip01','company','kolf']
LABELS=['Chip 1','Company merged','KOLF2.1J Strong']
TITLES=['Chip 1 (KO)','Company merged (KO)','KOLF2.1J Strong (CRISPRi)']
KINDS=['within_guide','same_target','different_target','NTC']
COLORS=['#496E8B','#A47C55','#4E8575']
STYLE={'font.family':'Arial','font.size':9.5,'axes.labelsize':10,
       'axes.titlesize':10.5,'axes.titleweight':'normal','axes.linewidth':.7,
       'axes.spines.top':False,'axes.spines.right':False,'axes.unicode_minus':False,
       'xtick.labelsize':9,'ytick.labelsize':9,'xtick.major.width':.7,
       'ytick.major.width':.7,'xtick.major.size':3,'ytick.major.size':3,
       'legend.fontsize':9,'legend.frameon':False,'svg.fonttype':'none',
       'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white'}

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def numeric_sha(values):return hashlib.sha256(np.asarray(values,dtype='<f8').tobytes()).hexdigest()

def assert_drawn_points(ax,expected):
    drawn=np.vstack([np.asarray(c.get_offsets(),float) for c in ax.collections])
    expected=np.asarray(expected,float)
    assert drawn.shape==expected.shape
    a=np.lexsort((drawn[:,1],drawn[:,0]));b=np.lexsort((expected[:,1],expected[:,0]))
    assert np.array_equal(drawn[a],expected[b]),'Scatter offsets do not equal the source values'

def configure(font):
    if font:font_manager.fontManager.addfont(str(font))
    plt.rcParams.update(STYLE)

def save(fig,out,name):
    for ext in ['png','svg']:fig.savefig(out/(name+'.'+ext),dpi=300,bbox_inches='tight',pad_inches=.06)
    plt.close(fig)

def variation(root,expected_summary,out):
    points=validate_points(root,10);old=pd.read_csv(expected_summary)
    fig,ax=plt.subplots(figsize=(10.8,6.0));fig.subplots_adjust(left=.09,right=.99,bottom=.13,top=.94)
    proof=[];ymax=0
    for j,(name,color) in enumerate(zip(DATASETS,COLORS)):
        for i,kind in enumerate(KINDS):
            values=points[name].loc[points[name].comparison.eq(kind),'mean_absolute_difference'].to_numpy(float)
            previous=old.loc[old.dataset.eq(name)&old.comparison.eq(kind)].iloc[0]
            assert len(values)==int(previous.pairs)
            assert np.allclose([np.median(values),np.quantile(values,.1),np.quantile(values,.9)],previous[['median','p10','p90']].to_numpy(float),rtol=1e-12,atol=1e-12)
            position=i+(j-1)*.24
            box=ax.boxplot([values],positions=[position],widths=.19,patch_artist=True,manage_ticks=False,showfliers=False,
                 medianprops=dict(color=color,lw=1.25),whiskerprops=dict(color=color,lw=.75),capprops=dict(color=color,lw=.75))
            box['boxes'][0].set(facecolor=matplotlib.colors.to_rgba(color,.12),edgecolor=color,linewidth=.75)
            jitter=np.random.default_rng(90210+j*10+i).uniform(-.065,.065,len(values))
            ax.scatter(position+jitter,values,s=4.5,alpha=.18,color=color,edgecolors='none',zorder=2)
            ymax=max(ymax,float(values.max()))
            proof.append(dict(dataset=name,comparison=kind,pairs=len(values),values_sha256=numeric_sha(values),median=float(np.median(values))))
    handles=[Line2D([],[],marker='o',markersize=4,linestyle='none',color=c,label=l) for c,l in zip(COLORS,LABELS)]
    ax.legend(handles=handles,loc='upper right',ncol=3,columnspacing=1.3,handletextpad=.4,borderaxespad=.3)
    ax.set_xticks(range(4));ax.set_xticklabels(['Within sgRNA','Same target','Different targets','NTC'])
    ax.set_xlim(-.55,3.55);ax.set_ylim(0,ymax*1.1)
    ax.set_ylabel('Mean absolute RNA-response difference\n(log2 units)')
    ax.yaxis.set_major_locator(ticker.MultipleLocator(.1));ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
    ax.set_axisbelow(True);ax.grid(axis='y',color='#ededed',linewidth=.45)
    save(fig,out,'01_组内与组间RNA差异')
    return dict(raw_sources={d:sha(root/d/'variation_points.csv') for d in DATASETS},data=proof,
                RNA_features=2272,cells_per_group=10,point_values_and_pair_counts_unchanged=True,jitter_seed_and_width_unchanged=True)

def accepted_targets(root,expected_summary):
    previous=pd.read_csv(expected_summary);selected={};proof=[]
    for d in DATASETS:
        p,good=select_points(pd.read_csv(root/d/'guide_target_RNA.csv'))
        old=previous.loc[previous.dataset.eq(d)].iloc[0]
        assert len(p)==int(old.plotted_guides) and int(good.sum())==int(old.guides_NTC_detection_ge10)
        assert np.isclose(100*p.loc[good,'RNA_remaining_ratio'].median(),old.median_RNA_remaining_percent,atol=1e-11,rtol=1e-12)
        selected[d]=(p,good)
        proof.append(dict(dataset=d,source_sha256=sha(root/d/'guide_target_RNA.csv'),plotted_guides=len(p),
                          targets=int(p.target_gene_id.nunique()),NTC_detection_ge10=int(good.sum()),NTC_detection_lt10=int((~good).sum()),
                          above200=int(p.RNA_remaining_ratio.gt(2).sum()),
                          original_x_sha256=numeric_sha(p.NTC_mean_per_10000),original_y_ratio_sha256=numeric_sha(p.RNA_remaining_ratio),
                          guide_axis_sha256=hashlib.sha256('\n'.join(p.guide.astype(str)).encode()).hexdigest(),
                          positive_points=int(p.RNA_remaining_ratio.gt(0).sum()),zero_points=int(p.RNA_remaining_ratio.eq(0).sum())))
    return selected,proof


def target(root,expected_summary,out):
    """Linear percentage scatter with a separate, true-valued high-value panel."""
    selected,proof=accepted_targets(root,expected_summary)
    all_x=pd.concat([p.NTC_mean_per_10000 for p,_ in selected.values()])
    lo=10**np.floor(np.log10(all_x.min()));hi=10**np.ceil(np.log10(all_x.max()))
    ymax=max(100*p.RNA_remaining_ratio.max() for p,_ in selected.values())
    tail_max=float(np.ceil(ymax/500)*500)
    fig=plt.figure(figsize=(12,6.1));axes=[];tops=[]
    width=.28;gap=.045;bottom=.13;height=width*12/6.1
    for i,(d,title) in enumerate(zip(DATASETS,TITLES)):
        ax=fig.add_axes([.065+i*(width+gap),bottom,width,height]);top=fig.add_axes([.065+i*(width+gap),.735,width,.16])
        axes.append(ax);tops.append(top);p,good=selected[d]
        for mask,color,alpha in [(~good,'#c5c5c5',.55),(good,'#3a5364',.7)]:
            q=p.loc[mask];values=100*q.RNA_remaining_ratio
            for area,use in [(ax,values.le(200)),(top,values.gt(200))]:
                area.scatter(q.loc[use,'NTC_mean_per_10000'],values.loc[use],s=7.5,color=color,alpha=alpha,linewidths=0,zorder=3)
        assert_drawn_points(SimpleNamespace(collections=list(ax.collections)+list(top.collections)),
                            np.column_stack((p.NTC_mean_per_10000,100*p.RNA_remaining_ratio)))
        for area in [ax,top]:
            area.set_xscale('log');area.set_xlim(lo,hi);area.set_yscale('linear');area.grid(False)
            area.set_xticks([10.**e for e in [-4,-2,0,2] if lo<=10.**e<=hi])
            area.xaxis.set_major_formatter(ticker.LogFormatterMathtext(10));area.xaxis.set_minor_locator(ticker.NullLocator())
            area.tick_params(axis='y',labelleft=i==0)
        ax.set_ylim(-5,205);ax.set_yticks([0,50,100,150,200]);ax.axhline(100,color='#777777',lw=.75,ls='--',zorder=1)
        top.set_ylim(200,tail_max);top.set_yticks([1000,3000]);top.tick_params(axis='x',bottom=False,labelbottom=False)
        top.spines['bottom'].set_visible(False);top.set_title(title,pad=8)
        # Mark the change between the two linear ranges; no point is moved.
        mark=.012
        ax.plot([-mark,mark],[1-mark,1+mark],transform=ax.transAxes,color='#555555',lw=.65,clip_on=False)
        top.plot([-mark,mark],[-mark,mark],transform=top.transAxes,color='#555555',lw=.65,clip_on=False)
    axes[0].set_ylabel('Target RNA remaining (%)')
    handles=[Line2D([],[],marker='o',markersize=3.5,linestyle='none',color='#3a5364',label='NTC detection ≥10%'),
             Line2D([],[],marker='o',markersize=3.5,linestyle='none',color='#c5c5c5',label='NTC detection <10%')]
    axes[-1].legend(handles=handles,loc='upper right',fontsize=8.5,handletextpad=.3,borderaxespad=.3)
    fig.supxlabel('Mean normalized target RNA in NTC (UMI/10k)',fontsize=10,y=.035)
    fig.canvas.draw();ratios=[float(ax.get_window_extent().height/ax.get_window_extent().width) for ax in axes]
    assert np.allclose(ratios,1,atol=1e-4)
    save(fig,out,'02_靶RNA剩余量_线性散点')
    return dict(data=proof,display='linear_scatter_with_separate_high_value_panel',main_y_limits=[-5,205],tail_y_limits=[200,tail_max],
                all_drawn_offsets_equal_original_values=True,all_zeros_at_true_zero=True,no_pseudocount_added=True,
                main_panel_height_to_width=ratios,x_limits=[float(lo),float(hi)],values_capped=False)


def ecdf(values):
    """Exact right-continuous empirical CDF, retaining ties and observed zeros."""
    values=np.asarray(values,float)
    assert len(values) and np.isfinite(values).all() and (values>=0).all()
    x,counts=np.unique(values,return_counts=True)
    y=100*np.cumsum(counts)/len(values)
    for threshold in [0.,50.,100.,200.]:
        idx=np.searchsorted(x,threshold,side='right')-1
        drawn=0. if idx<0 else y[idx]
        assert np.isclose(drawn,100*np.count_nonzero(values<=threshold)/len(values),atol=1e-12)
    assert y[-1]==100 and counts.sum()==len(values)
    return x,y


def target_cdf(root,expected_summary,out):
    selected,proof=accepted_targets(root,expected_summary)
    fig,axes=plt.subplots(1,2,figsize=(10.8,4.8),sharey=True)
    fig.subplots_adjust(left=.07,right=.985,bottom=.17,top=.87,wspace=.13)
    curves=[]
    for i,(ax,title,measurable) in enumerate(zip(axes,['All plotted sgRNAs','NTC detection ≥10%'],[False,True])):
        for d,label,color in zip(DATASETS,TITLES,COLORS):
            p,good=selected[d];use=p.loc[good] if measurable else p
            values=100*use.RNA_remaining_ratio.to_numpy(float);x,y=ecdf(values)
            xx=np.r_[-3,x,max(200.,float(x[-1]))+1];yy=np.r_[0,y,100]
            line,=ax.step(xx,yy,where='post',color=color,lw=1.45,label=f'{label} (n={len(values):,})')
            assert np.array_equal(line.get_xdata(),xx) and np.array_equal(line.get_ydata(),yy)
            curves.append(dict(dataset=d,cohort='NTC_detection_ge10' if measurable else 'all_plotted',guides=len(values),
                               original_remaining_percent_sha256=numeric_sha(values),zeros=int(np.count_nonzero(values==0)),
                               percent_at_or_below={str(t):float(100*np.count_nonzero(values<=t)/len(values)) for t in [0,50,100,200]},
                               above200=int(np.count_nonzero(values>200)),exact_ecdf_checked=True))
        ax.set_title(title,pad=9);ax.set_xlim(-3,200);ax.set_ylim(0,102)
        ax.set_xticks([0,50,100,150,200]);ax.set_yticks([0,25,50,75,100]);ax.set_axisbelow(True)
        ax.grid(axis='y',color='#ededed',linewidth=.45);ax.axvline(100,color='#999999',lw=.7,ls='--',zorder=0)
        ax.legend(loc='lower right',fontsize=8.0,handlelength=1.6,handletextpad=.5,borderaxespad=.5)
        ax.set_xlabel('Target RNA remaining (%)')
    axes[0].set_ylabel('sgRNAs at or below this level (%)')
    save(fig,out,'02_靶RNA剩余量_累计比例')
    return dict(data=proof,curves=curves,display='exact_ecdf_linear_percentage_axes',shown_x_range=[0,200],
                above200_retained_in_denominators=True,above200_not_moved_to_200=True,
                left_cohort_all_original_points=True,right_cohort_existing_detection_stratum=True,no_pseudocount_added=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--variation',type=Path);p.add_argument('--variation-summary',type=Path)
    p.add_argument('--target',type=Path,required=True);p.add_argument('--target-summary',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--font',type=Path)
    p.add_argument('--only-target',action='store_true');p.add_argument('--target-style',choices=['linear','cdf','both'],default='linear')
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True);configure(a.font)
    proof=dict(presentation_only=True,new_statistics=False,style=STYLE,source_code_sha256=sha(__file__))
    if not a.only_target:
        assert a.variation is not None and a.variation_summary is not None
        proof['variation']=variation(a.variation,a.variation_summary,a.output)
    if a.target_style in ['linear','both']:proof['target_linear']=target(a.target,a.target_summary,a.output)
    if a.target_style in ['cdf','both']:proof['target_cdf']=target_cdf(a.target,a.target_summary,a.output)
    (a.output/'numeric_invariance.json').write_text(json.dumps(proof,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps(dict(presentation_only=True,target_styles=a.target_style,exports=['PNG 300 dpi','SVG']),indent=2))
