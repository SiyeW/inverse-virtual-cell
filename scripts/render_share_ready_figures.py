"""Restrained, share-ready versions of the two accepted result figures.

This module only redraws accepted tables. It never reads count matrices,
changes cohorts, estimates effects, or performs a new statistical test.
"""
from pathlib import Path
import argparse,json,hashlib
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

def configure(font):
    if font:font_manager.fontManager.addfont(str(font))
    plt.rcParams.update(STYLE)

def save(fig,out,name):
    for ext in ['png','svg']:fig.savefig(out/(name+'.'+ext),dpi=300,bbox_inches='tight',pad_inches=.06)
    plt.close(fig)

def variation(root,expected_summary,out):
    points=validate_points(root,10);old=pd.read_csv(expected_summary)
    fig,ax=plt.subplots(figsize=(12,3.8));fig.subplots_adjust(left=.075,right=.995,bottom=.20,top=.91)
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

def target(root,expected_summary,out):
    tables={d:pd.read_csv(root/d/'guide_target_RNA.csv') for d in DATASETS};previous=pd.read_csv(expected_summary)
    selected={d:select_points(tables[d]) for d in DATASETS};all_x=pd.concat([v[0].NTC_mean_per_10000 for v in selected.values()])
    lo=10**np.floor(np.log10(all_x.min()));hi=10**np.ceil(np.log10(all_x.max()))
    fig,axes=plt.subplots(1,3,figsize=(12,3.65),sharex=True,sharey=True)
    fig.subplots_adjust(left=.065,right=.995,bottom=.22,top=.87,wspace=.12)
    proof=[]
    for d,title,ax in zip(DATASETS,TITLES,axes):
        p,good=selected[d];old=previous.loc[previous.dataset.eq(d)].iloc[0]
        assert len(p)==int(old.plotted_guides) and int(good.sum())==int(old.guides_NTC_detection_ge10)
        assert np.isclose(100*p.loc[good,'RNA_remaining_ratio'].median(),old.median_RNA_remaining_percent,atol=1e-11,rtol=1e-12)
        for mask,color,alpha in [(~good,'#c5c5c5',.55),(good,'#3a5364',.7)]:
            q=p.loc[mask];values=100*q.RNA_remaining_ratio;inside=values.le(200)
            ax.scatter(q.loc[inside,'NTC_mean_per_10000'],values.loc[inside],s=7.5,color=color,alpha=alpha,linewidths=0,zorder=3)
            ax.scatter(q.loc[~inside,'NTC_mean_per_10000'],np.full((~inside).sum(),200.),s=9,color=color,alpha=alpha,marker='^',linewidths=0,zorder=3)
        ax.axhline(100,color='#777777',lw=.75,ls='--',zorder=1)
        ax.set_title(title,pad=8);ax.set_xscale('log');ax.set_xlim(lo,hi);ax.set_ylim(-4,205)
        ticks=[10.**e for e in [-4,-2,0,2] if lo<=10.**e<=hi]
        ax.set_xticks(ticks);ax.xaxis.set_major_formatter(ticker.LogFormatterMathtext(10));ax.xaxis.set_minor_locator(ticker.NullLocator())
        ax.set_yticks([0,50,100,150,200]);ax.grid(False)
        clipped=int(p.RNA_remaining_ratio.gt(2).sum());assert clipped==int(old.above200)
        proof.append(dict(dataset=d,source_sha256=sha(root/d/'guide_target_RNA.csv'),plotted_guides=len(p),targets=int(p.target_gene_id.nunique()),
                          NTC_detection_ge10=int(good.sum()),NTC_detection_lt10=int((~good).sum()),above200=clipped,
                          original_x_sha256=numeric_sha(p.NTC_mean_per_10000),original_y_ratio_sha256=numeric_sha(p.RNA_remaining_ratio),
                          guide_axis_sha256=hashlib.sha256('\n'.join(p.guide.astype(str)).encode()).hexdigest()))
    axes[0].set_ylabel('Target RNA remaining (%)')
    fig.supxlabel('Mean normalized target RNA in NTC (UMI/10k)',fontsize=10,y=.035)
    handles=[Line2D([],[],marker='o',markersize=3.5,linestyle='none',color='#3a5364',label='NTC detection ≥10%'),
             Line2D([],[],marker='o',markersize=3.5,linestyle='none',color='#c5c5c5',label='NTC detection <10%')]
    axes[-1].legend(handles=handles,loc='upper right',fontsize=8.5,handletextpad=.3,borderaxespad=.3)
    save(fig,out,'02_靶RNA剩余量')
    return dict(data=proof,x_limits=[float(lo),float(hi)],y_display_limits=[-4,205],upper_display_cap_percent=200,
                all_original_guide_points_retained=True,original_remaining_ratios_unchanged=True,uniform_point_size=True,
                positive_controls_retained_without_decorative_overlay=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--variation',type=Path,required=True);p.add_argument('--variation-summary',type=Path,required=True)
    p.add_argument('--target',type=Path,required=True);p.add_argument('--target-summary',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--font',type=Path)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True);configure(a.font)
    proof=dict(variation=variation(a.variation,a.variation_summary,a.output),target=target(a.target,a.target_summary,a.output),
               presentation_only=True,new_statistics=False,style=STYLE,source_code_sha256=sha(__file__))
    (a.output/'numeric_invariance.json').write_text(json.dumps(proof,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps(dict(presentation_only=True,variation_pairs=[x['pairs'] for x in proof['variation']['data'] if x['comparison']=='same_target'],
                         target_guides=[x['plotted_guides'] for x in proof['target']['data']],exports=['PNG 300 dpi','SVG']),indent=2))
