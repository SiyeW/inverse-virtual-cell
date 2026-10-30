"""Show the complete first-chip guide set, keeping selected NTCs visible."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator,ScalarFormatter

def plot(root,output):
    target=pd.read_csv(root/'guide_screen_summary.csv');nt=pd.read_csv(root/'NTC_guide_summary.csv')
    nt['role']='non_targeting_control';group=pd.concat([target,nt],ignore_index=True)
    energy=pd.read_csv(root/'energy_by_guide.csv')
    assert len(target)==368 and len(nt)==97 and group.guide.is_unique and energy.guide.is_unique
    group=group.merge(energy[['guide','records','energy_unbiased']],on='guide',how='left',validate='one_to_one',suffixes=('_R',''))
    assert group.energy_unbiased.notna().all() and (group.records>=10).all()
    assert np.array_equal(group.RNA_association_candidate.to_numpy(),(group.guide_q_BY<.05).to_numpy())
    assert np.all(group.multiplicity_family_guides==465)
    plt.rcParams.update({'font.family':'Arial','font.size':10,'axes.labelsize':10,'axes.titlesize':11,
        'axes.linewidth':.7,'xtick.major.width':.7,'ytick.major.width':.7,'svg.fonttype':'none'})
    fig,axes=plt.subplots(1,2,figsize=(8.2,3.65),sharex=True,sharey=True)
    all_y=group.energy_unbiased;span=max(float(all_y.max()-all_y.min()),.02)
    lo=min(0,float(all_y.min()))-.06*span;hi=float(all_y.max())+.10*span
    largest=int(group.records.max())*1.15
    for ax,is_nt,color,title in zip(axes,[False,True],['#B97F53','#447D9E'],['Targeting sgRNAs (n = 368)','NTC sgRNAs (n = 97)']):
        sub=group[group.role.eq('non_targeting_control')==is_nt];passed=sub.RNA_association_candidate
        ax.scatter(sub.loc[~passed,'records'],sub.loc[~passed,'energy_unbiased'],s=18,c='#B8BEC3',alpha=.7,lw=0,label='Does not pass')
        ax.scatter(sub.loc[passed,'records'],sub.loc[passed,'energy_unbiased'],s=25,c=color,alpha=.95,lw=.35,edgecolors='white',label='Passes RNA-association rule')
        ax.axhline(0,color='#A0A4A7',lw=.7,zorder=0)
        ax.set_xscale('log');ax.set_xlim(9,largest);ax.set_ylim(lo,hi)
        ax.set_xticks([v for v in [10,20,50,100,200,500] if v<largest]);ax.xaxis.set_major_formatter(ScalarFormatter());ax.xaxis.set_minor_locator(NullLocator())
        ax.spines[['top','right']].set_visible(False);ax.set_title(title,pad=10)
        ax.set_xlabel('Cell barcodes per sgRNA')
        ax.text(.98,.96,f'{int(passed.sum())} / {len(sub)} pass',ha='right',va='top',transform=ax.transAxes,fontsize=10)
        ax.legend(loc='upper left',bbox_to_anchor=(0,-.25),frameon=False,fontsize=8.5,handletextpad=.45)
    axes[0].set_ylabel('Difference from NTCs\n(E-distance estimate)')
    fig.suptitle('Exploratory RNA-response screening: first chip',fontsize=12,y=1.01)
    fig.subplots_adjust(left=.10,right=.99,bottom=.28,top=.88,wspace=.10)
    output.parent.mkdir(parents=True,exist_ok=True)
    for ext in ['png','svg']:fig.savefig(output.with_suffix('.'+ext),dpi=300,bbox_inches='tight',facecolor='white')
    plt.close(fig)
    proof=dict(input_files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [root/'guide_screen_summary.csv',root/'NTC_guide_summary.csv',root/'energy_by_guide.csv']},
        guide_points=465,targeting_guides=368,NTC_guides=97,targeting_passing=int(target.RNA_association_candidate.sum()),
        NTC_passing=int(nt.RNA_association_candidate.sum()),fixed_RNA_panel=2000,company_results_substituted=False,
        not_a_DNA_or_functional_KO_call=True)
    output.with_suffix('.proof.json').write_text(json.dumps(proof,indent=2));print(json.dumps(proof))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();plot(a.input,a.output)
