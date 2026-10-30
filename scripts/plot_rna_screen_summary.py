"""Concise, complete-cohort screen counts and all available anchor comparisons."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

def plot(screen,profiles,output):
    target=pd.read_csv(screen/'guide_screen_summary.csv');nt=pd.read_csv(screen/'NTC_guide_summary.csv')
    pairs=pd.read_csv(profiles/'candidate_same_target_pairs.csv');coverage=pd.read_csv(profiles/'candidate_target_coverage.csv')
    assert len(target)==368 and len(nt)==97 and len(pairs)==7 and len(coverage)==6
    assert target.RNA_association_candidate.sum()==6 and nt.RNA_association_candidate.sum()==5
    assert pairs.target.nunique()==4 and coverage.compared_pairs.eq(0).sum()==2
    assert np.isfinite(pairs[['same_target_r','different_target_median_r']].to_numpy()).all()
    plt.rcParams.update({'font.family':'Arial','font.size':10,'axes.labelsize':10,'axes.titlesize':11,
        'axes.linewidth':.7,'xtick.major.width':.7,'ytick.major.width':.7,'svg.fonttype':'none'})
    fig,axes=plt.subplots(1,2,figsize=(8.1,3.65),gridspec_kw={'width_ratios':[.85,1.15]})
    counts=np.array([6,5]);totals=np.array([368,97]);rates=100*counts/totals
    ax=axes[0];ax.bar([0,1],rates,width=.52,color=['#B47C51','#447D9E'],edgecolor='none')
    for i,(n,total,rate) in enumerate(zip(counts,totals,rates)):
        ax.text(i,rate+.16,f'{n} / {total}',ha='center',va='bottom',fontsize=10)
    ax.set_xticks([0,1],['Targeting','NTC']);ax.set_ylim(0,6.3);ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.set_ylabel('sgRNAs retained (%)');ax.set_title('sgRNAs passing RNA-expression tests',pad=12)
    ax=axes[1];lo=-.06;hi=.18
    assert pairs.same_target_r.between(lo,hi).all() and pairs.different_target_median_r.between(lo,hi).all()
    ax.plot([lo,hi],[lo,hi],color='#A8ACAF',lw=.8,ls='--',zorder=0)
    ax.axhline(0,color='#DFE2E4',lw=.6,zorder=0);ax.axvline(0,color='#DFE2E4',lw=.6,zorder=0)
    ax.scatter(pairs.same_target_r,pairs.different_target_median_r,s=34,color='#B47C51',edgecolors='white',lw=.5)
    ax.set_xlim(lo,hi);ax.set_ylim(lo,hi);ax.xaxis.set_major_locator(MultipleLocator(.05));ax.yaxis.set_major_locator(MultipleLocator(.05))
    ax.set_xlabel('r: same-target sgRNA pair');ax.set_ylabel('Median r: different-target sgRNAs')
    ax.set_title('Same-target follow-up',pad=12)
    ax.text(.97,.96,'All 7 pairs from 4 targets',ha='right',va='top',transform=ax.transAxes,fontsize=9)
    for ax in axes:ax.spines[['top','right']].set_visible(False)
    fig.suptitle('First-chip exploratory RNA screening',fontsize=12,y=1.01)
    fig.subplots_adjust(left=.08,right=.99,bottom=.21,top=.85,wspace=.55)
    fig.text(.99,.01,'Pearson r ranges from −1 to 1; correlation axes are zoomed.',ha='right',va='bottom',fontsize=8,color='#55595C')
    output.parent.mkdir(parents=True,exist_ok=True)
    for ext in ['png','svg']:fig.savefig(output.with_suffix('.'+ext),dpi=300,bbox_inches='tight',facecolor='white')
    plt.close(fig)
    files=[screen/'guide_screen_summary.csv',screen/'NTC_guide_summary.csv',profiles/'candidate_same_target_pairs.csv',profiles/'candidate_target_coverage.csv']
    proof=dict(input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
        targeting=368,targeting_retained=6,NTC=97,NTC_retained=5,retained_targets=6,
        targets_with_eligible_peer=4,all_anchor_pairs_shown=7,no_second_eligible_guide=2,
        guide_differential_expression_tests_repeated=False,not_KO_efficiency=True,company_data_substituted=False)
    output.with_suffix('.proof.json').write_text(json.dumps(proof,indent=2));print(json.dumps(proof))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--screen',type=Path,required=True);p.add_argument('--profiles',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();plot(a.screen,a.profiles,a.output)
