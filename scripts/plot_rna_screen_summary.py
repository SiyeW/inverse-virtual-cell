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
    background=pd.read_csv(profiles/'different_target_background.csv')
    assert len(target)==368 and len(nt)==97 and len(pairs)==7 and len(coverage)==6
    assert target.RNA_association_candidate.sum()==6 and nt.RNA_association_candidate.sum()==5
    assert pairs.target.nunique()==4 and coverage.compared_pairs.eq(0).sum()==2
    assert np.isfinite(pairs[['same_target_r','different_target_median_r']].to_numpy()).all()
    plt.rcParams.update({'font.family':'Arial','font.size':10,'axes.labelsize':10,'axes.titlesize':11,
        'axes.linewidth':.7,'xtick.major.width':.7,'ytick.major.width':.7,'svg.fonttype':'none'})
    fig,axes=plt.subplots(1,2,figsize=(9.4,4.8),gridspec_kw={'width_ratios':[.78,1.22]})
    counts=np.array([6,5]);totals=np.array([368,97]);rates=100*counts/totals
    ax=axes[0];ax.bar([0,1],rates,width=.52,color=['#B47C51','#447D9E'],edgecolor='none')
    for i,(n,total,rate) in enumerate(zip(counts,totals,rates)):
        ax.text(i,rate+.16,f'{n} / {total}',ha='center',va='bottom',fontsize=10)
    ax.set_xticks([0,1],['Targeting','NTC']);ax.set_ylim(0,6.3);ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.set_ylabel('sgRNAs retained (%)');ax.set_title('SCEPTRE RNA tests → sgRNA selection\nAdded aggregation: Bonferroni + BY',pad=12,fontsize=10)
    ax=axes[1];labels=[]
    for i,row in pairs.iterrows():
        bg=background[(background.anchor==row.anchor)&(background.other_same_target==row.other_same_target)]
        assert len(bg)==row.background_guides
        ax.scatter(bg.different_target_r,np.full(len(bg),i),s=12,c='#A6ACB1',alpha=.65,lw=0,
            label='Other-target comparisons' if i==0 else None)
        ax.scatter([row.same_target_r],[i],s=37,color='#B47C51',edgecolors='white',lw=.5,zorder=3,
            label='Same-target pair' if i==0 else None)
        short=row.target.split('_chr')[0]
        labels.append(f'{short}: sgRNA {row.anchor.rsplit("-",1)[-1]} vs {row.other_same_target.rsplit("-",1)[-1]}')
    ax.axvline(0,color='#DFE2E4',lw=.6,zorder=0)
    lo=min(float(background.different_target_r.min()),float(pairs.same_target_r.min()),0)-.035
    hi=max(float(background.different_target_r.max()),float(pairs.same_target_r.max()),0)+.035
    ax.set_xlim(lo,hi);ax.set_ylim(len(pairs)-.5,-.6);ax.set_yticks(range(len(pairs)),labels,fontsize=9)
    ax.xaxis.set_major_locator(MultipleLocator(.1));ax.set_xlabel('Correlation of RNA changes (Pearson r)')
    rmin=int(pairs.RNAs_compared.min());rmax=int(pairs.RNAs_compared.max())
    ax.set_title(f'7 sgRNA pairs from all 4 comparable targets\nEach pair compares {rmin:,}–{rmax:,} RNA changes',pad=12,fontsize=10)
    ax.legend(loc='upper left',bbox_to_anchor=(0,-.17),frameon=False,fontsize=8.5,ncol=1,handletextpad=.45)
    for ax in axes:ax.spines[['top','right']].set_visible(False)
    fig.suptitle('First-chip exploratory RNA screening',fontsize=12,y=1.01)
    fig.subplots_adjust(left=.07,right=.99,bottom=.25,top=.79,wspace=.82)
    fig.text(.99,.01,'Pearson r ranges from −1 to 1; the displayed range is zoomed.',ha='right',va='bottom',fontsize=8,color='#55595C')
    output.parent.mkdir(parents=True,exist_ok=True)
    for ext in ['png','svg']:fig.savefig(output.with_suffix('.'+ext),dpi=300,bbox_inches='tight',facecolor='white')
    plt.close(fig)
    files=[screen/'guide_screen_summary.csv',screen/'NTC_guide_summary.csv',profiles/'candidate_same_target_pairs.csv',profiles/'candidate_target_coverage.csv',profiles/'different_target_background.csv']
    proof=dict(input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
        targeting=368,targeting_retained=6,NTC=97,NTC_retained=5,retained_targets=6,
        targets_with_eligible_peer=4,all_anchor_pairs_shown=7,no_second_eligible_guide=2,RNAs_per_pair_range=[rmin,rmax],
        guide_differential_expression_tests_repeated=False,not_KO_efficiency=True,company_data_substituted=False)
    output.with_suffix('.proof.json').write_text(json.dumps(proof,indent=2));print(json.dumps(proof))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--screen',type=Path,required=True);p.add_argument('--profiles',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();plot(a.screen,a.profiles,a.output)
