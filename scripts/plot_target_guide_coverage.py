"""Standalone version of the accepted group-meeting figure5 left panel."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

def plot(source,output):
    d=pd.read_csv(source)
    assert len(d)==260 and d.target.is_unique and d.eligible_guides.notna().all()
    counts=[int(d.eligible_guides.eq(v).sum()) for v in [0,1,2]]+[int(d.eligible_guides.ge(3).sum())]
    assert counts==[9,47,59,145]
    plt.rcParams.update({'font.family':'Arial','font.size':10,'axes.linewidth':.7,'svg.fonttype':'none'})
    fig,ax=plt.subplots(figsize=(4.9,3.8))
    ax.bar(range(4),counts,width=.65,color='#447D9E',edgecolor='none')
    for i,n in enumerate(counts):ax.text(i,n+4,str(n),ha='center',va='bottom',fontsize=11)
    ax.set_xticks(range(4),['0','1','2','≥3']);ax.set_ylim(0,165);ax.yaxis.set_major_locator(MultipleLocator(50))
    ax.set_xlabel('sgRNAs per target with ≥10 records');ax.set_ylabel('Number of targets')
    ax.set_title('Company merged records: target coverage',fontsize=11,pad=14)
    ax.spines[['top','right']].set_visible(False)
    fig.text(.5,.015,'260 designed targets; original counts unchanged.',ha='center',fontsize=8.5,color='#55595C')
    fig.subplots_adjust(left=.14,right=.98,bottom=.19,top=.85)
    output.parent.mkdir(parents=True,exist_ok=True)
    for ext in ['png','svg']:fig.savefig(output.with_suffix('.'+ext),dpi=300,bbox_inches='tight',facecolor='white')
    plt.close(fig)
    output.with_suffix('.proof.json').write_text(json.dumps(dict(source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        targets=260,counts_0_1_2_3plus=counts,threshold='10 QC-passed single-guide records per sgRNA',
        original_figure='20261024_company_merged/figures/05_guide_agreement.png left panel',scientific_values_unchanged=True),indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();plot(a.input,a.output)
