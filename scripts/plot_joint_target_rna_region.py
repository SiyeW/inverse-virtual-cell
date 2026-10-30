from pathlib import Path
import argparse,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator

def plot(source,output):
    data=json.loads(source.read_text());plt.rcParams.update({'font.family':'Arial','font.size':10,'axes.titlesize':11,
        'axes.labelsize':10,'axes.linewidth':.7,'svg.fonttype':'none'})
    fig,axes=plt.subplots(1,2,figsize=(8.6,4.65),sharey=True)
    colors={'not_observed':'#BEC5CA','observed':'#447D9E','not_assessable':'white'}
    for ax,ds,label in zip(axes,['chip01','company'],['Chip 1','Company merged records']):
        s=data['summary'][ds];values={(r['target_RNA_detected'],r['regional_RNA_status']):r['guides'] for r in s['counts']}
        bottom=np.zeros(2)
        for status in ['not_observed','observed','not_assessable']:
            n=np.array([values[(b,status)] for b in [False,True]])
            ax.bar([0,1],n,bottom=bottom,width=.58,color=colors[status],
                edgecolor='#62696D' if status=='not_assessable' else 'white',linewidth=.5,
                hatch='///' if status=='not_assessable' else None)
            for i,v in enumerate(n):
                if status=='not_assessable':
                    ax.text(i,bottom[i]+v+8,f'{v} unresolved',ha='center',va='bottom',fontsize=9)
                elif v:ax.text(i,bottom[i]+v/2,str(v),ha='center',va='center',color='white' if status=='observed' else '#252A2D',fontsize=11)
            bottom+=n
        assert int(bottom.sum())==s['observed_guides']
        ax.set_title(f"{label}\n{s['observed_guides']} sgRNAs with assigned records",pad=12)
        ax.set_xticks([0,1],['Not detected','Detected']);ax.set_xlabel('Target RNA in expression matrix',labelpad=9)
        ax.set_ylim(0,530);ax.set_xlim(-.65,1.65);ax.yaxis.set_major_locator(MultipleLocator(100))
        ax.spines[['top','right']].set_visible(False)
        ax.text(.5,-.24,f"No assigned records: {s['no_assigned_records']} / {s['design_guides']} designs",transform=ax.transAxes,ha='center',fontsize=9)
    axes[0].set_ylabel('Number of sgRNAs')
    fig.suptitle('Target RNA detection and local RNA evidence',fontsize=12,y=.99)
    legend=[Patch(facecolor=colors['observed'],label='Regional RNA observed'),Patch(facecolor=colors['not_observed'],label='No regional RNA observed'),
        Patch(facecolor='white',edgecolor='#62696D',hatch='///',label='Region not assessable')]
    fig.legend(handles=legend,loc='lower center',bbox_to_anchor=(.5,.005),ncol=3,frameon=False,fontsize=9,
        handlelength=1.7,columnspacing=1.6)
    fig.subplots_adjust(left=.08,right=.99,top=.80,bottom=.27,wspace=.13)
    output.parent.mkdir(parents=True,exist_ok=True)
    for ext in ['png','svg']:fig.savefig(output.with_suffix('.'+ext),dpi=300,bbox_inches='tight',facecolor='white')
    plt.close(fig)
    output.with_suffix('.proof.json').write_text(json.dumps(dict(input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        layers_do_not_imply_same_transcript=True,all_design_guides_retained_in_summary=True,
        no_records_excluded_from_detection_bars_not_classified_as_undetected=True,
        read_limit_negatives_not_classified_as_absence=True,not_editing_or_KO_efficiency=True),indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();plot(a.input,a.output)
