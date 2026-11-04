"""Render all 104-target RNA cases with coding labels on unchanged PCA axes."""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from label_coding_rna_records import display_category

COLORS={"Insufficient site information":"#C7CCCF","Reference RNA at frameshift sites":"#4C8C72",
        "Other variant RNA":"#C4925E","Internal frameshift RNA":"#B34B46"}
BLUE="#427F9C"


def save(fig,path,svg=False):
    path.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(path.with_suffix('.png'),dpi=190,bbox_inches='tight',facecolor='white')
    if svg:fig.savefig(path.with_suffix('.svg'),bbox_inches='tight',facecolor='white')
    plt.close(fig)


def target_display(row):
    name=str(row['target'])
    if name.startswith('PB.'):
        parts=str(row['target_gene_id']).split('_')
        return parts[0]+' (ORF '+parts[2]+')'
    return name


def plot_case(dataset,row,arrays,lookup,path):
    xy=arrays['target_PCA'][:,:2];ntc=arrays['NTC_PCA'][:,:2]
    records=[lookup.get((str(g),str(b)),{}) for g,b in zip(arrays['guide'],arrays['cell_barcode'])]
    categories=np.array([display_category(r) for r in records])
    strong=np.array([r.get('frameshift_distinct_UMI_lower_bound',0)>=2 for r in records])
    fig,axes=plt.subplots(1,3,figsize=(13.1,4.2),sharex=True,sharey=True)
    for ax in axes:
        ax.scatter(ntc[:,0],ntc[:,1],s=2,c=BLUE,alpha=.025,rasterized=True)
        ax.set_xlabel(f"PC1 ({arrays['variance_ratio'][0]*100:.1f}%)")
    axes[0].set_ylabel(f"PC2 ({arrays['variance_ratio'][1]*100:.1f}%)")
    guides=sorted(np.unique(arrays['guide']))
    palette=['#A8794E','#548F78','#9076A7','#C5786C','#8E925B','#D0AA55']
    for i,guide in enumerate(guides):
        take=arrays['guide']==guide
        axes[0].scatter(xy[take,0],xy[take,1],s=12,c=palette[i%len(palette)],alpha=.75,
                        label='sgRNA '+guide.rsplit('-',1)[-1]+f' ({take.sum()})',rasterized=True)
    handles=[Line2D([],[],marker='o',color=BLUE,ls='',markersize=4,label=f'NTC ({len(ntc):,})')]
    handles+=axes[0].get_legend_handles_labels()[0]
    axes[0].legend(handles=handles,frameon=False,fontsize=7.5,loc='best')
    for category,color in COLORS.items():
        take=categories==category
        if take.any():
            text={'Reference RNA at frameshift sites':'Site-reference RNA',
                  'Insufficient site information':'Insufficient site information'}.get(category,category)
            axes[1].scatter(xy[take,0],xy[take,1],s=10 if category=='Insufficient site information' else 24,
                            c=color,alpha=.55 if category=='Insufficient site information' else .9,
                            label=text+f' ({take.sum()})',rasterized=True)
    if strong.any():
        axes[1].scatter(xy[strong,0],xy[strong,1],s=50,facecolors='none',edgecolors='black',linewidths=.9,
                        label='Same frameshift: ≥2 distinct UMIs')
    axes[1].legend(frameon=False,fontsize=7,loc='best')
    sc=axes[2].scatter(xy[:,0],xy[:,1],s=11,c=np.log10(arrays['total_RNA']),cmap='viridis',alpha=.8,rasterized=True)
    fig.colorbar(sc,ax=axes[2],fraction=.045,pad=.025).set_label('log10(total RNA UMI)',fontsize=8)
    for ax,title in zip(axes,['sgRNA identity','RNA sequence evidence in the target ORF','RNA measurement amount']):
        ax.set_title(title,fontsize=10)
    fit=str(arrays['PCA_fit_status'])
    if fit in ['insufficient_target_records','constant_target_RNA']:
        for ax in axes:
            ax.clear();ax.set_axis_off();ax.text(.5,.5,'Target-only PCA not estimable\n'+fit.replace('_',' '),
                                               ha='center',va='center',transform=ax.transAxes)
    cohort='Company merged records' if dataset=='company' else 'Chip 1 (single capture)'
    fig.suptitle(f"{cohort} — {target_display(row)}\n{len(xy):,} target records; {len(guides)} sgRNAs; existing coordinates retained",fontsize=12,y=1.04)
    fig.text(.02,-.025,'All records retained; all NTC projected. Green: reference RNA at frameshift sites. Start/stop-region changes are included with other variants.',fontsize=8)
    fig.tight_layout(w_pad=1);save(fig,path)
    expected={
        'frameshift_records':int((categories=='Internal frameshift RNA').sum()),
        'site_reference_records':int((categories=='Reference RNA at frameshift sites').sum()),
        'other_variant_records':int((categories=='Other variant RNA').sum()),
        'insufficient_site_records':int((categories=='Insufficient site information').sum()),
    }
    if any(expected[k]!=int(row[k]) for k in expected):raise ValueError('Statistical and displayed record counts differ')
    return dict(dataset=dataset,target_gene_id=row['target_gene_id'],records=len(xy),**expected,
                target_PCA_sha256=hashlib.sha256(np.ascontiguousarray(arrays['target_PCA']).tobytes()).hexdigest(),
                NTC_PCA_sha256=hashlib.sha256(np.ascontiguousarray(arrays['NTC_PCA']).tobytes()).hexdigest())


def overview(cases,tests,out):
    fig,axes=plt.subplots(1,2,figsize=(11.6,4.7))
    for did,color,marker,label in [('company','#B27E52','o','Company records'),('chip01',BLUE,'^','Chip 1')]:
        data=cases[cases.dataset.eq(did)]
        axes[0].scatter(data.site_reference_records,data.frameshift_records,c=color,marker=marker,s=26,alpha=.65,label=label)
    axes[0].set(xlabel='Records with reference RNA at frameshift sites',ylabel='Records with internal frameshift RNA',
                title='A  Sequence evidence across all 104 targets')
    axes[0].legend(frameon=False,fontsize=8)
    colors=['#CBD0D3','#A9C5D3','#447F9B','#B34B46']
    names=['Too few sequence-informative records','No within-sgRNA comparison','Compared: q ≥ 0.05','Compared: q < 0.05']
    for i,did in enumerate(['company','chip01']):
        data=tests[tests.dataset.eq(did)&tests.comparison.eq('internal_frameshift')]
        if len(data)!=104:raise ValueError('Incomplete target overview')
        n=[int(data.status.eq('too_few_observed_records').sum()),int(data.status.eq('no_within_guide_comparison').sum()),
           int((data.status.eq('tested')&data.q_raw.ge(.05)).sum()),int((data.status.eq('tested')&data.q_raw.lt(.05)).sum())]
        if sum(n)!=104:raise ValueError('Undefined outcome in all-target overview')
        bottom=0
        for count,color in zip(n,colors):
            axes[1].bar(i,count,bottom=bottom,width=.55,color=color)
            if count>=3:axes[1].text(i,bottom+count/2,str(count),ha='center',va='center',color='white' if color in colors[2:] else '#262626')
            elif count:axes[1].text(i+.32,bottom+count/2,str(count),fontsize=9,va='center')
            bottom+=count
        axes[1].text(i,108,f'{n[3]} with q < 0.05',ha='center',fontsize=9,color=colors[3])
    axes[1].set(xticks=[0,1],xticklabels=['Company records','Chip 1'],ylabel='Number of targets',ylim=(0,116),
                title='B  Full-RNA association, conditional on sgRNA')
    axes[1].legend(handles=[Patch(facecolor=c,label=t) for c,t in zip(colors,names)],frameon=False,fontsize=7.5,
                   loc='upper center',bbox_to_anchor=(.5,-.12),ncol=2)
    fig.text(.02,-.12,'All 104 targets/source and all PCA records retained. Missing site evidence is not evidence of failed editing. Tests are RNA associations, not KO calls.',fontsize=8)
    fig.tight_layout(w_pad=2);save(fig,out/'13_target_RNA_state_association',True)


def ecdf(values):
    x=np.sort(np.asarray(values,float));x=x[np.isfinite(x)]
    return x,np.arange(1,len(x)+1)/len(x) if len(x) else np.array([])


def profile_figure(tables,out):
    if not any(len(t) and t.same_target.any() for t in tables.values()):return False
    fig,axes=plt.subplots(1,2,figsize=(10.1,4.1),sharex=True,sharey=True)
    for ax,did,title in zip(axes,['company','chip01'],['Company records','Chip 1']):
        data=tables[did]
        for same,color,label in [(False,'#9AA3A9','Different targets'),(True,'#B34B46','Same target')]:
            values=data.loc[data.same_target.eq(same),'correlation_raw'].dropna().to_numpy()
            x,y=ecdf(values)
            if len(x):
                ax.step(np.r_[-1,x,1],np.r_[0,y,1],where='post',color=color,lw=1.8,label=f'{label} (n = {len(x)})')
            else:ax.plot([],[],color=color,lw=1.8,label=f'{label} (n = 0)')
        ax.axvline(0,color='#B8BEC2',lw=.7,ls='--');ax.set_xlim(-1,1);ax.set_ylim(0,1.03)
        ax.set_title(title);ax.set_xlabel('Correlation of within-sgRNA RNA changes (Pearson r)')
        ax.legend(frameon=False,fontsize=8,loc='upper left');ax.grid(axis='y',color='#E6E8E9',lw=.6);ax.set_axisbelow(True)
    axes[0].set_ylabel('Cumulative proportion of guide pairs')
    fig.suptitle('Across-guide similarity of RNA changes',fontsize=13,y=1.02)
    fig.text(.02,-.035,'All eligible guides and all pairings; no response cutoff or random subset. Each guide uses all internal-frameshift and site-reference records.',fontsize=8)
    fig.text(.02,-.085,'Some means use one record per evidence class. Pair correlations are descriptive and are not independent replicate tests.',fontsize=8,color='#555B60')
    fig.tight_layout(w_pad=2);save(fig,out/'19_cross_guide_RNA_changes',True);return True


def gallery(cases,out):
    entries=cases.to_dict('records')
    for row in entries:row['image']=row['dataset']+'/'+row['target_gene_id']+'.png'
    text='''<!doctype html><html lang="en"><meta charset="utf-8"><title>RNA states and coding evidence: all 104 targets</title>
<style>body{font:15px system-ui,sans-serif;color:#252525;margin:24px;background:white}header,main{max-width:1500px;margin:auto}h1{font-size:24px}p{line-height:1.5}input,select,button{font:inherit;padding:7px;margin:5px 0}#target{width:100%}#figure{width:100%;height:auto;display:block;margin-top:15px;border:1px solid #eee}#details{overflow-wrap:anywhere}button{margin-right:10px}</style>
<header><h1>RNA states and coding evidence: all 104 targets</h1><p>Existing target-only PCA coordinates; all records retained. Red indicates internal frameshift RNA, green reference RNA observed at the same-guide candidate sites, orange other variants, grey insufficient site information. These are RNA observations, not DNA genotypes or functional KO calls.</p>
<select id="cohort"><option value="company">Company merged barcode records</option><option value="chip01">Chip 1 — single capture</option></select> <input id="search" placeholder="Find a target or Gene ID"><select id="target"></select><button id="previous">Previous</button><button id="next">Next</button><p id="details"></p></header><main><img id="figure" alt="Target-level RNA state and coding evidence"></main>
<script>const all=__DATA__;let shown=[];const cohort=document.getElementById('cohort'),list=document.getElementById('target');
function fill(){const term=document.getElementById('search').value.toLowerCase();shown=all.filter(x=>x.dataset===cohort.value&&(x.target+' '+x.target_gene_id).toLowerCase().includes(term)).sort((a,b)=>a.target.localeCompare(b.target));list.replaceChildren();for(const x of shown){const o=document.createElement('option');o.value=x.target_gene_id;o.textContent=x.target;list.append(o)}if(shown.length)show(shown[0].target_gene_id)}
function show(id){const x=shown.find(x=>x.target_gene_id===id);if(!x)return;list.value=id;document.getElementById('figure').src=x.image;document.getElementById('details').textContent=x.target_gene_id+' · '+x.RNA_records+' records: '+x.frameshift_records+' internal frameshift RNA, '+x.site_reference_records+' site-reference RNA, '+x.other_variant_records+' other variant RNA, '+x.insufficient_site_records+' insufficient site information. '+x.frameshift_two_UMI_records+' records support one internal frameshift with ≥2 distinct UMIs.'}
cohort.onchange=fill;document.getElementById('search').oninput=fill;list.onchange=()=>show(list.value);document.getElementById('previous').onclick=()=>{let i=shown.findIndex(x=>x.target_gene_id===list.value);show(shown[(i-1+shown.length)%shown.length]?.target_gene_id)};document.getElementById('next').onclick=()=>{let i=shown.findIndex(x=>x.target_gene_id===list.value);show(shown[(i+1)%shown.length]?.target_gene_id)};fill();show('PB.26158.49_chr2_1_2485_32_236_uORF_ATG');</script></html>'''.replace('__DATA__',json.dumps(entries,ensure_ascii=False).replace('</','<\\/'))
    (out/'case_figures/index.html').write_text(text,encoding='utf-8')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--config',type=Path,required=True)
    cfg=json.loads(parser.parse_args().config.read_text(encoding='utf-8'))
    root=Path(cfg['output']);out=root/'figures';out.mkdir(exist_ok=False)
    if json.loads((root/'status.json').read_text())['stage']!='statistics_complete_pending_figures':
        raise ValueError('Statistics are not complete')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none',
                         'axes.spines.top':False,'axes.spines.right':False})
    labels=json.loads(Path(cfg['labels']).read_text(encoding='utf-8'))
    if hashlib.sha256(Path(cfg['labels']).read_bytes()).hexdigest()!=cfg['labels_sha256']:raise ValueError('Label packet changed')
    receipt=json.loads((root/'PCA_input_receipt.json').read_text())['sha256_during_required_read']
    cases=pd.read_csv(root/'case_summary.csv');proof=[];tests=[];pairs={}
    for dataset in ['company','chip01']:
        rows=cases[cases.dataset.eq(dataset)]
        if len(rows)!=104:raise ValueError('All 104 target plots are required')
        lookup={(r['guide'],r['cell_barcode']):r for r in labels[dataset]['records']}
        for index,(_,row) in enumerate(rows.iterrows(),start=1):
            pca=Path(cfg['PCA_root'])/'RNA'/dataset/(row.target_gene_id+'.npz');raw=pca.read_bytes()
            if hashlib.sha256(raw).hexdigest()!=receipt[dataset+'/'+row.target_gene_id]:raise ValueError('Existing PCA coordinates changed')
            with np.load(io.BytesIO(raw),allow_pickle=False) as arrays:
                proof.append(plot_case(dataset,row,arrays,lookup,out/'case_figures'/dataset/str(row.target_gene_id)))
            if index%10==0 or index==len(rows):
                (root/'status.json').write_text(json.dumps(dict(stage='figures',dataset=dataset,completed=index,total=len(rows))),encoding='utf-8')
        tests.append(pd.read_csv(root/(dataset+'_association.csv')))
        pairs[dataset]=pd.read_csv(root/(dataset+'_guide_pairs.csv'))
        print(json.dumps(dict(stage='figures',dataset=dataset,completed=len(rows))),flush=True)
    overview(cases,pd.concat(tests,ignore_index=True),out)
    pair_plot=profile_figure(pairs,out);gallery(cases,out)
    (out/'plot_proof.json').write_text(json.dumps(dict(code_commit=cfg['commit'],
        all_208_cases_rendered=True,case_count_by_dataset={d:sum(r['dataset']==d for r in proof) for d in ['company','chip01']},
        unchanged_PCA_source_verified_for_every_case=True,paired_guide_profiles_plot=pair_plot,
        no_BAM_or_original_count_container_read=True,cases=proof),ensure_ascii=False,indent=2),encoding='utf-8')
    (root/'status.json').write_text(json.dumps(dict(stage='complete_pending_delivery',cases=208,pair_plot=pair_plot),indent=2),encoding='utf-8')


if __name__=='__main__':main()
