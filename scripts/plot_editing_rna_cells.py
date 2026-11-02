"""Scientific figures for all predeclared guides/targets, with no response selection."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from explore_editing_rna_cells import add_same_indel_support

BLUE = '#377b9d'
RED = '#c53932'
GRAY = '#b9bfc2'
ORANGE = '#bd8857'
GREEN = '#447d64'


def style():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'svg.fonttype': 'none', 'savefig.facecolor': 'white'})


def save(fig, root, stem, svg=False):
    root.mkdir(parents=True, exist_ok=True)
    fig.savefig(root / (stem + '.png'), dpi=190, bbox_inches='tight')
    if svg:
        fig.savefig(root / (stem + '.svg'), bbox_inches='tight')
    plt.close(fig)


def join_qc_counts(data, table):
    """Use actual cohort-specific QC denominators, never company counts for chip1."""
    index = table.set_index(['dataset', 'guide'], verify_integrity=True)
    changes = []
    for dataset in ['company', 'chip01']:
        for row in data[dataset]['guides']:
            actual = int(index.loc[(dataset, row['guide']), 'records'])
            if dataset == 'company' and row['assigned_QC_records'] != actual:
                raise ValueError('Company QC denominator unexpectedly changed')
            if row['assigned_QC_records'] != actual:
                changes.append(dict(dataset=dataset, guide=row['guide'], previous=row['assigned_QC_records'], actual=actual))
            row['assigned_QC_records'] = actual
            if row['observed_records'] > actual:
                raise ValueError('More RNA-supported records than assigned QC records')
    return changes


def primary_summary(data, root):
    a = data['company']; guides = [r for r in a['guides'] if r['primary_197']]
    assert len(guides) == 197
    names = {r['guide'] for r in guides}
    cells = [r for r in a['cells'] if r['guide'] in names]
    variants = [r for r in a['variants'] if r['guide'] in names]
    indels = [r for r in variants if r['kind'] in ['INS', 'DEL']]
    records = sum(r['assigned_QC_records'] for r in guides)
    observed = sum(r['observed_records'] for r in guides)
    counts = [records, observed, sum(r['cut_ref_records'] for r in guides),
              sum(r['indel_records'] for r in guides), sum(r['same_indel_ge2_records'] for r in guides)]
    fig, ax = plt.subplots(2, 2, figsize=(11.2, 8.2))
    values = np.asarray([r['observed_records'] for r in guides])
    bins = [(1,1),(2,4),(5,9),(10,19),(20,10**9)]
    y = [int(((values >= lo) & (values <= hi)).sum()) for lo,hi in bins]
    ax[0,0].bar(np.arange(5), y, color=BLUE, width=.65)
    for i,n in enumerate(y): ax[0,0].text(i,n+1,str(n),ha='center')
    ax[0,0].set(xticks=range(5),xticklabels=['1','2–4','5–9','10–19','20+'],
                xlabel='Barcode records with local RNA, per sgRNA',ylabel='Number of sgRNAs',
                title='A  Local coverage within the 197 sgRNAs')
    ax[0,0].set_ylim(0,max(y)*1.2)
    labels=['All assigned QC records','Any local RNA evidence','Reference RNA across cut ±3 bp',
            'Indel RNA support','Same indel: ≥2 distinct UMIs']
    yy=np.arange(5)[::-1]
    ax[0,1].scatter(counts,yy,s=48,c=[GRAY,BLUE,GREEN,ORANGE,RED],zorder=3)
    for x,z in zip(counts,yy):ax[0,1].annotate(f'{x:,}',(x,z),xytext=(7,0),textcoords='offset points',va='center')
    ax[0,1].set_xscale('log');ax[0,1].set_xlim(10,70000)
    ax[0,1].set(yticks=yy,yticklabels=labels,xlabel='Number of barcode records (log scale)',
                title='B  Barcode-level RNA observations (overlapping groups)')
    ax[0,1].tick_params(axis='y',length=0,labelsize=9)
    for supported,color,label in [(False,ORANGE,'1 RNA molecule'),(True,BLUE,'≥2 RNA molecules')]:
        rows=[r for r in indels if (r['guide_alt']>=2)==supported]
        grouped=Counter((r['cut_offset'],r['length_change']) for r in rows)
        ax[1,0].scatter([k[0] for k in grouped],[k[1] for k in grouped],
                        s=[22+13*(n-1) for n in grouped.values()],alpha=.73,c=color,label=label)
    nt=[r for r in indels if r['NTC_alt']>0]
    ax[1,0].scatter([r['cut_offset'] for r in nt],[r['length_change'] for r in nt],
                    facecolors='none',edgecolors='black',s=75,linewidths=.8,label='Also observed in NTC')
    ax[1,0].axvline(0,color='#555555',lw=.8);ax[1,0].axhline(0,color='#b7b7b7',lw=.6)
    ax[1,0].set(xlabel='Indel start relative to predicted cut (bp)',ylabel='Net sequence length change (bp)',
                title=f'C  All {len(indels)} guide–indel observations')
    ax[1,0].set_xticks([-18,-15,-12,-9,-6,-3,0,3])
    ax[1,0].legend(frameon=False,fontsize=8,loc='upper left')
    remainder=Counter(abs(r['length_change'])%3 for r in indels)
    vals=[remainder[0],remainder[1]+remainder[2]]
    ax[1,1].bar(range(2),vals,color=[GRAY,BLUE],width=.55)
    for i,n in enumerate(vals):ax[1,1].text(i,n+1,str(n),ha='center')
    ax[1,1].set(xticks=range(2),xticklabels=['Multiple of 3 bp','Not a multiple of 3 bp'],xlabel='Absolute indel length',
                ylabel='Guide–indel observations',title='D  Length pattern; reading-frame effect unassigned')
    ax[1,1].set_ylim(0,max(vals)*1.22)
    fig.suptitle('RNA sequence evidence: 197 sgRNAs, 104 targets',fontsize=15,y=1.01)
    fig.text(.02,-.025,'Company merged barcode records. Counts are not DNA editing or functional KO rates. Correct target ORFs remain unverified.',fontsize=9)
    fig.tight_layout(w_pad=2.5,h_pad=2.5)
    save(fig,root,'12_site_RNA_evidence_197',True)
    return dict(primary_sgRNAs=197,targets=104,assigned_records=records,local_RNA_records=observed,
                cut_reference_RNA_records=counts[2],indel_RNA_records=counts[3],indel_ge2_records=counts[4],
                guide_indel_observations=len(indels),indel_supported_guides=len({r['guide'] for r in indels}),
                indel_ge2_observations=sum(r['guide_alt']>=2 for r in indels),
                indel_also_NTC_observations=len(nt),length_remainder=dict(remainder))


def evidence_arrays(cells, guide, barcode):
    evidence={(r['guide'],r['cell_barcode']):r for r in cells}
    rows=[evidence.get((str(g),str(b)),{}) for g,b in zip(guide,barcode)]
    indel=np.array([r.get('max_same_indel_families',0) for r in rows]);ref=np.array([r.get('cut_ref_families',0) for r in rows])
    snv=np.array([r.get('SNV_families',0) for r in rows]);amb=np.array([r.get('ambiguous_families',0) for r in rows])
    category=np.full(len(rows),'No site information',dtype=object)
    category[(ref>0)&(amb==0)]='Reference-window RNA'
    category[snv>0]='Substitution RNA'
    category[indel>0]='Indel RNA'
    return category,indel,ref


def pca_case(dataset, row, arrays, cells, root):
    xy=arrays['target_PCA'][:,:2];ntc=arrays['NTC_PCA'][:,:2]
    categories,indel,ref=evidence_arrays(cells,arrays['guide'],arrays['cell_barcode'])
    fig,axes=plt.subplots(1,3,figsize=(13.1,4.2),sharex=True,sharey=True)
    for ax in axes:
        ax.scatter(ntc[:,0],ntc[:,1],s=2,c=BLUE,alpha=.025,rasterized=True)
        ax.set_xlabel(f"PC1 ({arrays['variance_ratio'][0]*100:.1f}%)")
    axes[0].set_ylabel(f"PC2 ({arrays['variance_ratio'][1]*100:.1f}%)")
    guide_names=sorted(np.unique(arrays['guide']))
    palette=['#a36f3c','#4d8b72','#8b6ba4','#c56d62','#85834d','#d29d43']
    for i,g in enumerate(guide_names):
        mask=arrays['guide']==g
        axes[0].scatter(xy[mask,0],xy[mask,1],s=12,c=palette[i%len(palette)],alpha=.75,label=g.rsplit('-',1)[-1],rasterized=True)
    handles=[Line2D([],[],marker='o',color=BLUE,ls='',markersize=4,label=f'NTC ({len(ntc):,})')]
    h,l=axes[0].get_legend_handles_labels();axes[0].legend(handles=handles+h,labels=[f'NTC ({len(ntc):,})']+['sgRNA '+x for x in l],frameon=False,fontsize=8,loc='best')
    colors={'No site information':GRAY,'Reference-window RNA':GREEN,'Substitution RNA':ORANGE,'Indel RNA':RED}
    for category,color in colors.items():
        mask=categories==category
        if mask.any():axes[1].scatter(xy[mask,0],xy[mask,1],s=10 if category=='No site information' else 24,c=color,
                                     alpha=.55 if category=='No site information' else .9,label=f'{category} ({mask.sum()})',rasterized=True)
    lookup={(r['guide'],r['cell_barcode']):r for r in cells}
    strong=np.array([lookup.get((str(g),str(b)),{}).get('same_indel_distinct_UMI_lower_bound',0)>=2 for g,b in zip(arrays['guide'],arrays['cell_barcode'])])
    if strong.any():axes[1].scatter(xy[strong,0],xy[strong,1],s=50,facecolors='none',edgecolors='black',linewidths=.9,label='Indel: ≥2 RNA molecules')
    axes[1].legend(frameon=False,fontsize=7,loc='best')
    sc=axes[2].scatter(xy[:,0],xy[:,1],s=11,c=np.log10(arrays['total_RNA']),cmap='viridis',alpha=.8,rasterized=True)
    cb=fig.colorbar(sc,ax=axes[2],fraction=.045,pad=.025);cb.set_label('log10(total RNA UMI)',fontsize=8)
    for ax,title in zip(axes,['sgRNA identity','Observed site-RNA evidence','RNA measurement amount']):ax.set_title(title,fontsize=10)
    cohort='Company merged records' if dataset=='company' else 'Chip 1 (single capture)'
    fig.suptitle(f"{cohort} — {row['target']}\n{len(xy):,} target records; {len(guide_names)} sgRNAs; same coordinates in all three panels",fontsize=12,y=1.04)
    fig.text(.02,-.025,'All available QC records retained. Site-RNA labels are observations, not DNA genotypes or functional KO calls.',fontsize=8)
    fig.tight_layout(w_pad=1)
    save(fig,root/dataset,str(row['target_gene_id']))
    return dict(dataset=dataset,target_gene_id=row['target_gene_id'],target=row['target'],RNA_records=len(xy),
                indel_RNA_records=int((indel>0).sum()),indel_ge2_records=int(strong.sum()),
                reference_window_only_records=int((categories=='Reference-window RNA').sum()),
                no_site_information_records=int((categories=='No site information').sum()),
                indel_sgRNAs=len(set(arrays['guide'][indel>0])),
                RNA_amount_PC1_correlation=float(np.corrcoef(np.log1p(arrays['total_RNA']),xy[:,0])[0,1]) if len(xy)>2 else None)


def variant_diagram(row, display, variants, root, slug):
    variants=[v for v in variants if v['guide']==row['guide']]
    variants=sorted(variants,key=lambda v:(v['pos0'],v['kind'],v['alt']))
    start,end=display['start0'],display['end0'];reference=display['reference']
    n=max(1,len(variants));fig,ax=plt.subplots(figsize=(13.2,max(3.3,1.9+.42*n)))
    ax.set_xlim(-15,len(reference)+24);ax.set_ylim(-n-.9,2.1)
    ax.axvspan(row['start0']-start-.5,row['end0']-start-.5,color='#edf0f1',zorder=0)
    for i,base in enumerate(reference):ax.text(i,1,base,ha='center',va='center',family='monospace',fontsize=9,color='#555555')
    ax.text(-1.8,1,'Reference',ha='right',va='center',fontsize=9)
    if row['cut0'] is not None:
        ax.axvline(row['cut0']-start-.5,color=RED,lw=.8,ls='--');ax.text(row['cut0']-start,1.6,'Predicted cut',color=RED,ha='center',fontsize=8)
    ax.text(len(reference)+1,1,'sgRNA RNA\nALT / informative',va='center',fontsize=8)
    ax.text(len(reference)+13,1,'NTC RNA\nALT / informative',va='center',fontsize=8)
    for j,v in enumerate(variants):
        y=-j;pos=v['pos0']-start
        ax.text(-1.8,y,f"{v['kind']} {v['length_change']:+d} bp" if v['kind']!='SNV' else f"{v['ref']}>{v['alt']}",ha='right',va='center',fontsize=9)
        for i,base in enumerate(reference):ax.text(i,y,base,ha='center',va='center',family='monospace',fontsize=9,color='#b7b7b7')
        if v['kind']=='DEL':
            a=pos+1;b=pos+len(v['ref'])
            ax.plot([max(0,a)-.4,min(len(reference),b)-.6],[y,y],color=RED,lw=2.2)
            if b>len(reference):ax.text(len(reference)-.3,y,'→',color=RED,va='center')
        elif v['kind']=='INS':
            ax.text(pos+.5,y,'▼',color=RED,ha='center',va='center',fontsize=10)
            ax.text(pos+.5,y+.28,'+'+v['alt'][len(v['ref']):],color=RED,ha='center',fontsize=8,family='monospace')
        else:ax.text(pos,y,v['alt'],ha='center',va='center',family='monospace',fontsize=9,color=RED,fontweight='bold')
        for prefix,x in [('guide',len(reference)+1),('NTC',len(reference)+13)]:
            denominator=sum(v[prefix+'_'+k] for k in ['alt','ref','other'])
            ax.text(x,y,f"{v[prefix+'_alt']} / {denominator}",va='center',fontsize=9)
    if not variants:ax.text(20,-.2,'No sequence variant observed in the available RNA molecules',ha='center',fontsize=10)
    ax.set_axis_off()
    fig.suptitle(f"{row['guide']}\n{row['chrom']}:{row['start0']+1:,}–{row['end0']:,}; sgRNA region shaded",fontsize=12,y=1.02)
    fig.text(.03,-.015,'Individual variant representations, not complete DNA alleles. Counts are UMI-family observations; ambiguous families excluded from denominators.',fontsize=8)
    save(fig,root,slug)


def original_read_alignment(row, display, root, slug):
    """Display actual stored RNA reads, without inventing full alleles or coverage."""
    from analyze_guided_rna_alleles import Read
    molecules=display['molecules']
    # Several pages preserve all displayed families without unreadably tall PNGs.
    start,end=display['start0'],display['end0'];reference=display['reference']
    def visible(item):
        read=Read(item['SAM'],'display');position=read.start
        for operation,length in read.ops:
            if operation in 'MDN=X':
                if position<end and position+length>start:return True
                position+=length
            elif operation=='I' and start<=position<=end:return True
        return False
    ordered=sorted([m for m in molecules if visible(m)],key=lambda r:(r['cell_barcode'],r['UMI'],r['GX']))
    record_numbers={b:i+1 for i,b in enumerate(sorted({m['cell_barcode'] for m in ordered}))}
    molecule_numbers={}
    for b in record_numbers:
        for i,m in enumerate([m for m in ordered if m['cell_barcode']==b]):molecule_numbers[(b,m['UMI'],m['GX'])]=i+1
    for page,offset in enumerate(range(0,max(1,len(ordered)),45)):
        current=ordered[offset:offset+45]
        fig,ax=plt.subplots(figsize=(12.5,max(3.0,1.5+.22*len(current))))
        ax.set_xlim(-21,len(reference)+1);ax.set_ylim(-len(current)-.7,2.3)
        ax.axvspan(row['start0']-start-.5,row['end0']-start-.5,color='#f1f3f4',zorder=0)
        if row['cut0'] is not None:ax.axvline(row['cut0']-start-.5,color=RED,lw=.7,ls='--',zorder=1)
        for i,base in enumerate(reference):ax.text(i,1.2,base,ha='center',va='center',family='monospace',fontsize=9)
        ax.text(-1.5,1.2,'Reference',ha='right',va='center',fontsize=9)
        for j,item in enumerate(current):
            read=Read(item['SAM'],'display');r=read.start;q=0;y=-j
            label=f"record {record_numbers[item['cell_barcode']]} / RNA {molecule_numbers[(item['cell_barcode'],item['UMI'],item['GX'])]}"
            ax.text(-1.5,y,label,ha='right',va='center',family='monospace',fontsize=7)
            for operation,length in read.ops:
                if operation in 'M=X':
                    a,b=max(start,r),min(end,r+length)
                    if b>a:
                        ax.add_patch(Rectangle((a-start-.45,y-.12),b-a-.1,.24,facecolor='#c3c9cd',edgecolor='none',zorder=2))
                        for pos in range(a,b):
                            ix=q+pos-r;base=read.seq[ix]
                            if base!=reference[pos-start]:
                                ax.add_patch(Rectangle((pos-start-.45,y-.18),.9,.36,facecolor='white',edgecolor='none',zorder=3))
                                color={'A':'#2b8761','C':'#287ba0','G':'#b18020','T':'#c63730'}.get(base,'#666666') if read.qualities[ix]>=20 else '#c2c2c2'
                                ax.text(pos-start,y,base,ha='center',va='center',family='monospace',fontsize=8,color=color,zorder=4)
                    r+=length;q+=length
                elif operation=='D':
                    a,b=max(start,r),min(end,r+length)
                    if b>a:ax.plot([a-start-.45,b-start-.55],[y,y],color='#202020',lw=1,zorder=4)
                    r+=length
                elif operation=='N':
                    a,b=max(start,r),min(end,r+length)
                    if b>a:ax.plot([a-start-.45,b-start-.55],[y,y],color=BLUE,lw=.7,ls=':',zorder=4)
                    r+=length
                elif operation=='I':
                    if start<=r<=end:
                        ax.text(r-start-.5,y,'▼',color='#8b4b92',ha='center',va='center',fontsize=8,zorder=5)
                        ax.text(r-start-.5,y+.21,'+'+str(length),color='#8b4b92',ha='center',fontsize=6,zorder=5)
                    q+=length
                elif operation=='S':q+=length
        if not current:ax.text(20,0,'No representative molecule available',ha='center')
        ax.set_axis_off()
        fig.suptitle(f"Original RNA alignments — {row['guide']}\n{row['chrom']}:{start+1:,}–{end:,}; {len(ordered)} representatives overlap this view; page {page+1}",fontsize=11,y=.99)
        fig.text(.03,-.01,'One original representative read per UMI family. Grey: aligned bases; letters: mismatches; black line: deletion; purple: insertion; blue dotted line: splice.',fontsize=7)
        fig.subplots_adjust(left=.015,right=.995,top=.92,bottom=.04)
        save(fig,root,slug+f'_reads_{page+1}')


def association_overview(case_rows, association, root):
    cases=pd.DataFrame(case_rows);fig,ax=plt.subplots(1,2,figsize=(11.2,4.4))
    for dataset,color,marker,label in [('company',ORANGE,'o','Company records'),('chip01',BLUE,'^','Chip 1')]:
        c=cases[cases.dataset.eq(dataset)]
        ax[0].scatter(c.reference_window_only_records,c.indel_RNA_records,c=color,marker=marker,s=24,alpha=.65,label=label)
        a=association[association.dataset.eq(dataset)&association.comparison.eq('any_indel_RNA')]
        a=a.dropna(subset=['q_raw','q_RNA_amount_residual'])
        if len(a):ax[1].scatter(-np.log10(a.q_raw),-np.log10(a.q_RNA_amount_residual),c=color,marker=marker,s=30,alpha=.8,label=f'{label} ({len(a)} tests)')
    ax[0].set(xlabel='Records with reference-window RNA only',ylabel='Records supporting an RNA indel',title='A  Site-RNA evidence across all 104 targets')
    ax[0].legend(frameon=False,fontsize=8)
    threshold=-np.log10(.05)
    ax[1].axvline(threshold,color='#777777',ls='--',lw=.8);ax[1].axhline(threshold,color='#777777',ls='--',lw=.8)
    ax[1].set(xlabel='RNA-state association: −log10(q)',ylabel='After RNA-amount adjustment: −log10(q)',title='B  Exploratory association, within-sgRNA permutations')
    ax[1].legend(frameon=False,fontsize=8)
    fig.text(.02,-.02,'PCA displays use all available records. Tests compare observed RNA-evidence groups, not KO versus wild type. Dashed lines: q = 0.05.',fontsize=8)
    fig.tight_layout(w_pad=2)
    save(fig,root,'13_target_RNA_state_association',True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);args=p.parse_args();cfg=json.loads(args.config.read_text())
    style();out=Path(cfg['figure_output']);out.mkdir(parents=True,exist_ok=True)
    # Export contains the two sequentially completed sequence cohorts and no source BAM.
    data={d:add_same_indel_support(json.loads((Path(cfg['sequence_output'])/'sequence'/(d+'.json')).read_text())) for d in ['company','chip01']}
    corrections=join_qc_counts(data,pd.read_csv(cfg['joint_table']))
    stats=primary_summary(data,out)
    pd.DataFrame(data['company']['guides']).to_csv(out/'guide_evidence_all_peers.csv',index=False)
    allcases=[]
    for dataset in ['company','chip01']:
        table=pd.read_csv(Path(cfg['RNA_output'])/'RNA'/(dataset+'_targets.csv'))
        for _,row in table.iterrows():
            file=Path(cfg['RNA_output'])/'RNA'/dataset/(row.target_gene_id+'.npz')
            if file.exists():
                with np.load(file,allow_pickle=False) as arrays:allcases.append(pca_case(dataset,row,arrays,data[dataset]['cells'],out/'case_figures'))
        print(json.dumps(dict(stage='PCA_figures',dataset=dataset,figures=len([x for x in allcases if x['dataset']==dataset]))),flush=True)
    cases=pd.DataFrame(allcases);cases.to_csv(out/'case_evidence_summary.csv',index=False)
    association=pd.concat([pd.read_csv(Path(cfg['association_output'])/'association'/(d+'.csv')) for d in ['company','chip01']],ignore_index=True)
    association_overview(allcases,association,out)
    primary=sorted([r for r in data['company']['guides'] if r['primary_197']],key=lambda r:r['guide'])
    display={r['guide']:r for r in data['company']['display']};site_index=[]
    for i,row in enumerate(primary):
        slug=f'G{i+1:03d}';variant_diagram(row,display[row['guide']],data['company']['variants'],out/'site_figures',slug)
        original_read_alignment(row,display[row['guide']],out/'site_figures',slug)
        site_index.append(dict(guide=row['guide'],target_gene_id=row['target_gene_id'],file='site_figures/'+slug+'.png'))
    atomic=lambda p,x:Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf-8')
    atomic(out/'site_index.json',site_index)
    atomic(out/'plot_proof.json',dict(config=cfg,primary_summary=stats,QC_metadata_reconciled=corrections,
                                    case_count_by_dataset=cases.groupby('dataset').size().to_dict(),
                                    site_figures=len(site_index),classification='RNA evidence only; no DNA genotype or KO call'))


if __name__=='__main__':main()
