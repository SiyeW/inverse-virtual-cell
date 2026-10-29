"""Full design-denominator overview, not selected-target illustrations or KO rate."""
from pathlib import Path
from collections import Counter,defaultdict
import argparse,csv,json,hashlib

LABELS=[
    '无合格单sgRNA记录','靶向位置不能唯一确定','读取受限，尚无可用覆盖',
    '未见可用的靶向区域RNA分子','有覆盖，未见插入/缺失',
    '插入/缺失：仅1个RNA分子支持','插入/缺失：至少2个RNA分子支持']
COLORS=['#e3e5e7','#c9cdd1','#dbb982','#abb5c0','#738b9c','#91b8ca','#28708f']


def cohort_summary(design,cohort):
    targets={r['guide'] for r in design if r['source_control_type'] in ['targeting','positive_control']}
    candidate=cohort['candidate'];rows=cohort['guide_rows'];variants=cohort['variant_rows']
    expected={g['guide'] for r in candidate['regions'] for g in r['guides']}
    if len(targets)!=candidate['all_target_guides'] or not expected<=targets:raise ValueError('Whole design identity differs')
    if len(rows)!=len(expected) or {r['guide'] for r in rows}!=expected:raise ValueError('Incomplete or duplicate analyzed guide inventory')
    if len(expected)!=candidate['uniquely_placed_observed']:raise ValueError('Candidate cardinality differs')
    if any(r['guide'] not in expected for r in variants):raise ValueError('Foreign variant guide')
    observed=candidate['observed_target_guides']
    if not len(expected)<=observed<=len(targets):raise ValueError('Impossible source denominator')
    by=defaultdict(list)
    for r in variants:
        if r['kind'] in ['INS','DEL'] and r['guide_alt']>0:by[r['guide']].append(r)
    counts=[len(targets)-observed,observed-len(expected),0,0,0,0,0]
    capped_by_category=[0]*len(counts)
    for r in rows:
        evidence=max([x['guide_alt'] for x in by[r['guide']]],default=0)
        if evidence>=2:index=6
        elif evidence==1:index=5
        elif r['own_RNA_families_observed']>0:index=4
        elif r['scan_capped']:index=2
        else:index=3
        counts[index]+=1
        if r['scan_capped']:capped_by_category[index]+=1
    assert sum(counts)==len(targets)
    points=[];without_ntc=0;NTC_only=0
    for r in variants:
        if r['kind'] not in ['INS','DEL']:continue
        if r['guide_alt']==0:
            if r['NTC_alt']>0:NTC_only+=1
            continue
        own=sum(r['guide_'+k] for k in ['alt','ref','other'])
        ntc=sum(r['NTC_'+k] for k in ['alt','ref','other'])
        if not own:raise ValueError('ALT support without informative RNA family')
        if not ntc:without_ntc+=1;continue
        points.append(dict(guide=r['guide'],chrom=r['chrom'],pos0=r['pos0'],ref=r['ref'],alt=r['alt'],
                           guide_alt=r['guide_alt'],NTC_alt=r['NTC_alt'],guide_n=own,NTC_n=ntc,
                           guide_percent=100*r['guide_alt']/own,NTC_percent=100*r['NTC_alt']/ntc,
                           scan_capped=r['scan_capped'],bcftools_candidate_seen=r['bcftools_candidate_seen']))
    return dict(label=cohort['label'],dataset_id=cohort['dataset_id'],design_n=len(targets),observed_n=observed,
                analyzed_n=len(expected),counts=counts,capped_by_category=capped_by_category,
                points=points,indel_rows_without_NTC_opportunity=without_ntc,
                NTC_only_indel_rows=NTC_only,SNV_rows=sum(r['kind']=='SNV' for r in variants),
                functional_KO_efficiency_assessed=False)


def plot(summaries,directory):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','Noto Sans CJK SC','DejaVu Sans'],
                         'axes.unicode_minus':False,'font.size':11,'svg.fonttype':'none'})
    fig,axs=plt.subplots(2,2,figsize=(13.0,9.5),gridspec_kw={'height_ratios':[1,1.35]})
    for col,s in enumerate(summaries):
        ax=axs[0,col];ys=np.arange(len(LABELS))
        ax.barh(ys,s['counts'],color=COLORS,height=.63)
        # Hatching marks limited sampling inside the observed categories too.
        for y,(total,limited) in enumerate(zip(s['counts'],s['capped_by_category'])):
            if limited:ax.barh(y,limited,left=total-limited,color='none',edgecolor='#655e52',hatch='///',height=.63,linewidth=.4)
            ax.text(total+6,y,str(total),va='center',fontsize=10)
        ax.set_yticks(ys,LABELS);ax.invert_yaxis();ax.set_xlim(0,max(s['counts'])*1.22+15)
        ax.set_xlabel('sgRNA条数（每条只计入一类）');ax.set_title(s['label']+'：全部'+str(s['design_n'])+'条靶向sgRNA',pad=14)
        ax.spines[['top','right']].set_visible(False)
        ax=axs[1,col];points=s['points']
        for threshold,color,label in [(1,'#b6bec5','仅1个RNA分子支持'),(2,'#28708f','至少2个RNA分子支持')]:
            selected=[p for p in points if (p['guide_alt']==1 if threshold==1 else p['guide_alt']>=2)]
            if selected:ax.scatter([p['NTC_percent'] for p in selected],[p['guide_percent'] for p in selected],s=23 if threshold==1 else 34,color=color,alpha=.72,edgecolors='white',linewidths=.4,label=label)
        ax.plot([0,100],[0,100],color='#858585',lw=.8,ls='--',zorder=0)
        ax.set(xlim=(-3,103),ylim=(-3,103),xlabel='NTC中支持同一插入/缺失的RNA分子（%）',ylabel='带对应sgRNA的细胞中支持该变化的RNA分子（%）')
        ax.set_title('靶向区域实际观察到的插入/缺失',pad=12)
        ax.set_aspect('equal',adjustable='box');ax.spines[['top','right']].set_visible(False)
        if points:ax.legend(loc='lower right',frameon=False,fontsize=9)
        else:ax.text(.5,.5,'没有同时具有自身及NTC可观察分子的插入/缺失点',ha='center',va='center',wrap=True,transform=ax.transAxes,fontsize=10)
        ax.text(.02,.98,f"{len(points)}个sgRNA–序列变化组合\n另有{s['indel_rows_without_NTC_opportunity']}个因NTC无覆盖未作点",va='top',transform=ax.transAxes,fontsize=9)
    fig.subplots_adjust(left=.22,right=.98,top=.94,bottom=.15,wspace=.73,hspace=.58)
    fig.text(.5,.066,'每个散点是一条sgRNA靶向区域内的一种插入/缺失；上方表示在对应sgRNA细胞中更常见。\nRNA分子按细胞条码和UMI去重；斜线部分表示读取触及上限。未检出不代表未编辑；图中比例不是DNA或功能KO效率。',ha='center',va='center',fontsize=10,linespacing=1.6)
    directory.mkdir(parents=True,exist_ok=True)
    for suffix in ['png','svg']:fig.savefig(directory/('全体sgRNA序列证据.'+suffix),dpi=300,facecolor='white')
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--design',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    data=json.loads(a.input.read_text());design=list(csv.DictReader(a.design.open()))
    if [c['dataset_id'] for c in data['cohorts']]!=['chip01','company']:raise ValueError('Expected exactly both complete cohorts')
    summaries=[cohort_summary(design,c) for c in data['cohorts']]
    plot(summaries,a.output)
    (a.output/'whole_cohort_plot_data.json').write_text(json.dumps(dict(input_SHA256=hashlib.sha256(a.input.read_bytes()).hexdigest(),design_SHA256=hashlib.sha256(a.design.read_bytes()).hexdigest(),summaries=summaries),ensure_ascii=False,indent=2))
