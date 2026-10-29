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
        if not ntc:without_ntc+=1
        points.append(dict(guide=r['guide'],chrom=r['chrom'],pos0=r['pos0'],ref=r['ref'],alt=r['alt'],
                           guide_alt=r['guide_alt'],NTC_alt=r['NTC_alt'],guide_n=own,NTC_n=ntc,
                           guide_percent=100*r['guide_alt']/own,NTC_percent=100*r['NTC_alt']/ntc if ntc else None,
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
    shared_bar_limit=max(max(s['counts']) for s in summaries)*1.22+15
    shared_count_limit=max([p['guide_n'] for s in summaries for p in s['points']],default=10)*1.5
    for col,s in enumerate(summaries):
        ax=axs[0,col];ys=np.arange(len(LABELS))
        ax.barh(ys,s['counts'],color=COLORS,height=.63)
        # Hatching marks limited sampling inside the observed categories too.
        for y,(total,limited) in enumerate(zip(s['counts'],s['capped_by_category'])):
            if limited:ax.barh(y,limited,left=total-limited,color='none',edgecolor='#655e52',hatch='///',height=.63,linewidth=.4)
            ax.text(total+6,y,str(total),va='center',fontsize=10)
        ax.set_yticks(ys,LABELS);ax.invert_yaxis();ax.set_xlim(0,shared_bar_limit)
        ax.set_xlabel('sgRNA条数（每条只计入一类）');ax.set_title(s['label']+'：全部'+str(s['design_n'])+'条靶向sgRNA',pad=14)
        ax.spines[['top','right']].set_visible(False)
        ax=axs[1,col];points=s['points']
        # Company results show NTC fractions are zero for most events; plotting
        # those fractions against own fractions collapses the data to a line.
        # Counts show both the available evidence and ALT support, without
        # turning 1/1 observed RNA families into an apparent 100% KO rate.
        totals=Counter((p['guide_n'],p['guide_alt']) for p in points)
        called=Counter((p['guide_n'],p['guide_alt']) for p in points if p['bcftools_candidate_seen'])
        assert sum(totals.values())==len(points)
        for grouped,color,label in [(totals,'#aeb8c0','比对读段中观察到'),(called,'#28708f','变异检测工具也检出')]:
            if grouped:ax.scatter([p[0] for p in grouped],[p[1] for p in grouped],s=[22*n for n in grouped.values()],color=color,alpha=.76,edgecolors='white',linewidths=.4,label=label)
        for (depth,alt),n in totals.items():
            if n>=5:ax.text(depth,alt,str(n),ha='center',va='center',fontsize=8,color='#243440')
        limit=shared_count_limit
        if points:
            ax.plot([.75,limit],[.75,limit],color='#858585',lw=.8,ls='--',zorder=0)
            ax.set_xscale('log');ax.set_yscale('log');ax.set_xlim(.75,limit);ax.set_ylim(.75,limit)
            ticks=[x for x in [1,2,5,10,20,50,100,200,500,1000,2000,5000,10000] if x<=limit]
            if len(ticks)>8:ticks=[x for x in ticks if x in [1,5,10,50,100,500,1000,5000,10000]]
            ax.set_xticks(ticks,[str(x) for x in ticks]);ax.set_yticks(ticks,[str(x) for x in ticks])
            ax.minorticks_off();ax.legend(loc='upper left',bbox_to_anchor=(0,.84),frameon=False,fontsize=9,markerscale=.7)
        else:ax.text(.5,.5,'本次未观察到由对应sgRNA记录支持的插入/缺失',ha='center',va='center',wrap=True,transform=ax.transAxes,fontsize=10)
        ax.set_xlabel('该位置可检查的RNA分子数');ax.set_ylabel('其中支持该插入/缺失的RNA分子数')
        ax.set_title('每种插入/缺失有多少RNA分子支持？',pad=12)
        ax.set_aspect('equal',adjustable='box');ax.spines[['top','right']].set_visible(False)
        ntc_also=sum(p['NTC_alt']>0 for p in points);ntc_seen=len(points)-s['indel_rows_without_NTC_opportunity']
        ax.text(.02,.98,f"{len(points)}个序列变化组合，涉及{len({p['guide'] for p in points})}条sgRNA\nNTC同位点有覆盖：{ntc_seen}；也见同种变化：{ntc_also}",va='top',transform=ax.transAxes,fontsize=9)
    fig.subplots_adjust(left=.22,right=.98,top=.94,bottom=.15,wspace=.73,hspace=.58)
    fig.text(.5,.066,'每个组合是一条sgRNA靶向区域内的一种插入/缺失；重叠组合合并成圆，圆面积表示组合数，数字标出≥5处的重叠数。\n对角线表示该位置可检查的分子都支持这种变化；两列同尺度。分子按条码/UMI/基因标记去PCR重复，可来自同一条码。\n斜线表示读取触及上限；序列变化仍是RNA候选证据，未检出不代表未编辑，不能据此给出DNA或功能KO效率。',ha='center',va='center',fontsize=9.3,linespacing=1.6)
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
