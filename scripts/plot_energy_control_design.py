"""Show the completed matched control experiments without efficacy labels."""
from pathlib import Path
import argparse,json


def draw(summary,output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    r=json.loads(Path(summary).read_text());c=r['cases'];out=Path(output);out.mkdir(exist_ok=True)
    assert r['all_cases_fixed50'] and not r['targeting_rows_analyzed']
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','DejaVu Sans'],'axes.unicode_minus':False,'font.size':11,'svg.fonttype':'none'})
    fig,axes=plt.subplots(1,2,figsize=(9,5.6),sharey=True);fig.subplots_adjust(left=.09,right=.97,bottom=.25,top=.79,wspace=.28)
    cases=[['original_shared_real_query','real_independent_A','real_independent_B'],['pooled_shared_A','pooled_shared_B','pooled_independent']]
    labels=[['原来共用的\n小份NTC参照','每组独立抽\nNTC参照A','每组独立抽\nNTC参照B'],['共用原参照\n随机分组1','共用原参照\n随机分组2','独立参照\n随机分组']]
    colors=[['#9b9fa3','#3d799d','#3d799d'],['#9b9fa3','#9b9fa3','#3d799d']]
    for ax,names,ticks,palette,title in zip(axes,cases,labels,colors,['保留真实NTC sgRNA分组','NTC打散后随机组成细胞组']):
        values=[c[n]['raw_P_lt05'] for n in names];assert all(c[n]['tests']==99 for n in names)
        ax.bar(range(3),values,color=palette,width=.62)
        for i,v in enumerate(values):ax.text(i,v+.6,str(v),ha='center',va='bottom',fontsize=12)
        ax.set_xticks(range(3),ticks,fontsize=9.5);ax.set_ylim(0,40);ax.set_title(title,pad=14,fontsize=12);ax.spines[['top','right']].set_visible(False)
    axes[0].set_ylabel('检出RNA差异的比较数（每种方案99次）')
    fig.suptitle('NTC参照检查：参照抽样与真实分组都会影响结果',fontsize=15,y=.96)
    fig.text(.5,.09,'每次都比较50个细胞与50个NTC；统计条件均为原始P＜0.05。\n灰色沿用同一小份参照，蓝色为各组独立抽参照；计数不表示KO或独立生物重复。',ha='center',va='center',fontsize=9.2,linespacing=1.8)
    for ext in ['png','svg']:fig.savefig(out/f'NTC参照最终检查.{ext}',dpi=300,facecolor='white')
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--summary',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();draw(a.summary,a.output)
