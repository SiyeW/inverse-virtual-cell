"""Direct RNA abundance figure; preserve every value in the companion table."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

DATASETS=['chip01','company','kolf']
TITLES=['我们的首芯片（KO）','公司合并矩阵（KO）','论文 KOLF Strong（CRISPRi预选）']


def select_points(table):
    if table.guide.duplicated().any():raise ValueError('A guide has more than one direct RNA result')
    mapped=table.mapping_status.eq('exact_catalogue_id')
    valid=table.loc[mapped&table.cells.ge(10)&table.NTC_mean_per_10000.gt(0)].copy()
    good=valid.NTC_detected_fraction.ge(.1)
    if valid.RNA_remaining_ratio.isna().any() or (valid.RNA_remaining_ratio<0).any():raise ValueError('Invalid remaining RNA ratio')
    if not np.allclose(valid.RNA_remaining_ratio,valid.mean_per_10000/valid.NTC_mean_per_10000,rtol=1e-12,atol=1e-12):raise ValueError('RNA ratio differs from stored means')
    return valid,good


def render(root,out,font=None):
    if font:font_manager.fontManager.addfont(font);plt.rcParams['font.family']=font_manager.FontProperties(fname=font).get_name()
    else:plt.rcParams['font.family']='Microsoft YaHei'
    plt.rcParams['axes.unicode_minus']=False
    out.mkdir(parents=True,exist_ok=True)
    inputs={name:pd.read_csv(root/name/'guide_target_RNA.csv') for name in DATASETS}
    all_x=pd.concat([select_points(t)[0].NTC_mean_per_10000 for t in inputs.values()])
    xlimits=(10**np.floor(np.log10(all_x.min())),10**np.ceil(np.log10(all_x.max())))
    fig,axes=plt.subplots(1,3,figsize=(16,6.2),sharey=True)
    rows=[];provenance={}
    for name,title,ax in zip(DATASETS,TITLES,axes):
        path=root/name/'guide_target_RNA.csv';table=inputs[name]
        p,good=select_points(table);strong=p.loc[good]
        for mask,color,label,zorder in [(~good,'#b5b5b5','NTC中检出不足10%',2),(good,'#327dad','NTC中检出至少10%',3)]:
            q=p.loc[mask];y=q.RNA_remaining_ratio*100
            sizes=12+14*np.log10(q.cells)
            within=y.le(200)
            ax.scatter(q.loc[within,'NTC_mean_per_10000'],y.loc[within],s=sizes.loc[within],c=color,alpha=.55,linewidths=.3,edgecolors='white',label=label,zorder=zorder)
            ax.scatter(q.loc[~within,'NTC_mean_per_10000'],np.full((~within).sum(),199.),s=sizes.loc[~within],marker='^',c=color,alpha=.75,linewidths=.3,edgecolors='white',zorder=zorder)
        pc=p.loc[p.source_control_type.eq('positive_control')]
        ax.scatter(pc.NTC_mean_per_10000,np.minimum(pc.RNA_remaining_ratio*100,199),s=60,marker='D',facecolors='none',edgecolors='#cf7924',linewidths=1,label='项目文库标注正控',zorder=4)
        ax.axhline(100,c='#333333',lw=1.1,zorder=1);ax.axhline(50,c='#999999',lw=.7,ls='--',zorder=1)
        ax.set_xscale('log');ax.set_ylim(-3,207);ax.set_xlim(*xlimits)
        ax.set_yticks([0,25,50,75,100,125,150,175,200]);ax.set_title(title,fontsize=12,pad=13)
        ax.set_xlabel('NTC中靶RNA的平均量\n（每10,000个总RNA UMI）',fontsize=11)
        ax.grid(axis='y',alpha=.15);ax.spines[['top','right']].set_visible(False)
        top=int((p.RNA_remaining_ratio>2).sum())
        med=100*strong.RNA_remaining_ratio.median() if len(strong) else np.nan
        ax.text(.98,.965,f'{len(p)}条sgRNA，{p.target_gene_id.nunique()}个靶RNA\n蓝点中位剩余量：{med:.1f}%',transform=ax.transAxes,ha='right',va='top',fontsize=10,bbox=dict(facecolor='white',alpha=.87,edgecolor='none'),zorder=5)
        excluded=dict(fewer_than_10=int(table.cells.lt(10).sum()),reference_not_detected=int((table.cells.ge(10)&table.NTC_mean_per_10000.eq(0)).sum()),mapping_missing=int(table.mapping_status.eq('missing_from_full_reference').sum()))
        ax.text(.025,-.24,f'未作点：少于10个细胞 {excluded["fewer_than_10"]}条；NTC未检出 {excluded["reference_not_detected"]}条\n超过200% {top}条，用顶部三角表示；完整值见表。',transform=ax.transAxes,fontsize=9,va='top')
        provenance[name]=dict(input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),plotted_guides=len(p),overflow_gt200=top,excluded=excluded,xlimits=xlimits,horizontal_clipping=False)
        rows.append(dict(dataset=name,plotted_guides=len(p),plotted_targets=p.target_gene_id.nunique(),reliably_detected_definition='NTC nonzero fraction >=0.10; descriptive screen, not assay success cutoff',guides_NTC_detection_ge10=len(strong),median_RNA_remaining_percent=med,above200=top,**excluded))
    axes[0].set_ylabel('sgRNA细胞组的靶RNA均值 / NTC均值（%）',fontsize=11)
    axes[1].legend(loc='lower right',fontsize=9,framealpha=.95)
    fig.suptitle('靶RNA剩下多少？100%表示与NTC相同，50%表示只有一半',fontsize=16,y=.98)
    fig.text(.5,.025,'每点一条sgRNA；使用该sgRNA全部合格细胞。三份数据主图都使用混合NTC，未给论文单独作批次校正。\n灰点的对照表达稀少，比例容易受抽样影响。RNA下降或未检出均不能换算为DNA敲除率；KO与CRISPRi不能按此图直接评优劣。',ha='center',fontsize=10)
    fig.subplots_adjust(left=.07,right=.985,top=.87,bottom=.30,wspace=.15)
    for suffix in ['png','svg']:fig.savefig(out/('靶RNA剩余量三份数据同图.'+suffix),dpi=190)
    plt.close(fig)
    pd.DataFrame(rows).to_csv(out.parent/'figure_summary.csv',index=False)
    (out.parent/'plot_provenance.json').write_text(json.dumps(provenance,indent=2),encoding='utf8')
    return rows


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--font');a=p.parse_args()
    print(json.dumps(render(a.input,a.output,a.font),indent=2))
