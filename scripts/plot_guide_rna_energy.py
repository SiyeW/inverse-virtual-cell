"""Plot all per-guide tests and annotate the unchanged unsupervised UMAP."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import pandas as pd


def draw(table,summary_file,points_file,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','Noto Sans CJK SC','DejaVu Sans'],'axes.unicode_minus':False,'font.size':12,'svg.fonttype':'none'})
    frame=pd.read_csv(table);summary=json.loads(Path(summary_file).read_text());out=Path(out);out.mkdir(parents=True,exist_ok=True)
    if len(frame)!=785 or not frame.guide.is_unique or summary['tested_guides']!=785:raise ValueError('Whole guide inventory not complete')
    if summary['calibration']['review_required']:raise ValueError('NTC diagnostic needs review before candidates may be plotted')
    target=frame.role.isin(['targeting','positive_control']);ntc=frame.role.eq('non_targeting_control');inter=frame.role.eq('intergenic_control')
    sig=frame.significant_adj.astype(str).eq('True')
    if not (frame.pvalue_adj.gt(0)&frame.pvalue_adj.le(1)).all():raise ValueError('Invalid adjusted P values')
    fig,ax=plt.subplots(figsize=(8.2,6.4));fig.subplots_adjust(left=.13,right=.96,bottom=.22,top=.88)
    y=-np.log10(frame.pvalue_adj.to_numpy())
    for mask,color,label in [(target,'#c58d59',f'靶向sgRNA（{int(target.sum())}条）'),(ntc,'#3d799d',f'NTC sgRNA（{int(ntc.sum())}条）'),(inter,'#8d969e',f'基因间区对照（{int(inter.sum())}条）')]:
        ax.scatter(frame.loc[mask,'energy_unbiased'],y[mask],s=26,c=color,alpha=.7,linewidths=.35,edgecolors='white',label=label)
    ax.axvline(0,color='#a6aaae',lw=.7);ax.axhline(-np.log10(.05),color='#555d63',ls='--',lw=.9)
    ax.text(.98,-np.log10(.05),' 校正P＝0.05',transform=ax.get_yaxis_transform(),ha='right',va='bottom',fontsize=10,color='#51595f')
    ax.set_xlabel('整体RNA差异（样本量校正E-distance）');ax.set_ylabel('统计支持（−log₁₀ 校正P）')
    ax.set_title('每条sgRNA与NTC的整体RNA差异',pad=16);ax.spines[['top','right']].set_visible(False)
    ax.legend(frameon=False,loc='upper left',fontsize=10)
    d=summary['calibration'];n_hit=int((target&sig).sum());n_ctrl=int((ntc&sig).sum())
    fig.text(.5,.09,f'每点一条sgRNA；靶向组中{n_hit}条、NTC中{n_ctrl}条满足本次校正P＜0.05。\nNTC诊断：原始P＜0.05为{d["raw_p_lt05_count"]}/106；未见明显异常不等于已证明全库误选率。\n负E值是有限样本估计的波动，不表示RNA下降；本图只检验RNA状态差异，不是DNA或功能KO效率。',ha='center',va='center',fontsize=9.5,linespacing=1.7)
    for ext in ['png','svg']:fig.savefig(out/('逐sgRNA整体RNA差异.'+ext),dpi=300,facecolor='white')
    plt.close(fig)
    with np.load(points_file,allow_pickle=False) as z:xy=z['UMAP'];roles=z['role'];labels=z['guide'];barcodes=z['cell_barcode']
    if len(set(barcodes))!=71664 or not set(labels)<=set(frame.guide):raise ValueError('UMAP and guide tests do not share full cell axes')
    retained=set(frame.loc[target&sig,'guide']);marked=np.isin(labels,list(retained));isntc=roles=='non_targeting_control'
    colors=np.full(len(xy),'#c4c7ca',dtype='U7');colors[isntc]='#447d9e';colors[marked]='#bf7541'
    order=np.random.default_rng(20261030).permutation(len(xy))
    fig,ax=plt.subplots(figsize=(7.6,6.8));fig.subplots_adjust(left=.12,right=.97,top=.88,bottom=.20)
    ax.scatter(xy[order,0],xy[order,1],c=colors[order],s=2,alpha=.5,linewidths=0,rasterized=True)
    ax.set_aspect('equal',adjustable='box');ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2');ax.spines[['top','right']].set_visible(False)
    ax.set_title('同一RNA状态图：标出差异sgRNA所在的记录',pad=16)
    legend=[('#bf7541',f'携带差异sgRNA（{int(marked.sum()):,}条记录）'),('#447d9e',f'NTC（{int(isntc.sum()):,}条记录）'),('#c4c7ca','其它记录')]
    ax.legend(handles=[Line2D([],[],marker='o',ls='',color=c,markersize=6,label=l) for c,l in legend],frameon=False,fontsize=9.5,loc='upper right')
    fig.text(.5,.085,'沿用未按响应选择的71,664条记录及原UMAP坐标，没有重新拟合来拉开各组。\n橙色表示“携带通过组级差异条件的sgRNA”，不表示该细胞已被判定响应或KO。',ha='center',va='center',fontsize=9.5,linespacing=1.7)
    for ext in ['png','svg']:fig.savefig(out/('差异sgRNA的UMAP位置.'+ext),dpi=300,facecolor='white')
    plt.close(fig)
    proof=dict(records=71664,tested_guides=785,targeting_groups=int(target.sum()),NTC_groups=int(ntc.sum()),targeting_significant=n_hit,NTC_significant=n_ctrl,marked_records=int(marked.sum()),
        points_SHA256=hashlib.sha256(Path(points_file).read_bytes()).hexdigest(),UMAP_coordinates_SHA256=hashlib.sha256(xy.tobytes()).hexdigest(),same_UMAP_no_refit=True,
        source_test_table_SHA256=hashlib.sha256(Path(table).read_bytes()).hexdigest(),screening_unit='sgRNA group, not individual cell or biological replicate')
    (out/'plot_proof.json').write_text(json.dumps(proof,ensure_ascii=False,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--table',type=Path,required=True);p.add_argument('--summary',type=Path,required=True);p.add_argument('--points',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();draw(a.table,a.summary,a.points,a.output)
