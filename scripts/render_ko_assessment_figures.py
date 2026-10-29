"""Render accepted Mixscape and RNA-observability tables; never reread counts.

These figures describe model RNA labels and observed RNA alignment coverage.
Neither is a plot of DNA editing or functional knockout efficiency.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle, Patch
import numpy as np
import pandas as pd


BLUE = '#32769B'
INK = '#25313A'
GRAY = '#AAB3BA'
PALE = '#EDF0F2'
AMBER = '#BB8E53'
ROLE = ['targeting', 'positive_control']


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_tables(root):
    expected = {
        'official_mixscape_guides_chip01.csv': 'a58699c03c2dc27f34b816de9f9f448f53669814d3b1b6fefe8799919055a4fe',
        'official_mixscape_guides_company.csv': 'cc92052638c5ec2254097cdf7286aa699a77b1298e163a56af845ad99e214002',
        'guide_sequence_observability.csv': '96c8522f99b034b0f96ed8b9a6ad8b6e032ef40847519af79eb2d7f1bb39524d',
    }
    result = {}
    for name, sha in expected.items():
        path = root / 'tables' / name
        if digest(path) != sha:
            raise ValueError(f'Accepted input table differs: {name}')
        result[name] = pd.read_csv(path)
    return result, expected


def setup(font):
    font_manager.fontManager.addfont(str(font))
    family = font_manager.FontProperties(fname=str(font)).get_name()
    plt.rcParams.update({
        'font.family': family, 'font.size': 10.5, 'axes.titlesize': 12,
        'axes.labelsize': 11, 'axes.edgecolor': '#6D777E',
        'text.color': INK, 'axes.labelcolor': INK, 'xtick.color': INK,
        'ytick.color': INK, 'axes.unicode_minus': False,
        'svg.fonttype': 'none', 'savefig.facecolor': 'white',
        'axes.spines.top': False, 'axes.spines.right': False,
        'figure.facecolor': 'white', 'axes.facecolor': 'white',
    })


def save(fig, output, name):
    paths = []
    for suffix in ['png', 'svg']:
        path = output / f'{name}.{suffix}'
        fig.savefig(path, dpi=300, facecolor='white')
        paths.append(path)
    plt.close(fig)
    return paths


def rna_response(tables, output):
    frames = {key: tables[f'official_mixscape_guides_{key}.csv'] for key in ['chip01', 'company']}
    cohort = {}
    for key, g in frames.items():
        if not g.guide.is_unique or len(g) != 1196:
            raise ValueError('Incomplete design-guide table')
        t = g[g.source_control_type.isin(ROLE)]
        q = t.groupby('target_gene_name')[['n_QC_records', 'n_model_evaluated', 'n_model_RNA_response']].sum()
        cohort[key] = dict(observed_targets=int(q.n_QC_records.gt(0).sum()),
                           evaluated_targets=int(q.n_model_evaluated.gt(0).sum()),
                           RNA_response_targets=int(q.n_model_RNA_response.gt(0).sum()),
                           RNA_response_guides=int(t.n_model_RNA_response.gt(0).sum()))
    g = frames['company']
    positive_targets = sorted(g.loc[g.n_model_RNA_response.gt(0), 'target_gene_name'].unique())
    assert len(positive_targets) == 7 and cohort['company']['RNA_response_guides'] == 20
    assert cohort['chip01']['RNA_response_targets'] == 0
    short = {t: t.split('_chr')[0] for t in positive_targets}
    assert len(set(short.values())) == len(short)
    values = np.full((len(positive_targets), 5), np.nan)
    annotations = {}
    source_rows = []
    for i, target in enumerate(positive_targets):
        subset = g[g.target_gene_name.eq(target)]
        for _, row in subset.iterrows():
            suffix = int(row.guide.rsplit('-', 1)[1])
            if not 1 <= suffix <= 5:
                raise ValueError('Unexpected guide suffix for this result')
            j = suffix - 1
            n = int(row.n_model_evaluated)
            if n > 0:
                percentage = 100 * int(row.n_model_RNA_response) / n
                assert np.isclose(percentage/100, row.model_RNA_response_fraction)
                values[i, j] = percentage
                # Keep a visible nonzero fractional percentage, e.g. 0.7%.
                label = f'{percentage:.0f}%' if np.isclose(percentage, round(percentage)) else f'{percentage:.1f}%'
                annotations[(i,j)] = (label, f'n={n}', 'observed')
            else:
                annotations[(i,j)] = ('无记录', '', 'no_records')
            source_rows.append(dict(target=target, display_target=short[target], guide=row.guide,
                                    n_records=n, n_RNA_response=int(row.n_model_RNA_response),
                                    RNA_response_percent=values[i,j] if n else None))
    fig = plt.figure(figsize=(10.8, 6.8))
    ax = fig.add_axes([.23, .16, .59, .65])
    cmap = matplotlib.colormaps['Blues']
    norm = Normalize(0,100)
    for i in range(len(positive_targets)):
        for j in range(5):
            item = annotations.get((i,j))
            if item is None:
                ax.add_patch(Rectangle((j-.5,i-.5),1,1,facecolor='white',edgecolor='#E2E7EA',lw=.7,hatch='///'))
                ax.text(j,i,'未设计',ha='center',va='center',fontsize=9,color='#899299')
            elif item[2] == 'no_records':
                ax.add_patch(Rectangle((j-.5,i-.5),1,1,facecolor=PALE,edgecolor='white',lw=1.8))
                ax.text(j,i,'无记录',ha='center',va='center',fontsize=9,color='#7D878F')
            else:
                value=values[i,j]
                ax.add_patch(Rectangle((j-.5,i-.5),1,1,facecolor=cmap(norm(value)),edgecolor='white',lw=1.8))
                color='white' if value >= 55 else INK
                ax.text(j,i-.13,item[0],ha='center',va='center',fontsize=11,fontweight='normal',color=color)
                ax.text(j,i+.18,item[1],ha='center',va='center',fontsize=8.7,color=color)
    ax.set(xlim=(-.5,4.5),ylim=(len(positive_targets)-.5,-.5))
    ax.set_xticks(range(5),[f'sgRNA {i}' for i in range(1,6)])
    ax.set_yticks(range(len(positive_targets)),[short[t] for t in positive_targets])
    ax.tick_params(axis='both',length=0,pad=9)
    for spine in ax.spines.values():spine.set_visible(False)
    cax=fig.add_axes([.855,.19,.02,.58])
    bar=fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,ticks=[0,25,50,75,100])
    bar.set_label('被模型识别为RNA发生变化的记录比例 (%)',labelpad=11)
    bar.outline.set_visible(False)
    fig.text(.06,.955,'Mixscape结果：有多少记录被识别出RNA变化？',fontsize=15,fontweight='normal',va='top')
    fig.text(.06,.89,f"出现模型RNA响应的靶点：首芯片 0 / {cohort['chip01']['evaluated_targets']}　　公司合并 7 / {cohort['company']['evaluated_targets']}",fontsize=11)
    fig.text(.06,.835,'下图展开公司这7个靶点的全部设计sgRNA，包含0%和无记录的情况。',fontsize=10,color='#57626A')
    fig.text(.06,.067,'每格为模型判定比例，n为参与计算的记录数；PB仅省略坐标后缀。这不是DNA编辑率或功能KO效率。',fontsize=9.2,color='#57626A')
    pd.DataFrame(source_rows).to_csv(output/'RNA_heatmap_display_values.csv',index=False)
    return save(fig,output,'03_Mixscape_RNA响应'), cohort


def sequence_coverage(tables, output):
    table=tables['guide_sequence_observability.csv']
    assert len(table)==1196 and table.guide.is_unique
    t=table[table.source_control_type.isin(ROLE) & table.n_assigned_QC_records.gt(0)].copy()
    assert len(t)==781
    own=t.own_guide_distinct_CB_UB_families_spanning.fillna(0).astype(int)
    capped=t.scan_capped.fillna(False).astype(bool)
    positive=own.gt(0)
    missing_umi=t.coverage_status.eq('own_assigned_reads_span_but_UMI_identity_missing')
    unresolved=t.mapping_status.ne('unique_exact_genomic_placement')
    limited=capped & ~positive
    absent=~positive & ~missing_umi & ~unresolved & ~limited
    masks=[positive,absent,limited,missing_umi,unresolved]
    assert np.equal(np.column_stack([x.to_numpy() for x in masks]).sum(1),1).all()
    counts=[int(x.sum()) for x in masks]
    assert counts==[212,512,46,7,4]
    labels=['有RNA片段跨过目标区域','未观察到这样的RNA片段','读取受限，目前无法判断','有片段，但缺少UMI信息','靶向序列有多个匹配位置']
    exact_depth=own[positive & ~capped]
    bound_depth=own[positive & capped]
    bins=[0,1,2,4,9,19,49,np.inf]
    depth_labels=['1','2','3–4','5–9','10–19','20–49','≥50']
    exact_counts=pd.cut(exact_depth,bins=bins,labels=depth_labels).value_counts(sort=False).to_numpy()
    bound_counts=pd.cut(bound_depth,bins=bins,labels=depth_labels).value_counts(sort=False).to_numpy()
    assert int(exact_counts.sum()+bound_counts.sum())==212
    assert int(own[positive].median())==2 and int(own.ge(10).sum())==17
    fig=plt.figure(figsize=(11.5,5.8))
    gs=fig.add_gridspec(1,2,left=.27,right=.965,bottom=.22,top=.79,width_ratios=[1.08,1],wspace=.32)
    ax=fig.add_subplot(gs[0]);bx=fig.add_subplot(gs[1])
    colors=[BLUE,GRAY,AMBER,'#D6DCE0','#D6DCE0']
    y=np.arange(len(labels))
    ax.barh(y,counts,height=.58,color=colors,edgecolor='none')
    ax.set_yticks(y,labels);ax.invert_yaxis();ax.set_xlim(0,570)
    ax.set_xticks([0,200,400]);ax.set_xlabel('sgRNA数量')
    ax.set_title('A  781条有细胞记录的靶向sgRNA',loc='left',pad=16,fontsize=11.5)
    ax.tick_params(axis='y',length=0,pad=8)
    for i,n in enumerate(counts):ax.text(n+9,i,str(n),ha='left',va='center',fontsize=10.5)
    ax.spines['left'].set_visible(False)
    ax.grid(axis='x',color='#E9ECEF',lw=.6);ax.set_axisbelow(True)
    x=np.arange(len(depth_labels))
    bx.bar(x,exact_counts,width=.72,color=BLUE,label='未触及读取上限')
    bx.bar(x,bound_counts,width=.72,bottom=exact_counts,color='#C6DCE8',edgecolor=BLUE,lw=.6,hatch='///',label='已触及上限：数量为下限')
    bx.set_xticks(x,depth_labels);bx.set_xlabel('每条sgRNA的UMI计数\n（同一细胞条码内去重）')
    bx.set_ylabel('sgRNA数量');bx.set_ylim(0,112)
    bx.set_title('B  这212条，各有多少个UMI？',loc='left',pad=16,fontsize=11.5)
    bx.grid(axis='y',color='#E9ECEF',lw=.6);bx.set_axisbelow(True)
    for i,n in enumerate(exact_counts+bound_counts):bx.text(i,n+2,str(int(n)),ha='center',va='bottom',fontsize=9.5)
    bx.legend(frameon=False,fontsize=8.5,loc='upper right')
    fig.text(.06,.955,'RNA测序片段覆盖了sgRNA对应的区域吗？',fontsize=15,fontweight='normal',va='top')
    fig.text(.06,.88,'公司合并数据　｜　212条有跨区域RNA片段；这212条的UMI计数中位数为2，17条达到10。',fontsize=10.8)
    fig.text(.06,.104,'只统计跨过完整20nt区域、且来自该sgRNA对应细胞记录的读段；同一细胞条码内相同UMI只计一次。',fontsize=9.2,color='#57626A')
    fig.text(.06,.054,'212条中18条触及读取上限，计数仅为已观察到的下限。有RNA覆盖不代表已经发生编辑。',fontsize=9.2,color='#57626A')
    category_table=pd.DataFrame(dict(category=labels,n_guides=counts))
    category_table.to_csv(output/'sequence_display_categories.csv',index=False)
    pd.DataFrame(dict(evidence_bin=depth_labels,uncapped=exact_counts,capped_lower_bound=bound_counts)).to_csv(output/'sequence_display_depth.csv',index=False)
    details=dict(observed_targeting_guides=781,coverage_categories=dict(zip(labels,counts)),
                 positive_guides=212,positive_capped=18,positive_median_families=2,
                 positive_ge10=17,depth_bins=depth_labels,uncapped_depth=exact_counts.tolist(),
                 capped_lower_bound_depth=bound_counts.tolist())
    return save(fig,output,'04_靶向区域测序覆盖'), details


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--input-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--font',type=Path,default=Path('C:/Windows/Fonts/msyh.ttc'))
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    setup(args.font);tables,hashes=source_tables(args.input_root)
    paths,cohort=rna_response(tables,args.output)
    new,coverage=sequence_coverage(tables,args.output);paths+=new
    proof=dict(plot_source_sha256=digest(__file__),input_SHA256=hashes,
               reuses_saved_results=True,raw_counts_reread=False,new_scientific_run=False,
               RNA_cohorts=cohort,sequence_coverage=coverage,
               files={p.name:dict(bytes=p.stat().st_size,sha256=digest(p)) for p in paths})
    (args.output/'plot_proof.json').write_text(json.dumps(proof,indent=2,ensure_ascii=False),encoding='utf8')
    print(json.dumps(dict(files=[str(p) for p in paths],RNA_cohorts=cohort,coverage_categories=coverage['coverage_categories']),ensure_ascii=False,indent=2))


if __name__=='__main__':main()
