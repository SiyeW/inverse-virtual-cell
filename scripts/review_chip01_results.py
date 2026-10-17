"""Describe the completed first pass without changing tests or forcing groups."""
from pathlib import Path
from datetime import datetime
import argparse,json,subprocess
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage,leaves_list
from scipy.spatial.distance import pdist
from analyze_expression_structure import sha256,atomic_json


def descriptive_similarity(values,cells,min_cells=50,max_genes=2000):
    """Coverage-limited descriptive ordering only; no significance or group calls."""
    values=np.asarray(values,dtype=float);cells=np.asarray(cells)
    if values.ndim!=2 or cells.shape!=(len(values),) or not np.isfinite(values).all():raise ValueError('Invalid profiles')
    rows=np.flatnonzero((cells>=min_cells)&(values.std(1)>1e-8))
    if len(rows)<2:return rows,np.empty((0,0)),np.array([],dtype=int)
    variance=values[rows].var(0);columns=np.flatnonzero(variance>1e-12)
    columns=columns[np.argsort(variance[columns],kind='stable')[-min(max_genes,len(columns)):]]
    if len(columns)<2:return rows,np.empty((0,0)),columns
    x=values[rows][:,columns];d=pdist(x,'correlation')
    if not np.isfinite(d).all():raise ValueError('Nonfinite descriptive distance')
    order=leaves_list(linkage(np.clip(d,0,2),'average'))
    return rows[order],np.corrcoef(x[order]),columns


def run(root,output):
    if output.exists():raise FileExistsError('Preserve prior review')
    effects=root/'effects';function=root/'functional_analysis';structure=root/'structure_v2'
    for folder in [effects,function,structure]:
        if json.loads((folder/'status.json').read_text())['stage']!='complete':raise ValueError('Stage incomplete')
    e=pd.read_csv(effects/'energy_results.csv');d=pd.read_csv(effects/'de_summary.csv')
    groups=pd.read_csv(function/'perturbation_groups.csv');s=pd.read_csv(function/'control_sensitivity_summary.csv')
    profiles=np.load(effects/'perturbation_expression_profiles.npz')
    if list(profiles['targets'])!=groups.target.tolist():raise ValueError('Profile/target alignment mismatch')
    output.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(12,9),layout='constrained')
    ax=axes[0,0];im=ax.scatter(e.cells,e.energy_full,c=e.bh_q,vmin=0,vmax=1,cmap='viridis',s=18)
    ax.set(xscale='log',xlabel='Cells per perturbation',ylabel='Full-data energy distance',title='Global effect and coverage');fig.colorbar(im,ax=ax,label='BH q')
    ax=axes[0,1];x=-np.log10(e.bh_q.clip(lower=1e-12));y=-np.log10(e.ntc_sensitivity_bh_q.clip(lower=1e-12))
    ax.scatter(x,y,s=16);threshold=-np.log10(.05);ax.axhline(threshold,color='red',ls='--');ax.axvline(threshold,color='red',ls='--')
    ax.set(xlabel='Primary -log10(BH q)',ylabel='NTC sensitivity -log10(BH q)',title='Same significance rule in both analyses')
    ax=axes[1,0];ax.scatter(d.cells,d.exploratory_de_genes,s=18,c=np.where(d.cells<20,'orange','steelblue'))
    ax.set(xscale='log',xlabel='Cells per perturbation',ylabel='Exploratory DE genes',title='Gene-level changes; orange means <20 cells')
    ax=axes[1,1];ax.hist(s.jaccard,bins=np.linspace(0,1,21),color='steelblue')
    ax.set(xlabel='DE gene-set Jaccard',ylabel='Perturbations',title='Control-exclusion sensitivity (not doublet validation)')
    fig.savefig(output/'effect_and_sensitivity_overview.png',dpi=170);plt.close(fig)
    order,corr,columns=descriptive_similarity(profiles['values'],groups.cells.to_numpy())
    mapping=groups.iloc[order][['target_code','target','cells','energy_bh_q']]
    mapping.to_csv(output/'descriptive_similarity_order.csv',index=False)
    if corr.size:
        fig,ax=plt.subplots(figsize=(12,10),layout='constrained');im=ax.imshow(corr,vmin=-1,vmax=1,cmap='RdBu_r')
        step=max(1,len(order)//40);ticks=np.arange(0,len(order),step);labels=mapping.target_code.to_numpy()[ticks]
        ax.set_xticks(ticks,labels,rotation=90,fontsize=6);ax.set_yticks(ticks,labels,fontsize=6)
        ax.set_title('Descriptive response similarity: >=50 cells, NOT validated functional groups');fig.colorbar(im,ax=ax)
        fig.savefig(output/'descriptive_response_similarity.png',dpi=170);plt.close(fig)
    policy=json.loads((structure/'qc_policy.json').read_text());scope=json.loads((function/'methods_and_scope.json').read_text())
    n_sig=int((e.bh_q<.05).sum());n_alt=int((e.ntc_sensitivity_bh_q<.05).sum())
    n_de=int((d.exploratory_de_genes>0).sum());n_groups=int(groups.loc[groups.expression_effect_group!='unassigned','expression_effect_group'].nunique())
    # These sources were already hash-accepted by the dependent functional stage.
    sources={str(p.relative_to(root)):sha256(p) for p in [effects/'output_manifest.json',function/'output_manifest.json',structure/'output_manifest.json']}
    pathway_hits=pd.read_csv(function/'pathway_enrichment_hits.csv.gz')
    result=dict(time=datetime.now().astimezone().isoformat(),code_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
      primary_energy_significant=n_sig,sensitivity_energy_significant=n_alt,energy_tests=len(e),de_contrasts=len(d),de_contrasts_with_hits=n_de,
      minimum_primary_energy_bh_q=float(e.bh_q.min()),median_de_genes=float(d.exploratory_de_genes.median()),maximum_de_genes=int(d.exploratory_de_genes.max()),
      reported_primary_expression_groups=n_groups,descriptive_similarity_targets=len(order),descriptive_features=len(columns),
      source_manifests_sha256=sources,primary_statistics_or_filters_changed=False,descriptive_view='>=50 cells, all eligible regardless of significance; no new group labels or p-values')
    atomic_json(output/'review_summary.json',result)
    lines=['# 第一张芯片：结果复核与解读','',
      '## 这次到底得到什么','',
      f'本轮在 {policy["main_cells"]:,} 个单sgRNA、RNA质控通过细胞中，完成 {len(d)} 个扰动的差异表达与 {len(e)} 个扰动的整体表达变化比较。',
      f'整体效应检验：主分析 {n_sig}/{len(e)}、排除8条对照后的敏感性 {n_alt}/{len(e)} 达到BH q<0.05。主分析最小q为 {e.bh_q.min():.4g}。不显著不等于证明没有效应，也不能据此判定DNA敲除失败。',
      f'基因层面：{n_de}/{len(d)} 个扰动出现至少一个符合本轮展示条件的基因（每扰动全基因BH q<0.05且绝对log2FC>=0.25）；每扰动中位数 {d.exploratory_de_genes.median():.0f} 个、最多 {d.exploratory_de_genes.max()} 个。它们是单芯片探索性候选，尤其少于20细胞的结果需要谨慎。','',
      '## 为什么没有报告显著功能组','',
      '主功能分组的入口要求整体效应q<0.05、至少20细胞且至少5个差异基因。'+('入口为空，因此没有可报告的主功能组：不是所有扰动之间的相关性都为零，也不是程序没有执行。' if n_sig==0 else f'本轮可报告的描述性主分组数为{n_groups}，仍不当作独立实验验证。'),
      f'另提供 {len(order)} 个至少50细胞扰动的描述性响应谱相关热图与层次排序，用于查看模式。这张图没有增加新的显著性检验或功能组标签，不能据图宣称找到已验证功能模块；不为得到正结果而放宽主流程门槛。',
      ('本轮通路候选文件为空，应如实报告未获得满足预设q/重叠门槛的通路，不强行命名功能。' if pathway_hits.empty else f'本轮有{len(pathway_hits)}个达到预设门槛的通路条目，按扰动/方向查表，不当作靶点功能被实验验证。'),'',
      '## 质控、覆盖度与公司指标一起看','',
      '公司RNA文库中位数14925 UMI、4947基因/细胞；有至少一个可靠protospacer的记录比例约71.63%。这些反映可用表达信息与guide识别覆盖，不是编辑效率。',
      '公司报告RNA测序饱和度约13.54%、guide约99.08%；两者按独立捕获分子的重复程度计数，不是读到多少种设计guide。高guide饱和度不能证明DNA敲除成功。',
      'RNA QC通过22745/25234；本轮单sgRNA候选14289，其中3851对照。至少20/50/100细胞的非对照靶点为183/78/15：覆盖度会限制发现弱效应的能力，但不将它当作已确诊的唯一原因。',
      '低UMI阈值4263、低基因数2350、线粒体比例上限8%，沿用已记录的基线规则；不按不可靠双细胞自动线硬删，不声称残余双细胞已经清零。','',
      '## 直接看哪些文件','',
      '- [效应、覆盖度与对照敏感性图](effect_and_sensitivity_overview.png)',
      '- [描述性响应谱相似性图](descriptive_response_similarity.png)；编号对照见 descriptive_similarity_order.csv。',
      '- [统计与功能初稿](../functional_analysis/04_standard_analysis_report.md)',
      '- [整体效应完整表](../effects/energy_results.csv)、[逐扰动差异基因索引](../effects/de_summary.csv)',
      '- [标准细胞降维图](../structure_v2/expression_umap.png)。细胞簇不是功能组，也没有据此猜细胞类型。','',
      '## 本轮边界','',
      '没有独立生物学重复；DNA indel/frameshift未评估；自定义ORF不猜宿主。这里只支持当前芯片、当前方案下的探索性结论，不等于所有芯片都无效，也不为阴性结果自行追加实验、改模型或无限换工具。',
      '原始输入不变，方法和来源见既有参数文件及本目录review_summary.json。公司指标来源为20261012_baseline/tables/company_metrics_summary.csv，覆盖度/QC见同批报告；这些是历史交付指标，不伪称本轮重新测序。']
    (output/'05_interpretation.md').write_text('\n'.join(lines)+'\n',encoding='utf8')
    atomic_json(output/'output_manifest.json',{p.name:dict(bytes=p.stat().st_size,sha256=sha256(p)) for p in output.iterdir() if p.name!='output_manifest.json'})
    atomic_json(output/'status.json',dict(stage='complete',code_commit=result['code_commit'],source_modified=False))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.root,a.output)
