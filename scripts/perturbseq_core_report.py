"""Portable Markdown report and read-only result acceptance for core-analysis v1."""
import json
from pathlib import Path, PurePosixPath
import re

import numpy as np
import pandas as pd

from perturbseq_core_io import sha256


def write_report(output, summary, roles, support, pairs):
    s = summary
    role_names = dict(non_targeting_control='真正NTC', intergenic_control='基因间区对照',
                      positive_control='来源标注正控', targeting='普通targeting')
    role_rows = '\n'.join(f'| {role_names[r.source_control_type]} | {r.library_guides} | {r.cells} | {r.zero_guides} | {r.eligible_guides} |'
                          for r in roles.itertuples())
    values = pairs.pearson.dropna() if 'pearson' in pairs else pd.Series(dtype=float)
    median = f'{values.median():.4f}' if len(values) else '未定义（没有有效guide对）'
    flag = '**工程前缀测试，不作为本芯片科学结果。**\n\n' if not s['scientific_result'] else ''
    text = f'''# {s['chip_id']}：核心质控与RNA响应分析

{flag}本报告保留完整细胞分母、零覆盖guide和未评价靶点。它描述实验覆盖、细胞质量及RNA响应，不输出已经验证的功能命中。

## 1. 细胞与guide分配

![细胞去向及质控交集](figures/01_technical_qc.png)

输入{s['source_cells']:,}个细胞记录；RNA质控后{s['rna_qc_cells']:,}个，RNA质控与单guide条件同时满足{s['analysis_cells']:,}个。
两项条件是交集，不是可以相加的流水阶段。未可靠分配和多guide仍在全量分母中；本流程不重新筛细胞，也不把高RNA计数直接当双细胞删除。

## 2. 四类guide覆盖与缺失

![按来源角色分层的覆盖](figures/02_four_type_coverage.png)

| 来源角色 | 文库guide数 | 质控后单guide细胞 | 零覆盖guide | 达到描述入口guide |
|---|---:|---:|---:|---:|
{role_rows}

曲线和零覆盖比例回答“设计的guide实际支持了多少细胞”，不是起始库到终点的丰度变化。没有时间点基线时，低覆盖不能直接叫耗竭。来源标注正控不等于本实验RNA表型已通过验证；缺少某类时如实保留为0，不制造该类背景。

## 3. NTC与targeting的细胞质量

![四类细胞RNA及guide质量](figures/03_four_type_cell_quality.png)

按角色直接对比RNA UMI、检测特征数、线粒体比例和guide UMI，并展示RNA质控保留率。箱线图只在显示时省略离群点，没有据此增加细胞过滤。相近的技术指标不能单独证明敲除成功；明显偏移也需要结合覆盖和实验设计解释。

## 4. NTC背景下的RNA log2FC

![等细胞数RNA效应分布与幅度](figures/05_matched_rna_effects.png)

本次固定深度归一化到{s['normalization_target_umi']:g} UMI，使用{s['response_features']:,}/{s['source_rna_features']:,}个RNA特征。
RNA log2FC = log2((guide组归一化均值+1)/(NTC参考均值+1))。参考采用全部{s['reference_cells']:,}个来源确认的NTC细胞（{s['reference_guides']}条guide），排除基因间区对照；评价NTC guide时，移除它自己的**全部**细胞，不让它进入自身参考。

每条合格guide固定抽{s['matched_cells']}个细胞，重复{s['matched_repeats']}次。左图用第一次固定抽样，RNA特征先在每条guide内等权；右图用每条guide多次抽样的RMS中位数，避免简单把样本数差别当效应强弱。
注意中心、宽度和尾部，不强求NTC服从正态分布。RNA特征、guide和重复抽样均不是独立生物学重复；RMS不是p值或命中阈值。唯一、明确对应的直接靶RNA从下游响应比较中屏蔽；同名多RNA或自定义ORF不猜测宿主基因。

## 5. RNA表达是否偏离参考对角线

![等细胞数RNA表达均值对比](figures/06_rna_mean_scatter.png)

横纵轴都是相同RNA特征的log1p归一化均值，分别来自不重叠NTC参考和固定大小组。NTC/基因间区示例按最高覆盖预选，标注正控逐靶点展示，不按显著性挑例子。
这是RNA表达比较，**不是Day 0/Day 10 sgRNA丰度图**。大多数点靠近对角线说明这些均值接近，不等于统计检验校准良好或DNA编辑率已知。

## 6. 同靶点多sgRNA支持

![靶点可评价覆盖及guide对响应一致性](figures/07_target_support.png)

全部{s['targets']}个靶点中，{s['targets_assessable']}个具备可描述的多guide响应，余下保留具体未评价原因。有效guide对的RNA log2FC相关中位数为{median}。
相关性衡量不同guide的RNA响应是否相似，不把同一芯片的guide当生物学重复，也不将低相关或覆盖不足直接解释为靶点无功能。

[逐靶点支持图册](逐靶点支持图册.md)按全部靶点名称顺序编号，展示全部可评估靶点和标注正控；未出图的靶点仍在[完整靶点表](tables/target_support.csv)。图册散点使用细胞最多的两条合格guide，相关矩阵对角线只是自身相关，不能作为重复验证。去掉最大guide后的结果是依赖性敏感性检查，不是独立重复。

## 7. 讨论与当前能回答的问题

本次完成覆盖、技术QC、控制背景RNA效应及多guide支持的有限描述。先依据覆盖确定哪些靶点有评价条件，再结合RNA幅度及多guide方向解释；不单靠曲线外观给整批“PASS”，不自动挑选“最佳阳性”。

原AD/MW渐近检验已在稀疏并列秩NTC比较中出现失准，默认核心流程不再运行这条显著性及后续功能聚类路线。旧诊断和历史结果保留供追查，不代表可靠功能命中，也未被本流程重新校准。能量距离、双细胞诊断和UMAP旧结果不是本次自动重跑的产物；需要时按路线图的独立QC/统计阶段处理，不能复用首芯片数值替代新芯片结果。

DNA indel/frameshift效率未从RNA/guide数据评估；正控预期RNA表型和独立生物学重复设计需实验依据。已有批次/细胞系/时间点等元数据保留在输入和配置，不把芯片自动当独立生物学重复。

技术复查：[本次参数](config.resolved.json)、[摘要与版本](summary.json)、[完整输出哈希](output_manifest.json)、[逐guide响应](tables/guide_effect_metadata.csv)、[固定RNA特征](tables/response_features.csv)。这些是复现凭据，不替代科学限制。
'''
    (output/'分析报告.md').write_text(text, encoding='utf8')


def verify_output(output):
    """Accept saved artifacts without reopening RNA source or computing new effects."""
    output = Path(output).resolve()
    if (output/'owner.lock').exists():
        raise ValueError('Run still owned; inspect writer before accepting')
    status = json.loads((output/'status.json').read_text(encoding='utf8'))
    if status['stage'] != 'complete':
        raise ValueError('Run did not complete; no automatic resume')
    manifest = json.loads((output/'output_manifest.json').read_text(encoding='utf8'))
    for rel, record in manifest.items():
        p = PurePosixPath(rel)
        path = output/rel
        if p.is_absolute() or '..' in p.parts or path.is_symlink() or not path.resolve().is_relative_to(output):
            raise ValueError('Unsafe result manifest path')
        if path.stat().st_size != record['bytes'] or sha256(path) != record['sha256']:
            raise ValueError('Result differs: ' + rel)
    actual = {p.relative_to(output).as_posix() for p in output.rglob('*') if p.is_file()}
    if actual != set(manifest) | {'output_manifest.json', 'status.json'}:
        raise ValueError('Unexpected or missing artifacts; inspect before acceptance')
    prepared = status.get('phase') == 'prepare'
    names = ['reports/01_coverage_and_stage1_qc.md'] if prepared else ['分析报告.md', '逐靶点支持图册.md']
    for name in names:
        text = (output/name).read_text(encoding='utf8')
        for match in re.finditer(r'\]\(([^)]+)\)', text):
            rel = match.group(1)
            if rel not in manifest and rel != 'output_manifest.json':
                raise ValueError('Report link not in verified manifest: ' + rel)
    from PIL import Image
    for rel in manifest:
        if rel.endswith('.png'):
            with Image.open(output/rel) as image:
                image.verify()
    if prepared:
        return dict(accepted=True, phase='prepare', verified_files=len(manifest),
            raw_source_reopened=False, analysis_rerun=False, fresh_source_hash=False,
            output_manifest_sha256=sha256(output/'output_manifest.json'),
            method='Saved-file SHA and images; initial H5AD writer recorded reopened row-total checks')
    cells = pd.read_csv(output/'tables/analysis_cells.csv.gz', index_col=0)
    overlay = pd.read_csv(output/'tables/control_annotations.csv')
    records = pd.read_csv(output/'tables/guide_effect_metadata.csv')
    selected = json.loads((output/'matched_cell_selections.json').read_text(encoding='utf8'))
    config = json.loads((output/'config.resolved.json').read_text(encoding='utf8'))
    summary = json.loads((output/'summary.json').read_text(encoding='utf8'))
    if cells.index.has_duplicates or cells.source_row.duplicated().any():
        raise ValueError('Duplicate analysis cell identity')
    if int(overlay.single_guide_after_stage1_qc.sum()) != len(cells):
        raise ValueError('Guide coverage denominator differs')
    source_rows = cells.set_index('source_row').guide_id.astype(str)
    for choice in selected:
        rows = choice['source_rows']
        if len(rows) != config['matched_cells'] or len(set(rows)) != len(rows) or not source_rows.loc[rows].eq(choice['guide']).all():
            raise ValueError('Invalid deterministic matched selection')
    if len(selected) != len(records)*config['matched_repeats'] or not records.reference_overlap_cells.eq(0).all():
        raise ValueError('Selection coverage/reference overlap differs')
    roles = overlay.set_index('guide').source_control_type
    ntc = cells.guide_id.map(roles).eq('non_targeting_control')
    with np.load(output/'tables/rna_effect_profiles.npz') as profiles:
        if not np.array_equal(profiles['guides'], records.guide.astype(str)):
            raise ValueError('Response row identity differs')
        reconstructed = np.log2((profiles['means']+1)/(profiles['reference_means']+1))
        effects = profiles['log2fc']
        if effects.shape != reconstructed.shape:
            raise ValueError('Response shape differs')
        np.testing.assert_allclose(effects[np.isfinite(effects)], reconstructed[np.isfinite(effects)])
        for row in records.itertuples():
            expected = int((ntc & cells.guide_id.ne(row.guide)).sum())
            if row.reference_cells != expected:
                raise ValueError('Whole-guide NTC leave-out reference count differs')
    if summary['validated_functional_hits'] or summary['inferential_de'] != 'not_run':
        raise ValueError('Unexpected inference in descriptive core output')
    return dict(accepted=True, run_id=summary['run_id'], scientific_result=summary['scientific_result'],
        verified_files=len(manifest), output_manifest_sha256=sha256(output/'output_manifest.json'),
        raw_source_reopened=False, analysis_rerun=False, fresh_source_hash=False,
        method='Saved-file SHA, coverage, effect formula, sampling, NTC reference and image-link audit')
