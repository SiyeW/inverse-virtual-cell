# 总体与芯片间比较

入口为 scripts/run_core_comparison.py，读取各芯片已经验收的核心结果，生成一份嵌图报告和可追查的表。它不读取原始FASTQ、BAM或细胞表达计数矩阵。首次采用的历史c82bcf4报告格式可直接读取，不需要重算已有芯片。

## 比较什么

1. 各芯片输入细胞与RNA质控后单guide细胞的覆盖。
2. NTC与targeting的RNA变化幅度：从每芯片已经保存的固定等细胞数抽样中，取共同RNA条目后重新计算RMS log2FC。每guide贡献一个数值，NTC参考始终在芯片内构造并完整排除自身guide。
3. 同一guide跨芯片的变化模式：在共同RNA条目上计算全部合格细胞平均log2FC的相关系数，分别展示NTC和targeting的中位数及可计算guide数量。常量向量记为未定义，不填成零。

共同RNA范围用实际gene ID的交集建立并保存，按ID对齐顺序；不能按数组列号直接比较。各输入必须采用相同固定归一化深度、匹配细胞数以及完整guide/target/control角色目录。原方法刻意排除guide直接靶点的RNA，因此profiles中相应NaN必须与direct_feature_index逐行一致；计算时保留这种排除，并报告实际条目数，其它缺失或无穷值拒绝。不同原始特征范围的旧RMS数值不直接混用。报告是芯片分层的总体对照，不产生混合细胞矩阵、独立生物重复推断、p值或功能hit。

## 输入和验收

配置见 configs/core_comparison.example.json。expected_chips必须完整列出本次要报告的芯片，输入不能缺项或重复。每个输入通过已验收output_manifest.json的SHA绑定；程序复核计算所需的summary、guide元数据、RNA profiles、完整guide目录和target支持表，并拒绝未完成、仍有owner.lock或只完成工程检查的正式输入。

正式运行还需要一份小的计数/参考可比性审查记录，包含compatible、basis和input_manifests（chip_id到本次输入manifest SHA的完整映射），以其SHA写入配置。basis应说明实际核对的参考、计数及实验来源；相同gene ID清单本身不证明参考序列/外显子注释完全等价。这个记录不能靠将布尔值改为true来代替科学核对。

    python scripts/run_core_comparison.py run --config /private/comparison.json
    python scripts/run_core_comparison.py verify --output /results/study/joint/core_comparison_01

--smoke只用于工程检查，输出明确标为非科学结果。运行采用全新外部目录，拒绝覆盖任何旧输出或与输入重叠；失败保留现场，不自动重启。输出包括总体Markdown、三张图、逐芯片/guide/靶点表、共同RNA清单、配置、版本、状态和逐文件SHA清单。verify只验收这些输出。

数据齐备时先完成每张芯片核心分析及验收，再执行本入口。代码与合成测试完成不等于真实总体分析已经完成。
