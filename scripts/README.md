# 脚本入口

**实验芯片默认分析：[`run_core_analysis.py`](run_core_analysis.py)**。步骤、输入契约、换芯片与查证方法见[核心路线](../docs/CORE_ANALYSIS_ROADMAP.md)。

| 用途 | 入口 |
|---|---|
| 多芯片目录/输入回执/前置检查及无BAM配置，绝不启动作业 | `manage_core_study.py init / check / render-cellranger`，见[运行与查证](../docs/CORE_STUDY_OPERATIONS.md) |
| sgRNA计数定义与可选靶点表达标注分开，保留原始靶点/对照身份 | `prepare_guide_counting_reference.py`，见[输入契约](../docs/GUIDE_COUNTING_REFERENCE.md)；不改RNA参考或把缺失靶点改成NTC |
| 接收回执约束下抽取成对FASTQ前缀，仅作小规模工程检查 | `prepare_fastq_prefix_smoke.py`；限制读取条数、逐条配对、输出SHA，不代表正式芯片结果 |
| 10x完整输入→全细胞原始counts与RNA QC | `run_core_analysis.py prepare`（调用`analyze_guide_capture_qc.py`） |
| 四类覆盖/细胞质量、NTC RNA效应、等n比较、多guide支持及报告 | `run_core_analysis.py run` |
| 已有核心结果只读验收，不重跑 | `run_core_analysis.py verify` |
| 已验独立芯片结果→共同RNA范围、NTC背景和跨芯片一致性报告 | `run_core_comparison.py run / verify`，见[总体比较](../docs/CORE_COMPARISON.md) |
| 单独表达结构/双细胞工作流审计 | `analyze_expression_structure.py`、`audit_scrublet_workflow.py`，不是自动删细胞 |
| 原版notebook入口 | `run_analysis_notebook.py`，保留旧回放；新批次优先用prepare |
| 旧试验/旧统计路线 | 已移入[archive](../archive/20261023_legacy_chip_analysis/README.md)，只显式查证/历史回放 |

其余脚本不串入实验QC主线：`preprocess_*`、`prepare_10x_guide_capture.py`、`extract_sparse_h5ad_smoke_subset.py`及`train_contrastivevi.py`/`inverse_retrieval.py`属于独立模型数据/研究路线；`copy_verified_delivery.py`、`finalize_verified_delivery.py`、`audit_finalized_delivery.py`、`migrate_verified_tree.py`属于数据接收/迁移。保留不等于本轮需要运行。
