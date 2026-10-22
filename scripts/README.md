# 脚本入口

**实验芯片默认分析：[`run_core_analysis.py`](run_core_analysis.py)**。步骤、输入契约、换芯片与查证方法见[核心路线](../docs/CORE_ANALYSIS_ROADMAP.md)。

| 用途 | 入口 |
|---|---|
| 多芯片目录/输入回执/前置检查及无BAM配置，绝不启动作业 | `manage_core_study.py init / check / render-cellranger`，见[运行与查证](../docs/CORE_STUDY_OPERATIONS.md) |
| 10x完整输入→全细胞原始counts与RNA QC | `run_core_analysis.py prepare`（调用`analyze_guide_capture_qc.py`） |
| 四类覆盖/细胞质量、NTC RNA效应、等n比较、多guide支持及报告 | `run_core_analysis.py run` |
| 已有核心结果只读验收，不重跑 | `run_core_analysis.py verify` |
| 单独表达结构/双细胞工作流审计 | `analyze_expression_structure.py`、`audit_scrublet_workflow.py`，不是自动删细胞 |
| 原版notebook入口 | `run_analysis_notebook.py`，保留旧回放；新批次优先用prepare |
| 旧试验/旧统计路线 | 已移入[archive](../archive/20261023_legacy_chip_analysis/README.md)，只显式查证/历史回放 |

其余脚本不串入实验QC主线：`preprocess_*`、`prepare_10x_guide_capture.py`、`extract_sparse_h5ad_smoke_subset.py`及`train_contrastivevi.py`/`inverse_retrieval.py`属于独立模型数据/研究路线；`copy_verified_delivery.py`、`finalize_verified_delivery.py`、`audit_finalized_delivery.py`、`migrate_verified_tree.py`属于数据接收/迁移。保留不等于本轮需要运行。
