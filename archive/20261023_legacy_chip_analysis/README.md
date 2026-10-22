# 历史芯片路线（非默认）

2026-10-23从活动入口剥离；共17个脚本、2个配置文件，205,881字节，移动前后逐文件SHA256一致。没有删除源数据、计算产物或旧结果。
原位置、归档位置和完整摘要见[archive_manifest.json](archive_manifest.json)，原版本为`c82bcf418731122de2f3244a4bbb1baed6225d7d`。
`.gitattributes`对本次归档载荷关闭换行转换，跨平台检出仍能对上本机移动时的原字节SHA；原历史Git提交自己的换行/内容另以该提交为准。

当前分析从[核心路线](../../docs/CORE_ANALYSIS_ROADMAP.md)和`python scripts/run_core_analysis.py`开始，不自动进入本目录。

| 类别 | 归档脚本（scripts/下） | 为什么不再默认运行 |
|---|---|---|
| 早期整套及功能聚类 | `run_chip01_standard.py`、`finish_chip01_analysis.py`、`review_chip01_results.py` | 已被后续审计/报告替代；弱分组不是验证模块 |
| 文献路线和数学定义域修正 | `run_chip01_literature.py`、`finish_literature_analysis.py`、`correct_literature_de_domain.py` | 原AD/MW推断的NTC校准未通过；修正数学定义域不等于修好概率模型 |
| 重复控制/双细胞诊断 | `diagnose_control_guides.py`、`review_control_covariates.py`、`diagnose_doublets.py` | 历史比较和重复尝试；保留原证据，不常规扫参/复读 |
| 前版补充QC | `run_chip01_qc_supplement.py`、`chip01_qc_supplement.py` | 旧z响应、旧参考和原检验的诊断；不作为新版描述性默认 |
| 已交付报告对齐原版 | `run_chip01_report_alignment.py`、`chip01_report_alignment.py`、`plot_chip01_report_alignment.py` | 这版不是废结果；保留精确历史实现。通用部分已抽到核心，避免再带着首芯片硬编码和失准检验一起执行 |
| 历史数值库/小测试 | `perturbation_statistics.py`、`replogle2022_statistics.py`、`smoke_perturbation_statistics.py` | 保留复现和回归用途；数值单测过关不等于真实稀疏数据的p值有效 |
| 文献配置 | `configs/replogle2022.json`、`configs/replogle2022.md` | 历史方法参数，不是新芯片默认政策 |

## 怎么查证/回放

1. 查已完成运行自己的summary、原始配置、输入/输出SHA与固定commit，先区分“只检查产物”和“重新计算”。
2. **精确重算历史运行**应在独立worktree检出该运行的原commit，使用原已验输入和一个全新输出；旧运行目录保护不动。当前归档只是代码的原字节保存，并不是自动移植全部旧路径/运行环境。
3. 只查看本分支归档CLI帮助，可显式使用：

```bash
python scripts/replay_archived_analysis.py run_chip01_report_alignment.py -- --help
```

适配器只补齐归档及活动工具的导入路径，不会自动恢复旧输入/参数、下载数据或启动任何队列；执行计算仍须明确给参数和正确commit。被归档的配置不会自动回填到原`configs/`位置。

`tests/conftest.py`仅在回归测试时加入本目录；当前核心入口不会导入这里。原数值/历史行为测试保留。活动`analyze_expression_structure.py`仍供历史脚本调用，完整历史版本以原commit为准。

## 没有被误归档的东西

原始交付、H5AD原始计数、已验主报告和图册、H/server数据归档未动。表达结构和适当的双细胞评估仍是可选QC；模型研究、通用预处理、迁移和传输验收独立保留，不因为本任务只整理实验QC而删除。
