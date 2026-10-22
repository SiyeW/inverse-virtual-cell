# Perturb-seq 核心分析路线

这是实验芯片的**当前默认入口**，不是AI训练路线。保留已采用的覆盖、质控、RNA效应和多sgRNA支持分析；旧统计试验放入[历史归档](../archive/20261023_legacy_chip_analysis/README.md)，不再默认串起来运行。

## 一条主线

| 顺序 | 要回答的问题 | 核心输出 |
|---|---|---|
| 0. 输入与细胞QC | 读到什么数据？哪些细胞有可靠guide？ | 全量稀疏原始counts、逐细胞QC、完整分母、输入SHA |
| 1. 文库与控制角色 | NTC、基因间区、标注正控、普通targeting各有多少支持？ | 全文库guide覆盖、零覆盖、四类细胞数 |
| 2. 技术质量 | 四类RNA及guide质量是否偏移？ | RNA UMI/特征数/线粒体比例/guide UMI及保留率 |
| 3. RNA效应背景 | NTC是否稳定，targeting是否偏离？ | 留掉自身guide的RNA log2FC、等n分布/RMS、表达对角线图 |
| 4. 多sgRNA支持 | 同一靶点不同guide是否一致，是否被单条guide主导？ | 全靶点支持表、guide对相关、去最大guide敏感性、顺序图册 |
| 5. 报告与查证 | 图的含义、限制和来源是否能查清？ | 嵌图Markdown主报告/图册、固定参数/抽样、输出SHA、只读验收 |

**默认不做**：原AD/MW显著性命中、据此筛选的功能聚类、为旧弱分组强加机制、DNA编辑率推断、模型训练。已发现原渐近检验在稀疏并列秩NTC中失准；图像看起来像正态、均值靠近对角线或数值实现相符，都不使其重新有效。

## 核心文件

- [`scripts/run_core_analysis.py`](../scripts/run_core_analysis.py)：唯一主命令，`prepare` / `run` / `verify`。
- [`perturbseq_core_io.py`](../scripts/perturbseq_core_io.py)：输入身份、原始CSR计数、guide目录和RNA特征契约。
- [`perturbseq_core_effects.py`](../scripts/perturbseq_core_effects.py)：从已接受流程抽出的描述性效应与靶点支持算法。
- [`perturbseq_core_plots.py`](../scripts/perturbseq_core_plots.py)、[`perturbseq_core_report.py`](../scripts/perturbseq_core_report.py)：静态图、顺序图册、自动图文报告及只读验收。
- [`requirements-core-analysis.txt`](../requirements-core-analysis.txt)：CPU分析依赖，不需要模型或GPU环境。每次运行另记录实际Python/包版本；不是声称已锁定所有平台的二进制环境。

## 输入约定

从已验收的完整Cell Ranger矩阵、guide calls和feature reference开始，或直接使用**已经包含全部细胞分母**的RNA-QC H5AD。只含单guide的H5AD不能冒充完整分母。

H5AD的`X`必须为原始非负整数、RNA-only CSR；不能传入归一化或log矩阵。RNA feature ID和cell barcode唯一。`obs`保留`guide_id`、`perturbed_gene`、`is_control`、`guide_status`、`qc_pass_stage1`、`main_analysis_stage1`及RNA/guide技术指标；三个策略字段是真布尔值，不接受字符串自动转真。

`main_analysis_stage1 = qc_pass_stage1 AND single_guide`，核心流程核对而不重选。批次、细胞系、时间点、真实生物学重复等元数据保留；芯片号不自动等于生物学重复。

guide目录是覆盖**全部设计guide（包括零覆盖）**的CSV：

```csv
guide,target_gene_name,source_control_type
NTC_001,Non-Targeting,non_targeting_control
IG_001,Non-Targeting,intergenic_control
POS_001,GENE_A,positive_control
KO_001,GENE_B,targeting
```

这只是格式示例，不是可直接使用的生物学注释。必须用实验设计/已核对目录映射真实guide ID和靶点；如来自另一个测序参考，应先按已记录序列方向完成身份对应，不能猜名称、逐次尝试方向挑最高匹配。当前CSV中的`target_gene_name`须匹配源H5AD的粗粒度靶点，阴性角色对应`control_target`且源细胞记为`ctrl`。一个靶点不允许冲突角色。标注正控不自动变为本RNA实验的已验证正控。

## 怎样运行

把配置复制到仓库外，填写真实路径、已验SHA、维度/大小、不同的输出目录，以及已推送的固定Git提交。正式运行要求干净工作树；禁止覆盖已有运行目录，失败后保留现场，不自动重试。

### 0. 若已有合格的全细胞H5AD，跳过这步

使用[`core_prepare.example.json`](../configs/core_prepare.example.json)，四个输入都要填SHA。路径用绝对路径，输出父目录须已经存在。

```bash
python scripts/run_core_analysis.py prepare --config /records/chipXX/prepare.json
python scripts/run_core_analysis.py verify --output /results/chipXX/NEW_BASELINE
```

调用既有`analyze_guide_capture_qc.py`，全量分块测量，保留原始计数和所有guide分配；记录RNA阈值并逐行回读计数总和。默认低UMI/低特征采用log1p下侧5 MAD并设500 UMI/200特征下限；线粒体上限为max(8%, 中位数+3 MAD)。高计数只标记，不自动删细胞。它不是FASTQ比对或最终双细胞排除。

输出文件名随`chip_id`变化，不再固定写`chip01`。guide catalogue的四类细分仍须提供，不能仅凭Cell Ranger粗标签发明。

### 1. 小前缀验证，再用不同目录正式分析

使用[`core_analysis.example.json`](../configs/core_analysis.example.json)。`source`的大小/维度/SHA、`catalogue`的SHA均须填写；`output`直接位于`allowed_output_parent`下且必须不存在。

```bash
python scripts/run_core_analysis.py run --config /records/chipXX/smoke.json --smoke
python scripts/run_core_analysis.py verify --output /results/chipXX/SMOKE_ONLY
# 验证后，用另一个output和run_id：
python scripts/run_core_analysis.py run --config /records/chipXX/formal.json
python scripts/run_core_analysis.py verify --output /results/chipXX/FORMAL_RUN
```

`--smoke`只读前`smoke_prefix_rows`行CSR和必要元数据，不重读全文件hash；输出明确标为非科学结果，阈值仍按配置。前缀没有足够NTC或合格guide时会停，不自动扩大前缀、降阈值或重跑。

正式运行完整验证一次源SHA；全程源只读，结束检查身份和小输入SHA。`verify`只查结果，不重新打开原始矩阵、重算效应或触发恢复；若用于`prepare`结果，会读输出H5AD做SHA，耗时按输出大小计。验收核对哈希/精确文件覆盖、图和链接；核心结果另核对抽样、效应公式、全guide留出和覆盖分母。这不是第二次独立原始数据重算。

### 2. 阅读与复核

正式`run`目录的最小交付：

```text
分析报告.md                   # 按科学顺序嵌入6张核心总览图，并讲解边界
逐靶点支持图册.md             # 所有可评价靶点和标注正控，按全靶点顺序编号
figures/                      # 含targets/，相对链接，可整目录携带
tables/                       # 全guide/全target、未评价原因、细胞/特征身份
matched_cell_selections.json   # 固定种子与每次抽样的原始行号
config.resolved.json
summary.json                  # 队列、参数、版本与代码提交
output_manifest.json          # 每个静态产物的SHA256
status.json                   # complete才可接受；owner.lock完成后释放
```

所有目标都留在表里；不能出图的目标不假装“无效应”。图册不按显著性挑图、不根据结果调阈值。

## 固定的方法和需要随芯片变化的参数

| 保持不变的规则 | 随芯片确认/记录的量 |
|---|---|
| 原始counts不改，完整分母/零覆盖保留 | 输入SHA/大小/维度、实验身份及全量guide目录 |
| 全部真正NTC为参考，NTC guide排除自身全部细胞 | 实际NTC数量；至少两条被观测NTC且留出后至少两细胞 |
| RNA log2FC为均值比，伪计数固定1 | 归一化尺度：新芯片默认NTC细胞UMI中位数；回放用原固定值 |
| 相同细胞数比较，guide不是生物学重复 | 默认每guide>=10细胞、抽10细胞×8次；显式改动另记run_id |
| 直接靶RNA只在全参考中唯一时屏蔽 | 标准基因/自定义ORF身份；不猜宿主基因 |
| label-blind固定特征族，不按p或强效应选RNA | 新芯片`pooled_raw_mean > 0.5`；或已验`frozen_ids`CSV |
| 缺覆盖和未评价明确保留，不自动判hit | 每靶点guide覆盖、多guide相关及正控预期证据 |

`pooled_raw_mean`在该芯片全体RNA-QC单guide细胞中选择RNA，没有guide/靶点显著性选择。它是新芯片的明确参数化描述策略，**不是首芯片旧5,531特征族的同义替代**。如需严格回放，使用`features: {"mode":"frozen_ids","path":"...csv","sha256":"..."}`，CSV列为`gene_id`、按源顺序且无重复；归一化也填原值。输出自动保存实际特征族。

## 首芯片复现与其它芯片迁移

首芯片已交付13图主报告和113图图册保持不动。其历史整套计算固定在`c82bcf418731122de2f3244a4bbb1baed6225d7d`，包含参考敏感性和原检验失准诊断；当前默认核心不是声称重跑了这些历史诊断。精确重现历史版本：在独立Git worktree检出该提交，使用原哈希输入/原完整配置和全新输出，不从当前归档目录盲启动旧任务。

当前核心的首芯片描述性回放采用原15,259 UMI、原5,531特征CSV、同一细胞策略、NTC全guide留出、原种子和抽样数；不需要重新生成AD核心对照或校准表。实验私有配置/guide目录/特征清单保存在项目D记录区，不提交Git。

迁移新芯片只更换输入身份、完整guide角色映射、实验信息和新输出位置，确认归一化/特征策略后小测试。**不要复制首芯片的细胞数、5,531基因数、旧2,505核心细胞或75条核心guide作为新芯片常量**。本次整理不等于其它芯片已经计算完成。

## 按需工具与历史路线的区别

- **表达结构/双细胞质控仍有用**：`analyze_expression_structure.py`、`audit_scrublet_workflow.py`保留。当前是读取既有baseline身份/来源记录的独立诊断工具，不是这次默认核心自动生成的UMAP或删细胞规则；新芯片需先确认输入布局和阈值诊断，不能复用首芯片图或不可靠自动截点。
- **旧AD/MW校准反例要留，但不常规重复**：归档保留完整代码和测试。未来引入推断统计时，先通过适配稀疏度/样本量/并列秩的NTC校准，再谈命中与功能分组；这不是本轮新增统计方法任务。
- **旧能量距离保留作历史结果**：不因旧DE失准就宣布它也必错；但旧参考/队列结果不能当作新芯片自动重算。单独提出该统计任务时，使用其明确表示/零假设/完整参考和成本门禁。
- **模型研究与接收工具分开**：KOLF/Adamson、contrastiveVI、GenePT及传输验收脚本不是本芯片QC的一环，也不因“不在核心”就删除。它们原位独立保留。

## 验证与维护

```bash
python -m pytest -q tests/test_core_analysis.py tests/test_guide_capture_qc.py
python -m pytest -q tests
```

覆盖：已接受算法数值一致、新芯片命名、固定特征身份、四类及缺少可选控制的完整报告、零/单guide未评价、原始计数不变、输入错配/旧输出拒绝、产物篡改拒绝、归档原字节检查。真实原始数据没有因为整理代码而重新计算；后续正式迁移仍需执行对应芯片的输入和小样本门禁。
