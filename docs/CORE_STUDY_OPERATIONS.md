# 多芯片核心分析：运行与查证

本文规定如何组织一批分析，不另建科学方法。科学步骤与解释边界见[核心分析路线](CORE_ANALYSIS_ROADMAP.md)；已有单芯片入口仍是`run_core_analysis.py prepare/run/verify`。批次准备入口为[`manage_core_study.py`](../scripts/manage_core_study.py)，只建记录、检查条件、生成配置，**绝不自动启动Cell Ranger或科学分析**。

## 交付与顺序

每个新芯片各有独立计数、完整细胞分母、核心分析和嵌图报告；一份总体报告补充芯片间质量、NTC背景及RNA效应一致性。已有芯片结果可复用，但必须先核对实验关系、参考/基因身份、guide目录和计数口径。不同芯片不自动视作独立生物学重复。

```text
每张芯片的RNA + CRISPR FASTQ（含同一捕获的全部补测）
  → Cell Ranger：不输出新BAM
  → 上游结果验收 → prepare → 工程前缀验证 → 正式core → 图文验收
                              各芯片验收齐全
                                → 保留芯片来源的联合分析/总体报告
```

联合不是再比对一次FASTQ，也不是对各芯片矩阵求平均。不能用来源不明的合并矩阵替代独立计数，不能以聚类或sgRNA猜芯片。总体比较由[run_core_comparison.py及其说明](CORE_COMPARISON.md)负责，输入为已验收的独立核心结果；批次准备入口仍只组织/检查条件。真实总体结果必须实际执行和验收后才能称完成。

## 唯一工作布局

输入保留原址，通过绝对路径、完整哈希回执和当前元数据引用；不搬家、不重命名、不为每张复制整批FASTQ。已接收结果、历史报告、公共数据和模型研究不属于可清理工作区。

```text
<project>/lab/deliveries/…                 受保护交付，永不写入
<project>/lab/analyses/<old-chip>/…        已验结果，永不覆盖
<project>/lab/studies/<study-id>_r1/
  study.json                             固定范围、输入、代码和运行条件
  records/
    input_inventory.json                 内容回执 + 当前元数据，不冒充新hash
    preflight.initial.json               已通过项、缺口、时间和测量边界
    workflow.json                        阶段依赖及拟交付范围
    setup_manifest.json                  此版本准备文件的逐文件SHA
  configs/
    <chip>.cellranger.plan.json           草案；未就绪时不给可执行CSV
    cellranger/<chip>.csv                 全部门禁满足后才能生成，不自动执行
  runs/<chip>/<stage>_<attempt>/          每阶段独立新目录，含日志/参数/回执
  joint/<scope>_<attempt>/                联合计算与跨芯片对照，保留来源
  reports/<scope>_<accepted-run>/         已验读者报告、图册、表和凭据
<reserved-native-work>/<study>/<chip>/    Cell Ranger运行中间文件；另行实查/预留
```

- **代码**在专用Git clone/worktree中，执行固定且已推用户origin的提交；不在数据盘里改运行中代码。
- **D记录区**保存轻量私有参数、检查结果、部署/验收凭据；不放大数据或环境。
- **服务器data**保存输入交付和正式结果。**本地H**按验收后的发布清单归档；不把每个中间文件再镜像一次。
- 临时空间和最终存储分开预留。通用厂家空间要求不是某一芯片的实测峰值；串行可复用临时空间，但已完成结果会累计。共享文件系统显示的空闲总量不是私人配额或已验证的高性能工作区。
- 不在home、/tmp或未预留的共享目录启动大量写入。CIFS的inode审计政策不代表Cell Ranger所有文件操作已被验证；运行工作文件系统需要独立的小文件行为探针。

## 冻结输入与配置

[`core_study.example.json`](../configs/core_study.example.json)是格式示例，路径/哈希/主机与数量需换成真实凭据，不能直接开跑。运行条件不齐可以生成**明确标为未执行的准备记录**。软件、参考或参数就绪后建立新的setup版本；不改旧`study.json`或重写旧检查结论。

必需证据：

| 范围 | 固定内容 | 验收含义 |
|---|---|---|
| 交付 | manifest/ledger的SHA、逐FASTQ源/目标完整SHA+MD5、公司MD5、数量/字节 | 已接收内容身份；本次metadata检查不是新全量读取 |
| 实验 | 每芯片RNA/CRISPR/补测对应，独立捕获关系，时间点/细胞系/真实重复 | 实验者的明确说明可作为来源证据；不因缺少机器可读表而反复要求确认，也不把芯片当生物重复 |
| 参考 | 所选参考的序列/注释身份、文件SHA清单、feature reference SHA | 评估自定义基因定义是否保留；不必拘泥某个供应方目录，但不能静默丢掉研究靶点或凭目录名认定等价 |
| 代码/软件 | 固定Git提交、环境/平台、Cell Ranger版本、命令及资源上限 | 老系统需先确认二进制兼容；不升级共享系统来适配 |
| 方法 | QC规则、四类guide目录、归一化/特征策略、匹配细胞数和种子 | 不照搬旧芯片细胞数或特征数，不自动恢复旧显著性方法 |

FASTQ逐项检查源/目标内容回执一致、vendor值、当前大小/mtime/device、R1/R2配对、芯片总数量与字节。已明确CIFS目标仅inode变化时单独审计，完整内容证明保留；其它变化停止并报告恢复成本，不重拷、重hash后才通知。

## 准备与检查

在真实执行主机的干净固定代码中运行；不在Windows解析服务器路径来创建目录。

```bash
python scripts/manage_core_study.py init --config /private/study.setup.r1.json
python scripts/manage_core_study.py check --study /results/STUDY_r1
# 只有软件、完整参考、实验捕获对应和工作空间等门禁都通过才会写CSV：
python scripts/manage_core_study.py render-cellranger --study /results/STUDY_READY_r2
```

`init`会拒绝已有目录、源/结果重叠、错误主机、脏代码或不符提交。`check`不会读取FASTQ正文，不复制/恢复数据，也不安装或启动程序。已有准备记录逐文件检查SHA，不能改旧配置来绕过门禁。`render-cellranger`即使成功也只写配置，绝不是生产任务启动或科学验收。输入完整性问题抛错停留；缺软件/参考/空间等正常准备缺口返回明确blockers。

正式运行另外需要：实际目录上的小探针、资源/总剩余工作复核、每输出一个写者、独立tmux与日志/状态、启动后父链/命令/创建时间与实际进展核对。先小前缀工程测试，再用全新目录正式运行；工程通过不构成生物学结论。失败保留全部现场，禁止自动重启/覆盖。

## 上游保留与验收

新Cell Ranger配置采用`[gene-expression]`下`create-bam,false`，不影响必要比对/计数。原始FASTQ、公司已有BAM不删；以后逐read检查可能需要从FASTQ重新处理，不能声称新无BAM结果具有旧BAM的全部明细。

除原始/过滤矩阵和guide归属外，保留guide UMI/阈值/feature reference、`molecule_info.h5`或相应分子汇总、原QC/web summary、config/命令/版本/日志及验收清单。不要只留两张表后丢掉复核所需信息。不在本轮自动删除运行中间文件；可回收工作目录只在生产写者退出、保留物逐SHA验收并经明确白名单后处理。

上游验收检查：运行确为该芯片、RNA/CRISPR和全部补测覆盖、无新BAM、矩阵/guide细胞对应、完整参考/guide身份、成功状态与无写者、输出逐SHA。随后才把四个合格输入交给核心prepare。core的独立验收仍按原脚本执行，不因批次状态替代文件和图文验收。

## 联合分析门禁

1. 各芯片独立核心结果已验收，包括完整分母和未评价原因；旧芯片的源counts保持原样。
2. 用`chip_id + original_cell_barcode`作跨芯片唯一身份，同时保留原barcode、原source_row、样本和实验元数据；不得把旧行号直接拼成全局行号。
3. 核对共同基因ID及顺序、完整guide角色目录、QC/归一化和用于横向比较的同一特征族；不把每张各选特征后的相关系数当同空间可比。
4. NTC参考先在芯片内构造，NTC guide仍排除自身全部细胞。总体效果必须同时给分芯片支持与贡献，避免某芯片样本量或技术差异主导；不盲目跨芯片池化背景。
5. 不把chip数当生物重复数、不新增p值/hit门禁、不在原始counts上做不可逆批次校正。不适合直接池化时在总体报告中并列/分层汇总，明确适用队列。

## 报告、发表口径与查证

科学报告按问题顺序嵌入真实图，紧随解释结论与限制；技术部署/传输记录另存。每个图能回到表、参数、输入身份和代码提交。完整靶点表保留“覆盖不足/未评价”，标注正控不冒充已验证正控，描述性结果不叫功能hit。

对外说“完成”前必须逐一核对：每芯片报告、总体报告、实际图、精确覆盖的文件SHA、输入保护、无活动写者及本地归档状态。准备文件、目录、测试通过和单张报告都不等于整批完成。

验证命令：

```bash
python -m pytest -q tests/test_core_study.py tests/test_core_analysis.py tests/test_guide_capture_qc.py
python -m pytest -q tests
```

检查器测试覆盖错误主机、源/输出重叠、丢弃来源、漏配对/漏文件、回执差异、vendor不一致、元数据漂移、仅inode审计、受保护输入不变、禁止BAM配置、拒绝重复初始化和准备记录篡改；不以这些工程测试代替真实数据分析。
