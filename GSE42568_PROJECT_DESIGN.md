# GSE42568 乳腺癌与正常乳腺组织表达谱项目设计

## 1. 架构结论与适用边界

本文件定义一个可直接进入实现阶段、但本轮不包含任何分析代码或生物学结果的项目规范。目标数据集为 GEO GSE42568，核心比较为 104 个乳腺癌组织样本与 17 个正常乳腺组织样本在 Affymetrix Human Genome U133 Plus 2.0 Array（GPL570）上的表达差异。GEO 存档明确给出 121 个样本、平台、样本级来源、处理方法及 log2 GC-RMA 信号；原始 CEL 文件亦可获得。[^1][^2][^3]

架构决策如下。

1. **标准版本使用 GEO 发布的 processed series matrix，不重新从 CEL 起步。**该矩阵已经是 log2 GC-RMA 表达值；重复取对数、再次分位数标准化或把它当作 RNA-seq count 使用，都会改变已发布数据的含义。
2. **主要估计对象是本队列中的癌组织—正常组织平均表达差异。**这是非配对、横断面的 bulk-tissue association，不是肿瘤发生的因果效应，也不是诊断、预后或治疗反应模型。
3. **主要差异表达模型纳入 processing-run proxy 与 tissue group。**样本名/CEL 文件名显示 13 个处理日期层级；其中 6 个层级同时含癌和正常样本，因此组别效应原则上可估计，但存在明显的部分混杂，必须用 group-only 与 mixed-run 子集进行敏感性验证。
4. **正常组只能称为 normal breast tissue。**GEO 没有说明它们是健康供者、癌旁组织、减乳组织或与肿瘤配对的组织；正常样本年龄全部缺失，也没有 pairing identifier。因此正文不得使用 “healthy controls”“matched normal” 或 “adjacent normal”。
5. **基因级主要结果采用预先规定、与结局标签无关的代表探针规则，并保留探针级与另一种聚合方法作敏感性分析。**不得按最小 P 值、最大 |t| 或最大组间差异挑探针。
6. **多重检验以 BH FDR 0.05 为统计发现标准；TREAT 的 1.2 倍最小效应阈值用于优先级清单。**完整结果表永远保留全部被检验基因，不能只交付“显著基因”。
7. **GO ORA 与 Hallmark GSEA 可进入标准包；KEGG 必须经过商业许可闸门。**在未确认项目用途对应的 KEGG 商业许可前，不调用 KEGG REST、不缓存/再分发 KEGG 内容，也不把 KEGG 图作为默认客户产物。可用 Reactome 作为许可清晰的路径数据库替代。[^16][^17]
8. **任何样本剔除都必须基于预先定义的技术证据，而不是 PCA 图“看起来不顺眼”。**仅凭肿瘤组更分散、聚类不完全分开或某样本影响显著，均不足以删除样本。

### 1.1 本轮架构核验快照

为避免设计建立在网页摘要的想象上，本轮对 GEO 官方 matrix 与 121 条样本记录做了只读核验。以下数值是 **2026-09-11 获取版本的输入审计快照**，不是差异表达结果：

| 项目 | 核验值 | 设计含义 |
|---|---:|---|
| 样本数 | 121 | 104 cancer，17 normal |
| 探针集数 | 54,675 | 与每条 GSM 的 row count 一致 |
| 缺失/非有限值 | 0 | 可进入下游 QC，但不代表样本质量合格 |
| 重复 GSM / 重复探针 ID | 0 / 0 | 主键完整 |
| 全矩阵范围 | 2.3128–16.0469 | 与 log2 微阵列信号一致 |
| 样本中位数范围 | 2.7262–3.2495 | 后续仍须画分布图并做鲁棒异常审计 |
| 样本 IQR 范围 | 2.5565–3.4813 | 后续仍须逐样本检查 |
| 常数探针集 | 96 | 不进入建模；进入过滤审计表 |
| `AFFX-` 控制探针集 | 62 | 用于技术审计，不进入基因推断 |
| series matrix SHA-256 | `43b93f72f1d81838dcd2a5b6ccc3bc1875e46230dad0d7346236773a99e32b76` | 仅标识本次取得的压缩文件版本 |
| series matrix 压缩大小 | 22,565,943 bytes | 与当次服务器响应一致 |

矩阵上游文件可能重新压缩或更新，所以哈希不是永久常量。正式运行必须同时保存 accession、直接 URL、HTTP Last-Modified、下载时间、字节数和 SHA-256；如果任一项与锁定清单不同，就视为一个新的数据版本，而不是悄悄覆盖旧结果。

## 2. 科学问题、估计对象与不作出的结论

### 2.1 主要问题

在 GSE42568 的 121 份 bulk breast-tissue expression profiles 中，癌组织相对于正常乳腺组织，哪些基因的平均 log2 表达存在稳定差异？哪些 GO biological processes 与经许可的 pathway/gene-set programs 在该对比中富集？

主要对比固定为：

> cancer − normal

因此：

- 正的 log2 fold change 表示癌组织中表达更高；
- 负的 log2 fold change 表示癌组织中表达更低；
- 线性模型效应是调整 processing-run proxy 后的组均值差；
- 结果描述限于该存档队列和该技术平台。

### 2.2 统计估计对象

主要 estimand 是所有通过预定质量闸门的样本中，控制离散 processing-run proxy 后，cancer 与 normal 的 gene-level log2 expression mean difference。它不等价于：

- 同一患者肿瘤与癌旁组织的配对差异；
- 年龄匹配后的差异；
- 纯肿瘤细胞与纯正常上皮细胞的差异；
- 癌症导致的因果改变；
- 肿瘤亚型、分级、ER、淋巴结、复发或生存的效应；
- 可在临床个体上直接使用的分类器或 biomarker panel。

### 2.3 明确不纳入主要模型的变量

年龄、肿瘤大小、grade、ER、node、术后治疗、复发和总生存信息不进入 cancer-normal 主要设计：年龄在 17 个正常样本中全部缺失；其余变量只对肿瘤有定义或发生在取材以后。把这些变量直接加进主要模型会删除所有正常样本、产生不可估计设计，或错误控制肿瘤后的变量。

## 3. 实验、数据、样本与平台的可核查事实

### 3.1 GEO 与原论文的一致信息

GSE42568 的 GEO 标题为 “Breast Cancer Gene Expression Analysis”。GEO 描述 104 个乳腺癌活检样本，取材发生在 tamoxifen 或 chemotherapy 之前，另有 17 个正常乳腺组织样本；所有样本均使用 GPL570。原论文将这 104 个原发乳腺癌样本描述为 1993–1997 年在 Dublin 的 St Vincent’s University Hospital 收集，组织在切除后 30 分钟内 snap-freeze，并以与研究中其它阵列相同的 GCRMA 流程处理。[^1][^4][^5]

原论文还整合了多个公开队列并使用 ComBat 处理跨研究批次；这 **不能** 推导出 GEO 中这个单队列 deposited matrix 已经经过 ComBat，也不能成为本项目再次对单队列盲目 ComBat 的理由。GEO 样本记录对本 matrix 的明确表述是 GC robust multi-array average，数值为 “Log2 GC-RMA signal intensity”。

### 3.2 样本数量与元数据完整度

| 字段 | Cancer，n=104 | Normal，n=17 | 在主要模型中的角色 |
|---|---:|---:|---|
| tissue group / source | 104 完整 | 17 完整 | 主要解释变量 |
| GPL570 / row count | 全部一致 | 全部一致 | 输入一致性闸门 |
| age | 104 完整 | 0 完整 | 不调整；只做肿瘤描述 |
| tumor size | 104 完整 | 不适用 | 不调整 |
| grade | 11 G1 / 40 G2 / 53 G3 | 不适用 | 不调整 |
| ER | 67 positive / 34 negative / 3 missing | 不适用 | 不调整 |
| node | 45 negative / 59 positive | 不适用 | 不调整 |
| relapse-free status | 56 censored / 48 event | 不适用 | 不属于当前问题 |
| overall-survival status | 69 censored / 35 death | 不适用 | 不属于当前问题 |
| pairing identifier | 未提供 | 未提供 | 设计为非配对 |
| normal provenance | 不适用 | 未充分说明 | 核心局限 |

样本级年龄核算为 31.06–89.93 岁，平均 58.83 岁；27 个肿瘤样本小于 50 岁、77 个不小于 50 岁。GEO series-level 文字写成 20 与 77，二者合计也不是 104。正式元数据表应保留这一差异，分析描述以可逐条追溯的 sample-level 字段为准，并在数据问题日志中记录，而不是静默“修正”网页。

GEO 还给出肿瘤大小 0.6–8 cm、样本级平均 2.7875 cm，以及 series-level 的组织学汇总；但组织学没有逐样本映射，不能将汇总计数反推给具体 GSM。

### 3.3 平台与可解释层级

GPL570 是 Affymetrix Human Genome U133 Plus 2.0 Array，属于 in situ oligonucleotide commercial array。官方平台记录称其覆盖完整 U133 probe-set collection，并增加约 9,921 个 probe sets；GSE42568 实际矩阵有 54,675 行。[^3]

平台行是 **probe set**，不是天然的一行一个当前基因。一个基因可能对应多个探针集，探针集可能跨转录本、交叉杂交或映射到多个基因，旧 probe annotation 也可能与当前基因命名不一致。因此：

- 探针级信号可以建模，但不能直接把 54,675 行称为 54,675 个基因；
- symbol 只作展示标签，稳定键使用 Entrez Gene ID；
- 任何 pathway universe 必须从真正进入 gene-level 检验的唯一基因构造；
- top-gene 结论必须同时呈现代表探针、构成探针数和方向一致性。

### 3.4 processing-run proxy 的证据

样本标题和原始 CEL 文件名中的日期后缀形成以下分层。3 个无日期后缀的 `T` 样本经有限 CEL header text audit 均显示 2004-11-26；一个标记 2005-01-25 的样本 header 显示次日扫描。因此该字段只能命名为 `processing_run_proxy`，不能声称是取材日或精确 scan date。

| processing-run proxy | Normal | Cancer | 是否可在层内比较组别 |
|---|---:|---:|---|
| 2004-11-26（3 个无后缀 T 样本） | 0 | 3 | 否 |
| 2004-12-03 | 0 | 8 | 否 |
| 2004-12-07 | 0 | 8 | 否 |
| 2004-12-08 | 0 | 8 | 否 |
| 2004-12-09 | 0 | 11 | 否 |
| 2004-12-10 | 0 | 10 | 否 |
| 2004-12-14 | 3 | 12 | 是 |
| 2004-12-15 | 3 | 10 | 是 |
| 2004-12-20 | 0 | 11 | 否 |
| 2004-12-21 | 1 | 10 | 是 |
| 2004-12-22 | 1 | 9 | 是 |
| 2005-01-25 | 8 | 3 | 是 |
| 2005-01-26 | 1 | 1 | 是 |
| **合计** | **17** | **104** | 6 个 mixed levels |

所有 normal 都位于 mixed levels，而 59 个 cancer 位于 cancer-only levels。组别与 run 不是完全混杂，所以主要效应在代数上可能估计；但估计在很大程度上由 62 个 mixed-run 样本（17 normal、45 cancer）提供信息。这个结构必须进入报告正文，不能只在 supplement 里出现。

样本名前缀 `N/P` 只出现在 normal，`S/T` 只出现在 cancer，与 tissue group 完全重合。它不能作为第二个“批次”协变量，否则会让组别不可估计；前缀也不能被解释成患者类型或取材来源，除非取得新的原始记录。

## 4. Processed GEO 数据是否适合直接分析

### 4.1 可用性判断

在本项目的目标——建立一个公开、可复现、范围受控的 cancer-normal case study——下，processed matrix **有条件适合**作为标准版输入，理由是：

- GEO 逐样本明确标注 GC-RMA，并明确数值是 log2 signal；
- 121 个样本全部为 GPL570，行数一致；
- 官方 matrix 维度、ID、有限值和数值范围通过本轮基础审计；
- 原论文报告使用 GCRMA 与 quantile normalization，和 GEO 存档方向一致；
- 主要问题是均值差异，不要求重新定义 probe-level background model。GCRMA 本身包含针对 GC affinity 的背景校正与多阵列汇总流程。[^6]

### 4.2 标准版不得做的“二次处理”

- 不再次 `log2(x + c)`；输入已经是 log2。
- 不将数值取整，不用 DESeq2/edgeR 等 count model。
- 不默认再次 quantile-normalize。GCRMA 在 probe level 正规化后再汇总，最终 probe-set 分布不必逐分位点完全相同；“分布不完全重合”不是二次正规化的充分理由。
- 不对整个矩阵预先 ComBat 再送入差异表达。主要模型用设计矩阵表示 run，既保留不确定性，又避免把感兴趣组别作为未知批次的一部分移除。
- 不在组间分别做 normalization，也不使用 cancer/normal 标签挑高变探针或决定过滤。

### 4.3 processed-only 能与不能验证的内容

标准版可以验证：文件完整性、样本和平台一致性、signal distributions、样本相关性、PCA、run/group 结构、可估计性、gene annotation、差异表达与 enrichment 稳定性。

标准版不能诚实声称验证：

- probe-level model 的 NUSE/RLE；
- RNA degradation slope；
- array image spatial artifacts、bubble、scratch 或 scanner saturation；
- CEL-level background correction 是否最优；
- alternative CDF/Brainarray 重汇总后的稳健性。

如果 processed QC 出现无法解释的严重异常，应触发 raw-CEL extension，而不是在 processed matrix 上叠加无法追溯的“修复”。

## 5. 数据冻结、命名与可追溯性

### 5.1 输入清单

标准运行必须锁定并生成 `data/metadata/source_manifest.tsv`，至少包含：

| 字段 | 含义 |
|---|---|
| `accession` | GSE42568 / GPL570 / 每个 GSM |
| `asset_role` | matrix、series metadata、sample metadata、platform annotation |
| `source_url` | 官方直接 URL，不是搜索结果页 |
| `retrieved_at_utc` | ISO 8601 UTC 时间 |
| `http_last_modified` | 上游返回值，若无则明确 NA |
| `content_length_bytes` | 下载字节数 |
| `sha256` | 原始取得文件的 SHA-256 |
| `compression_test` | gzip 完整性结果 |
| `parser_version` | 解析工具及版本 |
| `license_or_terms_url` | 数据/注释/基因集条款 |
| `notes` | 上游异常、重定向或版本说明 |

GEO 提供多种程序化下载路径；项目只使用 NCBI/GEO 官方 HTTPS/FTP 端点，并对下载后的压缩文件做 checksum 与解压完整性验证。[^2]

### 5.2 样本清单是唯一事实表

`data/metadata/sample_manifest.tsv` 一行一个 GSM，必须含：

- `gsm_id`、原始 title、原始 source name；
- 规范化后的 `tissue_group`，只允许 `cancer` / `normal`；
- `platform_id`、row count；
- `processing_run_proxy` 与推导来源；
- age、tumor size、grade、ER、node、RFS、OS 原始值与规范化值；
- 每个字段的 missingness reason：`not_reported`、`not_applicable`、`parse_failed` 不可混成空字符串；
- 是否进入 primary、每个 sensitivity analysis，以及排除原因；
- 原始文本字段的 provenance，确保标准化不覆盖原始记录。

任何 figure/table 的样本顺序都由这张 manifest 控制。表达矩阵列名必须与 manifest 的 GSM 集合一一相等；只允许显式重排，不能依赖文件当前列顺序。

### 5.3 数据版本变化策略

若官方文件哈希改变：

1. 保留旧 manifest 与旧报告；
2. 新建 dataset version，而不是覆盖；
3. 比较 gzip 解压后的文本哈希、维度、header、数值和 metadata；
4. 若只是压缩容器变化但解压内容相同，记录为 transport-only change；
5. 若内容变化，重新运行全部 pipeline，并在 changelog 中列出对结果的影响。

## 6. 分析总体设计

### 6.1 分析总体与样本包含规则

Primary population 为 121 个 deposited samples，条件是每个样本通过文件、身份、平台和最低技术质量闸门。组别不做下采样；104:17 的不平衡本身不会使 limma 失效，而随机丢弃 tumor 会损失精度并引入不必要的随机性。

样本排除只允许三类原因：

1. 文件或 sample ID 无法与官方记录匹配；
2. 数值损坏、平台/维度不一致等不可修复的技术失败；
3. 至少两个相互独立的异常指标同时触发，且有分布、相关性、raw header 或原始记录的佐证，由预先指定 reviewer 签字确认。

病理异质性、极端但可信的生物表达、对结果影响较大，均不是自动排除理由。所有排除都必须保留 primary-with-all 与 exclusion sensitivity 的并列结果。

### 6.2 主要线性模型

对每个 gene-level expression feature 拟合：

`expression ~ processing_run_proxy + tissue_group`

并计算 `cancer - normal` contrast。run 作为固定效应；不使用 patient random effect，因为没有配对或重复测量标识。假设 13 个 run levels 均保留且没有其它缺失，设计有截距、12 个 run 系数和 1 个 group 系数，预期 residual degrees of freedom 为 107。

运行前必须：

- 对设计矩阵做 QR rank check，要求 rank 等于列数；
- 输出 alias/estimability 报告；
- 报告 group coefficient 的标准误、设计矩阵 condition diagnostics 和每个 run 的组别构成；
- 若 group 不可估计，立即停止，不得通过删除列或把 run 数值化来强行运行。

### 6.3 limma 经验贝叶斯策略

标准版使用 limma linear model 与 moderated statistics。主要推断使用 `trend=TRUE` 与 `robust=TRUE`：前者允许 residual variance 随平均表达水平平滑变化，后者降低少量异常方差基因对经验贝叶斯估计的影响。是否启用 trend 必须由 mean–variance diagnostic 支持，同时保留无 trend 的结果作技术敏感性比较。limma 文档也明确区分在设计矩阵中调整批次与只为绘图使用的 `removeBatchEffect`。[^7]

不将 visualization-only batch-adjusted matrix 用于正式 DE。若 cancer 与 normal 的 variance trend 明显不同，则增加预先声明的 heteroskedastic sensitivity（例如 limma 的 group-aware vooma/voomaLmFit 方案）；它不替代主要模型，也不能在看到期望结果后才决定启用。

### 6.4 四个预先定义的分析集

| ID | 样本/模型 | 目的 | 主要判读 |
|---|---|---|---|
| P1 | 全部合格样本；run + group | 主要结果 | 正式 effect、CI、FDR、TREAT |
| S1 | 全部合格样本；group only | 检查 run 调整影响 | 与 P1 比较 effect/rank，不作为替代主结论 |
| S2 | 6 个 mixed runs，n=62；run + group | 让组别比较完全来自层内重叠 | 看方向与效应一致性；不要求相同显著基因数 |
| S3 | 排除经签字确认的技术失败；run + group | 检查异常样本影响 | 若无样本排除，则明确标记 not applicable |
| S4 | P1 的 `trend=FALSE` 与必要时 group-aware variance model | 检查均值—方差建模选择 | 比较 effect、SE 和 top-rank 稳定性 |

S2 样本量较小，显著性下降是预期的；稳定性判断侧重 effect size、direction 和 rank correlation，而不是要求每个 P1 FDR hit 在 S2 仍达到 0.05。

## 7. QC 方案与闸门

### 7.1 输入与结构 QC

以下为硬性 assert，任何一项失败都停止：

- gzip/文本可完整读取，SHA-256 与 manifest 一致；
- matrix 恰有 121 个唯一 GSM 列与 54,675 个唯一 probe-set 行；若上游版本变化导致不同维度，先进入版本审计，不把旧数字硬编码为永恒事实；
- matrix sample IDs 与 sample manifest 集合完全一致；
- 121 个 GSM 全部为 GPL570；
- 所有表达项可转换为有限 numeric，无 `NA`、`Inf`、`-Inf` 或非数值 token；
- `tissue_group` 计数为 104/17，且映射来源可回溯；
- processed scale 与 GEO 的 log2 声明一致；如出现百级/千级值或其它强烈反证，停止确认，而不是自动转换。

### 7.2 样本分布 QC

必须生成每样本 density、boxplot、median、IQR、MAD、1%/99% quantile 与 dynamic range。每个标量指标用 median/MAD 计算 robust z；`|robust z| > 3` 只标记为 warning，不单独触发删除。

判读原则：

- GC-RMA 后箱线图不需要像素级重合；
- group-wide shift 可能是真实组织组成差异，不是自动批次；
- 某一 run 的所有样本同步偏移提示技术效应；
- 单个样本同时在位置、离散度和相关性上异常，才需要深入核对 raw CEL/记录。

### 7.3 PCA QC

PCA 默认在通过 annotation/filter 的 gene-level matrix 上执行：按全部样本计算 MAD，选择 top 5,000 variable genes；对每个 gene 做中心化，但不按 gene 标准差缩放。这样避免低表达噪声因为 unit-variance scaling 获得与稳定高信号相同权重。

必须输出：

- PC1–PC2、PC1–PC3，分别按 tissue group 和 run 着色；
- 每个 PC 的 explained variance；
- top 1,000、top 5,000 与全部 eligible genes 的 PCA sensitivity；
- group、run 与前 10 PCs 的关联/partial variance summary；
- candidate outliers 的 robust PCA distance。

对前若干 PCs（累计解释至少 50% variance，最多 10 PCs）计算 robust distance；超过基于稳健协方差的 0.999 reference contour 只记为一个异常信号。PCA 分开不是分析成功条件，PCA 不分开也不是失败条件。

### 7.4 样本相关性与聚类 QC

使用全部 eligible gene-level rows 计算 Pearson sample correlation；Spearman 作为敏感性。聚类距离固定为 `1 - correlation`，average linkage。必须交付：

- 121×121 correlation heatmap，按 group/run 注释；
- 每个样本 median correlation 与 network connectivity；
- 低相关样本的 robust z；
- blind-to-group 的样本树，以及带 metadata annotation 的解释版。

最低相关性不设置一个跨平台通用的任意绝对阈值；异常判定使用队列内 robust lower-tail、PCA 和 distribution 的联合证据。若某样本仅因肿瘤生物学异质性而相关性较低，则保留。

### 7.5 批次与可估计性 QC

必须在任何 DE 之前生成：

- `tissue_group × processing_run_proxy` 交叉表；
- run/group 在 PCA 上的分布；
- 每个 PC 的 group 与 run 解释度；
- full-rank 和 contrast estimability 检查；
- P1、S1、S2 的 group effect rank/effect concordance。

若 run 与 group 实际完全混杂、group coefficient 不可估计，项目必须停在“描述性 QC”，不能发布 adjusted DE。若可估计但结果严重依赖 cancer-only runs，则降级结论，并以 S2 作为必要的证据边界。

### 7.6 样本异常升级规则

候选样本只有满足下列四类中的至少两类，才进入人工 exclusion review：

1. distribution metric robust z 超过 3；
2. median correlation 或 connectivity lower-tail robust z 低于 -3；
3. robust PCA distance 超过 0.999 contour；
4. 原始记录/CEL header/CEL-level extension 提供独立技术异常证据。

最终删除还需：明确理由、reviewer、时间、受影响分析、with/without 对比和决策日志。仅凭组标签、聚类位置或影响力不满足条件。

## 8. 探针注释、重复探针与 gene-level 规则

### 8.1 注释版本冻结

主键使用 Affymetrix probe-set ID，映射目标使用 Entrez Gene ID。实现时锁定 Bioconductor release 与完整 session info；在本设计日期可锁定 Bioconductor 3.23、`hgu133plus2.db` 3.13.0、相应 `AnnotationDbi`/`org.Hs.eg.db`/`GO.db` 版本。`hgu133plus2.db` 页面虽然持续发布包版本，其底层 Entrez 注释来源日期较早，因此 top findings 还需要与当次 NCBI Gene、HGNC 或 Ensembl 当前记录交叉核验，并把变化记录下来。[^8][^9]

不得只存 gene symbol：symbol 会更名、复用或丢失。输出至少保留 `probe_id`、`entrez_id`、`symbol_at_run`、`gene_name_at_run`、annotation package/version/source date。

### 8.2 映射审计分类

所有 54,675 个 probe sets 必须进入 mapping audit，并互斥归类：

- `control_affx`：`AFFX-` controls；
- `constant`：跨样本没有 variance；
- `unmapped`：没有当前 Entrez mapping；
- `one_to_one`：唯一 Entrez；
- `one_to_many`：一个 probe set 映射多个 Entrez；
- `many_probes_per_gene`：唯一映射，但一个 Entrez 有多个 probes；
- `cross_hybridizing`：`_x_at`；
- `unknown_suffix`：不符合已审查 suffix 规则。

`AnnotationDbi::select` 的 one-to-many 展开是数据事实，不是让程序“取第一行”的许可。每类必须报告数量、百分比及前若干实例。[^8]

### 8.3 Affymetrix suffix 处理

根据 Affymetrix 官方 probe-set 命名说明，`_at` 通常有较高唯一性；`_a_at` 可识别同一基因的 alternative transcripts；`_s_at` 可能共享于多个 transcripts；`_x_at` 的交叉杂交行为更不可预测。[^10]

主要 gene-level 候选限定为：

1. 非 `AFFX-`；
2. 非 constant；
3. 唯一映射到一个 Entrez；
4. 非 `_x_at`；
5. suffix 已知并可分级。

`_s_at` 不因 suffix 单独删除，但排在 `_at` 与 `_a_at` 后；one-to-many probes 不进入主要 gene-level test。它们仍保留在探针级 supplement，避免把不确定性藏起来。

### 8.4 一个基因多个探针的主要规则

对每个 Entrez，在完全不使用 tissue label、P value、fold change 或 t statistic 的前提下选择一个代表探针：

1. specificity tier 优先顺序：unique `_at` > `_a_at` > `_s_at`；
2. 在最高可用 tier 内，选择全部 121 样本 pooled median expression 最高的 probe set；
3. 完全相同时按 probe-set ID 字典序决定；
4. 把所有候选、tier、median 与最终选择写入 `probe_selection_audit.tsv`。

这样做避免把同一基因多次计入 FDR 和 enrichment，也避免按最显著探针 cherry-pick。代价是可能只代表某个转录本，所以必须同时做两层敏感性分析：

- **Probe-level layer**：所有合格、唯一 Entrez 的 probe sets 分别建模，结果只作技术 supplement；
- **Median-collapse layer**：每个样本内对同一 Entrez 的全部合格非 `_x_at` probes 取 median，再按同一模型分析。

headline gene 必须满足：代表探针与 median-collapse 方向一致；若有多个 constituent probes，报告同向比例及是否有强烈反向探针。方向不一致的基因可以保留在完整表，但不得作为无保留的单基因主结论。

### 8.5 注释质量闸门

下列情况触发停止并调查：

- 注释包/platform 不匹配；
- 去除 controls/constant 后，可进入主要 gene-level selection 的 probes 少于 60%；
- 主要 gene universe 少于 10,000 个唯一 Entrez；
- top genes 大量缺少当前名称或映射发生重大冲突；
- representative 与 median-collapse 对主要结论普遍反向。

60% 和 10,000 是本项目预先设定的最低可用性边界，不是 GPL570 的普适生物学定律。实际比例必须透明报告；即使通过，也不能隐去 unmapped/multi-mapped rows。

## 9. 差异表达分析规范

### 9.1 过滤

主要 gene-level test 的过滤只依据以下 outcome-blind 条件：control、constant、annotation ambiguity、cross-hybridization tier 与代表探针规则。默认不设“至少一组平均表达多少”的 group-aware 阈值，因为 processed GCRMA 没有统一 detection-call 标准，而用组别均值过滤可能改变检验空间。

如执行低信号敏感性过滤，阈值必须仅由 pooled sample distribution 或 Affymetrix present/absent 信息（若另行取得）定义，提前写入配置，并将过滤前后结果都交付。

### 9.2 主推断与效应阈值

每个基因至少输出：

- adjusted cancer-normal log2FC；
- fold change（`2^log2FC`，并明确方向）；
- ordinary SE、moderated t、raw P、BH adjusted P；
- 95% confidence interval；
- average expression；
- TREAT statistic/P/FDR；
- representative probe、probe tier、构成 probes 数与方向一致性；
- P1/S1/S2/S3/S4 的 effect 与 direction；
- annotation version 与结果状态。

统计发现层使用 moderated test 的 BH FDR < 0.05。优先级层使用 TREAT 检验真实绝对效应超过 `log2(1.2) = 0.2630`，并要求 TREAT BH FDR < 0.05。TREAT 直接检验超过最小效应，而不是先按 P 值显著、再事后用 fold-change 过滤；limma 手册也将 1.1、1.2、1.5 倍列为常见最小 fold-change 量级。[^7]

选择 1.2 倍的理由是：本项目是 bulk microarray discovery case，标准版需要对稳定但中等效应保持敏感，同时避免把任意非零差异都列为“优先”。必须同时报告预先规定的 1.1、1.5 倍 TREAT sensitivity counts；不得因为 1.2 倍结果太多或太少而换阈值。

报告中把三层概念分开：

- **tested**：进入模型的全部唯一基因；
- **statistically discovered**：moderated BH FDR < 0.05；
- **priority effect**：TREAT 1.2-fold BH FDR < 0.05，且通过 probe/stability review。

不得把 `P < 0.05` 未校正结果称为显著；不得把 `|log2FC|` 排名前列但 CI 很宽的基因包装成稳健结果。

### 9.3 多重检验

gene-level P values 在全部 tested genes 上用 Benjamini–Hochberg 控制 FDR。Probe-level supplement 在其自身完整 probe universe 单独校正，不能把其 adjusted P 与 gene-level adjusted P 混在同一列。每个 enrichment collection/ontology/direction 也在其预定义 family 内校正，并明确 family 边界。

### 9.4 不使用的模型

- 不使用配对 t-test；无 pairing ID。
- 不使用普通逐基因 t-test 替代 limma；样本不平衡且需要经验贝叶斯稳定方差。
- 不使用 DESeq2/edgeR；输入不是 counts。
- 不把 run 编成连续日期数；它是离散 processing proxy。
- 不把 sample prefix 与 run 同时强塞入设计；prefix 与 group 完全混杂。
- 不在 DE 前用 `removeBatchEffect` 生成的新矩阵作为 response。
- 不在看到 PCA/volcano 后改 contrast、过滤、TREAT threshold 或协变量。

## 10. 可视化规范

所有图必须有固定数据来源、排序、颜色、单位、样本数、版本和 caption。颜色需色盲友好；tissue group 在所有图保持同一映射，run 使用独立 palette。图中的 row z-score、batch adjustment 或 top-gene selection 必须在图注醒目标注，不能让展示变换被误认为统计模型输入。

### 10.1 标准图清单

| Figure ID | 图 | 数据与参数 | 回答的问题 |
|---|---|---|---|
| F01 | 数据与样本流图 | GEO → 121 samples → QC → mapped genes → tested genes | 数据如何进入结果 |
| F02 | metadata completeness | 字段×组缺失热图；run×group 柱图 | 哪些协变量可用、是否混杂 |
| F03 | per-sample density/boxplot | deposited log2 GC-RMA values | 分布/尺度/单样本异常 |
| F04 | QC metric panel | median/IQR/MAD/correlation/connectivity | 异常证据是否一致 |
| F05 | primary PCA | top 5,000 MAD genes，center yes，scale no | 主要变异与 group/run |
| F06 | visualization-only adjusted PCA | 仅为展示去 run、显式保留 group | run 调整后结构；不作 DE 输入 |
| F07 | sample correlation heatmap | Pearson，distance 1-r，average linkage | 样本相似性、异常与 run block |
| F08 | design/stability panel | P1 vs S1/S2 effect scatter、rank correlation | 结果对部分混杂是否敏感 |
| F09 | MA plot | P1 all tested genes；priority 标色 | effect 是否依赖平均信号 |
| F10 | volcano | x=adjusted log2FC，y=-log10 raw P；FDR/TREAT 分层 | effect 与证据的关系 |
| F11 | priority-gene heatmap | 最多 50 genes，见下方规则 | 样本层表达模式，不作为检验 |
| F12 | top-effect forest plot | log2FC 与 95% CI，显示敏感性 effect | effect 大小与不确定性 |
| F13 | GO ORA dotplot | up/down 分开，GeneRatio 与 FDR | 过度代表的 biological processes |
| F14 | GSEA NES plot | Hallmark/许可合格集合，NES/FDR | 全排序 program-level shift |
| F15 | selected enrichment curves | 预先规则选择正/负各最多 3 个 | leading edge 与 rank distribution |

### 10.2 Heatmap 的反误导规则

Priority heatmap 最多显示 50 个 TREAT priority genes：优先各取 25 个 up/down；一侧不足时不为凑数放入不合格基因。排序先按 TREAT FDR、再按 |moderated t|、最后 Entrez ID。行内 z-score 仅为可视化，截断到 ±2.5；聚类基于 z-scored rows，column annotations 至少含 group 与 run。图注必须写明：颜色不是原始 log2FC，DE 来自未做 row z-score 的模型。

若少于 10 个 priority genes，heatmap 标记 “not generated: insufficient priority genes”，不放宽门槛。完整表达与效应在表格中交付。

### 10.3 Volcano 与标签规则

Volcano 标签使用固定、可审计的规则：优先 TREAT FDR 最小的 up/down 各最多 10 个，再按 |effect| 与 ID tie-break；不得人工只标熟悉癌基因。y 轴超过绘图上限的点只做视觉截断，原始 P 保留在表中；零 P 若来自浮点下溢，用最小正值作显示并标注，不篡改统计值。

## 11. GO、KEGG、Reactome 与 GSEA

### 11.1 统一的 universe 原则

任何富集分析的背景都不是“全人类基因组”，也不是 GPL570 的全部历史注释。资源 R 的正式 universe 定义为：

> 通过本项目 gene-level QC、实际进入 P1 检验、且在资源 R 中有至少一个可用注释的唯一 Entrez genes。

这能控制“芯片能测到什么”和“本项目实际检验了什么”的选择机制。GO 官方指南也明确建议在实验只可能选择某个基因子集时提供自定义 reference list。[^11]

所有基因集在与 universe 交集后再应用 size threshold；报告原始 size、intersected size、命中数。输入 ID 和输出 ID 都保存，禁止只交 symbol。

### 11.2 GO over-representation analysis

标准版 primary ORA 为 GO Biological Process；Cellular Component 与 Molecular Function 进入 supplement。参数固定为：

- up 与 down 分开；输入为 TREAT 1.2-fold BH FDR < 0.05 的 priority genes；
- one-sided hypergeometric/Fisher over-representation；
- gene-set size 10–500（与 universe 相交后）；
- BH correction 在 ontology × direction 的预定义 family 内分别进行；
- report FDR < 0.05，同时交付全部被检验 terms；
- gene ratios 同时给出 numerator/denominator，避免只给小数；
- 保留 term ID、term name、gene count、universe count、P、FDR、gene IDs、annotation snapshot。

若某方向可映射 priority genes 少于 10 个，则该方向 ORA 标为 `not_evaluable`，不降低 DE 或 enrichment threshold。`clusterProfiler` 可作为实现工具，但 universe、ID 类型和多重检验 family 必须由项目显式控制，不能依赖不透明默认值。[^12]

GO terms 高度重叠。主报告的 term 精简只能用于展示：按 semantic similarity/overlap 聚类，每簇选 FDR 最小且名称最具体的代表 term；完整未经精简的表必须保留，且不得把精简后的 term 数重新当成多重检验基数。

### 11.3 KEGG 商业许可闸门

KEGG 的官方 REST 页面和 legal/subscribe 页面将数据使用，特别是非学术/商业用途，与许可条件绑定。这个项目明确用于商业生物信息服务作品集，因此默认政策是：**未有书面许可或客户可验证订阅时，不运行 KEGG，不缓存 KEGG pathway definitions，不再分发 KEGG maps。**[^16]

若 licence gate 通过，才执行 KEGG ORA：

- 锁定 organism=`hsa`、取得日期、API/数据库 release、许可证明的位置；
- 使用与 KEGG annotation 相交后的 tested Entrez universe；
- up/down 分开、size 10–500、BH FDR 0.05；
- 完整表包含 pathway ID、名称、命中/背景、P/FDR、gene IDs；
- 客户包是否可包含 pathway 图片由许可逐项决定。

若 gate 未通过，结果表中写 `KEGG_NOT_RUN_LICENSE_GATE`，并在标准版使用 Reactome pathway ORA。Reactome 官方将其数据以 CC0 提供，但引用、版本和 logo/软件条款仍需遵守。[^17]

### 11.4 GSEA / preranked enrichment

GSEA 不预先筛选 DE genes。标准 rank vector 包含全部 P1 tested unique Entrez genes，排名统计量为 P1 的 cancer-normal moderated t statistic：正值向 cancer-up，负值向 normal-up。选择 t 而不是 log2FC，是因为它同时包含效应方向与估计精度，且与主要模型一致。经典 GSEA 的目的正是检测协同但可能单基因效应不强的 gene sets。[^13][^14]

标准 primary collection 为当次许可确认可用于该商业交付的 MSigDB Hallmark snapshot；其 50 个集合经设计用于降低冗余。必须保存 collection release（设计时可见版本为 2026.1.Hs）、下载/包版本、source gene-set licence metadata 与 attribution。若某个集合有额外限制，就从商业包排除并记录；不能因为 `msigdbr` 能下载就推定所有上游内容均可任意商用。[^15]

若 Hallmark licence gate 未通过，标准替代为 GO BP 或 Reactome 的 preranked GSEA，并清楚改名，不能仍称 Hallmark。

`fgseaMultilevel` 参数预先固定：

- statistic：P1 moderated t；
- stable order：t 从大到小，完全相同再按 Entrez ID；不加入随机 jitter；
- minSize=15、maxSize=500，均在 universe intersection 后计算；
- score exponent=1；
- BH correction 在每个预定义 collection 内进行；
- 报告 NES、P、FDR、set size、leading-edge IDs、方向与 rank ties 比例；
- exact package/version、random seed 和 session info 写入记录，即便 multilevel 算法多数步骤是确定的。[^13]

如果 ties 占 tested genes 超过 1%，暂停解释并确认 statistic 精度与排序；必要时用平均秩敏感性分析，但不人为添加可改变结果的随机噪声。

### 11.5 Enrichment 的解释限制

- ORA/GSEA 说明基因集在统计排名中的富集，不证明 pathway activation、因果机制或药物靶点有效。
- 多个相似 GO terms 不是多个独立生物学发现。
- bulk tissue 富集可由细胞组成差异驱动。
- `normal-up` 应解释为在正常组织中相对更高，而不是“抑癌通路被关闭”。
- leading edge 用于解释集合驱动基因，不是经过独立 FDR 控制的新 biomarker list。

## 12. 参数选择、敏感性与稳定性判据

### 12.1 关键参数登记表

| 决策 | 主设置 | 为什么 | 预先敏感性 |
|---|---|---|---|
| 数据层 | deposited log2 GC-RMA | 范围受控、GEO 明确处理 | raw CEL 为扩展 |
| 组别方向 | cancer-normal | 单一、可解释 contrast | 不反复切换符号 |
| 批次 | 13-level run proxy fixed effect | 有部分但非完全混杂 | group-only；mixed runs |
| PCA genes | top 5,000 pooled MAD | 聚焦有信息变异，不用标签 | 1,000 / all eligible |
| PCA scaling | center yes，scale no | 避免放大低变异噪声 | 无需把 scale=yes 当主图 |
| gene key | Entrez | 稳定且 enrichment 兼容 | symbols 仅展示 |
| duplicate probes | specificity + pooled median | outcome-blind、可审计 | median collapse；probe-level |
| DE | limma trend+robust | 微阵列连续 log2 值、稳定方差 | no-trend；必要时 group-aware variance |
| multiplicity | BH FDR 0.05 | 控制探索性发现比例 | 完整 P/FDR 分布 |
| minimum effect | TREAT 1.2-fold | 直接检验最小相关效应 | 1.1 / 1.5-fold |
| ORA set size | 10–500 | 排除过小不稳定与过泛集合 | 完整 size 审计 |
| GSEA set size | 15–500 | 与常见 preranked 实践一致 | collection-specific audit |
| GSEA rank | moderated t | 方向+精度，与主模型同源 | log2FC rank supplement 可选 |

这些选择全部在看到 DE/enrichment 结果之前冻结。任何变更通过 `decision_log.md` 新增记录，说明触发证据、是否改变主要结果、批准人和版本号；不直接改掉旧配置。

### 12.2 稳定性交通灯

对 P1 与 S1、P1 与 S2，在共同 tested genes 上计算 Spearman rank correlation、Pearson effect correlation、sign concordance；对 P1 priority genes 单独报告 sign concordance 和 effect attenuation。S2 不以“仍显著”作为必要条件。

| 状态 | 预设证据 | 发布规则 |
|---|---|---|
| Green | effect Spearman ≥0.90，且 P1 priority sign concordance ≥95% | 可按主要模型发布，并保留局限 |
| Amber | Spearman 0.80–<0.90，或 sign concordance 80–<95% | 可交付，但 headline 限于跨分析一致 genes/pathways；正文突出批次敏感性 |
| Red | Spearman <0.80，或 sign concordance <80%，或大量效应反转 | 不发布稳健 cancer-normal gene list；升级 raw-CEL/设计调查 |

数值是本项目的 operational acceptance criteria，不宣称为领域统一标准。还必须展示 scatter plots 和 CI，避免单一汇总阈值掩盖局部反转。

### 12.3 样本影响与 leave-one-run-out

增加 leave-one-run-out diagnostic：每次移除一个 processing-run level，重算 group effect 的简化 sensitivity summary。它只用于识别单个 run 对 top genes/pathways 的支配程度，不作为 13 次独立发现分析。若任一 mixed run 的移除导致大量 headline effects 反向，状态至少降为 Amber；若 2005-01-25（含 8/17 normals）移除后总体结论崩溃，必须在摘要中明确正常组支撑过度集中。

### 12.4 生物学 sanity checks

Sanity check 只用于发现错误，不用于筛选“看起来合理”的结果：

- effect sign 是否与 contrast 定义一致；
- top gene 的 sample-level expression 是否支持模型估计，而非单个样本驱动；
- gene symbol、Entrez 与 probe annotation 是否一致；
- top genes 的多个 probes 是否方向一致；
- GO/GSEA leading edges 是否由极少数重复/高相关 genes 支配；
- adjustment 前后是否出现无法解释的大量方向翻转；
- 正常样本是否几乎全部由一个 run 决定。

已知乳腺癌 marker 可以作为后验 plausibility check，但不得用来决定阈值、删除样本或只汇报符合预期的 pathways。若经典 marker 未出现，先查 annotation、方向、批次、组织组成和 power，不能把它自动判定为 pipeline 错误。

## 13. 失败条件、暂停条件与升级路径

### 13.1 Hard stop

出现以下任一项，标准分析停止，不发布 DE/enrichment：

- 官方 accession、样本数、平台或矩阵列无法一致映射；
- 文件损坏、hash 未记录、非有限值或无法确认 log2 scale；
- tissue group 或 cancer-normal contrast 无法从原始记录唯一建立；
- 设计矩阵 rank deficient，group coefficient 不可估计；
- 注释/platform 错配，或 gene-level coverage 低于第 8.5 节边界；
- 严重异常样本无法通过现有 processed 数据解释，且足以实质改变结果；
- P1 与核心敏感性达到 Red；
- gene-set 数据许可不允许预定用途；该项只停止受影响的 enrichment 模块，不必停止合法的 DE/GO 模块。

### 13.2 Pause and review

以下为人工复核，不自动失败：

- 某样本触发至少两个 outlier signals；
- run adjustment 使大量 priority genes 反向；
- 代表探针与 median collapse 大范围不一致；
- age summary 与 GEO series text 不一致；
- 上游文件哈希变化；
- GSEA ties >1%；
- GO term 数异常多且高度冗余；
- priority gene 少于 10 或多到几乎覆盖整个 universe。

每个 pause 必须产生 issue record；resolution 可以是继续、降级、排除并做敏感性、或升级 raw-CEL extension。

### 13.3 Raw-CEL 升级触发

满足以下任一条件，应报价并启动扩展版，而不是在标准版偷偷增加工作：

- 无法解释的样本分布/相关性异常；
- 需要确认 scan/image/RNA degradation/NUSE/RLE；
- deposited preprocessing provenance 不足以支撑主要结论；
- 需要 alternative CDF/Brainarray；
- 客户要求从原始文件重建全部表达值；
- 主要结果对 run proxy 极度敏感，且 raw header/quality metrics 可能澄清。

## 14. 可复现工程与技术规范

### 14.1 推荐仓库结构

项目根目录由仓库所在位置决定；所有实现均使用仓库相对路径，运行时可从激活环境的 `PATH` 或可选的 `BIOINFO_ENV_PREFIX` 解析。建议结构如下：

```text
repository-root/
├── README.md
├── LICENSE
├── CITATION.cff
├── CHANGELOG.md
├── config/
│   ├── analysis.yml
│   ├── gene_set_licenses.yml
│   ├── environment.yml
│   └── conda-explicit-linux-64.txt
├── data/
│   ├── metadata/
│   │   ├── source_manifest.tsv
│   │   ├── sample_manifest.tsv
│   │   ├── field_dictionary.tsv
│   │   └── metadata_issues.tsv
│   ├── raw/                  # gitignored; downloaded processed inputs/CEL extension
│   └── derived/              # gitignored; deterministic intermediates
├── R/
│   ├── 00_validate_config.R
│   ├── 01_acquire_validate.R
│   ├── 02_build_metadata.R
│   ├── 03_processed_qc.R
│   ├── 04_annotate_collapse.R
│   ├── 05_fit_de.R
│   ├── 06_sensitivity.R
│   ├── 07_enrichment.R
│   ├── 08_figures_tables.R
│   └── 09_acceptance_checks.R
├── tests/
│   ├── testthat/
│   └── fixtures/             # synthetic/minimal, no production claims
├── reports/
│   ├── client_report.qmd
│   ├── methods_appendix.qmd
│   └── qc_appendix.qmd
├── results/
│   ├── tables/
│   ├── figures/
│   ├── logs/
│   └── checksums/
├── docs/
│   ├── GSE42568_PROJECT_DESIGN.md
│   ├── decision_log.md
│   └── interpretation_language.md
├── _targets.R
└── .gitignore
```

本轮只创建本设计文件；上述其它文件是实现合同，不是已存在声明。

### 14.2 工作流引擎与环境

实现建议使用 R + `targets` 构建显式 DAG，并以 pinned Conda environment 作为唯一环境锁；每个 target 只写入声明的 derived/result 路径。标准运行记录：

- OS、R、Bioconductor 与全部 package versions；
- locale、timezone、BLAS、线程数；
- Git commit、配置哈希、输入哈希、gene-set snapshot；
- random seeds；
- 开始/结束 UTC、命令退出状态；
- warning/error；
- 所有最终产物 SHA-256。

原始/processed matrix 和大体积派生文件不提交 Git；提交下载清单、checksum、解析规则与小型 synthetic fixtures。代码 LICENSE 不自动覆盖 GEO 数据、注释或基因集；每个外部资产单独记录条款。

### 14.3 阶段输入输出合同

| Stage | 必需输入 | 必需输出 | 通过条件 |
|---|---|---|---|
| A Acquire | accession + URLs | locked source files、manifest、checksums | 官方来源、完整、hash 固定 |
| B Metadata | series/sample records | sample manifest、dictionary、issues | 121 GSM 一一映射，group 104/17 |
| C Matrix | series matrix | immutable parsed matrix metadata | 维度/ID/finite/scale 通过 |
| D QC | matrix + manifest | QC tables/figures/flag ledger | 无未解决 hard stop |
| E Annotation | matrix + pinned DB | mapping audit、representative map、gene matrix | coverage 与 ambiguity gate 通过 |
| F DE | gene matrix + design | full P1/S1–S4 results、model diagnostics | contrast estimable、统计列完整 |
| G Enrichment | full tested genes + licensed sets | ORA/GSEA full tables + plots | universe/版本/licence 全记录 |
| H Report | validated artifacts | HTML/PDF/Markdown + machine-readable package | acceptance checks 全通过 |

每个 stage 都必须验证 schema、行数、唯一键、NA、排序和上游 hash。下游不直接读取“目录里最新文件”，而读取 manifest/target 指向的确定对象。

### 14.4 配置中必须显式声明的值

- accession、platform、expected sample groups；
- contrast direction；
- processing-run proxy 解析表与例外；
- annotation package/source/version；
- probe suffix tier 与 representative rule；
- limma trend/robust/TREAT 设置；
- FDR family 与 threshold；
- PCA/heatmap/top-label rules；
- GO/GSEA set-size 与 universe；
- KEGG/Hallmark licence gate 状态；
- outlier flag 与 exclusion policy；
- sensitivity traffic-light thresholds；
- report language与允许的 interpretation phrases。

配置改变必须使下游 targets 失效并重新计算，不能只重画一张图。

### 14.5 自动化测试

至少覆盖：

- download/hash mismatch 会失败；
- GSM 顺序改变会按 ID 对齐而非错位；
- 重复/缺失 ID 会失败；
- normal 被误标成 healthy 或 paired 的 schema 检查；
- log2 matrix 被再次 log 的 guard；
- one-to-many annotation 不会自动取第一行；
- representative probe 选择不读取 group/DE columns；
- contrast sign 的小型已知 fixture；
- universe 只含 tested-and-annotated genes；
- license gate 未通过时 KEGG target 不运行；
- 报告表与 figure 引用的 artifact hashes 一致；
- 重跑在相同环境下产生相同 key tables/checksums，允许仅时间戳不同。

CI 只运行 synthetic/minimal fixture 与静态检查；全量生物结果由受控 runner 生成。CI 通过不等于数据质量/生物学结论通过。

## 15. 标准客户交付包与 GitHub 作品集

### 15.1 标准版包含

- 数据来源、版本、checksum、sample manifest 与字段字典；
- processed matrix 完整性与样本级 QC；
- processing-run/group confounding audit；
- pinned annotation、mapping audit 与 probe-collapse sensitivity；
- P1 主 DE、S1–S4 敏感性和 leave-one-run-out diagnostic；
- 全部 DE table、priority table、probe-level supplement；
- PCA、correlation、QC、MA、volcano、heatmap、effect forest figures；
- GO BP ORA（CC/MF supplement）；
- 经许可的 Hallmark GSEA，或明确命名的 GO/Reactome GSEA fallback；
- KEGG license gate；许可通过才运行，否则 Reactome ORA 替代；
- client report、methods appendix、QC appendix；
- session info、config、logs、checksums、decision log、limitations；
- 可从官方输入重建结果的仓库结构与运行说明。

标准版只包含一个预定义主要对比 `cancer-normal`、一次预注册参数集和问题修复后的最终重跑。结果导向的反复换阈值、挑亚组或增加 endpoint 不属于标准版。

### 15.2 GitHub 公开作品集应呈现

README 的顺序建议为：问题与边界 → 数据 provenance → 设计中的最大风险（normal provenance、run confounding）→ 可复现命令 → QC gate → 主要结果导航 → limitations → 如何复用模板。README 不把“发现最多 DEG”当卖点，而强调为什么该结果可审计、何时应停止、哪些内容没有被过度承诺。

公开仓库不要提交大 matrix/CEL；提供官方下载步骤、hash manifest 和合法的小型 fixture。每个展示图链接到生成它的 table/config，release 附带最终 checksums。代码 LICENSE、GEO data attribution、MSigDB/KEGG/Reactome/GO 条款分别说明。

### 15.3 客户版验收标准

交付被视为完成必须同时满足：

1. 报告中样本数、组别、方向、平台与 manifest 一致；
2. 所有输入与最终产物有 checksum；
3. 所有 hard-stop checks 通过或项目以明确 failure report 结束；
4. 主设计与敏感性结果可从完整表复核；
5. 每张主图可定位到 source table、参数和脚本 stage；
6. 所有 gene/pathway 结果有稳定 ID、版本与 universe；
7. licensing gates 有可审计记录；
8. limitations 包含 normal provenance、年龄不可调整、非配对、bulk composition 与 run partial confounding；
9. 未交付任何未经范围批准的临床预测/因果主张；
10. 在锁定环境中 clean rebuild 成功，key output checksums 一致。

如果科学质量闸门失败，合格交付可以是“不可支持预定结论”的审计报告；不能为了形成漂亮作品集而绕过失败条件。

## 16. 扩展版清单与范围控制

以下均不属于首个标准 case，需单独目标、报价和验收：

| 扩展 | 新增价值 | 新增风险/需求 |
|---|---|---|
| Raw CEL 重处理 | NUSE/RLE、RNA degradation、image/spatial、统一 preprocessing | 约 935 MB CEL、更多计算与方法选择 |
| Alternative CDF/Brainarray | 更现代的一基因一 probe-set 映射 | 必须从 CEL 重汇总，结果不可与 deposited matrix 混称同一版本 |
| PAM50/分子亚型 | 解释肿瘤异质性 | normal 不参与同一问题；需验证 centroid/platform scaling |
| Tumor-only ER/grade/node | 临床病理关联 | 新 contrasts、多重检验、混杂设计 |
| Survival analysis | RFS/OS 关联 | censoring、临床协变量、模型假设、过拟合 |
| Deconvolution | 评估 immune/stromal composition | reference/平台偏差；不能等同真实细胞比例 |
| WGCNA/network | 共表达模块 | 样本量、批次与参数敏感，解释成本高 |
| 外部验证/meta-analysis | 评估跨队列复现 | 平台 harmonization、队列定义、批次与许可 |
| ML classifier | 预测任务 | 独立 test set、nested CV、calibration；本数据不应自证临床性能 |
| 组织来源追溯 | 澄清 17 normals | 需联系作者、伦理/样本记录，不一定可获得 |

标准版中不得“顺手”做 survival、classifier、WGCNA、药物预测或因果网络。这些会改变 estimand、验证要求和客户承诺。

## 17. 混杂因素、数据陷阱与解释语言

### 17.1 核心不可消除局限

1. **Normal provenance 未报告充分。**不能区分健康、癌旁或其它来源。
2. **年龄无法在组间调整。**104 tumors 有年龄，17 normals 全缺失；年龄可能是 residual confounder。
3. **Processing run 与 group 部分混杂。**所有 normal 集中在 6 个 mixed runs，且 8/17 normal 位于单一 2005-01-25 层级。
4. **Bulk tissue composition。**差异可能反映 epithelial、immune、stromal、adipose 比例，而非同一细胞类型内调控。
5. **非配对设计。**个体间差异进入组间比较。
6. **老平台与动态注释。**probe specificity 与基因映射不等同于现代 RNA-seq transcript quantification。
7. **存档选择过程不完全。**只分析 deposited 104/17，不能为未公开的纳入/排除机制补故事。
8. **横断面观察性数据。**无法证明表达改变导致癌症，也无法直接证明临床效用。

### 17.2 常见技术陷阱

- 把 processed values 当 raw intensity 再 log/normalize；
- 以 sample title 前缀作为可调整 batch；
- 用癌/正常标签选择 variable genes 之后再 PCA，造成监督泄漏；
- 从每个基因挑最显著 probe；
- 用全人类基因组作为 ORA universe；
- 将同一基因的多个 probes 重复送入 enrichment；
- 只按 nominal P < 0.05 汇报；
- 在 heatmap 的 z-score 上计算 fold change；
- 因 PCA 没有完美分组而删样本或追加 normalization；
- 把 ComBat-adjusted display matrix 当正式 inferential response；
- 以显著性 retention 评价 n=62 sensitivity，忽略 power 下降；
- 将 KEGG/MSigDB 包可安装误认为商业再分发自动获准；
- 从 GO/GSEA 推导 pathway 被激活、药物会有效或机制已证实。

### 17.3 允许与禁止的措辞

| 推荐措辞 | 禁止/需证据升级的措辞 |
|---|---|
| “在 GSE42568 中，癌组织相对正常乳腺组织的表达差异” | “癌症导致该基因改变” |
| “与 cancer status 相关” | “诊断 biomarker” |
| “gene set 在排序中富集” | “通路被激活/关闭” |
| “normal breast tissue；来源未充分说明” | “healthy/matched/adjacent controls” |
| “调整 processing-run proxy 后” | “消除了所有批次效应” |
| “在主要与敏感性模型中方向一致” | “已在独立队列验证” |
| “bulk-tissue signal” | “tumor-cell-intrinsic mechanism” |

## 18. 完整表格、图和报告文件清单

### 18.1 必交表格

| Table ID | 文件名 | 内容 |
|---|---|---|
| T01 | `source_manifest.tsv` | 来源、时间、大小、hash、条款 |
| T02 | `sample_manifest.tsv` | 121 样本、group/run/covariates/包含状态 |
| T03 | `metadata_issues.tsv` | 年龄汇总冲突、缺失、推导字段与 resolution |
| T04 | `matrix_integrity.tsv` | 维度、ID、finite、range、constant/control |
| T05 | `sample_qc_metrics.tsv` | distribution/PCA/correlation/connectivity flags |
| T06 | `sample_exclusion_ledger.tsv` | 证据、reviewer、with/without 状态 |
| T07 | `design_diagnostics.tsv` | run×group、rank、estimability、df |
| T08 | `probe_mapping_audit.tsv` | 全 probes 映射类别与 suffix |
| T09 | `probe_selection_audit.tsv` | 每个 Entrez candidates 与代表探针规则 |
| T10 | `de_primary_all_genes.tsv` | P1 全部 tested genes |
| T11 | `de_priority_genes.tsv` | TREAT priority + probe/stability flags |
| T12 | `de_probe_level_supplement.tsv.gz` | 合格 probes 的技术 supplement |
| T13 | `de_sensitivity_summary.tsv` | P1/S1–S4 effect、rank、sign |
| T14 | `leave_one_run_out_summary.tsv` | 每次去一 run 的 influence |
| T15 | `go_ora_all_terms.tsv` | BP/CC/MF、up/down 全部 terms |
| T16 | `pathway_ora_all_terms.tsv` | KEGG（有许可）或 Reactome，状态字段明确 |
| T17 | `gsea_all_gene_sets.tsv` | 全 collections、NES/FDR/leading edge |
| T18 | `gene_set_manifest.tsv` | collection release、license、universe、size |
| T19 | `acceptance_checks.tsv` | 每个硬/软闸门及证据路径 |
| T20 | `artifact_checksums.tsv` | 最终产物 hash |

任何表如果模块不运行，也要生成 status row，说明 `not_applicable`、`not_evaluable` 或 `license_gate`；不允许文件静默缺失。

### 18.2 必交图

主图 F01–F15 的 PDF 与 300-dpi PNG；用于网页的 SVG 只在不包含不适宜嵌入的字体/敏感元数据时生成。每图配同名 caption text 和 source-data TSV。完整 QC supplement 可额外包含 PCA sensitivity、Spearman heatmap、probe concordance 与 enrichment redundancy map。

### 18.3 必交报告

- `client_report.html` 与 `client_report.pdf`：决策、主要结果、限制、可行动结论；
- `methods_appendix.html/pdf`：全部参数、模型、版本、universe、许可；
- `qc_appendix.html/pdf`：样本/矩阵/annotation/design/stability 审计；
- `machine_readable_manifest.json`：所有 artifact、hash、schema/version；
- `sessionInfo.txt`、`pipeline.log`、`decision_log.md`；
- GitHub `README.md`：公开复现入口，不替代正式报告。

## 19. 实施顺序与审查点

1. **Freeze**：获取官方 matrix/metadata/platform，写 checksum 与 terms。
2. **Identity review**：确认 121 GSM、104/17、GPL570、log2 GCRMA、字段缺失。
3. **Design review**：冻结 contrast、run proxy、非配对状态、不可调整变量。
4. **Processed QC**：分布、PCA、correlation、outlier ledger；决定继续或 raw extension。
5. **Annotation review**：冻结映射、suffix tiers、coverage 和 representative selection。
6. **Model dry run**：只检查 design rank、schema、contrast sign，不看/不解释 top genes。
7. **Locked full run**：P1、S1–S4、leave-one-run-out；生成不可变表。
8. **Statistical review**：FDR/TREAT、probe concordance、stability traffic light。
9. **License review**：GO/Hallmark/KEGG/Reactome 逐项 gate 后运行 enrichment。
10. **Interpretation review**：每个 headline 回到 sample-level、probe 和 sensitivity 证据。
11. **Reproducibility review**：clean rebuild、hash、报告链接、CI fixture。
12. **Release**：客户包与 GitHub release 分开，分别检查允许公开的资产。

每个 review point 都能停止项目。实施团队不能把全部脚本跑完以后再补 provenance、批次设计或 licensing。

## 20. 需求覆盖矩阵

| 用户要求 | 本文件位置 |
|---|---|
| 1. 实验、数据、样本、平台核验 | 第 1.1、3、5 节 |
| 2. GEO processed 数据是否适合 | 第 4 节 |
| 3. cancer vs normal 设计 | 第 2、6 节 |
| 4. QC 方案 | 第 7、13 节 |
| 5. 注释与重复探针 | 第 8 节 |
| 6. 差异表达 | 第 9 节 |
| 7. PCA/相关性/heatmap/volcano | 第 7、10 节 |
| 8. GO/KEGG/GSEA 与 universe/rank | 第 11 节 |
| 9. 参数理由 | 第 12.1 节 |
| 10. sanity/failure criteria | 第 12.2–12.4、13 节 |
| 11. 文件/表/图清单 | 第 14、18 节 |
| 12. 混杂与陷阱 | 第 17 节 |
| 13. 标准版与扩展版 | 第 15、16 节 |
| 14. 技术规范 | 第 14、19 节 |

## 21. 最终放行决策模板

正式结果报告首页必须给出以下一个结论，不留模糊状态：

- **PASS — standard release**：所有 hard stops 通过，稳定性 Green，许可完整；
- **PASS WITH LIMITATIONS**：所有 hard stops 通过，稳定性 Amber；仅发布跨敏感性一致结论；
- **HOLD — raw-CEL or provenance review required**：processed 数据不足以解决关键质量/混杂问题；
- **FAIL — estimand not supportable**：组别不可估计、身份/平台错误或关键许可/数据问题无法解决。

即使 PASS，结论上限仍是“该队列中的稳健表达关联”。商业价值来自可追溯、可复算、知道何时不该下结论，而不是把一个公共微阵列数据集包装成临床产品。

## Sources

[^1]: NCBI Gene Expression Omnibus. “GSE42568: Breast Cancer Gene Expression Analysis.” Submitted 2012-11-27; public 2013-05-26; updated 2019-03-25. Accessed 2026-09-11. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42568

[^2]: NCBI Gene Expression Omnibus. “Downloading GEO Data.” NCBI. Accessed 2026-09-11. https://www.ncbi.nlm.nih.gov/geo/info/download.html

[^3]: NCBI Gene Expression Omnibus. “GPL570: [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array.” Platform record; annotation update shown by GEO. Accessed 2026-09-11. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570

[^4]: Clarke C, Madden SF, Doolan P, et al. “Correlating transcriptional networks to breast cancer survival: a large-scale coexpression analysis.” *Carcinogenesis*. 2013;34(10):2300–2308. doi:10.1093/carcin/bgt208. Full-text repository copy accessed 2026-09-11. https://www.tara.tcd.ie/bitstreams/6671e241-a8e6-486b-88e1-944a49fd00f6/download

[^5]: PubMed. “Correlating transcriptional networks to breast cancer survival: a large-scale coexpression analysis.” PMID 23740839. U.S. National Library of Medicine. Accessed 2026-09-11. https://pubmed.ncbi.nlm.nih.gov/23740839/

[^6]: Wu Z, Irizarry RA, Gentleman R, Murillo FM, Spencer F. “A Model-Based Background Adjustment for Oligonucleotide Expression Arrays” / Bioconductor `gcrma` vignette. Bioconductor. Accessed 2026-09-11. https://bioconductor.org/packages/release/bioc/vignettes/gcrma/inst/doc/gcrma2.0.pdf

[^7]: Ritchie ME, Phipson B, Wu D, et al.; limma authors. “limma: Linear Models for Microarray Data — User’s Guide.” Bioconductor package manual, release documentation accessed 2026-09-11. https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf

[^8]: Bioconductor. “AnnotationDbi: Manipulation of SQLite-based annotations in Bioconductor.” Package manual, release documentation accessed 2026-09-11. https://bioconductor.org/packages/release/bioc/manuals/AnnotationDbi/man/AnnotationDbi.pdf

[^9]: Bioconductor. “hgu133plus2.db: Affymetrix Human Genome U133 Plus 2.0 Array annotation data.” Package page/manual, release documentation accessed 2026-09-11. https://bioconductor.org/packages/release/data/annotation/html/hgu133plus2.db.html

[^10]: Affymetrix / Thermo Fisher Scientific. “GeneChip Expression Analysis: Data Analysis Fundamentals” (probe-set suffix and specificity guidance). Accessed 2026-09-11. https://documents.thermofisher.com/TFS-Assets/LSG/manuals/data_analysis_fundamentals_manual.pdf

[^11]: Gene Ontology Consortium. “GO Enrichment Analysis.” Gene Ontology Resource. Accessed 2026-09-11. https://www.geneontology.org/docs/go-enrichment-analysis/

[^12]: Bioconductor. “clusterProfiler: A universal enrichment tool for interpreting omics data.” Package manual, release documentation accessed 2026-09-11. https://bioconductor.org/packages/release/bioc/manuals/clusterProfiler/man/clusterProfiler.pdf

[^13]: Bioconductor. “fgsea: Fast Gene Set Enrichment Analysis.” Package manual, release documentation accessed 2026-09-11. https://bioconductor.org/packages/release/bioc/manuals/fgsea/man/fgsea.pdf

[^14]: Subramanian A, Tamayo P, Mootha VK, et al. “Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles.” *Proceedings of the National Academy of Sciences*. 2005;102(43):15545–15550. https://pmc.ncbi.nlm.nih.gov/articles/PMC1239896/

[^15]: Broad Institute / UC San Diego. “Molecular Signatures Database (MSigDB): Collections and License Terms.” Release information and terms accessed 2026-09-11. https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp and https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp

[^16]: Kanehisa Laboratories. “KEGG API” and “KEGG Copyright and Licensing.” Kyoto Encyclopedia of Genes and Genomes. Accessed 2026-09-11. https://www.kegg.jp/kegg/rest/ and https://www.kegg.jp/kegg/legal.html

[^17]: Reactome. “Reactome License Agreement.” Reactome Knowledgebase; CC0 data terms. Accessed 2026-09-11. https://reactome.org/license
