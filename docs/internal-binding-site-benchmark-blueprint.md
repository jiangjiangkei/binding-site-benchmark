# Binding Site Prediction Benchmark 内部开发蓝图

## 1. 目标

本项目的目标不是再做一个“方法列表”，而是建立一个可以稳定复跑、可扩展、可比较的 **protein-small molecule binding site prediction benchmark**。  
v1 先解决两件事：

1. 建一个可信、可冻结、低泄漏的测试集。
2. 复现一组当前常用且有代表性的 binding site prediction 方法，并用统一协议比较它们。

这份文档按“先做什么、产出什么、如何验收”来写，默认读者是内部开发同事。

## 2. v1 范围定义

### 2.1 纳入范围

- 任务类型：protein-small molecule binding site prediction
- 输入：单个蛋白结构
- 输出：binding site center、binding residues、或 pocket volume/voxels
- 结构来源优先级：
  - `v1-core`: 实验结构为主，优先 holo 结构
  - `v1-extended`: apo 结构、AlphaFold2 或其他 predicted structures
- 评测对象：可本地运行、可批处理、公开代码的方法

### 2.2 暂不纳入

- protein-protein binding site prediction
- peptide、DNA/RNA、金属离子主导的位点预测
- docking、pose prediction、binding affinity prediction
- 只能通过网页服务器运行、无法批量复现的方法

### 2.3 关键假设

为了把项目快速做实，**v1 只做 small-molecule binding site benchmark**。  
这是当前最成熟、方法最集中、数据也最容易标准化的子任务。

## 3. 为什么这样切

当前 ligand binding site prediction 的公开方法已经非常多。2024 年的系统比较工作集中评测了 13 个主流方法，涵盖了 `fpocket`、`P2Rank`、`DeepPocket`、`PUResNet`、`GrASP`、`IF-SitePred`、`VN-EGNN` 等，并明确指出：

- 现在的主流方法已经从几何法转向机器学习和深度学习。
- benchmark 最大的问题之一不是“没有方法”，而是 **测试集定义不一致、泄漏控制不足、输出格式不统一**。
- `top-N+2 recall` 值得作为通用主指标之一。

因此，v1 的核心不是追求方法越多越好，而是优先把 **数据定义、评测协议、复现接口** 做扎实。

## 4. v1 成功标准

项目达到以下条件，就可以认为 v1 成功：

### 4.1 数据侧

- 产出一个冻结版测试集：`test_core_v1`
- 每个 target 都有完整元数据：
  - `pdb_id`
  - `chain_id`
  - `uniprot_id`
  - `ligand_id`
  - `release_date`
  - `structure_type`
  - `site_residue_set`
  - `site_center`
  - `source_db`
- 完成训练泄漏检查并生成审计报告

### 4.2 方法侧

- 至少复现 4 个方法并跑通全量 benchmark
- 每个方法都有固定版本、固定环境、固定输入输出适配器
- 每个方法都能在 20 个样本的 smoke set 上稳定运行

### 4.3 评测侧

- 有统一的 residue-level prediction schema
- 有统一的 residue-level evaluation script
- 能输出：
  - residue-level AUROC / AUPRC
  - residue-level precision / recall / F1
  - residue-level top-(N+2) recall
  - top1 / top3 / top5 success
  - DCA / DCC
  - runtime
  - failure rate

## 5. 方法清单与优先级

### 5.1 方法分层

#### P0：必须先跑通的稳定基线

- `fpocket`
- `P2Rank`
- `P2Rank + rescoring / PRANK 路线`
- `DeepPocket`

这 4 个方法是 v1 的主干，原因很直接：

- 社区使用广
- 有公开代码
- 可本地批量运行
- 在近期比较工作里仍然是标准参照系

#### P1：第二批主流深度学习方法

- `PUResNetV2.0`
- `GrASP`
- `IF-SitePred`
- `ProtCross`

这几类方法代表当前更“现代”的深度学习路线，但依赖更复杂，建议在 P0 跑稳后再接入。

#### P2：前沿/实验性方法

- `VN-EGNN`

它代表最新一类等变 GNN 路线，可以作为扩展项，但不应阻塞 v1 交付。

### 5.2 纳入规则

方法是否进入 benchmark，统一看下面 5 条：

1. 有公开实现和明确论文出处。
2. 能在本地运行，不依赖人工网页操作。
3. 可以接受批量输入。
4. 输出可被规范化为统一 schema。
5. 许可证和模型权重可用于内部研究复现。

## 6. 测试集蓝图

## 6.1 总体策略

v1 不建议直接从某一篇论文的测试集原样照搬。  
推荐做法是：

1. 用 `BioLiP` 作为 biologically relevant ligand 的主要来源。
2. 用 `PDB` 原始结构作为结构输入来源。
3. 用 `PLINDER` 的思路做泄漏控制、apo/predicted 扩展和后续泛化子集设计。
4. 参考 `LIGYSIS` 的经验，避免只看单个 asymmetric unit 和单个复合物视角。

### 6.2 推荐的数据分层

#### `test_core_v1`

这是必须先做的主测试集，特点如下：

- holo 结构
- 实验结构
- small-molecule ligand
- 高质量结构
- 尽量单链可解释
- 尽量避免跨链共同形成的口袋

#### `test_apo_v1`

- 与 `test_core_v1` 中一部分 target 配对的 apo 结构
- 用于测试方法对构象变化的鲁棒性

#### `test_pred_v1`

- 与 `test_core_v1` 配对的 AlphaFold2 或其他 predicted structures
- 用于测试方法迁移到非实验结构时的性能

结论很明确：  
**v1 先冻结 `test_core_v1`，`apo` 和 `predicted` 作为扩展集并行准备。**

## 6.3 `test_core_v1` 具体筛选规则

### A. 结构筛选

- 只保留实验结构
- 优先 `X-ray`
- 分辨率建议 `<= 2.8 Å`
- 删除明显缺失关键口袋区域的结构
- 删除主链严重断裂、残基编号异常、无法标准解析的结构

### B. 配体筛选

- 仅保留 biologically relevant small molecules
- 排除：
  - buffer/additive
  - glycerol、ethylene glycol 一类结晶添加剂
  - 单纯金属离子
  - DNA/RNA
  - 长肽
  - 共价修饰如果无法和普通 pocket 统一评测则先剔除

### C. 位点筛选

- 只保留可清晰映射到蛋白表面 pocket 的位点
- 暂不纳入需要多条链共同定义且边界模糊的口袋
- 一个 target 最好只保留 1 个主位点用于 `core` 评测
- 多位点蛋白可以进入扩展集，但不建议先放进 v1-core

### D. 冗余控制

- 按 `UniProt` 聚合蛋白
- 先做 sequence clustering，建议 `30%` 或 `40%` identity 阈值
- 再做 pocket-level 去冗余
- 同一蛋白的多个高度相似复合物，只保留一个代表结构

### E. 时间切分与泄漏控制

这是最关键的部分。

因为很多方法的训练集来自 `scPDB`、`PDBbind`、`HOLO4K`、`COACH420` 或其衍生版本，v1 必须显式防止“测试集本来就在别人训练里”。

建议做两层泄漏控制：

1. **时间切分**
   - 测试集优先选择较新的 PDB release
   - 建议从 `2023-01-01` 或更晚开始筛候选
2. **相似性切分**
   - 对照公开训练源做 sequence-level overlap 检查
   - 能做的话再补 pocket similarity / ligand similarity 审计

### F. 推荐规模

`test_core_v1` 建议先做：

- `300-500` 个蛋白 target

这个规模的优点：

- 足够支撑稳定排序
- 不至于把方法复现和批处理压垮
- 手工抽查仍然可做

## 6.4 标注定义

每个 reference site 至少要有 3 层标签。v1 的主评估口径固定为
**residue-level binding residues**；site center / DCA / DCC 只作为兼容
pocket-center 方法的辅助视角。

### 1. Ligand 原子集合

- 直接保存参考 ligand 的 heavy atoms
- 用于 DCA 评测

### 2. Binding residues

- 定义为与 ligand 任一原子距离 `<= 4.0 Å` 或 `<= 4.5 Å` 的残基
- 阈值必须写死在配置里，不要在不同实验间漂移

### 3. Site center

建议同时保存两种中心：

- `ligand_centroid`
- `binding_residue_centroid`

主评测先统一使用一种，另一种保留用于敏感性分析。  
如果要和多数 pocket-center literature 对齐，优先保留 ligand-centered 的定义。

## 6.5 人工质控

自动流程之后必须加人工质控，不然测试集会被脏样本拖垮。

建议：

- 从候选集中随机抽 `50` 个样本
- 每个样本人工确认：
  - ligand 是否生物学相关
  - pocket 是否合理
  - site residues 是否明显错位
  - target 是否存在多链依赖问题

如果 50 个样本里错误率高于 `10%`，就不能冻结测试集，必须回到规则层调整。

## 7. 方法复现蓝图

## 7.1 复现原则

每个方法只做两件事：

1. 跑出作者默认推荐的结果。
2. 统一转成 benchmark 自己的标准输出。

v1 不要一开始就做大规模调参。  
先跑默认配置，后面再做 `default track` 和 `tuned track` 两条赛道。

## 7.2 复现顺序

### 第一批

1. `fpocket`
2. `P2Rank`
3. `P2Rank-rescore / PRANK 路线`
4. `DeepPocket`

### 第二批

5. `PUResNetV2.0`
6. `GrASP`
7. `IF-SitePred`

### 第三批

8. `VN-EGNN`

## 7.3 每个方法都要补齐的元信息

每接入一个方法，都必须记录：

- 论文链接
- 代码仓库链接
- commit hash / release tag
- 模型权重来源和校验值
- 运行环境
- CPU/GPU 需求
- 输入格式
- 原始输出格式
- 已知失败模式
- 官方默认参数

## 7.4 方法适配层设计

仓库里不要把方法逻辑和 benchmark 逻辑混在一起。  
建议每个方法都走统一的 adapter 模式。

### 标准输入

每个方法 adapter 只接收：

- `target_id`
- `input_structure_path`
- `output_dir`
- `config_path`

### 标准输出

统一输出成一个 `json` 或 `csv`，至少包含：

- `target_id`
- `method`
- `pocket_rank`
- `score`
- `center_x`
- `center_y`
- `center_z`
- `residue_ids`
- `raw_output_path`
- `runtime_sec`
- `status`

### 状态定义

- `success`
- `partial_success`
- `failed_parse`
- `failed_runtime`
- `failed_oom`
- `failed_unknown`

## 7.5 为什么一定要做标准化输出

因为这些方法原生输出差异很大：

- 有的给 pocket center
- 有的给 residues
- 有的给 voxel segmentation
- 有的不给 rank
- 有的不给 score

如果不做统一 prediction schema，后面所有评测代码都会变成一次性脚本。

## 7.6 环境管理建议

不要试图把所有方法塞进一个 conda 环境。  
推荐：

- 每个方法一个独立环境
- benchmark 主分析一个独立环境

建议目录：

- `envs/fpocket.yml`
- `envs/p2rank.yml`
- `envs/deeppocket.yml`
- `envs/benchmark.yml`

如果某些方法依赖太重，直接上 Docker/Singularity，不要硬拼。

## 7.7 烟雾测试

每个方法接入后先跑 `smoke set`，不要直接跑全量。

`smoke set` 建议：

- `20` 个 target
- 覆盖：
  - 小蛋白
  - 大蛋白
  - 明显 pocket
  - 边界模糊 pocket
  - 不同配体大小

只有满足下面条件，才允许进全量：

- 成功率 `>= 90%`
- 输出可以完整被 parser 解析
- 平均 runtime 在可接受范围内

## 8. 评测协议

## 8.1 主指标

v1 建议主指标如下：

1. `Top1 DCA success rate`
2. `Top3 DCA success rate`
3. `Top5 DCA success rate`
4. `Top-N+2 recall`

其中 `Top-N+2 recall` 建议作为全局主指标之一。  
这里的 `N` 是 reference sites 数量，`+2` 是为了减轻不同方法输出 pocket 数不同带来的不公平。

## 8.2 辅指标

- `DCC success rate`
- residue-level `precision / recall / F1`
- pocket-level overlap
- 平均预测 pocket 数
- 冗余 pocket 比例
- runtime
- peak memory
- failure rate

## 8.3 为什么主指标优先选 DCA/DCC

因为不同方法输出粒度差别太大。  
center-based metrics 是最容易统一的主评测层。

对于输出 residues/segmentation 的方法，可以额外做 residue-level 和 volume-level 评测，但不能让这些指标成为唯一主指标。

## 8.4 评测赛道

建议从一开始就定义两条赛道：

### `default-track`

- 完全使用作者默认参数
- 不允许针对测试集调参

### `tuned-track`

- 允许有限调参
- 但必须记录全部改动

v1 报告优先展示 `default-track`。

## 9. 推荐仓库结构

```text
binding-site-benchmark/
├── docs/
│   └── internal-binding-site-benchmark-blueprint.md
├── configs/
│   ├── dataset/
│   ├── methods/
│   └── evaluation/
├── data/
│   ├── raw/
│   ├── interim/
│   ├── processed/
│   ├── splits/
│   └── manifests/
├── envs/
├── methods/
│   ├── fpocket/
│   ├── p2rank/
│   ├── deeppocket/
│   ├── puresnet/
│   ├── grasp/
│   ├── if_sitepred/
│   └── vnegnn/
├── src/
│   ├── dataset/
│   ├── runners/
│   ├── parsers/
│   ├── evaluation/
│   └── utils/
├── results/
│   ├── smoke/
│   ├── benchmark/
│   └── reports/
└── scripts/
```

## 10. 具体开发顺序

## 10.1 Phase 0：把骨架先搭起来

目标：先让仓库具备“能长项目”的形状。

交付物：

- 目录骨架
- 配置文件模板
- `target manifest` schema
- `prediction schema`
- `evaluation config`

验收标准：

- 任意新方法都能按统一目录接入
- 任意 target 都能被唯一标识

## 10.2 Phase 1：测试集构建

### Step 1. 拉原始数据

输入源：

- `BioLiP`
- `PDB`
- 可选：
  - `PLINDER`
  - `scPDB`
  - `PDBbind`

交付物：

- `data/raw/biolip/`
- `data/raw/pdb/`
- `data/raw/metadata/`

### Step 2. 生成候选 target 清单

做法：

- 先按 ligand 类型过滤
- 再按结构质量过滤
- 再按解析成功率过滤

交付物：

- `data/interim/candidate_targets.csv`

### Step 3. 做去冗余和泄漏检查

做法：

- sequence clustering
- 对照常见训练集 overlap 审计
- 时间切分过滤

交付物：

- `data/interim/nonredundant_targets.csv`
- `results/reports/leakage_audit_v1.md`

### Step 4. 生成 site labels

做法：

- ligand atom extraction
- binding residue assignment
- site center calculation

交付物：

- `data/processed/test_core_v1_targets.csv`
- `data/processed/test_core_v1_sites.json`

### Step 5. 手工抽查并冻结

交付物：

- `results/reports/manual_qc_v1.md`
- `data/splits/test_core_v1.txt`

## 10.3 Phase 2：方法复现

### Step 1. 先做 P0

顺序固定：

1. `fpocket`
2. `P2Rank`
3. `P2Rank/PRANK rescoring`
4. `DeepPocket`

每接一个方法都做：

- 环境文件
- 运行脚本
- 输出 parser
- smoke test

### Step 2. 再做 P1

- `PUResNetV2.0`
- `GrASP`
- `IF-SitePred`

### Step 3. 最后做 P2

- `VN-EGNN`

## 10.4 Phase 3：统一评测与报告

交付物：

- benchmark 主结果表
- per-method 明细
- runtime 报告
- failure analysis
- 数据集卡片
- 方法卡片

## 11. 建议的脚本与文件命名

为了让整个项目后面不乱，建议一开始就把命名约束住。

### 数据侧脚本

- `scripts/build_candidate_targets.py`
- `scripts/filter_biolip_ligands.py`
- `scripts/remove_redundancy.py`
- `scripts/audit_leakage.py`
- `scripts/build_site_labels.py`
- `scripts/freeze_test_split.py`

### 方法侧脚本

- `scripts/run_fpocket.py`
- `scripts/run_p2rank.py`
- `scripts/run_deeppocket.py`
- `scripts/run_puresnet.py`
- `scripts/run_grasp.py`
- `scripts/run_if_sitepred.py`
- `scripts/run_vnegnn.py`

### 评测脚本

- `scripts/normalize_predictions.py`
- `scripts/evaluate_predictions.py`
- `scripts/aggregate_results.py`
- `scripts/generate_report.py`

## 12. 8 周可执行里程碑

## Week 1

- 冻结项目范围
- 建仓库骨架
- 明确测试集标签 schema
- 列出 P0/P1/P2 方法清单

## Week 2

- 拉取 BioLiP/PDB 元数据
- 搭候选 target 生成脚本
- 出第一版候选集

## Week 3

- 做 sequence 去冗余
- 做时间切分和训练泄漏审计
- 输出第一版 `test_core_v1` 候选

## Week 4

- 生成 site labels
- 做人工抽查
- 冻结 `test_core_v1`

## Week 5

- 接入 `fpocket`
- 接入 `P2Rank`
- 建立统一 parser 与 smoke benchmark

## Week 6

- 接入 `DeepPocket`
- 接入 `P2Rank/PRANK rescoring`
- 输出第一版 baseline 结果

## Week 7

- 接入 `PUResNetV2.0` / `GrASP` / `IF-SitePred` 中至少 1-2 个
- 跑第一轮全量 benchmark

## Week 8

- 汇总结果
- 做 failure analysis
- 形成内部 benchmark v1 报告

## 13. 主要风险与规避

### 风险 1：测试集泄漏

这是最大风险。  
规避方式：

- 时间切分
- sequence overlap 审计
- 对照公开训练源

### 风险 2：配体标签脏

很多 PDB 配体并不是真正生物学相关。  
规避方式：

- 优先用 BioLiP
- 人工抽查
- 保留排除名单

### 风险 3：方法复现成本高于预期

尤其是深度学习方法。  
规避方式：

- 先做 P0
- 每个方法独立环境
- 先 smoke 后 full run

### 风险 4：不同方法输出不可比

规避方式：

- 一开始就定义统一 schema
- 主指标以 center-based metrics 为主

### 风险 5：项目过早追求“最全”

规避方式：

- v1 只求稳跑通
- 先 4 个方法，再扩到 7-8 个

## 14. 建议的第一版交付顺序

如果现在就开始干，最合理的顺序是：

1. 先把仓库骨架和 schema 定下来。
2. 先把 `test_core_v1` 做出来，不要先沉迷装方法。
3. 先接 `fpocket + P2Rank`，建立最小闭环。
4. 然后接 `DeepPocket` 和 `PRANK`。
5. 等评测闭环稳定后，再拉 `PUResNetV2.0 / GrASP / IF-SitePred`。

一句话概括：  
**这个项目的关键不是“模型越新越好”，而是“测试集可信、复现实验可重复、评测协议统一”。**

## 15. 参考依据

- BioLiP 提供 biologically relevant ligand-protein interaction 标注，适合作为 v1 数据源基础。
- 2024 年的 comparative evaluation 给出了当前主流方法集合，也强调了 benchmark 协议和冗余控制的重要性。
- PLINDER 强调了大规模、低泄漏 split、apo/predicted structure 扩展和真实泛化评测的重要性。
- P2Rank、DeepPocket、PUResNet、GrASP、IF-SitePred、VN-EGNN 可以作为当前方法谱系的代表。
