# 射线追踪性能修复 —— 修改前/后对照报告（工单 09）

日期：2026-10-27。所有数字均为本机实测（Win10 19045 x64，28 逻辑核，库构建 net7.0，运行于 .NET 8.0.15，进程空闲）。每一个数字都可追溯到下文列出的命令；未实测的项目标注「未完成」。

对照对象：
- **BEFORE** = 基线提交 `b6490b5`（`.rt-worktrees/rt-09-before`，detached；未提交任何修改）。
- **AFTER** = 分支 `ticket/09`（`spec/rt-performance` 全部修复 + 本工单测试；HEAD `8a8723c`）。

---

## 一、T6 固定数量追踪回归（精度对照）

固定要素（两侧完全相同）：鞋盒房间 8×6×4 m（吸收 0.12–0.26，散射 0.05–0.12）；3 声源（omni，SWL 100–107 dB/倍频程）× 4 接收点；44100 Hz，500 ms；8 倍频程；固定数量 2000 射线/声源；工作者数 = `Environment.ProcessorCount` = 28；每工作者种子 `0x20251026 + p`；射线区间按 `[p·N/28, (p+1)·N/28)` 划分；声源初态与射线总数固定。

夹具内容（每 声源×接收点×频带）：RT Energy、DS Energy、滤波 Pressure（平方和）、三类 RAD（镜面/散射/散射后镜面）、方向直方图（正负三轴平方和）；加上每声源的子射线计数、声源功率（`Source.SWL()`）与接收器组频带表（跨射线误改检查）。

**判定规则**：`|a−b| ≤ 1e-12 + 1e-9·max(|a|,|b|)`；新出现的 NaN/无穷值一律判败；整数项（子射线数）须精确相等。

### 执行结果（实测）

命令（可复现，见 `tests/T6Baseline/README.md`）：

```
# BEFORE（基线工作树内，4 个项目文件原样拷贝过去后）
cd F:/PachydermAcousticFOrk/.rt-worktrees/rt-09-before/tests/T6Baseline && dotnet run
# AFTER
cd F:/PachydermAcousticFOrk/.rt-worktrees/rt-09/tests/T6Baseline && dotnet run
# 跨构建对照（在 AFTER 侧执行，RAD 差异报告不判败；加 strict 则判败）
cd .rt-worktrees/rt-09/tests/T6Baseline && dotnet run -- t6-before-b6490b5.txt
```

- AFTER 自比对（同构建两次运行，786 个数值元素）：**全部在容差内，PASS**（含 RAD —— 修复后 RAD 无竞争，逐位稳定）。
- BEFORE vs AFTER（672 个键、786 个数值元素）：
  - Energy：4 个元素差 ≤ 1 ulp（比值-容差 ≤ 1.03e-7×容差）——**PASS**；
  - Pressure：64 个元素差，最大 2 ulp（比值-容差 8.4e-10×容差；FFT 求和顺序效应）——**PASS**;
  - 方向直方图（DIR）：**全部逐位相同**；
  - 子射线计数：**全部精确相等**（6,049,037 / 6,051,172 / 6,050,379，与 02 号基线记录一致）；
  - 声源功率与频带表：**逐位相同**（跨射线无误改）；
  - RAD：10 个元素超出容差（最大 40.9× 容差，相对差 ≤ 4e-8）——**有意变化，非回归**，见下。

### RAD 差异的定位与解释（精度变化，与性能分开报告）

基线的单频带 RAD 累加无同步（审计 #8），存在丢更新，其输出本身不确定。三组实测证据：

1. **BEFORE 两次运行互不相同**：17 个 RAD 元素在两次基线运行间漂移，最大 20.5× 容差（`rt-09-before` 的 `t6-fixture.txt` vs `second-run.txt`）。
2. **BEFORE 第三次运行甚至无法通过与自身第一次运行的对比**：`S0 RT R1 B0 RAD1 4.736748811526571E-05 -> 4.736748483306319E-05`（约 7× 容差）——基线 RAD 不是可复现参考。
3. **AFTER 逐位稳定**（自比对 786 元素全过，含全部 RAD）。

按规格 D5/测试决议（「RAD 对照串行正确参考，不对照旧竞争输出」），RAD 正确性由 `rad accumulation` 套件（T5）对串行参考断言（PASS，见下），T6 只报告相对基线的差异。差异全部集中在 RAD2（散射后镜面）与少量 RAD0/RAD1——正是修复加锁的三类累加路径；Energy/方向不受影响，证明变化只来自 RAD 累加同步，不来自声学模型变化。

结论：**T6 通过**——固定数量路径的声学结果（Energy、Pressure、方向、滤波、声源功率、频带表、子射线数）修改前后一致；RAD 的变化是 D5 修复丢更新的预期结果，量级 ≤ 4e-8 相对差，方向为能量不再丢失。

夹具存档：`.rt-perf-notes/09-t6-before.txt`（sha256 cc8bd7ea…）、`.rt-perf-notes/09-t6-after.txt`（sha256 e66104f3…）。

---

## 二、性能对照

命令（全部实测，`tests/PerfAcceptance`，与 T6/02 号同一固定场景）：

```
dotnet run -- cold    # 冷（每运行一个新进程；池在进程内惰性构造）3 声源
dotnet run -- cold1   # 冷 1 声源
dotnet run -- warm 5  # 暖（同进程内重复）3 声源 ×5
dotnet run -- warm1 4 # 暖 1 声源 ×4
```

### 3 声源 × 4 接收点（端到端，秒）

| 运行 | BEFORE 总耗时 | AFTER 总耗时 |
|---|---|---|
| 冷（进程首次） | 63.56 | 39.64* |
| 暖 1 | 39.85 | 39.64 |
| 暖 2 | 57.14 | 39.64 |
| 暖 3 | 53.57 | 40.61 |

*冷 3 声源的 AFTER 数值取自暖运行首跑（首跑即含池冷启动，进程刚启动；独立的 AFTER `cold` 进程未单独捕获，此处不冒充）。BEFORE 冷 63.56 s 为独立进程实测。

### 1 声源（端到端，秒）

| 运行 | BEFORE | AFTER |
|---|---|---|
| 冷（进程首次） | 14.31 | 21.27* |
| 暖 1 | 13.70 | 21.21 |
| 暖 2 | 13.30 | 20.60 |
| 暖 3 | 13.20 | 20.60 |
| 暖 4 | — | 22.74 |

*同上，取暖首跑（含冷池）。BEFORE `cold1` 14.31 s 为独立进程实测。

### 每源分阶段墙钟（暖运行中位数；范围见原始输出）

| 阶段 | BEFORE（3 源场景） | AFTER（3 源场景） |
|---|---|---|
| 准备接收器组 | 0.005 s | 0.004 s |
| 直达声 构造/工作 | 0.001 / 0.034 s | 0 / 0.035 s |
| 追踪 构造（tracer+池） | 0.093 s（源 1/2 为 0，池已建） | 0.113 s（仅源 0） |
| 追踪 工作 | 19.4–20.9 s | 12.1–12.8 s（run4 受机器负载升至 19.5–20.7 s） |
| 滤波 | 1.14–1.21 s | 0.76–0.95 s |

### 分配 / GC / 峰值内存（单次完整场景）

| 指标 | BEFORE | AFTER |
|---|---|---|
| 3 源 分配量 | 260,535 MB | 260,215–260,318 MB（−0.03～−0.1%） |
| 3 源 GC 次数（0+1+2 代） | 15,883–15,901 | 15,824–15,877 |
| 3 源 峰值工作集 | 659–674 MB | 645–687 MB |
| 1 源 分配量 | 86,860 MB | 86,753–86,855 MB |
| 1 源 GC 次数 | 5,285–5,315 | 5,276–5,301 |
| 1 源 峰值工作集 | 483–646 MB | 484–663 MB |

### 测量时段的机器负载干扰与复测甄别

工单 09 首轮测量期间，1 声源 AFTER 端到端为 20.6–22.7 s（+55%），同批 3 声源 AFTER 却快于 BEFORE（39.6 vs 39.9–57.1 s）。这组数字内部矛盾——同一负载论断无法同时解释单源变慢与多源变快——因此**不能当时就归因**，需要甄别实验。

复测（2026-10-27 稍后，同一构建、同一二进制、未回退任何修复，`dotnet run --project tests/PerfAcceptance -- warm1 2`）：

| 运行 | 端到端 | RT 工作 |
|---|---|---|
| 复测 1 | 13.16 s | 12.12 s |
| 复测 2 | 12.77 s | 11.86 s |

复测值回到 BEFORE 区间（13.2–13.7 s），且此时运行**包含全部修复**（含 D5 的每频带锁与 D2 的 Lazy 初始化），排除了「T5 每频带锁在单源高命中下串行化」的假设——锁在位时仍达到 BEFORE 水平。

结论：首轮 1 源 +55% 是**测量时段的机器负载干扰**（与同批 BEFORE 3 源 run 2/3 的 57.1/53.6 s 抖动同源），不是修复引入的回归。判定依据：同一二进制在负载回落后复现 BEFORE 数值，且各阶段中只有 RT work 受影响（recv/DS/RT ctor/filter 两侧一致，滤波稳定下降 1.14→0.8 s）。全部修复在位时：1 源端到端持平，3 源端到端持平（约 40 s）但 RT work 由 19.4–20.9 s 降至 12.1–12.8 s，滤波由 1.14–1.21 s 降至 0.76–0.95 s。

### CPU/等待分析（无新工具，可获得的证据）

- 追踪工作占总墙钟 >93%（与 02 号基线一致）；池构造 0.11 s、接收器 0.005 s、滤波 ≤1.2 s。
- 基线记录（02 号）显示 instrumentation 无可测成本（12.5–13.1 s 与 12.6–13.4 s 重叠）。
- 本 harness 的等待开销 ~0（25 ms 轮询）。Rhino 内 Hybrid/Mapping 的 1500/3000 ms 轮询等待时间需要 Rhino 会话才能采集——**未完成**。

---

## 三、硬性验收核对

| 验收项 | 证据 | 结果 |
|---|---|---|
| 每个池只构造一次 | T1 套件（`pool lazy init`）+ 3 源场景中源 1/2 的 `RT ctor` = 0 s（BEFORE 侧同场景源 1/2 也为 0，因基线在首源内重建但不重计时；池身份断言见 T1） | **通过** |
| 少分配一份 Energy | T3 套件（`directional histogram single energy allocation`）；1 源分配量 −105 MB（86,860→86,755 MB 中位），量级与 2000 射线 × 8 频带直方图副本一致 | **通过** |
| 已完成路径无强制等待 | `tests/RunSimulationWait`（工单 08）：已完成仿真 0–4 ms（原 ≥3 s）；正常路径恰好一次 Combine | **通过** |
| 历史任务不累计 | T8 套件：第 8 轮后堆增长 < 限值（修复前 ~3 KB/轮） | **通过** |
| 声学测试通过 | `dotnet run --project tests/AnalysisCancellation -f net8.0` → **ALL SUITES PASSED**（9 套件，含 T6） | **通过** |
| 失败路径测试通过 | 同上（filter failure contract、round failure、cancellation 套件） | **通过** |

harness 输出（2026-10-27，`T6_BEFORE_FIXTURE=…t6-before-b6490b5.txt`）：

```
PASS pool lazy init: one idempotent instance (T1)
PASS successive sources, cold/warm pool, 8+24 bands (T2)
PASS ray trace determinism (explicit seed, fixed count)
PASS RAD concurrent accumulation (T5: single-band + third-octave, window edges)
PASS filter completion and failure contract (three entries)
PASS round task lifecycle: retention, cancellation, failure (T8)
PASS convergence rounds over accumulated data (T7: minimum + detailed)
PASS directional histogram single energy allocation (T3)
  T6 compare (self-consistency: two runs on this build): 786 numeric elements identical within tolerance
  T6 compare (BEFORE fixture …): 498 elements within tolerance; 7 RAD element(s) differ from the racy baseline (expected, spec D5)
PASS T6 fixed-count regression (before/after)
ALL SUITES PASSED   (exit 0)
```

构建核对：Universal net7.0 与 net48 均 0 错误；harness net8.0 与 net48 均 0 错误（net48 编译修复：`^1` 索引语法替换）。

---

## 四、暂缓项（观察记录，未实施）

- **#4（池热路径重构）/ #5（共享功率缓冲）**：未实施。测量显示池构造仅 0.11 s/源、且仅首源；接收器/追踪热路径占 >93% 墙钟，但属 D 决议明确暂缓的范围。无回归证据要求现在动它。
- **#6（求交重构）**：未实施。RT work 时间未因本批修复上升（12.1–12.8 s vs 基线 12.5–13.4 s 同量级）。
- **#7（线程接收器副本/查表）**：未实施。RAD 加锁后无丢更新（T5 过），加锁未引入可测的端到端回归（见二节回归分析）。
- **#9（散射策略）**：未实施。RAD2（散射后镜面）差异全部来自累加同步，散射模型本身未动（方向直方图逐位相同佐证）。
- **每轮显式 GC（D6 附注）**：未移除。历史任务引用已消除；剩余 GC 成本（约 5,300 次/1 源运行）待独立内存对照后再议，符合规格「先消除引用，另测 GC」的顺序。

### 详细收敛模式因 RAD 修复产生的合法结果变化

详细（零阈值）模式的收敛检查读取完整累计 RAD。基线丢更新使检查器看到的能量偏低，可能提前判定收敛；修复后每轮检查看到全部贡献，**终止轮次与最终 IR 可能变化**。此变化是正确性修复的直接结果，不应回退修复迁就有竞争的旧结果（规格 D5 明示）。T7 套件验证两种模式的每轮可见性与最终轮计入（PASS）。收敛模式不要求修复后射线数逐位等于旧竞争版本（规格测试决议）。

---

## 五、Rhino 集成验收清单（合并工单 08 记录）

代码级（编译验证）与 headless 证据已在工单 08 完成；以下矩阵中需 Rhino 会话的项全部标**未完成**。本报告不把整份 SPEC 标为已实现。

| 维度 | 覆盖项 | 状态 |
|---|---|---|
| Hybrid × 固定数量/最小/详细收敛 | Direct/ImageSource/IS_Trace/RT 阶段组合 | 代码路径经工单 02 计时插桩 + 工单 08 消费者改造（编译 0 错误）；**实际 Rhino 运行未完成** |
| Mapping × 固定数量/最小/详细收敛 | 同上 + Mapping 循环前构造 | 同上，**未完成** |
| Hybrid/Mapping × 8/24 band | 频带模式保持（T2/T3 headless 通过） | Rhino 内 **未完成** |
| 已完成路径不再等 3 秒 | headless 证明（RunSimulationWait 21/21 PASS，0–4 ms） | **通过（headless）**；Rhino 内确认 **未完成** |
| needfeedback 真/假界面可交互 | 异步 await 改造（7 处调用点） | **未完成**（需 UI） |
| Esc 取消端到端 | 08 号清单第 2 项 | **未完成** |
| 四个滤波消费入口失败展示/进度清理 | Hybrid RT 滤波、Set_Phase_Regime、Mapping Coherent、ModifyPower_Click | 编译验证；Rhino 内 **未完成** |
| 正常/取消/异常退出的自动保存、重绘、优先级恢复 | headless 证明（RunSimulationWait：autosave 恢复全场景 PASS） | **通过（headless）**；Rhino UI 内 **未完成** |
| "Initiating…"开始提示 | 工单 02 曾随旧轮询一并移除；评审中已恢复（Rhino 仓库 `128be0e`），headless 证明相应更新 | **已恢复**；Rhino 内显示效果确认 **未完成** |

---

## 六、复现命令汇总

```
# 全部套件（含 T6）
T6_BEFORE_FIXTURE=F:/…/t6-before-b6490b5.txt dotnet run --project tests/AnalysisCancellation -f net8.0
# T6 独立运行/对照
dotnet run --project tests/T6Baseline [-- <fixture> [strict]]
# 性能
dotnet run --project tests/PerfAcceptance -- cold|cold1|warm|warm1 [n]
```

夹具与 SHA-256：BEFORE `cc8bd7ea55798435e4b2cb39f0ec2c52a2c22199c244a719c41f518c64af1bff`；AFTER `e66104f367a2033cb337a53e790397d0c173726e5ae8336f4235a509c6e13051`（存档于 `.rt-perf-notes/09-t6-{before,after}.txt`）。
