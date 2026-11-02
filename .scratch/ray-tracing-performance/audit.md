# 射线追踪性能审查 — 多声源场景的成本与正确性边界

方法：静态源码核查；未构建、未运行 Rhino 仿真、未做 profile。
范围：`PachydermAcoustic_Universal/`、`Hare/`、`PachydermAcoustic_Rhinoceros/`。
状态：审计已修订，代码尚未修改。实施要求见 [SPEC](spec.md)。

## 修订说明

保留原报告 #1–#10 编号，便于追溯。此次修订撤回“已满核所以声源并行无收益”“几何无分配”“RAD 只影响面板”“恢复注释即可安全并行”“可降低一个数量级”等未经证明或不成立的结论；纠正内存口径，并补充数组所有权、异步完成和收敛依赖。

下文区分三件事：**源码事实**、**待测性能影响**、**修改的正确性条件**。条目顺序不代表已测得的收益或瓶颈排序。源码位置均相对 `F:/PachydermAcousticFOrk/`，仅反映本轮核查时的代码。

## 结论先说

- 多声源会放大接收器构造、对象池重复创建、追踪和滤波成本；各部分占比尚未测量。
- Hybrid 与 Mapping 的声源计算循环都是串行，但这不能证明串行就是首要瓶颈，也不能证明并行化没有收益。
- 固定射线数模式按 `Environment.ProcessorCount` 投递任务，不等于 CPU 已满载；收敛模式按 `Source.Revolution_Period()` 投递每轮任务，基类默认值为 80，并非按核数投递。
- 应优先修复明确的正确性问题和重复工作，再根据分阶段测量决定是否改并行度、缓冲所有权或声学算法。

证据：`PachydermAcoustic_Rhinoceros/Pach_Hybrid_Control.cs:1188–1299`；`PachydermAcoustic_Rhinoceros/Pach_Mapping_Control.cs:799–876`；`PachydermAcoustic_Universal/Pachyderm_Acoustic_Universal/Classes_Complete_Trace - TaskAll.cs:218–247、297–333`；`PachydermAcoustic_Universal/Pachyderm_Acoustic_Universal/Classes_Sources.cs:166`。

## 成本放大结构

Hybrid 每个声源依次创建接收器、执行所选仿真阶段，再生成射线结果的滤波器；完成的接收器结果保留在声源结果数组中。Mapping 在声源计算循环之前构造各声源的映射接收器，不能统一描述成“每进入声源循环才分配银行”。

| 阶段 | 已确认的工作 | 尚待测量 |
|---|---|---|
| 接收器构造 | 分配能量、压力、方向和 RAD 数据；方向直方图有重复 Energy 分配 | 构造耗时、累计分配、存活量及峰值 |
| 追踪器构造／启动 | 重复创建静态对象池；启动追踪任务 | 冷启动与后续声源的成本差异 |
| 追踪 | 求交、分裂、衰减、接收器累加、同步和临时分配 | CPU、锁等待、分配与 GC 的占比 |
| 滤波 | 接收器逐个生成压力及方向滤波器 | FFT 槽位争用、耗时与临时内存 |
| UI 等待 | 内层轮询及合并后的同步等待，外层还有二次轮询 | 计算结束到结果可见的延迟 |

## #1 对象池初始化不幂等

**事实：** `BroadRayPool.Initialize()` 每次创建新池，默认分配 200,000 个 BroadRay；`OctaveRayPool.Initialize()` 也会重建池，默认容量 5,000。SplitRayTracer 构造调用两者，IS_Trace 启动还会调用 BroadRayPool 初始化。这里是对象池，不是线程池。两个 `initialized` 字段仅写未读。

**影响：** 重复构造与旧池的回收成本真实存在，但尚不能称为第一大瓶颈。

**修改边界：** 初始化应具备线程安全的一次创建／安全发布语义，`Initialize` 和 `Instance` 必须走同一个入口；一行空值检查不够，构造函数也不应提前发布半初始化实例。改变初始化不等于证明租借算法线程安全，更不等于可以开启声源级并行。

`GeometryPool` 会缓存池引用，但当前三个仓库未找到其消费者；“旧池被永久保留”只能写为该类型被访问后的条件性风险，不能断言当前路径已泄漏。

证据：`PachydermAcoustic_Universal/Pachyderm_Acoustic_Universal/Classes_Rays.cs:486–642`；`PachydermAcoustic_Universal/Pachyderm_Acoustic_Universal/Classes_Complete_Trace - TaskAll.cs:76–77`；`PachydermAcoustic_Universal/Pachyderm_Acoustic_Universal/Classes_IS_Trace.cs:67`。

## #2 Directional_Histogram 重复分配 Energy

**事实：** 两个方向直方图构造函数在基类已分配 Energy 后再次替换它。第一份数组立即失去该字段的引用，成为可回收垃圾；不是整个接收器的全部数据都分配两遍。

**建议：** 删除子类重复的 Energy 分配，保留基类初始化及全部方向数据初始化。验收覆盖两种构造方式和 8／24 个频带。

证据：`PachydermAcoustic_Universal/Pachyderm_Acoustic_Universal/Classes_Receivers.cs:1122–1152、1444–1499`。

### 接收器内存：累计分配不等于存活量

以下仅适用于含方向直方图的球形接收器，假设 44,100 Hz、1 秒、8 band；MB 使用十进制，暂不计对象头、数组引用、滤波结果和 FFT 临时缓存。

| 数据 | 构造累计分配 MB | 构造结束后可达的数据 MB |
|---|---:|---:|
| 基类 Energy，随后被替换 | 2.8224 | 0 |
| 子类最终持有的 Energy | 2.8224 | 2.8224 |
| Pressure | 2.8224 | 2.8224 |
| 正负三轴方向数据 | 16.9344 | 16.9344 |
| 三类 Reflection_Analysis_Data | 8.4672 | 8.4672 |
| **合计** | **33.8688** | **31.0464** |

24 band 分别为 **101.6064 MB 累计分配**、**93.1392 MB 可达数据**。修复 #2 可降低累计构造分配，不会把已经可回收的旧数组再从稳定存活量中扣一次。实际驻留、峰值和回收时间必须测量；“数十 GB”还取决于声源数、接收点数、时长和频带数，不能直接认定为第二大成因。

Mapping 的相应构造实参为 1,000 Hz，且存在不同接收器配置，不能套用上述 44.1 kHz 内存表。RAD 三类均被写入，详细收敛还会读取，不能为了省内存直接删维度。

证据：`PachydermAcoustic_Universal/Pachyderm_Acoustic_Universal/Classes_Receivers.cs:722–740、803、842、855–864`；`PachydermAcoustic_Rhinoceros/Pach_Mapping_Control.cs:799`。

## #3 Create_Filter 当前串行，但不能直接恢复注释中的并行版本

**事实：** Receiver_Bank 按接收器串行调用滤波；该方法声明为 `async void`，当前没有有效 await，另有未使用的 Task 数组和被注释的并行代码。

**正确性条件：** 注释中的 lambda 捕获 for 循环变量 `rec`；恢复 await 后，外层 Task.Run 也无法等待 async void 的真正完成。因此“取消注释即可并行”是错误建议。

**最小修复方向：** 批量 void 入口先保持串行、同步完成的契约，去掉多余 async 及死代码，保证完成或抛错后才返回，并在失败时复位 CreatingFilter。若将来确需异步并行，必须统一可等待契约、逐个迁移调用方，并解决索引捕获、进度线程归属、失败传播与峰值内存问题。

另外两个返回滤波数组的重载也设置 CreatingFilter，却未在返回前复位；完成状态修复需一并覆盖，保持它们的签名及返回值语义。

**性能限制：** ETCToFilter 内部调用 FFT 时使用固定槽位 0，底层按槽位加锁。接收器级并行是否能加速尚未证明，不应作为默认修复。

证据：`PachydermAcoustic_Universal/Pachyderm_Acoustic_Universal/Classes_Receivers.cs:574–610、1719`；`PachydermAcoustic_Rhinoceros/Pach_Hybrid_Control.cs:1285`；`PachydermAcoustic_Universal/Pachyderm_Acoustic_Universal/Pach_SP.cs:801、846、1815`。

## #4 BroadRayPool.pull 存在重复分配，但不能盲目原地清零或替换锁

**事实：** pull 为 Energy、Decimation 和空 Freq_Bands 重新分配数组。两个 new_BroadRay 重载随后都替换 Energy；其中一个创建完整频带表，另一个直接引用调用方的 `_Octaves`，不是两个重载都创建新频带表。

**所有权约束：** Energy 引用调用方给出的数组，Decimation 会由追踪器替换，Freq_Bands 可能与追踪器共享。原地清空这些数组可能修改不属于池的状态。要减少分配，必须先明确“谁拥有可写缓冲、活多久、何时可复用”，并覆盖 8／24 band 及父子射线重叠存活。

**租借约束：** 当前 BroadRayPool 通过环形索引选择对象，release 仅释放信号量、不携带对象身份；信号量数量本身不能证明任意归还顺序下不会再次出租仍在使用的对象。不能把 lock 简单换成 Interlocked 取模，就宣称整个对象池正确或更快。

本条的热路径重构需另行验证所有权和租约，不列入无条件的小修复；#1 的初始化修复不包含该重构。

证据：`PachydermAcoustic_Universal/Pachyderm_Acoustic_Universal/Classes_Rays.cs:524–599`；`PachydermAcoustic_Universal/Pachyderm_Acoustic_Universal/Classes_Complete_Trace - TaskAll.cs:357–368`。

## #5 DirPower.Clone 不能直接改成共享 SourcePower

**事实：** 基类 DirPower 为每次调用返回功率副本，Directions 将它交给 new_BroadRay，后者不做防御性复制。射线会屏蔽频带、吸收能量及分裂，原地修改 Energy。

**结论：** Clone 有当前所有权契约上的必要性。直接返回共享 SourcePower 会写坏声源及其他射线；按 threadid 复用也必须证明重叠活跃射线不会共享可写状态。基类另一处 Clone 属于功率访问方法，不能都计成每射线开销。

保留现有独立功率副本；如果 profile 证明值得改，再设计填充独占缓冲的契约，并覆盖方向性声源、线／面声源、8／24 band 和所有生成入口。

证据：`PachydermAcoustic_Universal/Pachyderm_Acoustic_Universal/Classes_Sources.cs:131–145、221–232、649–655`；`PachydermAcoustic_Universal/Pachyderm_Acoustic_Universal/Classes_Rays.cs:472–481、531、546`；`PachydermAcoustic_Universal/Pachyderm_Acoustic_Universal/Classes_Complete_Trace - TaskAll.cs:359`。

## #6 求交仍有分配，不能认定几何不是瓶颈

**事实：** Polygon_Scene.shoot 的普通命中分支创建三个 List，失败分支也创建三个，是互斥路径，不是失败再额外创建三个。X_Event_NH 分支复用事件已有列表，但不等于整次查询无分配。

Voxel_Grid 已使用 DDA、FirstTest 和最近命中剪枝，查询临时状态按线程保存；但 X_Event 是 class，Shoot 入口和成功路径有显式 new，场景包装层也预先 new 了事件。不能据此写“无分配、已证明线程安全、瓶颈不在这里”。

**修改边界：** 输出列表会逃逸给调用方，移入 ThreadLocal 缓冲前必须验证生命周期、重入及跨查询覆盖问题；场景也不只使用 Voxel_Grid。先测量，再决定是否改变返回契约，而不是直接复用输出容器。

证据：`PachydermAcoustic_Universal/Pachyderm_Acoustic_Universal/Classes_Scenes.cs:796–825`；`Hare/Voxel_Grid.cs:330–373`。

## #7 接收器衰减和锁是候选热点，线程副本不能只在结束时合并

**事实：** 宽带命中按选中频带计算衰减，完整频带时为 8／24 次 Math.Pow；直方图累加还计算压力并同步。Combine_Clones 已支持合并 Energy、Pressure、方向和三类 RAD，但 SplitRayTracer 当前未以线程银行副本方式使用它。

**正确性条件：** 追踪器每轮任务结束后立即检查收敛；检查器持有主银行的 Energy，详细模式还读取 RAD。完整副本若只在最终结束时合并，最小收敛模式可能一直看不到新能量；只拆 RAD 也会让详细模式读取旧数据。

将来若使用副本，必须在每轮检查前提供一致、累计且不重复计数的合并视图，不能简单把每轮的累计副本反复相加；也不能用共享写引用假装合并。接收器分块还可能增加追踪工作，不能只看到内存下降。

距离量化查表会引入误差；副本会放大内存。两者均需 profile、误差／内存预算和独立设计，不能承诺线性扩展。

证据：`PachydermAcoustic_Universal/Pachyderm_Acoustic_Universal/Classes_Receivers.cs:467–481、777–806、1170–1182`；`PachydermAcoustic_Universal/Pachyderm_Acoustic_Universal/Classes_Complete_Trace - TaskAll.cs:331、792、859–933、1097–1115`。

## #8 单频带 CheckRay 写 RAD 缺少同步，是正确性缺陷

**事实：** BandCount 为 1 时，三类 RAD 都存在未加锁的读改写；宽带和多分量分支则使用对应频带数组作为锁。应补齐一致的累加同步和有效样本边界，避免嵌套持有直方图锁与 RAD 锁。

**影响：** 不只是反射分析面板。详细收敛模式（Raycount = 0）的 SpecMin_Convergence_Check 读取 RAD，决定窗口及最低射线数；丢失更新可能改变停止时机和最终 IR。固定射线数模式及最小收敛模式（Raycount < 0）需区别对待。

验收应以串行正确累加／明确的收敛条件为参考，不以有竞争的旧多线程结果为金标准。修复本条无需先引入 #7 的高内存副本方案。

证据：`PachydermAcoustic_Universal/Pachyderm_Acoustic_Universal/Classes_Receivers.cs:832–864`；`PachydermAcoustic_Universal/Pachyderm_Acoustic_Universal/Classes_Complete_Trace - TaskAll.cs:182–183、859–933`。

## #9 Scatter_Early 会产生多个子射线，但收益与增长率尚未证明

**事实：** 每个非零能量的几何八度组都会生成散射分支，没有额外的正能量阈值；存在透射时，每次早期调用最多产生 8 × 3 = 24 条 OctaveRay 子射线。此处 8 个几何组不因第三倍频程直接变成 24 个组。

**限定：** 这些子射线后续进入 Scatter_Late，而非递归照搬整个宽带 Scatter_Early。Late 也可能分裂，因此整体数量取决于路径、系数、阈值与截止时间；不能由“单次最多 24 条”推出 24 的指数树或确定的指数增长。

能量阈值、概率抽样／俄罗斯轮盘等会改变采样及误差特性，需要单独的声学验证与兼容策略。“可降低一个数量级”没有实测依据，不属于本 SPEC 的实现范围。

证据：`PachydermAcoustic_Universal/Pachyderm_Acoustic_Universal/Classes_Materials.cs:2995–3070、3129–3156`；`PachydermAcoustic_Universal/Pachyderm_Acoustic_Universal/Classes_Complete_Trace - TaskAll.cs:416、422–427、473–475`。

## #10 Calculate_Conv 保留历史任务并每轮显式 GC

**事实：** 两种收敛模式都复用持续增长的 MainRays，每轮追加 Revolution_Period 个任务，再等待整个历史集合；每轮调用 GC.Collect。固定射线数模式不执行该循环。

**建议：** 先把本轮任务集合限制在本轮，等待全部完成后再检查收敛，保持种子、射线调度输入、累计数据和取消边界。不要遗留仍运行的任务或丢失异常。

显式 GC 是否移除，应在清理历史引用后比较 GC 暂停、总耗时与内存峰值；不能未经测量就承诺删除一定更快、更省内存。本条属于任务生命周期／内存管理，不属于改变散射模型来减少射线数。

证据：`PachydermAcoustic_Universal/Pachyderm_Acoustic_Universal/Classes_Complete_Trace - TaskAll.cs:213、297–336`。

## UI 等待：独立于射线算法的可修复成本

RunSimulation 首次等待使用 do 循环；needfeedback 为 false 时仍同步 Sleep。Combine_ThreadLocal_Results 后的第二个 do 循环无条件先 Sleep 3 秒，当前已完成的仿真也会付出这段等待；在 Rhino UI 上下文中会阻塞界面。

Hybrid／Mapping 的外层 GetAwaiter 加 Delay 轮询又引入最多约一个轮询间隔的残余等待，并非每源固定 1.5／3 秒；线程调度拥塞还可能额外延迟。应直接 await 现有 Task，内层先检查状态、确需等待时再异步延迟，保留结果合并和 finally 中的环境恢复。

证据：`PachydermAcoustic_Rhinoceros/Rhino_Pach_Tools.cs:37–109`；`PachydermAcoustic_Rhinoceros/Pach_Hybrid_Control.cs:1194–1195、1212–1213、1238–1240、1275–1277`；`PachydermAcoustic_Rhinoceros/Pach_Mapping_Control.cs:825–826、875–876`。

## 测量与实施顺序

1. 保存可重复的基线：固定场景、声源／接收点、时长、采样率、频带、模式、射线数或每轮样本规则；区分冷／暖池，记录平台及运行时。
2. 独立修复 #8 的 RAD 正确性，以及 #3 的同步完成／失败契约；这不是通过删声学工作换速度。
3. 实施 #1 的安全初始化、#2 的重复分配删除、#10 的本轮任务生命周期和 UI 等待修复，并分别验证结果与阶段指标。
4. 用 profile 决定是否另行推进 #4／#5／#6／#7／#9；不按未经测量的“最大收益”承诺排序。

现有日志可辅助观察，但不能独立完成上述诊断：`_ts` 从 Begin 起计算追踪时间，不包含此前接收器／池构造和之后滤波；活跃分布式路径中 `_rayTotal` 主要在 End_Trace 入口递增，途中逃逸的子射线不都计入，不能称完整生成数或全部求交工作量。

至少分别记录每源准备／构造、追踪、滤波、等待和端到端墙钟时间；嵌套阶段不得重复相加。分配量、可达内存、进程峰值和 GC 暂停分别报告。跨 await 的全流程不能拿某个线程的分配字节数冒充总分配。

证据：`PachydermAcoustic_Rhinoceros/Pach_Hybrid_Control.cs:1285–1289`；`PachydermAcoustic_Universal/Pachyderm_Acoustic_Universal/Classes_Complete_Trace - TaskAll.cs:221、242、440–445、491–494`。

## 本次没有验证的事项

- 没有构建、仿真、profile 或硬件测量；所有加速幅度、瓶颈排序、峰值和跨平台表现均未验证。
- 没有全面审查 Direct_Sound、ImageSourceData 和 IS_Trace 内部的优化空间；只核对相关入口、池调用和等待关系。
- 没有证明全部几何实现、对象池租约、FFT 或多个并发仿真之间的全局线程安全；不得由局部修复推导这些保证。
