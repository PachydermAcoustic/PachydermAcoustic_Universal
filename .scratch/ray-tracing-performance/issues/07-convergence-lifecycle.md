# 07: 收敛轮任务生命周期（D6／审计 #10）

**What to build:** Calculate_Conv 的任务集合仅表示当前轮：不再持有和重复等待以前各轮的任务；当前轮任务全部结束后才检查收敛、开始下一轮或完成退出。保持每轮 Revolution_Period、随机种子输入规则、声学运算、累计结果和既有取消边界。任务失败不得被当成成功收敛，不得留下继续写接收器的未观察任务；失败路径结算本轮已启动工作并走现有失败清理。收敛检查仍在本轮全部任务完成后读取主银行的完整累计 Energy/RAD。本期不默认删除或重新配置每轮显式 GC——先消除历史任务引用并测量剩余 GC 成本。

**Blocked by:** 01（断言工程）, 04（收敛检查读取 RAD，需先修复累加正确性）

**Status:** done

- [x] 任务集合仅含当前轮任务；不保留历史轮任务；当前轮全部结束后才检查收敛/开下一轮/退出
- [x] 每轮 Revolution_Period、随机种子输入规则、声学运算、累计结果、取消边界保持
- [x] 任务失败不当作成功收敛；失败路径结算本轮已启动工作、走现有失败清理，无未观察任务继续写接收器
- [x] T7 断言通过：最小/详细两种模式分别跨越多轮，每轮检查可见上一轮累计数据；最后一轮所有贡献计入；有效能量夹具按原收敛条件结束；零能量/无命中用例按原规则处理
- [x] T8 断言通过：多轮后无按轮数持续增长的历史任务保留；取消及受控任务失败不提前宣布成功
- [x] 保留每轮显式 GC；记录清理历史引用后的 GC 成本观察供后续独立决策

## Comments

**实现（Classes_Complete_Trace - TaskAll.cs）：** Calculate_Conv 的轮循环内改为局部 List<ValueTask> RoundRays（容量 Revolution_Period），仅承载本轮任务；逐个 await 结算整轮（含故障任务，首个异常被捕获），随后置空引用，不再保留或重复等待历史轮任务。原字段 MainRays 删除；Begin 固定射线数路径改用局部 FixedRays（同一逻辑，一次性等待后 conclude），与收敛路径生命周期互不共享。每轮任务故障时：整轮先完全结算（无未观察任务继续写接收器），首个异常存入新增公共属性 SplitRayTracer.Failure，立即 break——不执行收敛检查（故障不会被当作成功收敛），不再启动后续轮，随后照常 Conclude_Simulation 与 _ts 记录。取消边界不变（轮内 cast_ray_distributed 首行检查 + 检查后 IsCancellationRequested break）。每轮 GC.Collect() 保留未动。

**T7/T8（tests/AnalysisCancellation/ConvergenceLifecycleSuite.cs）：** T7 用 5x4x3 m 鞋盒房间（Polygon_Scene，α=0.3、s=0.15，需显式 partition()），最小/详细两模式各跑至原收敛条件结束（修复后实测最小 ~35 轮、详细 ~120 轮）；断言轮数从 _currentRay / Revolution_Period 推导（详细模式 SpecMin 零数据早退不 Fill，Payloads 少计）、每轮 Fill 时刻银行累计严格递增、末轮贡献计入银行、最小模式最后三轮 Fill 载荷 conv1/conv2/convinf 均 ≤0.1 阈值。T8 断言零能量最小模式恰跑 25 轮后取消（不自动收敛）、银行恰为零；保留断言改为 round 8→25 窗口 <20 KB（修复前同窗口实测 +44.7/+45.5 KB，必定失败）；详细零能量模式用 3 s 定时取消（其检查零数据时在 Fill 前早退，无轮数钩子），断言停止且不宣告成功；真实能量第 3 轮取消断言恰 3x80 条射线（本轮全部结算）；FlakyScene（Empty_Scene 子类，List-shoot 第 51 次起抛 InvalidOperationException）注入故障：断言 Failure 为注入异常、达到 Stopped、无检查执行、返回后 150 ms 内无 shoot 活动/银行写入。固定射线数模式回归由既有 RayTraceDeterminismSuite 覆盖（Begin 局部化后仍 PASS）。

**GC 观察记录（供后续独立决策，未删除 GC.Collect）：** 修复后（历史任务引用已消除），鞋盒场景每轮显式 GC.Collect() 实测 12–15 ms，live 集合 ~74 MB（harness 全部套件共享进程），且 120+ 轮恒定不随轮数增长——成本由 live 集合主导，与每轮追踪分配无关。是否移除仍需按 spec 用内存峰值对照独立决定；本票只提供观察。
