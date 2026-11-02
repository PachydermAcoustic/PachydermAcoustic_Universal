# 08: Rhino 等待流程（D7／审计 UI 等待章节）

**What to build:** RunSimulation 保留现有 Task 返回签名，但启动后及结果合并后先检查现有状态协议，确需等待才异步延迟——两个阶段和 needfeedback 的两个取值都不使用阻塞式休眠。已完成状态不付出首次或合并后的固定延迟；内部仍在运行时保留 3,000 ms 轮询间隔，不引入新的配置或完成通知框架。正常完成路径保留一次 Combine_ThreadLocal_Results 及其顺序；合并启动的后续工作继续按状态协议等待；失败/取消不在仍有写入者时强行合并。Hybrid/Mapping 中对 RunSimulation 结果的消费者直接 await，去掉外层 GetAwaiter 加 Delay 的二次轮询。正常、取消及异常退出保留自动保存、重绘、进程优先级等现有 finally 恢复，UI 更新保持在正确线程，不新增取消 API。

**Blocked by:** 03（相同文件；直接 await Task 是滤波后台调用契约的一部分）

**Status:** done

- [x] RunSimulation 保留 Task 返回签名；启动后与合并后先检查状态协议，确需等待才异步延迟；两阶段 × needfeedback 两值均无阻塞式休眠
- [x] 已完成状态不付出首次或合并后的固定延迟（含现有合并前固定 3 秒）；仍在运行时保留 3,000 ms 轮询
- [x] 正常完成路径保留一次 Combine_ThreadLocal_Results 及其顺序；失败/取消不在仍有写入者的情况下强行合并，保留各仿真类型安全退出语义
- [x] Hybrid/Mapping 消费者直接 await RunSimulation 的 Task；删除外层 GetAwaiter+Delay 二次轮询；其后读取结果与启动下一阶段的顺序不变
- [x] 正常/取消/异常退出保留自动保存、重绘、进程优先级等 finally 恢复；UI 更新在正确线程；无新取消 API
- [ ] Rhino 集成清单执行并记录（headless 部分已完成并记录于 .rt-perf-notes/08-rhino-checklist.md 与 report.md 第五节；Rhino 会话内交互项未完成）：已完成仿真与无后续工作的合并路径不再额外等待 3 秒；确有合并后工作时正确等待；needfeedback 真/假时界面均可继续交互

----

**2026-10-26 (ticket 08 implementation):** Rhino commit `4ef9f72` (branch `ticket/08`, repo `PachydermAcoustic_Rhinoceros`); headless proof in `tests/RunSimulationWait/` (this repo, links the production `RunSimulation.cs` against Rhino stubs). Ran: Rhino plugin builds net7.0 + net48 (0 errors), `dotnet run --project tests/RunSimulationWait` (21/21 PASS, exit 0), shared harness `tests/AnalysisCancellation` still green. Proven headlessly: completed sims pay 0–4 ms (no first/post-combine fixed delay, both needfeedback values), 500 ms running work waits exactly one 3,000 ms poll in phase 1 and in the combine-started phase, combine runs exactly once, failure/cancel propagate before any combine, autosave restored on every exit. Still needs a live Rhino session (see `.rt-perf-notes/08-rhino-checklist.md`): UI interactivity during waits, Esc-cancel UX, progress-box/error-dialog behavior of the reworked filter call sites.

Note on 重绘 (redraw): nothing in the Rhino plugin disables `RedrawEnabled` anymore (repo-wide grep finds no reference), so the pre-spec finally's redraw restore was vestigial; RunSimulation's finally keeps the live restores (autosave, process priority) and does not add a dead redraw call.
