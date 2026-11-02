# 02: 每源阶段计时基线（D1）

**What to build:** Hybrid 与 Mapping 主流程记录每个声源的 准备/构造（区分接收器与追踪器/池成本）、追踪、滤波、等待和端到端墙钟时间，形成可重复基线。单调计时；说明阶段边界与嵌套关系；未启用的阶段标为未执行而不是伪造零耗时。Mapping 循环前的构造成本归属到相应声源或明确列为批次准备。沿用现有日志渠道与开关，不新增日志服务、配置或依赖。

**Blocked by:** None (can start immediately)

**Status:** done

- [x] Hybrid 与 Mapping 每源摘要含准备/构造、追踪、滤波、等待、端到端墙钟；准备/构造区分接收器与追踪器/池成本
- [x] Mapping 循环前构造成本按声源归属或明确列为批次准备，无遗漏
- [x] 单调时钟计时；嵌套阶段不重复相加；工作时长与包含该工作的 await 总时间分别记录
- [x] 未启用阶段标记为未执行；不新增逐射线计数/日志热点
- [x] 保留现有公开计数字段含义；不修改功率谱、共享频带表或声学抽样规则（结果不变）
- [x] 产出一份修改前的分阶段基线记录（固定场景），供 09 号工单对照

## Comments

Implemented on branch `ticket/02` (Universal: 384db2d, Rhino: d8e4344).

- `Simulation_Type._workClock` (additive Stopwatch) is restarted by `RCPachTools.RunSimulation` around Begin()..Combine_ThreadLocal_Results; Rhino-side `Source_Stage_Timing` prints one per-source summary per source via the existing `RhinoApp.WriteLine` channel (prepare receivers, direct/image/is-trace/ray-trace with construct-vs-work-vs-await-wait split, filter span, end-to-end wall; stages that did not run are reported as not run).
- Hybrid and Mapping pre-loop scene construction/partitioning get an explicit batch-prep line; Mapping pre-loop map-mesh+receiver construction is attributed per source and labelled nested (do not add to batch total).
- Headless BEFORE baseline (fixed scenario, 4 warm runs) + AFTER comparison: `.rt-perf-notes/02-baseline.md`; fixture `.rt-perf-notes/02-fixture-before.txt` (sha256 1f0e5fab...c70eb), reproduced by `dotnet run --project tests/PhaseTimingBaseline`. Acoustic fixture identical within spec tolerance => results unchanged.
- Rhino-session numbers (UI wait, image source, IS-trace, Hybrid/Mapping batch prep in Rhino) remain incomplete: instrumentation prints them, but no Rhino run has been captured yet.
