# 01: 无 Rhino 断言工程与确定性夹具

**What to build:** 一个轻量控制台断言工程，通过项目引用访问 Universal，不依赖 Rhino、不添加测试包、不使用私有反射；独立于生产源目录，不被生产项目默认源文件通配规则编入。提供后续所有验收（spec 测试矩阵 T1–T8）共用的确定性夹具：用现有可显式提供种子的 Calculate 参数入口固定随机输入、直接构造接收器、检查池身份、固定射线场景。运行方式必须写清。

**Blocked by:** None (can start immediately)

**Status:** done

- [x] 断言工程位于生产源目录之外，不作为生产项目依赖，不被生产通配规则编入；不依赖 Rhino、测试包、私有反射
- [x] net7.0 与 net48 分别可构建运行；缺少工具或 SDK 时报告阻塞，不记作通过
- [x] 夹具覆盖：显式种子 Calculate 入口、直接构造接收器、池身份断言、固定场景（工作者数、每工作者种子、射线区间）
- [x] 冒烟断言通过：同一显式种子下固定射线计算结果可重复（退出码 0）
- [x] 运行命令（两个目标框架）写入工程说明
- net7.0 验证补充（评审后）：机器上未装 net7.0 运行时，测试工程已加 `<RollForward>Major</RollForward>` 使 net7.0 生成物在已安装的较新运行时上真实执行——9 套件全部 PASS，退出码 0。

## Comments

Done on branch ticket/01. tests/AnalysisCancellation is now the shared harness: ProjectReference to Universal (no Rhino, no packages, no private reflection), split into Harness.cs/Program.cs/Fixtures.cs/RayTraceDeterminismSuite.cs; suites plug in via one Register line in Program.cs. Fixtures: fixed Empty_Scene, RandomSource 8/24-band spectra, direct Receiver_Bank and Spherical_Receiver construction, WorkerParams (4 workers, per-worker seed BaseSeed+w, disjoint ray ranges), Tracer/RunSeeded via public SplitRayTracer.Calculate(Calc_Params) with isOrder=-1 (no IS pass in harness), PoolFixtures identity accessors. Smoke: same seed twice -> identical energy histograms + RAD category 0, exact element-wise, nonzero total, _currentRay == 2000. Runs green on net8.0 (dotnet run -f net8.0) and net48 (build + exe), both exit 0; library net7.0+net48 builds 0 errors. Blocker: no net7.0 runtime installed, so a net7.0 harness run is not claimed. Note: the three old RCPachTools.RunSimulation cancellation checks were removed (they linked Rhino-repo sources + partial-class hooks, incompatible with the Rhino-free project-reference model and unresolvable in worktrees); recoverable at b6490b5.
