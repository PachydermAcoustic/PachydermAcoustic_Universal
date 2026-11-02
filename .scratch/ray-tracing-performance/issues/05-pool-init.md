# 05: 对象池幂等安全初始化（D2／审计 #1）

**What to build:** BroadRayPool 与 OctaveRayPool 各只有一个线程安全的惰性创建入口（可用标准库一次初始化机制）：Initialize 仅确保实例存在，Instance 返回同一实例，构造期间不向外发布半初始化实例。重复 Initialize 不替换已有池、不隐式重置容量/位置/信号量/在用对象。保留默认容量 200,000/5,000 与现有租借、生成、归还入口签名；删除仅写未读的初始化标志。核查 SplitRayTracer、IS_Trace、声源生成和 GeometryPool 等消费者，保持调用兼容。不改环形出租与匿名归还算法，不开启多仿真实例并发。

**Blocked by:** 01（断言工程承载 T1/T2 断言）

**Status:** done

- [x] BroadRayPool 与 OctaveRayPool 各一个线程安全惰性创建入口；Initialize 仅确保实例存在；Instance 返回同一实例；构造期间不发布实例
- [x] 重复 Initialize 不替换池、不重置容量、位置、信号量或在用对象；持有射线时再次初始化不修改该射线
- [x] 保留默认容量与租借/生成/归还入口签名；不引入重置接口；删除仅写未读的初始化标志；SplitRayTracer、IS_Trace、声源生成、GeometryPool 消费者兼容
- [x] T1 断言通过：首次并发访问 Initialize/Instance 与反复调用得到相同、完整实例；全部借用配对归还；若同一射线在归还前再次被出租，报告为独立正确性阻塞，不暗改租约算法
- [x] T2 断言通过：至少两个不同功率谱的声源依次计算，冷/暖池结果正确；先前结果不被下一声源覆盖；8/24 band 均通过
- [x] 通过连续多源及冷/暖运行验证；若回归暴露在用射线被再次出租，暂停相关合并并单列阻塞

---

**Implementation note (2026-10-26, branch `ticket/05`):** both pools now create themselves through a single `Lazy<T>` (`LazyThreadSafetyMode.ExecutionAndPublication`); `Initialize()` only forces materialization and never replaces an existing pool; the constructors no longer self-publish (`instance = this` removed), so no half-constructed instance is ever visible. Write-only `initialized` flags and the unused `insLock` were deleted; capacities (200,000 / 5,000), ring-rent and anonymous-release algorithms, and all borrow/generate/release signatures unchanged. Consumers verified compatible: `SplitRayTracer` ctor and `IS_Trace` (their `Initialize()` calls are now idempotent no-ops when a pool exists), all `BroadRayPool.Instance.new_*` source/direction overloads, and `GeometryPool` — its type-init capture now always yields the live instance because the instance is never replaced. Assertions live in `tests/AnalysisCancellation/PoolInitSuite.cs` (T1 lifecycle + lease canary, T2 successive sources 8/24 bands cold/warm); red was confirmed against the pre-fix code (concurrent first access produced multiple instances). No in-use-ray re-rental observed at existing concurrency (8 workers × 250 pull/release cycles) — lease algorithm untouched as required.
