using RayTraceTests;

// Shared deterministic-assertion harness (ticket 01). Register one line per suite here;
// suite logic lives in its own file so parallel tickets don't edit the same code.
// PoolInitSuite runs first so T1 exercises cold (first-ever) pool access.
PoolInitSuite.Register();
RayTraceDeterminismSuite.Register();
RadAccumulationSuite.Register();
FilterContractSuite.Register();
ConvergenceLifecycleSuite.Register();
HistogramEnergySuite.Register();

// T6 fixed-count BEFORE/AFTER regression (ticket 09): dumps the fixture twice and
// compares; also compares against the BEFORE fixture when T6_BEFORE_FIXTURE is set
// (see tests/T6Baseline/README.md).
T6RegressionSuite.Register();

return await Harness.RunAll();
