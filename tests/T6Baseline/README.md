# T6 fixed-count BEFORE/AFTER regression (ticket 09)

Deterministic fixed-count comparison: the SAME fixed scenario runs against the BASELINE
build and the FIXED build, and the dumped fixtures are compared element-wise.

Fixed per spec T6: scene (shoebox 8×6×4 m), worker count (`Environment.ProcessorCount`),
per-worker seeds (`0x20251026 + p`), disjoint ray ranges, source initial state, total ray
count (2000 rays/source). 3 sources × 4 receivers × 8 octave bands, 44100 Hz, 500 ms.

Dumped per source/receiver/band: RT Energy, DS Energy, filter Pressure (sum of squares),
all three RAD categories, direction histograms (+/− axis, sum of squared components),
sub-ray counts; plus source power (`Source.SWL()`) and the receiver-bank SWL/band table
after all rays (cross-ray corruption check).

Verdict rule: `|a−b| ≤ 1e-12 + 1e-9·max(|a|,|b|)`; any newly appearing NaN/Infinity fails;
integer tokens (sub-ray counts) must match exactly.

## Layout

- `T6Scenario.cs`, `T6Shared.cs` — the scenario + comparator. Compiled into BOTH
  `tests/T6Baseline` (standalone runner) and `tests/AnalysisCancellation`
  (`T6RegressionSuite`, one `Register()` line in `Program.cs`).
- The baseline worktree `rt-09-before` (detached at `b6490b5`) gets a verbatim copy of
  this project; because it compiles the same files against the pre-fix library, its
  fixture is the BEFORE side. Nothing is committed there.

## Reproduce (exact commands, Windows)

```
# AFTER side (fixed build, branch ticket/09)
cd F:/PachydermAcousticFOrk/.rt-worktrees/rt-09/tests/T6Baseline
dotnet run                       # records t6-fixture.txt (or compares if present)
dotnet run -- t6-fixture.txt strict   # same-build compare, RAD deltas fail

# BEFORE side (baseline commit b6490b5)
#   copy the four project files into the baseline worktree's tests/ and run there:
cd F:/PachydermAcousticFOrk/.rt-worktrees/rt-09-before/tests/T6Baseline
dotnet run                       # records the BEFORE t6-fixture.txt

# Cross-build compare (run in the AFTER worktree, referencing the BEFORE fixture):
cd F:/PachydermAcousticFOrk/.rt-worktrees/rt-09/tests/T6Baseline
dotnet run -- <path-to-before>/t6-fixture.txt        # RAD reported, not failed
dotnet run -- <path-to-before>/t6-fixture.txt strict # RAD deltas also fail

# Or through the shared harness:
T6_BEFORE_FIXTURE=<path-to-before>/t6-fixture.txt dotnet run --project tests/AnalysisCancellation -f net8.0
```

## Expected verdict (measured 2026-10-27, 28 workers)

- directions, source power, band table, sub-ray counts: **bitwise identical**;
  Energy within 1 ulp; Pressure within 2 ulp (FFT summation order).
- RAD categories differ on a minority of elements: the baseline RAD accumulation is
  racy (lost updates, spec D5). Two BEFORE runs differ from each other by up to
  ~20× the tolerance on RAD, so the baseline RAD values are not a valid reference
  (spec: “RAD 对照串行正确参考，不对照旧竞争输出”). The fixed build is stable
  run-to-run; its RAD is checked against a serial reference by the
  `rad accumulation` suite in `tests/AnalysisCancellation` instead.

If the compare fails on non-RAD elements, that is a real regression: localize it with
the printed key (`S<src> RT R<rec> B<band> E|P|DIR ...`).
