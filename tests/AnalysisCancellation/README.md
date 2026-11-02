# AnalysisCancellation — shared deterministic-assertion harness (ticket 01)

Plain console assertion app for the ray-tracing performance work (spec
`.scratch/ray-tracing-performance/spec.md` § Testing Decisions, matrix T1–T8).
Tickets 03–07 add their suites here.

Constraints (kept): outside the production source directory, accesses Universal
**via project reference only**, no Rhino dependency, no test packages, no
private reflection. Exit code `0` = all suites passed; `1` = at least one
failure. Each suite prints one `PASS <name>` or `FAIL <name>: <error>` line.

## Run

net8.0 (primary; the net7.0-targeted Universal library runs on the installed
net8 runtime):

```
dotnet run --project tests/AnalysisCancellation -f net8.0
```

net48:

```
dotnet build tests/AnalysisCancellation -f net48
tests/AnalysisCancellation/bin/Debug/net48/AnalysisCancellation.exe
```

Both commands verified on this machine: `PASS ray trace determinism (explicit
seed, fixed count)` / `ALL SUITES PASSED`, exit code 0.

**net7.0 blocker:** no net7.0 runtime is installed here (runtimes 8/9/10 only),
so a net7.0-flavored run of this harness is not executable on this machine and
is not claimed as passing. The referenced Universal library itself still builds
for `net7.0;net48`.

## Layout

- `Program.cs` — entry point; registers suites (one line per suite).
- `Harness.cs` — `Suite`/`Assert`/`RunAll` core. To add a suite (tickets 03–07):
  create `XSuite.cs` with `public static void Register()`, add one
  `XSuite.Register();` line to `Program.cs`. No other file changes.
- `Fixtures.cs` — shared deterministic fixtures (public production surface only):
  - `Fixtures.Scene()` — fixed `Empty_Scene` (no surfaces; each ray's fate is
    fixed by its seed alone).
  - `Fixtures.Source()/ThirdOctaveSource()` — `RandomSource` with fixed 8/24-band
    spectra; `dbOffset` gives a second distinct spectrum for T2.
  - `Fixtures.ReceiverBank()` — direct `Receiver_Bank` construction; `Fixtures.Receiver()`
    — direct single `Spherical_Receiver` construction.
  - `Fixtures.WorkerParams()` — fixed worker count, per-worker seed
    (`BaseSeed + worker`), disjoint ray ranges (T6 split).
  - `Fixtures.Tracer()/RunSeeded()` — fresh `SplitRayTracer` plus the
    explicit-seed, fixed-count entry `SplitRayTracer.Calculate(Calc_Params)`,
    awaited per worker sequentially so repeat runs are bit-identical.
  - `PoolFixtures` — pool identity accessors (`BroadRayPool.Instance`,
    `OctaveRayPool.Instance`, `GeometryPool.MainPool/OctPool`) for T1; they
    assert nothing themselves.
- `RayTraceDeterminismSuite.cs` — ticket 01 smoke: same explicit seed + fixed
  ray count ⇒ identical receiver energy histograms and RAD category 0, exact
  element-wise; `_currentRay` reaches the fixed total; recorded energy > 0 so
  the comparison is not vacuous.

## History note

The three `RCPachTools.RunSimulation` cancellation checks that previously lived
in `Program.cs` (commit b6490b5) were removed by ticket 01: they compiled Rhino
-repo source files (`PachydermAcoustic_Rhinoceros/RunSimulation.cs`) plus
partial-class hooks into this project, which is incompatible with the mandated
Rhino-free project-reference model (and did not resolve inside worktrees).
They remain recoverable from b6490b5 and belong with the Rhino-repo suites.
