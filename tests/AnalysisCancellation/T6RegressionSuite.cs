using RayTraceTests;

/// <summary>
/// T6 fixed-count BEFORE/AFTER regression (ticket 09, spec Testing Decisions T6).
/// Runs the shared fixed-count scenario from <c>tests/T6Baseline</c> (same compiled code
/// the baseline worktree runs via its own copy of that project), then:
///  - always: dumps the fixture twice in-process and compares the two runs element-wise
///    (tolerance 1e-12 + 1e-9*max, no new NaN/Infinity), proving the fixed-count path is
///    repeatable on this build — including RAD, which on the fixed build is race-free;
///  - when <c>T6_BEFORE_FIXTURE</c> points at the fixture recorded from baseline commit
///    b6490b5 (see tests/T6Baseline/README.md): compares BEFORE vs AFTER. Energy, Pressure,
///    direction, source power and band-table elements must match within tolerance; RAD
///    elements are reported, not failed: the baseline RAD accumulation is racy (spec D5)
///    and its output is not a valid reference (spec: “RAD 对照串行正确参考，不对照旧竞争输出”).
/// </summary>
internal static class T6RegressionSuite
{
    public static void Register() => Harness.Suite("T6 fixed-count regression (before/after)", Run);

    private static async Task Run()
    {
        string first = await T6Scenario.Run();
        string second = await T6Scenario.Run();
        string a = Path.Combine(Path.GetTempPath(), $"t6-a-{Guid.NewGuid():N}.txt");
        string b = Path.Combine(Path.GetTempPath(), $"t6-b-{Guid.NewGuid():N}.txt");
        File.WriteAllText(a, first);
        File.WriteAllText(b, second);
        T6Compare.Files(a, b, "self-consistency: two runs on this build");

        string before = Environment.GetEnvironmentVariable("T6_BEFORE_FIXTURE");
        if (before != null)
        {
            Harness.Assert(File.Exists(before), $"T6_BEFORE_FIXTURE points to a missing file: {before}");
            T6Compare.Files(before, a, $"BEFORE fixture {Path.GetFileName(before)} (RAD reported, not failed: baseline RAD is the racy pre-D5 output)", radReportOnly: true);
        }
        else
        {
            Console.WriteLine("  T6: T6_BEFORE_FIXTURE not set; BEFORE/AFTER comparison not executed in this run (see tests/T6Baseline/README.md)");
        }
    }
}
