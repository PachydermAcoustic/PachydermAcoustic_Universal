namespace RayTraceTests;

/// <summary>
/// Tiny assertion harness shared by all suites (ticket 01). Plain console app:
/// each suite prints one PASS line on success, one FAIL line with the first
/// error otherwise; the process exit code is 0 only when every suite passed.
/// </summary>
internal static class Harness
{
    private static readonly List<(string Name, Func<Task> Body)> suites = new();
    private static int failures;

    /// <summary>Registers a suite for <see cref="RunAll"/>. Tickets 03–07 add suites here, one line each.</summary>
    public static void Suite(string name, Func<Task> body) => suites.Add((name, body));

    public static void Assert(bool condition, string message)
    {
        if (!condition) throw new Exception(message);
    }

    public static async Task<int> RunAll()
    {
        foreach (var (name, body) in suites)
        {
            try
            {
                await body();
                Console.WriteLine($"PASS {name}");
            }
            catch (Exception error)
            {
                failures++;
                Console.WriteLine($"FAIL {name}: {error.Message}");
            }
        }
        Console.WriteLine(failures == 0 ? "ALL SUITES PASSED" : $"{failures} SUITE(S) FAILED");
        return failures == 0 ? 0 : 1;
    }
}
