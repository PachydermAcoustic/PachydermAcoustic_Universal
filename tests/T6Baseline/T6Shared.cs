// Helpers shared verbatim between the T6 runner projects (ticket 09). Lives in the
// shared tests tree so the baseline worktree copy compiles the identical file.
using Hare.Geometry;
using Pachyderm_Acoustic;

internal static class T6Shared
{
    public const int BASE_SEED = 0x20251026;
    public const int RAY_COUNT = 2000;

    public static void Assert(bool condition, string message)
    {
        if (!condition) throw new Exception(message);
    }

    public static Point[] Face(params double[] xyz)
    {
        // 4 corners given as x,y,z triplets: c0,c1,c2 (cross product of c1-c0 and c2-c0 is the interior normal), c3.
        return new Point[] { new Point(xyz[0], xyz[1], xyz[2]), new Point(xyz[3], xyz[4], xyz[5]), new Point(xyz[6], xyz[7], xyz[8]), new Point(xyz[9], xyz[10], xyz[11]) };
    }

    public static async Task Run_To_Completion(Simulation_Type sim)
    {
        // Replicates RCPachTools.RunSimulation coordination without Rhino.
        sim.Begin();
        while (sim.ThreadState() != System.Threading.ThreadState.Stopped) await Task.Delay(25);
        sim.Combine_ThreadLocal_Results();
        while (sim.ThreadState() == System.Threading.ThreadState.Running) await Task.Delay(25);
    }
}
