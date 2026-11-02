// Fixed multi-source scenario with per-source stage timing (ticket 09, Deliverable B).
// Same scene/seeds as tests/PhaseTimingBaseline and tests/T6Baseline; parameterized by
// source count (1 or 3). Referenced by tests/PerfAcceptance/Program.cs.
using Pachyderm_Acoustic;
using Pachyderm_Acoustic.Environment;
using System.Diagnostics;
using Hare.Geometry;
using Material = Pachyderm_Acoustic.Environment.Material;

namespace RayTraceTests;

internal static class Perf
{
    internal static class Scenario
    {
        const int RECEIVERS = 4, SAMPLE_RATE = 44100, CUTOFF_MS = 500, BANDS = 8, RAY_COUNT = 2000;
        static readonly Point[] ReceiverPts = new Point[]
        {
            new Point(2, 4, 1.2), new Point(4, 3, 1.2), new Point(6, 4, 1.2), new Point(4, 1.5, 1.2),
        };

        /// <summary>Runs the fixed scenario; returns per-source stage lines.</summary>
        public static async Task<List<string>> Run(int sourceCount)
        {
            double[] swl = new double[8] { 100, 101, 102, 103, 104, 105, 106, 107 };
            Point[] srcPts = new Point[] { new Point(2, 2, 1.7), new Point(6, 2, 1.7), new Point(4, 5, 1.7) };
            var sources = new RandomSource[sourceCount];
            for (int i = 0; i < sourceCount; i++) sources[i] = new RandomSource(swl, srcPts[i], i, false);

            var scene = new Polygon_Scene(20, 50, 101325, 0, false, true);
            scene.Construct(
                new Point[][][] { new Point[][] { T6Shared.Face(0,0,0, 8,0,0, 8,6,0, 0,6,0), T6Shared.Face(0,0,4, 0,6,4, 8,6,4, 8,0,4),
                                                  T6Shared.Face(0,0,0, 0,0,4, 8,0,4, 8,0,0), T6Shared.Face(0,6,0, 8,6,0, 8,6,4, 0,6,4),
                                                  T6Shared.Face(0,0,0, 0,6,0, 0,6,4, 0,0,4), T6Shared.Face(8,0,0, 8,0,4, 8,6,4, 8,6,0) } },
                new List<Material> { new Basic_Material(new double[] { 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26 }) },
                new List<Scattering> { new Lambert_Scattering(new double[] { 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12 }) },
                new List<double[]> { new double[8] });
            var partitionPts = new List<Point>(ReceiverPts);
            foreach (var s in sources) s.AppendPts(ref partitionPts);
            scene.partition(partitionPts);

            var lines = new List<string>();
            int workers = Environment.ProcessorCount;
            for (int s = 0; s < sourceCount; s++)
            {
                var e2e = Stopwatch.StartNew();
                var sw = Stopwatch.StartNew();
                Receiver_Bank rec = new Receiver_Bank(ReceiverPts, sources[s], scene, SAMPLE_RATE, CUTOFF_MS, Receiver_Bank.Type.Stationary, false);
                sw.Stop(); var prepareRecv = sw.Elapsed;

                sw = Stopwatch.StartNew();
                Direct_Sound DS = new Direct_Sound(sources[s], rec, scene, Enumerable.Range(0, BANDS).ToArray());
                sw.Stop(); var dsCtor = sw.Elapsed;
                sw = Stopwatch.StartNew();
                await T6Shared.Run_To_Completion(DS);
                sw.Stop(); var dsWork = sw.Elapsed;

                sw = Stopwatch.StartNew();
                SplitRayTracer RT = new SplitRayTracer(sources[s], rec, scene, (CUTOFF_MS / 1000.0) * scene.Sound_speed(0), new int[] { 0, BANDS - 1 }, 0, RAY_COUNT, null);
                sw.Stop(); var rtCtor = sw.Elapsed;
                sw = Stopwatch.StartNew();
                var tasks = new List<Task>();
                for (int p = 0; p < workers; p++)
                {
                    int start = (int)Math.Floor((double)p * RAY_COUNT / workers);
                    int end = p == workers - 1 ? RAY_COUNT : (int)Math.Floor((double)(p + 1) * RAY_COUNT / workers);
                    Calc_Params cp = new Calc_Params(scene, start, end, p, T6Shared.BASE_SEED + p);
                    tasks.Add(Task.Run(async () => await RT.Calculate(cp)));
                }
                await Task.WhenAll(tasks);
                RT.Combine_ThreadLocal_Results();
                sw.Stop(); var rtWork = sw.Elapsed;

                sw = Stopwatch.StartNew();
                DS.Create_Filter();
                rec.Create_Filter(null);
                sw.Stop(); var filter = sw.Elapsed;
                e2e.Stop();

                lines.Add($"  src {s}: e2e {e2e.Elapsed.TotalSeconds:0.###} s | recv {prepareRecv.TotalSeconds:0.###} | DS ctor {dsCtor.TotalSeconds:0.###} work {dsWork.TotalSeconds:0.###} | RT ctor {rtCtor.TotalSeconds:0.###} work {rtWork.TotalSeconds:0.###} | filter {filter.TotalSeconds:0.###}");
            }
            return lines;
        }
    }
}
