//'Pachyderm-Acoustic: Geometrical Acoustics for Rhinoceros (GPL)
//'
//'This file is part of Pachyderm-Acoustic.
//'
//'Pachyderm-Acoustic is free software; you can redistribute it and/or modify
//'it under the terms of the GNU General Public License as published
//'by the Free Software Foundation; either version 3 of the License, or
//'(at your option) any later version.
//'
//'Pachyderm-Acoustic is distributed in the hope that it will be useful,
//'but WITHOUT ANY WARRANTY; without even the implied warranty of
//'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//'GNU General Public License for more details.
//'
//'You should have received a copy of the GNU General Public License
//'along with Pachyderm-Acoustic; if not, write to the Free Software
//'Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using Hare.Geometry;
using Polygon_Scene = Pachyderm_Acoustic.Environment.Polygon_Scene;

namespace Pachyderm_Acoustic
{
    namespace Numeric
    {
        namespace TimeDomain
        {
            public partial class Acoustic_Compact_FDTD
            {
                /// <summary>
                /// Fitted (cut) boundary cells (Bilbao, Hamilton, Botts and Savioja, TASLP 2016, Section III-C): each cell keeps only
                /// its fluid volume and open face areas, and walls act over their true area. Polygon_Scene only. The scheme then runs
                /// at cT/h = Cut_Cell_Courant, with every cell satisfying eq. (38).
                /// </summary>
                public static bool Use_Cut_Cells = false;

                /// <summary>cT/h with cut cells, at most 1/sqrt3. Cells are merged until sum sigma / nu &lt;= 4 / (cT/h)^2.</summary>
                public static double Cut_Cell_Courant = 0.4;

                /// <summary>Gives cell pieces on the far side of a wall their own unknowns instead of merging them into a neighbour.</summary>
                public static bool Cut_Cell_Separate_Regions = true;

                /// <summary>Gaps between polygons narrower than twice this (metres) are treated as closed; at most a quarter of a node spacing.</summary>
                public static double Cut_Cell_Seam_Tolerance = 0.01;

                /// <summary>What the cut-cell pass did, or null if it did not run.</summary>
                public Cut_Cell_Summary Cut_Cells { get; private set; }

                /// <summary>True when the cut-cell pass will run for this scene (it needs the polygons of a Polygon_Scene).</summary>
                private bool Cut_Cells_Apply { get { return Use_Cut_Cells && Rm is Polygon_Scene; } }

                private double Courant_Time_Step(double dy_, double c)
                {
                    return Cut_Cells_Apply ? 2 * dy_ * Math.Min(Cut_Cell_Courant, 1 / Math.Sqrt(3)) / c : dy_ * Math.Sqrt(2) / c;
                }

                /// <summary>Interior coupling (cT/h)^2 / 2 for the time step in use.</summary>
                private void Set_Interior_Coupling(double c)
                {
                    double r = c * dt / (2 * dy);
                    RDD_Node.Coupling = 0.5 * r * r;
                }

                private void Fit_Boundary_Cells()
                {
                    Cut_Cells = null;
                    if (!Cut_Cells_Apply) return;
                    Cut_Cell_Builder builder = new Cut_Cell_Builder(this, (Polygon_Scene)Rm);
                    Cut_Cells = builder.Build();
                }

                public sealed class Cut_Cell_Summary
                {
                    /// <summary>Cells cut by a single plane, treated exactly.</summary>
                    public int ExactCells;
                    /// <summary>Cells at edges, corners or polygon borders, whose pieces, volumes and face areas are found by sampling.</summary>
                    public int SampledCells;
                    /// <summary>Cell pieces merged into a neighbour: pieces on the far side of the wall from their node, and cells too small for eq. (38).</summary>
                    public int MergedPieces;
                    /// <summary>Volume (in cells) of pieces with no neighbour to merge into, left out.</summary>
                    public double DroppedVolume;
                    /// <summary>Largest sum of open face fractions over volume fraction after merging.</summary>
                    public double MaxStiffness;
                    /// <summary>Wall area given to cut cells, counting each side of a wall (physical units).</summary>
                    public double WallArea;
                    /// <summary>Volume (in cells) of pieces on the far side of a wall from their node, merged into a neighbour.</summary>
                    public double FarSideVolume;
                    /// <summary>Merges made for eq. (38), and the volume (in cells) of the groups merged.</summary>
                    public int StiffnessMerges;
                    public double StiffnessMergedVolume;
                    /// <summary>Far-side pieces kept as cells of their own (Cut_Cell_Separate_Regions), and their volume in cells.</summary>
                    public int SeparateRegions;
                    public double SeparateRegionVolume;
                    /// <summary>The stiffness limit for the time step in use, 4 / (cT/h)^2.</summary>
                    public double StiffnessLimit;
                    /// <summary>Wall area given to cut cells per material, counting each side of a wall (physical units).</summary>
                    public Dictionary<Environment.Material, double> WallAreaByMaterial = new Dictionary<Environment.Material, double>();
                    /// <summary>Wall area found in cells but given to no piece (should be zero).</summary>
                    public double UnassignedWallArea;
                    /// <summary>Time taken by the cut-cell pass.</summary>
                    public double Seconds;
                    /// <summary>Time per phase (s).</summary>
                    public string Phases;
                    public override string ToString()
                    {
                        return $"{ExactCells} plane-cut cells, {SampledCells} sampled cells, {MergedPieces} pieces merged ({FarSideVolume:F1} cells far-side, {StiffnessMerges} for stiffness, {StiffnessMergedVolume:F1} cells), {SeparateRegions} far-side regions kept ({SeparateRegionVolume:F1} cells), {DroppedVolume:F3} cells of volume dropped, max stiffness {MaxStiffness:F3} (limit {StiffnessLimit:F2}), wall area {WallArea:F1}, {Seconds:F1} s";
                    }
                }

                /// <summary>
                /// Builds the cut cells after the staircase mesher. Works in ideal lattice coordinates, in units of dy with x scaled
                /// by sqrt2 dy / dx, where neighbours are at distance 2 and each cell is a regular rhombic dodecahedron.
                /// </summary>
                private sealed class Cut_Cell_Builder
                {
                    struct V3
                    {
                        public double X, Y, Z;
                        public V3(double x, double y, double z) { X = x; Y = y; Z = z; }
                        public static V3 operator +(V3 a, V3 b) { return new V3(a.X + b.X, a.Y + b.Y, a.Z + b.Z); }
                        public static V3 operator -(V3 a, V3 b) { return new V3(a.X - b.X, a.Y - b.Y, a.Z - b.Z); }
                        public static V3 operator *(double s, V3 a) { return new V3(s * a.X, s * a.Y, s * a.Z); }
                        public double Dot(V3 b) { return X * b.X + Y * b.Y + Z * b.Z; }
                        public V3 Cross(V3 b) { return new V3(Y * b.Z - Z * b.Y, Z * b.X - X * b.Z, X * b.Y - Y * b.X); }
                        public double Length { get { return Math.Sqrt(X * X + Y * Y + Z * Z); } }
                    }

                    static readonly double R2 = Math.Sqrt(2);
                    const double FaceTolerance = 1e-7;   // ideal units (a node spacing is 2)
                    // Link k of an RDD node, as RDD_Node.Link_Nodes numbers them; Opp[k] is the link back.
                    static readonly V3[] E =
                    {
                        new V3(R2, 1, 1), new V3(R2, -1, -1), new V3(R2, -1, 1), new V3(R2, 1, -1),
                        new V3(-R2, 1, 1), new V3(-R2, -1, -1), new V3(-R2, -1, 1), new V3(-R2, 1, -1),
                        new V3(0, 2, 0), new V3(0, -2, 0), new V3(0, 0, 2), new V3(0, 0, -2)
                    };
                    static readonly int[] Opp = { 5, 4, 7, 6, 1, 0, 3, 2, 9, 8, 11, 10 };
                    static readonly V3[][] Faces;       // cell faces about the origin, outward-oriented
                    static readonly double FaceArea, CellVolume;

                    static Cut_Cell_Builder()
                    {
                        // Vertices of the cell { q : E_k . q <= 2 } from triple plane intersections.
                        List<V3> verts = new List<V3>();
                        for (int i = 0; i < 12; i++)
                            for (int j = i + 1; j < 12; j++)
                                for (int k = j + 1; k < 12; k++)
                                {
                                    if (!Solve3(E[i], E[j], E[k], 2, 2, 2, out V3 v)) continue;
                                    bool inside = true;
                                    for (int m = 0; m < 12 && inside; m++) if (E[m].Dot(v) > 2 + 1e-9) inside = false;
                                    if (inside && !verts.Any(u => (u - v).Length < 1e-9)) verts.Add(v);
                                }
                        Faces = new V3[12][];
                        for (int k = 0; k < 12; k++)
                        {
                            List<V3> on = verts.Where(v => Math.Abs(E[k].Dot(v) - 2) < 1e-9).ToList();
                            V3 c = new V3(on.Average(v => v.X), on.Average(v => v.Y), on.Average(v => v.Z));
                            V3 n = (1.0 / E[k].Length) * E[k];
                            V3 u = on[0] - c; u = (1.0 / u.Length) * u;
                            V3 w = n.Cross(u);
                            Faces[k] = on.OrderBy(v => Math.Atan2(w.Dot(v - c), u.Dot(v - c))).ToArray();
                        }
                        FaceArea = VectorArea(Faces[0]).Length;
                        double vol = 0;
                        for (int k = 0; k < 12; k++) vol += Centre(Faces[k]).Dot(VectorArea(Faces[k]));
                        CellVolume = vol / 3;

                        double step = 2 * R2 / SampleCount;
                        List<V3> sampleList = new List<V3>();
                        SampleIndex = new int[SampleCount, SampleCount, SampleCount];
                        for (int i = 0; i < SampleCount; i++)
                            for (int j = 0; j < SampleCount; j++)
                                for (int k = 0; k < SampleCount; k++)
                                {
                                    V3 g = new V3(-R2 + (i + 0.5) * step, -R2 + (j + 0.5) * step, -R2 + (k + 0.5) * step);
                                    bool ok = true;
                                    for (int m = 0; m < 12 && ok; m++) if (E[m].Dot(g) > 2) ok = false;
                                    SampleIndex[i, j, k] = ok ? sampleList.Count : -1;
                                    if (ok) sampleList.Add(g);
                                }
                        Samples = sampleList.ToArray();
                        SamplesByRadius = Enumerable.Range(0, Samples.Length).OrderBy(i => Samples[i].Length).ToArray();
                        FaceSamples = new V3[12][];
                        FaceSampleNearest = new int[12][][];
                        for (int f = 0; f < 12; f++)
                        {
                            V3 v0 = Faces[f][0], e1 = Faces[f][1] - v0, e2 = Faces[f][3] - v0;
                            FaceSamples[f] = new V3[FaceSampleCount * FaceSampleCount];
                            FaceSampleNearest[f] = new int[FaceSamples[f].Length][];
                            for (int a = 0; a < FaceSampleCount; a++)
                                for (int b = 0; b < FaceSampleCount; b++)
                                {
                                    V3 q = v0 + ((a + 0.5) / FaceSampleCount) * e1 + ((b + 0.5) / FaceSampleCount) * e2;
                                    int idx = a * FaceSampleCount + b;
                                    FaceSamples[f][idx] = q;
                                    FaceSampleNearest[f][idx] = Enumerable.Range(0, Samples.Length).OrderBy(i => (Samples[i] - q).Length).Take(12).ToArray();
                                }
                        }
                    }

                    static bool Solve3(V3 a, V3 b, V3 c, double ra, double rb, double rc, out V3 x)
                    {
                        double det = a.Dot(b.Cross(c));
                        x = new V3();
                        if (Math.Abs(det) < 1e-12) return false;
                        x = (1.0 / det) * (ra * b.Cross(c) + rb * c.Cross(a) + rc * a.Cross(b));
                        return true;
                    }

                    static V3 VectorArea(IList<V3> p)
                    {
                        V3 a = new V3();
                        for (int i = 0; i < p.Count; i++) a = a + p[i].Cross(p[(i + 1) % p.Count]);
                        return 0.5 * a;
                    }

                    static V3 Centre(IList<V3> p)
                    {
                        V3 c = new V3();
                        foreach (V3 v in p) c = c + v;
                        return (1.0 / p.Count) * c;
                    }

                    /// <summary>Sutherland-Hodgman: the part of convex polygon p with n . q &gt;= d.</summary>
                    static List<V3> Clip(IList<V3> p, V3 n, double d)
                    {
                        List<V3> o = new List<V3>(p.Count + 2);
                        for (int i = 0; i < p.Count; i++)
                        {
                            V3 a = p[i], b = p[(i + 1) % p.Count];
                            double fa = n.Dot(a) - d, fb = n.Dot(b) - d;
                            if (fa >= 0) o.Add(a);
                            if ((fa > 0 && fb < 0) || (fa < 0 && fb > 0)) o.Add(a + (fa / (fa - fb)) * (b - a));
                        }
                        return o;
                    }

                    sealed class Tri
                    {
                        public V3 A, B, C, N;     // ideal coordinates; N unit normal
                        public double D;          // N . q = D on the plane
                        public V3 M0, M1, M2;     // in-plane unit normals of edges AB, BC, CA, pointing into the triangle
                        public Environment.Material Mat;
                    }

                    sealed class Hit
                    {
                        public int Tri;
                        public double IdealArea, PhysicalArea;
                        public V3 Centroid;                    // of the triangle's part inside the cell
                        public bool SharedFace;                // lies in a face shared with a neighbour that counts its area
                    }

                    enum Kind { Untouched, Exact, Sampled }

                    // Cells cut by more than one plane are measured on SampleCount^3 interior points and FaceSampleCount^2 points per face.
                    const int SampleCount = 24, FaceSampleCount = 8;
                    static readonly V3[] Samples;
                    static readonly int[,,] SampleIndex;
                    static readonly int[] SamplesByRadius;
                    /// <summary>Index offsets within three sample spacings, nearest first.</summary>
                    static readonly (int di, int dj, int dk)[] ShellOffsets =
                        (from i in Enumerable.Range(-3, 7) from j in Enumerable.Range(-3, 7) from k in Enumerable.Range(-3, 7)
                         where i * i + j * j + k * k <= 9
                         orderby i * i + j * j + k * k
                         select (i, j, k)).ToArray();
                    static readonly V3[][] FaceSamples;
                    static readonly int[][][] FaceSampleNearest;

                    sealed class Cell
                    {
                        public int X, Y, Z;
                        public long Key;
                        public V3 C;
                        public List<Hit> Hits = new List<Hit>();
                        public Kind Kind = Kind.Untouched;
                        public V3 N; public double D;          // plane, for exact cells
                        public double NuPlus;                  // volume fraction on N . q >= D
                        public double[] SigmaPlus;             // open fraction of each face on that side
                        public int NodeSide;                   // +1 or -1
                        public Dictionary<Environment.Material, double> Wall; // physical wall area per material (exact cells)
                        public int[] Frag = { -1, -1, -1 };     // fragment index by side + 1 (0 = whole cell)
                        public int[] SampleFragment;           // sampled cells: fragment of each sample, -1 where none
                        public Pieces Pending;                 // sampled cells: pieces found, not yet numbered
                    }

                    /// <summary>The pieces of a sampled cell, numbered from 0 within the cell.</summary>
                    sealed class Pieces
                    {
                        public int[] OfSample;
                        public double[] Nu;
                        public bool[] Primary;
                        public Dictionary<Environment.Material, double>[] Wall;
                    }

                    sealed class Fragment
                    {
                        public Cell Cell;
                        public int Side;                       // 0 whole, +1 / -1 for the sides of an exact cell, 2 sampled piece
                        public double Nu;
                        public bool Primary;                   // contains the cell's node
                        public Dictionary<Environment.Material, double> Wall = new Dictionary<Environment.Material, double>();
                    }

                    readonly Acoustic_Compact_FDTD fvm;
                    readonly Polygon_Scene scene;
                    readonly Node[][][] F;
                    readonly Point min;
                    readonly double dx, dy, dz, c, fs, lambda, w0, v0Physical, Stiffness_Limit;
                    readonly long yStride, zStride;
                    readonly Dictionary<long, Cell> cells = new Dictionary<long, Cell>();
                    readonly List<Fragment> frags = new List<Fragment>();
                    readonly List<(int a, int b, double sigma)> couplings = new List<(int, int, double)>();
                    List<Tri> tris;
                    /// <summary>Reach of a polygon beyond its edges, in ideal units.</summary>
                    readonly double seam;

                    public Cut_Cell_Builder(Acoustic_Compact_FDTD fvm, Polygon_Scene scene)
                    {
                        this.fvm = fvm; this.scene = scene;
                        F = fvm.PFrame; min = fvm.Bounds.Min_PT;
                        dx = fvm.dx; dy = fvm.dy; dz = fvm.dz;
                        c = scene.Sound_speed(0);
                        fs = 1.0 / fvm.dt;
                        lambda = c * fvm.dt / dx;
                        w0 = RDD_Node.Coupling;
                        // Eq. (38): sum of open face fractions over volume fraction at most 4 / (cT/h)^2.
                        double r = c * fvm.dt / (2 * dy);
                        Stiffness_Limit = 4 / (r * r);
                        v0Physical = 4 * dx * dy * dz;
                        seam = Math.Max(0, Math.Min(Cut_Cell_Seam_Tolerance, 0.25 * dy)) / dy;
                        int ymax = 0, zmax = 0;
                        foreach (Node[][] a in F) { ymax = Math.Max(ymax, a.Length); foreach (Node[] b in a) zmax = Math.Max(zmax, b.Length); }
                        zStride = zmax + 1; yStride = (ymax + 1) * zStride;
                    }

                    long KeyOf(int x, int y, int z) { return x * yStride + y * zStride + z; }

                    V3 Ideal(Point p) { return new V3((p.x - min.x) * R2 / dx, (p.y - min.y) / dy, (p.z - min.z) / dz); }

                    static V3 CentreOf(int x, int y, int z)
                    {
                        int mod = x % 2;
                        return new V3((x - 0.5) * R2, 2 * (y + 0.5 - 0.5 * mod), 2 * (z + 0.5 - 0.5 * mod));
                    }

                    double PhysicalArea(IList<V3> p)
                    {
                        V3 a = new V3();
                        for (int i = 0; i < p.Count; i++)
                        {
                            V3 u = new V3(p[i].X * dx / R2, p[i].Y * dy, p[i].Z * dz);
                            V3 v = new V3(p[(i + 1) % p.Count].X * dx / R2, p[(i + 1) % p.Count].Y * dy, p[(i + 1) % p.Count].Z * dz);
                            a = a + u.Cross(v);
                        }
                        return 0.5 * a.Length;
                    }

                    bool InGrid(int x, int y, int z) { return x >= 0 && x < F.Length && y >= 0 && y < F[x].Length && z >= 0 && z < F[x][y].Length; }

                    /// <summary>The node at the far end of link k, as RDD_Node.Link_Nodes connects it.</summary>
                    bool Neighbour(int x, int y, int z, int k, out int nx, out int ny, out int nz)
                    {
                        int mod = x % 2;
                        nx = x; ny = y; nz = z;
                        if (k < 8)
                        {
                            nx = k < 4 ? x + 1 : x - 1;
                            int kk = k % 4;
                            if (mod == 0) { ny = y + (kk == 0 || kk == 3 ? 1 : 0); nz = z + (kk == 0 || kk == 2 ? 1 : 0); }
                            else { ny = y - (kk == 1 || kk == 2 ? 1 : 0); nz = z - (kk == 1 || kk == 3 ? 1 : 0); }
                        }
                        else if (k == 8) ny = y + 1;
                        else if (k == 9) ny = y - 1;
                        else if (k == 10) nz = z + 1;
                        else nz = z - 1;
                        return InGrid(nx, ny, nz);
                    }

                    Cell GetCell(int x, int y, int z)
                    {
                        long key = KeyOf(x, y, z);
                        if (!cells.TryGetValue(key, out Cell cell))
                        {
                            cell = new Cell { X = x, Y = y, Z = z, Key = key, C = CentreOf(x, y, z) };
                            cells.Add(key, cell);
                        }
                        return cell;
                    }

                    public Cut_Cell_Summary Build()
                    {
                        System.Diagnostics.Stopwatch clock = System.Diagnostics.Stopwatch.StartNew();
                        Cut_Cell_Summary summary = new Cut_Cell_Summary { StiffnessLimit = Stiffness_Limit };
                        List<string> laps = new List<string>(); double lastLap = 0;
                        Action<string> lap = name => { double now = clock.Elapsed.TotalSeconds; laps.Add($"{name} {now - lastLap:F1}"); lastLap = now; };
                        CollectTriangles();
                        FindCutCells();
                        lap("triangles+cells");
                        // Fit each material's admittance once, in parallel.
                        System.Threading.Tasks.Parallel.ForEach(tris.Select(t => t.Mat).Where(m => m != null).Distinct().ToList(),
                            m => new DIF_IWB_2p(m, fs, 0, lambda));
                        lap("material fits");
                        List<Cell> hit = cells.Values.ToList();
                        System.Threading.Tasks.Parallel.ForEach(hit, Classify);
                        lap("classify");
                        // Find the pieces of sampled cells in parallel; number them in order below.
                        System.Threading.Tasks.Parallel.ForEach(hit.Where(cl => cl.Kind == Kind.Sampled).ToList(), cl => cl.Pending = ComputePieces(cl));
                        lap("sampled pieces");
                        foreach (Cell cell in cells.Values)
                        {
                            if (cell.Kind == Kind.Exact) summary.ExactCells++;
                            else if (cell.Kind == Kind.Sampled) summary.SampledCells++;
                        }
                        List<Cell> geometry = cells.Values.Where(cl => cl.Kind != Kind.Untouched).ToList();
                        if (geometry.Count == 0) { summary.Seconds = clock.Elapsed.TotalSeconds; return summary; }

                        // Affected: cells with geometry and their neighbours. Everything else stays as the mesher built it.
                        HashSet<long> affected = new HashSet<long>();
                        foreach (Cell cl in geometry)
                        {
                            affected.Add(cl.Key);
                            for (int k = 0; k < 12; k++) if (Neighbour(cl.X, cl.Y, cl.Z, k, out int nx, out int ny, out int nz)) { GetCell(nx, ny, nz); affected.Add(KeyOf(nx, ny, nz)); }
                        }
                        foreach (Cell cl in cells.Values.ToList()) MakeFragments(cl);
                        List<long> keys = affected.ToList();
                        // Every cell a coupling can reach exists before the couplings are found in parallel.
                        foreach (long key in keys)
                        {
                            Cell a = cells[key];
                            for (int k = 0; k < 12; k++)
                            {
                                if (!Neighbour(a.X, a.Y, a.Z, k, out int nx, out int ny, out int nz)) continue;
                                Cell b = GetCell(nx, ny, nz);
                                if (b.Frag[1] < 0 && b.Kind == Kind.Untouched) MakeFragments(b);
                            }
                        }
                        lap("fragments");
                        AssignOrphanWalls();
                        lap("orphan walls");
                        List<(int a, int b, double sigma)>[] found = new List<(int, int, double)>[keys.Count];
                        System.Threading.Tasks.Parallel.For(0, keys.Count, i =>
                        {
                            List<(int, int, double)> local = new List<(int, int, double)>();
                            long key = keys[i];
                            Cell a = cells[key];
                            for (int k = 0; k < 12; k++)
                            {
                                if (!Neighbour(a.X, a.Y, a.Z, k, out int nx, out int ny, out int nz)) continue;
                                long kb = KeyOf(nx, ny, nz);
                                if (affected.Contains(kb) && kb < key) continue;
                                Couple(a, cells[kb], k, local);
                            }
                            found[i] = local;
                        });
                        foreach (List<(int, int, double)> l in found) couplings.AddRange(l);
                        lap("couplings");
                        Group(summary, affected);
                        lap("group+nodes");
                        summary.Phases = string.Join(", ", laps) + " [group: " + GroupLaps + "]";
                        summary.WallArea = wallTotal;
                        summary.WallAreaByMaterial = wallByMaterial;
                        summary.UnassignedWallArea = unassigned;
                        summary.Seconds = clock.Elapsed.TotalSeconds;
                        return summary;
                    }


                    void CollectTriangles()
                    {
                        tris = new List<Tri>();
                        for (int i = 0; i < scene.Count(); i++)
                        {
                            Point[] p = scene.polygon_unrounded(i);
                            if (p == null || p.Length < 3) continue;
                            Environment.Material mat = scene.AbsorptionValue[i];
                            for (int j = 1; j + 1 < p.Length; j++)
                            {
                                Tri t = new Tri { A = Ideal(p[0]), B = Ideal(p[j]), C = Ideal(p[j + 1]), Mat = mat };
                                V3 n = (t.B - t.A).Cross(t.C - t.A);
                                double len = n.Length;
                                if (!(len > 1e-12)) continue;
                                t.N = (1.0 / len) * n; t.D = t.N.Dot(t.A);
                                t.M0 = Inward(t.N, t.A, t.B, t.C); t.M1 = Inward(t.N, t.B, t.C, t.A); t.M2 = Inward(t.N, t.C, t.A, t.B);
                                tris.Add(t);
                            }
                        }
                    }

                    /// <summary>Unit vector in the plane, perpendicular to edge a-b, pointing towards c.</summary>
                    static V3 Inward(V3 n, V3 a, V3 b, V3 c)
                    {
                        V3 m = n.Cross(b - a);
                        double len = m.Length;
                        if (!(len > 0)) return new V3();
                        m = (1.0 / len) * m;
                        return m.Dot(c - a) >= 0 ? m : -1.0 * m;
                    }

                    void FindCutCells()
                    {
                        ConcurrentDictionary<long, ConcurrentBag<Hit>> found = new ConcurrentDictionary<long, ConcurrentBag<Hit>>();
                        double reach = R2 + 1e-9;     // the cell's circumradius
                        System.Threading.Tasks.Parallel.For(0, tris.Count, ti =>
                        {
                            Tri t = tris[ti];
                            double x0 = Math.Min(t.A.X, Math.Min(t.B.X, t.C.X)) - reach, x1 = Math.Max(t.A.X, Math.Max(t.B.X, t.C.X)) + reach;
                            double y0 = Math.Min(t.A.Y, Math.Min(t.B.Y, t.C.Y)) - reach, y1 = Math.Max(t.A.Y, Math.Max(t.B.Y, t.C.Y)) + reach;
                            double z0 = Math.Min(t.A.Z, Math.Min(t.B.Z, t.C.Z)) - reach, z1 = Math.Max(t.A.Z, Math.Max(t.B.Z, t.C.Z)) + reach;
                            int xa = Math.Max(0, (int)Math.Floor(x0 / R2 + 0.5)), xb = Math.Min(F.Length - 1, (int)Math.Ceiling(x1 / R2 + 0.5));
                            V3[] poly = { t.A, t.B, t.C };
                            for (int x = xa; x <= xb; x++)
                            {
                                int mod = x % 2;
                                int ya = Math.Max(0, (int)Math.Floor(y0 / 2 - 0.5 + 0.5 * mod)), yb = Math.Min(F[x].Length - 1, (int)Math.Ceiling(y1 / 2 - 0.5 + 0.5 * mod));
                                for (int y = ya; y <= yb; y++)
                                {
                                    int za = Math.Max(0, (int)Math.Floor(z0 / 2 - 0.5 + 0.5 * mod)), zb = Math.Min(F[x][y].Length - 1, (int)Math.Ceiling(z1 / 2 - 0.5 + 0.5 * mod));
                                    for (int z = za; z <= zb; z++)
                                    {
                                        V3 cc = CentreOf(x, y, z);
                                        if (Math.Abs(t.N.Dot(cc) - t.D) > reach) continue;
                                        // A little tolerance, so a polygon lying in a cell face is found by the cells on both sides.
                                        List<V3> clipped = new List<V3>(poly);
                                        for (int k = 0; k < 12 && clipped.Count >= 3; k++) clipped = Clip(clipped, -1.0 * E[k], -2 - FaceTolerance - E[k].Dot(cc));
                                        if (clipped.Count < 3) continue;
                                        double area = VectorArea(clipped).Length;
                                        if (!(area > 1e-9 * FaceArea)) continue;
                                        // In a face both cells see it; its area is counted once, by the cell with the smaller key
                                        // (the other side's fluid gets it as the neighbour across that face).
                                        bool shared = false;
                                        for (int k = 0; k < 12 && !shared; k++)
                                        {
                                            if (!clipped.All(v => Math.Abs(E[k].Dot(v - cc) - 2) < 10 * FaceTolerance)) continue;
                                            if (Neighbour(x, y, z, k, out int nx, out int ny, out int nz) && KeyOf(nx, ny, nz) < KeyOf(x, y, z)) shared = true;
                                        }
                                        found.GetOrAdd(KeyOf(x, y, z), _ => new ConcurrentBag<Hit>()).Add(new Hit { Tri = ti, IdealArea = area, PhysicalArea = PhysicalArea(clipped), Centroid = Centre(clipped), SharedFace = shared });
                                    }
                                }
                            }
                        });
                        foreach (KeyValuePair<long, ConcurrentBag<Hit>> kv in found.OrderBy(kv => kv.Key))
                        {
                            long key = kv.Key;
                            int x = (int)(key / yStride), y = (int)(key % yStride / zStride), z = (int)(key % zStride);
                            Cell cell = GetCell(x, y, z);
                            cell.Hits.AddRange(kv.Value.OrderBy(h => h.Tri));
                        }
                    }

                    /// <summary>
                    /// Exact when every polygon in the cell lies in one plane and together they cover the plane's whole section
                    /// of the cell (no polygon edge ends inside it); otherwise sampled.
                    /// </summary>
                    void Classify(Cell cell)
                    {
                        if (cell.Hits.Count == 0) return;
                        double total = cell.Hits.Sum(h => h.IdealArea);
                        if (total < 1e-6 * FaceArea) { cell.Hits.Clear(); return; }
                        Tri t0 = tris[cell.Hits[0].Tri];
                        bool planar = true;
                        foreach (Hit h in cell.Hits)
                        {
                            Tri t = tris[h.Tri];
                            if (Math.Abs(Math.Abs(t.N.Dot(t0.N)) - 1) > 1e-9) { planar = false; break; }
                            foreach (V3 v in new[] { t.A, t.B, t.C }) if (Math.Abs(t0.N.Dot(v) - t0.D) > 1e-7) planar = false;
                            if (!planar) break;
                        }
                        cell.Kind = Kind.Sampled;
                        if (!planar) return;

                        V3 n = t0.N; double d = t0.D;
                        double vol = 0; V3 sum = new V3();
                        double[] sigma = new double[12];
                        for (int k = 0; k < 12; k++)
                        {
                            V3[] face = Faces[k].Select(v => v + cell.C).ToArray();
                            List<V3> part = Clip(face, n, d);
                            if (part.Count < 3) continue;
                            V3 a = VectorArea(part);
                            sigma[k] = Math.Min(1.0, a.Length / FaceArea);
                            sum = sum + a;
                            vol += Centre(part).Dot(a);
                        }
                        double cap = sum.Length;     // the plane's section closes the surface
                        if (Math.Abs(total - cap) > 1e-6 * Math.Max(cap, FaceArea)) return;   // a polygon ends inside the cell, or overlaps
                        vol = (vol - d * cap) / 3;  // cap: outward normal -n, on which q . (-n) = -d
                        cell.Kind = Kind.Exact;
                        cell.N = n; cell.D = d;
                        cell.NuPlus = Math.Max(0, Math.Min(1, vol / CellVolume));
                        cell.SigmaPlus = sigma;
                        cell.NodeSide = n.Dot(cell.C) - d >= 0 ? 1 : -1;
                        cell.Wall = new Dictionary<Environment.Material, double>();
                        foreach (Hit h in cell.Hits)
                        {
                            if (h.SharedFace) continue;
                            Environment.Material m = tris[h.Tri].Mat;
                            cell.Wall.TryGetValue(m, out double a0);
                            cell.Wall[m] = a0 + h.PhysicalArea;
                        }
                    }

                    void MakeFragments(Cell cell)
                    {
                        if (cell.Frag[0] >= 0 || cell.Frag[1] >= 0 || cell.Frag[2] >= 0 || cell.SampleFragment != null) return;
                        if (cell.Kind == Kind.Exact)
                        {
                            foreach (int side in new[] { 1, -1 })
                            {
                                double nu = side > 0 ? cell.NuPlus : 1 - cell.NuPlus;
                                if (nu < 1e-9)
                                {
                                    // No fluid on this side within the cell (the wall lies in one of its faces): the wall's
                                    // area on this side goes to the neighbour across that face.
                                    foreach (Hit h in cell.Hits)
                                        if (!h.SharedFace) orphans.Add((cell, h.Centroid + (0.05 * side) * cell.N, tris[h.Tri].Mat, h.PhysicalArea));
                                    continue;
                                }
                                cell.Frag[side + 1] = frags.Count;
                                // Each side of the wall is an absorbing face for the fluid on that side.
                                frags.Add(new Fragment { Cell = cell, Side = side, Nu = nu, Primary = side == cell.NodeSide, Wall = new Dictionary<Environment.Material, double>(cell.Wall) });
                            }
                        }
                        else if (cell.Kind == Kind.Sampled) SamplePieces(cell);
                        else
                        {
                            cell.Frag[1] = frags.Count;
                            frags.Add(new Fragment { Cell = cell, Side = 0, Nu = 1, Primary = true });
                        }
                    }

                    /// <summary>
                    /// Pieces of a cell cut by several polygons: samples joined where the segment between neighbouring samples
                    /// crosses no polygon. Each piece's volume is its share of the samples; each polygon's part in the cell is
                    /// a wall of the pieces on either side of it.
                    /// </summary>
                    void SamplePieces(Cell cell)
                    {
                        Pieces p = cell.Pending ?? ComputePieces(cell);
                        cell.Pending = null;
                        int first = frags.Count;
                        for (int j = 0; j < p.Nu.Length; j++)
                            frags.Add(new Fragment { Cell = cell, Side = 2, Nu = p.Nu[j], Primary = p.Primary[j], Wall = p.Wall[j] });
                        cell.SampleFragment = p.OfSample.Select(j => first + j).ToArray();
                    }

                    /// <summary>Finds the pieces of a sampled cell; touches no shared state.</summary>
                    Pieces ComputePieces(Cell cell)
                    {
                        int ns = Samples.Length;
                        int[] root = Enumerable.Range(0, ns).ToArray();
                        int R(int i) { while (root[i] != i) { root[i] = root[root[i]]; i = root[i]; } return i; }
                        for (int i = 0; i < SampleCount; i++)
                            for (int j = 0; j < SampleCount; j++)
                                for (int k = 0; k < SampleCount; k++)
                                {
                                    int a = SampleIndex[i, j, k];
                                    if (a < 0) continue;
                                    foreach ((int di, int dj, int dk) in new[] { (1, 0, 0), (0, 1, 0), (0, 0, 1) })
                                    {
                                        if (i + di >= SampleCount || j + dj >= SampleCount || k + dk >= SampleCount) continue;
                                        int b = SampleIndex[i + di, j + dj, k + dk];
                                        if (b < 0) continue;
                                        if (Crosses(cell.C + Samples[a], cell.C + Samples[b], cell.Hits)) continue;
                                        int ra = R(a), rb = R(b);
                                        if (ra != rb) root[ra] = rb;
                                    }
                                }
                        int nodeRoot = -1;
                        foreach (int i in SamplesByRadius)
                            if (!Crosses(cell.C, cell.C + Samples[i], cell.Hits)) { nodeRoot = R(i); break; }
                        Dictionary<int, int> index = new Dictionary<int, int>();
                        List<int> count = new List<int>();
                        List<bool> primary = new List<bool>();
                        int[] of = new int[ns];
                        for (int i = 0; i < ns; i++)
                        {
                            int r = R(i);
                            if (!index.TryGetValue(r, out int j)) { index[r] = j = count.Count; count.Add(0); primary.Add(r == nodeRoot); }
                            count[j]++;
                            of[i] = j;
                        }
                        Pieces p = new Pieces { OfSample = of, Nu = count.Select(c0 => (double)c0 / ns).ToArray(), Primary = primary.ToArray() };
                        p.Wall = Enumerable.Range(0, count.Count).Select(_ => new Dictionary<Environment.Material, double>()).ToArray();
                        foreach (Hit h in cell.Hits)
                        {
                            if (h.SharedFace) continue;
                            Tri t = tris[h.Tri];
                            foreach (int side in new[] { 1, -1 })
                            {
                                V3 q = h.Centroid + (0.05 * side) * t.N;
                                int f = -1;
                                foreach (int i in NearestSamples(cell.C, q))
                                    if (!Crosses(q, cell.C + Samples[i], cell.Hits)) { f = of[i]; break; }
                                if (f < 0)
                                {
                                    // The fluid on this side is a sliver too thin to hold a sample: its wall goes to the
                                    // neighbouring cell that fluid opens into (see AssignOrphanWalls).
                                    orphans.Add((cell, q, t.Mat, h.PhysicalArea));
                                    continue;
                                }
                                p.Wall[f].TryGetValue(t.Mat, out double a0);
                                p.Wall[f][t.Mat] = a0 + h.PhysicalArea;
                            }
                        }
                        return p;
                    }

                    /// <summary>
                    /// The fragment of a cell that contains point q (ideal coordinates), or -1. For a sampled cell: the piece of
                    /// the nearest sample q reaches without crossing a polygon (candidates: the given nearest list, else all).
                    /// </summary>
                    /// <summary>Samples in order of distance from q (relative to the cell's centre c): the nearby ones first via the sample lattice, then all the rest.</summary>
                    static IEnumerable<int> NearestSamples(V3 c, V3 q)
                    {
                        double step = 2 * R2 / SampleCount;
                        V3 l = q - c;
                        int i0 = Math.Max(0, Math.Min(SampleCount - 1, (int)Math.Floor((l.X + R2) / step)));
                        int j0 = Math.Max(0, Math.Min(SampleCount - 1, (int)Math.Floor((l.Y + R2) / step)));
                        int k0 = Math.Max(0, Math.Min(SampleCount - 1, (int)Math.Floor((l.Z + R2) / step)));
                        foreach ((int di, int dj, int dk) in ShellOffsets)
                        {
                            int i = i0 + di, j = j0 + dj, k = k0 + dk;
                            if (i < 0 || j < 0 || k < 0 || i >= SampleCount || j >= SampleCount || k >= SampleCount) continue;
                            int idx = SampleIndex[i, j, k];
                            if (idx >= 0) yield return idx;
                        }
                        foreach (int idx in Enumerable.Range(0, Samples.Length).OrderBy(n => (Samples[n] - l).Length)) yield return idx;
                    }

                    int PieceAt(Cell cell, V3 q, int[] nearest)
                    {
                        if (cell.Kind == Kind.Untouched) return cell.Frag[1];
                        if (cell.Kind == Kind.Exact) return cell.Frag[SideOf(cell, q, cell.C) + 1];
                        IEnumerable<int> order = nearest != null ? (IEnumerable<int>)nearest : NearestSamples(cell.C, q);
                        foreach (int i in order)
                            if (!Crosses(q, cell.C + Samples[i], cell.Hits)) return cell.SampleFragment[i];
                        return -1;
                    }

                    int SideOf(Cell cell, V3 q, V3 fallback)
                    {
                        double s = cell.N.Dot(q) - cell.D;
                        if (Math.Abs(s) < 1e-9) s = cell.N.Dot(fallback) - cell.D;
                        return s >= 0 ? 1 : -1;
                    }

                    int FragmentOf(Cell cell, int side) { return cell.Kind == Kind.Exact ? cell.Frag[side + 1] : cell.Frag[1]; }

                    static void AddCoupling(List<(int, int, double)> to, int fa, int fb, double sigma)
                    {
                        if (fa < 0 || fb < 0 || !(sigma > 1e-12)) return;
                        to.Add((fa, fb, Math.Min(1.0, sigma)));
                    }

                    /// <summary>The coupling(s) across the face between cell a and its neighbour b along link k.</summary>
                    void Couple(Cell a, Cell b, int k, List<(int, int, double)> to)
                    {
                        V3 m = 0.5 * (a.C + b.C);
                        if (a.Kind == Kind.Sampled || b.Kind == Kind.Sampled)
                        {
                            // Sample the shared face: each point couples the pieces of a and b that contain it.
                            Dictionary<(int, int), int> hitsOn = new Dictionary<(int, int), int>();
                            V3[] pts = FaceSamples[k];
                            for (int i = 0; i < pts.Length; i++)
                            {
                                V3 q = a.C + pts[i];
                                int pa = PieceAt(a, q, a.Kind == Kind.Sampled ? FaceSampleNearest[k][i] : null);
                                int pb = PieceAt(b, q, null);
                                if (pa < 0 || pb < 0) continue;
                                hitsOn.TryGetValue((pa, pb), out int c0); hitsOn[(pa, pb)] = c0 + 1;
                            }
                            foreach (KeyValuePair<(int, int), int> kv in hitsOn.OrderBy(kv => kv.Key)) AddCoupling(to, kv.Key.Item1, kv.Key.Item2, (double)kv.Value / pts.Length);
                            return;
                        }
                        if (a.Kind == Kind.Exact && b.Kind == Kind.Exact && SamePlane(a, b))
                        {
                            int o = a.N.Dot(b.N) > 0 ? 1 : -1;
                            double sp = a.SigmaPlus[k];
                            AddCoupling(to, a.Frag[2], b.Frag[o + 1], sp);
                            AddCoupling(to, a.Frag[0], b.Frag[-o + 1], 1 - sp);
                            return;
                        }
                        // A wall crossing this face would cross both cells, so here the face lies wholly on one side of each.
                        int fa =a.Kind == Kind.Exact ? a.Frag[SideOf(a, m, b.C) + 1] : a.Frag[1];
                        int fb = b.Kind == Kind.Exact ? b.Frag[SideOf(b, m, a.C) + 1] : b.Frag[1];
                        AddCoupling(to, fa, fb, 1);
                    }

                    static bool SamePlane(Cell a, Cell b)
                    {
                        double dot = a.N.Dot(b.N);
                        if (Math.Abs(Math.Abs(dot) - 1) > 1e-9) return false;
                        return Math.Abs(a.D - (dot > 0 ? b.D : -b.D)) < 1e-7;
                    }

                    bool Crosses(V3 p, V3 q, List<Hit> hits)
                    {
                        foreach (Hit h in hits) if (SegmentTriangle(p, q, tris[h.Tri])) return true;
                        return false;
                    }

                    /// <summary>
                    /// Whether segment p-q passes through the triangle, taken to reach `seam` beyond its edges so that a link
                    /// through a shared edge or corner of two triangles, or through a seam between two polygons, counts as crossing.
                    /// </summary>
                    bool SegmentTriangle(V3 p, V3 q, Tri t)
                    {
                        V3 dir = q - p;
                        double denom = t.N.Dot(dir);
                        if (Math.Abs(denom) < 1e-14) return false;
                        const double tol = 1e-9;
                        double tt = (t.D - t.N.Dot(p)) / denom;
                        if (tt < -tol || tt > 1 + tol) return false;
                        V3 x = p + tt * dir;
                        double reach = seam + tol;
                        return t.M0.Dot(x - t.A) >= -reach && t.M1.Dot(x - t.B) >= -reach && t.M2.Dot(x - t.C) >= -reach;
                    }

                    // Union-find over fragments.
                    int[] parent;
                    int Find(int i) { while (parent[i] != i) { parent[i] = parent[parent[i]]; i = parent[i]; } return i; }
                    void Union(int i, int j) { i = Find(i); j = Find(j); if (i != j) parent[i] = j; }

                    public string GroupLaps = "";
                    void Group(Cut_Cell_Summary summary, HashSet<long> affected)
                    {
                        System.Diagnostics.Stopwatch gclock = System.Diagnostics.Stopwatch.StartNew(); double glast = 0;
                        Action<string> glap = name => { double now = gclock.Elapsed.TotalSeconds; GroupLaps += $"{(GroupLaps.Length > 0 ? ", " : "")}{name} {now - glast:F1}"; glast = now; };
                        int n = frags.Count;
                        parent = Enumerable.Range(0, n).ToArray();
                        List<(int other, double sigma)>[] adj = new List<(int, double)>[n];
                        for (int i = 0; i < n; i++) adj[i] = new List<(int, double)>();
                        foreach ((int a, int b, double s) in couplings) { adj[a].Add((b, s)); adj[b].Add((a, s)); }

                        // A piece on the far side of the wall from its node: with separate regions it is a cell of its own (if it
                        // couples to anything); otherwise it joins the piece it shares the most open face with on its own side.
                        bool[] region = new bool[n];
                        if (Cut_Cell_Separate_Regions)
                            for (int i = 0; i < n; i++) region[i] = !frags[i].Primary && adj[i].Count > 0;
                        bool Owns(int i) { return frags[i].Primary || region[i]; }
                        for (int i = 0; i < n; i++)
                        {
                            if (Owns(i) || (Cut_Cell_Separate_Regions && !frags[i].Primary)) continue;
                            int best = -1; double bs = 0;
                            foreach ((int o, double s) in adj[i]) if (s > bs || (s == bs && best >= 0 && frags[o].Primary && !frags[best].Primary)) { bs = s; best = o; }
                            if (best >= 0) { Union(i, best); summary.MergedPieces++; summary.FarSideVolume += frags[i].Nu; }
                        }

                        // Groups without any node are isolated slivers: left out.
                        Dictionary<int, double> volume = new Dictionary<int, double>();
                        Dictionary<int, bool> hasNode = new Dictionary<int, bool>();
                        void Tally()
                        {
                            volume.Clear(); hasNode.Clear();
                            for (int i = 0; i < n; i++)
                            {
                                int g = Find(i);
                                volume.TryGetValue(g, out double v); volume[g] = v + frags[i].Nu;
                                hasNode.TryGetValue(g, out bool p); hasNode[g] = p || Owns(i);
                            }
                        }
                        Tally();

                        glap("far-side merge");
                        // Eq. (38): merge each group stiffer than the limit into the neighbour it shares the most face with. Only groups of
                        // affected cells absorb others.
                        bool[] unaffectedRoot = new bool[n];
                        for (int i = 0; i < n; i++) if (!affected.Contains(frags[i].Cell.Key)) unaffectedRoot[Find(i)] = true;
                        double[] extArr = new double[n];
                        bool[] stiffFlag = new bool[n];
                        int passes = 0;
                        for (int pass = 0; pass < 100; pass++)
                        {
                            passes++;
                            Array.Clear(extArr, 0, n);
                            foreach ((int a, int b, double s) in couplings)
                            {
                                int ga = Find(a), gb = Find(b);
                                if (ga == gb || !hasNode[ga] || !hasNode[gb]) continue;
                                extArr[ga] += s; extArr[gb] += s;
                            }
                            List<int> stiff = new List<int>();
                            foreach (int g in volume.Keys) if (extArr[g] > Stiffness_Limit * volume[g] * (1 + 1e-9)) { stiff.Add(g); stiffFlag[g] = true; }
                            if (stiff.Count == 0) break;
                            stiff = stiff.OrderByDescending(g => extArr[g] / volume[g]).ThenBy(g => g).ToList();
                            // Face shared with each neighbour, for the stiff groups only.
                            Dictionary<int, Dictionary<int, double>> shared = new Dictionary<int, Dictionary<int, double>>();
                            foreach ((int a, int b, double s) in couplings)
                            {
                                int ga = Find(a), gb = Find(b);
                                if (ga == gb || !hasNode[ga] || !hasNode[gb]) continue;
                                if (stiffFlag[ga]) { if (!shared.TryGetValue(ga, out var sa)) shared[ga] = sa = new Dictionary<int, double>(); sa.TryGetValue(gb, out double x1); sa[gb] = x1 + s; }
                                if (stiffFlag[gb]) { if (!shared.TryGetValue(gb, out var sb)) shared[gb] = sb = new Dictionary<int, double>(); sb.TryGetValue(ga, out double x2); sb[ga] = x2 + s; }
                            }
                            foreach (int g in stiff) stiffFlag[g] = false;
                            HashSet<int> touched = new HashSet<int>();
                            bool merged = false;
                            foreach (int g in stiff)
                            {
                                if (touched.Contains(g) || !shared.ContainsKey(g)) continue;
                                var options = shared[g].Where(kv => !unaffectedRoot[kv.Key]).OrderByDescending(kv => kv.Value).ThenBy(kv => kv.Key).ToList();
                                if (options.Count == 0) continue;
                                int target = options[0].Key;
                                if (touched.Contains(target)) continue;
                                summary.StiffnessMerges++; summary.StiffnessMergedVolume += volume[g];
                                bool unaff = unaffectedRoot[g] || unaffectedRoot[target];
                                Union(g, target);
                                unaffectedRoot[Find(target)] = unaff;
                                touched.Add(g); touched.Add(Find(target));
                                summary.MergedPieces++;
                                merged = true;
                            }
                            Tally();
                            if (!merged) break;
                        }
                        glap($"stiffness merges ({passes} passes)");

                        // Members, representatives and couplings per group.
                        Dictionary<int, List<int>> members = new Dictionary<int, List<int>>();
                        for (int i = 0; i < n; i++) { int g = Find(i); if (!members.TryGetValue(g, out var l)) members[g] = l = new List<int>(); l.Add(i); }
                        Dictionary<int, Dictionary<int, double>> link = new Dictionary<int, Dictionary<int, double>>();
                        foreach ((int a, int b, double s) in couplings)
                        {
                            int ga = Find(a), gb = Find(b);
                            if (ga == gb || !hasNode[ga] || !hasNode[gb]) continue;
                            if (!link.TryGetValue(ga, out var la)) link[ga] = la = new Dictionary<int, double>();
                            la.TryGetValue(gb, out double x1); la[gb] = x1 + s;
                            if (!link.TryGetValue(gb, out var lb)) link[gb] = lb = new Dictionary<int, double>();
                            lb.TryGetValue(ga, out double x2); lb[ga] = x2 + s;
                        }

                        Dictionary<int, Cell> rep = new Dictionary<int, Cell>();
                        HashSet<int> regionGroups = new HashSet<int>();
                        foreach (KeyValuePair<int, List<int>> g in members)
                        {
                            if (!hasNode[g.Key]) { summary.DroppedVolume += volume[g.Key]; continue; }
                            Fragment best = g.Value.Select(i => frags[i]).Where(f => f.Primary).OrderByDescending(f => f.Nu).ThenBy(f => f.Cell.Kind == Kind.Untouched ? 0 : 1).FirstOrDefault();
                            if (best == null)
                            {
                                // A far-side region with no node of its own: its unknown is hosted by a node in the grid.
                                best = g.Value.Select(i => frags[i]).OrderByDescending(f => f.Nu).First();
                                regionGroups.Add(g.Key);
                                summary.SeparateRegions++; summary.SeparateRegionVolume += volume[g.Key];
                            }
                            rep[g.Key] = best.Cell;
                            double ext = link.TryGetValue(g.Key, out var lg) ? lg.Values.Sum() : 0;
                            if (g.Value.Any(i => frags[i].Cell.Kind != Kind.Untouched)) summary.MaxStiffness = Math.Max(summary.MaxStiffness, ext / volume[g.Key]);
                        }

                        glap("members+links");
                        // New nodes. A group that is a single, whole, unaffected cell keeps its original node.
                        Dictionary<long, Node> replaced = new Dictionary<long, Node>();
                        Dictionary<int, Cut_Cell_Node> repNode = new Dictionary<int, Cut_Cell_Node>();
                        bool Original(int g) { return members[g].Count == 1 && !affected.Contains(frags[members[g][0]].Cell.Key); }
                        foreach (int g in rep.Keys)
                        {
                            if (Original(g)) continue;
                            Cell r = rep[g];
                            Cut_Cell_Node node = new Cut_Cell_Node(RDD_Location(min, r.X, r.Y, r.Z, dx, dy, dz));
                            repNode[g] = node;
                            if (!regionGroups.Contains(g)) replaced[r.Key] = node;
                        }
                        foreach (int g in rep.Keys)
                        {
                            if (Original(g)) continue;
                            foreach (int i in members[g])
                            {
                                Fragment f = frags[i];
                                if (!f.Primary || f.Cell == rep[g]) continue;
                                replaced[f.Cell.Key] = new Slave_Node(RDD_Location(min, f.Cell.X, f.Cell.Y, f.Cell.Z, dx, dy, dz), repNode[g]);
                            }
                        }
                        foreach (long key in affected)
                        {
                            Cell cell = cells[key];
                            if (replaced.ContainsKey(key)) continue;
                            int any = cell.Frag.Where(i => i >= 0 && rep.ContainsKey(Find(i)) && !Original(Find(i))).DefaultIfEmpty(-1).First();
                            if (any >= 0) replaced[key] = new Slave_Node(RDD_Location(min, cell.X, cell.Y, cell.Z, dx, dy, dz), repNode[Find(any)]);
                        }
                        // Regions ride along with the node in their cell's slot, or take the slot if nothing else has.
                        foreach (int g in regionGroups)
                        {
                            long key = rep[g].Key;
                            if (!replaced.TryGetValue(key, out Node host)) replaced[key] = repNode[g];
                            else if (host is Cut_Cell_Node cn) cn.AddCompanion(repNode[g]);
                            else if (host is Slave_Node sn) sn.AddCompanion(repNode[g]);
                            else throw new InvalidOperationException("cut cells: no host for a far-side region");
                        }
                        foreach (int g in rep.Keys)
                        {
                            if (Original(g)) continue;
                            List<Node> nb = new List<Node>(); List<double> w = new List<double>();
                            if (link.TryGetValue(g, out var lg))
                                foreach (KeyValuePair<int, double> kv in lg)
                                {
                                    if (!rep.ContainsKey(kv.Key)) continue;
                                    Cell other = rep[kv.Key];
                                    Node target = Original(kv.Key) ? F[other.X][other.Y][other.Z] : repNode[kv.Key];
                                    nb.Add(target); w.Add(w0 * kv.Value / volume[g]);
                                }
                            repNode[g].Set(nb.ToArray(), w.ToArray(), Filters(members[g], volume[g]), volume[g]);
                        }
                        glap("nodes+filters");
                        foreach (KeyValuePair<long, Node> kv in replaced)
                        {
                            long key = kv.Key;
                            int x = (int)(key / yStride), y = (int)(key % yStride / zStride), z = (int)(key % zStride);
                            F[x][y][z] = kv.Value;
                        }
                        // Unaffected neighbours keep their links; point them at the new nodes.
                        foreach (KeyValuePair<long, Node> kv in replaced)
                        {
                            long key = kv.Key;
                            int x = (int)(key / yStride), y = (int)(key % yStride / zStride), z = (int)(key % zStride);
                            for (int k = 0; k < 12; k++)
                            {
                                if (!Neighbour(x, y, z, k, out int nx, out int ny, out int nz)) continue;
                                if (affected.Contains(KeyOf(nx, ny, nz))) continue;
                                if (F[nx][ny][nz] is RDD_Node rn && rn.Links2 != null) rn.Links2[Opp[k]] = kv.Value;
                            }
                        }
                    }

                    /// <summary>
                    /// Boundary filters of a group, one per material. The scheme's wall term is lambda M sum (v[n+1] - v[n-1])
                    /// with lambda M = c T S / (2 V), so M = dx S / (2 V) for wall area S in a group of volume V.
                    /// </summary>
                    double wallTotal;
                    readonly Dictionary<Environment.Material, double> wallByMaterial = new Dictionary<Environment.Material, double>();
                    double unassigned;
                    readonly ConcurrentBag<(Cell cell, V3 q, Environment.Material mat, double area)> orphans = new ConcurrentBag<(Cell, V3, Environment.Material, double)>();

                    /// <summary>Gives the wall of a sliver too thin to hold a sample to the piece across the face it lies against.</summary>
                    void AssignOrphanWalls()
                    {
                        foreach ((Cell cell, V3 q, Environment.Material mat, double area) in orphans.OrderBy(o => o.cell.Key).ThenBy(o => o.q.X).ThenBy(o => o.q.Y).ThenBy(o => o.q.Z))
                        {
                            int kBest = 0; double best = double.NegativeInfinity;
                            for (int k = 0; k < 12; k++) { double e = E[k].Dot(q - cell.C); if (e > best) { best = e; kBest = k; } }
                            int f = -1;
                            if (Neighbour(cell.X, cell.Y, cell.Z, kBest, out int nx, out int ny, out int nz) && cells.TryGetValue(KeyOf(nx, ny, nz), out Cell nb))
                            {
                                V3 across = q + (Math.Max(0, (2 - best) / 4) + 0.01) * E[kBest];
                                f = PieceAt(nb, across, null);
                            }
                            if (f < 0) { unassigned += area; continue; }
                            frags[f].Wall.TryGetValue(mat, out double a0);
                            frags[f].Wall[mat] = a0 + area;
                        }
                    }

                    IIR_DIF[] Filters(List<int> group, double volumeInCells)
                    {
                        Dictionary<Environment.Material, double> weight = new Dictionary<Environment.Material, double>();
                        double groupVolume = volumeInCells * v0Physical;
                        foreach (int i in group)
                            foreach (KeyValuePair<Environment.Material, double> w in frags[i].Wall)
                            {
                                wallTotal += w.Value;
                                wallByMaterial.TryGetValue(w.Key, out double t0); wallByMaterial[w.Key] = t0 + w.Value;
                                weight.TryGetValue(w.Key, out double m0);
                                weight[w.Key] = m0 + dx * w.Value / (2 * groupVolume);
                            }
                        return weight.Select(w => new DIF_IWB_2p(w.Key, fs, w.Value, lambda)).ToArray();
                    }
                }
            }

            /// <summary>
            /// A (possibly merged) cut cell: p[n+1] = sum_k W_k p_k + (2 - sum_k W_k) p[n] - p[n-1] - wall term, with
            /// W_k = c^2 T^2 S_k / (V h) for each open face area S_k, cell volume V and link length h.
            /// </summary>
            public sealed class Cut_Cell_Node : Acoustic_Compact_FDTD.Node
            {
                Acoustic_Compact_FDTD.Node[] neighbours = new Acoustic_Compact_FDTD.Node[0];
                double[] W = new double[0];
                double centre;
                Acoustic_Compact_FDTD.IIR_DIF[] filters = new Acoustic_Compact_FDTD.IIR_DIF[0];
                double denom = 1, kappa;

                public Cut_Cell_Node(Point loc) : base(loc) { }

                /// <summary>Fluid volume of the (merged) cell, in full cells.</summary>
                public double Volume { get; private set; }
                /// <summary>Sum of the boundary filters' a_b (lambda M b0).</summary>
                public double Kappa { get { return kappa; } }

                internal void Set(Acoustic_Compact_FDTD.Node[] nb, double[] w, Acoustic_Compact_FDTD.IIR_DIF[] f, double volume)
                {
                    neighbours = nb; W = w; filters = f; Volume = volume;
                    centre = 2 - w.Sum();
                    kappa = f.Sum(x => x.a_b);
                    denom = 1 / (1 + kappa);
                }

                public int NeighbourCount { get { return neighbours.Length; } }
                /// <summary>The cells this one is coupled to and the weight of each coupling (for diagnostics).</summary>
                public IReadOnlyList<Acoustic_Compact_FDTD.Node> Neighbours { get { return neighbours; } }
                public IReadOnlyList<double> Weights { get { return W; } }
                public double CouplingSum { get { return W.Sum(); } }

                public override void Link_Nodes(ref Acoustic_Compact_FDTD.Node[][][] Frame, int x, int y, int z) { }

                public override void UpdateP()
                {
                    double s = centre * Pn - Pn_1;
                    for (int i = 0; i < neighbours.Length; i++) s += W[i] * neighbours[i].P;
                    double mem = 0;
                    for (int i = 0; i < filters.Length; i++) mem += filters[i].g_b_term();
                    Pnf = (s + kappa * Pn_1 + mem) * denom;
                    for (int i = 0; i < filters.Length; i++) filters[i].Update(Pnf, Pn_1);
                    for (int i = 0; i < companions.Length; i++) companions[i].UpdateP();
                }

                public override void UpdateT()
                {
                    base.UpdateT();
                    for (int i = 0; i < companions.Length; i++) companions[i].UpdateT();
                }

                public override void reset()
                {
                    base.reset();
                    for (int i = 0; i < companions.Length; i++) companions[i].reset();
                }

                Acoustic_Compact_FDTD.Node[] companions = new Acoustic_Compact_FDTD.Node[0];
                /// <summary>Cells not in the grid (far-side regions) that this node updates along with itself.</summary>
                internal void AddCompanion(Acoustic_Compact_FDTD.Node n) { companions = companions.Append(n).ToArray(); }
                public int CompanionCount { get { return companions.Length; } }
                /// <summary>The far-side cells this node updates along with itself.</summary>
                public IReadOnlyList<Acoustic_Compact_FDTD.Node> Companions { get { return companions; } }
            }

            /// <summary>A cell merged into another: its pressure is the group's.</summary>
            public sealed class Slave_Node : Acoustic_Compact_FDTD.Node
            {
                readonly Acoustic_Compact_FDTD.Node master;
                public Slave_Node(Point loc, Acoustic_Compact_FDTD.Node master) : base(loc) { this.master = master; }
                /// <summary>The node of the group this cell was merged into.</summary>
                public Acoustic_Compact_FDTD.Node Master { get { return master; } }
                public override void Link_Nodes(ref Acoustic_Compact_FDTD.Node[][][] Frame, int x, int y, int z) { }
                public override void UpdateP() { for (int i = 0; i < companions.Length; i++) companions[i].UpdateP(); }
                public override void UpdateT()
                {
                    Pn_1 = Pn;
                    Pn = master.Pnf * Attenuation;
                    for (int i = 0; i < companions.Length; i++) companions[i].UpdateT();
                }
                public override void reset()
                {
                    base.reset();
                    for (int i = 0; i < companions.Length; i++) companions[i].reset();
                }
                Acoustic_Compact_FDTD.Node[] companions = new Acoustic_Compact_FDTD.Node[0];
                internal void AddCompanion(Acoustic_Compact_FDTD.Node n) { companions = companions.Append(n).ToArray(); }
                /// <summary>The far-side cells this node updates along with itself.</summary>
                public IReadOnlyList<Acoustic_Compact_FDTD.Node> Companions { get { return companions; } }
            }
        }
    }
}
