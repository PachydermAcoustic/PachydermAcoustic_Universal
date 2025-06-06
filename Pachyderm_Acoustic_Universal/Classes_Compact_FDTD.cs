﻿using System;
using System.Linq;
using System.Collections.Generic;
using Pachyderm_Acoustic.Environment;
using Hare.Geometry;
using System.Threading;
using Eto.Forms;

namespace Pachyderm_Acoustic
{
    namespace Numeric
    {
        namespace TimeDomain
        {
            public partial class Acoustic_Compact_FDTD : Simulation_Type
            {
                Scene Rm;
                public double tmax;       //seconds
                private double rho0 = 1.21;    //density of air kg/m^3
                public double dt;
                public double dx, dy, dz;
                public readonly int xDim, yDim, zDim;
                public Node[][][] PFrame;            //pressure scalar field initialisation
                Vector[] Dir = new Vector[13];
                Point[] Orig = new Point[13];
                double fmax;
                protected AABB Bounds;
                protected AABB Bounds_Inner;
                PML Layers;
                public Signal_Driver_Compact SD;
                public Microphone_Compact Mic;
                public int no_of_Layers;
                public int n;
                double time_ms;
                int threadct = System.Environment.ProcessorCount;
                double floorplaneoverride = 0;

                //PML Settings
                int PML_LayerNo = 50;
                double PML_MaxAtten = 0.13;
                
                public enum GridType
                {
                    Freefield,
                    ScatteringLab,
                    TransparencyLab,
                    Terrain
                }

                public static Hare.Geometry.Point RDD_Location(Hare.Geometry.Point MinPt, int x, int y, int z, double dx, double dy, double dz)
                {
                    int mod = x % 2;
                    return new Point(MinPt.x + (((double)x - 0.5) * dx), MinPt.y + 2 * (((double)y + (0.5 - 0.5 * mod)) * dy), MinPt.z + 2 * (((double)z + (0.5 - 0.5 * mod)) * dz));
                }

                public int[] RDD_Location(Hare.Geometry.Point p)
                {
                    int[] loc = new int[3];
                    loc[0] = (int)Math.Floor((p.x - Bounds.Min_PT.x) / dx);
                    loc[1] = (int)Math.Floor((p.y - Bounds.Min_PT.y) / dy);// (4 * dx / Utilities.Numerics.rt2));
                    loc[2] = (int)Math.Floor((p.z - Bounds.Min_PT.z) / dz);// (2 * dx / Utilities.Numerics.rt2));

                    //p -= Bounds.Min_PT;
                    //p.x -= 05 * dx;
                    //loc[0] = (int)Math.Round(p.x / dx);
                    //int mod = loc[0] % 2;
                    //p.y -= (0.5 - 0.5 * mod) * 2 * dy;
                    //loc[1] =  (int)Math.Round(p.y / (2 * dy));
                    //p.z -= (0.5 - 0.5 * mod) * 2 * dz;
                    //loc[2] = (int)Math.Round(p.z / (2 * dz));
                    return loc;
                }

                public Hare.Geometry.Point RDD_Location(int x, int y, int z)
                {
                    int mod = x % 2;
                    return new Point(Bounds.Min_PT.x + (((double)x - 0.5) * dx), Bounds.Min_PT.y + 2 * (((double)y + (0.5 - 0.5 * mod)) * dy), Bounds.Min_PT.z + 2 * (((double)z + (0.5 - 0.5 * mod)) * dz));
                }

                public Acoustic_Compact_FDTD(Scene Rm_in, ref Signal_Driver_Compact S_in, ref Microphone_Compact M_in, double fmax_in, double tmax_ms_in, GridType GT, Point SampleOrigin, double mindimx, double mindimy, double mindimz, bool PML = true, double floorplane = 0)
                {
                    floorplaneoverride = floorplane;
                    Rm = Rm_in;
                    SD = S_in;
                    Mic = M_in;
                    fmax = fmax_in;
                    tmax = tmax_ms_in;
                    Rm.partition();
                    if (GT == GridType.Freefield)
                    {
                        Build_Freefield_FVM13(ref xDim, ref yDim, ref zDim, PML, mindimx, mindimy, mindimz);
                        SD.Connect_Grid_Freefield(PFrame, Bounds, dx, dy, dz, tmax, dt, no_of_Layers);
                        Mic.Connect_Grid_UniqueOnly_Freefield(PFrame, Bounds, dx, tmax, dt, no_of_Layers);
                    }
                    else if (GT == GridType.ScatteringLab)
                    {
                        Build_ScatteringLaboratory_FVM13(ref xDim, ref yDim, ref zDim, PML, SampleOrigin, mindimx, mindimy, mindimz);
                        SD.Connect_Grid_Laboratory(PFrame, Bounds, Bounds_Inner, dx, dy, dz, tmax, dt, no_of_Layers);
                        Mic.Connect_Grid_Hemisphere_Laboratory(PFrame, Bounds, SampleOrigin, 1, dx, tmax, dt, no_of_Layers);
                    }
                    else if (GT == GridType.TransparencyLab)
                    {
                        //TODO: Build a custom lab with freefield condition at bottom boundary...
                        Build_TransparencyLaboratory_FVM13(ref xDim, ref yDim, ref zDim, PML, SampleOrigin, mindimx, mindimy, mindimz);
                        SD.Connect_Grid_Laboratory(PFrame, Bounds, Bounds_Inner, dx, dy, dz, tmax, dt, no_of_Layers);
                        Mic.Connect_Grid_UniqueOnly_Laboratory(PFrame, Bounds, dx, tmax, dt, no_of_Layers);
                    }
                    else if (GT == GridType.Terrain)
                    {
                        Hare.Geometry.Point Origin = (Rm_in.Max() + Rm.Min()) / 2;
                        Build_ScatteringLaboratory_FVM13(ref xDim, ref yDim, ref zDim, PML, Origin, mindimx, mindimy, mindimz);
                        SD.Connect_Grid_Laboratory(PFrame, Bounds, Bounds, dx, dy, dz, tmax, dt, no_of_Layers);
                        Mic.Connect_Grid_UniqueOnly_Laboratory(PFrame, Bounds, dx, tmax, dt, no_of_Layers);
                    }
                }

                public void Build_TransparencyLaboratory_FVM13(ref int xDim, ref int yDim, ref int zDim, bool PML, Point SampleOrigin, double xmin, double ymin, double zmin)
                {
                    double rt2 = Math.Sqrt(2);
                    double rt3 = Math.Sqrt(3);

                    double dydz = Rm.Sound_speed(0) / fmax * .1;
                    dx = 2 * dydz / Math.Sqrt(2);

                    Bounds = new AABB(Rm.Min() - new Vector(.05 * dx, .05 * dydz, .05 * dydz), Rm.Max() + new Point(.05 * dx, .05 * dydz, .05 * dydz));
                    Bounds_Inner = new AABB(Bounds.Min_PT.x, Bounds.Min_PT.y, Bounds.Min_PT.z, Bounds.Max_PT.x, Bounds.Max_PT.y, Bounds.Max_PT.z);

                    no_of_Layers = 0;
                    double max_Layer = 0;

                    if (PML)
                    {
                        no_of_Layers = PML_LayerNo;
                        max_Layer = PML_MaxAtten;
                    }

                    //double x_length = Bounds.X_Length(), y_length = Bounds.Y_Length(), z_length = Bounds.Z_Length();
                    double x_length, y_length, z_length;
                    Point MinPt = new Hare.Geometry.Point(SampleOrigin.x, SampleOrigin.y, SampleOrigin.z);
                    //if (x_length < xmin)
                    //{
                    x_length = xmin;
                    MinPt.x -= x_length / 2;
                    //}
                    //if (y_length < ymin)
                    //{
                    y_length = ymin;
                    MinPt.y -= y_length / 2;
                    //}
                    //if (z_length < zmin)
                    //{
                    z_length = zmin;
                    MinPt.z -= z_length /2;
                    //}

                    x_length += (no_of_Layers * 4 + 1) * dx;
                    y_length += (no_of_Layers * 4 + 1) * dydz;
                    z_length += (no_of_Layers * 4 + 1) * dydz;

                    //estimated distance between nodes
                    xDim = (int)Math.Ceiling(x_length / (dx));                                //set number of nodes in x direction
                    dx = x_length / xDim;                                                   //refined distance between nodes
                    yDim = (int)Math.Ceiling(y_length / dydz);                                //set number of nodes in y direction
                    dy = y_length / yDim;
                    zDim = (int)Math.Ceiling(z_length / dydz);                                //set number of nodes in z direction
                    dz = z_length / zDim;

                    dt = dy * rt2 / (Rm.Sound_speed(0));                           //set time step small enough to satisfy courrant condition
                    dxrt2 = dx * rt2;
                    dxrt3 = dx * rt3;

                    Dir = new Vector[13]{
                        new Vector(-1, 0, 0),
                        new Vector(0, -1, 0),
                        new Vector(0, 0, -1),

                        new Vector(-1/rt2,-1/rt2,double.Epsilon),
                        new Vector(1/rt2, -1/rt2,double.Epsilon),
                        new Vector(-1/rt2, double.Epsilon, 1/rt2),
                        new Vector(double.Epsilon, -1/rt2, 1/rt2),
                        new Vector(1/rt2, double.Epsilon, 1/rt2),
                        new Vector(double.Epsilon, 1/rt2, 1/rt2),

                        new Vector(-dx, -dy, dz),
                        new Vector(dx, -dy, dz),
                        new Vector(dx, dy, dz),
                        new Vector(-dx, dy, dz)
                    };

                    foreach (Vector V in Dir) V.Normalize();

                    PFrame = new Node[xDim][][];// yDim, zDim];                               //pressure scalar field initialisation

                    List<Bound_Node_RDD> Bound = new List<Bound_Node_RDD>();

                    int xDimt = xDim;
                    int yDimt = yDim;
                    int zDimt = zDim;

                    MinPt.x -= 2 * dx * no_of_Layers;
                    MinPt.y -= 2 * dy * no_of_Layers; 
                    MinPt.z -= 2 * dz * no_of_Layers;
                    Bounds = new AABB(MinPt, MinPt + new Point(x_length, y_length, z_length));

                    //System.Threading.Tasks.Parallel.For(0, xDim, (x) =>
                    for (int x = 0; x < PFrame.Length; x++)
                    {
                        int mod = x % 2;
                        PFrame[x] = new Node[(int)(Math.Floor((double)yDimt / 2) + yDimt % 2 * mod)][];
                        Random Rnd = new Random(x);
                        for (int y = 0; y < PFrame[x].Length; y++)
                        {
                            PFrame[x][y] = new Node[(int)(Math.Floor((double)zDimt / 2) + yDimt % 2 * mod)];
                            for (int z = 0; z < PFrame[x][y].Length; z++)
                            {
                                List<Environment.Material> abs;
                                List<Bound_Node.Boundary> BDir;
                                Point Loc = Acoustic_Compact_FDTD.RDD_Location(MinPt, x, y, z, dx, dy, dz); //new Point(MinPt.x + 2 * (((double)x - 0.5) * dx), MinPt.y + 2 * (((double)y + (0.5 - 0.5 * mod)) * dy), MinPt.z + 2 * (((double)z + (0.5 - 0.5 * mod)) * dz));
                                if (!Intersect_13Pt(Loc, SD.frequency, out BDir, out abs, ref Rnd))
                                {
                                    PFrame[x][y][z] = new RDD_Node(Loc);//, rho0, dt, dx, Rm.Sound_speed, new int[] { x, y, z });
                                }
                                else
                                {
                                    PFrame[x][y][z] = new Bound_Node_RDD_MaterialFilter(Loc, rho0, dt, dx, Rm.Sound_speed(0), new int[] { x, y, z }, abs, BDir); // abs,
                                    Bound.Add(PFrame[x][y][z] as Bound_Node_RDD_MaterialFilter);
                                }
                            }
                        }
                    }//);

                    //Node.Attenuation = Math.Sqrt(Math.Pow(10, -.1 * Rm.AttenuationPureTone(Bounds.Min_PT, SD.frequency) * dt));//PFrame[0][0][0].Pt, SD.frequency

                    bool failed = false;
                    //Make Mesh Templates:
                    //Build_Mesh_Sections();
                    //for (int x = 0; x < PFrame.Length; x++)
                    System.Threading.Tasks.Parallel.For(0, xDim, (x) =>
                    {
                        for (int y = 0; y < PFrame[x].Length; y++)
                        {
                            for (int z = 0; z < PFrame[x][y].Length; z++)
                            {
                                PFrame[x][y][z].Link_Nodes(ref PFrame, x, y, z);
                            }
                        }
                    });

                    yDim = PFrame[0].Length;
                    zDim = PFrame[0][0].Length;

                    for (int i = 0; i < Bound.Count; i++) if (Bound[i] != null) Bound[i].Complete_Boundary();
                    if (failed) return;

                    //Set up PML...
                    Layers = new Acoustic_Compact_FDTD.PML(no_of_Layers, max_Layer, PFrame, false);
                }

                public void Build_ScatteringLaboratory_FVM13(ref int xDim, ref int yDim, ref int zDim, bool PML, Point SampleOrigin, double xmin, double ymin, double zmin)
                {
                    double rt2 = Math.Sqrt(2);
                    double rt3 = Math.Sqrt(3);

                    double dydz = Rm.Sound_speed(0) / fmax * .1;
                    dx = 2 * dydz / Math.Sqrt(2);

                    Bounds = new AABB(Rm.Min() - new Vector(.025 * dx, .05 * dydz, .05 * dydz), Rm.Max() + new Point(.025 * dx, .05 * dydz, .05 * dydz));
                    Bounds.Min_PT.z = floorplaneoverride;
                    Bounds_Inner = new AABB(Bounds.Min_PT.x, Bounds.Min_PT.y, Bounds.Min_PT.z, Bounds.Max_PT.x, Bounds.Max_PT.y, Bounds.Max_PT.z);

                    no_of_Layers = 0;
                    double max_Layer = 0;

                    if (PML)
                    {
                        no_of_Layers = PML_LayerNo;
                        max_Layer = PML_MaxAtten;
                    }

                    double x_length = Bounds.X_Length(), y_length = Bounds.Y_Length(), z_length = Bounds.Z_Length();
                    Point MinPt = new Hare.Geometry.Point(SampleOrigin.x, SampleOrigin.y, 0 + floorplaneoverride);

                    //MinPt.x -= (x_length / 4);// + no_of_Layers / 2 * dx);
                    //MinPt.y -= (y_length / 2);// + no_of_Layers / 2 * dydz);

                    List<Bound_Node_RDD> Bound = new List<Bound_Node_RDD>();

                    if (x_length < xmin)
                    {
                        x_length = xmin;
                        MinPt.x -= (xmin - x_length) / 2;
                    }
                    if (y_length < ymin)
                    {
                        y_length = ymin;
                        MinPt.y -= (ymin - y_length) / 2;
                    }
                    if (z_length < zmin - floorplaneoverride)
                    {
                        z_length = zmin - floorplaneoverride;
                        MinPt.z -= (zmin - z_length) / 2;
                    }

                    x_length += (no_of_Layers * 4 + 1) * dx;
                    y_length += (no_of_Layers * 4 + 1) * dydz;
                    z_length += (no_of_Layers * 2 + 1) * dydz;

                    //estimated distance between nodes
                    xDim = (int)Math.Ceiling(x_length / dx);                                //set number of nodes in x direction
                    dx = x_length / xDim;                                                   //refined distance between nodes
                    yDim = (int)Math.Ceiling(y_length / dydz);                                //set number of nodes in y direction
                    dy = y_length / yDim;
                    zDim = (int)Math.Ceiling(z_length / dydz);                                //set number of nodes in z direction
                    dz = z_length / zDim;

                    dt = dy * rt2 / (Rm.Sound_speed(0));                           //set time step small enough to satisfy courrant condition
                    dxrt2 = dx * rt2;
                    dxrt3 = dx * rt3;

                    ////MinPt -= new Point(dx * no_of_Layers + 1, dy * no_of_Layers + 1, 0);
                    //x_length += (no_of_Layers * 2) * dx;
                    //y_length += (no_of_Layers * 2 * rt2) * dy;
                    //z_length += (no_of_Layers * rt2) * dz;

                    //xDim = (int)Math.Ceiling(x_length / dx);                                //set number of nodes in x direction
                    //yDim = (int)Math.Ceiling(y_length / dy);                                //set number of nodes in y direction
                    //zDim = (int)Math.Ceiling(z_length / dz);                                //set number of nodes in z direction

                    //estimated distance between nodes
                    dt = dy * rt2 / (Rm.Sound_speed(0));                           //set time step small enough to satisfy courrant condition
                    dxrt2 = dx * rt2;
                    dxrt3 = dx * rt3;

                    Dir = new Vector[13]{
                        new Vector(-1, 0, 0),
                        new Vector(0, -1, 0),
                        new Vector(0, 0, -1),

                        new Vector(-1/rt2,-1/rt2,double.Epsilon),
                        new Vector(1/rt2, -1/rt2,double.Epsilon),
                        new Vector(-1/rt2, double.Epsilon, 1/rt2),
                        new Vector(double.Epsilon, -1/rt2, 1/rt2),
                        new Vector(1/rt2, double.Epsilon, 1/rt2),
                        new Vector(double.Epsilon, 1/rt2, 1/rt2),

                        new Vector(-dx, -dy, dz),
                        new Vector(dx, -dy, dz),
                        new Vector(dx, dy, dz),
                        new Vector(-dx, dy, dz)
                    };

                    foreach (Vector V in Dir) V.Normalize();

                    int xDimt = xDim;
                    int yDimt = yDim / 2;
                    int zDimt = zDim / 2;

                    PFrame = new Node[xDimt][][];// yDim, zDim];                               //pressure scalar field initialisation

                    int threadct = System.Environment.ProcessorCount;
                    List<System.Threading.Thread> T = new List<System.Threading.Thread>();

                    floorplaneoverride = Math.Round(floorplaneoverride / dydz) * dydz;

                    MinPt.x -= x_length / 2;
                    MinPt.y -= y_length / 2;
                    MinPt.z = 0 + floorplaneoverride;
                    Bounds = new AABB(MinPt, MinPt + new Point(x_length, y_length, z_length));

                    //System.Threading.Tasks.Parallel.For(0, xDim, (x) =>
                    //for (int x = 0; x < PFrame.Length; x++)
                    for (int p = 0; p < threadct; p++)
                    {
                        T.Add(new System.Threading.Thread(proc_id =>
                        {
                            for (int x = (int)Math.Floor((double)((int)proc_id  * PFrame.Length) / threadct); x < Math.Floor((double)(((int)proc_id + 1) * PFrame.Length) / threadct); x++)
                            {
                                int mod = x % 2;
                                PFrame[x] = new Node[(int)(Math.Floor((double)yDimt) + yDimt % 2 * mod)][];
                                Random Rnd = new Random(x);
                                for (int y = 0; y < PFrame[x].Length; y++)
                                {
                                    PFrame[x][y] = new Node[(int)(Math.Floor((double)zDimt) + yDimt % 2 * mod)];
                                    for (int z = 0; z < PFrame[x][y].Length; z++)
                                    {
                                        Point Loc = Acoustic_Compact_FDTD.RDD_Location(MinPt, x, y, z, dx, dy, dz); //new Point(MinPt.x + 2 * (((double)x - 0.5) * dx), MinPt.y + 2 * (((double)y + (0.5 - 0.5 * mod)) * dy), MinPt.z + 2 * (((double)z + (0.5 - 0.5 * mod)) * dz));
                                        if (Rm is Empty_Scene)
                                        {
                                            List<Bound_Node.Boundary> BDir = new List<Bound_Node.Boundary>();
                                            BDir.Add(Bound_Node.Boundary.AZNeg);
                                            if (z == 0)
                                            {
                                                PFrame[x][y][z] = new Bound_Node_RDD(Loc, rho0, dt, dx, Rm.Sound_speed(0), new int[] { x, y, z }, BDir); //abs,
                                                Bound.Add(PFrame[x][y][z] as Bound_Node_RDD);
                                            }
                                            else
                                            {
                                                PFrame[x][y][z] = new RDD_Node(Loc);//, rho0, dt, dx, Rm.Sound_speed, new int[] { x, y, z });
                                            }
                                        }
                                        else
                                        {
                                            List<Environment.Material> abs;
                                            List<Bound_Node.Boundary> BDir;

                                            if (!Intersect_13Pt(Loc, SD.frequency, out BDir, out abs, ref Rnd))
                                            {
                                                PFrame[x][y][z] = new RDD_Node(Loc);//, rho0, dt, dx, Rm.Sound_speed, new int[] { x, y, z });
                                            }
                                            else
                                            {
                                                PFrame[x][y][z] =
                                                    new Bound_Node_RDD(Loc, rho0, dt, dx, Rm.Sound_speed(0), new int[] { x, y, z }, BDir); //abs,
                                                Bound.Add(PFrame[x][y][z] as Bound_Node_RDD);
                                            }
                                        }
                                    }
                                }
                            }
                        }));
                        T[p].Start(p);
                    }//);

                    foreach (System.Threading.Thread t in T) t.Join();

                    //Node.Attenuation = Math.Sqrt(Math.Pow(10, -.1 * Rm.AttenuationPureTone(Bounds.Min_PT, SD.frequency) * dt));//PFrame[0][0][0].Pt, SD.frequency

                    bool failed = false;
                    //Make Mesh Templates:
                    //Build_Mesh_Sections();
                    //for (int x = 0; x < PFrame.Length; x++)
                    //System.Threading.Tasks.Parallel.For(0, xDim, (x) =>
                    //{

                    T = new List<System.Threading.Thread>();

                    for (int p = 0; p < threadct; p++)
                    {
                        T.Add(new System.Threading.Thread(proc_id =>
                        {
                            for (int x = (int)Math.Floor((double)((int)proc_id * PFrame.Length) / threadct); x < Math.Floor((double)(((int)proc_id + 1) * PFrame.Length) / threadct); x++)
                            {
                                for (int y = 0; y < PFrame[x].Length; y++)
                                {
                                    for (int z = 0; z < PFrame[x][y].Length; z++)
                                    {
                                        PFrame[x][y][z].Link_Nodes(ref PFrame, x, y, z);
                                    }
                                }
                            }
                        }));
                        T[p].Start(p);
                    }

                    foreach (System.Threading.Thread t in T) t.Join();

                    yDim = PFrame[0].Length;
                    zDim = PFrame[0][0].Length;

                    for (int i = 0; i < Bound.Count; i++) if (Bound[i] != null) Bound[i].Complete_Boundary();
                    if (failed) return;

                    //Set up PML...
                    Layers = new Acoustic_Compact_FDTD.PML(no_of_Layers, max_Layer, PFrame, false);
                }

                public void Build_Freefield_FVM13(ref int xDim, ref int yDim, ref int zDim, bool PML, double xmin, double ymin, double zmin)
                {
                    double rt2 = Math.Sqrt(2);
                    double rt3 = Math.Sqrt(3);

                    double dydz = Rm.Sound_speed(0) / fmax * .1;
                    dx = 2 * dydz / Math.Sqrt(2);

                    Bounds = new AABB(Rm.Min() - new Vector(.05 * dx, .05 * dydz, .05 * dydz), Rm.Max() + new Point(.05 * dx, .05 * dydz, .05 * dydz));
                    Bounds_Inner = new AABB(Bounds.Min_PT.x, Bounds.Min_PT.y, Bounds.Min_PT.z, Bounds.Max_PT.x, Bounds.Max_PT.y, Bounds.Max_PT.z);

                    no_of_Layers = 0;
                    double max_Layer = 0;

                    if (PML)
                    {
                        no_of_Layers = PML_LayerNo;
                        max_Layer = PML_MaxAtten;
                    }

                    double x_length = Bounds.X_Length(), y_length = Bounds.Y_Length(), z_length = Bounds.Z_Length();
                    Point MinPt = Bounds.Min_PT;
                    if (x_length < xmin)
                    {
                        x_length = xmin;
                        MinPt.x -= (xmin - x_length) / 2;
                    }
                    if (y_length < ymin)
                    {
                        y_length = ymin;
                        MinPt.y -= (ymin - y_length) / 2;
                    }
                    if (z_length < zmin)
                    {
                        z_length = zmin;
                        MinPt.z -= (zmin - z_length) / 2;
                    }

                    x_length += (no_of_Layers * 4 + 1) * dx;
                    y_length += (no_of_Layers * 4 + 1) * dydz;
                    z_length += (no_of_Layers * 4 + 1) * dydz;

                    //estimated distance between nodes
                    xDim = (int)Math.Ceiling(x_length / (dx));                                //set number of nodes in x direction
                    dx = x_length / xDim;                                                   //refined distance between nodes
                    yDim = (int)Math.Ceiling(y_length / dydz);                                //set number of nodes in y direction
                    dy = y_length / yDim;
                    zDim = (int)Math.Ceiling(z_length / dydz);                                //set number of nodes in z direction
                    dz = z_length / zDim;

                    dt = dy * rt2 / (Rm.Sound_speed(0));                           //set time step small enough to satisfy courrant condition
                    dxrt2 = dx * rt2;
                    dxrt3 = dx * rt3;

                    Dir = new Vector[13]{
                        new Vector(-1, 0, 0),
                        new Vector(0, -1, 0),
                        new Vector(0, 0, -1),

                        new Vector(-1/rt2,-1/rt2,double.Epsilon),
                        new Vector(1/rt2, -1/rt2,double.Epsilon),
                        new Vector(-1/rt2, double.Epsilon, 1/rt2),
                        new Vector(double.Epsilon, -1/rt2, 1/rt2),
                        new Vector(1/rt2, double.Epsilon, 1/rt2),
                        new Vector(double.Epsilon, 1/rt2, 1/rt2),

                        new Vector(-dx, -dy, dz),
                        new Vector(dx, -dy, dz),
                        new Vector(dx, dy, dz),
                        new Vector(-dx, dy, dz)
                    };

                    foreach (Vector V in Dir) V.Normalize();
                    
                    PFrame = new Node[xDim][][];// yDim, zDim];                               //pressure scalar field initialisation

                    List<Bound_Node_RDD> Bound = new List<Bound_Node_RDD>();

                    int xDimt = xDim;
                    int yDimt = yDim;
                    int zDimt = zDim;

                    MinPt.x -= 2 * dx * no_of_Layers;
                    MinPt.y -= 2 * dy * no_of_Layers;
                    MinPt.z -= 2 * dz * no_of_Layers;
                    Bounds = new AABB(MinPt, MinPt + new Point(x_length, y_length, z_length));

                    //System.Threading.Tasks.Parallel.For(0, xDim, (x) =>
                    for (int x = 0; x < PFrame.Length; x++)
                    {
                        int mod = x % 2;
                        PFrame[x] = new Node[(int)(Math.Floor((double)yDimt / 2) + yDimt % 2 * mod)][];
                        Random Rnd = new Random(x);
                        for (int y = 0; y < PFrame[x].Length; y++)
                        {
                            PFrame[x][y] = new Node[(int)(Math.Floor((double)zDimt / 2) + yDimt % 2 * mod)];
                            for (int z = 0; z < PFrame[x][y].Length; z++)
                            {
                                List<Environment.Material> abs;
                                List<Bound_Node.Boundary> BDir;
                                Point Loc = Acoustic_Compact_FDTD.RDD_Location(MinPt, x, y, z, dx, dy, dz); //new Point(MinPt.x + 2 * (((double)x - 0.5) * dx), MinPt.y + 2 * (((double)y + (0.5 - 0.5 * mod)) * dy), MinPt.z + 2 * (((double)z + (0.5 - 0.5 * mod)) * dz));
                                if (!Intersect_13Pt(Loc, SD.frequency, out BDir, out abs, ref Rnd))
                                {
                                    PFrame[x][y][z] = new RDD_Node(Loc);//, rho0, dt, dx, Rm.Sound_speed, new int[] { x, y, z });
                                }
                                else
                                {
                                    PFrame[x][y][z] =
                                        new Bound_Node_RDD(Loc, rho0, dt, dx, Rm.Sound_speed(0), new int[] { x, y, z }, BDir); // abs,
                                    Bound.Add(PFrame[x][y][z] as Bound_Node_RDD);
                                }
                            }
                        }
                    }//);

                    //Node.Attenuation = Math.Sqrt(Math.Pow(10, -.1 * Rm.AttenuationPureTone(Bounds.Min_PT, SD.frequency) * dt)); //PFrame[0][0][0].Pt, SD.frequency

                    bool failed = false;
                    //Make Mesh Templates:
                    //Build_Mesh_Sections();
                    //for (int x = 0; x < PFrame.Length; x++)
                    System.Threading.Tasks.Parallel.For(0, xDim, (X) =>
                    {
                        for (int y = 0; y < PFrame[X].Length; y++)
                        {
                            for (int z = 0; z < PFrame[X][y].Length; z++)
                            {
                                PFrame[X][y][z].Link_Nodes(ref PFrame, X, y, z);
                            }
                        }
                    });

                    yDim = PFrame[0].Length;
                    zDim = PFrame[0][0].Length;

                    for (int i = 0; i < Bound.Count; i++) if (Bound[i] != null) Bound[i].Complete_Boundary();
                    if (failed) return;

                    //Set up PML...
                    Layers = new Acoustic_Compact_FDTD.PML(no_of_Layers, max_Layer, PFrame, true);
                }

                public void RuntoCompletion()
                {
                    for (time_ms = 0; time_ms < tmax / 1000; time_ms += dt)
                    {
                        Increment();
                    }
                }

                //public void Build_Interp(ref int xDim, ref int yDim, ref int zDim)
                //{
                //    dx = Rm.Sound_speed(0) / fmax * .1;

                //    Bounds = new AABB(Rm.Min() - new Point(.05 * dx, .05 * dx, .05 * dx), Rm.Max() + new Point(.05 * dx, .05 * dx, .05 * dx));
                //    double x_length = Bounds.X_Length();
                //    double y_length = Bounds.Y_Length();
                //    double z_length = Bounds.Z_Length();

                //    //estimated distance between nodes
                //    xDim = (int)Math.Ceiling(x_length / dx);                                //set number of nodes in x direction
                //    dx = x_length / xDim;                                                       //refined distance between nodes
                //    yDim = (int)Math.Ceiling(y_length / dx);                                //set number of nodes in y direction
                //    double dy = y_length / yDim;
                //    zDim = (int)Math.Ceiling(z_length / dx);                              //set number of nodes in z direction
                //    double dz = z_length / zDim;

                //    double rt2 = Math.Sqrt(2);
                //    double rt3 = Math.Sqrt(3);
                //    dxrt2 = dx * rt2;
                //    dxrt3 = dx * rt3;

                //    Dir = new Vector[13]{
                //        new Vector(-1,double.Epsilon,double.Epsilon),
                //        new Vector(double.Epsilon,-1,double.Epsilon),
                //        new Vector(double.Epsilon,double.Epsilon,-1),

                //        new Vector(-1/rt2,-1/rt2,double.Epsilon),
                //        new Vector(1/rt2, -1/rt2,double.Epsilon),

                //        new Vector(-1/rt2, double.Epsilon, 1/rt2),
                //        new Vector(double.Epsilon, -1/rt2, 1/rt2),
                //        new Vector(1/rt2, double.Epsilon, 1/rt2),
                //        new Vector(double.Epsilon, 1/rt2, 1/rt2),

                //        new Vector(-1/rt3,-1/rt3,1/rt3),
                //        new Vector(1/rt3,-1/rt3,1/rt3),
                //        new Vector(1/rt3,1/rt3,1/rt3),
                //        new Vector(-1/rt3,1/rt3,1/rt3)
                //    };

                //    foreach (Vector V in Dir) V.Normalize();

                //    PFrame = new P_Node[xDim][][];// yDim, zDim];                               //pressure scalar field initialisation

                //    //System.Threading.Tasks.Parallel.For(0, xDim, (x) =>
                //    for (int x = 0; x < xDim; x++)
                //    {
                //        PFrame[x] = new P_Node[yDim][];
                //        Random Rnd = new Random(x);
                //        for (int y = 0; y < yDim; y++)
                //        {
                //            PFrame[x][y] = new P_Node[zDim];
                //            for (int z = 0; z < zDim; z++)
                //            {
                //                List<double[]> abs;
                //                List<Bound_Node.Boundary> BDir;
                //                Point Loc = new Point(Bounds.Min_PT.x + (((double)x + 0.5) * dx), Bounds.Min_PT.y + (((double)y + 0.5) * dy), Bounds.Min_PT.z + (((double)z + 0.5) * dz));
                //                if (!Intersect_26Pt(Loc, out BDir, out abs, ref Rnd))
                //                {
                //                    PFrame[x][y][z] = new P_Node(Loc);//, rho0, dt, dx, Rm.Sound_speed, new int[] { x, y, z });
                //                }
                //                else
                //                {
                //                    PFrame[x][y][z] =
                //                        new Bound_Node_IWB(Loc, rho0, dt, dx, Rm.Sound_speed(0), new int[] { x, y, z }, abs, BDir);
                //                }
                //            }
                //        }
                //    }//);

                //    bool failed = false;
                //    //Make Mesh Templates:
                //    Build_Mesh_Sections();
                //    for (int x = 0; x < xDim; x++)
                //    //System.Threading.Tasks.Parallel.For(0, xDim, (x) =>
                //    {
                //        for (int y = 0; y < yDim; y++)
                //        {
                //            for (int z = 0; z < zDim; z++)
                //            {
                //                try
                //                {
                //                    PFrame[x][y][z].Link_Nodes(ref PFrame, x, y, z);
                //                }
                //                catch
                //                {
                //                    //Display faulty voxels in a display conduit here...
                //                    UI.CellConduit.Instance.Add(PFrame[x][y][z], x, y, z, dx, Bounds.Min_PT);
                //                    failed = true;
                //                }
                //            }
                //        }
                //    }//);

                //    //dt = dx / (Math.Sqrt(1.0 / 3.0) * (Rm.Sound_speed(0)));                           //set time step small enough to satisfy courrant condition
                //    //dt = dx / (Rm.Sound_speed(0));                           //set time step small enough to satisfy courrant condition
                //    //dt = dx / (Math.Sqrt(.75) * (Rm.Sound_speed));                           //set time step small enough to satisfy courrant condition
                //    dt = Math.Sqrt(dx * dx + dy * dy + dz * dz) / Rm.Sound_speed(0);
                //    if (failed) return;

                //    SD.Connect_Grid_Freefield(PFrame, Bounds, dx, dy, dz, tmax, dt, 0);
                //    Mic.Connect_Grid_Freefield(PFrame, Bounds, dx, tmax, dt, 0);
                //}

                public double SampleFrequency
                {
                    get { return 1f / dt; }
                }

                public void reset(double frequency, Signal_Driver_Compact.Signal_Type s)
                {
                    SD.reset(frequency, s);
                    Mic.reset();
                    time_ms = 0;
                    n = 0;
                    for (int X = 0; X < PFrame.Length; X++)
                        for (int Y = 0; Y < PFrame[X].Length; Y++)
                            for (int Z = 0; Z < PFrame[X][Y].Length; Z++)
                                PFrame[X][Y][Z].reset();
                }

                public double Increment()
                {
                    //Rhino.RhinoApp.CommandPrompt = string.Format("Running {0} Hz., {1} ms.", SD.frequency, Math.Round(time_ms * 1000));
                    List<Hare.Geometry.Point> Pts = new List<Hare.Geometry.Point>();
                    List<double> Pressure = new List<double>();

                    System.Threading.Tasks.Parallel.For(0, threadct * 10, proc_id =>
                    //for (int p = 0; p < threadct; p++)
                    //{
                    {
                        for (int x = (int)Math.Floor((double)((int)proc_id * PFrame.Length) / (10 * threadct)); x < Math.Floor((double)(((int)proc_id + 1) * PFrame.Length) / (10 * threadct)); x++)
                        {
                            for (int y = 0; y < PFrame[x].Length; y++)
                            {
                                for (int z = 0; z < PFrame[x][y].Length; z++)
                                {
                                    PFrame[x][y][z].UpdateP();
                                }
                            }
                        }
                    });

                    //for (int p = 0; p < threadct; p++)
                    //{
                    //inc2.Add(new System.Threading.Thread(
                    System.Threading.Tasks.Parallel.For(0, threadct * 10, proc_id =>
                        {
                            for (int x = (int)Math.Floor((double)((int)proc_id * PFrame.Length) / (threadct * 10)); x < Math.Floor((double)(((int)proc_id + 1) * PFrame.Length) / (threadct * 10)); x++)
                            {
                                {
                                    for (int y = 0; y < PFrame[x].Length; y++)
                                    {
                                        for (int z = 0; z < PFrame[x][y].Length; z++)
                                        {
                                            PFrame[x][y][z].UpdateT();
                                        }
                                    }
                                }
                            }
                        }
                    );

                    //for (int i = 0; i < threadct; i++) inc1[i].Start(i);
                    //foreach (System.Threading.Thread i in inc1) i.Join();
                    //for (int i = 0; i < threadct; i++) inc2[i].Start(i);
                    //foreach (System.Threading.Thread i in inc2) i.Join();

                    Layers.Attenuate();

                    SD.Drive(n);
                    Mic.Record(n);
                    n += 1;
                    time_ms = n * dt;
                    return time_ms;
                }

                /////////////////
                //Display Methods : Get Points and Pressure for display output
                /////////////////
                public void Pressure_Points(ref List<List<Hare.Geometry.Point>> Pts, ref List<List<double>> Pressure, int[] X, int[] Y, int[] Z, double Low_P, bool Volume, bool Vectored, bool Colored)
                {
                    Pts = new List<List<Hare.Geometry.Point>>();
                    Pressure = new List<List<double>>();
                    if (Volume)
                    {
                        for (int x = 0; x < xDim; x++)
                        {
                            List<double> P = new List<double>();
                            List<Hare.Geometry.Point> PtList = new List<Hare.Geometry.Point>();
                            //System.Threading.Tasks.Parallel.For(0, yDim, (y) =>
                            for (int y = 0; y < yDim; y++)
                            {
                                for (int z = 0; z < zDim; z++)
                                {
                                    if (PFrame[x][y][z].P < Low_P) continue;
                                    if (Colored) P.Add(PFrame[x][y][z].P);
                                    //if (Vectored)
                                    //{
                                    //    PtList.Add((new Rhino.Geometry.Point3d(PFrame[x][y][z].Pt.x, PFrame[x][y][z].Pt.y, PFrame[x][y][z].Pt.z) + PFrame[x][y][z].VelocityDirection() * 3 * dx * (Magnitude ? PFrame[x][y][z].P: PFrame[x][y][z].P) ));
                                    //}
                                    //else
                                    //{
                                    PtList.Add(Acoustic_Compact_FDTD.RDD_Location(Bounds.Min_PT, x, y, z, dx, dy, dz));//new Rhino.Geometry.Point3d(PFrame[x][y][z].Pt.x, PFrame[x][y][z].Pt.y, PFrame[x][y][z].Pt.z)
                                    //}
                                }
                            }//);
                            Pressure.Add(P);
                            Pts.Add(PtList);
                        }
                        return;
                    }
                    if (X != null || X.Length > 0)
                    {
                        foreach (int x in X)
                        {
                            List<double> P = new List<double>();
                            List<Hare.Geometry.Point> PtList = new List<Hare.Geometry.Point>();
                            //System.Threading.Tasks.Parallel.For(0, yDim, (y) =>
                            for (int y = 0; y < PFrame[x].Length; y++)
                            {
                                for (int z = 0; z < PFrame[x][y % 2].Length; z++)
                                {
                                    lock (Pts)
                                    {
                                        if (Colored) P.Add(PFrame[x][y][z].P);
                                        //if (Vectored)
                                        //{
                                        //    PtList.Add((new Rhino.Geometry.Point3d(PFrame[x][y][z].Pt.x, PFrame[x][y][z].Pt.y, PFrame[x][y][z].Pt.z) + PFrame[x][y][z].VelocityDirection() * PFrame[x][y][z].P * 3 * dx));
                                        //}
                                        //else
                                        //{
                                        PtList.Add(Acoustic_Compact_FDTD.RDD_Location(Bounds.Min_PT, x, y, z, dx, dy, dz));//new Rhino.Geometry.Point3d(PFrame[x][y][z].Pt.x, PFrame[x][y][z].Pt.y, PFrame[x][y][z].Pt.z)
                                        //}
                                    }
                                }
                            }//);
                            Pressure.Add(P);
                            Pts.Add(PtList);
                        }
                    }
                    if (Y != null || Y.Length > 0)
                    {
                        foreach (int y in Y)
                        {
                            List<double> P = new List<double>();
                            List<Hare.Geometry.Point> PtList = new List<Hare.Geometry.Point>();
                            //System.Threading.Tasks.Parallel.For(0, xDim, (x) =>
                            for (int x = 0; x < PFrame.Length; x++)
                            {
                                for (int z = 0; z < PFrame[x][y % 2].Length; z++)
                                {
                                    lock (Pts)
                                    {
                                        if (Colored) P.Add(PFrame[x][y][z].P);
                                        //if (Vectored)
                                        //{
                                        //    PtList.Add((new Rhino.Geometry.Point3d(PFrame[x][y][z].Pt.x, PFrame[x][y][z].Pt.y, PFrame[x][y][z].Pt.z) + PFrame[x][y][z].VelocityDirection() * PFrame[x][y][z].P * 3 * dx));
                                        //}
                                        //else
                                        //{
                                        PtList.Add(Acoustic_Compact_FDTD.RDD_Location(Bounds.Min_PT, x, y, z, dx, dy, dz));//new Rhino.Geometry.Point3d(PFrame[x][y][z].Pt.x, PFrame[x][y][z].Pt.y, PFrame[x][y][z].Pt.z)
                                        //}
                                    }
                                }
                            }//);
                            Pressure.Add(P);
                            Pts.Add(PtList);
                        }
                    }
                    if (Z != null || Z.Length > 0)
                    {
                        foreach (int z in Z)
                        {
                            List<double> P = new List<double>();
                            List<Hare.Geometry.Point> PtList = new List<Hare.Geometry.Point>();
                            //System.Threading.Tasks.Parallel.For(0, xDim, (x) =>
                            for (int x = 0; x < PFrame.Length; x++)
                            {
                                for (int y = 0; y < PFrame[x].Length; y++)
                                {
                                    lock (Pts)
                                    {
                                        try { if (Colored) P.Add(PFrame[x][y][z].P); }
                                        catch 
                                        {
                                            if (x > PFrame.Length - 1) x = PFrame.Length - 1;
                                            else if (y > PFrame[x].Length - 1)  y = PFrame[x].Length - 1;
                                            if (Colored) P.Add(PFrame[x][y][z].P);
                                        }
                                        //if (Vectored)
                                        //{
                                        //    PtList.Add((new Rhino.Geometry.Point3d(PFrame[x][y][z].Pt.x, PFrame[x][y][z].Pt.y, PFrame[x][y][z].Pt.z) + PFrame[x][y][z].VelocityDirection() * PFrame[x][y][z].P * 3 * dx));
                                        //}
                                        //else
                                        //{
                                        PtList.Add(Acoustic_Compact_FDTD.RDD_Location(Bounds.Min_PT, x, y, z, dx, dy, dz));//new Rhino.Geometry.Point3d(PFrame[x][y][z].Pt.x, PFrame[x][y][z].Pt.y, PFrame[x][y][z].Pt.z)
                                        //}
                                    }
                                }
                            }//);
                            Pressure.Add(P);
                            Pts.Add(PtList);
                        }
                    }
                }

                double dxrt2;
                double dxrt3;

                //public bool Intersect_26Pt(Point Center, out List<Bound_Node_IWB.Boundary> Bout, out List<double[]> alpha, ref Random Rnd)
                //{
                //    Bout = new List<Bound_Node.Boundary>();
                //    alpha = new List<double[]>();

                //    Center += new Point((Rnd.NextDouble() - .5) * 1E-6, (Rnd.NextDouble() - .5) * 1E-6, (Rnd.NextDouble() - .5) * 1E-6);

                //    //double u, v;

                //    X_Event XPt = new X_Event();
                //    //new Vector(-1,0,0),
                //    if (Rm.shoot(new Ray(Center, Dir[0], 0, Rnd.Next()), 0, out XPt) && XPt.t < dx)
                //    {
                //        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center + Dir[0] * dx)));
                //        Bout.Add(Bound_Node.Boundary.AXNeg);
                //        //TODO: Intelligently Assign Absorption
                //        alpha.Add(Rm.AbsorptionValue[XPt.Poly_id].Coefficient_A_Broad());
                //    }
                //    XPt = new X_Event();
                //    if (Rm.shoot(new Ray(Center, Dir[0] * -1, 0, Rnd.Next()), 0, out XPt) && XPt.t < dx)
                //    {
                //        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center - Dir[0] * dx)));
                //        Bout.Add(Bound_Node.Boundary.AXPos);
                //        //TODO: Intelligently Assign Absorption
                //        alpha.Add(Rm.AbsorptionValue[XPt.Poly_id].Coefficient_A_Broad());
                //    }
                //    //new Vector(0,-1,0),
                //    XPt = new X_Event();
                //    if (Rm.shoot(new Ray(Center, Dir[1], 0, Rnd.Next()), 0, out XPt) && XPt.t < dx)
                //    {
                //        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center + Dir[1] * dx)));
                //        Bout.Add(Bound_Node.Boundary.AYNeg);
                //        alpha.Add(Rm.AbsorptionValue[XPt.Poly_id].Coefficient_A_Broad());
                //    }
                //    XPt = new X_Event();
                //    if (Rm.shoot(new Ray(Center, Dir[1] * -1, 0, Rnd.Next()), 0, out XPt) && XPt.t < dx)
                //    {
                //        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center - Dir[1] * dx)));
                //        Bout.Add(Bound_Node.Boundary.AYPos);
                //        alpha.Add(Rm.AbsorptionValue[XPt.Poly_id].Coefficient_A_Broad());
                //    }
                //    //new Vector(0,0,-1),
                //    XPt = new X_Event();
                //    if (Rm.shoot(new Ray(Center, Dir[2], 0, Rnd.Next()), 0, out XPt) && XPt.t < dx)
                //    {
                //        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center + Dir[2] * dx)));
                //        Bout.Add(Bound_Node.Boundary.AZNeg);
                //        alpha.Add(Rm.AbsorptionValue[XPt.Poly_id].Coefficient_A_Broad());
                //    }
                //    XPt = new X_Event();
                //    if (Rm.shoot(new Ray(Center, Dir[2] * -1, 0, Rnd.Next()), 0, out XPt) && XPt.t < dx)
                //    {
                //        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center - Dir[2] * dx)));
                //        Bout.Add(Bound_Node.Boundary.AZPos);
                //        alpha.Add(Rm.AbsorptionValue[XPt.Poly_id].Coefficient_A_Broad());
                //    }
                //    //new Vector(-1/rt2,-1/rt2,0),
                //    XPt = new X_Event();
                //    if (Rm.shoot(new Ray(Center, Dir[3], 0, Rnd.Next()), 0, out XPt) && XPt.t < dxrt2)
                //    {
                //        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center + Dir[3] * dx)));
                //        Bout.Add(Bound_Node.Boundary.SDXNegYNeg);
                //        alpha.Add(Rm.AbsorptionValue[XPt.Poly_id].Coefficient_A_Broad());
                //    }
                //    XPt = new X_Event();
                //    if (Rm.shoot(new Ray(Center, Dir[3] * -1, 0, Rnd.Next()), 0, out XPt) && XPt.t < dxrt2)
                //    {
                //        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center - Dir[3] * dx)));
                //        Bout.Add(Bound_Node.Boundary.SDXPosYPos);
                //        alpha.Add(Rm.AbsorptionValue[XPt.Poly_id].Coefficient_A_Broad());
                //    }
                //    //new Vector(1/rt2, -1/rt2,0),
                //    XPt = new X_Event();
                //    if (Rm.shoot(new Ray(Center, Dir[4], 0, Rnd.Next()), 0, out XPt) && XPt.t < dxrt2)
                //    {
                //        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center + Dir[4] * dx)));
                //        Bout.Add(Bound_Node.Boundary.SDXPosYNeg);
                //        alpha.Add(Rm.AbsorptionValue[XPt.Poly_id].Coefficient_A_Broad());
                //    }
                //    XPt = new X_Event();
                //    if (Rm.shoot(new Ray(Center, Dir[4] * -1, 0, Rnd.Next()), 0, out XPt) && XPt.t < dxrt2)
                //    {
                //        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center - Dir[4] * dx)));
                //        Bout.Add(Bound_Node.Boundary.SDXNegYPos);
                //        alpha.Add(Rm.AbsorptionValue[XPt.Poly_id].Coefficient_A_Broad());
                //    }
                //    //new Vector(-1/rt2, 0, 1/rt2),
                //    XPt = new X_Event();
                //    if (Rm.shoot(new Ray(Center, Dir[5], 0, Rnd.Next()), 0, out XPt) && XPt.t < dxrt2)
                //    {
                //        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center + Dir[5] * dx)));
                //        Bout.Add(Bound_Node.Boundary.SDXNegZPos);
                //        alpha.Add(Rm.AbsorptionValue[XPt.Poly_id].Coefficient_A_Broad());
                //    }
                //    XPt = new X_Event();
                //    if (Rm.shoot(new Ray(Center, Dir[5] * -1, 0, Rnd.Next()), 0, out XPt) && XPt.t < dxrt2)
                //    {
                //        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center - Dir[5] * dx)));
                //        Bout.Add(Bound_Node.Boundary.SDXPosZNeg);
                //        alpha.Add(Rm.AbsorptionValue[XPt.Poly_id].Coefficient_A_Broad());
                //    }
                //    //new Vector(0, -1/rt2, 1/rt2),
                //    XPt = new X_Event();
                //    if (Rm.shoot(new Ray(Center, Dir[6], 0, Rnd.Next()), 0, out XPt) && XPt.t < dxrt2)
                //    {
                //        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center + Dir[6] * dx)));
                //        Bout.Add(Bound_Node.Boundary.SDYNegZPos);
                //        alpha.Add(Rm.AbsorptionValue[XPt.Poly_id].Coefficient_A_Broad());
                //    }
                //    XPt = new X_Event();
                //    if (Rm.shoot(new Ray(Center, Dir[6] * -1, 0, Rnd.Next()), 0, out XPt) && XPt.t < dxrt2)
                //    {
                //        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center - Dir[6] * dx)));
                //        Bout.Add(Bound_Node.Boundary.SDYPosZNeg);
                //        alpha.Add(Rm.AbsorptionValue[XPt.Poly_id].Coefficient_A_Broad());
                //    }
                //    //new Vector(1/rt2, 0, 1/rt2),
                //    XPt = new X_Event();
                //    if (Rm.shoot(new Ray(Center, Dir[7], 0, Rnd.Next()), 0, out XPt) && XPt.t < dxrt2)
                //    {
                //        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center + Dir[7] * dx)));
                //        Bout.Add(Bound_Node.Boundary.SDXPosZPos);
                //        alpha.Add(Rm.AbsorptionValue[XPt.Poly_id].Coefficient_A_Broad());
                //    }
                //    XPt = new X_Event();
                //    if (Rm.shoot(new Ray(Center, Dir[7] * -1, 0, Rnd.Next()), 0, out XPt) && XPt.t < dxrt2)
                //    {
                //        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center - Dir[7] * dx)));
                //        Bout.Add(Bound_Node.Boundary.SDXNegZNeg);
                //        alpha.Add(Rm.AbsorptionValue[XPt.Poly_id].Coefficient_A_Broad());
                //    }
                //    //new Vector(0, 1/rt2, 1/rt2),
                //    XPt = new X_Event();
                //    if (Rm.shoot(new Ray(Center, Dir[8], 0, Rnd.Next()), 0, out XPt) && XPt.t < dxrt2)
                //    {
                //        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center + Dir[8] * dx)));
                //        Bout.Add(Bound_Node.Boundary.SDYPosZPos);
                //        alpha.Add(Rm.AbsorptionValue[XPt.Poly_id].Coefficient_A_Broad());
                //    }
                //    XPt = new X_Event();
                //    if (Rm.shoot(new Ray(Center, Dir[8] * -1, 0, Rnd.Next()), 0, out XPt) && XPt.t < dxrt2)
                //    {
                //        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center - Dir[8] * dx)));
                //        Bout.Add(Bound_Node.Boundary.SDYNegZNeg);
                //        alpha.Add(Rm.AbsorptionValue[XPt.Poly_id].Coefficient_A_Broad());
                //    }
                //    //new Vector(-1/rt3,-1/rt3,1/rt3),
                //    XPt = new X_Event();
                //    if (Rm.shoot(new Ray(Center, Dir[9], 0, Rnd.Next()), 0, out XPt) && XPt.t < dxrt2)
                //    {
                //        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center + Dir[9] * dx)));
                //        Bout.Add(Bound_Node.Boundary.DXNegYNegZPos);
                //        alpha.Add(Rm.AbsorptionValue[XPt.Poly_id].Coefficient_A_Broad());
                //    }
                //    XPt = new X_Event();
                //    if (Rm.shoot(new Ray(Center, Dir[9] * -1, 0, Rnd.Next()), 0, out XPt) && XPt.t < dxrt2)
                //    {
                //        Bout.Add(Bound_Node.Boundary.DXPosYPosZNeg);
                //        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center - Dir[9] * dx)));
                //        alpha.Add(Rm.AbsorptionValue[XPt.Poly_id].Coefficient_A_Broad());
                //    }
                //    //new Vector(1/rt3,-1/rt3,1/rt3),
                //    XPt = new X_Event();
                //    if (Rm.shoot(new Ray(Center, Dir[10], 0, Rnd.Next()), 0, out XPt) && XPt.t < dxrt2)
                //    {
                //        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center + Dir[10] * dx)));
                //        Bout.Add(Bound_Node.Boundary.DXPosYNegZPos);
                //        alpha.Add(Rm.AbsorptionValue[XPt.Poly_id].Coefficient_A_Broad());
                //    }
                //    XPt = new X_Event();
                //    if (Rm.shoot(new Ray(Center, Dir[10] * -1, 0, Rnd.Next()), 0, out XPt) && XPt.t < dxrt2)
                //    {
                //        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center - Dir[10] * dx)));
                //        Bout.Add(Bound_Node.Boundary.DXNegYPosZNeg);
                //        alpha.Add(Rm.AbsorptionValue[XPt.Poly_id].Coefficient_A_Broad());
                //    }
                //    //new Vector(1/rt3,1/rt3,1/rt3),
                //    XPt = new X_Event();
                //    if (Rm.shoot(new Ray(Center + new Point(0, 1E-6, -1E-6), Dir[11], 0, Rnd.Next()), 0, out XPt) && XPt.t < dxrt2)
                //    {
                //        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center + Dir[11] * dx)));
                //        Bout.Add(Bound_Node.Boundary.DXPosYPosZPos);
                //        alpha.Add(Rm.AbsorptionValue[XPt.Poly_id].Coefficient_A_Broad());
                //    }
                //    XPt = new X_Event();
                //    if (Rm.shoot(new Ray(Center, Dir[11] * -1, 0, Rnd.Next()), 0, out XPt) && XPt.t < dxrt2)
                //    {
                //        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center - Dir[11] * dx)));
                //        Bout.Add(Bound_Node.Boundary.DXNegYNegZNeg);
                //        alpha.Add(Rm.AbsorptionValue[XPt.Poly_id].Coefficient_A_Broad());
                //    }
                //    //new Vector(-1/rt3,1/rt3,1/rt3)
                //    XPt = new X_Event();
                //    if (Rm.shoot(new Ray(Center, Dir[12], 0, Rnd.Next()), 0, out XPt) && XPt.t < dxrt2)
                //    {
                //        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center + Dir[12] * dx)));
                //        Bout.Add(Bound_Node.Boundary.DXNegYPosZPos);
                //        alpha.Add(Rm.AbsorptionValue[XPt.Poly_id].Coefficient_A_Broad());
                //    }
                //    XPt = new X_Event();
                //    if (Rm.shoot(new Ray(Center, Dir[12] * -1, 0, Rnd.Next()), 0, out XPt) && XPt.t < dxrt2)
                //    {
                //        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center - Dir[12] * dx)));
                //        Bout.Add(Bound_Node.Boundary.DXPosYNegZNeg);
                //        alpha.Add(Rm.AbsorptionValue[XPt.Poly_id].Coefficient_A_Broad());
                //    }

                //    return Bout.Count > 0;
                //}

                public bool Intersect_13Pt(Point Center, double frequency, out List<Bound_Node.Boundary> Bout, out List<Environment.Material> mat, ref Random Rnd)
                {
                    Bout = new List<Bound_Node.Boundary>();
                    mat = new List<Environment.Material>();

                    Center += new Point((Rnd.NextDouble() - .5) * 1E-6, (Rnd.NextDouble() - .5) * 1E-6, (Rnd.NextDouble() - .5) * 1E-6);

                    X_Event XPt = new X_Event();

                    double dx2 = 2 * dy + double.Epsilon;

                    XPt = new X_Event();
                    //new Vector(0, -1, 0)
                    if (Rm.shoot(new Ray(Center, -1 * Dir[1], 0, Rnd.Next()), 0, out XPt) && XPt.t < dx2)
                    {
                        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center - Dir[1] * 2 * dx)));
                        Bout.Add(Bound_Node.Boundary.AYPos);
                        mat.Add(Rm.AbsorptionValue[XPt.Poly_id]);
                    }
                    XPt = new X_Event();
                    if (Rm.shoot(new Ray(Center, Dir[1], 0, Rnd.Next()), 0, out XPt) && XPt.t < dx2)
                    {
                        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center + Dir[1] * 2 * dx)));
                        Bout.Add(Bound_Node.Boundary.AYNeg);
                        mat.Add(Rm.AbsorptionValue[XPt.Poly_id]);
                    }

                    XPt = new X_Event();
                    //new Vector(0, 0, -1),
                    if (Rm.shoot(new Ray(Center, -1 * Dir[2], 0, Rnd.Next()), 0, out XPt) && XPt.t < dx2)
                    {
                        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center - Dir[2] * 2 * dx)));
                        Bout.Add(Bound_Node.Boundary.AZPos);
                        mat.Add(Rm.AbsorptionValue[XPt.Poly_id]);
                    }
                    XPt = new X_Event();
                    if (Rm.shoot(new Ray(Center, Dir[2], 0, Rnd.Next()), 0, out XPt) && XPt.t < dx2)
                    {
                        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center + Dir[2] * 2 * dx)));
                        Bout.Add(Bound_Node.Boundary.AZNeg);
                        mat.Add(Rm.AbsorptionValue[XPt.Poly_id]);
                    }

                    XPt = new X_Event();
                    //new Vector(-dx, -dy, dz),
                    if (Rm.shoot(new Ray(Center, Dir[9], 0, Rnd.Next()), 0, out XPt) && XPt.t < dx2)
                    {
                        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center + Dir[9] * 2 * dx)));
                        Bout.Add(Bound_Node.Boundary.DXNegYNegZPos);
                        mat.Add(Rm.AbsorptionValue[XPt.Poly_id]);
                    }
                    XPt = new X_Event();
                    if (Rm.shoot(new Ray(Center, Dir[9] * -1, 0, Rnd.Next()), 0, out XPt) && XPt.t < dx2)
                    {
                        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center - Dir[9] * 2 * dx)));
                        Bout.Add(Bound_Node.Boundary.DXPosYPosZNeg);
                        mat.Add(Rm.AbsorptionValue[XPt.Poly_id]);
                    }

                    XPt = new X_Event();
                    //new Vector(dx, -dy, dz),
                    if (Rm.shoot(new Ray(Center, Dir[10], 0, Rnd.Next()), 0, out XPt) && XPt.t < dx2)
                    {
                        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center + Dir[10] * 2 * dx)));
                        Bout.Add(Bound_Node.Boundary.DXPosYNegZPos);
                        mat.Add(Rm.AbsorptionValue[XPt.Poly_id]);
                    }
                    XPt = new X_Event();
                    if (Rm.shoot(new Ray(Center, Dir[10] * -1, 0, Rnd.Next()), 0, out XPt) && XPt.t < dx2)
                    {
                        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center - Dir[10] * 2 * dx)));
                        Bout.Add(Bound_Node.Boundary.DXNegYPosZNeg);
                        mat.Add(Rm.AbsorptionValue[XPt.Poly_id]);
                    }

                    XPt = new X_Event();
                    //new Vector(dx, dy, dz),
                    if (Rm.shoot(new Ray(Center, Dir[11], 0, Rnd.Next()), 0, out XPt) && XPt.t < dx2)
                    {
                        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center + Dir[11] * 2 * dx)));
                        Bout.Add(Bound_Node.Boundary.DXPosYPosZPos);
                        mat.Add(Rm.AbsorptionValue[XPt.Poly_id]);
                    }
                    XPt = new X_Event();
                    if (Rm.shoot(new Ray(Center, Dir[11] * -1, 0, Rnd.Next()), 0, out XPt) && XPt.t < dx2)
                    {
                        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center - Dir[11] * 2 * dx)));
                        Bout.Add(Bound_Node.Boundary.DXNegYNegZNeg);
                        mat.Add(Rm.AbsorptionValue[XPt.Poly_id]);
                    }

                    XPt = new X_Event();
                    //new Vector(-dx, dy, dz)
                    if (Rm.shoot(new Ray(Center, Dir[12], 0, Rnd.Next()), 0, out XPt) && XPt.t < dx2)
                    {
                        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center + Dir[12] * 2 * dx)));
                        Bout.Add(Bound_Node.Boundary.DXNegYPosZPos);
                        mat.Add(Rm.AbsorptionValue[XPt.Poly_id]);
                    }
                    XPt = new X_Event();
                    if (Rm.shoot(new Ray(Center, Dir[12] * -1, 0, Rnd.Next()), 0, out XPt) && XPt.t < dx2)
                    {
                        //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(new Rhino.Geometry.Line(Utilities.PachTools.HPttoRPt(Center), Utilities.PachTools.HPttoRPt(Center - Dir[12] * 2 * dx)));
                        Bout.Add(Bound_Node.Boundary.DXPosYNegZNeg);
                        mat.Add(Rm.AbsorptionValue[XPt.Poly_id]);
                    }
                    return Bout.Count > 0;
                }

                public override string ProgressMsg()
                {
                    throw new NotImplementedException();
                }

                public override System.Threading.ThreadState ThreadState()
                {
                    throw new NotImplementedException();
                }

                public override void Combine_ThreadLocal_Results()
                {
                    throw new NotImplementedException();
                }

                public override void Begin()
                {
                    throw new NotImplementedException();
                }

                public override string Sim_Type()
                {
                    return "Finite Difference Time Domain";
                }

                public double P(int x, int y, int z)
                {
                    return this.PFrame[x][y][z].P;
                }

                public abstract class Node
                {
                    //public Point Pt;
                    public double Pnf, Pn, Pn_1;

                    public abstract void Link_Nodes(ref Node[][][] Frame, int x, int y, int z);
                    public abstract void UpdateP();
                    public static double Attenuation = 1;

                    public Node(Point loc)
                    {
                        //Pt = loc;
                    }

                    public virtual void UpdateT()
                    {
                        Pn_1 = Pn;
                        Pn = Pnf * Attenuation;
                    }

                    public void reset()
                    {
                        Pnf = 0; Pn = 0; Pn_1 = 0;
                    }

                    /// <summary>
                    /// Call after UpdateT, to levy additional attenuation...
                    /// </summary>
                    /// <param name="coefficient"></param>
                    public void Attenuate(double coefficient)
                    {
                        Pn *= coefficient;
                    }

                    //public Rhino.Geometry.Vector3d VelocityDirection()
                    //{
                    //    Rhino.Geometry.Vector3d V = new Rhino.Geometry.Vector3d(X, Y, Z);
                    //    V.Unitize();
                    //    return V;
                    //}

                    public double P
                    {
                        get
                        {
                            return Pn;
                        }
                    }
                }

                public class P_Node : Node
                {
                    protected double X, Y, Z;
                    protected Node Xpos_Link;
                    protected Node Ypos_Link;
                    protected Node Zpos_Link;
                    protected Node Xneg_Link;
                    protected Node Yneg_Link;
                    protected Node Zneg_Link;

                    protected Node[] Links2;
                    protected Node[] Links3;

                    public P_Node(Point loc)
                        : base(loc)
                    { }

                    public override void Link_Nodes(ref Node[][][] Frame, int x, int y, int z)
                    {
                        int xdim = Frame.Length - 1;
                        int ydim = Frame[0].Length - 1;
                        int zdim = Frame[0][0].Length - 1;

                        if (x < xdim)
                        {
                            Xpos_Link = Frame[x + 1][y][z] as P_Node;
                        }
                        else Xpos_Link = new Null_Node();
                        //FrameX[x + 1, y, z].id.AddRange(id);
                        if (y < ydim)
                        {
                            Ypos_Link = Frame[x][y + 1][z] as P_Node;
                        }
                        else Ypos_Link = new Null_Node();
                        //FrameY[x, y + 1, z].id.AddRange(id);
                        if (z < zdim)
                        {
                            Zpos_Link = Frame[x][y][z + 1] as P_Node;
                        }
                        else Zpos_Link = new Null_Node();
                        //FrameZ[x, y, z + 1].id.AddRange(id);
                        if (x > 0)
                        {
                            Xneg_Link = Frame[x - 1][y][z] as P_Node;
                        }
                        else Xneg_Link = new Null_Node();
                        //FrameX[x, y, z].id.AddRange(id);
                        if (y > 0)
                        {
                            Yneg_Link = Frame[x][y - 1][z] as P_Node;
                        }
                        else Yneg_Link = new Null_Node();
                        //FrameY[x, y, z].id.AddRange(id);
                        if (z > 0)
                        {
                            Zneg_Link = Frame[x][y][z - 1] as P_Node;
                        }
                        else Zneg_Link = new Null_Node();
                        //FrameZ[x, y, z].id.AddRange(id);

                        Links2 = new P_Node[12];
                        Links3 = new P_Node[8];

                        if (x < xdim && y < ydim) Links2[0] = Frame[x + 1][y + 1][z] as P_Node; else Links2[0] = new Null_Node();
                        if (x > 0 && y > 0) Links2[1] = Frame[x - 1][y - 1][z] as P_Node; else Links2[1] = new Null_Node();
                        if (x < xdim && y > 0) Links2[2] = Frame[x + 1][y - 1][z] as P_Node; else Links2[2] = new Null_Node();
                        if (x > 0 && y < ydim) Links2[3] = Frame[x - 1][y + 1][z] as P_Node; else Links2[3] = new Null_Node();
                        if (x < xdim && z < zdim) Links2[4] = Frame[x + 1][y][z + 1] as P_Node; else Links2[4] = new Null_Node();
                        if (x > 0 && z > 0) Links2[5] = Frame[x - 1][y][z - 1] as P_Node; else Links2[5] = new Null_Node();
                        if (x > 0 && z < zdim) Links2[6] = Frame[x - 1][y][z + 1] as P_Node; else Links2[6] = new Null_Node();
                        if (x < xdim && z > 0) Links2[7] = Frame[x + 1][y][z - 1] as P_Node; else Links2[7] = new Null_Node();
                        if (y < ydim && z < zdim) Links2[8] = Frame[x][y + 1][z + 1] as P_Node; else Links2[8] = new Null_Node();
                        if (y > 0 && z > 0) Links2[9] = Frame[x][y - 1][z - 1] as P_Node; else Links2[9] = new Null_Node();
                        if (y > 0 && z < zdim) Links2[10] = Frame[x][y - 1][z + 1] as P_Node; else Links2[10] = new Null_Node();
                        if (y < ydim && z > 0) Links2[11] = Frame[x][y + 1][z - 1] as P_Node; else Links2[11] = new Null_Node();

                        if (x < xdim && y < ydim && z < zdim) Links3[0] = Frame[x + 1][y + 1][z + 1] as P_Node; else Links3[0] = new Null_Node();
                        if (x > 0 && y > 0 && z > 0) Links3[1] = Frame[x - 1][y - 1][z - 1] as P_Node; else Links3[1] = new Null_Node();
                        if (x > 0 && y > 0 && z < zdim) Links3[2] = Frame[x - 1][y - 1][z + 1] as P_Node; else Links3[2] = new Null_Node();
                        if (x < xdim && y < ydim && z > 0) Links3[3] = Frame[x + 1][y + 1][z - 1] as P_Node; else Links3[3] = new Null_Node();
                        if (x > 0 && y < ydim && z < zdim) Links3[4] = Frame[x - 1][y + 1][z + 1] as P_Node; else Links3[4] = new Null_Node();
                        if (x < xdim && y > 0 && z > 0) Links3[5] = Frame[x + 1][y - 1][z - 1] as P_Node; else Links3[5] = new Null_Node();
                        if (x > 0 && y < ydim && z > 0) Links3[6] = Frame[x - 1][y + 1][z - 1] as P_Node; else Links3[6] = new Null_Node();
                        if (x < xdim && y > 0 && z < zdim) Links3[7] = Frame[x + 1][y - 1][z + 1] as P_Node; else Links3[7] = new Null_Node();
                    }

                    //public virtual void UpdateSLF()
                    //{
                    //    X = Xpos_Link.P + Xneg_Link.P;
                    //    Y = Ypos_Link.P + Yneg_Link.P;
                    //    Z = Zpos_Link.P + Zneg_Link.P;
                    //    Pnf = Courrant2 * (X + Y + Z) - Pn_1;
                    //}

                    //public virtual void UpdateOCTA()
                    //{
                    //    double p3 = 0;
                    //    foreach (P_Node node in Links3) p3 += node.P;
                    //    Pnf += p3 / 4 - Pn_1; ;
                    //}

                    //public virtual void UpdateCCP()
                    //{
                    //    double p2 = 0;
                    //    foreach (P_Node node in Links2) p2 += node.P;
                    //    Pnf += p2 / 4 - Pn - Pn_1;
                    //}

                    public override void UpdateP()
                    {
                        X = Xpos_Link.P + Xneg_Link.P;
                        Y = Ypos_Link.P + Yneg_Link.P;
                        Z = Zpos_Link.P + Zneg_Link.P;
                        ////IWB
                        //Pnf = 0.25 * (X + Y + Z) - 1.5 * Pn - Pn_1;
                        //double p2 = 0;
                        //foreach (P_Node node in Links2) p2 += node.P;
                        //Pnf += p2 * 0.125;
                        //double p3 = 0;
                        //foreach (P_Node node in Links3) p3 += node.P;
                        //Pnf += p3 * 0.0625;
                        ////IISO2
                        Pnf = (15.0 / 48.0) * (X + Y + Z) - (9.0 / 8.0) * Pn - Pn_1;
                        double p2 = 0;
                        foreach (P_Node node in Links2) p2 += node.P;
                        Pnf += p2 * (3.0 / 32.0);
                        double p3 = 0;
                        foreach (P_Node node in Links3) p3 += node.P;
                        Pnf += p3 * (1.0 / 64.0);
                        Pnf -= (9d / 8d) * Pn - Pn_1;
                        ////IISO
                        //Pnf = (.25) * (X + Y + Z) - Pn - Pn_1;
                        //double p2 = 0;
                        //foreach (P_Node node in Links2) p2 += node.P;
                        //Pnf += p2 * (.125);
                        ////CCP
                        //Pnf = -Pn - Pn_1;//-Pn - Pn_1;
                        //double p2 = 0;
                        //foreach (P_Node node in Links2) p2 += node.P;
                        //Pnf += p2 * (.25);
                        ////OCTA
                        //Pnf = -Pn_1;
                        //double p3 = 0;
                        //foreach (P_Node node in Links3) p3 += node.P;
                        //Pnf += p3 * (.25);
                        ////SLF
                        //Pnf = (1/3) * (X + Y + Z) - (9.0/8.0) * Pn - Pn_1;
                    }
                }

                public class RDD_Node : Node
                {
                    public Node[] Links2;
                    public RDD_Node(Point loc)//, double rho0, double dt, double dx, double c, int[] id_in)
                    : base(loc)
                    { }

                    public override void Link_Nodes(ref Node[][][] Frame, int x, int y, int z)
                    {
                        int mod = x % 2;

                        int xdim = Frame.Length - 1;
                        int ydim = Frame[mod].Length - 1;
                        int zdim = Frame[mod][0].Length - 1;

                        Links2 = new Node[12];

                        /*
                        [0] = x+y+z+
                        [1] = x+y-z-
                        [2] = x+y-z+
                        [3] = x+y+z-
                        [4] = x-y+z+
                        [5] = x-y-z-
                        [6] = x-y-z+
                        [7] = x-y+z-
                        [8] = y+
                        [9] = y-
                        [10] = z+
                        [11] = z-
                        */

                        if (mod == 0)
                        {
                            if (x < xdim)
                            {
                                if (y < ydim && z < zdim) Links2[0] = Frame[x + 1][y + 1][z + 1]; else Links2[0] = new Null_Node();
                                if (y > 0 && z > 0) Links2[1] = Frame[x + 1][y][z]; else Links2[1] = new Null_Node();
                                if (y > 0 && z < zdim) Links2[2] = Frame[x + 1][y][z + 1]; else Links2[2] = new Null_Node();
                                if (y < ydim && z > 0) Links2[3] = Frame[x + 1][y + 1][z]; else Links2[3] = new Null_Node();
                            }
                            else for (int i = 0; i < 4; i++) Links2[i] = new Null_Node();
                            if (x > 0)
                            {
                                if (y < ydim && z < zdim) Links2[4] = Frame[x - 1][y + 1][z + 1]; else Links2[4] = new Null_Node();
                                if (y > 0 && z > 0) Links2[5] = Frame[x - 1][y][z]; else Links2[5] = new Null_Node();
                                if (y > 0 && z < zdim) Links2[6] = Frame[x - 1][y][z + 1]; else Links2[6] = new Null_Node();
                                if (y < ydim && z > 0) Links2[7] = Frame[x - 1][y + 1][z]; else Links2[7] = new Null_Node();
                            }
                            else for (int i = 4; i < 8; i++) Links2[i] = new Null_Node();
                        }
                        else
                        {
                            if (x < xdim)
                            {
                                if (y < ydim && z < zdim) Links2[0] = Frame[x + 1][y][z]; else Links2[0] = new Null_Node();
                                if (y > 0 && z > 0) Links2[1] = Frame[x + 1][y - 1][z - 1]; else Links2[1] = new Null_Node();
                                if (y > 0 && z < zdim) Links2[2] = Frame[x + 1][y - 1][z]; else Links2[2] = new Null_Node();
                                if (y < ydim && z > 0) Links2[3] = Frame[x + 1][y][z - 1]; else Links2[3] = new Null_Node();
                            }
                            else for (int i = 0; i < 4; i++) Links2[i] = new Null_Node();
                            if (x > 0)
                            {
                                if (y < ydim && z < zdim) Links2[4] = Frame[x - 1][y][z]; else Links2[4] = new Null_Node();
                                if (y > 0 && z > 0) Links2[5] = Frame[x - 1][y - 1][z - 1]; else Links2[5] = new Null_Node();
                                if (y > 0 && z < zdim) Links2[6] = Frame[x - 1][y - 1][z]; else Links2[6] = new Null_Node();
                                if (y < ydim && z > 0) Links2[7] = Frame[x - 1][y][z - 1]; else Links2[7] = new Null_Node();
                            }
                            else for (int i = 4; i < 8; i++) Links2[i] = new Null_Node();
                        }

                        if (y < ydim) Links2[8] = Frame[x][y + 1][z]; else Links2[8] = new Null_Node();
                        if (y > 0) Links2[9] = Frame[x][y - 1][z]; else Links2[9] = new Null_Node();
                        if (z < zdim) Links2[10] = Frame[x][y][z + 1]; else Links2[10] = new Null_Node();
                        if (z > 0) Links2[11] = Frame[x][y][z - 1]; else Links2[11] = new Null_Node();

                        //foreach (Node n in Links2) Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(Utilities.PachTools.HPttoRPt(Pt), Utilities.PachTools.HPttoRPt(n.Pt));
                    }

                    public override void UpdateP()
                    {
                        double p2 = 0;
                        foreach (Node node in Links2) p2 += node.P;
                        Pnf = p2 * 0.25 - Pn - Pn_1;
                    }
                }

                public class Null_Node : Node
                {
                    public Null_Node()//Point loc, double rho0, double dt, double dx, double c, int[] id_in)
                        : base(new Point())//, rho0, dt, dx, c, id_in)
                    {
                        Instance = this;
                    }

                    public static Null_Node Instance
                    {
                        get;
                        private set;
                    }

                    public override void UpdateP()
                    {
                        Pnf = 0;
                        Pn = 0;
                        Pn_1 = 0;
                    }

                    public override void UpdateT()
                    {
                    }

                    public override void Link_Nodes(ref Node[][][] Frame, int x, int y, int z)
                    {
                    }
                }

                private class PML
                {
                    List<Node>[] Layers;
                    double[] Coefficients;
                    int layerCt;

                    public PML(int no_of_Layers, double max_Coefficient, Node[][][] PFrame, bool freefield)
                    {
                        Layers = new List<Node>[no_of_Layers];
                        Coefficients = new double[no_of_Layers];
                        layerCt = no_of_Layers;

                        for (int i = 0; i < layerCt; i++)
                        {
                            Coefficients[i] = 1 - max_Coefficient * (layerCt - i) / layerCt;
                            Layers[i] = new List<Node>();

                            for (int y = 0; y < PFrame[i].Length - i; y++) for (int z = 0; z < PFrame[i][y].Length - i; z++) Layers[i].Add(PFrame[i][y][z]);
                            int last = PFrame.Length - 1 - i;
                            for (int y = 0; y < PFrame[last].Length - i; y++) for (int z = 0; z < PFrame[last][y].Length - i; z++) Layers[i].Add(PFrame[last][y][z]);
                            int end = PFrame.Length - i - 1;
                            for (int x = i + 1; x < end; x++)
                            {
                                for (int j = i; j < PFrame[x][i].Length; j++) Layers[i].Add(PFrame[x][i][j]);
                                for (int j = i; j < PFrame[x][PFrame[x].Length - i - 1].Length; j++) Layers[i].Add(PFrame[x][PFrame[x].Length - i - 1][j]);
                                for (int j = i + 1; j < PFrame[x].Length - 1; j++)
                                {
                                    if (freefield) Layers[i].Add(PFrame[x][j][i]);
                                    Layers[i].Add(PFrame[x][j][PFrame[x][j].Length - 1 - (int)Math.Round(i/Utilities.Numerics.rt2)]);
                                }
                            }
                        }
                    }

                    public void Attenuate()
                    {
                        for (int i = 0; i < layerCt; i++)
                        {
                            foreach (Node n in Layers[i])
                            {
                                n.Attenuate(Coefficients[i]);
                            }
                        }
                    }
                }
            }

            public class Signal_Driver_Compact
            {
                List<Acoustic_Compact_FDTD.Node> SrcNode;
                public List<int> X, Y, Z;
                double[][] signal;
                double f;
                List<Hare.Geometry.Point> Loc;
                List<double[]> SWL;
                Signal_Type S;
                double w;
                double tmax;
                double dt;
                public List<int> delays;
                public int delayshortcut = 0;

                public enum Signal_Type
                {
                    Dirac_Pulse,
                    Sine_Tone,
                    Gaussian_Pulse,
                    Sine_Pulse,
                    SteadyState_Noise,
                    SS_Noise_Pulse,
                    Spectrum
                }

                public Signal_Driver_Compact(Signal_Type S_in, double freq, double w_in, Source[] Loc_in)
                {
                    Loc = new List<Point>();
                    SWL = new List<double[]>();
                    delays = new List<int>();
                    S = S_in;
                    f = freq;
                    w = w_in;
                    for (int i = 0; i < Loc_in.Length; i++)
                    {
                        if (Loc_in[i].Type() == "Line Source")
                        {
                            Loc.AddRange((Loc_in[i] as LineSource).Samples);
                            for (int j = 0; j < (Loc_in[i] as LineSource).Samples.Length; j++) SWL.Add((Loc_in[i] as LineSource).Power);
                        }
                        else
                        { 
                            Loc.Add(Loc_in[i].Origin);
                            SWL.Add(Loc_in[i].SoundPower);
                        }
                    }
                }

                public void Connect_Grid_Freefield(Acoustic_Compact_FDTD.Node[][][] Frame, AABB Bounds, double dx, double dy, double dz, double _tmax, double _dt, int no_of_Layers)
                {
                    tmax = _tmax / 1000;
                    dt = _dt;

                    SrcNode = new List<Acoustic_Compact_FDTD.Node>();
                    X = new List<int>();
                    Y = new List<int>();
                    Z = new List<int>();

                    List<Point> TempLOC = Loc;
                    List<double[]> TempSWL = SWL;
                    Loc = new List<Point>();
                    SWL = new List<double[]>();

                    for (int i = 0; i < TempLOC.Count; i++)
                    {
                        X.Add((int)Math.Floor((TempLOC[i].x - Bounds.Min_PT.x) / (dx)));// + (int)(no_of_Layers) / Utilities.Numerics.rt2);
                        if (X[X.Count-1] >= Frame.Length || X[X.Count-1] < 0) continue;
                        Y.Add((int)Math.Floor((TempLOC[i].y - Bounds.Min_PT.y) / (dx * Utilities.Numerics.rt2)));// + no_of_Layers / 2;
                        if (Y[Y.Count-1] >= Frame[X[i]].Length || Y[Y.Count-1] < 0) continue;
                        Z.Add((int)Math.Floor((TempLOC[i].z - Bounds.Min_PT.z) / (dx * Utilities.Numerics.rt2)));// + no_of_Layers / 2;
                        if (Z[Z.Count-1] >= Frame[X[i]][Y[i]].Length || Z[Z.Count-1] < 0) continue;
                        SrcNode.Add(Frame[X[X.Count-1]][Y[Y.Count-1]][Z[Z.Count-1]]);
                        Loc.Add(TempLOC[i]);
                        SWL.Add(TempSWL[i]);
                        delays.Add(0);
                    }

                    Generate_Signal();
                }

                public void Connect_Grid_Laboratory(Acoustic_Compact_FDTD.Node[][][] Frame, AABB Bounds, AABB Inner_Bounds, double dx, double dy, double dz, double _tmax, double _dt, int no_of_Layers)
                {
                    tmax = _tmax / 1000;
                    dt = _dt;

                    SrcNode = new List<Acoustic_Compact_FDTD.Node>();
                    X = new List<int>();
                    Y = new List<int>();
                    Z = new List<int>();

                    List<double[]> TempSWL = SWL;
                    SWL = new List<double[]>();

                    int[] minIndices = new int[3];
                    int[] maxIndices = new int[3];
                    minIndices[0] = (int)Math.Floor((Inner_Bounds.Min_PT.x - Bounds.Min_PT.x) / dx);
                    minIndices[1] = (int)Math.Floor((Inner_Bounds.Min_PT.y - Bounds.Min_PT.y) / (dx * Utilities.Numerics.rt2));
                    minIndices[2] = (int)Math.Floor((Inner_Bounds.Min_PT.z - Bounds.Min_PT.z) / (dx * Utilities.Numerics.rt2));
                    maxIndices[0] = (int)Math.Floor((Inner_Bounds.Max_PT.x - Bounds.Min_PT.x) / dx);
                    maxIndices[1] = (int)Math.Floor((Inner_Bounds.Max_PT.y - Bounds.Min_PT.y) / (dx * Utilities.Numerics.rt2));
                    maxIndices[2] = (int)Math.Floor((Inner_Bounds.Max_PT.z - Bounds.Min_PT.z) / (dx * Utilities.Numerics.rt2));

                    for (int i = 0; i < Loc.Count; i++)
                    {
                        if (Inner_Bounds.IsPointInBox(Loc[i].x, Loc[i].y, Loc[i].z))
                        {
                            X.Add((int)Math.Floor((Loc[i].x - Bounds.Min_PT.x) / (dx)));
                            if (X[X.Count - 1] >= Frame.Length || X[X.Count - 1] < 0) continue;
                            Y.Add((int)Math.Floor((Loc[i].y - Bounds.Min_PT.y) / (dx * Utilities.Numerics.rt2)));
                            if (Y[Y.Count - 1] >= Frame[X[i]].Length || Y[Y.Count - 1] < 0) continue;
                            Z.Add((int)Math.Floor((Loc[i].z - Bounds.Min_PT.z) / (dx * Utilities.Numerics.rt2)));
                            if (Z[Z.Count - 1] >= Frame[X[i]][Y[i]].Length || Z[Z.Count - 1] < 0) continue;
                            SrcNode.Add(Frame[X[X.Count - 1]][Y[Y.Count - 1]][Z[Z.Count - 1]]);
                            delays.Add(0);
                            SWL.Add(TempSWL[i]);
                        }
                        else
                        {
                            // Source is outside - find all visible faces and iterate only over valid indices on those faces

                            // Check X-Min face (left face)
                            if (Loc[i].x < Inner_Bounds.Min_PT.x)
                            {
                                int x = minIndices[0]; // Fixed x at minimum face

                                // Loop only through valid y and z on this face
                                for (int y = Math.Max(0, minIndices[1]); y <= Math.Min(maxIndices[1], Frame[x].Length - 1); y++)
                                {
                                    for (int z = Math.Max(0, minIndices[2]); z <= Math.Min(maxIndices[2], Frame[x][y].Length - 1); z++)
                                    {
                                        X.Add(x);
                                        Y.Add(y);
                                        Z.Add(z);
                                        SrcNode.Add(Frame[x][y][z]);
                                        double d = (Loc[i] - Acoustic_Compact_FDTD.RDD_Location(Bounds.Min_PT, x, y, z, dx, dx * Utilities.Numerics.rt2, dx * Utilities.Numerics.rt2)).Length();
                                        double poweradj = Math.Log10(4 * Math.PI * d * d);
                                        SWL.Add(new double[8]{TempSWL[i][0] + poweradj, TempSWL[i][1] + poweradj, TempSWL[i][2] + poweradj, TempSWL[i][3] + poweradj, TempSWL[i][4] + poweradj, TempSWL[i][5] + poweradj, TempSWL[i][6] + poweradj, TempSWL[i][7] + poweradj });
                                        delays.Add((int)Math.Round((d / 343)/dt)); // Store delay for this source
                                    }
                                }
                            }

                            // Check X-Max face (right face)
                            if (Loc[i].x > Inner_Bounds.Max_PT.x)
                            {
                                int x = maxIndices[0]; // Fixed x at maximum face

                                // Loop only through valid y and z on this face
                                for (int y = Math.Max(0, minIndices[1]); y <= Math.Min(maxIndices[1], Frame[x].Length - 1); y++)
                                {
                                    for (int z = Math.Max(0, minIndices[2]); z <= Math.Min(maxIndices[2], Frame[x][y].Length - 1); z++)
                                    {
                                        X.Add(x);
                                        Y.Add(y);
                                        Z.Add(z);
                                        SrcNode.Add(Frame[x][y][z]);
                                        double d = (Loc[i] - Acoustic_Compact_FDTD.RDD_Location(Bounds.Min_PT, x, y, z, dx, dx * Utilities.Numerics.rt2, dx * Utilities.Numerics.rt2)).Length();
                                        double poweradj = Math.Log10(4 * Math.PI * d * d);
                                        SWL.Add(new double[8] { TempSWL[i][0] + poweradj, TempSWL[i][1] + poweradj, TempSWL[i][2] + poweradj, TempSWL[i][3] + poweradj, TempSWL[i][4] + poweradj, TempSWL[i][5] + poweradj, TempSWL[i][6] + poweradj, TempSWL[i][7] + poweradj });
                                        delays.Add((int)Math.Round((d / 343) / dt)); // Store delay for this source
                                    }
                                }
                            }

                            // Check Y-Min face (front face)
                            if (Loc[i].y < Inner_Bounds.Min_PT.y)
                            {
                                int y = minIndices[1]; // Fixed y at minimum face

                                // Loop through x-range, but skip corners already handled by x-faces
                                for (int x = Math.Max(0, minIndices[0] + 1); x <= Math.Min(maxIndices[0] - 1, Frame.Length - 1); x++)
                                {
                                    if (y >= Frame[x].Length) continue; // Skip invalid indices

                                    for (int z = Math.Max(0, minIndices[2]); z <= Math.Min(maxIndices[2], Frame[x][y].Length - 1); z++)
                                    {
                                        X.Add(x);
                                        Y.Add(y);
                                        Z.Add(z);
                                        SrcNode.Add(Frame[x][y][z]);
                                        double d = (Loc[i] - Acoustic_Compact_FDTD.RDD_Location(Bounds.Min_PT, x, y, z, dx, dx * Utilities.Numerics.rt2, dx * Utilities.Numerics.rt2)).Length();
                                        double poweradj = Math.Log10(4 * Math.PI * d * d);
                                        SWL.Add(new double[8] { TempSWL[i][0] + poweradj, TempSWL[i][1] + poweradj, TempSWL[i][2] + poweradj, TempSWL[i][3] + poweradj, TempSWL[i][4] + poweradj, TempSWL[i][5] + poweradj, TempSWL[i][6] + poweradj, TempSWL[i][7] + poweradj });
                                        delays.Add((int)Math.Round((d / 343) / dt)); // Store delay for this source
                                    }
                                }
                            }

                            // Check Y-Max face (back face)
                            if (Loc[i].y > Inner_Bounds.Max_PT.y)
                            {
                                int y = maxIndices[1]; // Fixed y at maximum face

                                // Loop through x-range, but skip corners already handled by x-faces
                                for (int x = Math.Max(0, minIndices[0] + 1); x <= Math.Min(maxIndices[0] - 1, Frame.Length - 1); x++)
                                {
                                    if (y >= Frame[x].Length) continue; // Skip invalid indices

                                    for (int z = Math.Max(0, minIndices[2]); z <= Math.Min(maxIndices[2], Frame[x][y].Length - 1); z++)
                                    {
                                        X.Add(x);
                                        Y.Add(y);
                                        Z.Add(z);
                                        SrcNode.Add(Frame[x][y][z]);
                                        double d = (Loc[i] - Acoustic_Compact_FDTD.RDD_Location(Bounds.Min_PT, x, y, z, dx, dx * Utilities.Numerics.rt2, dx * Utilities.Numerics.rt2)).Length();
                                        double poweradj = Math.Log10(4 * Math.PI * d * d);
                                        SWL.Add(new double[8] { TempSWL[i][0] + poweradj, TempSWL[i][1] + poweradj, TempSWL[i][2] + poweradj, TempSWL[i][3] + poweradj, TempSWL[i][4] + poweradj, TempSWL[i][5] + poweradj, TempSWL[i][6] + poweradj, TempSWL[i][7] + poweradj });
                                        delays.Add((int)Math.Round((d / 343) / dt)); // Store delay for this source
                                    }
                                }
                            }

                            // Check Z-Min face (bottom face)
                            if (Loc[i].z < Inner_Bounds.Min_PT.z)
                            {
                                int z = minIndices[2]; // Fixed z at minimum face

                                // Loop through x and y ranges, but skip edges already handled by x and y faces
                                for (int x = Math.Max(0, minIndices[0] + 1); x <= Math.Min(maxIndices[0] - 1, Frame.Length - 1); x++)
                                {
                                    for (int y = Math.Max(0, minIndices[1] + 1); y <= Math.Min(maxIndices[1] - 1, Frame[x].Length - 1); y++)
                                    {
                                        if (z >= Frame[x][y].Length) continue; // Skip invalid indices

                                        X.Add(x);
                                        Y.Add(y);
                                        Z.Add(z);
                                        SrcNode.Add(Frame[x][y][z]);
                                        double d = (Loc[i] - Acoustic_Compact_FDTD.RDD_Location(Bounds.Min_PT, x, y, z, dx, dx * Utilities.Numerics.rt2, dx * Utilities.Numerics.rt2)).Length();
                                        double poweradj = Math.Log10(4 * Math.PI * d * d);
                                        SWL.Add(new double[8] { TempSWL[i][0] + poweradj, TempSWL[i][1] + poweradj, TempSWL[i][2] + poweradj, TempSWL[i][3] + poweradj, TempSWL[i][4] + poweradj, TempSWL[i][5] + poweradj, TempSWL[i][6] + poweradj, TempSWL[i][7] + poweradj });
                                        delays.Add((int)Math.Round((d / 343) / dt)); // Store delay for this source
                                    }
                                }
                            }

                            // Check Z-Max face (top face)
                            if (Loc[i].z > Inner_Bounds.Max_PT.z)
                            {
                                int z = maxIndices[2]; // Fixed z at maximum face

                                // Loop through x and y ranges, but skip edges already handled by x and y faces
                                for (int x = Math.Max(0, minIndices[0] + 1); x <= Math.Min(maxIndices[0] - 1, Frame.Length - 1); x++)
                                {
                                    for (int y = Math.Max(0, minIndices[1] + 1); y <= Math.Min(maxIndices[1] - 1, Frame[x].Length - 1); y++)
                                    {
                                        if (z >= Frame[x][y].Length) continue; // Skip invalid indices

                                        X.Add(x);
                                        Y.Add(y);
                                        Z.Add(z);
                                        SrcNode.Add(Frame[x][y][z]);
                                        double d = (Loc[i] - Acoustic_Compact_FDTD.RDD_Location(Bounds.Min_PT, x, y, z, dx, dx * Utilities.Numerics.rt2, dx * Utilities.Numerics.rt2)).Length();
                                        double poweradj = Math.Log10(4 * Math.PI * d * d);
                                        SWL.Add(new double[8] { TempSWL[i][0] + poweradj, TempSWL[i][1] + poweradj, TempSWL[i][2] + poweradj, TempSWL[i][3] + poweradj, TempSWL[i][4] + poweradj, TempSWL[i][5] + poweradj, TempSWL[i][6] + poweradj, TempSWL[i][7] + poweradj });
                                        delays.Add((int)Math.Round((d / 343) / dt)); // Store delay for this source
                                    }
                                }
                            }
                        }
                    }

                    delayshortcut = delays.Min();

                    Generate_Signal();
                }

                public void reset(double frequency, Signal_Type s)
                {
                    f = frequency;
                    S = s;
                    Generate_Signal();
                }

                private void Generate_Signal()
                {
                    signal = new double[Loc.Count][];
                    Random R = new Random((int)System.DateTime.Now.Ticks);
                    double[] noise = new double[(int)Math.Ceiling(tmax / dt)];
                    for (int n = 0; n < tmax / dt; n++) noise[n] = R.NextDouble();
                    
                    for(int i = 0; i < SWL.Count; i++)
                    {
                        double p = 0;
                        for (int j = 0; j < SWL[i].Length; j++) SWL[i][j] = Pachyderm_Acoustic.Utilities.AcousticalMath.Pressure_SPL(SWL[i][j]);
                    }

                    double f2pi = f * 2 * Math.PI;

                    for (int i = 0; i < Loc.Count; i++)
                    {
                        R = new Random((int)System.DateTime.Now.Ticks);
                        signal[i] = new double[(int)Math.Ceiling(tmax / dt)];

                        switch (S)
                        {
                            case Signal_Type.Dirac_Pulse:
                                signal[i][1] = 1;
                                break;
                            case Signal_Type.Sine_Tone:
                                double sum_e = 0;
                                for (int n = 0; n < tmax / dt; n++)
                                {
                                    double ph = f2pi * n * dt;
                                    signal[i][n] = Math.Sin(ph);
                                    sum_e += signal[i][n] * signal[i][n];
                                }
                                for (int j = 0; j < signal[i].Length; j++) signal[i][j] *= tmax * 20 / Math.Sqrt(sum_e);
                                break;
                            case Signal_Type.Gaussian_Pulse:
                                double sum = 0;
                                for (int n = 0; n < tmax / dt; n++)
                                {
                                    signal[i][n] = Math.Exp(-.5 * Math.Pow((double)n / w, 2));
                                    sum += signal[i][n];
                                }
                                for (int n = 0; n < tmax / dt; n++)
                                {
                                    signal[i][n] /= sum;
                                }
                                break;
                            case Signal_Type.Sine_Pulse:
                                //for (int n = 0; n < tmax / dt; n++) signal[n] = Math.Exp(-.5 * Math.Pow((double)n / 1, 2)) * Math.Sin(f2pi * n * dt);
                                double kk = Math.PI / 60 / 343 / dt;
                                double offset = Math.PI / kk / 343;
                                signal[i] = new double[60];
                                double param = 2 * (0.371 * 60 - 8.1);
                                double sumsig = 0;
                                for (int n = 0; n < 60; n++)
                                {
                                    double t = n * dt;
                                    signal[i][n] = (t - offset / 2) * Math.Pow(Math.Sin(kk * 343 * t), param);
                                    sumsig += signal[i][n] * signal[i][n];
                                }
                                for (int j = 0; j < signal[i].Length; j++) signal[i][j] = Math.Sqrt(signal[i][j] * signal[i][j] / sumsig) * ((signal[i][j] < 0) ? -1 : 1);
                                break;
                            case Signal_Type.SteadyState_Noise:
                                for (int n = 0; n < tmax / dt; n++) signal[i][n] = R.NextDouble();
                                break;
                            case Signal_Type.SS_Noise_Pulse:
                                for (int n = 0; n < tmax / dt; n++) signal[i][n] = Math.Exp(-.5 * Math.Pow((double)n / w, 2) * R.NextDouble());
                                break;
                            case Signal_Type.Spectrum:
                                signal[i] = Pachyderm_Acoustic.Audio.Pach_SP.Filter2Signal(noise, SWL[i], (int)(1.0 / dt), 0);
                                double[][] filt = new double[8][];
                                for(int o = 0; o < 8; o++)
                                {
                                    filt[o] = Pachyderm_Acoustic.Audio.Pach_SP.FIR_Bandpass(signal[i], o, (int)(1.0 / dt), 0);
                                    double p = 0;
                                    for (int j = 0; j < filt[o].Length; j++) p += Math.Abs(filt[o][j]);

                                    p =Pachyderm_Acoustic.Utilities.AcousticalMath.SPL_Pressure(p);
                                }

                                break;
                        }
                    }
                }

                public void Drive(int t)
                {
                    t += delayshortcut;
                    for (int i = 0; i < SrcNode.Count; i++)
                    {
                        if (t < signal[0].Length + delays[i])
                        {
                            if (t < delays[i])
                            {
                                SrcNode[i].Pn = 0;
                            }
                            else if (t < signal[0].Length + delays[i])
                            {
                                SrcNode[i].Pn = signal[0][t - delays[i]] * SWL[i][4];
                            }
                            else
                            {
                                SrcNode[i].Pn = 0;
                            }
                        }
                        else SrcNode[i].Pn = 0;
                    }
                }

                public System.Numerics.Complex[] Frequency_Response(int length)
                {
                    if (signal.Length > length) length = signal.Length;
                    double[] signal_clone = new double[length];
                    Array.Copy(signal, signal_clone, signal.Length);
                    return Audio.Pach_SP.FFT_General(signal_clone, 0);
                }

                public double frequency
                {
                    get
                    {
                        return f;
                    }
                    set
                    {
                        f = value;
                        Generate_Signal();
                    }
                }
            }

            public class Microphone_Compact
            {
                double[][] recording;
                List<double[][]> storage = new List<double[][]>();
                Acoustic_Compact_FDTD.Node[] RecNode;
                Hare.Geometry.Point[] Loc;
                public int[] X, Y, Z;
                double tmax;
                double Sample_Freq;
                int no_of_samples;
                public Topology grid_template;

                public Microphone_Compact()
                {
                    Loc = new Point[0];
                    Random R = new Random();
                    RecNode = new Acoustic_Compact_FDTD.Node[Loc.Length];
                }

                public Microphone_Compact(Hare.Geometry.Point[] Loc_in)
                {
                    Loc = Loc_in;
                    Random R = new Random();
                    RecNode = new Acoustic_Compact_FDTD.Node[Loc.Length];
                }

                public void Connect_Grid_Freefield(Acoustic_Compact_FDTD.Node[][][] Frame, AABB Bounds, double dx, double _tmax, double dt, int no_of_Layers)
                {
                    Sample_Freq = 1 / dt;
                    no_of_samples = (int)Math.Ceiling(_tmax / dt / 1000);
                    tmax = _tmax;
                    recording = new double[Loc.Length][];
                    X = new int[Loc.Length];
                    Y = new int[Loc.Length];
                    Z = new int[Loc.Length];

                    for (int i = 0; i < Loc.Length; i++)
                    {
                        recording[i] = new double[no_of_samples];
                        X[i] = (int)Math.Floor((Loc[i].x - Bounds.Min_PT.x) / (dx));// + (int)(no_of_Layers) / Utilities.Numerics.rt2);
                        Y[i] = (int)Math.Floor((Loc[i].y - Bounds.Min_PT.y) / (dx * Utilities.Numerics.rt2));// + no_of_Layers / 2;
                        Z[i] = (int)Math.Floor((Loc[i].z - Bounds.Min_PT.z) / (dx * Utilities.Numerics.rt2));// + no_of_Layers / 2;
                        RecNode[i] = Frame[X[i]][Y[i]][Z[i]];
                    }
                }

                public void Connect_Grid_Laboratory(Acoustic_Compact_FDTD.Node[][][] Frame, AABB Bounds, double dx, double _tmax, double dt, int no_of_Layers)
                {
                    no_of_samples = (int)Math.Ceiling(_tmax / dt / 1000);
                    tmax = _tmax;
                    recording = new double[Loc.Length][];
                    X = new int[Loc.Length];
                    Y = new int[Loc.Length];
                    Z = new int[Loc.Length];

                    for (int i = 0; i < Loc.Length; i++)
                    {
                        recording[i] = new double[no_of_samples];
                        X[i] = (int)Math.Floor((Loc[i].x - Bounds.Min_PT.x) / (dx));// + (int)(no_of_Layers) / Utilities.Numerics.rt2);
                        Y[i] = (int)Math.Floor((Loc[i].y - Bounds.Min_PT.y) / (dx * Utilities.Numerics.rt2));// + no_of_Layers / 2;
                        Z[i] = (int)Math.Floor((Loc[i].z - Bounds.Min_PT.z) / (dx * Utilities.Numerics.rt2));// + no_of_Layers / 2;
                        RecNode[i] = Frame[X[i]][Y[i]][Z[i]];
                    }
                }

                public void Connect_Grid_UniqueOnly_Freefield(Acoustic_Compact_FDTD.Node[][][] Frame, AABB Bounds, double dx, double _tmax, double dt, int no_of_Layers)
                {
                    no_of_samples = (int)Math.Ceiling(_tmax / dt / 1000);
                    tmax = _tmax;

                    bool[,,] Check_Frame = new bool[Frame.Length, Frame[0].Length, Frame[0][0].Length];
                    List<Acoustic_Compact_FDTD.Node> Rec_Unique = new List<Acoustic_Compact_FDTD.Node>();
                    X = new int[Loc.Length];
                    Y = new int[Loc.Length];
                    Z = new int[Loc.Length];

                    int j = 0;
                    for (int i = 0; i < Loc.Length; i++)
                    {
                        int x = (int)Math.Floor((Loc[i].x - Bounds.Min_PT.x) / (dx));// + (int)(no_of_Layers) / Utilities.Numerics.rt2);
                        int y = (int)Math.Floor((Loc[i].y - Bounds.Min_PT.y) / (dx * Utilities.Numerics.rt2));// + no_of_Layers / 2;
                        int z = (int)Math.Floor((Loc[i].z - Bounds.Min_PT.z) / (dx * Utilities.Numerics.rt2));// + no_of_Layers / 2;
                        //Rhino.RhinoDoc.ActiveDoc.Objects.AddPoint(Utilities.PachTools.HPttoRPt(Acoustic_Compact_FDTD.RDD_Location(Bounds.Min_PT, X, Y, Z, dx, dx * Utilities.Numerics.rt2, dx * Utilities.Numerics.rt2)));
                        if (!Check_Frame[x, y, z])
                        {
                            X[j] = x; Y[j] = y; Z[j] = z;
                            j++;
                            Check_Frame[X[i], Y[i], Z[i]] = true;
                            Rec_Unique.Add(Frame[X[i]][Y[i]][Z[i]]);
                        }
                    }
                    Array.Resize(ref X, j);
                    Array.Resize(ref Y, j);
                    Array.Resize(ref Z, j);

                    recording = new double[Rec_Unique.Count][];
                    for (int i = 0; i < recording.Length; i++) recording[i] = new double[no_of_samples];

                    RecNode = Rec_Unique.ToArray();
                }

                //public void Connect_Grid_UniqueOnly_Laboratory(Acoustic_Compact_FDTD.Node[][][] Frame, AABB Bounds, double dx, double _tmax, double dt, int no_of_Layers)
                //{
                //    no_of_samples = (int)Math.Ceiling(_tmax / dt / 1000);
                //    tmax = _tmax;

                //    bool[,,] Check_Frame = new bool[Frame.Length, Frame[0].Length, Frame[0][0].Length];
                //    List<Acoustic_Compact_FDTD.Node> Rec_Unique = new List<Acoustic_Compact_FDTD.Node>();
                //    X = new int[Loc.Length];
                //    Y = new int[Loc.Length];
                //    Z = new int[Loc.Length];

                //    for (int i = 0; i < Loc.Length; i++)
                //    {
                //        X[i] = (int)Math.Floor((Loc[i].x - Bounds.Min_PT.x) / (dx));// + (int)(no_of_Layers) / Utilities.Numerics.rt2);
                //        Y[i] = (int)Math.Floor((Loc[i].y - Bounds.Min_PT.y) / (2 * dx / Utilities.Numerics.rt2));// + no_of_Layers / 2;
                //        Z[i] = (int)Math.Floor((Loc[i].z - Bounds.Min_PT.z) / (2 * dx / Utilities.Numerics.rt2));// + no_of_Layers / 2;
                //        //Rhino.RhinoDoc.ActiveDoc.Objects.AddPoint(Utilities.PachTools.HPttoRPt(Acoustic_Compact_FDTD.RDD_Location(Bounds.Min_PT, X, Y, Z, dx, dx * Utilities.Numerics.rt2, dx * Utilities.Numerics.rt2)));
                //        if (!Check_Frame[X[i], Y[i], Z[i]])
                //        {
                //            Check_Frame[X[i], Y[i], Z[i]] = true;
                //            Rec_Unique.Add(Frame[X[i]][Y[i]][Z[i]]);
                //        }
                //    }

                //    recording = new double[Rec_Unique.Count][];
                //    for (int i = 0; i < recording.Length; i++) recording[i] = new double[no_of_samples];

                //    RecNode = Rec_Unique.ToArray();
                //}
                public void Connect_Grid_Hemisphere_Laboratory(Acoustic_Compact_FDTD.Node[][][] Frame, AABB Bounds, Point center, double radius, double dx, double _tmax, double dt, int no_of_Layers)
                {
                    List<Point> Loc_List = new List<Point>();
                    no_of_samples = (int)Math.Ceiling(_tmax / dt / 1000);
                    tmax = _tmax;

                    List<Point> ptlist = new List<Point>();
                    List<Acoustic_Compact_FDTD.Node> Rec = new List<Acoustic_Compact_FDTD.Node>();
                    List<int> Xl = new List<int>();
                    List<int> Yl = new List<int>();
                    List<int> Zl = new List<int>();

                    // Revise this appraoch given the known coordinates on the surface of a sphere... this appraoch is ineffective...
                    //for (int i = 0; i < Frame.Length; i++)
                    //{
                    //    for (int j = 0; j < Frame[i].Length; j++)
                    //    {
                    //        for (int k = 0; k < Frame[i][j].Length; k++)
                    //        {
                    //            Point p = Acoustic_Compact_FDTD.RDD_Location(Bounds.Min_PT, i, j, k, dx, 4 * dx / Utilities.Numerics.rt2, 4 * dx / Utilities.Numerics.rt2);// /2
                    //            if (p.z < center.z) continue;
                    //            double dist = (p - center).Length();
                    //            if (dist < radius + dx / 2 && dist > radius - dx / 2)
                    //            {
                    //                Xl.Add(i);
                    //                Yl.Add(j);
                    //                Zl.Add(k);
                    //                ptlist.Add(p);
                    //                Rec.Add(Frame[i][j][k]);
                    //            }
                    //        }
                    //    }
                    //}

                    grid_template = Utilities.Geometry.GeoHemiSphere(5, radius);

                    for (int i = 0; i < grid_template.Vertex_Count; i++)
                    {
                        Point p = grid_template[i];            
                        int x = (int)Math.Floor((p.x - Bounds.Min_PT.x) / (dx));
                        int y = (int)Math.Floor((p.y - Bounds.Min_PT.y) / (dx * Utilities.Numerics.rt2));
                        int z = (int)Math.Max(Math.Floor((p.z + center.z - Bounds.Min_PT.z) / (dx * Utilities.Numerics.rt2)),0);

                        Xl.Add(x);
                        Yl.Add(y);
                        Zl.Add(z);
                        ptlist.Add(p);
                        Rec.Add(Frame[x][y][z]);
                    }

                    //for (int phi = 0; phi < 72; phi++)
                    //{
                    //    for (int theta = 0; theta < 20; theta++)
                    //    {
                    //        Point p = radius * new Point(Math.Cos(phi * Math.PI / 36) * Math.Cos(theta * Math.PI / 36), Math.Sin(phi * Math.PI / 36) * Math.Cos(theta * Math.PI / 36), Math.Sin(theta * Math.PI / 36));
                    //        int x = (int)Math.Floor((p.x - Bounds.Min_PT.x) / (dx));
                    //        int y = (int)Math.Floor((p.y - Bounds.Min_PT.y) / (dx * Utilities.Numerics.rt2));
                    //        int z = (int)Math.Floor((p.z - Bounds.Min_PT.z) / (dx * Utilities.Numerics.rt2));

                    //        Xl.Add(x);
                    //        Yl.Add(y);
                    //        Zl.Add(z);
                    //        ptlist.Add(p);
                    //        Rec.Add(Frame[x][y][z]);
                    //    }
                    //}

                    Loc = ptlist.ToArray();
                    X = Xl.ToArray();
                    Y = Yl.ToArray();
                    Z = Zl.ToArray();

                    recording = new double[Rec.Count][];
                    for (int i = 0; i < recording.Length; i++) recording[i] = new double[no_of_samples];

                    RecNode = Rec.ToArray();
                }

                public void Connect_Grid_UniqueOnly_Laboratory(Acoustic_Compact_FDTD.Node[][][] Frame, AABB Bounds, double dx, double _tmax, double dt, int no_of_Layers)
                {
                    no_of_samples = (int)Math.Ceiling(_tmax / dt / 1000);
                    tmax = _tmax;

                    bool[,,] Check_Frame = new bool[Frame.Length, Frame[0].Length, Frame[0][0].Length];
                    List<Acoustic_Compact_FDTD.Node> Rec_Unique = new List<Acoustic_Compact_FDTD.Node>();
                    X = new int[Loc.Length];
                    Y = new int[Loc.Length];
                    Z = new int[Loc.Length];

                    int j = 0;
                    for (int i = 0; i < Loc.Length; i++)
                    {
                        int x = (int)Math.Floor((Loc[i].x - Bounds.Min_PT.x) / (dx));
                        int y = (int)Math.Floor((Loc[i].y - Bounds.Min_PT.y) / (dx * Utilities.Numerics.rt2));
                        int z = (int)Math.Floor((Loc[i].z - Bounds.Min_PT.z) / (dx * Utilities.Numerics.rt2));

                        Point p = Acoustic_Compact_FDTD.RDD_Location(Bounds.Min_PT, x, y, z, dx, 2 * dx / Utilities.Numerics.rt2, 2 * dx / Utilities.Numerics.rt2);// /2

                        if (!Check_Frame[x, y, z])
                        {
                            X[j] = x; Y[j] = y; Z[j] = z;
                            j++;
                            Check_Frame[X[i], Y[i], Z[i]] = true;
                            Rec_Unique.Add(Frame[X[i]][Y[i]][Z[i]]);
                        }
                    }
                    Array.Resize(ref X, j);
                    Array.Resize(ref Y, j);
                    Array.Resize(ref Z, j);

                    recording = new double[Rec_Unique.Count][];
                    for (int i = 0; i < recording.Length; i++) recording[i] = new double[no_of_samples];

                    RecNode = Rec_Unique.ToArray();
                }

                public void Record(int n)
                {
                    for (int i = 0; i < RecNode.Length; i++)
                    {
                        recording[i][n] = RecNode[i].P;
                    }
                }

                public void reset()
                {
                    storage.Add(recording);
                    recording = new double[Loc.Length][];
                    for (int i = 0; i < RecNode.Length; i++)
                    {
                        recording[i] = new double[no_of_samples];
                    }
                }

                public double RMS(int index, int omitsamples = 0)
                {
                    double sum = 0;
                    if (index > storage[0].Length - 1) throw new Exception("Index out of range");
                    for (int i = omitsamples; i < storage[0][index].Length; i++) sum += storage[0][index][i] * storage[0][index][i];
                    return Math.Sqrt(sum);
                }

                public double RMS(int index, int octave_id, int omitsamples = 0)
                {
                    if (octave_id > 7) return RMS(index);
                    double sum = 0;
                    if (index > storage[0].Length - 1) throw new Exception("Index out of range");

                    double[] filtered = Pachyderm_Acoustic.Audio.Pach_SP.FIR_Bandpass(Recordings()[0][index], octave_id, (int)Math.Round(this.Sample_Freq), 0);
                    for (int i = omitsamples; i < filtered.Length; i++) sum += filtered[i] * filtered[i];

                    return Math.Sqrt(sum);
                }

                public List<double[][]> Recordings()
                {
                    return storage;
                }

                public double[][] Recordings_Current()
                {
                    return recording;
                }

                public double[] Recordings(int mic_id, int omit = 0, int storage_id = 0)
                {
                    if (omit == 0) return storage[0][mic_id];
                    int l = storage[storage_id][mic_id].Count() - omit;
                    double[] range = new double[l];
                    range = storage[storage_id][mic_id].ToList().GetRange(omit, l).ToArray();
                    return range;
                }
            }
        }
    }
}