using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Pachyderm_Acoustic
{
    public class Sphere_Plot
    {
        Hare.Geometry.Vector[] Vertices;
        List<int[]> Faces;
        Hare.Geometry.Point Ctr;

        public Sphere_Plot(List<Hare.Geometry.Point> pts, Hare.Geometry.Point Center, double max_dist)
        {
            Ctr = Center;
            Faces = new List<int[]>();
            //Put data set in spherical coordinates...
            Vertices = new Hare.Geometry.Vector[pts.Count];
            for (int i = 0; i < pts.Count; i++)
            {
                Vertices[i] = pts[i] - Center;
                Vertices[i].Normalize();
            }

            double md2 = max_dist * max_dist * 25;

            for (int i = 0; i < Vertices.Length; i++)
            {
                for (int j = i+1; j < Vertices.Length; j++)
                {
                    Hare.Geometry.Vector vd = Vertices[i] - Vertices[j];
                    //if (i == j) continue;
                    if (vd.x * vd.x + vd.y * vd.y + vd.z * vd.z > md2)
                        continue;
                    for (int k = j+1; k < Vertices.Length; k++)
                    {
                        //if (i != k && j != k)
                        //{
                            vd = Vertices[i] - Vertices[k];
                            if (vd.x * vd.x + vd.y * vd.y + vd.z * vd.z < md2)
                            {
                                Hare.Geometry.Point ctr;
                                double r2;
                                CircumCircleRadius(Vertices[i], Vertices[j], Vertices[k], out r2, out ctr);
                                if (r2 < max_dist * 10)
                                {
                                    r2 *= r2;
                                    int[] f = new int[3] { i, j, k };
                                    if (check_circumsphere(ctr, r2, f)) Faces.Add(f);
                                }
                            }
                        //}
                    }
                }
            }

            for (int i = 0; i < pts.Count; i++)
            {
                Vertices[i].Normalize();
            }
        }

        public bool check_circumsphere(Hare.Geometry.Point ctr, double r2, int[] exclusive)
        {
             //Parallel.For(0, Vertices.Length, i =>
             for(int i = 0; i < Vertices.Length; i++)
             {
                 if (i != exclusive[0] && i != exclusive[1] && i != exclusive[2])
                 {
                     Hare.Geometry.Vector vd = Vertices[i] - ctr;
                     if (vd.x * vd.x + vd.y * vd.y + vd.z * vd.z < r2)
                        return false;
                 }
             }//) ;
            return true;
        }

        public void CircumCircleRadius(Hare.Geometry.Point a, Hare.Geometry.Point b, Hare.Geometry.Point c, out double radius, out Hare.Geometry.Point ctr)
        {
            Hare.Geometry.Point acmid = (a + c) / 2;
            Hare.Geometry.Point abmid = (a + b) / 2;
            Hare.Geometry.Vector ac = c - a;
            Hare.Geometry.Vector ab = a - b;
            Hare.Geometry.Vector abXac = Hare.Geometry.Hare_math.Cross(ab, ac);
            Hare.Geometry.Vector Pab = Hare.Geometry.Hare_math.Cross(abXac, ab);
            Hare.Geometry.Vector Pac = Hare.Geometry.Hare_math.Cross(abXac, ac);

            // this is the vector from a TO the circumsphere center
            radius = (Pab.x * acmid.y - Pab.x * abmid.y + Pab.y * abmid.x - Pab.y * acmid.x) / (Pac.x * Pab.y - Pab.x * Pac.y);
            ctr = acmid + Pac * radius;
            radius = (ctr - a).Length();
        }

        public Hare.Geometry.Topology Output(IEnumerable<double> magnitude, double Min, double Max, double Diameter)
        {
            if (magnitude.Count() != Vertices.Length) throw new Exception("Invalid data input to spherical plot...");
            Hare.Geometry.Point[] points = new Hare.Geometry.Point[Vertices.Length];
            for(int i = 0; i < magnitude.Count(); i++)
            {
                double mag = (magnitude.ElementAt(i));
                if (double.IsInfinity(mag)) mag = 0;
                mag = Math.Max(mag, Min);
                mag = Math.Min(mag, Max);
                mag -= Min;
                mag /= (Max - Min) * Diameter;
                points[i] = mag * Vertices[i];
            }
            Hare.Geometry.Topology T = new Hare.Geometry.Topology();
            T.Set_Topology(points, Faces.ToArray());
            T.Finish_Topology();
            return T;
        }
    }
}