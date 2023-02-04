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

            int mod = Vertices.Length - 1;

            for (int i = 0; i < Vertices.Length; i++)
            {
                if ((i+1) % 20 == 0) continue;
                Faces.Add(new int[3] {i % mod, (i + 1) % mod, (i + 20) % mod});
                Faces.Add(new int[3] { (i + 1) % mod, (i + 21) % mod, (i + 20) % mod });

                //for (int j = i + 1; j < Vertices.Length; j++)
                //{
                //    Hare.Geometry.Vector vd = Vertices[i] - Vertices[j];
                //    //if (i == j) continue;
                //    if (vd.x * vd.x + vd.y * vd.y + vd.z * vd.z > md2)
                //        continue;
                //    for (int k = j + 1; k < Vertices.Length; k++)
                //    {
                //        vd = Vertices[i] - Vertices[k];
                //        if (vd.x * vd.x + vd.y * vd.y + vd.z * vd.z < md2)
                //        {
                //            Hare.Geometry.Point ctr;
                //            double r2;
                //            Utilities.PachTools.CircumCircleRadius(Vertices[i], Vertices[j], Vertices[k], out r2, out ctr);
                //            if (r2 < max_dist * 10)
                //            {
                //                r2 *= r2;
                //                int[] f = new int[3] { i, j, k };
                //                //Parallel.For(0, Vertices.Length, i =>
                //                for (int x = 0; x < Vertices.Length; x++)
                //                {
                //                    if (x != f[0] && x != f[1] && x != f[2])
                //                    {
                //                        if (Utilities.PachTools.check_circumsphere(ctr, Vertices[x], r2)) Faces.Add(f);
                //                    }
                //                }//) ;
                //            }
                //        }
                //    }
                //}
            }

            for (int i = 0; i < pts.Count; i++)
            {
                Vertices[i].Normalize();
            }
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
                points[i] = mag * Vertices[i] + Ctr;
            }
            Hare.Geometry.Topology T = new Hare.Geometry.Topology();
            T.Set_Topology(points, Faces.ToArray());
            T.Finish_Topology();
            return T;
        }
    }
}