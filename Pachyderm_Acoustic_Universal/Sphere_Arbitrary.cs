using Hare.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Pachyderm_Acoustic
{
    public class Sphere_Plot
    {
        Hare.Geometry.Point Ctr;
        Topology hemisphere;

        public Sphere_Plot(Hare.Geometry.Point Center)
        {
            hemisphere = Utilities.Geometry.GeoHemiSphere(5, 1);
            Ctr = Center;
        }

        public Hare.Geometry.Topology Output(IEnumerable<double> magnitude, double Min, double Max, double Diameter)
        {
            if (magnitude.Count() != hemisphere.Vertex_Count) throw new Exception("Invalid data input to spherical plot...");
            Hare.Geometry.Point[] points = new Hare.Geometry.Point[hemisphere.Vertex_Count];
            for(int i = 0; i < magnitude.Count(); i++)
            {
                double mag = (magnitude.ElementAt(i));
                if (double.IsInfinity(mag)) mag = 0;
                mag = Math.Max(mag, Min);
                mag = Math.Min(mag, Max);
                mag -= Min;
                mag /= (Max - Min) * Diameter;
                points[i] = mag * hemisphere[i] + Ctr;
            }
            Hare.Geometry.Topology T = Utilities.Geometry.GeoHemiSphere(5,1);

            for (int i = 0; i < T.Vertex_Count; i++) T.Set_Vertex(i, points[i]);
            T.Finish_Topology();
            return T;
        }
    }
}