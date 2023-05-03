//'Pachyderm-Acoustic: Geometrical Acoustics for Rhinoceros (GPL) by Arthur van der Harten 
//' 
//'This file is part of Pachyderm-Acoustic. 
//' 
//'Copyright (c) 2008-2023, Arthur van der Harten 
//'Pachyderm-Acoustic is free software; you can redistribute it and/or modify 
//'it under the terms of the GNU General Public License as published 
//'by the Free Software Foundation; either version 3 of the License, or 
//'(at your option) any later version. 
//'Pachyderm-Acoustic is distributed in the hope that it will be useful, 
//'but WITHOUT ANY WARRANTY; without even the implied warranty of 
//'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
//'GNU General Public License for more details. 
//' 
//'You should have received a copy of the GNU General Public 
//'License along with Pachyderm-Acoustic; if not, write to the Free Software 
//'Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. 

using System;
using System.Collections.Generic;
using Hare.Geometry;

namespace Pachyderm_Acoustic
{
    public class Balloon
    {
        public Topology[] m_HareMesh;
        public Topology m_DisplayMesh;
        public float CurrentAlt, CurrentAzi, CurrentAxi;
        protected Hare.Geometry.Point CurrentPos = new Hare.Geometry.Point(0, 0, 0);
        protected string SWL;
        protected int Type;
        public string[] Code;

        public Balloon()
        { }

        public Balloon(string[] ballooncodes, Hare.Geometry.Point Center) //, double[] SWL
        {
            string[][][] balloons = new string[8][][];
            Code = ballooncodes;

            int minstable = 8;

            if (ballooncodes.Length == 8)
            {
                for (int oct = 0; oct < 8; oct++)
                {
                    if (ballooncodes[oct] == null) continue;
                    minstable = Math.Min(minstable, oct);
                    string[] lines = ballooncodes[oct].Split(new char[1] { ';' }, StringSplitOptions.RemoveEmptyEntries);
                    balloons[oct] = new string[lines.Length][];
                    for (int i = 0; i < lines.Length; i++)
                    {
                        balloons[oct][i] = lines[i].Split(new char[1] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                    }
                }
            }
            int umax, vmax;

            if (balloons[minstable][0].Length == 19) //[36].Count == 36 || balloons[4].Count == 19 || balloons[4].Count == 10)
            {
                Type = 0;
                umax = 19;
                vmax = 36;
            }
            else if (balloons[minstable][0].Length == 37)//[72].Count == 72 || balloons[4].Count == 37 || balloons[4].Count == 19)
            {
                Type = 1;
                umax = 37;
                vmax = 72;
            }
            else throw new Exception("Balloon Resolution - the number of samples in the balloon is non-standard. This is not supported in Pachyderm at this time...");

            for(int oct = 0; oct < 8; oct++)
            {
                if (balloons[oct] == null)
                {
                    balloons[oct] = new string[vmax][];
                    for(int v = 0; v < vmax; v++)
                    {
                        balloons[oct][v] = new string[umax];
                        for (int u = 0; u < umax; u++) balloons[oct][v][u] = "40";
                    }
                }
            }

            m_HareMesh = new Topology[8];
            for (int oct = 0; oct < 8; oct++)
            {
                Vector[,] Magnitude = new Vector[umax, vmax];
                double[,] mag = new double[umax, vmax];
                double Theta, Phi;
                //int idx = 0;

                for (int v = 0; v < balloons[oct].Length; v++)
                {
                    for (int u = 0; u < balloons[oct][v].Length; u++)
                    {
                            Theta = u * System.Math.PI / (umax - 1);
                            Phi = 2 * v * System.Math.PI / (vmax) + System.Math.PI / 2;
                            Magnitude[u, v] = new Vector(Math.Sin(Theta) * Math.Cos(Phi), Math.Cos(Theta), Math.Sin(Theta) * Math.Sin(Phi));
                            Magnitude[u, v].Normalize();
                            
                            double swlmag = (-double.Parse(balloons[oct][v][u]) + 60);
                            if (swlmag < 0) swlmag = 0;
                            Magnitude[u, v] *= swlmag;
                            mag[u, v] = swlmag;

                            if (balloons[oct].Length == vmax / 4 + 1)//Quarter data
                            {
                                int v_2 = vmax / 2 + 1;
                                Magnitude[(vmax - u) % vmax, v] = Magnitude[u, v];
                                Magnitude[v_2 - u, v] = Magnitude[u, v];
                                Magnitude[v_2 + u, v] = Magnitude[u, v];
                                mag[(vmax - u) % 36, v] = mag[u, v];
                                mag[v_2 - u, v] = mag[u, v];
                                mag[v_2 + u, v] = mag[u, v];
                            }
                            if (balloons[oct].Length == vmax / 2 + 1)//Half data
                            {
                                Magnitude[(vmax - u) % vmax, v] = Magnitude[u, v];
                                mag[(vmax - u) % vmax, v] = mag[u, v];
                            }
                    }
                }

                if (Code[oct] == null)
                {
                    Code[oct] = "";
                    for (int v = 0; v < balloons[oct].Length; v++)
                    {
                        for (int u = 0; u < balloons[oct][v].Length; u++)
                        {
                            Code[oct] += "0.0 ";
                        }
                        Code[oct] += ';';
                    }
                }
                //Create a mesh of the points...
                List<Hare.Geometry.Point[]> list = new List<Hare.Geometry.Point[]>();
                double Minx = double.MaxValue, Miny = double.MaxValue, Minz = double.MaxValue, Maxx = double.MinValue, Maxy = double.MinValue, Maxz = double.MinValue;

                for (int u = 0; u < umax - 1; u++)
                {
                    for (int v = 0; v < vmax; v++)
                    {
                        Hare.Geometry.Point[] Poly = new Hare.Geometry.Point[4];
                        Poly[0] = new Hare.Geometry.Point(Magnitude[u, v].x, Magnitude[u, v].y, Magnitude[u, v].z);
                        Poly[1] = new Hare.Geometry.Point(Magnitude[u, (v + 1) % vmax].x, Magnitude[u, (v + 1) % vmax].y, Magnitude[u, (v + 1) % vmax].z);
                        Poly[2] = new Hare.Geometry.Point(Magnitude[u + 1, (v + 1) % vmax].x, Magnitude[u + 1, (v + 1) % vmax].y, Magnitude[u + 1, (v + 1) % vmax].z);
                        Poly[3] = new Hare.Geometry.Point(Magnitude[u + 1, v].x, Magnitude[u + 1, v].y, Magnitude[u + 1, v].z);
                        list.Add(Poly);

                        foreach (Hare.Geometry.Point p in Poly)
                        {
                            if (p.x < Minx) Minx = p.x;
                            if (p.y < Miny) Miny = p.y;
                            if (p.z < Minz) Minz = p.z;
                            if (p.x > Maxx) Maxx = p.x;
                            if (p.y > Maxy) Maxy = p.y;
                            if (p.z > Maxz) Maxz = p.z;
                        }
                    }
                }

                m_HareMesh[oct] = new Topology(new Hare.Geometry.Point(Minx, Miny, Minz), new Hare.Geometry.Point(Maxx, Maxy, Maxz));
                foreach (Hare.Geometry.Point[] Poly in list) m_HareMesh[oct].Add_Polygon(Poly);
                if (oct == 4)
                {
                    m_DisplayMesh = new Topology(new Hare.Geometry.Point(Minx, Miny, Minz), new Hare.Geometry.Point(Maxx, Maxy, Maxz));
                    foreach (Hare.Geometry.Point[] Poly in list) m_DisplayMesh.Add_Polygon(Poly);
                    for (int i = 0; i < m_DisplayMesh.Vertex_Count; i++) m_DisplayMesh.Set_Vertex(i, m_DisplayMesh[i] / 90);
                }
            }
            Update_Position(Center);
        }

        public virtual Topology[] Balloons(double[] spl_values)
        {
            Topology[] Balloon = new Topology[8];

            int umax, vmax;
            switch (Type)
            {
                case 0:
                    umax = 19;
                    vmax = 36;
                    break;
                case 1:
                    umax = 37;
                    vmax = 72;
                    break;
                default:
                    throw new Exception("Balloon type not valid for this version of CLF.");
            }

            for (int oct = 0; oct < 8; oct++)
            {
                string[] values;
                if (Code[oct] != "")
                {
                    values = Code[oct].Split(';');
                }
                else
                {
                    values = new string[umax * vmax];
                    for (int i = 0; i < values.Length; i++)
                    {
                        values[i] = "0";
                    }
                }

                Vector[,] Magnitude = new Vector[umax, vmax];
                double Theta, Phi;
                for (int v = 0; v < vmax; v++)
                {
                    string[] valuesplit = values[v].Split(new char[1] { ' ' },  StringSplitOptions.RemoveEmptyEntries);
                    for (int u = 0; u < umax; u++)
                    {
                        if (double.IsInfinity(spl_values[oct]) || double.IsNaN(spl_values[oct]))
                        {
                            Magnitude[u, v] = new Vector();
                        }
                        else
                        {
                            Theta = u * System.Math.PI / (umax - 1);
                            Phi = 2 * v * System.Math.PI / (vmax) + System.Math.PI / 2;
                            Magnitude[u, v] = new Vector(Math.Sin(Theta) * Math.Cos(Phi), Math.Cos(Theta), Math.Sin(Theta) * Math.Sin(Phi));
                            Magnitude[u, v].Normalize();
                            Magnitude[u, v] *= (-double.Parse(valuesplit[u]) + spl_values[oct]);
                        }
                    }
                }

                for (int u = 0; u < umax; u++)
                {
                    for (int v = 0; v < vmax; v++)
                    {
                        double roll = CurrentAxi * Math.PI / 180;
                        double yaw = CurrentAlt * Math.PI / 180;
                        double pitch = CurrentAzi * Math.PI / 180;

                        double x = Magnitude[u, v].x;
                        Magnitude[u, v].x = x * Math.Cos(roll) - Magnitude[u, v].z * Math.Sin(roll);
                        Magnitude[u, v].z = x * Math.Sin(roll) + Magnitude[u, v].z * Math.Cos(roll);
                        double y = Magnitude[u, v].y;
                        Magnitude[u, v].y = y * Math.Cos(yaw) - Magnitude[u, v].z * Math.Sin(yaw);
                        Magnitude[u, v].z = y * Math.Sin(yaw) + Magnitude[u, v].z * Math.Cos(yaw);
                        x = Magnitude[u, v].x;
                        Magnitude[u, v].x = x * Math.Cos(pitch) - Magnitude[u, v].y * Math.Sin(pitch);
                        Magnitude[u, v].y = x * Math.Sin(pitch) + Magnitude[u, v].y * Math.Cos(pitch);
                    }
                }

                List<Hare.Geometry.Point[]> list = new List<Hare.Geometry.Point[]>();
                //Create a mesh of the points...
                for (int u = 0; u < umax - 1; u++)
                {
                    for (int v = 0; v < vmax; v++)
                    {
                        Hare.Geometry.Point[] Poly = new Hare.Geometry.Point[3];
                        Poly[0] = new Hare.Geometry.Point(Magnitude[u, v].x, Magnitude[u, v].y, Magnitude[u, v].z);
                        Poly[1] = new Hare.Geometry.Point(Magnitude[u, (v + 1) % vmax].x, Magnitude[u, (v + 1) % vmax].y, Magnitude[u, (v + 1) % vmax].z);
                        Poly[2] = new Hare.Geometry.Point(Magnitude[u + 1, v].x, Magnitude[u + 1, v].y, Magnitude[u + 1, v].z);
                        list.Add(Poly);
                        //Balloon[oct - 1].Add_Polygon(Poly);

                        Poly = new Hare.Geometry.Point[3];
                        Poly[0] = new Hare.Geometry.Point(Magnitude[u, (v + 1) % vmax].x, Magnitude[u, (v + 1) % vmax].y, Magnitude[u, (v + 1) % vmax].z);
                        Poly[1] = new Hare.Geometry.Point(Magnitude[u + 1, (v + 1) % vmax].x, Magnitude[u + 1, (v + 1) % vmax].y, Magnitude[u + 1, (v + 1) % vmax].z);
                        Poly[2] = new Hare.Geometry.Point(Magnitude[u + 1, v].x, Magnitude[u + 1, v].y, Magnitude[u + 1, v].z);
                        list.Add(Poly);
                        //Balloon[oct - 1].Add_Polygon(Poly);
                    }
                }

                double Minx = double.MaxValue, Miny = double.MaxValue, Minz = double.MaxValue, Maxx = double.MinValue, Maxy = double.MinValue, Maxz = double.MinValue;
                foreach (Hare.Geometry.Point[] p in list)
                {
                    foreach (Hare.Geometry.Point p0 in p)
                    {
                        if (p0.x < Minx) Minx = p0.x;
                        if (p0.y < Miny) Miny = p0.y;
                        if (p0.z < Minz) Minz = p0.z;
                        if (p0.x > Maxx) Maxx = p0.x;
                        if (p0.y > Maxy) Maxy = p0.y;
                        if (p0.z > Maxz) Maxz = p0.z;
                    }
                }

                Balloon[oct] = new Topology(new Hare.Geometry.Point(Minx, Miny, Minz), new Hare.Geometry.Point(Maxx, Maxy, Maxz));
                foreach (Hare.Geometry.Point[] p in list) Balloon[oct].Add_Polygon(p);
            }
            Balloon[0].Finish_Topology(new List<Point>());
            Balloon[1].Finish_Topology(new List<Point>());
            Balloon[2].Finish_Topology(new List<Point>());
            Balloon[3].Finish_Topology(new List<Point>());
            Balloon[4].Finish_Topology(new List<Point>());
            Balloon[5].Finish_Topology(new List<Point>());
            Balloon[6].Finish_Topology(new List<Point>());
            Balloon[7].Finish_Topology(new List<Point>());
            return Balloon;
        }

        public static string[] Read_Generic(string Path)
        {
            System.IO.StreamReader SR = new System.IO.StreamReader(Path);
            List<string[]>[] balloons = new List<string[]>[8];
            string[] Code = new string[10];
            string[] BANDS = new string[0];
            int[] BandIDS = new int[0];
            while (!SR.EndOfStream)
            {
                string line = SR.ReadLine();
                int oct = -1;
                switch (line)
                {
                    case "63":
                        oct = 0;
                        break;
                    case "125":
                        oct = 1;
                        break;
                    case "250":
                        oct = 2;
                        break;
                    case "500":
                        oct = 3;
                        break;
                    case "1k":
                        oct = 4;
                        break;
                    case "2k":
                        oct = 5;
                        break;
                    case "4k":
                        oct = 6;
                        break;
                    case "8k":
                        oct = 7;
                        break;
                    case "":
                        continue;
                    default:
                        line = line.Replace('<', ' ');
                        line = line.Replace('>', ' ');
                        string[] divided = line.Split(new char[1] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                        if (divided[0] == "OCTAVEBANDS")
                        {
                            BANDS = new string[divided.Length - 3];
                            Array.ConstrainedCopy(divided, 2, BANDS, 0, BANDS.Length);
                            BandIDS = new int[BANDS.Length];
                            for (int i = 0; i < BANDS.Length; i++)
                            {
                                BandIDS[i] = Utilities.PachTools.OctaveStr2Int(BANDS[i]);
                            }
                        }
                        if (BANDS.Length == 0) continue;

                        int codeindex;
                        if (divided[0] == "NOMSPL_AT_1M")
                        {
                            codeindex = 8;
                        }
                        else if (divided[0] == "MAXSPL_AT_1M")
                        {
                            codeindex = 9;
                        }
                        else
                        {
                            continue;
                        }

                        if (oct > 0) continue;
                        Code[codeindex] = "";
                        for (int testoct = 0; testoct < 8; testoct++)
                        {
                            bool band0 = true;
                            int bandoct = 0;
                            for (bandoct = 0; bandoct < BandIDS.Length; bandoct++)
                            {
                                if (testoct == BandIDS[bandoct])
                                {
                                    Code[codeindex] += double.Parse(divided[bandoct + 2]) + 11 + ";";
                                    band0 = false;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
                                    break;
                                }
                            }
                            if (band0) Code[codeindex] += "0;";
                        }
                        break;             
                }

                if (oct < 0) continue;

                string balloonline = SR.ReadLine();
                while (balloonline != "")
                {
                    Code[oct] += balloonline + ";";
                    balloonline = SR.ReadLine();
                }
            }

            return Code;
            ///
            /// Balloon63
            /// Balloon125
            /// Balloon250
            /// Balloon500
            /// Balloon1k
            /// Balloon2k
            /// Balloon4k
            /// Balloon8k
            /// Nominal SWL
            /// Max SWL
            ///
        }
        
        public void Update_Position()
        {
            Update_Position(CurrentPos);
        }

        public void Update_Position(Hare.Geometry.Point Center)
        {
            for (int i = 0; i < m_HareMesh[4].Vertex_Count; i++)
            {
                Point Po = m_HareMesh[4][i];
                Point P = new Point(Po.x, Po.y, Po.z)/90; //new Point((float)(Po.x + Center.x - CurrentPos.x), (float)(Po.y + Center.y - CurrentPos.y), (float)(Po.z + Center.z - CurrentPos.z));
                m_DisplayMesh.Set_Vertex(i, P);
            }
            Update_Aim();
            for (int i = 0; i < m_HareMesh[4].Vertex_Count; i++)
            {
                Point Po = m_DisplayMesh[i];
                Point P = new Point((float)(Po.x + Center.x), (float)(Po.y + Center.y), (float)(Po.z + Center.z)); //new Point((float)(Po.x + Center.x - CurrentPos.x), (float)(Po.y + Center.y - CurrentPos.y), (float)(Po.z + Center.z - CurrentPos.z));
                m_DisplayMesh.Set_Vertex(i, P);
            }
            CurrentPos = Center;
        }

        public void Update_Aim()
        {
            for (int i = 0; i < m_DisplayMesh.Vertex_Count; i++)
            {
                Point P = new Point((float)m_DisplayMesh[i].x, m_DisplayMesh[i].y, m_DisplayMesh[i].z);//(m_DisplayMesh[i].x - CurrentPos.x), (float)(m_DisplayMesh[i].y - CurrentPos.y), (float)(m_DisplayMesh[i].z - CurrentPos.z));
                double roll = CurrentAxi * Math.PI / 180;
                double yaw = CurrentAlt * Math.PI / 180;
                double pitch = CurrentAzi * Math.PI / 180;

                double x = m_DisplayMesh[i].x;
                P = new Point((float)(x * Math.Cos(roll) - P.z * Math.Sin(roll)), P.y, (float)(x * Math.Sin(roll) + P.z * Math.Cos(roll)));
                double y = P.y;
                P = new Point(P.x, (float)(y * Math.Cos(yaw) - P.z * Math.Sin(yaw)), (float)(y * Math.Sin(yaw) + P.z * Math.Cos(yaw)));
                x = P.x;
                P = new Point((float)(x * Math.Cos(pitch) - P.y * Math.Sin(pitch)), (float)(x * Math.Sin(pitch) + P.y * Math.Cos(pitch)), P.z);

                m_DisplayMesh.Set_Vertex(i, P);// new Point((float)(P.x + CurrentPos.x), (float)(P.y + CurrentPos.y), (float)(P.z + CurrentPos.z)));
            }
        }
    }

    public class Speaker_Balloon: Balloon
    {
        /// <summary>
        /// Original mesh used for the final source object.
        /// </summary>
        public string[] code;

        public Speaker_Balloon(string[] Ballooncode_in, string SWL_in, int Type_in, Hare.Geometry.Point Center)//int sym_in,
        {
            code = Ballooncode_in;
            Type = Type_in;
            if (SWL_in != "" || SWL_in != null)
            {
                SWL = SWL_in;
            }
            else
            {
                SWL_in = "120; 120; 120; 120; 120; 120; 120; 120";
            }
            string[] swl = SWL.Split(';');
            double[] swl_values = new double[8];

            Vector upper = new Vector(0, 1, 0);

            int umax, vmax;
            switch (Type)
            {
                case 0:
                    umax = 19;
                    vmax = 36;
                    break;
                case 1:
                    umax = 37;
                    vmax = 72;
                    break;
                default:
                    throw new Exception();
            }

            m_HareMesh = new Topology[8];
            for (int oct = 1; oct < 9; oct++)
            {
                string[] values;
                swl_values[oct - 1] = double.Parse(swl[oct - 1]);
                if (code[oct - 1] != "")
                {
                    values = code[oct - 1].Split(';');
                }
                else
                {
                    values = new string[umax * vmax];
                    for (int i = 0; i < values.Length; i++)
                    {
                        values[i] = "0";
                    }
                }

                Vector[,] Magnitude = new Vector[umax, vmax];
                double Theta, Phi;
                int idx = 0;

                for (int v = 0; v < vmax; v++)
                {
                    for (int u = 0; u < umax; u++)
                    {
                        if (double.IsInfinity(swl_values[oct - 1]) || double.IsNaN(swl_values[oct - 1]))
                        {
                            Magnitude[u, v] = new Vector();
                        }
                        else
                        {
                            Theta = u * System.Math.PI / (umax - 1);
                            Phi = 2 * v * System.Math.PI / (vmax) + System.Math.PI / 2;
                            Magnitude[u, v] = new Vector(Math.Sin(Theta) * Math.Cos(Phi), Math.Cos(Theta), Math.Sin(Theta) * Math.Sin(Phi));
                            Magnitude[u, v].Normalize();
                            double swlmag = (double.Parse(values[idx]) + 60);
                            if (swlmag < 0) swlmag = 0;
                            Magnitude[u, v] *= swlmag;
                        }
                        idx++;
                    }
                }

                //Create a mesh of the points...
                List<Hare.Geometry.Point[]> list = new List<Hare.Geometry.Point[]>();
                double Minx = double.MaxValue, Miny = double.MaxValue, Minz = double.MaxValue, Maxx = double.MinValue, Maxy = double.MinValue, Maxz = double.MinValue;

                for (int u = 0; u < umax - 1; u++)
                {
                    for (int v = 0; v < vmax; v++)
                    {
                        Hare.Geometry.Point[] Poly = new Hare.Geometry.Point[4];
                        Poly[0] = Magnitude[u, v];
                        Poly[1] = Magnitude[u, (v + 1) % vmax];
                        Poly[2] = Magnitude[u + 1, (v + 1) % vmax];
                        Poly[3] = Magnitude[u + 1, v];
                        list.Add(Poly);

                        foreach (Hare.Geometry.Point p in Poly)
                        {
                            if (p.x < Minx) Minx = p.x;
                            if (p.y < Miny) Miny = p.y;
                            if (p.z < Minz) Minz = p.z;
                            if (p.x > Maxx) Maxx = p.x;
                            if (p.y > Maxy) Maxy = p.y;
                            if (p.z > Maxz) Maxz = p.z;
                        }
                    }
                }

                m_HareMesh[oct - 1] = new Topology(new Hare.Geometry.Point(Minx, Miny, Minz), new Hare.Geometry.Point(Maxx, Maxy, Maxz));
                foreach (Hare.Geometry.Point[] Poly in list) m_HareMesh[oct - 1].Add_Polygon(Poly);
                if (oct == 4)
                {
                    m_DisplayMesh = new Topology(new Hare.Geometry.Point(Minx, Miny, Minz), new Hare.Geometry.Point(Maxx, Maxy, Maxz));
                    foreach (Hare.Geometry.Point[] Poly in list) m_DisplayMesh.Add_Polygon(Poly);
                }
            }
            Update_Position(Center);
        }

        public override Topology[] Balloons(double[] spl_values)
        {
            Topology[] Balloon = new Topology[8];

            int umax, vmax;
            switch (Type)
            {
                case 0:
                    umax = 19;
                    vmax = 36;
                    break;
                case 1:
                    umax = 37;
                    vmax = 72;
                    break;
                default:
                    throw new Exception("Balloon type not valid for this version of CLF.");
            }

            for (int oct = 1; oct < 9; oct++)
            {
                string[] values;
                if (code[oct - 1] != "")
                {
                    values = code[oct - 1].Split(';');
                }
                else
                {
                    values = new string[umax * vmax];
                    for (int i = 0; i < values.Length; i++)
                    {
                        values[i] = "0";
                    }
                }

                Vector[,] Magnitude = new Vector[umax, vmax];
                double Theta, Phi;
                int idx = 0;
                for (int v = 0; v < vmax; v++)
                {
                    for (int u = 0; u < umax; u++)
                    {
                        if (double.IsInfinity(spl_values[oct - 1]) || double.IsNaN(spl_values[oct - 1]))
                        {
                            Magnitude[u, v] = new Vector();
                        }
                        else
                        {
                            Theta = u * System.Math.PI / (umax - 1);
                            Phi = 2 * v * System.Math.PI / (vmax) + System.Math.PI / 2;
                            Magnitude[u, v] = new Vector(Math.Sin(Theta) * Math.Cos(Phi), Math.Cos(Theta), Math.Sin(Theta) * Math.Sin(Phi));
                            Magnitude[u, v].Normalize();
                            Magnitude[u, v] *= (double.Parse(values[idx]) + spl_values[oct - 1]);
                        }
                        idx++;
                    }
                }

                for (int u = 0; u < umax; u++)
                {
                    for (int v = 0; v < vmax; v++)
                    {
                        double roll = CurrentAxi * Math.PI / 180;
                        double yaw = CurrentAlt * Math.PI / 180;
                        double pitch = CurrentAzi * Math.PI / 180;

                        double x = Magnitude[u, v].x;
                        Magnitude[u, v].x = x * Math.Cos(roll) - Magnitude[u, v].z * Math.Sin(roll);
                        Magnitude[u, v].z = x * Math.Sin(roll) + Magnitude[u, v].z * Math.Cos(roll);
                        double y = Magnitude[u, v].y;
                        Magnitude[u, v].y = y * Math.Cos(yaw) - Magnitude[u, v].z * Math.Sin(yaw);
                        Magnitude[u, v].z = y * Math.Sin(yaw) + Magnitude[u, v].z * Math.Cos(yaw);
                        x = Magnitude[u, v].x;
                        Magnitude[u, v].x = x * Math.Cos(pitch) - Magnitude[u, v].y * Math.Sin(pitch);
                        Magnitude[u, v].y = x * Math.Sin(pitch) + Magnitude[u, v].y * Math.Cos(pitch);
                    }
                }

                List<Hare.Geometry.Point[]> list = new List<Hare.Geometry.Point[]>();
                //Create a mesh of the points...
                for (int u = 0; u < umax - 1; u++)
                {
                    for (int v = 0; v < vmax; v++)
                    {
                        Hare.Geometry.Point[] Poly = new Hare.Geometry.Point[3];
                        Poly[0] = Magnitude[u, v];
                        Poly[1] = Magnitude[u, (v + 1) % vmax];
                        Poly[2] = Magnitude[u + 1, v];
                        list.Add(Poly);
                        //Balloon[oct - 1].Add_Polygon(Poly);

                        Poly = new Hare.Geometry.Point[3];
                        Poly[0] = Magnitude[u, (v + 1) % vmax];
                        Poly[1] = Magnitude[u + 1, (v + 1) % vmax];
                        Poly[2] = Magnitude[u + 1, v];
                        list.Add(Poly);
                        //Balloon[oct - 1].Add_Polygon(Poly);
                    }
                }

                double Minx = double.MaxValue, Miny = double.MaxValue, Minz = double.MaxValue, Maxx = double.MinValue, Maxy = double.MinValue, Maxz = double.MinValue;
                foreach (Hare.Geometry.Point[] p in list)
                {
                    foreach (Hare.Geometry.Point p0 in p)
                    {
                        if (p0.x < Minx) Minx = p0.x;
                        if (p0.y < Miny) Miny = p0.y;
                        if (p0.z < Minz) Minz = p0.z;
                        if (p0.x > Maxx) Maxx = p0.x;
                        if (p0.y > Maxy) Maxy = p0.y;
                        if (p0.z > Maxz) Maxz = p0.z;
                    }
                }

                Balloon[oct - 1] = new Topology(new Hare.Geometry.Point(Minx, Miny, Minz), new Hare.Geometry.Point(Maxx, Maxy, Maxz));
                foreach (Hare.Geometry.Point[] p in list) Balloon[oct - 1].Add_Polygon(p);
            }
            Balloon[0].Finish_Topology(new List<Point>());
            Balloon[1].Finish_Topology(new List<Point>());
            Balloon[2].Finish_Topology(new List<Point>());
            Balloon[3].Finish_Topology(new List<Point>());
            Balloon[4].Finish_Topology(new List<Point>());
            Balloon[5].Finish_Topology(new List<Point>());
            Balloon[6].Finish_Topology(new List<Point>());
            Balloon[7].Finish_Topology(new List<Point>());
            return Balloon;
        }
    }
}