//'Pachyderm-Acoustic: Geometrical Acoustics for Rhinoceros (GPL)   
//' 
//'This file is part of Pachyderm-Acoustic. 
//' 
//'Copyright (c) 2008-2025, Open Research in Acoustical Science and Education, Inc. - a 501(c)3 nonprofit 
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

using Hare.Geometry;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;

namespace Pachyderm_Acoustic
{
    namespace Environment
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

                for (int oct = 0; oct < 8; oct++)
                {
                    if (balloons[oct] == null)
                    {
                        balloons[oct] = new string[vmax][];
                        for (int v = 0; v < vmax; v++)
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
                            Poly[0] = new Hare.Geometry.Point(Magnitude[u, v].dx, Magnitude[u, v].dy, Magnitude[u, v].dz);
                            Poly[1] = new Hare.Geometry.Point(Magnitude[u, (v + 1) % vmax].dx, Magnitude[u, (v + 1) % vmax].dy, Magnitude[u, (v + 1) % vmax].dz);
                            Poly[2] = new Hare.Geometry.Point(Magnitude[u + 1, (v + 1) % vmax].dx, Magnitude[u + 1, (v + 1) % vmax].dy, Magnitude[u + 1, (v + 1) % vmax].dz);
                            Poly[3] = new Hare.Geometry.Point(Magnitude[u + 1, v].dx, Magnitude[u + 1, v].dy, Magnitude[u + 1, v].dz);
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
                    //foreach (Hare.Geometry.Point[] Poly in list) m_HareMesh[oct].Add_Polygon(Poly);
                    m_HareMesh[oct].Build_Topology(list.ToArray());
                    m_HareMesh[oct].Finish_Topology();
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
                        string[] valuesplit = values[v].Split(new char[1] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
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

                            double x = Magnitude[u, v].dx;
                            Magnitude[u, v].dx = x * Math.Cos(roll) - Magnitude[u, v].dz * Math.Sin(roll);
                            Magnitude[u, v].dz = x * Math.Sin(roll) + Magnitude[u, v].dz * Math.Cos(roll);
                            double y = Magnitude[u, v].dy;
                            Magnitude[u, v].dy = y * Math.Cos(yaw) - Magnitude[u, v].dz * Math.Sin(yaw);
                            Magnitude[u, v].dz = y * Math.Sin(yaw) + Magnitude[u, v].dz * Math.Cos(yaw);
                            x = Magnitude[u, v].dx;
                            Magnitude[u, v].dx = x * Math.Cos(pitch) - Magnitude[u, v].dy * Math.Sin(pitch);
                            Magnitude[u, v].dy = x * Math.Sin(pitch) + Magnitude[u, v].dy * Math.Cos(pitch);
                        }
                    }

                    List<Hare.Geometry.Point[]> list = new List<Hare.Geometry.Point[]>();
                    //Create a mesh of the points...
                    for (int u = 0; u < umax - 1; u++)
                    {
                        for (int v = 0; v < vmax; v++)
                        {
                            Hare.Geometry.Point[] Poly = new Hare.Geometry.Point[3];
                            Poly[0] = new Hare.Geometry.Point(Magnitude[u, v].dx, Magnitude[u, v].dy, Magnitude[u, v].dz);
                            Poly[1] = new Hare.Geometry.Point(Magnitude[u, (v + 1) % vmax].dx, Magnitude[u, (v + 1) % vmax].dy, Magnitude[u, (v + 1) % vmax].dz);
                            Poly[2] = new Hare.Geometry.Point(Magnitude[u + 1, v].dx, Magnitude[u + 1, v].dy, Magnitude[u + 1, v].dz);
                            list.Add(Poly);
                            //Balloon[oct - 1].Add_Polygon(Poly);

                            Poly = new Hare.Geometry.Point[3];
                            Poly[0] = new Hare.Geometry.Point(Magnitude[u, (v + 1) % vmax].dx, Magnitude[u, (v + 1) % vmax].dy, Magnitude[u, (v + 1) % vmax].dz);
                            Poly[1] = new Hare.Geometry.Point(Magnitude[u + 1, (v + 1) % vmax].dx, Magnitude[u + 1, (v + 1) % vmax].dy, Magnitude[u + 1, (v + 1) % vmax].dz);
                            Poly[2] = new Hare.Geometry.Point(Magnitude[u + 1, v].dx, Magnitude[u + 1, v].dy, Magnitude[u + 1, v].dz);
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

            public static string[] Read_BRAS_csv(string Path)
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
                    Point P = new Point(Po.x, Po.y, Po.z) / 90; //new Point((float)(Po.x + Center.x - CurrentPos.x), (float)(Po.y + Center.y - CurrentPos.y), (float)(Po.z + Center.z - CurrentPos.z));
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

        public class Speaker_Balloon : Balloon
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
                            Poly[0] = new Point(Magnitude[u, v].dx, Magnitude[u, v].dy, Magnitude[u, v].dz);
                            Poly[1] = new Point(Magnitude[u, (v + 1) % vmax].dx, Magnitude[u, (v + 1) % vmax].dy, Magnitude[u, (v + 1) % vmax].dz);
                            Poly[2] = new Point(Magnitude[u + 1, (v + 1) % vmax].dx, Magnitude[u + 1, (v + 1) % vmax].dy, Magnitude[u + 1, (v + 1) % vmax].dz);
                            Poly[3] = new Point(Magnitude[u + 1, v].dx, Magnitude[u + 1, v].dy, Magnitude[u + 1, v].dz);
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

                            double x = Magnitude[u, v].dx;
                            Magnitude[u, v].dx = x * Math.Cos(roll) - Magnitude[u, v].dz * Math.Sin(roll);
                            Magnitude[u, v].dz = x * Math.Sin(roll) + Magnitude[u, v].dz * Math.Cos(roll);
                            double y = Magnitude[u, v].dy;
                            Magnitude[u, v].dy = y * Math.Cos(yaw) - Magnitude[u, v].dz * Math.Sin(yaw);
                            Magnitude[u, v].dz = y * Math.Sin(yaw) + Magnitude[u, v].dz * Math.Cos(yaw);
                            x = Magnitude[u, v].dx;
                            Magnitude[u, v].dx = x * Math.Cos(pitch) - Magnitude[u, v].dy * Math.Sin(pitch);
                            Magnitude[u, v].dy = x * Math.Sin(pitch) + Magnitude[u, v].dy * Math.Cos(pitch);
                        }
                    }

                    List<Hare.Geometry.Point[]> list = new List<Hare.Geometry.Point[]>();
                    //Create a mesh of the points...
                    for (int u = 0; u < umax - 1; u++)
                    {
                        for (int v = 0; v < vmax; v++)
                        {
                            Hare.Geometry.Point[] Poly = new Hare.Geometry.Point[3];
                            Poly[0] = new Point(Magnitude[u, v].dx, Magnitude[u, v].dy, Magnitude[u, v].dz);
                            Poly[1] = new Point(Magnitude[u, (v + 1) % vmax].dx, Magnitude[u, (v + 1) % vmax].dy, Magnitude[u, (v + 1) % vmax].dz);
                            Poly[2] = new Point(Magnitude[u + 1, v].dx, Magnitude[u + 1, v].dy, Magnitude[u + 1, v].dz);
                            list.Add(Poly);
                            //Balloon[oct - 1].Add_Polygon(Poly);

                            Poly = new Hare.Geometry.Point[3];
                            Poly[0] = new Point(Magnitude[u, (v + 1) % vmax].dx, Magnitude[u, (v + 1) % vmax].dy, Magnitude[u, (v + 1) % vmax].dz);
                            Poly[1] = new Point(Magnitude[u + 1, (v + 1) % vmax].dx, Magnitude[u + 1, (v + 1) % vmax].dy, Magnitude[u + 1, (v + 1) % vmax].dz);
                            Poly[2] = new Point(Magnitude[u + 1, v].dx, Magnitude[u + 1, v].dy, Magnitude[u + 1, v].dz);
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

        /// <summary>
        /// Independent cabinet display geometry for CLF loudspeaker sources.
        /// Stores mesh vertices/faces and independent line segments using Hare geometry only.
        /// </summary>
        [Serializable]
        public class Cabinet_Geometry
        {
            public List<Point> Vertices = new List<Point>();
            public List<int[]> Faces = new List<int[]>();
            public List<Cabinet_Line> Lines = new List<Cabinet_Line>();

            public bool HasMesh
            {
                get { return Vertices != null && Vertices.Count > 0 && Faces != null && Faces.Count > 0; }
            }

            public bool HasLines
            {
                get { return Lines != null && Lines.Count > 0; }
            }

            public bool HasGeometry
            {
                get { return HasMesh || HasLines; }
            }

            public Cabinet_Geometry()
            {
            }

            public Cabinet_Geometry Clone()
            {
                Cabinet_Geometry g = new Cabinet_Geometry();

                foreach (Point p in Vertices)
                {
                    g.Vertices.Add(new Point(p.x, p.y, p.z));
                }

                foreach (int[] f in Faces)
                {
                    g.Faces.Add((int[])f.Clone());
                }

                foreach (Cabinet_Line l in Lines)
                {
                    g.Lines.Add(new Cabinet_Line(
                        new Point(l.A.x, l.A.y, l.A.z),
                        new Point(l.B.x, l.B.y, l.B.z)));
                }

                return g;
            }

            /// <summary>
            /// Applies scale, aiming, and translation to a cabinet geometry object.
            /// This follows the same rotation order used by Balloon.Update_Aim:
            /// axial rotation, altitude, azimuth.
            /// </summary>
            public Cabinet_Geometry Transform(
                Point origin,
                double altitudeDeg,
                double azimuthDeg,
                double axialRotationDeg,
                double unitScale)
            {
                Cabinet_Geometry g = new Cabinet_Geometry();

                foreach (Point p in Vertices)
                {
                    g.Vertices.Add(TransformPoint(
                        p,
                        origin,
                        altitudeDeg,
                        azimuthDeg,
                        axialRotationDeg,
                        unitScale));
                }

                foreach (int[] f in Faces)
                {
                    g.Faces.Add((int[])f.Clone());
                }

                foreach (Cabinet_Line l in Lines)
                {
                    g.Lines.Add(new Cabinet_Line(
                        TransformPoint(l.A, origin, altitudeDeg, azimuthDeg, axialRotationDeg, unitScale),
                        TransformPoint(l.B, origin, altitudeDeg, azimuthDeg, axialRotationDeg, unitScale)));
                }

                return g;
            }

            /// <summary>
            /// Converts mesh faces to a Hare Topology.
            /// Line-only geometry is intentionally not included.
            /// </summary>
            public Topology ToTopology()
            {
                if (!HasMesh) return null;

                List<Point[]> polygons = new List<Point[]>();

                foreach (int[] face in Faces)
                {
                    if (face == null || face.Length < 3) continue;

                    Point[] poly = new Point[face.Length];
                    bool valid = true;

                    for (int i = 0; i < face.Length; i++)
                    {
                        int index = face[i];

                        if (index < 0 || index >= Vertices.Count)
                        {
                            valid = false;
                            break;
                        }

                        Point p = Vertices[index];
                        poly[i] = new Point(p.x, p.y, p.z);
                    }

                    if (valid) polygons.Add(poly);
                }

                if (polygons.Count == 0) return null;

                double minx = double.MaxValue;
                double miny = double.MaxValue;
                double minz = double.MaxValue;
                double maxx = double.MinValue;
                double maxy = double.MinValue;
                double maxz = double.MinValue;

                foreach (Point[] poly in polygons)
                {
                    foreach (Point p in poly)
                    {
                        if (p.x < minx) minx = p.x;
                        if (p.y < miny) miny = p.y;
                        if (p.z < minz) minz = p.z;

                        if (p.x > maxx) maxx = p.x;
                        if (p.y > maxy) maxy = p.y;
                        if (p.z > maxz) maxz = p.z;
                    }
                }

                Topology top = new Topology(
                    new Point(minx, miny, minz),
                    new Point(maxx, maxy, maxz));

                foreach (Point[] poly in polygons)
                {
                    top.Add_Polygon(poly);
                }

                top.Finish_Topology();

                return top;
            }

            private static Point TransformPoint(
                Point p,
                Point origin,
                double altitudeDeg,
                double azimuthDeg,
                double axialRotationDeg,
                double unitScale)
            {
                double x = p.x * unitScale;
                double y = p.y * unitScale;
                double z = p.z * unitScale;

                double roll = axialRotationDeg * Math.PI / 180.0;
                double yaw = altitudeDeg * Math.PI / 180.0;
                double pitch = azimuthDeg * Math.PI / 180.0;

                // Match Balloon.Update_Aim() convention.
                // 1. Axial rotation: x-z plane.
                double x0 = x;
                x = x0 * Math.Cos(roll) - z * Math.Sin(roll);
                z = x0 * Math.Sin(roll) + z * Math.Cos(roll);

                // 2. Altitude rotation: y-z plane.
                double y0 = y;
                y = y0 * Math.Cos(yaw) - z * Math.Sin(yaw);
                z = y0 * Math.Sin(yaw) + z * Math.Cos(yaw);

                // 3. Azimuth rotation: x-y plane.
                x0 = x;
                x = x0 * Math.Cos(pitch) - y * Math.Sin(pitch);
                y = x0 * Math.Sin(pitch) + y * Math.Cos(pitch);

                return new Point(
                    origin.x + x,
                    origin.y + y,
                    origin.z + z);
            }
        }

        [Serializable]
        public struct Cabinet_Line
        {
            public Point A;
            public Point B;

            public Cabinet_Line(Point a, Point b)
            {
                A = a;
                B = b;
            }
        }

        /// <summary>
        /// Parses the compact cabinet geometry notation generated by the CLF reader.
        /// This class is Rhino-independent.
        /// </summary>
        public static class Cabinet_Geometry_Parser
        {
            public static Cabinet_Geometry Parse(
                string pointNotation,
                string faceNotation,
                string lineNotation)
            {
                Cabinet_Geometry g = new Cabinet_Geometry();

                Dictionary<int, Point> pointMap = ParsePoints(pointNotation);

                // Normalize external point IDs to compact vertex indices.
                Dictionary<int, int> idToVertexIndex = new Dictionary<int, int>();

                foreach (int id in pointMap.Keys.OrderBy(k => k))
                {
                    idToVertexIndex[id] = g.Vertices.Count;

                    Point p = pointMap[id];
                    g.Vertices.Add(new Point(p.x, p.y, p.z));
                }

                foreach (int[] face in ParseFaces(faceNotation))
                {
                    List<int> remapped = new List<int>();

                    foreach (int id in face)
                    {
                        int vertexIndex;
                        if (idToVertexIndex.TryGetValue(id, out vertexIndex))
                        {
                            remapped.Add(vertexIndex);
                        }
                    }

                    if (remapped.Count >= 3)
                    {
                        g.Faces.Add(remapped.ToArray());
                    }
                }

                foreach (Cabinet_Line line in ParseLines(lineNotation))
                {
                    g.Lines.Add(line);
                }

                return g;
            }

            public static Cabinet_Geometry ParseAndTransform(
                string pointNotation,
                string faceNotation,
                string lineNotation,
                Point origin,
                double altitudeDeg,
                double azimuthDeg,
                double axialRotationDeg,
                double unitScale)
            {
                Cabinet_Geometry local = Parse(pointNotation, faceNotation, lineNotation);

                return local.Transform(
                    origin,
                    altitudeDeg,
                    azimuthDeg,
                    axialRotationDeg,
                    unitScale);
            }

            public static Cabinet_Geometry ParseFromSourceUserStrings(
                Func<string, string> getUserString,
                Point origin,
                double unitScale)
            {
                if (getUserString == null) return new Cabinet_Geometry();

                string p = getUserString("CLF_CabinetPoints");
                string f = getUserString("CLF_CabinetFaces");
                string l = getUserString("CLF_CabinetLines");

                double alt = 0;
                double azi = 0;
                double axial = 0;

                string aiming = getUserString("Aiming");
                ParseAiming(aiming, out alt, out azi, out axial);

                return ParseAndTransform(
                    p,
                    f,
                    l,
                    origin,
                    alt,
                    azi,
                    axial,
                    unitScale);
            }

            public static void ParseAiming(
                string aiming,
                out double altitudeDeg,
                out double azimuthDeg,
                out double axialRotationDeg)
            {
                altitudeDeg = 0;
                azimuthDeg = 0;
                axialRotationDeg = 0;

                if (string.IsNullOrWhiteSpace(aiming)) return;

                string[] parts = aiming.Split(';');

                if (parts.Length > 0) TryParse(parts[0], out altitudeDeg);
                if (parts.Length > 1) TryParse(parts[1], out azimuthDeg);
                if (parts.Length > 2) TryParse(parts[2], out axialRotationDeg);
            }

            private static Dictionary<int, Point> ParsePoints(string notation)
            {
                Dictionary<int, Point> points = new Dictionary<int, Point>();

                if (string.IsNullOrWhiteSpace(notation)) return points;
                if (!notation.StartsWith("P|")) return points;

                string body = notation.Substring(2);

                string[] entries = body.Split(
                    new char[] { ';' },
                    StringSplitOptions.RemoveEmptyEntries);

                foreach (string entry in entries)
                {
                    string[] idAndPoint = entry.Split(':');
                    if (idAndPoint.Length != 2) continue;

                    int id;
                    if (!int.TryParse(
                        idAndPoint[0],
                        NumberStyles.Integer,
                        CultureInfo.InvariantCulture,
                        out id))
                    {
                        continue;
                    }

                    string[] xyz = idAndPoint[1].Split(',');
                    if (xyz.Length != 3) continue;

                    double x, y, z;
                    if (!TryParse(xyz[0], out x)) continue;
                    if (!TryParse(xyz[1], out y)) continue;
                    if (!TryParse(xyz[2], out z)) continue;

                    points[id] = new Point(x, y, z);
                }

                return points;
            }

            private static List<int[]> ParseFaces(string notation)
            {
                List<int[]> faces = new List<int[]>();

                if (string.IsNullOrWhiteSpace(notation)) return faces;
                if (!notation.StartsWith("F|")) return faces;

                string body = notation.Substring(2);

                string[] entries = body.Split(
                    new char[] { ';' },
                    StringSplitOptions.RemoveEmptyEntries);

                foreach (string entry in entries)
                {
                    string[] idAndFace = entry.Split(':');
                    if (idAndFace.Length != 2) continue;

                    string[] ids = idAndFace[1].Split(
                        new char[] { ',' },
                        StringSplitOptions.RemoveEmptyEntries);

                    List<int> face = new List<int>();

                    foreach (string idString in ids)
                    {
                        int id;
                        if (int.TryParse(
                            idString,
                            NumberStyles.Integer,
                            CultureInfo.InvariantCulture,
                            out id))
                        {
                            face.Add(id);
                        }
                    }

                    if (face.Count >= 3)
                    {
                        faces.Add(face.ToArray());
                    }
                }

                return faces;
            }

            private static List<Cabinet_Line> ParseLines(string notation)
            {
                List<Cabinet_Line> lines = new List<Cabinet_Line>();

                if (string.IsNullOrWhiteSpace(notation)) return lines;
                if (!notation.StartsWith("L|")) return lines;

                string body = notation.Substring(2);

                string[] entries = body.Split(
                    new char[] { ';' },
                    StringSplitOptions.RemoveEmptyEntries);

                foreach (string entry in entries)
                {
                    string[] v = entry.Split(',');
                    if (v.Length != 6) continue;

                    double x1, y1, z1, x2, y2, z2;

                    if (!TryParse(v[0], out x1)) continue;
                    if (!TryParse(v[1], out y1)) continue;
                    if (!TryParse(v[2], out z1)) continue;
                    if (!TryParse(v[3], out x2)) continue;
                    if (!TryParse(v[4], out y2)) continue;
                    if (!TryParse(v[5], out z2)) continue;

                    lines.Add(new Cabinet_Line(
                        new Point(x1, y1, z1),
                        new Point(x2, y2, z2)));
                }

                return lines;
            }

            private static bool TryParse(string s, out double value)
            {
                return double.TryParse(s, NumberStyles.Float | NumberStyles.AllowThousands, CultureInfo.InvariantCulture, out value);
            }
        }
    }
}