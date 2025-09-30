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
using Pachyderm_Acoustic.Pach_Graphics;
using Pachyderm_Acoustic.Environment;
using Hare.Geometry;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Pachyderm_Acoustic
{
    namespace Utilities
    {
        public class FileIO
        {
            ///// <summary>
            ///// Writes a Pac1 file. [Useable by scripts and interface components alike.]
            ///// </summary>
            ///// <param name="Direct_Data">Array of Completed Direct Sound Simulations Required</param>
            ///// <param name="IS_Data">Array of Completed Image Source Simulations. Enter null if opted out.</param>
            ///// <param name="Receiver">Array of Completed Ray-Tracing Simulation Receivers. Enter null if opted out.</param>
            //public static void Write_Pac1(Direct_Sound[] Direct_Data, ImageSourceData[] IS_Data = null, Environment.Receiver_Bank[] Receiver = null)
            //{
            //    Eto.Forms.SaveFileDialog sf = new Eto.Forms.SaveFileDialog();
            //    sf.CurrentFilter = ".pac1";
            //    sf.Filters.Add("Pachyderm Ray Data file (*.pac1)|*.pac1|" + "All Files|");

            //    if (sf.ShowDialog() == Eto.Forms.DialogResult.Ok)
            //    {
            //        Write_Pac1(sf.FileName, Direct_Data, IS_Data, Receiver);
            //    }
            //}

            /// <summary>
            /// Writes a Pac1 file. [Useable by scripts and interface components alike.]
            /// </summary>
            /// <param name="filename">The location of the final saved file...</param>
            /// <param name="Direct_Data">Array of Completed Direct Sound Simulations Required</param>
            /// <param name="IS_Data">Array of Completed Image Source Simulations. Enter null if opted out.</param>
            /// <param name="Receiver">Array of Completed Ray-Tracing Simulation Receivers. Enter null if opted out.</param>
            public static void Write_Pac1(string filename, Direct_Sound[] Direct_Data, ImageSourceData[] IS_Data = null, Environment.Receiver_Bank[] Receiver = null)
            {
                if (Direct_Data == null && IS_Data == null && IS_Data == null && Receiver != null)
                {
                    throw new Exception("There is no simulated data to save.");
                }

                int maxAttempts = 3;
                int currentAttempt = 0;
                bool success = false;

                while (!success && currentAttempt < maxAttempts)
                {
                    System.IO.BinaryWriter sw = null;
                    try
                    {
                        sw = new System.IO.BinaryWriter(System.IO.File.Open(filename, System.IO.FileMode.Create));

                        //1. Date & Time
                        sw.Write(System.DateTime.Now.ToString());
                        //2. Library Version... if less than 1.1, assume only 1 source.
                        sw.Write(PachTools.Version());
                        //3. Cut off Time (seconds) and SampleRate
                        sw.Write((double)Receiver[0].CO_Time);
                        sw.Write(Receiver[0].SampleRate);
                        //3.1 Third Octave? (version >= 2.7) true for 1/3 octave, false for octave band.
                        if (double.Parse(PachTools.Version().Substring(0, 2)) >= 2.7) sw.Write(Direct_Data[0].Io[0].Length > 8);
                        //4.0 Source Count(int)
                        Hare.Geometry.Point[] SRC = new Hare.Geometry.Point[Direct_Data.Length];
                        for (int i = 0; i < Direct_Data.Length; i++) SRC[i] = Direct_Data[i].Src_Origin;

                        sw.Write(SRC.Length);
                        for (int i = 0; i < SRC.Length; i++)
                        {
                            //4.1 Source Location x (double)    
                            sw.Write(SRC[i].x);
                            //4.2 Source Location y (double)
                            sw.Write(SRC[i].y);
                            //4.3 Source Location z (double)
                            sw.Write(SRC[i].z);
                        }
                        //5. No of Receivers
                        sw.Write(Receiver[0].Rec_List.Length);

                        //6. Write the coordinates of each receiver point
                        //6b. Write the environmental characteristics at each receiver point (Rho * C); V2.0 only...
                        for (int q = 0; q < Receiver[0].Rec_List.Length; q++)
                        {
                            sw.Write(Receiver[0].Rec_List[q].Origin.x);
                            sw.Write(Receiver[0].Rec_List[q].Origin.y);
                            sw.Write(Receiver[0].Rec_List[q].Origin.z);
                            sw.Write(Receiver[0].Rec_List[q].Rho_C);
                        }

                        for (int s = 0; s < SRC.Length; s++)
                        {
                            if (Direct_Data != null)
                            {
                                //7. Write Direct Sound Data
                                Direct_Data[s].Write_Data(ref sw);
                            }

                            if (IS_Data != null && IS_Data[0] != null)
                            {
                                //8. Write Image Source Sound Data
                                IS_Data[s].Write_Data(ref sw);
                            }

                            if (Receiver != null)
                            {
                                //9. Write Ray Traced Sound Data
                                Receiver[s].Write_Data(ref sw);
                            }
                        }
                        sw.Write("End");
                        success = true;
                    }
                    catch (System.UnauthorizedAccessException ex)
                    {
                        // Permission issues - not worth retrying
                        Eto.Forms.MessageBox.Show("Access Denied",
                            $"Unable to save file due to permission issues.\n\nPlease ensure you have write access to the folder or try saving to a different location.\n\nError details: {ex.Message}");
                        break;
                    }
                    catch (System.IO.IOException ex)
                    {
                        // Potentially transient IO issues - try again
                        currentAttempt++;
                        if (currentAttempt >= maxAttempts)
                        {
                            Eto.Forms.MessageBox.Show("File Save Error",
                                $"Could not save the file after {maxAttempts} attempts.\n\nThe file may be in use by another application or the storage device may be full.\n\nError details: {ex.Message}");
                        }
                        else
                        {
                            // Wait a moment before retrying
                            System.Threading.Thread.Sleep(500);
                        }
                    }
                    catch (Exception ex)
                    {
                        // Other unexpected errors
                        Eto.Forms.MessageBox.Show("Unexpected Error",
                            $"An unexpected error occurred while saving the file:\n\n{ex.Message}\n\nPlease report this issue to the developers.");
                        break;
                    }
                    finally
                    {
                        // Ensure the writer is closed and disposed
                        if (sw != null)
                        {
                            sw.Close();
                            sw.Dispose();
                        }
                    }
                }
            }

            public static bool Read_Pac1(string filename, ref Direct_Sound[] Direct_Data, ref ImageSourceData[] IS_Data, ref Environment.Receiver_Bank[] Receiver, IProgressFeedback VB = null)
            {
                System.IO.BinaryReader sr = new System.IO.BinaryReader(System.IO.File.Open(filename, System.IO.FileMode.Open));
                try
                {
                    //1. Date & Time
                    string Savedate = sr.ReadString();
                    //2. Plugin Version
                    string Pach_version = sr.ReadString();
                    //3. Cut off Time and SampleRate
                    double CO_TIME = sr.ReadDouble();
                    int SampleRate = sr.ReadInt32();
                    //3.1 Third Octave? (version >= 2.7) true for 1/3 octave, false for octave band.
                    bool Third_Octave = double.Parse(Pach_version.Substring(0, 3)) >= 2.7 ? sr.ReadBoolean() : false;
                    //4. Source Count          
                    int SrcCt = 1;
                    string vers_shallow = Pach_version.Substring(0, 3);
                    if (double.Parse(vers_shallow) >= 1.1) SrcCt = sr.ReadInt32();
                    //4.1 Source Location x                   
                    //4.2 Source Location y
                    //4.3 Source Location z
                    Hare.Geometry.Point[] SrcPt = new Hare.Geometry.Point[SrcCt];
                    for (int s = 0; s < SrcCt; s++) SrcPt[s] = new Hare.Geometry.Point(sr.ReadDouble(), sr.ReadDouble(), sr.ReadDouble());
                    //5. No of Receivers
                    int Rec_Ct = sr.ReadInt32();
                    //6. Write the coordinates of each receiver point
                    //6b. Write the environmental characteristics at each receiver point (Rho * C); V2.0 only...
                    Hare.Geometry.Point[] Recs = new Hare.Geometry.Point[Rec_Ct];
                    double[] Rho_C = new double[Rec_Ct];
                    for (int q = 0; q < Rec_Ct; q++)
                    {
                        Recs[q] = new Hare.Geometry.Point(sr.ReadDouble(), sr.ReadDouble(), sr.ReadDouble());
                        if (double.Parse(Pach_version.Substring(0, 3)) >= 2.0) Rho_C[q] = sr.ReadDouble();
                        else Rho_C[q] = 400;
                    }

                    Direct_Data = new Direct_Sound[SrcCt];
                    IS_Data = new ImageSourceData[SrcCt];
                    Receiver = new Environment.Receiver_Bank[SrcCt];

                    int DDCT = 0;
                    int ISCT = 0;
                    int RTCT = 0;
                    do
                    {
                        string readin = sr.ReadString();
                        switch (readin)
                        {
                            case "Direct_Sound":
                            case "Direct_Sound w sourcedata":
                                //9. Read Direct Sound Data
                                Direct_Data[DDCT] = Direct_Sound.Read_Data(ref sr, Recs, SrcPt[DDCT], Rho_C, Third_Octave, Pach_version);
                                Direct_Data[DDCT].CO_Time = CO_TIME;
                                Direct_Data[DDCT].SampleFreq = (int)SampleRate;
                                DDCT++;
                                //float s1 = sr.ReadSingle();
                                //float s2 = sr.ReadSingle();
                                break;
                            case "Image-Source_Data":
                                //10. Read Image Source Sound Data
                                IS_Data[ISCT] = ImageSourceData.Read_Data(ref sr, Recs.Length, Direct_Data[DDCT - 1], false, Rho_C, ISCT, Pach_version, VB);
                                ISCT++;
                                break;
                            case "Ray-Traced_Data":
                                //11. Read Ray Traced Sound Data
                                Receiver[RTCT] = Environment.Receiver_Bank.Read_Data(ref sr, Direct_Data[RTCT].SWL, Rec_Ct, Recs, Rho_C, Direct_Data[RTCT].Delay_ms, ref SampleRate, Pach_version, VB);
                                RTCT++;
                                break;
                            case "End":
                            case "":
                                sr.Close();
                                return true;
                        }
                    } while (true);
                }
                catch (System.Exception X)
                {
                    sr.Close();
                    Eto.Forms.MessageBox.Show("File Read Failed...", String.Format("Results file was corrupt or incomplete. We apologize for this inconvenience. Please report this to the software author. It will be much appreciated. \r\n Exception Message: {0}. \r\n Method: {1}", X.Message, X.TargetSite));
                    return false;
                }
            }

            /// <summary>
            /// Writes pachyderm mapping file.
            /// </summary>
            /// <param name="filename">The location the new file is to be written to.</param>
            /// <param name="Rec_List">The list of receivers to be written.</param>
            public static void Write_pachm(string filename, PachMapReceiver[] Rec_List)
            {
                int maxAttempts = 3;
                int currentAttempt = 0;
                bool success = false;

                while (!success && currentAttempt < maxAttempts)
                {
                    System.IO.BinaryWriter sw = null;
                    try
                    {
                        sw = new System.IO.BinaryWriter(System.IO.File.Open(filename, System.IO.FileMode.Create));

                        //1. Write calculation type. (string)  
                        sw.Write(Rec_List[0].Data_Type());
                        Boolean Directional = Rec_List[0].Data_Type() == "Type;Map_Data";
                        //2. Write the number of samples in each histogram. (int)
                        sw.Write((UInt32)Rec_List[0].SampleCT);
                        //3. Write the sample rate. (int) 
                        sw.Write((UInt32)Rec_List[0].SampleRate);
                        //4. Write the number of Receivers (int)
                        int Rec_Ct = Rec_List[0].Rec_List.Length;
                        sw.Write((UInt32)Rec_Ct);
                        //4.5 Announce the Version
                        sw.Write("Version");
                        sw.Write(PachTools.Version());
                        //5. Announce that the following data pertains to the form of the analysis mesh. (string)
                        sw.Write("Mesh Information");
                        //6. Announce Mesh Vertices (string)
                        sw.Write("Mesh Vertices");
                        //Write the number of vertices & faces (int) (int)
                        sw.Write((UInt32)Rec_List[0].Map_Mesh.Vertex_Count);
                        sw.Write((UInt32)Rec_List[0].Map_Mesh.Polygon_Count);

                        for (int i = 0; i < Rec_List[0].Map_Mesh.Vertex_Count; i++)
                        {
                            //Write Vertex: (double) (double) (double)
                            sw.Write(Rec_List[0].Map_Mesh[i].x);
                            sw.Write(Rec_List[0].Map_Mesh[i].y);
                            sw.Write(Rec_List[0].Map_Mesh[i].z);
                        }
                        //7. Announce Mesh Faces (string)
                        sw.Write("Mesh Faces");
                        for (int i = 0; i < Rec_List[0].Map_Mesh.Polygon_Count; i++)
                        {
                            // Write mesh vertex indices: (int) (int) (int) (int)
                            sw.Write((UInt32)Rec_List[0].Map_Mesh.Polys[i].Points[0].index);
                            sw.Write((UInt32)Rec_List[0].Map_Mesh.Polys[i].Points[1].index);
                            sw.Write((UInt32)Rec_List[0].Map_Mesh.Polys[i].Points[2].index);
                            if (Rec_List[0].Map_Mesh.Polys[i].VertextCT > 3) sw.Write((UInt32)Rec_List[0].Map_Mesh.Polys[i].Points[3].index);
                            else sw.Write((UInt32)Rec_List[0].Map_Mesh.Polys[i].Points[2].index);
                        }
                        //7.5: Announce the number of sources.
                        sw.Write("SourceswLoc");
                        sw.Write(Rec_List.Length);
                        //7.5a: Announce the Type of Source
                        for (int i = 0; i < Rec_List.Length; i++)
                        {
                            ///////////////////////
                            sw.Write(Rec_List[i].Src.x);
                            sw.Write(Rec_List[i].Src.y);
                            sw.Write(Rec_List[i].Src.z);
                            ///////////////////////
                            sw.Write(Rec_List[i].SrcType);
                            sw.Write(Rec_List[i].delay_ms);//v.2.0.0.1
                        }

                        //8. Announce that the following data pertains to the receiver histograms (string)
                        sw.Write("Receiver Hit Data");
                        //8a. Announce whether or not data is linked to vertices rather than faces (bool)
                        sw.Write(Rec_List[0].Rec_Vertex);

                        for (int s = 0; s < Rec_List.Length; s++)
                        {
                            for (int i = 0; i < Rec_Ct; i++)
                            {
                                //Write Receiver Index (int)
                                sw.Write((UInt32)i);
                                //Write the direct sound arrival time.
                                sw.Write((Rec_List[s].Rec_List[i] as PachMapReceiver.Map_Receiver).Direct_Time);
                                //Write Impedance of Air
                                sw.Write(Rec_List[0].Rec_List[i].Rho_C);

                                for (int Octave = 0; Octave < 8; Octave++)
                                {
                                    //Write Octave (int)
                                    sw.Write((UInt32)Octave);
                                    double[] Hist = Rec_List[s].Rec_List[i].GetEnergyHistogram(Octave);
                                    for (int e = 0; e < Rec_List[s].SampleCT; e++)
                                    {
                                        //Write each energy value in the histogram (double)...
                                        sw.Write(Hist[e]);
                                        sw.Write(Rec_List[s].Rec_List[i].Pressure_rms(e, Octave));//v.2.5
                                                                                                  //Write each directional value in the histogram (double) (double) (double);
                                        if (Directional)
                                        {
                                            Hare.Geometry.Vector DirPos = Rec_List[s].Directions_Pos(Octave, e, i);
                                            Hare.Geometry.Vector DirNeg = Rec_List[s].Directions_Neg(Octave, e, i);
                                            sw.Write(DirPos.dx);
                                            sw.Write(DirPos.dy);
                                            sw.Write(DirPos.dz);
                                            sw.Write(DirNeg.dx);
                                            sw.Write(DirNeg.dy);
                                            sw.Write(DirNeg.dz);
                                        }
                                    }
                                }
                                sw.Write("End_Receiver_Hits");
                            }
                        }
                        sw.Write("End_of_File");
                        success = true;
                    }
                    catch (System.UnauthorizedAccessException ex)
                    {
                        // Permission issues - not worth retrying
                        Eto.Forms.MessageBox.Show("Access Denied",
                            $"Unable to save mapping file due to permission issues.\n\nPlease ensure you have write access to the folder or try saving to a different location.\n\nError details: {ex.Message}");
                        break;
                    }
                    catch (System.IO.IOException ex)
                    {
                        // Potentially transient IO issues - try again
                        currentAttempt++;
                        if (currentAttempt >= maxAttempts)
                        {
                            Eto.Forms.MessageBox.Show("File Save Error",
                                $"Could not save the mapping file after {maxAttempts} attempts.\n\nThe file may be in use by another application or the storage device may be full.\n\nError details: {ex.Message}");
                        }
                        else
                        {
                            // Wait a moment before retrying
                            System.Threading.Thread.Sleep(500);
                        }
                    }
                    catch (System.Exception ex)
                    {
                        // Other unexpected errors
                        Eto.Forms.MessageBox.Show("Unexpected Error",
                            $"An unexpected error occurred while saving the mapping file:\n\n{ex.Message}\n\nPlease report this issue to the developers.");
                        break;
                    }
                    finally
                    {
                        // Ensure the writer is closed and disposed
                        if (sw != null)
                        {
                            sw.Close();
                            sw.Dispose();
                        }
                    }
                }
            }

            /// <summary>
            /// reads a file and populates the map receiver instance.
            /// </summary>
            /// <returns></returns>
            public static bool Read_pachm(string path, ref PachMapReceiver[] Map)
            {
                System.IO.BinaryReader sr = new System.IO.BinaryReader(System.IO.File.Open(path, System.IO.FileMode.Open));
                //1. Write calculation type. (string)
                string CalcType = sr.ReadString();
                if (CalcType != "Type;Map_Data" && CalcType != "Type;Map_Data_NoDir") throw new Exception("Map Data File Expected");
                bool Directional = (CalcType == "Type;Map_Data");

                //2. Write the number of samples in each histogram. (int)
                int SampleCT = (int)sr.ReadUInt32();
                //3. Write the sample rate. (int) 
                int SampleRate = (int)sr.ReadUInt32();
                //4. Write the number of Receivers (int)
                int Rec_CT = (int)sr.ReadUInt32();
                //4.5 Write the version number
                double version = 1.1;
                double rev = 0;
                //5. Announce that the following data pertains to the form of the analysis mesh. (string)
                int s_ct = 1;
                Hare.Geometry.Topology Map_Mesh = new Hare.Geometry.Topology();
                Map = new PachMapReceiver[1];
                Map[0] = new PachMapReceiver();
                //double[] Rho_C = null;

                double[] delay;

                do
                {
                    switch (sr.ReadString())
                    {
                        case "Version":
                            //Pach1.7 = Versioning functionality added.
                            string v = sr.ReadString();
                            version = double.Parse(v.Substring(0, 3));
                            rev = int.Parse(v.Split(new char[1] { '.' })[3]);
                            break;
                        case "Mesh Information":
                            //6. Announce Mesh Vertices (string)
                            //Write the number of vertices & faces (int) (int)
                            if (sr.ReadString() != "Mesh Vertices") throw new Exception("Mesh Vertices Expected");

                            int VC = (int)sr.ReadUInt32();
                            int FC = (int)sr.ReadUInt32();

                            Hare.Geometry.Point[] Vertices = new Hare.Geometry.Point[VC];

                            for (int i = 0; i < VC; i++)
                            {
                                //Write Vertex: (double) (double) (double)
                                //Map_Mesh.Vertices.Add(new Rhino.Geometry.Point3d(sr.ReadSingle(), sr.ReadSingle(), sr.ReadSingle()));
                                Vertices[i] = new Hare.Geometry.Point(sr.ReadDouble(), sr.ReadDouble(), sr.ReadDouble());
                            }

                            //7. Announce Mesh Faces (string)
                            if (sr.ReadString() != "Mesh Faces") throw new Exception("Mesh Faces Expected");

                            Hare.Geometry.Point[][] Polys = new Hare.Geometry.Point[FC][];
                            for (int i = 0; i < FC; i++)
                            {
                                // Write mesh vertex indices: (int) (int) (int) (int)
                                int F0 = (int)sr.ReadUInt32(), F1 = (int)sr.ReadUInt32(), F2 = (int)sr.ReadUInt32(), F3 = (int)sr.ReadUInt32();

                                if (F2 == F3) Polys[i] = new Hare.Geometry.Point[3] { Vertices[F0], Vertices[F1], Vertices[F2] };// Map_Mesh.Faces.AddFace(F0, F1, f2);
                                else Polys[i] = new Hare.Geometry.Point[4] { Vertices[F0], Vertices[F1], Vertices[F2], Vertices[F3] };//Map_Mesh.Faces.AddFace(F0, F1, f2, F3);
                            }

                            Map_Mesh = new Hare.Geometry.Topology(Polys);
                            break;
                        case "Sources":
                            //7.5: Announce the number of sources.
                            s_ct = sr.ReadInt32();
                            delay = new double[s_ct];
                            Map = new PachMapReceiver[s_ct];
                            //7.5a Announce the type of source.

                            for (int s = 0; s < s_ct; s++)
                            {
                                Map[s] = new PachMapReceiver();
                                Map[s].CutOffTime = (double)SampleCT / (double)SampleRate;
                                Map[s].SampleCT = SampleCT;
                                Map[s].SampleRate = SampleRate;
                                Map[s].Map_Mesh = Map_Mesh;
                                Map[s].Rec_List = new PachMapReceiver.Map_Receiver[Rec_CT];
                                Map[s].SrcType = sr.ReadString();
                                //4.4 Source delay (ms)
                                if (version > 2.0 || (version == 2.0 && rev >= 1))
                                {
                                    delay[s] = sr.ReadDouble();
                                }
                            }
                            break;
                        case "SourceswLoc":
                            //7.5: Announce the number of sources.
                            s_ct = sr.ReadInt32();
                            delay = new double[s_ct];
                            Map = new PachMapReceiver[s_ct];
                            //7.5a Announce the type of source.

                            for (int s = 0; s < s_ct; s++)
                            {
                                Map[s] = new PachMapReceiver();
                                Map[s].CutOffTime = (double)SampleCT / (double)SampleRate * 1000;
                                Map[s].SampleCT = SampleCT;
                                Map[s].SampleRate = SampleRate;
                                Map[s].Map_Mesh = Map_Mesh;
                                Map[s].Rec_List = new PachMapReceiver.Map_Receiver[Rec_CT];
                                Map[s].Src = new Hare.Geometry.Point(sr.ReadDouble(), sr.ReadDouble(), sr.ReadDouble());
                                Map[s].SrcType = sr.ReadString();
                                //4.4 Source delay (ms)
                                if (version > 2.0 || (version == 2.0 && rev >= 1))
                                {
                                    delay[s] = sr.ReadDouble();
                                }
                            }
                            break;
                        case "Receiver Hit Data":
                            if (Map[0] == null)
                            {
                                Map = new PachMapReceiver[1];
                                Map[0] = new PachMapReceiver();
                                Map[0].CutOffTime = (double)SampleCT / (double)SampleRate;
                                Map[0].SampleRate = SampleRate;
                                Map[0].SampleCT = SampleCT;
                                Map[0].Map_Mesh = Map_Mesh;
                                Map[0].Rec_List = new PachMapReceiver.Map_Receiver[Rec_CT];
                                Map[0].SrcType = "Geodesic";
                            }

                            //8. Announce that the following data pertains to the receiver histograms (string)                        
                            //8a. Announce whether or not data is linked to vertices rather than faces (bool)
                            bool vert_Receiver = sr.ReadBoolean();
                            for (int s = 0; s < s_ct; s++)
                            {
                                Map[s].Rec_Vertex = vert_Receiver;
                                for (int i = 0; i < Map[s].Rec_List.Length; i++)
                                {
                                    //for version 1.7 and up, write direct sound arrival time.
                                    //Write Receiver Index (int)
                                    int j = (int)sr.ReadUInt32();
                                    //Write Direct Sound Arrival Time.
                                    double Direct_Time;
                                    if (version >= 1.7) Direct_Time = sr.ReadDouble(); else Direct_Time = (Map[s].Src - Map[s].Rec_List[i].Origin).Length() / 343f;
                                    //Write Impedance of Air
                                    double Rho_C = version >= 2.0 ? sr.ReadDouble() : 400;

                                    if (vert_Receiver)
                                    {
                                        Map[s].Rec_List[i] = new PachMapReceiver.Map_Receiver(Map_Mesh[i], Map[s].Src, Direct_Time, Rho_C, i, SampleRate, SampleCT, Directional);
                                    }
                                    else
                                    {
                                        Hare.Geometry.Point RecLoc = Map_Mesh.Polys[i].Points[0] + Map_Mesh.Polys[i].Points[1] + Map_Mesh.Polys[i].Points[2];
                                        if (Map_Mesh.Polys[i].VertextCT > 3)
                                        {
                                            RecLoc += Map_Mesh.Polys[i].Points[3];
                                            RecLoc /= 4;
                                        }
                                        else RecLoc /= 3;
                                        Map[s].Rec_List[i] = new PachMapReceiver.Map_Receiver(RecLoc, Map[s].Src, Direct_Time, Rho_C, i, SampleRate, SampleCT, Directional);
                                    }

                                    for (int Octave = 0; Octave < 8; Octave++)
                                    {
                                        //Write Octave (int)
                                        int Oct_out = (int)sr.ReadUInt32();
                                        if (Oct_out != Octave) throw new Exception(string.Format("Octave {0} Expected", Octave));
                                        double[] Hist = Map[s].Rec_List[i].GetEnergyHistogram(Octave);
                                        if (Directional)
                                        {
                                            if (version < 2.5)
                                            {
                                                for (int e = 0; e < SampleCT; e++)
                                                    Map[s].Rec_List[i].Combine_Sample(e, sr.ReadDouble(), 0, new Hare.Geometry.Vector(sr.ReadSingle(), sr.ReadSingle(), sr.ReadSingle()), new Hare.Geometry.Vector(sr.ReadSingle(), sr.ReadSingle(), sr.ReadSingle()), Octave);
                                            }
                                            else
                                            {
                                                for (int e = 0; e < SampleCT; e++)
                                                    Map[s].Rec_List[i].Combine_Sample(e, sr.ReadDouble(), sr.ReadDouble(), new Hare.Geometry.Vector(sr.ReadSingle(), sr.ReadSingle(), sr.ReadSingle()), new Hare.Geometry.Vector(sr.ReadSingle(), sr.ReadSingle(), sr.ReadSingle()), Octave);
                                            }
                                        }
                                        else
                                        {
                                            if (version < 2.5)
                                            {
                                                for (int e = 0; e < SampleCT; e++)
                                                    Map[s].Rec_List[i].Combine_Sample(e, sr.ReadDouble(), 0, new Hare.Geometry.Vector(0, 0, 0), new Hare.Geometry.Vector(0, 0, 0), Octave);
                                            }
                                            else
                                            {
                                                for (int e = 0; e < SampleCT; e++)
                                                    Map[s].Rec_List[i].Combine_Sample(e, sr.ReadDouble(), sr.ReadDouble(), new Hare.Geometry.Vector(0, 0, 0), new Hare.Geometry.Vector(0, 0, 0), Octave);
                                            }
                                        }
                                    }
                                    if (sr.ReadString() != "End_Receiver_Hits") throw new Exception("End of Receiver Hits Expected");
                                }
                            }
                            break;
                        case "End_of_File":
                            sr.Close();
                            return true;
                    }
                } while (true);
                throw new Exception("Unsuccessful Read");
            }

            /// <summary>
            /// Reads a 3D model from a text file and outputs the Polygon Scene object.
            /// </summary>
            /// <param name="filePath">Path to the model file</param>
            /// <returns>Parameters ready for the Construct method</returns>
            public static Polygon_Scene ReadModelFile(string filePath)
            {
                // Vertices need to be in a dictionary because they're referenced by ID
                Dictionary<int, Point> vertices = new Dictionary<int, Point>();
                
                // Simple lists for materials and their properties
                List<string> materialNames = new List<string>();
                List<double[]> absorptionCoeffs = new List<double[]>();
                List<double[]> scatteringCoeffs = new List<double[]>();
                List<double[]> transmissionCoeffs = new List<double[]>();
                
                // Lists to store object and polygon data
                List<int> objectMaterialIds = new List<int>();
                List<List<List<int>>> objectPolygonVertices = new List<List<List<int>>>();
                List<List<bool>> objectPolygonCurved = new List<List<bool>>();
                List<List<double[]>> objectPolygonKurvature = new List<List<double[]>>();
                List<List<Vector[]>> objectPolygonFrameAxes = new List<List<Vector[]>>();
                
                using (StreamReader reader = new StreamReader(filePath))
                {
                    string line;
                    string section = null;
                    int currentObjectIndex = -1;

                    while ((line = reader.ReadLine()) != null)
                    {
                        line = line.Trim();

                        // Skip empty lines and comments
                        if (string.IsNullOrWhiteSpace(line) || line.StartsWith("#"))
                            continue;

                        // Check for section headers
                        if (line.StartsWith("[VERTICES]", StringComparison.OrdinalIgnoreCase))
                        {
                            section = "VERTICES";
                            continue;
                        }
                        else if (line.StartsWith("[MATERIALS]", StringComparison.OrdinalIgnoreCase))
                        {
                            section = "MATERIALS";
                            continue;
                        }
                        else if (line.StartsWith("[OBJECT]", StringComparison.OrdinalIgnoreCase))
                        {
                            section = "OBJECT";
                            // Extract material ID from [OBJECT MaterialID]
                            string[] parts = line.Substring(7).Trim().Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                            if (parts.Length > 0)
                            {
                                int materialId = int.Parse(parts[0]) - 1; // Adjust to 0-based index
                                objectMaterialIds.Add(materialId);
                                objectPolygonVertices.Add(new List<List<int>>());
                                objectPolygonCurved.Add(new List<bool>());
                                objectPolygonKurvature.Add(new List<double[]>());
                                objectPolygonFrameAxes.Add(new List<Vector[]>());
                                currentObjectIndex = objectMaterialIds.Count - 1;
                            }
                            continue;
                        }

                        // Process each section
                        if (section == "VERTICES")
                        {
                            string[] parts = line.Split(new[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                            if (parts.Length >= 4)
                            {
                                int id = int.Parse(parts[0]);
                                double x = double.Parse(parts[1]);
                                double y = double.Parse(parts[2]);
                                double z = double.Parse(parts[3]);
                                vertices[id] = new Point(x, y, z);
                            }
                        }
                        else if (section == "MATERIALS")
                        {
                            string[] parts = line.Split(new[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                            if (parts.Length >= 2)
                            {
                                int id = int.Parse(parts[0]);
                                string name = parts[1];
                                
                                // Create arrays for coefficients
                                double[] abs = new double[8];
                                double[] scat = new double[8];
                                double[] trans = new double[8];
                                
                                // Read material properties until [END]
                                string propLine;
                                while ((propLine = reader.ReadLine()) != null)
                                {
                                    propLine = propLine.Trim();

                                    if (string.IsNullOrWhiteSpace(propLine) || propLine.StartsWith("#"))
                                        continue;

                                    if (propLine.StartsWith("[END]", StringComparison.OrdinalIgnoreCase))
                                        break;

                                    string[] propParts = propLine.Split(new[] { '=' }, 2, StringSplitOptions.RemoveEmptyEntries);
                                    if (propParts.Length == 2)
                                    {
                                        string propName = propParts[0].Trim().ToUpperInvariant();
                                        string propValue = propParts[1].Trim();

                                        switch (propName)
                                        {
                                            case "ABSORPTION":
                                                abs = propValue.Split(',').Select(s => double.Parse(s.Trim())).ToArray();
                                                break;
                                            case "SCATTERING":
                                                scat = propValue.Split(',').Select(s => double.Parse(s.Trim())).ToArray();
                                                break;
                                            case "TRANSMISSION":
                                                trans = propValue.Split(',').Select(s => double.Parse(s.Trim())).ToArray();
                                                break;
                                        }
                                    }
                                }
                                
                                // Ensure arrays are the right size
                                if (abs.Length != 8)
                                    Array.Resize(ref abs, 8);
                                if (scat.Length != 8)
                                    Array.Resize(ref scat, 8);
                                if (trans.Length != 8)
                                    Array.Resize(ref trans, 8);
                                
                                // Make sure id-1 is a valid index (materials are 1-based in the file)
                                while (materialNames.Count <= id - 1)
                                {
                                    materialNames.Add(string.Empty);
                                    absorptionCoeffs.Add(new double[8]);
                                    scatteringCoeffs.Add(new double[8]);
                                    transmissionCoeffs.Add(new double[8]);
                                }
                                
                                // Store the material data
                                materialNames[id - 1] = name;
                                absorptionCoeffs[id - 1] = abs;
                                scatteringCoeffs[id - 1] = scat;
                                transmissionCoeffs[id - 1] = trans;
                            }
                        }
                        else if (section == "OBJECT" && currentObjectIndex >= 0)
                        {
                            if (line.StartsWith("POLY", StringComparison.OrdinalIgnoreCase))
                            {
                                string data = line.Substring(4).Trim();
                                string[] parts = data.Split(new[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);

                                if (parts.Length >= 1)
                                {
                                    // Parse vertex IDs
                                    List<int> vertexIds = parts[0].Split(',').Select(int.Parse).ToList();
                                    bool isCurved = false;
                                    double[] kurvature = new double[2];
                                    Vector[] frameAxes = new Vector[2] { new Vector(), new Vector() };
                                    
                                    // Parse optional properties
                                    for (int i = 1; i < parts.Length; i++)
                                    {
                                        string prop = parts[i];
                                        if (prop.StartsWith("CURVE=", StringComparison.OrdinalIgnoreCase))
                                        {
                                            string curveData = prop.Substring(6);
                                            string[] values = curveData.Split(',');
                                            if (values.Length == 2)
                                            {
                                                isCurved = true;
                                                kurvature[0] = double.Parse(values[0]);
                                                kurvature[1] = double.Parse(values[1]);
                                            }
                                        }
                                        else if (prop.StartsWith("FRAME=", StringComparison.OrdinalIgnoreCase))
                                        {
                                            string frameData = prop.Substring(6);
                                            string[] values = frameData.Split(',');
                                            if (values.Length == 6)
                                            {
                                                isCurved = true;
                                                frameAxes[0] = new Vector(
                                                    double.Parse(values[0]),
                                                    double.Parse(values[1]),
                                                    double.Parse(values[2])
                                                );
                                                frameAxes[1] = new Vector(
                                                    double.Parse(values[3]),
                                                    double.Parse(values[4]),
                                                    double.Parse(values[5])
                                                );
                                            }
                                        }
                                    }
                                    
                                    // Validate curvature and frame axes consistency
                                    if ((kurvature[0] != 0 || kurvature[1] != 0) && 
                                        (frameAxes[0].IsZeroVector() || frameAxes[1].IsZeroVector()))
                                    {
                                        throw new FormatException("Polygon has curvature but missing frame axes");
                                    }

                                    if ((!frameAxes[0].IsZeroVector() || !frameAxes[1].IsZeroVector()) && 
                                        (kurvature[0] == 0 && kurvature[1] == 0))
                                    {
                                        throw new FormatException("Polygon has frame axes but missing curvature");
                                    }
                                    
                                    // Store the polygon data
                                    objectPolygonVertices[currentObjectIndex].Add(vertexIds);
                                    objectPolygonCurved[currentObjectIndex].Add(isCurved);
                                    objectPolygonKurvature[currentObjectIndex].Add(kurvature);
                                    objectPolygonFrameAxes[currentObjectIndex].Add(frameAxes);
                                }
                            }
                        }
                    }
                }
                
                // Create the output structures for Construct()
                int objectCount = objectMaterialIds.Count;
                
                Point[][][] model = new Point[objectCount][][];
                List<Environment.Material> materials = new List<Environment.Material>();
                List<Scattering> scatterings = new List<Scattering>();
                List<double[]> transmissions = new List<double[]>();
                bool[] isCurved_List = new bool[objectCount];
                double[][][] kurvatures = new double[objectCount][][];
                Vector[][][] frameAxes_List = new Vector[objectCount][][];
                
                // Build the model
                for (int objIndex = 0; objIndex < objectCount; objIndex++)
                {
                    int materialId = objectMaterialIds[objIndex];
                    
                    // Add material properties
                    materials.Add(new Environment.Basic_Material(absorptionCoeffs[materialId]));
                    scatterings.Add(new Lambert_Scattering(scatteringCoeffs[materialId]));
                    transmissions.Add(transmissionCoeffs[materialId]);
                    
                    // Check if any polygon in this object is curved
                    isCurved_List[objIndex] = objectPolygonCurved[objIndex].Any(c => c);

                    // Create polygon arrays
                    int polyCount = objectPolygonVertices[objIndex].Count;
                    model[objIndex] = new Point[polyCount][];
                    kurvatures[objIndex] = new double[polyCount][];
                    frameAxes_List[objIndex] = new Vector[polyCount][];
                    
                    for (int polyIndex = 0; polyIndex < polyCount; polyIndex++)
                    {
                        // Convert vertex IDs to Points
                        model[objIndex][polyIndex] = objectPolygonVertices[objIndex][polyIndex]
                            .Select(id => vertices[id])
                            .ToArray();
                        
                        // Set curvature and frame axes
                        kurvatures[objIndex][polyIndex] = objectPolygonKurvature[objIndex][polyIndex];
                        frameAxes_List[objIndex][polyIndex] = objectPolygonFrameAxes[objIndex][polyIndex];
                    }
                }

                // Create and construct the Polygon_Scene with standard defaults
                // Temperature: 20°C, Humidity: 50%, Pressure: 101325 Pa (standard atmospheric pressure)
                // Air_Choice: 0 (default), EdgeCorrection: true, IsAcoustic: true
                Polygon_Scene scene = new Polygon_Scene(20.0, 50.0, 101325.0, 0, true, false);
                scene.Construct(model, materials, scatterings, transmissions, isCurved_List, kurvatures, frameAxes_List);

                return scene;
            }
        }
    }
}