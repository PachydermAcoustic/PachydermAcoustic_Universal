using HDF.PInvoke;
using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Text;
using Hare.Geometry;
using Pachyderm_Acoustic.Utilities;

namespace Pachyderm_Acoustic
{
    namespace Audio
    {
        public class HRTF
        {
            public Topology T;
            public Voxel_Grid VG;
            Hare.Geometry.Point[][] P;
            int Fnum = 0;
            List<int>[] bins;
            List<Vector> emitterVectors;
            double[,] emitterdirs;
            double[][][] HRIR;
            int[] Face_Centroid_ID;
            Point origin;
            Random rand = new Random();
            int Fs;
            Vector[] PrincipalDirections;
            double[] Azi, Alt;

            public HRTF(string SofafilePath)
            {
                // Open the .SOFA file
                long fileId = H5F.open(SofafilePath, H5F.ACC_RDONLY);
                if (fileId < 0)
                {
                    Console.WriteLine("Error opening file.");
                    return;
                }

                try
                {
                    // Read the SOFA Conventions attribute
                    string Conv = ReadGlobalString(fileId, "Conventions");
                    string Conv_SOFA = ReadGlobalString(fileId, "SOFAConventions");
                    emitterdirs = ReadEmitterPosition(fileId);
                    //Read samplefrequency...
                    Fs = 48000;

                    //Create a reference sphere...
                    Pachyderm_Acoustic.Utilities.Geometry.GeoSphere(1);
                    PrincipalDirections = new Vector[T.Polygon_Count];
                    for (int i = 0; i < T.Polygon_Count; i++)
                    {
                        Vector d = new Vector(T.Polys[i].Centroid.x, T.Polys[i].Centroid.y, T.Polys[i].Centroid.z);
                        d.Normalize();
                        PrincipalDirections[i] = d;
                    }
                    double[,] view = ReadSourcePosition(fileId);

                    emitterVectors = new List<Vector>();
                    Alt = new double[T.Polygon_Count];
                    Azi = new double[T.Polygon_Count];
                    bins = new List<int>[T.Polygon_Count];
                    for (int i = 0; i < T.Polygon_Count; i++) bins[i] = new List<int>();

                    Face_Centroid_ID = new int[T.Polygon_Count];
                    for (int i = 0; i < T.Polygon_Count; i++) Face_Centroid_ID[i] = -1;
                    double[] mindist = new double[T.Polygon_Count];
                    for (int i = 0; i < T.Polygon_Count; i++) mindist[i] = double.MaxValue;

                    origin = new Point(0, 0, 0);
                    for (int i = 0; i < emitterdirs.GetLength(0); i++)
                    {
                        Vector d = new Vector(Math.Cos(Math.PI * emitterdirs[i, 0] / 180), Math.Sin(Math.PI * emitterdirs[i, 0] / 180), Math.Asin(Math.PI * emitterdirs[i, 1] / 180));
                        d.Normalize();
                        emitterVectors.Add(d);
                        X_Event X;
                        VG.Shoot(new Ray(origin, d, 0, rand.Next()), 0, out X);
                        if (X.Hit)
                        {
                            bins[X.Poly_id].Add(i);
                            double l = (T.Polygon_Centroid(X.Poly_id) - X.X_Point).Length();
                            if (l < mindist[X.Poly_id])
                            {
                                Azi[X.Poly_id] = emitterdirs[i, 0];
                                Alt[X.Poly_id] = emitterdirs[i, 1];
                                mindist[X.Poly_id] = l;
                                Face_Centroid_ID[X.Poly_id] = i;
                            }
                        }
                    }

                    // Read additional attributes or datasets here
                    // For example, reading the HRTF dataset
                    ReadHrtfDataset(fileId, "Data.IR");
                }
                finally
                {
                    H5F.close(fileId);
                }
            }

            public double[][] HeadRelatedIR(int face_id)
            {
                return HRIR[Face_Centroid_ID[face_id]];
            }

            public int SampleCt
            {
                get
                {
                    return HRIR[0][0].Length;
                }
            }

            public int DirsCT
            {
                get
                {
                    return T.Polys.Count;
                }
            }

            public Vector[] Directions
            {
                get { return PrincipalDirections; }
            }

            private double[,] ReadSourcePosition(long fileId)
            {
                double[,] retbuffer = null;
                string[] dId = new string[] { "ListenerView" };
                // Open the dataset
                foreach (string dataset in dId)
                {
                    long datasetId = H5D.open(fileId, dataset);

                    if (datasetId < 0)
                    {
                        Console.WriteLine("Error opening EmitterPosition dataset.");
                        return null;
                    }

                    try
                    {
                        // Get the datatype and dataspace of the dataset
                        long typeId = H5D.get_type(datasetId);
                        long spaceId = H5D.get_space(datasetId);

                        // Determine the number of elements (positions) in the dataspace
                        ulong[] dims = new ulong[2]; // Assuming EmitterPosition is 2D: [N x 3] for N emitters and 3 coordinates (x, y, z)
                        H5S.get_simple_extent_dims(spaceId, dims, null);

                        // Allocate a buffer for the double array
                        double[,] buffer = new double[dims[0], dims[1]];

                        // Read the dataset
                        GCHandle pinnedArray = GCHandle.Alloc(buffer, GCHandleType.Pinned);
                        try
                        {
                            IntPtr pointer = pinnedArray.AddrOfPinnedObject();
                            H5D.read(datasetId, typeId, H5S.ALL, H5S.ALL, H5P.DEFAULT, pointer);
                        }
                        finally
                        {
                            pinnedArray.Free();
                        }

                        if (buffer.GetLength(0) > 1) retbuffer = buffer; else continue;
                    }
                    finally
                    {
                        // Close resources
                        H5D.close(datasetId);
                    }
                    if (retbuffer != null) break;
                }

                return retbuffer;
            }

            private double[,] ReadEmitterPosition(long fileId)
            {
                double[,] retbuffer = null;
                string[] dId = new string[] { "EmitterPosition", "SourcePosition" };
                // Open the dataset
                foreach (string dataset in dId)
                {
                    long datasetId = H5D.open(fileId, dataset);

                    if (datasetId < 0)
                    {
                        Console.WriteLine("Error opening EmitterPosition dataset.");
                        return null;
                    }

                    try
                    {
                        // Get the datatype and dataspace of the dataset
                        long typeId = H5D.get_type(datasetId);
                        long spaceId = H5D.get_space(datasetId);

                        // Determine the number of elements (positions) in the dataspace
                        ulong[] dims = new ulong[2]; // Assuming EmitterPosition is 2D: [N x 3] for N emitters and 3 coordinates (x, y, z)
                        H5S.get_simple_extent_dims(spaceId, dims, null);

                        // Allocate a buffer for the double array
                        double[,] buffer = new double[dims[0], dims[1]];

                        // Read the dataset
                        GCHandle pinnedArray = GCHandle.Alloc(buffer, GCHandleType.Pinned);
                        try
                        {
                            IntPtr pointer = pinnedArray.AddrOfPinnedObject();
                            H5D.read(datasetId, typeId, H5S.ALL, H5S.ALL, H5P.DEFAULT, pointer);
                        }
                        finally
                        {
                            pinnedArray.Free();
                        }

                        if (buffer.GetLength(0) > 1) retbuffer = buffer; else continue;
                    }
                    finally
                    {
                        // Close resources
                        H5D.close(datasetId);
                    }
                    if (retbuffer != null) break;
                }

                return retbuffer;
            }

            static double[] ReadListenerView(long fileid)
            {
                // Open the dataset for the ListenerView field
                long datasetId = H5D.open(fileid, "ListenerView");
                if (datasetId < 0)
                {
                    Console.WriteLine("Failed to open dataset.");
                    H5F.close(fileid);
                    return null;
                }

                // Get the dataspace of the dataset
                long dataspaceId = H5D.get_space(datasetId);
                if (dataspaceId < 0)
                {
                    Console.WriteLine("Failed to get dataspace.");
                    H5D.close(datasetId);
                    H5F.close(fileid);
                    return null;
                }

                // Determine the size of the dataspace
                ulong[] dims = new ulong[3];
                H5S.get_simple_extent_dims(dataspaceId, dims, null);

                // Allocate array for reading data
                double[] listenerView = new double[dims[0]];

                // Read the data
                GCHandle hnd = GCHandle.Alloc(listenerView, GCHandleType.Pinned);
                H5D.read(datasetId, H5T.NATIVE_DOUBLE, H5S.ALL, H5S.ALL, H5P.DEFAULT, hnd.AddrOfPinnedObject());
                hnd.Free();

                // Output the ListenerView values
                Console.WriteLine("ListenerView:");
                foreach (var value in listenerView)
                {
                    Console.WriteLine(value);
                }

                // Cleanup
                H5S.close(dataspaceId);
                H5D.close(datasetId);
                H5F.close(fileid);
                return listenerView;
            }

            //private double[] ReadListinerView(long fileId)
            //{
            //    double[] retbuffer = null;
            //    // Open the dataset
            //        long datasetId = H5D.open(fileId, "ListenerView");

            //        if (datasetId < 0)
            //        {
            //            Console.WriteLine("Error opening EmitterPosition dataset.");
            //            return null;
            //        }

            //        try
            //        {
            //            // Get the datatype and dataspace of the dataset
            //            long typeId = H5D.get_type(datasetId);
            //            long spaceId = H5D.get_space(datasetId);

            //            // Determine the number of elements (positions) in the dataspace
            //            //int dims = H5S.(spaceId);

            //            // Allocate a buffer for the double array
            //            double[] buffer = new double[3];

            //            // Read the dataset
            //            GCHandle pinnedArray = GCHandle.Alloc(buffer, GCHandleType.Pinned);
            //            try
            //            {
            //                IntPtr pointer = pinnedArray.AddrOfPinnedObject();
            //                H5D.read(datasetId, typeId, H5S.ALL, H5S.ALL, H5P.DEFAULT, pointer);
            //            }
            //            finally
            //            {
            //                pinnedArray.Free();
            //            }
            //        }
            //        finally
            //        {
            //            // Close resources
            //            H5D.close(datasetId);
            //        }

            //    return retbuffer;
            //}

            private string ReadGlobalString(long fileId, string Field)
            {
                // Open the attribute
                long attrId = H5A.open_by_name(fileId, "/", Field, H5P.DEFAULT, H5P.DEFAULT);
                if (attrId < 0)
                {
                    H5F.close(fileId);
                    return null;
                }

                // Get the datatype of the attribute
                long typeId = H5A.get_type(attrId);
                if (typeId < 0)
                {
                    H5A.close(attrId);
                    H5F.close(fileId);
                    return null;
                }

                // Determine the size of the attribute
                IntPtr size = H5T.get_size(typeId);
                byte[] buffer = new byte[size.ToInt32()];

                // Read the attribute
                GCHandle pinnedArray = GCHandle.Alloc(buffer, GCHandleType.Pinned);
                IntPtr pointer = pinnedArray.AddrOfPinnedObject();
                if (H5A.read(attrId, typeId, pointer) < 0)
                {
                    Console.WriteLine("Error reading attribute.");
                    pinnedArray.Free();
                    H5T.close(typeId);
                    H5A.close(attrId);
                    H5F.close(fileId);
                    return null;
                }
                pinnedArray.Free();

                // Convert the byte array to a string
                string conventions = Encoding.ASCII.GetString(buffer).TrimEnd('\0');

                // Close resources
                H5T.close(typeId);
                H5A.close(attrId);

                return conventions;
            }

            //private double[] ReadListenerView(long fileId)
            //{
            //    // Open the attribute
            //    long attrId = H5A.open_by_name(fileId, "/", "ListenerView", H5P.DEFAULT, H5P.DEFAULT);
            //    if (attrId < 0)
            //    {
            //        H5F.close(fileId);
            //        return null;
            //    }

            //    // Get the datatype of the attribute
            //    long typeId = H5A.get_type(attrId);
            //    if (typeId < 0)
            //    {
            //        H5A.close(attrId);
            //        H5F.close(fileId);
            //        return null;
            //    }

            //    // Get the dataspace of the attribute
            //    long spaceId = H5A.get_space(attrId);
            //    if (spaceId < 0)
            //    {
            //        H5T.close(typeId);
            //        H5A.close(attrId);
            //        return null;
            //    }

            //    // Determine the number of elements in the dataspace
            //    ulong[] dims = new ulong[1];
            //    int rank = H5S.get_simple_extent_dims(spaceId, dims, null);
            //    if (rank < 0)
            //    {
            //        Console.WriteLine("Error getting dataspace dimensions.");
            //        H5S.close(spaceId);
            //        H5T.close(typeId);
            //        H5A.close(attrId);
            //        return null;
            //    }

            //    // Allocate a buffer for the double array
            //    double[] buffer = new double[dims[0]];

            //    // Read the attribute
            //    GCHandle pinnedArray = GCHandle.Alloc(buffer, GCHandleType.Pinned);
            //    IntPtr pointer = pinnedArray.AddrOfPinnedObject();
            //    if (H5A.read(attrId, typeId, pointer) < 0)
            //    {
            //        Console.WriteLine("Error reading attribute.");
            //        pinnedArray.Free();
            //        H5S.close(spaceId);
            //        H5T.close(typeId);
            //        H5A.close(attrId);
            //        return null;
            //    }
            //    pinnedArray.Free();

            //    // Close resources
            //    H5S.close(spaceId);
            //    H5T.close(typeId);
            //    H5A.close(attrId);

            //    return buffer;
            //}


            private float[] ReadSourcePositionUnits(long fileId)
            {
                // Open the attribute
                long attrId = H5A.open_by_name(fileId, "SourcePosition", "Units", H5P.DEFAULT, H5P.DEFAULT);
                if (attrId < 0)
                {
                    H5F.close(fileId);
                    return null;
                }

                // Get the datatype of the attribute
                long typeId = H5A.get_type(attrId);
                if (typeId < 0)
                {
                    H5A.close(attrId);
                    H5F.close(fileId);
                    return null;
                }

                // Get the dataspace of the attribute
                long spaceId = H5A.get_space(attrId);
                if (spaceId < 0)
                {
                    H5T.close(typeId);
                    H5A.close(attrId);
                    return null;
                }

                // Determine the number of elements in the dataspace
                ulong[] dims = new ulong[1];
                int rank = H5S.get_simple_extent_dims(spaceId, dims, null);
                if (rank < 0)
                {
                    Console.WriteLine("Error getting dataspace dimensions.");
                    H5S.close(spaceId);
                    H5T.close(typeId);
                    H5A.close(attrId);
                    return null;
                }

                // Allocate a buffer for the float array
                float[] buffer = new float[dims[0]];

                // Read the attribute
                GCHandle pinnedArray = GCHandle.Alloc(buffer, GCHandleType.Pinned);
                IntPtr pointer = pinnedArray.AddrOfPinnedObject();
                if (H5A.read(attrId, typeId, pointer) < 0)
                {
                    Console.WriteLine("Error reading attribute.");
                    pinnedArray.Free();
                    H5S.close(spaceId);
                    H5T.close(typeId);
                    H5A.close(attrId);
                    return null;
                }
                pinnedArray.Free();

                // Close resources
                H5S.close(spaceId);
                H5T.close(typeId);
                H5A.close(attrId);

                return buffer;
            }

            private void ReadHrtfDataset(long fileId, string datasetName)
            {
                long datasetId = H5D.open(fileId, datasetName);
                if (datasetId < 0)
                {
                    Console.WriteLine($"Error opening dataset {datasetName}.");
                    return;
                }

                // Allocate a buffer for the data
                double[,,] dataIR;

                try
                {
                    // Assuming the dataset is 3D: Measurements x Receivers x Samples
                    long dataspaceId = H5D.get_space(datasetId);
                    ulong[] dims = new ulong[3]; // Adjust based on actual dataset dimensions
                    float[] spherical = ReadSourcePositionUnits(fileId);

                    H5S.get_simple_extent_dims(dataspaceId, dims, null);

                    dataIR = new double[dims[0], dims[1], dims[2]];
                    GCHandle handle = GCHandle.Alloc(dataIR, GCHandleType.Pinned);
                    try
                    {
                        H5D.read(datasetId, H5T.NATIVE_FLOAT, H5S.ALL, H5S.ALL, H5P.DEFAULT, handle.AddrOfPinnedObject());
                    }
                    finally
                    {
                        handle.Free();
                    }

                    // Process the HRTF data as needed for integration with Pachyderm
                    // This might involve converting the data into a format Pachyderm can use for simulations
                }
                finally
                {
                    H5D.close(datasetId);
                }

                double NormL = 0;
                double NormR = 0;

                int idx = dataIR.GetLength(0);
                int channel = dataIR.GetLength(1);
                int samples = dataIR.GetLength(2);

                HRIR = new double[idx][][];
                for (int i = 0; i < idx; i++)
                {
                    HRIR[i] = new double[2][];
                    HRIR[i][0] = new double[dataIR.GetLength(2)];
                    HRIR[i][1] = new double[dataIR.GetLength(2)];
                    for (int k = 0; k < samples; k++)
                    {
                        NormL = Math.Max(NormL, Math.Abs(dataIR[i, 0, k]));
                        NormR = Math.Max(NormL, Math.Abs(dataIR[i, 1, k]));
                        HRIR[i][0][k] = dataIR[i, 0, k];
                        HRIR[i][1][k] = dataIR[i, 1, k];
                    }
                }

                for (int i = 0; i < idx; i++)
                {
                    for (int j = 0; j < channel; j++)
                    {
                        for (int k = 0; k < samples; k++)
                        {
                            HRIR[i][0][k] /= NormL;
                            HRIR[i][1][k] /= NormR;
                        }
                    }
                }
            }

            /// <summary>
            /// Finds the appropriate HRTF for the incoming direction of sound...
            /// </summary>
            /// <param name="Dir_incoming">The direction sound is coming from, considering 0,0,0 as origin, and the direction the listener perceives it from (the opposite of the direction sound is traveling.).</param>
            /// <param name="normalized">Indicate if the direction is normalized. A value of false will cause the vector to be normalized in this method.</param>
            /// <param name="samplingfrequency">The required sampling frequency of the output signal HRTF.../param>
            /// <returns></returns>
            public double[][] Lookup(Vector Dir_incoming, bool normalized, int samplingfrequency)
            {
                if (!normalized) Dir_incoming.Normalize();
                X_Event x;
                Ray r = new Ray(origin, Dir_incoming, 0, rand.Next());
                if (!VG.Shoot(r, 0, out x)) throw new Exception("Invalid directional input to HRTF lookup...");
                int[] pos = bins[x.Poly_id].ToArray();

                double min = double.MaxValue;
                int idx = 0;
                foreach (int p in pos)
                {
                    Vector dl = Dir_incoming - emitterVectors[p];
                    double l2 = dl.dx * dl.dx + dl.dy * dl.dy + dl.dz * dl.dz;
                    if (l2 < min)
                    {
                        min = l2;
                        idx = p;
                    }
                }

                double[][] hrir_out = new double[HRIR[idx].Length][];
                for (int i = 0; i < HRIR[idx].Length; i++)
                {
                    hrir_out[i] = new double[HRIR[idx][i].Length];
                    for (int j = 0; j < HRIR[idx][i].Length; j++) hrir_out[i][j] = HRIR[idx][i][j];
                    if (samplingfrequency != Fs) hrir_out[i] = Pachyderm_Acoustic.Audio.Pach_SP.Resample_Cubic(hrir_out[i], Fs, samplingfrequency, 0, true);
                }

                return hrir_out;
            }

            public double[][] Rose_To_Sphere(double[][] Vector_Rose)
            {
                double PI_180 = Math.PI / 180;
                double[][] sphere = new double[Directions.Length][];
                double power = Directions.Length / 2 - 1;

                int[] ids = new int[3];

                for (int i = 0; i < Directions.Length; i++)
                {
                    sphere[i] = new double[Vector_Rose[0].Length];

                    if (emitterdirs[i, 0] <= 90 || emitterdirs[i, 0] >= 270) { ids[0] = 0; }
                    else { ids[0] = 1; }
                    if (emitterdirs[i, 0] <= 180) { ids[1] = 2; }
                    else { ids[1] = 3; }
                    if (emitterdirs[i, 1] >= 0) { ids[2] = 4; }
                    else { ids[2] = 5; }

                    double rosemodx = Math.Pow(Math.Cos(PI_180 * emitterdirs[i, 0]), power);
                    double rosemody = Math.Pow(Math.Sin(PI_180 * emitterdirs[i, 1]), power);
                    double rosemodz = Math.Pow(Math.Sin(PI_180 * emitterdirs[i, 2]), power);

                    for (int j = 0; j < Vector_Rose[0].Length; j++)
                    {
                        sphere[i][j] = Vector_Rose[ids[0]][j] * rosemodx + Vector_Rose[ids[1]][j] * rosemody + Vector_Rose[ids[2]][j] * rosemodz;
                    }
                }

                return sphere;
            }

            public double[][] Filter(double[][] Vector_Rose)
            {
                double[][] Sphere = Rose_To_Sphere(Vector_Rose);
                double[][] signal = new double[2][];
                signal[0] = new double[Vector_Rose[0].Length + HRIR[0][0].Length];
                signal[1] = new double[Vector_Rose[0].Length + HRIR[0][0].Length];

                for (int i = 0; i < Sphere.Length; i++)
                {
                    double[] filtL = Pachyderm_Acoustic.Audio.Pach_SP.FFT_Convolution_double(Sphere[i], HRIR[Face_Centroid_ID[i]][0], 0);
                    double[] filtR = Pachyderm_Acoustic.Audio.Pach_SP.FFT_Convolution_double(Sphere[i], HRIR[Face_Centroid_ID[i]][1], 0);
                    for (int t = 0; t < filtL.Length; t++) { signal[0][t] += filtL[t]; signal[1][t] += filtR[t]; }
                }
                return signal;
            }

            public double[] Azimuth
            {
                get
                {
                    return Azi;
                }
            }
            public double[] Altitude
            {
                get
                {
                    return Alt;
                }
            }

            double[][] Loaded_Filter;
            double[][][] Loaded_HRIR;
            double[][] Translation;
            public void Load(Direct_Sound Direct, ImageSourceData ISData, Pachyderm_Acoustic.Environment.Receiver_Bank RTData, double CO_Time_ms, int Sampling_Frequency, int Rec_ID, bool Start_at_Zero, bool flat)
            {
                Loaded_Filter = new double[6][];
                Loaded_Filter[0] = Pachyderm_Acoustic.Utilities.IR_Construction.Aurfilter_Directional(Direct, ISData, RTData, CO_Time_ms, Sampling_Frequency, Rec_ID, Start_at_Zero, 0, 0, true, flat);
                Loaded_Filter[1] = Pachyderm_Acoustic.Utilities.IR_Construction.Aurfilter_Directional(Direct, ISData, RTData, CO_Time_ms, Sampling_Frequency, Rec_ID, Start_at_Zero, 0, 180, true, flat);
                Loaded_Filter[2] = Pachyderm_Acoustic.Utilities.IR_Construction.Aurfilter_Directional(Direct, ISData, RTData, CO_Time_ms, Sampling_Frequency, Rec_ID, Start_at_Zero, 0, 90, true, flat);
                Loaded_Filter[3] = Pachyderm_Acoustic.Utilities.IR_Construction.Aurfilter_Directional(Direct, ISData, RTData, CO_Time_ms, Sampling_Frequency, Rec_ID, Start_at_Zero, 0, 270, true, flat);
                Loaded_Filter[4] = Pachyderm_Acoustic.Utilities.IR_Construction.Aurfilter_Directional(Direct, ISData, RTData, CO_Time_ms, Sampling_Frequency, Rec_ID, Start_at_Zero, 90, 0, true, flat);
                Loaded_Filter[5] = Pachyderm_Acoustic.Utilities.IR_Construction.Aurfilter_Directional(Direct, ISData, RTData, CO_Time_ms, Sampling_Frequency, Rec_ID, Start_at_Zero, -90, 0, true, flat);

                double power = (Directions.Length / 2 - 1) / 2;

                Loaded_HRIR = new double[DirsCT][][];
                Translation = new double[DirsCT][];
                for (int i = 0; i < DirsCT; i++)
                {
                    if (Face_Centroid_ID[i] == -1) continue;
                    Loaded_HRIR[i] = new double[2][];
                    for (int j = 0; j < 2; j++)
                    {
                        Loaded_HRIR[i][j] = Pach_SP.Resample(HRIR[Face_Centroid_ID[i]][j], 48000, Sampling_Frequency, 0, true);
                    }
                    Translation[i] = new double[3];
                    Translation[i][0] = Math.Pow(Math.Abs(Directions[i].dx), power);
                    Translation[i][1] = Math.Pow(Math.Abs(Directions[i].dy), power);
                    Translation[i][2] = Math.Pow(Math.Abs(Directions[i].dz), power);
                }
            }

            public double[][] Binaural_IR(double _azi, double _alt)
            {
                double[][] filt = PachTools.Rotate_Vector_Rose(Loaded_Filter, -_azi, -_alt, true);
                double[][] Signal = new double[2][];
                Signal[0] = new double[filt.Length + 2048 + HRIR[0][0].Length];
                Signal[1] = new double[Signal[0].Length];

                for (int i = 0; i < DirsCT; i++)
                {
                    if (Face_Centroid_ID[i] == -1) continue;
                    double[] s = new double[filt.Length];
                    if (Directions[i].dx > 0) for (int t = 0; t < filt.Length; t++) s[t] += filt[t][0] * Translation[i][0];
                    else for (int t = 0; t < 0; t++) s[t] = filt[t][1] * Translation[i][0];
                    if (Directions[i].dy > 0) for (int t = 0; t < filt.Length; t++) s[t] += filt[t][2] * Translation[i][1];
                    else for (int t = 0; t < 0; t++) s[t] = filt[t][3] * Translation[i][1];
                    if (Directions[i].dz > 0) for (int t = 0; t < filt.Length; t++) s[t] += filt[t][4] * Translation[i][2];
                    else for (int t = 0; t < 0; t++) s[t] = filt[t][5] * Translation[i][2];

                    double[][] hs = new double[2][];
                    hs[0] = Pach_SP.FFT_Convolution_double(s, Loaded_HRIR[i][0], 0);
                    hs[1] = Pach_SP.FFT_Convolution_double(s, Loaded_HRIR[i][1], 0);

                    for (int t = 0; t < Signal[0].Length; t++)
                    {
                        Signal[0][t] += hs[0][t];
                        Signal[1][t] += hs[1][t];
                    }
                }

                return Signal;
            }
        }
    }
}