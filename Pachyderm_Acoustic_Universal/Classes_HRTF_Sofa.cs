using HDF.PInvoke;
using System;
using System.Linq;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Text;
using Hare.Geometry;
using Vector = Hare.Geometry.Vector;
using Pachyderm_Acoustic.Utilities;
using Pachyderm_Acoustic.Audio;

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
            float[][][] HRIR;
            int[] Face_Centroid_ID;
            Point origin;
            Random rand = new Random();
            private int Fs;
            Vector[] PrincipalDirections;
            double[] Azi, Alt;

            public bool AngularDistancePass { get; private set; }
            public bool GlobeCoveragePass { get; private set; }
            public string ValidationMessage { get; private set; }
            public bool ValidationPassed => AngularDistancePass && GlobeCoveragePass;

            public double AvgAngularDistanceFront { get; private set; }
            public double AvgAngularDistanceBack { get; private set; }
            public double MaxCoverageGap { get; private set; }

            double[][][] Loaded_HRIR;
            double[][][] Loaded_HRIR_44100;
            double[][] Translation;

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
                    emitterdirs = ReadSourcePosition(fileId);
                    Fs = ReadSamplingFrequency(fileId);

                    AngularDistancePass = CheckHrtfAngularDistance();
                    GlobeCoveragePass = CheckHrtfGlobeCoverage(15.0);

                    AvgAngularDistanceFront = GetAvgAngularDistanceFront();
                    AvgAngularDistanceBack = GetAvgAngularDistanceBack();
                    MaxCoverageGap = GetMaxCoverageGap();

                    if (!AngularDistancePass && !GlobeCoveragePass)
                    {
                        ValidationMessage = "HRTF failed both angular distance and globe coverage checks:\n" +
                                            $"Average angular distance in front hemisphere: {AvgAngularDistanceFront:F2}° (tolerance: 8°)\n" +
                                            $"Average angular distance in back hemisphere: {AvgAngularDistanceBack:F2}° (tolerance: 8°)\n" +
                                            $"Maximum coverage gap found: {MaxCoverageGap:F2}° (tolerance: 15°)\n" +
                                            "For more detailed diagnostics and a visualisation of the sampling density of your .SOFA file, please see: \n" +
                                            "https://github.com/domfrbassett/HRTF_Diagnostic_Tool";
                    }
                    else if (!AngularDistancePass)
                    {
                        ValidationMessage = "HRTF failed angular distance check.\n" +
                                            $"Average angular distance in front hemisphere: {AvgAngularDistanceFront:F2}° (tolerance: 8°)\n" +
                                            $"Average angular distance in back hemisphere: {AvgAngularDistanceBack:F2}° (tolerance: 8°)\n" +
                                            "For more detailed diagnostics and a visualisation of the sampling density of your .SOFA file, please see: \n" +
                                            "https://github.com/domfrbassett/HRTF_Diagnostic_Tool";
                    }
                    else if (!GlobeCoveragePass)
                    {
                        ValidationMessage = "HRTF failed globe coverage check.\n" +
                                            $"Maximum coverage gap found: {MaxCoverageGap:F2}° (tolerance: 15°)\n" +
                                            "For more detailed diagnostics and a visualisation of the sampling density of your .SOFA file, please see: \n" +
                                            "https://github.com/domfrbassett/HRTF_Diagnostic_Tool";
                    }

                    int order = GetRequiredSubdivisionOrder(emitterdirs.GetLength(0));

                    //Create a reference sphere...
                    VG = Pachyderm_Acoustic.Utilities.Geometry.GeoSphere(order);
                    T = VG.Model[0];
                    PrincipalDirections = new Vector[T.Polygon_Count];
                    for (int i = 0; i < T.Polygon_Count; i++)
                    {
                        Vector d = new Vector(T.Polys[i].Centroid.x, T.Polys[i].Centroid.y, T.Polys[i].Centroid.z);
                        d.Normalize();
                        PrincipalDirections[i] = d;
                    }

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
                        double azRad = Math.PI * emitterdirs[i, 0] / 180.0;
                        double elRad = Math.PI * emitterdirs[i, 1] / 180.0;
                        double x = Math.Cos(elRad) * Math.Cos(azRad);
                        double y = Math.Cos(elRad) * Math.Sin(azRad);
                        double z = Math.Sin(elRad);
                        Vector d = new Vector(x, y, z);
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

                    // Read HRIR data
                    ReadHrtfDataset(fileId, "Data.IR");

                    var validHRIRs = new List<double[][]>();
                    var validTranslations = new List<double[]>();
                    var validDirections = new List<Vector>();

                    int validFaces = Face_Centroid_ID.Count(id => id != -1);
                    double power = (validFaces / 2.0 - 1) / 2.0;

                    for (int i = 0; i < Directions.Length; i++)
                    {
                        if (Face_Centroid_ID[i] == -1) continue;
                        var hrir = Pach_SP_HRTF.ResampleHRIRWDL(HRIR[Face_Centroid_ID[i]], Fs, 44100);

                        var trans = new double[3];
                        trans[0] = Math.Pow(Math.Abs(Directions[i].dx), power);
                        trans[1] = Math.Pow(Math.Abs(Directions[i].dy), power);
                        trans[2] = Math.Pow(Math.Abs(Directions[i].dz), power);

                        validHRIRs.Add(hrir);
                        validTranslations.Add(trans);
                        validDirections.Add(Directions[i]);
                    }

                    Loaded_HRIR_44100 = validHRIRs.ToArray(); // Store HRIRs resampled to 44100 Hz - useful for many applications and can be retrieved quickly
                    Translation = validTranslations.ToArray();
                    Directions = validDirections.ToArray();
                }
                finally
                {
                    H5F.close(fileId);
                }
            }

            /// <summary>
            /// This method determines the required subdivision order based on the total number of sources. This assumes a fully covered sphere, which is rarely the case with real HRTF datasets. 
            /// However, it provides a reasonable estimate for the subdivision order needed to achieve a certain number of directions. In the future, this could be improved by analysing the actual distribution of sources in the dataset.
            /// 
            /// </summary>
            /// <param name="totalSources"></param>
            /// <returns></returns>
            public static int GetRequiredSubdivisionOrder(int totalSources)
            {
                if (totalSources <= 20) return 0;    // Order 0: 20 faces
                if (totalSources <= 80) return 1;    // Order 1: 80 faces
                if (totalSources <= 320) return 2;   // Order 2: 320 faces
                if (totalSources <= 1280) return 3;  // Order 3: 1280 faces
                return 4;                            // Order 4: 5120 faces
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
                set { PrincipalDirections = value; }
            }

            private bool CheckHrtfAngularDistance()
            {
                int nSources = emitterdirs.GetLength(0);

                // Convert spherical to cartesian unit vectors
                var vectors = new Vector[nSources];
                for (int i = 0; i < nSources; i++)
                {
                    vectors[i] = PachTools.SphericalToCartesian(emitterdirs[i, 0], emitterdirs[i, 1]);
                }

                // Separate front/back using X coordinate (positive = front, negative = back)
                var frontVectors = vectors.Where(v => v.dx >= 0).ToArray();
                var backVectors = vectors.Where(v => v.dx < 0).ToArray();

                // Fail if either hemisphere has no sources
                if (frontVectors.Length == 0 || backVectors.Length == 0)
                    return false;

                // Compute average minimum angular distance
                AvgAngularDistanceFront = AverageMinAngularDistance(frontVectors);
                AvgAngularDistanceBack = AverageMinAngularDistance(backVectors);

                // Check against threshold (8°)
                return AvgAngularDistanceFront <= 8 && AvgAngularDistanceBack <= 8;
            }

            private bool CheckHrtfGlobeCoverage(double maxAllowedGapDeg)
            {
                var emitters = Enumerable.Range(0, emitterdirs.GetLength(0))
                    .Select(i => PachTools.SphericalToCartesian(emitterdirs[i, 0], emitterdirs[i, 1]))
                    .ToArray();

                int N_az = 73; // Number of azimuthal divisions - approximately 5° spacing
                int N_el = 19; // Number of elevation divisions - approximately 5° spacing

                var refPoints = new List<Hare.Geometry.Vector>();

                for (int elIdx = 0; elIdx < N_el; elIdx++)
                {
                    double el = 90.0 * elIdx / (N_el - 1);
                    for (int azIdx = 0; azIdx < N_az; azIdx++)
                    {
                        double az = 360.0 * azIdx / N_az;
                        refPoints.Add(PachTools.SphericalToCartesian(az, el));
                    }
                }

                MaxCoverageGap = 0.0;
                for (int i = 0; i < refPoints.Count; i++)
                {
                    double minAngle = emitters.Min(emitter => AngularDistanceDegrees(refPoints[i], emitter));
                    if (minAngle > MaxCoverageGap) MaxCoverageGap = minAngle;

                    if (minAngle > maxAllowedGapDeg)
                    {
                        return false;
                    }
                }

                return true;
            }

            public static double AngularDistanceDegrees(Hare.Geometry.Vector v1, Hare.Geometry.Vector v2)
            {
                double dot = Hare.Geometry.Hare_math.Dot(v1, v2);
                dot = Math.Min(1.0, Math.Max(-1.0, dot)); // Clamp to [-1, 1] to avoid NaN
                return Math.Acos(dot) * 180.0 / Math.PI; // Convert radians to degrees
            }

            private double AverageMinAngularDistance(Hare.Geometry.Vector[] vectors)
            {
                int n = vectors.Length;
                double[] minAngles = new double[n];

                for (int i = 0; i < n; i++)
                {
                    double minAngle = 180.0;
                    for (int j = 0; j < n; j++)
                    {
                        if (i == j) continue;
                        double angle = AngularDistanceDegrees(vectors[i], vectors[j]);
                        if (angle < minAngle) minAngle = angle;
                    }
                    minAngles[i] = minAngle;
                }

                return minAngles.Average();
            }

            public double GetAvgAngularDistanceFront()
            {
                return AvgAngularDistanceFront;
            }

            public double GetAvgAngularDistanceBack()
            {
                return AvgAngularDistanceBack;
            }

            public double GetMaxCoverageGap()
            {
                return MaxCoverageGap;
            }

            public int ReadSamplingFrequency(long fileId)
            {
                long datasetId = H5D.open(fileId, "/Data.SamplingRate");
                if (datasetId < 0)
                {
                    throw new Exception("Failed to open dataset '/Data.SamplingRate'. This dataset is required to get the sampling frequency.");
                }

                try
                {
                    long typeId = H5D.get_type(datasetId);
                    if (typeId < 0)
                    {
                        throw new Exception("Failed to get datatype for '/Data.SamplingRate'.");
                    }

                    try
                    {
                        double[] Fs = new double[1];
                        GCHandle hnd = GCHandle.Alloc(Fs, GCHandleType.Pinned);

                        try
                        {
                            IntPtr ptr = hnd.AddrOfPinnedObject();
                            if (H5D.read(datasetId, typeId, H5S.ALL, H5S.ALL, H5P.DEFAULT, ptr) < 0)
                            {
                                throw new Exception("Failed to read '/Data.SamplingRate'.");
                            }

                            System.Diagnostics.Debug.WriteLine($"Input sampling rate: {(int)Math.Round(Fs[0])} Hz");
                            System.Diagnostics.Debug.WriteLine("Output sampling rate (target): 44100 Hz");
                            return (int)Math.Round(Fs[0]);
                        }
                        finally
                        {
                            hnd.Free();
                        }
                    }
                    finally
                    {
                        H5T.close(typeId);
                    }
                }
                finally
                {
                    H5D.close(datasetId);
                }
            }

            private double[,] ReadSourcePosition(long fileId)
            {
                string[] candidates = { "SourcePosition", "EmitterPosition" };

                foreach (var name in candidates)
                {
                    long datasetId = H5D.open(fileId, name);
                    if (datasetId < 0)
                    {
                        Console.WriteLine($"Dataset {name} not found.");
                        continue;
                    }

                    long spaceId = -1;

                    try
                    {
                        spaceId = H5D.get_space(datasetId);
                        if (spaceId < 0)
                        {
                            Console.WriteLine($"Error getting dataspace for {name}.");
                            continue;
                        }

                        ulong[] dims = new ulong[2];
                        if (H5S.get_simple_extent_dims(spaceId, dims, null) < 0)
                        {
                            Console.WriteLine($"Error getting dataset dimensions for {name}.");
                            continue;
                        }

                        // Accept either N x 3 or 3 x N
                        bool transposed = false;
                        int numSources;
                        if (dims[1] == 3)
                        {
                            numSources = (int)dims[0]; // N x 3
                        }
                        else if (dims[0] == 3)
                        {
                            transposed = true; // 3 x N
                            numSources = (int)dims[1];
                        }
                        else
                        {
                            Console.WriteLine($"Unexpected {name} format. Expected Nx3 or 3xN, got {dims[0]}x{dims[1]}.");
                            return null;
                        }

                        // Flat buffer
                        int len = numSources * 3;
                        double[] flat = new double[len];

                        // Pin the flat buffer
                        GCHandle h = GCHandle.Alloc(flat, GCHandleType.Pinned);
                        try
                        {
                            IntPtr ptr = h.AddrOfPinnedObject();

                            if (H5D.read(datasetId, H5T.NATIVE_DOUBLE, H5S.ALL, H5S.ALL, H5P.DEFAULT, ptr) < 0)
                            {
                                Console.WriteLine($"Failed to read {name} dataset.");
                                continue;
                            }
                        }
                        finally
                        {
                            h.Free();
                        }

                        // Build final result [numSources, 2] (azimuth, elevation)
                        double[,] result = new double[numSources, 2];
                        if (!transposed)
                        {
                            for (int i = 0; i < numSources; ++i)
                            {
                                result[i, 0] = flat[i * 3 + 0]; // azimuth
                                result[i, 1] = flat[i * 3 + 1]; // elevation
                            }
                        }
                        else
                        {
                            for (int i = 0; i < numSources; ++i)
                            {
                                result[i, 0] = flat[0 * numSources + i]; // azimuth row
                                result[i, 1] = flat[1 * numSources + i]; // elevation row
                            }
                        }

                        Console.WriteLine($"Read {numSources} source positions from {name}.");
                        return result;
                    }
                    finally
                    {
                        if (spaceId >= 0) H5S.close(spaceId);
                        if (datasetId >= 0) H5D.close(datasetId);
                    }
                }

                Console.WriteLine("No valid SourcePosition or EmitterPosition dataset found.");
                return null;
            }

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

                try
                {
                    long dataspaceId = H5D.get_space(datasetId);
                    ulong[] dims = new ulong[3];
                    H5S.get_simple_extent_dims(dataspaceId, dims, null);
                    int measurements = (int)dims[0];
                    int channels = (int)dims[1];
                    int samples = (int)dims[2];

                    float[,,] dataIR = new float[measurements, channels, samples];
                    GCHandle handle = GCHandle.Alloc(dataIR, GCHandleType.Pinned);
                    try
                    {
                        H5D.read(datasetId, H5T.NATIVE_FLOAT, H5S.ALL, H5S.ALL, H5P.DEFAULT, handle.AddrOfPinnedObject());
                    }
                    finally
                    {
                        handle.Free();
                    }

                    HRIR = new float[measurements][][];
                    for (int i = 0; i < measurements; i++)
                    {
                        HRIR[i] = new float[channels][];
                        for (int ch = 0; ch < channels; ch++)
                        {
                            HRIR[i][ch] = new float[samples];
                            for (int k = 0; k < samples; k++)
                                HRIR[i][ch][k] = dataIR[i, ch, k];
                        }
                    }
                }
                finally
                {
                    H5D.close(datasetId);
                }
            }

            double[][] Loaded_Filter;

            public void Load(Direct_Sound Direct, ImageSourceData ISData, Pachyderm_Acoustic.Environment.Receiver_Bank RTData, SystemResponseCompensation.SystemCompensationSettings sysCompSettings, double CO_Time_ms, int targetFs, int Rec_ID, bool Start_at_Zero, bool flat, bool auto)
            {
                Loaded_Filter = new double[6][];
                Loaded_Filter[0] = Pachyderm_Acoustic.Utilities.IR_Construction.Aurfilter_Directional(Direct, ISData, RTData, CO_Time_ms, targetFs, Rec_ID, Start_at_Zero, 0, 0, true, flat);
                Loaded_Filter[1] = Pachyderm_Acoustic.Utilities.IR_Construction.Aurfilter_Directional(Direct, ISData, RTData, CO_Time_ms, targetFs, Rec_ID, Start_at_Zero, 0, 180, true, flat);
                Loaded_Filter[2] = Pachyderm_Acoustic.Utilities.IR_Construction.Aurfilter_Directional(Direct, ISData, RTData, CO_Time_ms, targetFs, Rec_ID, Start_at_Zero, 0, 90, true, flat);
                Loaded_Filter[3] = Pachyderm_Acoustic.Utilities.IR_Construction.Aurfilter_Directional(Direct, ISData, RTData, CO_Time_ms, targetFs, Rec_ID, Start_at_Zero, 0, 270, true, flat);
                Loaded_Filter[4] = Pachyderm_Acoustic.Utilities.IR_Construction.Aurfilter_Directional(Direct, ISData, RTData, CO_Time_ms, targetFs, Rec_ID, Start_at_Zero, 90, 0, true, flat);
                Loaded_Filter[5] = Pachyderm_Acoustic.Utilities.IR_Construction.Aurfilter_Directional(Direct, ISData, RTData, CO_Time_ms, targetFs, Rec_ID, Start_at_Zero, -90, 0, true, flat);

                bool _sysCompApplied = false;

                if (targetFs != Fs)
                {
                    if (targetFs == 44100)
                    {
                        var validHRIRs = new List<double[][]>();

                        for (int i = 0; i < Loaded_HRIR_44100.Length; i++)
                        {
                            var hrirCopy = new double[Loaded_HRIR_44100[i].Length][];
                            for (int ch = 0; ch < Loaded_HRIR_44100[i].Length; ch++)
                            {
                                hrirCopy[ch] = new double[Loaded_HRIR_44100[i][ch].Length];
                                Array.Copy(Loaded_HRIR_44100[i][ch], hrirCopy[ch], Loaded_HRIR_44100[i][ch].Length);
                            }
                            validHRIRs.Add(hrirCopy);
                        }

                        Loaded_HRIR = validHRIRs.ToArray();
                    }
                    else
                    {
                        var validHRIRs = new List<double[][]>();
                        for (int i = 0; i < DirsCT; i++)
                        {
                            int faceId = Face_Centroid_ID[i];
                            if (faceId == -1) continue;
                            var hrirResampled = new double[2][];

                            hrirResampled = Pach_SP_HRTF.ResampleHRIRWDL(
                                HRIR[faceId],
                                Fs,
                                targetFs
                            );

                            validHRIRs.Add(hrirResampled);
                        }
                        Loaded_HRIR = validHRIRs.ToArray();
                    }
                }
                else
                {
                    var validHRIRs = new List<double[][]>();

                    for (int i = 0; i < DirsCT; i++)
                    {
                        int faceId = Face_Centroid_ID[i];
                        if (faceId == -1) continue;

                        var hrirCopy = new double[HRIR[faceId].Length][];

                        for (int ch = 0; ch < HRIR[faceId].Length; ch++)
                        {
                            hrirCopy[ch] = new double[HRIR[faceId][ch].Length];
                            for (int k = 0; k < HRIR[faceId][ch].Length; k++)
                            {
                                hrirCopy[ch][k] = HRIR[faceId][ch][k];
                            }
                        }

                        validHRIRs.Add(hrirCopy);
                    }

                    Loaded_HRIR = validHRIRs.ToArray();
                }

                if (!sysCompSettings.InputIsDTF)
                {
                    Pach_SP_HRTF.ApplySystemCompensation(Loaded_HRIR, Directions, Fs, 0, sysCompSettings, auto);
                }

                Pach_SP_HRTF.ShiftHRIRPairs(Loaded_HRIR, 0.8);
            }

            public double[][] Binaural_IR(double _azi, double _alt)
            {
                // Rotate the directional filters according to the given azimuth and elevation
                double[][] rotatedDirectionalFilters = PachTools.Rotate_Vector_Rose(Loaded_Filter, -_azi, -_alt, true);

                List<double[]> directionalSignalsList = new List<double[]>();

                // Build directional component signals from the rotated filters
                for (int i = 0; i < Directions.Length; i++)
                {
                    double[] directionalSignal = Pach_SP_HRTF.BuildDirectionalComponent(i, rotatedDirectionalFilters, Directions, Translation);
                    directionalSignalsList.Add(directionalSignal);
                }

                double[][] directionalSignals = directionalSignalsList.ToArray();
                double[] dryDirectionalSignal = Pach_SP_HRTF.BuildDrySignal(directionalSignals);
                double dryRMS = Pach_SP_HRTF.ComputeRMS(dryDirectionalSignal);

                // Convolve each directional signal with its corresponding HRIR
                double[][][] convolvedSignals = new double[directionalSignals.Length][][];
                for (int idx = 0; idx < directionalSignals.Length; idx++)
                {
                    convolvedSignals[idx] = new double[2][];
                    convolvedSignals[idx][0] = Pach_SP.FFT_Convolution_double(directionalSignals[idx], Loaded_HRIR[idx][0], 0);
                    convolvedSignals[idx][1] = Pach_SP.FFT_Convolution_double(directionalSignals[idx], Loaded_HRIR[idx][1], 0);
                }

                // Sum across directions to get the final binaural signal
                double[][] Signal = Pach_SP_HRTF.SumAcrossDirections(convolvedSignals);
                Pach_SP_HRTF.NormaliseStereoByDryRMS(Signal, dryRMS);

                return Signal;
            }
        }
    }
}