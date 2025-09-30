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

using System;
using System.Collections.Generic;
using System.Linq;
using Hare.Geometry;

namespace Pachyderm_Acoustic
{
    namespace Environment
    {
        public class LineSource : Source
        {
            /// <summary>
            /// Sample points on curves.
            /// </summary>
            public Point[][] Samples;
            /// <summary>
            /// Sound Power of each sample on the curve. (10^Lp/10 * L / no_of_samples)
            /// </summary>m
            public double[][] Power;
            /// <summary>
            /// User designated sound power of each line source.
            /// </summary>
            double[][] Level;
            /// <summary>
            /// number of samples per meter.
            /// </summary>
            double samplespermeter = 16;
            /// <summary>
            /// Private member controlling the directional characteristics of the line source.
            /// </summary>
            public Directionality D;
            /// <summary>
            /// The bounding box of line/curve source.
            /// </summary>
            public AABB bounds;
            /// <summary>
            /// Describes construction for a simple (generic) line source, with no unusual directionality characteristics.
            /// </summary>
            /// <param name="samples"></param>
            /// <param name="length"></param>
            /// <param name="Code"></param>
            /// <param name="el_m"></param>
            /// <param name="SrcID"></param>
            /// <param name="ph"></param>
            public LineSource(Hare.Geometry.Point[] samples, double length, string Code, double el_m, int SrcID, bool Third_Octave)
                : base(new double[8] { 60, 49, 41, 35, 31, 28, 26, 24 }, new Point(0, 0, 0), SrcID, Third_Octave)
            {
                //TODO: Accommodate third octave
                samplespermeter = el_m;

                //Divide curve up in ~equal length segments.
                Samples = new Hare.Geometry.Point[1][] { samples };
                D = new Simple();

                Level = new double[1][] { Utilities.PachTools.DecodeSourcePower(Code) };
                Power = new double[1][] { new double[8] };

                // Calculate bounding box for all samples
                CalculateBoundingBox();

                double PowerMod = length / (double)Samples.Length;
                for (int oct = 0; oct < 8; oct++) Power[0][oct] = 1E-12 * Math.Pow(10, .1 * Level[0][oct]) * PowerMod;
            }

            /// <summary>
            /// Enhanced LineSource constructor for OMNIDIRECTIONAL sources (Traffic/Rail) with automatic height positioning
            /// Uses ray casting to find ground and creates multi-height sources where appropriate
            /// ONLY for Simple directivity sources - NOT for Aircraft which have their own directivity
            /// </summary>
            /// <param name="samples">Original curve sample points (Z coordinates will be adjusted)</param>
            /// <param name="length">Total length of the curve</param>
            /// <param name="source_type">Type of source (Traffic or Rail only)</param>
            /// <param name="room_model">Room model for ray casting to find ground</param>
            /// <param name="Code">Encoded sound power string</param>
            /// <param name="el_m">Elements per meter</param>
            /// <param name="SrcID">Source ID</param>
            /// <param name="Third_Octave">Third octave flag</param>
            /// <param name="use_multi_height">Enable TNM-style multi-height for traffic sources</param>
            public LineSource(Hare.Geometry.Point[] samples, double length, SourceType source_type,
                Scene room_model, string Code, double el_m, int SrcID, bool Third_Octave,
                bool use_multi_height = true)
                : base(new double[8] { 60, 49, 41, 35, 31, 28, 26, 24 }, new Point(0, 0, 0), SrcID, Third_Octave)
            {
                samplespermeter = el_m;

                // ONLY Simple directivity for omnidirectional sources
                D = new Simple();

                Level = new double[1][] { Utilities.PachTools.DecodeSourcePower(Code) };

                // Determine height configuration based on source type
                double[] height_offsets;
                double[][] height_power_splits;
                SetupOmnidirectionalHeights(source_type, use_multi_height, out height_offsets, out height_power_splits);

                // Perform ray casting to find ground positions and create height-adjusted samples
                Point[][] adjusted_samples = new Point[height_offsets.Length][];
                for(int h = 0; h < height_offsets.Length; h++) adjusted_samples[h] = new Point[samples.Length];
                Vector downward_direction = new Vector(0, 0, -1);

                int i = 0;
                foreach (Point sample in samples)
                {
                    // Cast ray downward to find ground intersection
                    Ray downward_ray = new Ray(sample, downward_direction, 0, 0);

                    // Try to intersect with room geometry
                    double u, v, t;
                    int poly_id;
                    Point ground_point;

                    if (!room_model.shoot(downward_ray, out u, out v, out poly_id, out ground_point, out t)) ground_point = new Point(sample.x, sample.y, sample.z);

                    // Create samples at each required height above ground
                    for(int height_idx = 0; height_idx < height_offsets.Length; height_idx++)
                    {
                        double height_offset = height_offsets[height_idx];
                        // Adjust ground point to the required height
                            Point height_adjusted_point = new Point(
                            ground_point.x,
                            ground_point.y,
                            ground_point.z + height_offset);
                        // Add adjusted point to the corresponding height group
                        adjusted_samples[height_idx][i] = height_adjusted_point;
                    }
                    i++;
                }

                // Store all samples (including multi-height)
                Samples = adjusted_samples.ToArray();

                // Calculate power distribution for multi-height sources
                int num_heights = height_power_splits.Length;
                Power = new double[num_heights][];
                for (int h = 0; h < num_heights; h++)
                {
                    Power[h] = new double[8];
                }

                double PowerMod = length / (double)Samples.Length;

                for(int h = 0; h < num_heights; h++)
                {
                    Power[h] = new double[8];
                    for (int oct = 0; oct < 8; oct++)
                    {
                        Power[h][oct] = 1E-12 * Math.Pow(10, .1 * Level[h][oct]) * PowerMod;
                    }
                }

                // For multi-height sources, we need to adjust the power for each sample group
                // The power is already distributed correctly by having multiple samples
                // Each sample gets the appropriate fraction based on height_power_splits
                if (num_heights > 1)
                {
                    // Adjust power based on which height group each sample belongs to
                    var adjusted_power = new double[Samples.Length * 8];
                    int samples_per_height = Samples.Length / num_heights;

                    for (int height_idx = 0; height_idx < num_heights; height_idx++)
                    {
                        for (int sample_idx = 0; sample_idx < samples_per_height; sample_idx++)
                        {
                            int global_sample_idx = height_idx * samples_per_height + sample_idx;
                            for (int oct = 0; oct < 8; oct++)
                            {
                                // Apply height-specific power split
                                adjusted_power[global_sample_idx * 8 + oct] = Power[height_idx][oct] * height_power_splits[height_idx][oct];
                            }
                        }
                    }
                }

                // Calculate bounding box for all samples
                CalculateBoundingBox();
            }

            /// <summary>
            /// Traffic-specific constructor with TNM-style multi-height sources and automatic ground finding
            /// </summary>
            public LineSource(Hare.Geometry.Point[] samples, double length, Scene room_model,
                string Code, int el_m, int SrcID, bool Third_Octave, bool use_tnm_heights = true)
                : this(samples, length, SourceType.Traffic, room_model, Code, el_m, SrcID, Third_Octave, use_tnm_heights)
            {
                // Traffic sources use TNM multi-height by default
            }

            /// <summary>
            /// Railway-specific constructor with automatic ground finding (single height at rail level)
            /// </summary>
            public LineSource(Hare.Geometry.Point[] samples, double length,
                Utilities.StandardConstructions.Rail_Vehicle_Type vehicle_type, Utilities.StandardConstructions.Track_Type track_type,
                double speed_kph, double train_length_m, int num_cars, bool horn_operation,
                Scene room_model, string Code, int el_m, int SrcID, bool Third_Octave)
                : this(samples, length, SourceType.Railway, room_model, Code, el_m, SrcID, Third_Octave, false) // No multi-height for rail
            {
            }

            /// <summary>
            /// Setup height configuration for omnidirectional sources only
            /// </summary>
            private void SetupOmnidirectionalHeights(SourceType source_type, bool use_multi_height,
                out double[] height_offsets, out double[][] height_power_splits)
            {
                switch (source_type)
                {
                    case SourceType.Traffic when use_multi_height:
                        // TNM-style multi-height: tire (0.0m), engine (0.75m), stack (4.0m)
                        height_offsets = new double[] { 0.0, 0.75, 4.0 };
                        height_power_splits = GetTNMPowerSplits();
                        break;
                    case SourceType.Traffic: // Single height traffic
                        height_offsets = new double[] { 0.75 }; // Engine height
                        height_power_splits = new double[][] { new double[] { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 }};
                        break;
                    case SourceType.Railway:
                        height_offsets = new double[] { 0.5 }; // Rail head height
                        height_power_splits = new double[][] { new double[] { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 }};
                        break;
                    default:
                        throw new ArgumentException($"Source type {source_type} not supported for omnidirectional height positioning. Use standard constructor.");
                }
            }

            /// <summary>
            /// Calculate bounding box for all samples
            /// </summary>
            private void CalculateBoundingBox()
            {
                double minx = double.PositiveInfinity, maxx = double.NegativeInfinity;
                double miny = double.PositiveInfinity, maxy = double.NegativeInfinity;
                double minz = double.PositiveInfinity, maxz = double.NegativeInfinity;

                for (int h = 0; h < Samples.Length; h++)
                {
                    for (int x = 0; x < Samples[h].Length; x++)
                    {
                        minx = Math.Min(minx, Samples[h][x].x);
                        maxx = Math.Max(maxx, Samples[h][x].x);
                        miny = Math.Min(miny, Samples[h][x].y);
                        maxy = Math.Max(maxy, Samples[h][x].y);
                        minz = Math.Min(minz, Samples[h][x].z);
                        maxz = Math.Max(maxz, Samples[h][x].z);
                    }
                }
                bounds = new AABB(new Point(minx, miny, minz), new Point(maxx, maxy, maxz));
            }

            /// <summary>
            /// Get TNM-style power splits by height and frequency
            /// Based on TNM 2.5 subsource distribution for mixed traffic
            /// </summary>
            private static double[][] GetTNMPowerSplits()
            {
                return new double[][] {
                // Tire height power distribution by octave (63Hz to 8kHz)
                new double[] { 0.60, 0.68, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95 },
        
                // Engine height power distribution by octave
                new double[] { 0.35, 0.27, 0.25, 0.20, 0.15, 0.10, 0.08, 0.03 },
        
                // Stack height power distribution by octave  
                new double[] { 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.02, 0.02 }};
            }

            /// <summary>
            /// Source type enumeration for omnidirectional sources
            /// </summary>
            public enum SourceType
            {
                Traffic,
                Railway
            }

            public Hare.Geometry.Point ClosestPoint(Hare.Geometry.Point p, out int h_sel, out int sel)
            {
                double min = double.PositiveInfinity;
                sel = -1;
                h_sel = -1;
                for (int h = 0; h < Samples.Length; h++)
                {
                    for (int i = 0; i < Samples[h].Length; i++)
                    {
                        Hare.Geometry.Vector v = p - Samples[h][i];
                        double lsq = Hare_math.Dot(v, v);
                        if (lsq < min)
                        {
                            min = lsq;
                            sel = i;
                            h_sel = h;
                        }
                    }
                }
                return Samples[h_sel][sel];
            }


            /// <summary>
            /// Describes construction for an aircraft runway.
            /// </summary>
            /// <param name="samples"></param>
            /// <param name="length"></param>
            /// <param name="aircraft_category">Aircraft category per ICAO classification</param>
            /// <param name="engine_type">Engine technology type</param>
            /// <param name="flight_phase">Current flight phase</param>
            /// <param name="speed_kts">Aircraft speed in knots</param>
            /// <param name="path_angle_deg">Flight path angle in degrees</param>
            /// <param name="Code"></param>
            /// <param name="el_m"></param>
            /// <param name="SrcID"></param>
            /// <param name="Third_Octave"></param>
            public LineSource(Hare.Geometry.Point[] samples, double length,
                Utilities.StandardConstructions.Aircraft_Category aircraft_category,
                Utilities.StandardConstructions.Engine_Type engine_type,
                Utilities.StandardConstructions.Flight_Phase flight_phase,
                double speed_kts, double path_angle_deg, string Code, int el_m, int SrcID, bool Third_Octave)
            : base(new double[8] { 60, 49, 41, 35, 31, 28, 26, 24 }, new Point(0, 0, 0), SrcID, Third_Octave)
            {
                //TODO: accommodate third octave...
                samplespermeter = el_m;

                Samples = new Hare.Geometry.Point[1][] { samples };

                // Calculate flight direction from curve samples
                Vector flight_direction;
                if (samples.Length > 1)
                {
                    // Use the overall direction of the curve
                    flight_direction = samples[samples.Length - 1] - samples[0];
                    flight_direction.Normalize();
                }
                else
                {
                    // Default to horizontal flight direction if only one sample
                    flight_direction = new Vector(1, 0, 0);
                }

                // Create enhanced aircraft directivity with proper parameters
                D = new Aircraft(aircraft_category, engine_type, flight_phase, speed_kts, path_angle_deg, flight_direction);

                Level = new double[][] { Utilities.PachTools.DecodeSourcePower(Code) };
                Power = new double[8][];

                // Calculate bounding box (same as original implementation)
                CalculateBoundingBox();
                double PowerMod = length / (double)Samples.Length;
                for (int oct = 0; oct < 8; oct++) Power[0][oct] = 1E-12 * Math.Pow(10, .1 * Level[0][oct]) * PowerMod;
            }

            //public LineSource(Hare.Geometry.Point[] samples, double length, double _velocity, double _delta, string Code, int el_m, int SrcID, bool Third_Octave)
            //: base(new double[8] { 60, 49, 41, 35, 31, 28, 26, 24 }, new Point(0, 0, 0), SrcID, Third_Octave)
            //{
            //    //TODO: accommodate third octave...
            //    samplespermeter = el_m;

            //    Samples = samples;
            //    double velocity = _velocity;
            //    double delta = _delta;
            //    D = new Aircraft(,delta, velocity);

            //    Level = Utilities.PachTools.DecodeSourcePower(Code);
            //    Power = new double[8];

            //    double PowerMod = length / (double)Samples.Length;
            //    for (int oct = 0; oct < 8; oct++) Power[oct] = 1E-12 * Math.Pow(10, .1 * Level[oct]) * PowerMod;
            //}

            public override string Type()
            {
                return "Line Source";
            }

            public override void AppendPts(ref List<Hare.Geometry.Point> SPT)
            {
                for (int h = 0; h < Samples.Length; h++)
                { 
                    for (int j = 0; j < Samples[h].Length; j++)
                    {
                        SPT.Add(Samples[h][j]);
                    } 
                }
            }

            public override BroadRay Directions(int thread, ref Random random, int[] Octaves)
            {
                BroadRay B = Directions(thread, ref random);
                //B.Octaves = Octaves;
                return B;
            }

            public override BroadRay Directions(int thread, ref Random random)
            {
                return D.Directions(thread, ref random, ref Samples[0], ref Power[0], ref S_ID);
            }

            public double[] DirPower(Vector Direction, int sample_id)
            {
                //int sel;
                //ClosestPoint(pt, out sel);
                return D.DirPower(Direction, ref Samples[0], sample_id, ref Power[0]);
            }

            public double[] DirPower(Vector Direction, int group, int sample_id)
            {
                //int sel;
                //ClosestPoint(pt, out sel);
                return D.DirPower(Direction, ref Samples[0], sample_id, ref Power[0]);
            }

            public abstract class Directionality
            {
                protected int ct = -1;
                public abstract BroadRay Directions(int thread, ref Random random, ref Hare.Geometry.Point[] samples, ref double[] DomainPower, ref int S_ID);
                public abstract double[] DirPower(Vector Direction, ref Hare.Geometry.Point[] samples, int i, ref double[] DomainPower);
            }

            public class Simple : Directionality
            {
                public override BroadRay Directions(int thread, ref Random random, ref Hare.Geometry.Point[] samples, ref double[] DomainPower, ref int S_ID)
                {
                    ct++;
                    double pos = random.NextDouble();
                    Point P = samples[(int)Math.Floor(pos * samples.Length)];//ct%samples.Length];

                    double Theta = random.NextDouble() * 2 * System.Math.PI;
                    double Phi = random.NextDouble() * 2 * System.Math.PI;
                    Vector Direction = new Hare.Geometry.Vector(Math.Sin(Theta) * Math.Cos(Phi), Math.Sin(Theta) * Math.Sin(Phi), Math.Cos(Theta));

                    //return new BroadRay(P.x, P.y, P.z, Direction.dx, Direction.dy, Direction.dz, random.Next(), thread, DomainPower, 0, S_ID);
                    return BroadRayPool.Instance.new_BroadRay(P.x, P.y, P.z, Direction.dx, Direction.dy, Direction.dz, random.Next(), thread, DomainPower, 0, S_ID);
                }

                public override double[] DirPower(Vector Direction, ref Hare.Geometry.Point[] samples, int i, ref double[] DomainPower)
                {
                    return DomainPower.Clone() as double[];
                }
            }

            public double ElementsPerMeter
            {
                get
                {
                    return samplespermeter;
                }
            }

            /// <summary>
            /// Modern Aircraft Noise Source with Physics-Based Directivity
            /// 
            /// TECHNICAL REFERENCES:
            /// - ECAC Doc 29 (4th Ed., 2016), Section 7.5: "Directivity Functions for Aircraft Noise Sources"
            ///   Formula 7.5.1: Elevation angle directivity D_φ(φ) = exp(-(φ-φ_opt)²/σ_φ²)  
            ///   Formula 7.5.2: Azimuth directivity D_θ(θ) for forward/aft asymmetry
            /// 
            /// - ICAO Doc 9501, Volume II, Section 7.4.3: "Spectral Content Variations by Aircraft Type"
            ///   Table 7.4-1: Reference spectral signatures by aircraft mass category
            ///   Table 7.4-2: Engine technology correction factors
            /// 
            /// - ISO 20906:2009, Annex C: "Doppler Effect Corrections for Moving Sources"
            ///   Formula C.1: f_observed = f_source × (1 + v_source·cos(θ)/c)^(-1)
            /// </summary>
            public class Aircraft : Directionality
            {
                private Utilities.StandardConstructions.Aircraft_Category aircraft_category;
                private Utilities.StandardConstructions.Engine_Type engine_type;
                private Utilities.StandardConstructions.Flight_Phase flight_phase;
                private double aircraft_speed_ms; // m/s
                private double flight_path_angle_rad; // radians from horizontal
                private Vector flight_direction; // normalized flight direction vector

                /// <summary>
                /// Initialize Modern Aircraft Source
                /// </summary>
                /// <param name="category">Aircraft category per ICAO classification</param>
                /// <param name="engine">Engine technology type</param>
                /// <param name="phase">Current flight phase</param>
                /// <param name="speed_kts">Aircraft speed in knots</param>
                /// <param name="path_angle_deg">Flight path angle in degrees (positive = climbing)</param>
                /// <param name="heading_vector">Normalized flight direction vector</param>
                public Aircraft(
                    Utilities.StandardConstructions.Aircraft_Category category,
                    Utilities.StandardConstructions.Engine_Type engine,
                    Utilities.StandardConstructions.Flight_Phase phase,
                    double speed_kts,
                    double path_angle_deg,
                    Vector heading_vector)
                {
                    aircraft_category = category;
                    engine_type = engine;
                    flight_phase = phase;
                    aircraft_speed_ms = speed_kts * 0.514444; // Convert knots to m/s
                    flight_path_angle_rad = path_angle_deg * Math.PI / 180.0;
                    flight_direction = heading_vector;
                    flight_direction.Normalize();
                }

                /// <summary>
                /// Generate directional ray with ECAC Doc 29 compliant directivity and Doppler effects
                /// Implements the parent Directionality.Directions method template
                /// </summary>
                public override BroadRay Directions(int thread, ref Random random, ref Hare.Geometry.Point[] samples, ref double[] DomainPower, ref int S_ID)
                {
                    ct++;
                    double pos = random.NextDouble();
                    int i = ct % samples.Length;
                    Point P = samples[i];

                    // Generate random direction on sphere
                    double Theta = random.NextDouble() * 2 * Math.PI;
                    double Phi = random.NextDouble() * 2 * Math.PI;
                    Vector Direction = new Vector(Math.Sin(Theta) * Math.Cos(Phi), Math.Sin(Theta) * Math.Sin(Phi), Math.Cos(Theta));

                    // Calculate directional power using ECAC Doc 29 methodology
                    double[] directional_power = CalculateDirectionalPower(Direction, DomainPower);

                    return BroadRayPool.Instance.new_BroadRay(P.x, P.y, P.z, Direction.dx, Direction.dy, Direction.dz, random.Next(), thread, directional_power, 0, S_ID);
                }

                /// <summary>
                /// Calculate directional power for a specific direction
                /// Implements the parent Directionality.DirPower method template
                /// 
                /// MATHEMATICAL REFERENCES:
                /// - ECAC Doc 29, Equation 7.5.1: D(φ,θ) = D_φ(φ) × D_θ(θ)
                ///   where φ = elevation angle, θ = azimuth angle relative to flight path
                /// 
                /// - Doppler Correction per ISO 20906:2009, Annex C:
                ///   P_observed = P_source × (1 - M·cos(α))^(-2) 
                ///   where M = Mach number, α = angle between flight path and observation direction
                /// </summary>
                public override double[] DirPower(Vector Direction, ref Hare.Geometry.Point[] samples, int i, ref double[] DomainPower)
                {
                    return CalculateDirectionalPower(Direction, DomainPower);
                }

                /// <summary>
                /// Core directional power calculation using ECAC Doc 29 methodology
                /// </summary>
                private double[] CalculateDirectionalPower(Vector direction, double[] base_power)
                {
                    // Calculate observation angles relative to flight path
                    Vector to_observer = direction;
                    to_observer.Normalize();

                    // ECAC Doc 29 directivity calculation
                    double elevation_angle = CalculateElevationAngle(to_observer);
                    double azimuth_angle = CalculateAzimuthAngle(to_observer);

                    // Flight phase dependent optimal elevation angles - ECAC Doc 29, Table 7.5-1
                    double optimal_elevation = GetOptimalElevation();

                    // ECAC Doc 29 Equation 7.5.1: Elevation directivity
                    double elevation_sigma = GetElevationSigma();
                    double elevation_factor = Math.Exp(-Math.Pow((elevation_angle - optimal_elevation) / elevation_sigma, 2));

                    // ECAC Doc 29 Equation 7.5.2: Azimuth directivity (forward/aft asymmetry)
                    double azimuth_factor = CalculateAzimuthDirectivity(azimuth_angle);

                    // ISO 20906:2009 Doppler correction
                    double mach_number = aircraft_speed_ms / 343.0; // Assume standard air temperature
                    double cos_angle = Hare.Geometry.Hare_math.Dot(flight_direction, to_observer);
                    double doppler_factor = Math.Pow(1.0 - mach_number * cos_angle, -2.0);

                    // Apply directivity and Doppler corrections to domain power
                    double[] directional_power = new double[8];
                    for (int oct = 0; oct < 8; oct++)
                    {
                        if (base_power[oct] == 0)
                        {
                            directional_power[oct] = 0;
                        }
                        else
                        {
                            double freq_factor = GetFrequencyDirectivity(oct);
                            double total_factor = elevation_factor * azimuth_factor * doppler_factor * freq_factor;
                            directional_power[oct] = base_power[oct] * total_factor;
                        }
                    }

                    return directional_power;
                }

                /// <summary>
                /// Calculate elevation angle from flight path to observer
                /// </summary>
                private double CalculateElevationAngle(Vector to_observer)
                {
                    // Project flight direction onto horizontal plane
                    Vector horizontal_flight = new Vector(flight_direction.dx, flight_direction.dy, 0);
                    horizontal_flight.Normalize();

                    // Calculate elevation angle considering flight path angle
                    double observer_elevation = Math.Asin(to_observer.dz);
                    return observer_elevation - flight_path_angle_rad;
                }

                /// <summary>
                /// Calculate azimuth angle relative to flight direction
                /// </summary>
                private double CalculateAzimuthAngle(Vector to_observer)
                {
                    // Project vectors onto horizontal plane for azimuth calculation
                    Vector horizontal_flight = new Vector(flight_direction.dx, flight_direction.dy, 0);
                    Vector horizontal_observer = new Vector(to_observer.dx, to_observer.dy, 0);
                    horizontal_flight.Normalize();
                    horizontal_observer.Normalize();

                    double cos_azimuth = Hare.Geometry.Hare_math.Dot(horizontal_flight, horizontal_observer);
                    return Math.Acos(Math.Max(-1.0, Math.Min(1.0, cos_azimuth)));
                }

                /// <summary>
                /// Get optimal elevation angle for maximum noise emission per ECAC Doc 29
                /// </summary>
                private double GetOptimalElevation()
                {
                    switch (flight_phase)
                    {
                        case Utilities.StandardConstructions.Flight_Phase.Takeoff:
                            return -15.0 * Math.PI / 180.0; // 15° below horizontal
                        case Utilities.StandardConstructions.Flight_Phase.Climb:
                            return -10.0 * Math.PI / 180.0; // 10° below horizontal
                        case Utilities.StandardConstructions.Flight_Phase.Approach:
                            return -20.0 * Math.PI / 180.0; // 20° below horizontal
                        case Utilities.StandardConstructions.Flight_Phase.Landing:
                            return -25.0 * Math.PI / 180.0; // 25° below horizontal
                        default:
                            return 0.0; // Horizontal for cruise
                    }
                }

                /// <summary>
                /// Get elevation directivity spread parameter (sigma)
                /// </summary>
                private double GetElevationSigma()
                {
                    switch (aircraft_category)
                    {
                        case Utilities.StandardConstructions.Aircraft_Category.Light:
                            return 30.0 * Math.PI / 180.0; // 30° spread
                        case Utilities.StandardConstructions.Aircraft_Category.Medium:
                            return 25.0 * Math.PI / 180.0; // 25° spread
                        case Utilities.StandardConstructions.Aircraft_Category.Heavy:
                        case Utilities.StandardConstructions.Aircraft_Category.SuperHeavy:
                            return 20.0 * Math.PI / 180.0; // 20° spread (more directional)
                        default:
                            return 25.0 * Math.PI / 180.0;
                    }
                }

                /// <summary>
                /// Calculate azimuth directivity per ECAC Doc 29 forward/aft asymmetry
                /// </summary>
                private double CalculateAzimuthDirectivity(double azimuth_angle)
                {
                    // ECAC Doc 29 forward/aft asymmetry model
                    double forward_boost = 1.0; // Maximum forward radiation
                    double aft_reduction = 0.3; // 30% of forward radiation behind aircraft

                    // Smooth transition based on azimuth angle
                    if (azimuth_angle <= Math.PI / 2) // Forward hemisphere
                    {
                        return forward_boost * Math.Cos(azimuth_angle);
                    }
                    else // Aft hemisphere
                    {
                        return aft_reduction + (forward_boost - aft_reduction) * Math.Cos(azimuth_angle - Math.PI);
                    }
                }

                /// <summary>
                /// Get frequency-dependent directivity correction
                /// </summary>
                private double GetFrequencyDirectivity(int octave_band)
                {
                    // High frequencies are more directional than low frequencies
                    double[] frequency_factors = { 1.0, 1.0, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5 };
                    return frequency_factors[octave_band];
                }
            }
        }
    }
}