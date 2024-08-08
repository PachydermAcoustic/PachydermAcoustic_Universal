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
using System.IO;

namespace Pachyderm_Acoustic
{
    public partial class Pach_Properties
    {
        string SettingsPath;

        int GeometrySystem;
        int Priority;
        int Thread_Spec;
        int Spatial_Partition;

        int VGDomain;
        int OctDepth;
        string LibraryPath;
        bool Save;
        int FiltPhase;

        private static Pach_Properties instance;

        public static Pach_Properties Instance
        {
            get
            {
                if (instance == null)
                {
                    instance = new Pach_Properties();
                }
                return instance;
            }
        }

        private Pach_Properties()
        {
            SettingsPath = System.Environment.GetFolderPath(System.Environment.SpecialFolder.ApplicationData) + "\\Pachyderm";
            if (!System.IO.Directory.Exists(SettingsPath)) System.IO.Directory.CreateDirectory(SettingsPath);
            SettingsPath += "\\Pach_Settings.pset";

            //FileStream S = new FileStream(SettingsPath, FileMode.OpenOrCreate);
            //BinaryReader Reader = new BinaryReader(S);
            //BinaryWriter Writer = new BinaryWriter(S);

            Recover_Settings();
        }

        public Pach_Properties(string Path_to_Settings_Folder)
        {
            SettingsPath = Path_to_Settings_Folder;
            SettingsPath += "\\Pach_Settings.pset";
            Recover_Settings();
        }

        private void Recover_Settings()
        {
            FileStream S = new FileStream(SettingsPath, FileMode.OpenOrCreate);
            BinaryReader Reader = new BinaryReader(S);
            BinaryWriter Writer = new BinaryWriter(S);

            try
            {
                string p_version = Reader.ReadString();
                string newversion = Utilities.PachTools.Version();
                if (p_version != newversion) newSettings(ref Reader, ref Writer);

                //2. Geometry System(int)
                GeometrySystem = Reader.ReadInt32();

                //switch (Reader.ReadInt32())
                //{
                //    case 1:
                //        GEO_NURBS.Checked = true;
                //        break;
                //    case 2:
                //        GEO_MESH.Checked = true;
                //        break;
                //    case 3:
                //        GEO_MR_MESH.Checked = true;
                //        break;
                //}

                //3. Processing
                Priority = Reader.ReadInt32();
                //switch (Reader.ReadInt32())
                //{
                //    case 1:
                //        PR_Single.Checked = true;
                //        break;
                //    case 2:
                //        PR_MULTI_ALL.Checked = true;
                //        break;
                //    case 3:
                //        PR_MULTI_EXP.Checked = true;
                //        break;
                //}

                //4. ThreadCount(int)
                Thread_Spec = Reader.ReadInt32();
                //Thread_Spec.Value = Reader.ReadInt32();

                //5. Spatial Partition Selection (int)
                Spatial_Partition = Reader.ReadInt32();
                //switch (Reader.ReadInt32())
                //{
                //    case 1:
                //        VGSP_CHECK.Checked = true;
                //        break;
                //    case 2:
                //        OCT_CHECK.Checked = true;
                //        break;
                //}
                //6. Voxel Grid Domain(int)

                VGDomain = Reader.ReadInt32();
                //7. Octree Depth(int)
                //OCT_DEPTH.Value = Reader.ReadInt32();
                OctDepth = Reader.ReadInt32();
                //8. Material Library Path
                LibraryPath = Reader.ReadString();
                //9. Save Results after simulation?
                Save = Reader.ReadBoolean();
                //10. Save Filter Method
                FiltPhase = Reader.ReadInt32();
                //switch (Reader.ReadInt32())
                //{
                //    case 1:
                //        Filt_LinearPhase.Checked = true;
                //        Audio.Pach_SP.Filter = new Audio.Pach_SP.Linear_Phase_System();
                //        break;
                //    case 2:
                //        Filt_MinPhase.Checked = true;
                //        Audio.Pach_SP.Filter = new Audio.Pach_SP.Minimum_Phase_System();
                //        break;
                //}
                Reader.Close();
                Writer.Close();
            }
            catch (Exception x)
            {
                //if (x is EndOfStreamException)
                //{ 
                    newSettings(ref Reader, ref Writer); //}
                //else
                //{ throw x; }
            }
        }

        /// <summary>
        /// This method creates a new settings file for Pachyderm set to certain defaults.
        /// Called when previous methods find that there is no suitable settings file.
        /// </summary>
        /// <param name="Reader"></param>
        /// <param name="Writer"></param>
        private void newSettings(ref BinaryReader Reader, ref BinaryWriter Writer)
        {
            Reader.Close();
            Writer.Close();
            FileStream S = new FileStream(SettingsPath, FileMode.Create);
            Writer = new BinaryWriter(S);
            //1. Plugin Version(string)
            Writer.Write(Utilities.PachTools.Version());
            //2. Geometry System(int)
            GeometrySystem = 2;
            Writer.Write(2);
            //3. Priority(int)
            //0 = 'High'
            //1 = 'AboveNormal'
            //2 = 'Normal'
            Priority = 2;
            Writer.Write(2);
            //4. ThreadCount(int)
            Thread_Spec = System.Environment.ProcessorCount;
            Writer.Write(System.Environment.ProcessorCount);
            //5. Spatial Partition Selection (int)
            Spatial_Partition = 1;
            Writer.Write(1);
            //6. Voxel Grid Domain(int)
            VGDomain = 7;
            Writer.Write(7);
            //7. Octree Depth(int)
            OctDepth = 3;
            Writer.Write(3);
            //8. Material Library Path
            LibraryPath = System.Environment.GetFolderPath(System.Environment.SpecialFolder.ApplicationData) + "\\Pachyderm";
            //mlPath = System.IO.Path.GetDirectoryName(mlPath);
            Writer.Write(LibraryPath);
            //9. Save Results after simulation?
            Save = false;
            Writer.Write(false);
            //10. Save Filter Method
            FiltPhase = 1;
            { Writer.Write(1); }
            Reader.Close();
            Writer.Close();
        }

        /// <summary>
        /// Get the number of processors/threads to be used.
        /// </summary>
        /// <returns></returns>
        public int TaskPriority()
        {
            return Priority;
        }

        /// <summary>
        /// Returns user's choice for geometry system.
        /// </summary>
        /// <returns></returns>
        public int Geometry_Spec()
        {
            if (GeometrySystem == 1) return 0; //NURBS
            else if (GeometrySystem == 2) return 1; //MESH via Hare
            else if (GeometrySystem == 3) return 2; // Unimplemented Multi-resolution system.
            throw new Exception("Is there a new Geometrical Option that needs to be implemented?");
        }

        public int Priority_Choice
        {
            get 
            {
                return Priority; 
            }
            set 
            {
                Priority = value;
                Commit_Settings();
            }
        }

        /// <summary>
        /// Specifies the geometry system... 1 is NURBS (rhino specific). 2 is Mesh. 3 is Multi-resolution mesh (not yet implemented).
        /// </summary>
        public int Geometry_System
        {
            get 
            {
                return GeometrySystem; 
            }
            set
            {
                GeometrySystem = value;
                Commit_Settings();
            }
        }


        /// <summary>
        /// Returns the user selected spatial partition. 1 is for Voxel Grid. 2 is for Octree (not yet implemented.).
        /// </summary>
        /// <returns></returns>
        public int Spatial_Optimization
        {
            get
            {
                if (Spatial_Partition == 1) return 0;
                if (Spatial_Partition == 2) return 1;
                throw new Exception("Is there a new spatial partition to be implemented?");
            }
            set 
            {
                Spatial_Partition = value;
                Commit_Settings();
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public int ThreadCount
        {
            get
            {
                return Thread_Spec;
            }
            set
            {
                Thread_Spec = value;
                Commit_Settings();
            }
        }

        /// <summary>
        /// Filter Phase Regime. 
        /// </summary>
        public int FilterPhase
        {
            get
            {
                return FiltPhase;
            }
            set
            {
                FiltPhase = value;
                if (FiltPhase == 1) Audio.Pach_SP.Filter = new Audio.Pach_SP.Linear_Phase_System();
                else Audio.Pach_SP.Filter = new Audio.Pach_SP.Minimum_Phase_System();
                
                Commit_Settings();
            }
        }


        /// <summary>
        /// Returns the user selected domain of the voxel grid (how many times to subdivide voxels. Functionally 8^x voxels.)
        /// </summary>
        /// <returns></returns>
        public int VoxelGrid_Domain
        {
            get
            {
                return (int)VGDomain;
            }
            set
            {
                VGDomain = value;
                Commit_Settings();
            }
        }


        /// <summary>
        /// If octree is chosen, returns the depth of the octree (how many branches deep)
        /// </summary>
        /// <returns></returns>
        public int Oct_Depth
        {
            get
            {
                return (int)OctDepth;
            }
            set 
            {
                OctDepth = value;
                Commit_Settings();
            }
        }

        /// <summary>
        /// Gets the path to the material library.
        /// </summary>
        /// <returns></returns>
        public string Lib_Path
        {
            get
            {
                return LibraryPath;
            }
            set 
            {
                LibraryPath = value;
                Commit_Settings();
            }
        }

        /// <summary>
        /// Determines whether Pachyderm will save results automatically upon starting a simulation.
        /// </summary>
        public bool SaveResults
        {
            get
            {
                return Save;
            }
            set 
            {
                Save = value; 
                Commit_Settings();
            }
        }

        /// <summary>
        /// Records settings to the pachyderm settings file
        /// </summary>
        private void Commit_Settings()
        {
            File.Delete(SettingsPath);
            FileStream S = new FileStream(SettingsPath, FileMode.CreateNew);
            BinaryWriter Writer = new BinaryWriter(S);
            //1. Plugin Version(string)
            Writer.Write(Utilities.PachTools.Version());
            //2. Geometry System(int)
            //if (GEO_NURBS.Checked == true)
            //{ Writer.Write(1); }
            //else if (GEO_MESH.Checked == true)
            //{ Writer.Write(2); }
            //else
            //{ Writer.Write(3); }
            Writer.Write(GeometrySystem);

            //3. Processing(int)
            //if (PR_Single.Checked == true)
            //{ Writer.Write(1); }
            //else if (PR_MULTI_ALL.Checked == true)
            //{ Writer.Write(2); }
            //else
            //{ Writer.Write(3); }
            Writer.Write(Priority);

            //4. ThreadCount(int)
            //Writer.Write((int)Thread_Spec.Value);
            Writer.Write(Thread_Spec);

            //5. Spatial Partition Selection (int)
            //if (VGSP_CHECK.Checked == true)
            //{ Writer.Write(1); }
            //else
            //{ Writer.Write(2); }
            Writer.Write(Spatial_Partition);

            //6. Voxel Grid Domain(int)
            Writer.Write(VGDomain);
            //7. Octree Depth(int)
            Writer.Write(OctDepth);
            //8. Material Library Path
            Writer.Write(LibraryPath);
            //9. Save Results after simulation?
            Writer.Write(Save);

            //10. Save Filter Method
            if (FiltPhase == 1)
            {
                Audio.Pach_SP.Filter = new Audio.Pach_SP.Linear_Phase_System();
                Writer.Write(1);
            }
            else
            {
                Audio.Pach_SP.Filter = new Audio.Pach_SP.Minimum_Phase_System();
                Writer.Write(2);
            }
            Writer.Close();
        }
    }
}