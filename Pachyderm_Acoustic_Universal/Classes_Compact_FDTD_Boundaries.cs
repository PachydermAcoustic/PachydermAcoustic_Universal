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
using System.Runtime.Serialization;
using Hare.Geometry;

namespace Pachyderm_Acoustic
{
    namespace Numeric
    {
        namespace TimeDomain
        {
            public partial class Acoustic_Compact_FDTD : Simulation_Type
            {
                public abstract class Bound_Node : Node
                {
                    public Bound_Node(Point loc)
                    : base(loc)
                    { }

                    [Flags]
                    public enum Boundary
                    {
                        None = 0x0000000,
                        AXPos = 0x0000001,
                        AXNeg = 0x0000002,
                        AYPos = 0x0000004,
                        AYNeg = 0x0000008,
                        AZPos = 0x0000010,
                        AZNeg = 0x0000020,
                        SDXPosYPos = 0x0000040,
                        SDXPosYNeg = 0x0000080,
                        SDXNegYPos = 0x0000100,
                        SDXNegYNeg = 0x0000200,
                        SDXPosZPos = 0x0000400,
                        SDXPosZNeg = 0x0000800,
                        SDXNegZPos = 0x0001000,
                        SDXNegZNeg = 0x0002000,
                        SDYPosZPos = 0x0004000,
                        SDYPosZNeg = 0x0008000,
                        SDYNegZPos = 0x0010000,
                        SDYNegZNeg = 0x0020000,
                        DXPosYPosZPos = 0x0040000,
                        DXNegYPosZPos = 0x0080000,
                        DXPosYNegZPos = 0x0100000,
                        DXNegYNegZPos = 0x0200000,
                        DXPosYPosZNeg = 0x0400000,
                        DXNegYPosZNeg = 0x0800000,
                        DXPosYNegZNeg = 0x1000000,
                        DXNegYNegZNeg = 0x2000000
                    }
                    #region Quick Recognition Boundary Conditions
                    public static Boundary[] Face_Combos = new Boundary[10]
                    {
                    Boundary.AXNeg | Boundary.AXPos | Boundary.AYNeg | Boundary.AYPos | Boundary.AZNeg | Boundary.AZPos,
                    Boundary.AXNeg | Boundary.AXPos,
                    Boundary.AYNeg | Boundary.AYPos,
                    Boundary.AZNeg | Boundary.AZPos,
                    Boundary.AXNeg,
                    Boundary.AXPos,
                    Boundary.AYNeg,
                    Boundary.AYPos,
                    Boundary.AZNeg,
                    Boundary.AZPos
                    };

                    public static Boundary[] Wall = new Boundary[6]
                    {
                        Boundary.AZPos | Boundary.DXNegYNegZPos | Boundary.DXNegYPosZPos | Boundary.DXPosYNegZPos | Boundary.DXPosYPosZPos | Boundary.SDXNegZPos | Boundary.SDXPosZPos | Boundary.SDYNegZPos | Boundary.SDYPosZPos,
                        Boundary.AZNeg | Boundary.DXNegYNegZNeg | Boundary.DXNegYPosZNeg | Boundary.DXPosYNegZNeg | Boundary.DXPosYPosZNeg | Boundary.SDXNegZNeg | Boundary.SDXPosZNeg | Boundary.SDYNegZNeg | Boundary.SDYPosZNeg,
                        Boundary.AYPos | Boundary.DXNegYPosZNeg | Boundary.DXNegYPosZPos | Boundary.DXPosYPosZNeg | Boundary.DXPosYPosZPos | Boundary.SDXNegYPos | Boundary.SDXPosYPos | Boundary.SDYPosZNeg | Boundary.SDYPosZPos,
                        Boundary.AYNeg | Boundary.DXNegYNegZNeg | Boundary.DXNegYNegZPos | Boundary.DXPosYNegZNeg | Boundary.DXPosYNegZPos | Boundary.SDXNegYNeg | Boundary.SDXPosYNeg | Boundary.SDYNegZNeg | Boundary.SDYNegZPos,
                        Boundary.AXPos | Boundary.DXPosYNegZNeg | Boundary.DXPosYNegZPos | Boundary.DXPosYPosZNeg | Boundary.DXPosYPosZPos | Boundary.SDXPosYNeg | Boundary.SDXPosYPos | Boundary.SDXPosZNeg | Boundary.SDXPosZPos,
                        Boundary.AXNeg | Boundary.DXNegYNegZNeg | Boundary.DXNegYNegZPos | Boundary.DXNegYPosZNeg | Boundary.DXNegYPosZPos | Boundary.SDXNegYNeg | Boundary.SDXNegYPos | Boundary.SDXNegZNeg | Boundary.SDXNegZPos
                    };

                    public static Boundary[] Outer_Edges = new Boundary[12]
                    {
                        Boundary.AZPos | Boundary.DXNegYNegZPos | Boundary.DXNegYPosZPos | Boundary.DXPosYNegZPos | Boundary.DXPosYPosZPos | Boundary.SDXNegZPos | Boundary.SDXPosZPos | Boundary.SDYNegZPos | Boundary.SDYPosZPos |    Boundary.AYPos | Boundary.SDXNegYPos | Boundary.SDXPosYPos |    Boundary.SDYPosZNeg | Boundary.DXNegYPosZNeg | Boundary.DXPosYPosZNeg,
                        Boundary.AZPos | Boundary.DXNegYNegZPos | Boundary.DXNegYPosZPos | Boundary.DXPosYNegZPos | Boundary.DXPosYPosZPos | Boundary.SDXNegZPos | Boundary.SDXPosZPos | Boundary.SDYNegZPos | Boundary.SDYPosZPos |    Boundary.AYNeg | Boundary.SDXNegYNeg | Boundary.SDXPosYNeg |    Boundary.SDYNegZNeg | Boundary.DXNegYNegZNeg | Boundary.DXPosYNegZNeg,
                        Boundary.AZPos | Boundary.DXNegYNegZPos | Boundary.DXNegYPosZPos | Boundary.DXPosYNegZPos | Boundary.DXPosYPosZPos | Boundary.SDXNegZPos | Boundary.SDXPosZPos | Boundary.SDYNegZPos | Boundary.SDYPosZPos |    Boundary.AXPos | Boundary.SDXPosYNeg | Boundary.SDXPosYPos |    Boundary.SDXPosZNeg | Boundary.DXPosYNegZNeg | Boundary.DXPosYPosZNeg,
                        Boundary.AZPos | Boundary.DXNegYNegZPos | Boundary.DXNegYPosZPos | Boundary.DXPosYNegZPos | Boundary.DXPosYPosZPos | Boundary.SDXNegZPos | Boundary.SDXPosZPos | Boundary.SDYNegZPos | Boundary.SDYPosZPos |    Boundary.AXNeg | Boundary.SDXNegYNeg | Boundary.SDXNegYPos |    Boundary.SDXNegZNeg | Boundary.DXNegYNegZNeg | Boundary.DXNegYPosZNeg,
                        Boundary.AZNeg | Boundary.DXNegYNegZNeg | Boundary.DXNegYPosZNeg | Boundary.DXPosYNegZNeg | Boundary.DXPosYPosZNeg | Boundary.SDXNegZNeg | Boundary.SDXPosZNeg | Boundary.SDYNegZNeg | Boundary.SDYPosZNeg |    Boundary.AYPos | Boundary.SDXNegYPos | Boundary.SDXPosYPos |    Boundary.SDYPosZPos | Boundary.DXNegYPosZPos | Boundary.DXPosYPosZPos,
                        Boundary.AZNeg | Boundary.DXNegYNegZNeg | Boundary.DXNegYPosZNeg | Boundary.DXPosYNegZNeg | Boundary.DXPosYPosZNeg | Boundary.SDXNegZNeg | Boundary.SDXPosZNeg | Boundary.SDYNegZNeg | Boundary.SDYPosZNeg |    Boundary.AYNeg | Boundary.SDXNegYNeg | Boundary.SDXPosYNeg |    Boundary.SDYNegZPos | Boundary.DXNegYNegZPos | Boundary.DXPosYNegZPos,
                        Boundary.AZNeg | Boundary.DXNegYNegZNeg | Boundary.DXNegYPosZNeg | Boundary.DXPosYNegZNeg | Boundary.DXPosYPosZNeg | Boundary.SDXNegZNeg | Boundary.SDXPosZNeg | Boundary.SDYNegZNeg | Boundary.SDYPosZNeg |    Boundary.AXPos | Boundary.SDXPosYNeg | Boundary.SDXPosYPos |    Boundary.SDXPosZPos | Boundary.DXPosYNegZPos | Boundary.DXPosYPosZPos,
                        Boundary.AZNeg | Boundary.DXNegYNegZNeg | Boundary.DXNegYPosZNeg | Boundary.DXPosYNegZNeg | Boundary.DXPosYPosZNeg | Boundary.SDXNegZNeg | Boundary.SDXPosZNeg | Boundary.SDYNegZNeg | Boundary.SDYPosZNeg |    Boundary.AXNeg | Boundary.SDXNegYNeg | Boundary.SDXNegYPos |    Boundary.SDXNegZPos | Boundary.DXNegYNegZPos | Boundary.DXNegYPosZPos,
                        Boundary.DXNegYNegZPos | Boundary.SDXNegZPos | Boundary.DXNegYPosZPos | Boundary.SDYPosZPos | Boundary.DXPosYPosZPos | Boundary.SDXNegYNeg | Boundary.AXNeg | Boundary.SDXNegYPos | Boundary.AYPos | Boundary.SDXPosYPos | Boundary.DXNegYNegZNeg | Boundary.SDXNegZNeg | Boundary.DXNegYPosZNeg | Boundary.SDYPosZNeg | Boundary.DXPosYPosZNeg,
                        Boundary.DXNegYPosZPos | Boundary.SDYPosZPos | Boundary.DXPosYPosZPos | Boundary.SDXPosZPos | Boundary.DXPosYNegZPos | Boundary.SDXNegYPos | Boundary.AYPos | Boundary.SDXPosYPos | Boundary.AXPos | Boundary.SDXPosYNeg | Boundary.DXNegYPosZNeg | Boundary.SDYPosZNeg | Boundary.DXPosYPosZNeg | Boundary.SDXPosZNeg | Boundary.DXPosYNegZNeg,
                        Boundary.DXPosYPosZPos | Boundary.SDXPosZPos | Boundary.DXPosYNegZPos | Boundary.SDYNegZPos | Boundary.DXNegYNegZPos | Boundary.SDXPosYPos | Boundary.AXPos | Boundary.SDXPosYNeg | Boundary.AYNeg | Boundary.SDXNegYNeg | Boundary.DXPosYPosZNeg | Boundary.SDXPosZNeg | Boundary.DXPosYNegZNeg | Boundary.SDYNegZNeg | Boundary.DXNegYNegZNeg,
                        Boundary.DXPosYNegZPos | Boundary.SDYNegZPos | Boundary.DXNegYNegZPos | Boundary.SDXNegZPos | Boundary.DXNegYPosZPos | Boundary.SDXPosYNeg | Boundary.AYNeg | Boundary.SDXNegYNeg | Boundary.AXNeg | Boundary.SDXNegYPos | Boundary.DXPosYNegZNeg | Boundary.SDYNegZNeg | Boundary.DXNegYNegZNeg | Boundary.SDXNegZNeg | Boundary.DXNegYPosZNeg
                    };

                    public static Boundary[] OuterCorner = new Boundary[8]
                    {
                        Boundary.AXNeg | Boundary.AYPos | Boundary.AZNeg | Boundary.SDXPosYPos | Boundary.SDXNegYPos | Boundary.SDXNegYNeg | Boundary.SDXPosZNeg | Boundary.SDXNegZPos | Boundary.SDXNegZNeg | Boundary.SDYPosZPos | Boundary.SDYPosZNeg | Boundary.SDYNegZNeg | Boundary.DXPosYPosZPos | Boundary.DXNegYPosZPos | Boundary.DXNegYNegZPos | Boundary.DXPosYPosZNeg | Boundary.DXNegYPosZNeg | Boundary.DXPosYNegZNeg | Boundary.DXNegYNegZNeg,
                        Boundary.AXPos | Boundary.AYPos | Boundary.AZNeg | Boundary.SDXPosYPos | Boundary.SDXPosYNeg | Boundary.SDXNegYPos | Boundary.SDXPosZPos | Boundary.SDXPosZNeg | Boundary.SDXNegZNeg | Boundary.SDYPosZPos | Boundary.SDYPosZNeg | Boundary.SDYNegZNeg | Boundary.DXPosYPosZPos | Boundary.DXNegYPosZPos | Boundary.DXPosYNegZPos | Boundary.DXPosYPosZNeg | Boundary.DXNegYPosZNeg | Boundary.DXPosYNegZNeg | Boundary.DXNegYNegZNeg,
                        Boundary.AXPos | Boundary.AYNeg | Boundary.AZPos | Boundary.SDXPosYPos | Boundary.SDXPosYNeg | Boundary.SDXNegYNeg | Boundary.SDXPosZPos | Boundary.SDXPosZNeg | Boundary.SDXNegZPos | Boundary.SDYPosZPos | Boundary.SDYNegZPos | Boundary.SDYNegZNeg | Boundary.DXPosYPosZPos | Boundary.DXNegYPosZPos | Boundary.DXPosYNegZPos | Boundary.DXNegYNegZPos | Boundary.DXPosYPosZNeg | Boundary.DXPosYNegZNeg | Boundary.DXNegYNegZNeg,
                        Boundary.AXNeg | Boundary.AYNeg | Boundary.AZPos | Boundary.SDXPosYNeg | Boundary.SDXNegYPos | Boundary.SDXNegYNeg | Boundary.SDXPosZPos | Boundary.SDXNegZPos | Boundary.SDXNegZNeg | Boundary.SDYPosZPos | Boundary.SDYNegZPos | Boundary.SDYNegZNeg | Boundary.DXPosYPosZPos | Boundary.DXNegYPosZPos | Boundary.DXPosYNegZPos | Boundary.DXNegYNegZPos | Boundary.DXNegYPosZNeg | Boundary.DXPosYNegZNeg | Boundary.DXNegYNegZNeg,
                        Boundary.AXPos | Boundary.AYNeg | Boundary.AZNeg | Boundary.SDXPosYPos | Boundary.SDXPosYNeg | Boundary.SDXNegYNeg | Boundary.SDXPosZPos | Boundary.SDXPosZNeg | Boundary.SDXNegZNeg | Boundary.SDYPosZNeg | Boundary.SDYNegZPos | Boundary.SDYNegZNeg | Boundary.DXPosYPosZPos | Boundary.DXPosYNegZPos | Boundary.DXNegYNegZPos | Boundary.DXPosYPosZNeg | Boundary.DXNegYPosZNeg | Boundary.DXPosYNegZNeg | Boundary.DXNegYNegZNeg,
                        Boundary.AXNeg | Boundary.AYNeg | Boundary.AZNeg | Boundary.SDXPosYNeg | Boundary.SDXNegYPos | Boundary.SDXNegYNeg | Boundary.SDXPosZNeg | Boundary.SDXNegZPos | Boundary.SDXNegZNeg | Boundary.SDYPosZNeg | Boundary.SDYNegZPos | Boundary.SDYNegZNeg | Boundary.DXNegYPosZPos | Boundary.DXPosYNegZPos | Boundary.DXNegYNegZPos | Boundary.DXPosYPosZNeg | Boundary.DXNegYPosZNeg | Boundary.DXPosYNegZNeg | Boundary.DXNegYNegZNeg,
                        Boundary.AXNeg | Boundary.AYPos | Boundary.AZPos | Boundary.SDXPosYPos | Boundary.SDXNegYPos | Boundary.SDXNegYNeg | Boundary.SDXPosZPos | Boundary.SDXNegZPos | Boundary.SDXNegZNeg | Boundary.SDYPosZPos | Boundary.SDYPosZNeg | Boundary.SDYNegZPos | Boundary.DXPosYPosZPos | Boundary.DXNegYPosZPos | Boundary.DXPosYNegZPos | Boundary.DXNegYNegZPos | Boundary.DXPosYPosZNeg | Boundary.DXNegYPosZNeg | Boundary.DXNegYNegZNeg,
                        Boundary.AXPos | Boundary.AYPos | Boundary.AZPos | Boundary.SDXPosYPos | Boundary.SDXPosYNeg | Boundary.SDXNegYPos | Boundary.SDXPosZPos | Boundary.SDXPosZNeg | Boundary.SDXNegZPos | Boundary.SDYPosZPos | Boundary.SDYPosZNeg | Boundary.SDYNegZPos | Boundary.DXPosYPosZPos | Boundary.DXNegYPosZPos | Boundary.DXPosYNegZPos | Boundary.DXNegYNegZPos | Boundary.DXNegYPosZNeg | Boundary.DXPosYNegZNeg | Boundary.DXPosYPosZNeg
                    };

                    public static Boundary[] InnerCorner = new Boundary[12]
                    {
                        Boundary.DXNegYNegZNeg | Boundary.DXPosYPosZPos,
                        Boundary.DXNegYNegZPos | Boundary.DXPosYPosZNeg,
                        Boundary.DXNegYPosZPos | Boundary.DXPosYNegZNeg,
                        Boundary.DXNegYPosZNeg | Boundary.DXPosYNegZPos,
                        Boundary.DXPosYPosZPos,
                        Boundary.DXNegYNegZNeg,
                        Boundary.DXNegYNegZPos,
                        Boundary.DXPosYPosZNeg,
                        Boundary.DXNegYPosZPos,
                        Boundary.DXPosYNegZNeg,
                        Boundary.DXPosYNegZPos,
                        Boundary.DXNegYPosZNeg
                    };

                    public static Boundary[] InnerEdge = new Boundary[18]
                    {
                        Boundary.SDXNegYNeg | Boundary.SDXPosYPos,
                        Boundary.SDXNegYPos | Boundary.SDXPosYNeg,
                        Boundary.SDXNegZNeg | Boundary.SDXPosZPos,
                        Boundary.SDXNegZPos | Boundary.SDXPosZNeg,
                        Boundary.SDYPosZPos | Boundary.SDYNegZNeg,
                        Boundary.SDYPosZNeg | Boundary.SDYNegZPos,
                        Boundary.SDXPosYPos,
                        Boundary.SDXNegYNeg,
                        Boundary.SDXPosYNeg,
                        Boundary.SDXNegYPos,
                        Boundary.SDXPosZPos,
                        Boundary.SDXNegZNeg,
                        Boundary.SDXNegZPos,
                        Boundary.SDXPosZNeg,
                        Boundary.SDYPosZPos,
                        Boundary.SDYNegZNeg,
                        Boundary.SDYNegZPos,
                        Boundary.SDYPosZNeg
                    };

                    public static Boundary[] EdgeCombos = new Boundary[12]
                    {
                        Boundary.DXPosYPosZNeg | Boundary.DXPosYPosZPos | Boundary.SDXPosYPos,
                        Boundary.DXNegYNegZNeg | Boundary.DXNegYNegZPos | Boundary.SDXNegYNeg,
                        Boundary.DXPosYNegZNeg | Boundary.DXPosYNegZPos | Boundary.SDXPosYNeg,
                        Boundary.DXNegYPosZNeg | Boundary.DXNegYPosZPos | Boundary.SDXNegYPos,

                        Boundary.DXPosYNegZPos | Boundary.DXPosYPosZPos | Boundary.SDXPosZPos,
                        Boundary.DXNegYNegZNeg | Boundary.DXNegYPosZNeg | Boundary.SDXNegZNeg,
                        Boundary.DXNegYNegZPos | Boundary.DXNegYPosZPos | Boundary.SDXNegZPos,
                        Boundary.DXPosYPosZNeg | Boundary.DXPosYNegZNeg | Boundary.SDXPosZNeg,

                        Boundary.DXNegYPosZPos | Boundary.DXPosYPosZPos | Boundary.SDYPosZPos,
                        Boundary.DXPosYNegZNeg | Boundary.DXNegYNegZNeg | Boundary.SDYNegZNeg,
                        Boundary.DXNegYPosZNeg | Boundary.DXPosYPosZNeg | Boundary.SDYPosZNeg,
                        Boundary.DXNegYNegZPos | Boundary.DXPosYNegZPos | Boundary.SDYNegZPos
                    };
                #endregion
            }

                /// <summary>
                /// Default boundary node - This node is a dirichlet condition. If you need absorption, use Bound_Node_RDD_MaterialFilter.
                /// </summary>
                public class Bound_Node_RDD: RDD_Node
                {
                    public List<Bound_Node.Boundary> B_List;
                    protected Bound_Node.Boundary Flags;
                    protected int[] id;

                    public Bound_Node_RDD(Point loc, double rho0, double dt, double dx, double C, int[] id_in, List<Bound_Node.Boundary> B_in)
                    : base(loc)
                    {
                        id = id_in;
                        B_List = B_in;
                        Flags = Bound_Node.Boundary.None;
                        for (int i = 0; i < B_in.Count; i++) Flags |= B_in[i];
                    }

                    public void Complete_Boundary()
                    {
                        /*
                        [0] = x+y+z+
                        [1] = x+y-z-
                        [2] = x+y-z+
                        [3] = x+y+z-
                        */
                        foreach (Bound_Node.Boundary b in B_List)
                        {
                            if (b == Bound_Node.Boundary.DXPosYPosZPos)
                            {
                                if (Links2[0] is Null_Node) continue;
                                (Links2[0] as RDD_Node).Links2[5] = new Null_Node();
                                Links2[0] = new Null_Node();
                                //(Links2[0] as RDD_Node).Links2[5] = (Links2[0] as RDD_Node).Links2[0];
                                //Links2[0] = Links2[5];
                            }
                            else if (b == Bound_Node.Boundary.DXPosYNegZNeg)
                            {
                                if (Links2[1] is Null_Node) continue;
                                (Links2[1] as RDD_Node).Links2[4] = new Null_Node();
                                Links2[1] = new Null_Node();
                                //(Links2[1] as RDD_Node).Links2[4] = (Links2[1] as RDD_Node).Links2[1];
                                //Links2[1] = Links2[4];
                            }
                            else if (b == Bound_Node.Boundary.DXPosYNegZPos)
                            {
                                if (Links2[2] is Null_Node) continue;
                                (Links2[2] as RDD_Node).Links2[7] = new Null_Node();
                                Links2[2] = new Null_Node();
                                //(Links2[2] as RDD_Node).Links2[7] = (Links2[2] as RDD_Node).Links2[2];
                                //Links2[2] = Links2[7];
                            }
                            else if (b == Bound_Node.Boundary.DXPosYPosZNeg)
                            {
                                if (Links2[3] is Null_Node) continue;
                                (Links2[3] as RDD_Node).Links2[6] = new Null_Node();
                                Links2[3] = new Null_Node();
                                //(Links2[3] as RDD_Node).Links2[6] = (Links2[3] as RDD_Node).Links2[3];
                                //Links2[3] = Links2[6];
                            }
                            /*
                            [4] = x-y+z+
                            [5] = x-y-z-
                            [6] = x-y-z+
                            [7] = x-y+z-
                            */
                            else if (b == Bound_Node.Boundary.DXNegYPosZPos)
                            {
                                if (Links2[4] is Null_Node) continue;
                                (Links2[4] as RDD_Node).Links2[1] = new Null_Node();
                                Links2[4] = new Null_Node();
                                //(Links2[4] as RDD_Node).Links2[1] = (Links2[4] as RDD_Node).Links2[4];
                                //Links2[4] = Links2[1];
                            }
                            else if (b == Bound_Node.Boundary.DXNegYNegZNeg)
                            {
                                if (Links2[5] is Null_Node) continue;
                                (Links2[5] as RDD_Node).Links2[0] = new Null_Node();
                                Links2[5] = new Null_Node();
                                //(Links2[5] as RDD_Node).Links2[0] = (Links2[5] as RDD_Node).Links2[5];
                                //Links2[5] = Links2[0];
                            }
                            else if (b == Bound_Node.Boundary.DXNegYNegZPos)
                            {
                                if (Links2[6] is Null_Node) continue;
                                (Links2[6] as RDD_Node).Links2[3] = new Null_Node();
                                Links2[6] = new Null_Node();
                                //(Links2[6] as RDD_Node).Links2[3] = (Links2[6] as RDD_Node).Links2[6];
                                //Links2[6] = Links2[3];
                            }
                            else if (b == Bound_Node.Boundary.DXNegYPosZNeg)
                            {
                                if (Links2[7] is Null_Node) continue;
                                (Links2[7] as RDD_Node).Links2[2] = new Null_Node();
                                Links2[7] = new Null_Node();
                                //(Links2[7] as RDD_Node).Links2[2] = (Links2[7] as RDD_Node).Links2[7];
                                //Links2[7] = Links2[2];
                            }
                            /*
                            [8] = y+
                            [9] = y-    
                            [10] = z+
                            [11] = z-
                            */
                            else if (b == Bound_Node.Boundary.AYPos)
                            {
                                if (Links2[8] is Null_Node) continue;
                                (Links2[8] as RDD_Node).Links2[9] = new Null_Node();
                                Links2[8] = new Null_Node();
                                //    (Links2[8] as RDD_Node).Links2[9] = (Links2[8] as RDD_Node).Links2[8];
                                //    Links2[8] = Links2[9];
                            }
                            else if (b == Bound_Node.Boundary.AYNeg)
                            {
                                if (Links2[9] is Null_Node) continue;
                                (Links2[9] as RDD_Node).Links2[8] = new Null_Node();
                                Links2[9] = new Null_Node();
                                //(Links2[9] as RDD_Node).Links2[8] = (Links2[9] as RDD_Node).Links2[9];
                                //Links2[9] = Links2[8];
                            }
                            else if (b == Bound_Node.Boundary.AZPos)
                            {
                                if (Links2[10] is Null_Node) continue;
                                (Links2[10] as RDD_Node).Links2[11] = new Null_Node();
                                Links2[10] = new Null_Node();
                                //(Links2[10] as RDD_Node).Links2[11] = (Links2[10] as RDD_Node).Links2[10];
                                //Links2[10] = Links2[11];
                            }
                            else if (b == Bound_Node.Boundary.AZNeg)
                            {
                                if (Links2[11] is Null_Node) continue;
                                (Links2[11] as RDD_Node).Links2[10] = new Null_Node();
                                Links2[11] = new Null_Node();
                                //(Links2[11] as RDD_Node).Links2[10] = (Links2[11] as RDD_Node).Links2[11];
                                //Links2[11] = Links2[10];
                            }
                        }
                        /*
                        [0] = x+y+z+
                        [1] = x+y-z-
                        [2] = x+y-z+
                        [3] = x+y+z-
                        [4] = x-y+z+
                        [5] = x-y-z-
                        [6] = x-y-z+
                        [7] = x-y+z-
                        [8] = y+
                        [9] = y-
                        [10] = z+
                        [11] = z-
                        */
                    }
                }

                public class Bound_Node_RDD_MaterialFilter : Bound_Node_RDD
                {
                    private readonly IIR_DIF[] filters;
                    private readonly double denom;       // 1 / (1 + kappaTotal)
                    private readonly double kappaMinus1; // kappaTotal - 1
                    private int activeNeighborCount;     // number of non-null neighbors

                    public Bound_Node_RDD_MaterialFilter(Point loc, double rho0, double dt, double dx, double C, int[] id_in, List<Environment.Material> materials, List<Bound_Node.Boundary> B_in)
                    :base(loc, rho0, dt, dx, C, id_in, B_in)
                    {
                        const double M = 0.25;

                        double lambda = (dx > 0.0) ? (C * dt / dx) : 1.0;
                        if (lambda <= 0.0) lambda = 1.0;

                        filters = new IIR_DIF[materials.Count];

                        double ksum = 0.0;
                        double fs = 1.0 / dt;

                        for (int i = 0; i < filters.Length; i++)
                        {
                            filters[i] = new DIF_IWB_2p(materials[i], fs, M, lambda);
                            ksum += filters[i].a_b;
                        }

                        //ksum = Math.Min(ksum, 0.95);

                        denom = 1.0 / (1.0 + ksum);
                        kappaMinus1 = ksum - 1.0;
                    }

                    public new void Complete_Boundary()
                    {
                        base.Complete_Boundary();
                    }

                    public override void UpdateP()
                    {
                        double p2 = 0.0;
                        foreach (Node node in Links2)
                        {
                            if (node is Null_Node)
                            {
                                p2 += Pn;
                            }
                            else
                            {
                                p2 += node.P;
                            }
                        }

                        // Memory term from filters
                        double mem = 0.0;
                        for (int i = 0; i < filters.Length; i++) mem += filters[i].g_b_term();

                        // Boundary update:
                        // Pnf = (0.25*p2 - Pn + (κ - 1)*Pn_1 + mem) / (1 + κ)
                        Pnf = (0.25 * p2 - Pn + kappaMinus1 * Pn_1 + mem) * denom;

                        // Advance filter states
                        for (int i = 0; i < filters.Length; i++) filters[i].Update(Pnf, Pn_1);
                    }

                    public override void UpdateT()
                    {
                        base.UpdateT();
                    }
                }

                public abstract class IIR_DIF
                {
                    public double a_b { get; protected set; }

                    public abstract double g_b_term();                  // memory term using g^n
                    public abstract void Update(double Pnf, double Pn_1);// advance states using (Pnf - Pn_1)
                    public abstract void SetCourant(double lambda);     // allow updating if needed
                }

                public sealed class DIF_IWB_2p : IIR_DIF
                {
                    private readonly double[] b;   // admittance numerator (normalized so a[0]=1)
                    private readonly double[] a;   // admittance denominator
                    private readonly double[] s;   // DF-II-T state registers
                    private readonly double M;     // boundary mask/weight (e.g. 0.25 for CCP/RDD 12-pt face)
                    private double lambda;         // Courant number
                    private readonly double b0;    // b[0] cached

                    public DIF_IWB_2p(Environment.Material mat, double sampleFrequency, double weight, double courant)
                    {
                        M = weight;
                        lambda = (courant > 0.0) ? courant : 1.0;

                        double[] f_axis;
                        double[] aY, bY;

                        (aY, bY) = mat.Estimate_IIR_Coefficients(sampleFrequency, sampleFrequency * Math.Sqrt(2) / 10, out f_axis);

                        aY = aY ?? new double[0];
                        bY = bY ?? new double[0];
                        int order = Math.Max(2, Math.Max(aY.Length, bY.Length));
                        if (aY.Length < order) Array.Resize(ref aY, order);
                        if (bY.Length < order) Array.Resize(ref bY, order);

                        b = (double[])bY.Clone();
                        a = (double[])aY.Clone();

                        // Invalid leading coefficients are rejected by EnsurePassivity below.
                        double a0 = a[0];
                        for (int i = 0; i < order; i++)
                        {
                            b[i] /= a0;
                            a[i] /= a0;
                        }

                        EnsurePassivity(b, a, mat, sampleFrequency);

                        b0 = b[0];

                        a_b = lambda * b0 * M;

                        s = new double[order-1]; // initialized to 0
                    }

                    private static bool IsStable(double[] denominator)
                    {
                        // Schur step-down tests z^N + a[1] z^(N-1) + ... + a[N].
                        // Unlike numerator scaling, this rejects poles on/outside the unit circle.
                        double[] polynomial = (double[])denominator.Clone();
                        for (int n = polynomial.Length - 1; n > 0; n--)
                        {
                            double reflection = polynomial[n] / polynomial[0];
                            if (double.IsNaN(reflection) || Math.Abs(reflection) >= 1.0 - 1e-10) return false;
                            double scale = 1.0 - reflection * reflection;
                            double[] reduced = new double[n];
                            for (int k = 0; k < n; k++)
                                reduced[k] = (polynomial[k] - reflection * polynomial[n - k]) / scale;
                            polynomial = reduced;
                        }
                        return true;
                    }

                    private static void EnsurePassivity(double[] b, double[] a, Environment.Material mat, double fs)
                    {
                        bool valid = a[0] == 1.0 && b[0] >= 0.0;
                        for (int k = 0; k < a.Length; k++)
                            valid &= !double.IsNaN(a[k]) && !double.IsInfinity(a[k]) &&
                                     !double.IsNaN(b[k]) && !double.IsInfinity(b[k]);
                        valid = valid && IsStable(a);

                        if (valid)
                        {
                            // Re(Y) has the sign of Re(B * conjugate(A)). Represent the latter
                            // as a cosine polynomial; bound interpolation error between samples
                            // using |q''| <= sum(k^2 * |c[k]|), rather than trusting samples alone.
                            double[] cosine = new double[a.Length];
                            for (int j = 0; j < b.Length; j++)
                                for (int k = 0; k < a.Length; k++)
                                    cosine[Math.Abs(j - k)] += b[j] * a[k];
                            double curvature = 0.0, magnitude = 0.0;
                            for (int k = 0; k < cosine.Length; k++)
                            {
                                curvature += k * k * Math.Abs(cosine[k]);
                                magnitude += Math.Abs(cosine[k]);
                            }
                            const int nFreq = 4096;
                            double minimum = double.PositiveInfinity;
                            for (int i = 0; i <= nFreq; i++)
                            {
                                double w = Math.PI * i / nFreq, q = 0.0;
                                for (int k = 0; k < cosine.Length; k++) q += cosine[k] * Math.Cos(k * w);
                                minimum = Math.Min(minimum, q);
                            }
                            double errorBound = curvature * Math.Pow(Math.PI / nFreq, 2) / 8.0;
                            valid = !double.IsNaN(minimum) && !double.IsInfinity(magnitude) &&
                                    minimum >= errorBound + 1e-12 * magnitude;
                            // Exactly rigid (zero numerator) is passive, with stable state poles.
                            if (magnitude == 0.0) valid = true;
                        }
                        if (valid) return;

                        // Conservative first-order low-pass based on normal-incidence absorption.
                        // Keep an actual stable pole and state register, including for rigid walls.
                        double gain = 0.0;
                        var reflectionMagnitude = mat.Reflection_Narrow(fs * Math.Sqrt(2) / 20.0).Magnitude;
                        if (!double.IsNaN(reflectionMagnitude) && !double.IsInfinity(reflectionMagnitude))
                        {
                            reflectionMagnitude = Math.Min(1.0, reflectionMagnitude);
                            gain = (1.0 - reflectionMagnitude) / (1.0 + reflectionMagnitude);
                        }
                        Array.Clear(a, 0, a.Length);
                        Array.Clear(b, 0, b.Length);
                        a[0] = 1.0;
                        a[1] = -0.5;
                        b[0] = 0.5 * gain;
                        System.Diagnostics.Trace.TraceWarning("FVM material admittance failed stability/passivity validation; using a conservative first-order boundary filter.");
                    }

                    public override void SetCourant(double newLambda)
                    {
                        lambda = (newLambda > 0.0) ? newLambda : 1.0;
                        a_b = lambda * b0 * M;
                    }

                    public override double g_b_term()
                    {
                        // States already include lambda through x; admittance opposes the pressure change.
                        return s.Length == 0 ? 0.0 : -M * s[0];
                    }

                    public override void Update(double Pnf, double Pn_1)
                    {
                        if (s.Length == 0) return;
                        double x = lambda * (Pnf - Pn_1);
                        double y = b0 * x + s[0];
                        int N = s.Length;
                        for (int i = 0; i < N - 1; i++)
                        {
                            double bi = (i + 1 < b.Length) ? b[i + 1] : 0.0;
                            double ai = (i + 1 < a.Length) ? a[i + 1] : 0.0;
                            s[i] = s[i + 1] + (bi * x) - (ai * y);
                        }

                        int last = N - 1;
                        double bLast = (last + 1 < b.Length) ? b[last + 1] : 0.0;
                        double aLast = (last + 1 < a.Length) ? a[last + 1] : 0.0;
                        s[last] = (bLast * x) - (aLast * y);

                        if (double.IsNaN(s[0]) || double.IsInfinity(s[0]) || Math.Abs(s[0]) > 1e12)
                            Array.Clear(s, 0, s.Length);
                    }
                }
            }
        }
    }
}