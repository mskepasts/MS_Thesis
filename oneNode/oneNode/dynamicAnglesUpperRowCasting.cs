using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Rhino;
using Rhino.Geometry;

namespace oneNode
{
    public class Casting
    {
        public double TopAngle1 { get; set; }
        public double TopAngle2 { get; set; }
        public double SideAngle { get; set; }
        public Casting(double topAngle1, double topAngle2, double sideAngle)
        {
            TopAngle1 = topAngle1;
            TopAngle2 = topAngle2;
            SideAngle = sideAngle;
        }
        // Other properties, methods, events...
    }
    public class dynamicAnglesUpperRowCastings : GH_Component
    {
        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public dynamicAnglesUpperRowCastings()
            : base("dynamic angles upper Row Castings", "DA upper Row C",
            "Description",
            "Category", "Subcategory")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddLineParameter("TopLines", "TL", "Sorted top lines of the row", GH_ParamAccess.list);
            pManager.AddLineParameter("DiagonalLines", "DL", "Sorted bottom lines of the row", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Top Line Index", "TL_index", "Starting index of the top lines", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Diagonal Index", "DL_index", "Starting index of the diagonal lines", GH_ParamAccess.item);
            pManager.AddNumberParameter("Top Angles0", "TA0s", "All the top angles used to test", GH_ParamAccess.list);
            pManager.AddNumberParameter("Top Angles4", "TA4s", "All the top angles used to test", GH_ParamAccess.list);
            pManager.AddNumberParameter("Top Angles5", "TA5s", "All the top angles used to test", GH_ParamAccess.list);
            pManager.AddNumberParameter("Top Angles9", "TA9s", "All the top angles used to test", GH_ParamAccess.list);
            pManager.AddNumberParameter("Def Angles1", "DA1s", "All the top angles used to test", GH_ParamAccess.list);
            pManager.AddNumberParameter("Def Angles2", "DA2s", "All the top angles used to test", GH_ParamAccess.list);
            pManager.AddNumberParameter("Side Angles5", "SA5s", "All the side angles used to test", GH_ParamAccess.list);
            pManager.AddMeshParameter("Dome mesh", "DM", "Dome mesh", GH_ParamAccess.item);
            pManager.AddNumberParameter("Castings", "CM", "Castings", GH_ParamAccess.tree);
            pManager.AddIntegerParameter("Casting from 2 becoming 6", "C6", "Dome mesh", GH_ParamAccess.list, new List<int>());
            pManager.AddIntegerParameter("Casting from 1 becoming 7", "C7", "C7", GH_ParamAccess.list, new List<int>());
            pManager.AddIntegerParameter("Machine Number", "MN", "MN", GH_ParamAccess.item, 0);
            pManager.AddIntegerParameter("MoldCastings", "MCM", "Mold Castings", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.Register_LineParam("Original Lines", "oLines", "The lines coming into the component");
            pManager.Register_PointParam("Changed Points", "pOut", "The new points changed based on the angles");
            pManager.Register_LineParam("Changed Lines", "newLines", "The lines after angle has been changed");
            pManager.Register_GenericParam("AnglesTop", "aOutT", "The angles of the node going from top corner counterclockwise");
            pManager.Register_GenericParam("AnglesSide", "aOutS", "The angles of the node going from top corner counterclockwise");
            pManager.Register_GenericParam("Testing point", "pOutTest", "The new points changed based on the angles");
            pManager.Register_DoubleParam("Cost", "Cost", "Cost in order to change the node to fit the required angles");
            pManager.Register_DoubleParam("AngleCost", "ACost", "Cost in order to change the node to fit the required angles");
            pManager.Register_PlaneParam("side planes", "testSP", "sideplane test");
            pManager.Register_PlaneParam("top planes", "testTP", "topplane test");
            pManager.Register_MeshParam("New meshes", "testMesh", "mesh test");
            pManager.Register_DoubleParam("castings*dist", "CND", "casting Numbers * dist");
            pManager.Register_GenericParam("Castings", "C", "castings");
            pManager.Register_GenericParam("Machined Castings", "MC", "Machined Castings");
            pManager.Register_IntegerParam("Casting Build Up", "CBD", "casting Numbers Machined");
            pManager.Register_IntegerParam("Molded Castings", "MdC", "Molded Castings");
        }

        // Sort the lines based on the lines inputed
        public List<Point3d> GetEvalPoints(Line topLine, Line diagonalOne, Line diagonalTwo)
        {
            List<Point3d> points = new List<Point3d>();
            if (topLine.From == diagonalTwo.To)
            {
                points.Add(topLine.From);
                points.Add(diagonalTwo.From);
                points.Add(topLine.To);
            }
            else if (topLine.From == diagonalTwo.From)
            {
                points.Add(topLine.From);
                points.Add(diagonalTwo.To);
                points.Add(topLine.To);
            }
            else if (topLine.To == diagonalTwo.From)
            {
                points.Add(topLine.To);
                points.Add(diagonalTwo.To);
                points.Add(topLine.From);
            }
            else
            {
                points.Add(topLine.To);
                points.Add(diagonalTwo.From);
                points.Add(topLine.From);
            }

            if (points.Contains(diagonalOne.From))
                points.Add(diagonalOne.To);
            else
                points.Add(diagonalOne.From);

            return points;
        }

        // Change out point from last angle being found
        public List<Point3d> GetAlteredPoints(Point3d pNew, List<Point3d> p)
        {
            p[2] = pNew;
            return p;
        }

        // Find side plane angle
        public double FindingSideAngle(Point3d p, Point3d pOther, Plane basePlane)
        {
            Point3d holderPt = p;
            var xformPt = Rhino.Geometry.Transform.Translation(basePlane.ZAxis);
            holderPt.Transform(xformPt);
            Plane planeSide = new Plane();
            Plane.FitPlaneToPoints(new List<Point3d> { p, basePlane.ClosestPoint(pOther), holderPt }, out planeSide);
            planeSide.Origin = p;
            double angle = Vector3d.VectorAngle(pOther - p, basePlane.ClosestPoint(pOther) - p, planeSide);
            if (angle > Math.PI) { angle = 2 * Math.PI - angle; }
            return angle;
        }

        // Find side plane for debugging
        public Plane FindingSidePlane(Point3d p, Point3d pOther, Plane basePlane)
        {
            Point3d holderPt = p;
            var xformPt = Rhino.Geometry.Transform.Translation(basePlane.ZAxis);
            holderPt.Transform(xformPt);
            Plane planeSide = new Plane();
            Plane.FitPlaneToPoints(new List<Point3d> { p, basePlane.ClosestPoint(pOther), holderPt }, out planeSide);
            double angle = Vector3d.VectorAngle(pOther - p, basePlane.ClosestPoint(pOther) - p, planeSide);
            planeSide.Origin = p;
            if (angle > Math.PI) { angle = 2 * Math.PI - angle; }
            return planeSide;
        }

        // Find angles based on 4 joint
        public Point3d FindXandYandZ(double theta2, double thetaS2, double theta42, Plane plane2, Plane plane4, List<Point3d> p)
        {
            theta42 = RhinoMath.ToRadians(theta42);
            theta2 = RhinoMath.ToRadians(theta2);
            thetaS2 = RhinoMath.ToRadians(thetaS2);


            Line holderLn = new Line(p[1], plane2.ClosestPoint(p[3]));
            Plane planeS = new Plane();
            Plane.FitPlaneToPoints(new List<Point3d> { p[1], p[3], plane2.ClosestPoint(p[3]) }, out planeS);
            planeS.Origin = plane2.Origin;
            var xform2 = Rhino.Geometry.Transform.Rotation((thetaS2), planeS.ZAxis, planeS.Origin);
            var xform = Rhino.Geometry.Transform.Rotation((-theta2), plane2.ZAxis, plane2.Origin);
            xform *= xform2;
            holderLn.Transform(xform);


            Line holderLn4 = new Line(p[3], plane4.ClosestPoint(p[1]));
            var xform4 = Rhino.Geometry.Transform.Rotation((theta42), plane4.ZAxis, plane4.Origin);
            holderLn4.Transform(xform4);
            Plane planeInt = new Plane();
            Point3d holderPt = p[3];
            var xformPt = Rhino.Geometry.Transform.Translation(plane4.ZAxis);
            holderPt.Transform(xformPt);
            Plane.FitPlaneToPoints(new List<Point3d> { p[3], holderLn4.To, holderPt }, out planeInt);
            planeInt.Origin = plane4.Origin;
            BoundingBox box = new BoundingBox(new List<Point3d> { p[3], holderLn4.To, holderPt });
            holderLn.ExtendThroughBox(box);


            double t = 0.5;
            Rhino.Geometry.Intersect.Intersection.LinePlane(holderLn, planeInt, out t);



            return holderLn.PointAt(t);
        }

        // Find angles top (theta)
        public List<double> FindAnglesTop(List<Point3d> pts, List<Plane> plns)
        {
            // Node four
            double theta42 = Vector3d.VectorAngle(pts[3] - pts[1], pts[3] - pts[0], plns[3]);
            double theta43 = Vector3d.VectorAngle(pts[3] - pts[0], pts[3] - pts[2], plns[3]);

            // Node three
            double theta3 = Vector3d.VectorAngle(pts[2] - pts[3], pts[2] - pts[0], plns[2]);

            // Node two
            double theta2 = Vector3d.VectorAngle(pts[1] - pts[0], pts[1] - pts[3], plns[1]);

            // Node one
            double theta12 = Vector3d.VectorAngle(pts[0] - pts[3], pts[0] - pts[1], plns[0]);
            double theta13 = Vector3d.VectorAngle(pts[0] - pts[2], pts[0] - pts[3], plns[0]);

            List<double> angles = new List<double>();
            angles.Add(Math.Round(RhinoMath.ToDegrees(theta3), 2));
            angles.Add(Math.Round(RhinoMath.ToDegrees(theta13), 2));
            angles.Add(Math.Round(RhinoMath.ToDegrees(theta12), 2));
            angles.Add(Math.Round(RhinoMath.ToDegrees(theta2), 2));
            angles.Add(Math.Round(RhinoMath.ToDegrees(theta42), 2));
            angles.Add(Math.Round(RhinoMath.ToDegrees(theta43), 2));




            return angles;

        }

        // Find the side planes (debugging)
        public List<Plane> FindPlanesSide(List<Point3d> pts, List<Plane> plns)
        {

            List<Plane> planes = new List<Plane>();



            planes.Add((FindingSidePlane(pts[2], pts[3], plns[2])));
            planes.Add((FindingSidePlane(pts[2], pts[0], plns[2])));
            planes.Add((FindingSidePlane(pts[0], pts[2], plns[0])));
            planes.Add((FindingSidePlane(pts[0], pts[3], plns[0])));
            planes.Add((FindingSidePlane(pts[0], pts[1], plns[0])));
            planes.Add((FindingSidePlane(pts[1], pts[0], plns[1])));
            planes.Add((FindingSidePlane(pts[1], pts[3], plns[1])));
            planes.Add((FindingSidePlane(pts[3], pts[1], plns[3])));
            planes.Add((FindingSidePlane(pts[3], pts[0], plns[3])));
            planes.Add((FindingSidePlane(pts[3], pts[2], plns[3])));

            return planes;

        }

        // Find angles side (sigma)
        public List<double> FindAnglesSide(List<Point3d> pts, List<Plane> plns)
        {

            List<double> angles = new List<double>();



            angles.Add(Math.Round(RhinoMath.ToDegrees(FindingSideAngle(pts[2], pts[3], plns[2])), 2));
            angles.Add(Math.Round(RhinoMath.ToDegrees(FindingSideAngle(pts[2], pts[0], plns[2])), 2));
            angles.Add(Math.Round(RhinoMath.ToDegrees(FindingSideAngle(pts[0], pts[2], plns[0])), 2));
            angles.Add(Math.Round(RhinoMath.ToDegrees(FindingSideAngle(pts[0], pts[3], plns[0])), 2));
            angles.Add(Math.Round(RhinoMath.ToDegrees(FindingSideAngle(pts[0], pts[1], plns[0])), 2));
            angles.Add(Math.Round(RhinoMath.ToDegrees(FindingSideAngle(pts[1], pts[0], plns[1])), 2));
            angles.Add(Math.Round(RhinoMath.ToDegrees(FindingSideAngle(pts[1], pts[3], plns[1])), 2));
            angles.Add(Math.Round(RhinoMath.ToDegrees(FindingSideAngle(pts[3], pts[1], plns[3])), 2));
            angles.Add(Math.Round(RhinoMath.ToDegrees(FindingSideAngle(pts[3], pts[0], plns[3])), 2));
            angles.Add(Math.Round(RhinoMath.ToDegrees(FindingSideAngle(pts[3], pts[2], plns[3])), 2));

            return angles;

        }

        // Find x and y based on the angles given




        public double CheckAngleTolerance(double val)
        {
            double factor = 10.0;
            if (val / (1.0 / factor) < 20)
                val = val % (1.0 / factor);
            if (val < 0.01)
                val = 0.0;
            return val;
        }




        // Finding the cost of the function


        // Finding the top planes
        public List<Plane> GetTopPlanes(List<Point3d> pts, Mesh domeMesh)
        {
            List<Plane> planes = new List<Plane>();
            foreach (Point3d pt in pts)
            {
                if (domeMesh.ClosestMeshPoint(pt, 4) != null)
                {
                    Vector3d normal = domeMesh.NormalAt(domeMesh.ClosestMeshPoint(pt, 4));
                    planes.Add(new Plane(pt, normal));
                }
            }

            return planes;
        }


        // adds the casting to the list
        public Boolean AddCasting(double angle, double sideAngle, int num, List<Casting> castings)
        {
            if (num == 1)
                castings.Add(new Casting(angle, 60 - angle, sideAngle));
            else
                castings.Add(new Casting(60 - angle, angle, sideAngle));
            return true;
        }

        public double GetTopAngle(List<double> topAngles1, List<double> topAngles2, int casting1, int casting2)
        {
            return topAngles1[casting1] + topAngles2[casting2];
        }

        public Mesh FindNewMesh(Point3d pNew, List<Point3d> p, Mesh grid)
        {
            Mesh newGrid = new Mesh();
            var nodes = new Grasshopper.Kernel.Geometry.Node2List();

            List<Point3d> pts = new List<Point3d>();
            foreach (Point3d point in grid.Vertices.ToPoint3dArray())
            {
                if (point.Equals(p[0]))
                {
                    nodes.Append(new Grasshopper.Kernel.Geometry.Node2(pNew.X, pNew.Y));
                    pts.Add(pNew);
                }
                else
                {
                    nodes.Append(new Grasshopper.Kernel.Geometry.Node2(point.X, point.Y));
                    pts.Add((Point3d)point);
                }
            }
            //solve Delaunay
            var delMesh = new Mesh();
            var faces = new List<Grasshopper.Kernel.Geometry.Delaunay.Face>();

            faces = Grasshopper.Kernel.Geometry.Delaunay.Solver.Solve_Faces(nodes, 0);

            //output
            delMesh = Grasshopper.Kernel.Geometry.Delaunay.Solver.Solve_Mesh(nodes, 0, ref faces);
            for (int i = 0; i < pts.Count; i++)
            {
                delMesh.Vertices.SetVertex(i, pts[i]);
            }
            newGrid = delMesh;

            newGrid.RebuildNormals();


            return newGrid;
        }

        public List<List<double>> CastingsToNestedList(List<Casting> castings)
        {
            List<List<double>> result = new List<List<double>>();
            foreach (Casting casting in castings)
                result.Add(new List<double> { casting.TopAngle1, casting.TopAngle2, casting.SideAngle });
            return result;
        }
        public DataTree<double> CastingsToTree(List<Casting> castings)
        {
            DataTree<double> result = new DataTree<double>();
            int range = 0;
            foreach (Casting casting in castings)
            {
                GH_Path pth = new GH_Path(range);
                result.AddRange(new List<double> { casting.TopAngle1, casting.TopAngle2, casting.SideAngle }, pth);
                range++;
            }
            return result;
        }

        public Boolean CheckCastings(double ta, List<double> sa, List<Casting> castings, List<int> moldInd, List<int> castingIndex, List<int> currentIndex, List<Casting> machinedCastings, double reqAngle = -1, int reqNum = -1)
        {
            double tolerance = 0.03;
            double machining = 1.0;
            if (reqAngle == -1)
            {
                for (int i = 0; i < castings.Count; i++)
                    for (int j = 0; j < castings.Count; j++)
                    {
                        Casting casting1 = castings[i];
                        Casting casting2 = castings[j];
                        if (Math.Abs(casting1.SideAngle - sa[0]) < tolerance && Math.Abs(casting2.SideAngle - sa[1]) < tolerance
                            && Math.Abs(casting2.TopAngle2 - Math.Abs(casting1.TopAngle1 - ta)) < 2 * tolerance
                            && Math.Abs(casting1.TopAngle1 - Math.Abs(casting2.TopAngle2 - ta)) < 2 * tolerance)
                        {
                            castingIndex[currentIndex[0]] = i;
                            castingIndex[currentIndex[1]] = j;
                            return true;
                        }

                        if (Math.Abs(casting1.SideAngle - sa[0]) < tolerance && Math.Abs(casting2.SideAngle - sa[1]) < machining
                            && Math.Abs(casting2.TopAngle2 - Math.Abs(casting1.TopAngle1 - ta)) < machining
                            && Math.Abs(casting1.TopAngle1 - Math.Abs(casting2.TopAngle2 - ta)) < machining)
                        {
                            if (moldInd.Contains(i))
                            {
                                castingIndex[currentIndex[0]] = i;
                                castingIndex[currentIndex[1]] = castings.Count;
                                AddCasting(Math.Abs(casting1.TopAngle1 - ta), sa[1], 2, castings);
                                AddCasting(Math.Abs(casting1.TopAngle1 - ta), sa[1], 2, machinedCastings);
                                return true;
                            }
                        }
                        if (Math.Abs(casting1.SideAngle - sa[0]) < machining && Math.Abs(casting2.SideAngle - sa[1]) < tolerance
                            && Math.Abs(casting2.TopAngle2 - Math.Abs(casting1.TopAngle1 - ta)) < machining
                            && Math.Abs(casting1.TopAngle1 - Math.Abs(casting2.TopAngle2 - ta)) < machining)
                        {
                            if (moldInd.Contains(j))
                            {
                                castingIndex[currentIndex[1]] = j;
                                castingIndex[currentIndex[0]] = castings.Count;
                                AddCasting(Math.Abs(casting2.TopAngle2 - ta), sa[0], 1, castings);
                                AddCasting(Math.Abs(casting2.TopAngle2 - ta), sa[0], 1, machinedCastings);
                                return true;
                            }
                        }


                    }
                return false;
            }
            else
            {
                if (reqNum == 1)
                {
                    for (int i = 0; i < castings.Count; i++)
                    {
                        Casting casting = castings[i];
                        if (Math.Abs(casting.SideAngle - sa[0]) < tolerance
                                && Math.Abs(casting.TopAngle2 - Math.Abs(reqAngle - ta)) < 2 * tolerance)
                        {
                            castingIndex[currentIndex[0]] = i;
                            return true;
                        }
                        if (Math.Abs(casting.SideAngle - sa[0]) < machining
                               && Math.Abs(casting.TopAngle2 - Math.Abs(reqAngle - ta)) < 2 * machining)
                        {
                            if (moldInd.Contains(i))
                            {
                                castingIndex[currentIndex[0]] = castings.Count;
                                AddCasting(Math.Abs(reqAngle - ta), sa[0], 2, castings);
                                AddCasting(Math.Abs(reqAngle - ta), sa[0], 2, machinedCastings);
                                return true;
                            }
                        }
                    }
                    castingIndex[currentIndex[0]] = castings.Count;
                    moldInd.Add(castings.Count);
                    AddCasting(Math.Abs(reqAngle - ta), sa[0], 2, castings);
                   
                    return true;
                }

                else
                {
                    for (int i = 0; i < castings.Count; i++)
                    {
                        Casting casting = castings[i];
                        if (Math.Abs(casting.SideAngle - sa[0]) < tolerance
                                && Math.Abs(casting.TopAngle1 - Math.Abs(reqAngle - ta)) < 2 * tolerance)

                        {
                            castingIndex[currentIndex[0]] = i;
                            return true;
                        }
                        if (Math.Abs(casting.SideAngle - sa[0]) < machining
                               && Math.Abs(casting.TopAngle1 - Math.Abs(reqAngle - ta)) < 2 * machining)
                        {
                            
                            if (moldInd.Contains(i))
                            {
                                castingIndex[currentIndex[0]] = castings.Count;
                                AddCasting(Math.Abs(reqAngle - ta), sa[0], 1, castings);
                                AddCasting(Math.Abs(reqAngle - ta), sa[0], 1, machinedCastings);
                                return true;
                            }
                        }
                    }
                    castingIndex[currentIndex[0]] = castings.Count;
                    moldInd.Add(castings.Count);
                    AddCasting(Math.Abs(reqAngle - ta), sa[0], 1, castings);
         
                    return true;
                }
            }
        }



        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<Casting> castings = new List<Casting>();
            List<Casting> machinedCastings = new List<Casting>();
           
            List<double> topAngles0 = new List<double>();
            List<double> topAngles4 = new List<double>();
            List<double> topAngles5 = new List<double>();
            List<double> topAngles9 = new List<double>();
            List<double> sideAngles0 = new List<double>();
            List<double> sideAngles4 = new List<double>();
            List<double> sideAngles5 = new List<double>();
            List<double> sideAngles9 = new List<double>();
            List<double> defAngles1 = new List<double>();
            List<double> defAngles2 = new List<double>();
            List<int> castingIndex_7 = new List<int>();
            List<int> castingIndex_6 = new List<int>();
            Boolean failed = false;
            // Object test = null;
            // List<List<Object>> test = new List<List<Object>>();

            GH_Structure<Grasshopper.Kernel.Types.GH_Number> test = new GH_Structure<Grasshopper.Kernel.Types.GH_Number>();

            List<int> moldInd = new List<int>();
            int MN = 0;
            List<Line> topLines = new List<Line>();
            List<Line> diagonalLines = new List<Line>();

            int topLineIndexStart = 0;
            int diagonalLineIndexStart = 0;

            Mesh domeMesh = new Mesh();

            if (!DA.GetDataList(0, topLines)) { return; }
            if (!DA.GetDataList(1, diagonalLines)) { return; }
            if (!DA.GetData(2, ref topLineIndexStart)) { return; }
            if (!DA.GetData(3, ref diagonalLineIndexStart)) { return; }
            if (!DA.GetDataList(4, topAngles0)) { return; }
            if (!DA.GetDataList(5, topAngles4)) { return; }
            if (!DA.GetDataList(6, topAngles5)) { return; }
            if (!DA.GetDataList(7, topAngles9)) { return; }
            if (!DA.GetDataList(8, defAngles1)) { return; }
            if (!DA.GetDataList(9, defAngles2)) { return; }

            if (!DA.GetDataList(10, sideAngles5)) { return; }

            if (!DA.GetData(11, ref domeMesh)) { return; }
            if (!DA.GetDataTree(12, out test)) { }

            if (test.Branches.Count > 1)

                foreach (List<Grasshopper.Kernel.Types.GH_Number> obj in test.Branches)
                {




                    double topAngle1 = obj[0].Value;
                    double topAngle2 = obj[1].Value;
                    double sideAngle = obj[2].Value;
                    castings.Add(new Casting(topAngle1, topAngle2, sideAngle));
                }
            if (!DA.GetDataList(16, moldInd)) { }


            if (!DA.GetDataList(13, castingIndex_6)) { }
            if (!DA.GetDataList(14, castingIndex_7)) { }
            if (!DA.GetData(15, ref MN)) { }

            int tl_current = topLineIndexStart;
            int dl_current = diagonalLineIndexStart;
            int branch_current = 0;

            Point3d pLast = new Point3d();
            DataTree<Point3d> newP = new DataTree<Point3d>();
            DataTree<double> costs = new DataTree<double>();
            DataTree<double> angleCosts = new DataTree<double>();
            DataTree<Point3d> newP1 = new DataTree<Point3d>();
            DataTree<Line> originalLines = new DataTree<Line>();
            DataTree<Line> newLines = new DataTree<Line>();
            DataTree<double> topAnglesFinal = new DataTree<double>();
            DataTree<double> sideAnglesFinal = new DataTree<double>();
            DataTree<Plane> sidePlanes = new DataTree<Plane>();
            DataTree<Plane> topPlanes = new DataTree<Plane>();
            DataTree<int> castingPositions = new DataTree<int>();
            DataTree<Mesh> meshes = new DataTree<Mesh>();

            Mesh finalMesh = domeMesh;
            double value1 = 2.0;
            Casting nextCasting_9 = null;
            Casting nextCasting_0 = null;
            int nextCastingIndex_0 = -1;
            int nextCastingIndex_9 = -1;

            double dist = double.PositiveInfinity;

            for (int i = topLineIndexStart; i < topLines.Count - 1; i++)
            {
                double topAngle0 = topAngles0[topAngles0.Count - 1];
                if (branch_current < topAngles0.Count)
                    topAngle0 = topAngles0[branch_current];

                double topAngle4 = topAngles4[topAngles4.Count - 1];
                if (branch_current < topAngles4.Count)
                    topAngle4 = topAngles4[branch_current];

                double topAngle5 = topAngles5[topAngles5.Count - 1];
                if (branch_current < topAngles5.Count)
                    topAngle5 = topAngles5[branch_current];

                double topAngle9 = topAngles9[topAngles9.Count - 1];
                if (branch_current < topAngles9.Count)
                    topAngle9 = topAngles9[branch_current];



                double sideAngle5 = sideAngles5[sideAngles5.Count - 1];
                if (branch_current < sideAngles5.Count)
                    sideAngle5 = sideAngles5[branch_current];





                double topAngle1 = defAngles1[defAngles1.Count - 1];
                if (branch_current < defAngles1.Count)
                    topAngle1 = defAngles1[branch_current];

                double topAngle2 = defAngles2[defAngles2.Count - 1];
                if (branch_current < defAngles2.Count)
                    topAngle2 = defAngles2[branch_current];

                List<int> castingIndex = new List<int>();
                for (int j = 0; j < 10; j++)
                    castingIndex.Add(-1);
                GH_Path pth = new GH_Path(branch_current);
                domeMesh = finalMesh;
                List<Point3d> p = GetEvalPoints(topLines[tl_current], diagonalLines[dl_current], diagonalLines[dl_current + 1]);
                List<Plane> planes = GetTopPlanes(p, domeMesh);
                List<Plane> finalPlanes = GetTopPlanes(p, domeMesh);
                List<Plane> holderPlanes = new List<Plane>();
                double angleCost = double.PositiveInfinity;
                double cost = double.PositiveInfinity;

                List<int> bestFit = new List<int>();

                if (tl_current != topLineIndexStart)
                    p = GetAlteredPoints(pLast, p);

                List<Point3d> holderPoints = new List<Point3d>();
                Point3d p1New_current = p[0];
                Mesh holderMesh = new Mesh();



                p1New_current = FindXandYandZ(topAngle1, sideAngle5, topAngle2, planes[1], planes[3], p);
                for (int q = 0; q < 4; q++)
                    if (q == 0)
                        holderPoints.Add(p1New_current);
                    else
                        holderPoints.Add(p[q]);
                finalPlanes = GetTopPlanes(p, FindNewMesh(p1New_current, p, domeMesh));
                List<double> tas = FindAnglesTop(holderPoints, finalPlanes);
                List<double> sas = FindAnglesSide(holderPoints, finalPlanes);
                castingIndex[0] = nextCastingIndex_0;
                castingIndex[9] = nextCastingIndex_9;


                if (p1New_current.DistanceTo(p[0]) < dist)
                    dist = p1New_current.DistanceTo(p[0]);


                if (p1New_current.DistanceTo(p[0]) < value1)
                {
                    if (castingIndex_6.Count -  1 < branch_current || castingIndex_6.Count < 2)
                    {
                        // angle 65
                        if (castings.Count == 0)

                        {
                            moldInd.Add(castings.Count);
                            AddCasting(60 - topAngle5, sideAngle5, 2, castings);
                            
                            castingIndex[5] = 0;
                            CheckCastings(tas[3], new List<double> { sas[6] }, castings, moldInd, castingIndex, new List<int> { 6 }, machinedCastings, castings[castingIndex[5]].TopAngle2, 2);
                        }
                        

                        else if (!CheckCastings(tas[3], new List<double> { sas[6], sas[5] }, castings, moldInd, castingIndex, new List<int> { 6, 5 }, machinedCastings))
                        {
                            castingIndex[5] = castings.Count;
                            moldInd.Add(castings.Count);
                            AddCasting(60 - topAngle5, sideAngle5, 2, castings);
                            


                            CheckCastings(tas[3], new List<double> { sas[6] }, castings, moldInd, castingIndex, new List<int> { 6 }, machinedCastings, castings[castingIndex[5]].TopAngle2, 2);

                        }



                        // angle 10
                        if (nextCasting_0 == null)
                            if (!CheckCastings(tas[0], new List<double> { sas[0], sas[1] }, castings, moldInd, castingIndex, new List<int> { 0, 1 }, machinedCastings))
                            {
                                castingIndex[0] = castings.Count;
                                moldInd.Add(castings.Count);
                                AddCasting(60 - topAngle0, sas[0], 2, castings);
                                

                                CheckCastings(tas[0], new List<double> { sas[1] }, castings, moldInd, castingIndex, new List<int> { 1 }, machinedCastings, castings[castings.Count - 1].TopAngle2, 2);

                            }
                            else { }

                        else
                            CheckCastings(tas[0], new List<double> { sas[1] }, castings, moldInd, castingIndex, new List<int> { 1 }, machinedCastings, nextCasting_0.TopAngle2, 2);


                        // angle 98, angel 87

                        if (nextCasting_9 == null)
                            if (!CheckCastings(tas[5], new List<double> { sas[9], sas[8] }, castings, moldInd, castingIndex, new List<int> { 9, 8 }, machinedCastings))
                            {
                                castingIndex[9] = castings.Count;
                                moldInd.Add(castings.Count);
                                AddCasting(topAngle9, sas[9], 1, castings);
                                

                                CheckCastings(tas[5], new List<double> { sas[8] }, castings, moldInd, castingIndex, new List<int> { 8 }, machinedCastings, castings[castings.Count - 1].TopAngle1, 1);

                            }
                            else { }

                        else
                            CheckCastings(tas[5], new List<double> { sas[8] }, castings, moldInd, castingIndex, new List<int> { 8 }, machinedCastings, nextCasting_9.TopAngle1, 1);
                        CheckCastings(tas[4], new List<double> { sas[7] }, castings, moldInd, castingIndex, new List<int> { 7 }, machinedCastings, castings[castingIndex[8]].TopAngle1, 1);


                        // angle 43, angel 32
                        if (!CheckCastings(tas[2], new List<double> { sas[4], sas[3] }, castings, moldInd, castingIndex, new List<int> { 4, 3 }, machinedCastings))
                        {
                            castingIndex[4] = castings.Count;
                            moldInd.Add(castings.Count);
                            AddCasting(topAngle4, sas[4], 1, castings);
                            

                            CheckCastings(tas[2], new List<double> { sas[3] }, castings, moldInd, castingIndex,  new List<int> { 3 }, machinedCastings, castings[castingIndex[4]].TopAngle1, 1);
                        }

                        CheckCastings(tas[1], new List<double> { sas[2] }, castings, moldInd, castingIndex, new List<int> { 2 }, machinedCastings, castings[castingIndex[3]].TopAngle1, 1);
                    }
                    else
                    {


                        castingIndex[6] = castingIndex_6[branch_current];
                        castingIndex[7] = castingIndex_7[branch_current];
                        //56 because 6 is found from below
                        CheckCastings(tas[3], new List<double> { sas[5] }, castings, moldInd, castingIndex, new List<int> { 5 }, machinedCastings, castings[castingIndex_6[branch_current]].TopAngle1, 1);





                        // angle 10
                        if (nextCasting_0 == null)
                            if (!CheckCastings(tas[0], new List<double> { sas[0], sas[1] }, castings, moldInd, castingIndex, new List<int> { 0, 1 }, machinedCastings))
                            {
                                castingIndex[0] = castings.Count;
                                moldInd.Add(castings.Count);
                                AddCasting(60 - topAngle0, sas[0], 2, castings);
                                
                                CheckCastings(tas[0], new List<double> { sas[1] }, castings, moldInd, castingIndex, new List<int> { 1 }, machinedCastings, castings[castings.Count - 1].TopAngle2, 2);

                            }
                            else { }

                        else
                            CheckCastings(tas[0], new List<double> { sas[1] }, castings, moldInd, castingIndex, new List<int> { 1 }, machinedCastings, nextCasting_0.TopAngle2, 2);


                        // angle 98, angel 87

                        if (nextCasting_9 == null)
                            if (!CheckCastings(tas[5], new List<double> { sas[9], sas[8] }, castings, moldInd, castingIndex, new List<int> { 9, 8 }, machinedCastings))
                            {



                                CheckCastings(tas[4], new List<double> { sas[8] }, castings, moldInd, castingIndex, new List<int> { 8 }, machinedCastings, castings[castingIndex_7[branch_current]].TopAngle2, 2);

                            }
                            else { }

                        else
                            CheckCastings(tas[5], new List<double> { sas[8] }, castings, moldInd, castingIndex, new List<int> { 8 }, machinedCastings, nextCasting_9.TopAngle1, 1);


                        // angle 43, angel 32
                        if (!CheckCastings(tas[2], new List<double> { sas[4], sas[3] }, castings, moldInd, castingIndex, new List<int> { 4, 3 }, machinedCastings))
                        {
                            castingIndex[4] = castings.Count;
                            moldInd.Add(castings.Count);
                            AddCasting(topAngle4, sas[4], 1, castings);
                            
                            CheckCastings(tas[2], new List<double> { sas[3] }, castings, moldInd, castingIndex, new List<int> { 3 }, machinedCastings, castings[castingIndex[4]].TopAngle1, 1);
                        }

                        CheckCastings(tas[1], new List<double> { sas[2] }, castings, moldInd, castingIndex, new List<int> { 2 }, machinedCastings, castings[castingIndex[3]].TopAngle1, 1);
                    } 
                    // after the angles
                    finalMesh = FindNewMesh(p1New_current, p, domeMesh);

                    nextCastingIndex_0 = castingIndex[4];
                    nextCastingIndex_9 = castingIndex[5];

                    nextCasting_0 = castings[castingIndex[4]];
                    nextCasting_9 = castings[castingIndex[5]];

                    List<Line> holderOriginalLines = new List<Line>();
                    holderOriginalLines.Add(new Line(p[2], p[3]));
                    holderOriginalLines.Add(new Line(p[1], p[3]));
                    holderOriginalLines.Add(new Line(p[2], p[0]));
                    holderOriginalLines.Add(new Line(p[0], p[3]));
                    holderOriginalLines.Add(new Line(p[0], p[1]));
                    originalLines.AddRange(holderOriginalLines, pth);


                    List<Line> holderChangeLines = new List<Line>();
                    holderChangeLines.Add(new Line(p[2], p[3]));
                    holderChangeLines.Add(new Line(p[1], p[3]));
                    holderChangeLines.Add(new Line(p[2], p1New_current));
                    holderChangeLines.Add(new Line(p1New_current, p[3]));
                    holderChangeLines.Add(new Line(p1New_current, p[1]));
                    newLines.AddRange(holderChangeLines, pth);

                    List<Point3d> newP_current = new List<Point3d>();
                    newP_current.Add(p1New_current);
                    newP_current.Add(p[1]);
                    newP_current.Add(p[2]);
                    newP_current.Add(p[3]);
                    newP.AddRange(newP_current, pth);




                    topAnglesFinal.AddRange(FindAnglesTop(newP_current, finalPlanes), pth);
                    sideAnglesFinal.AddRange(FindAnglesSide(newP_current, finalPlanes), pth);
                    sidePlanes.AddRange(FindPlanesSide(newP_current, finalPlanes), pth);
                    costs.Add(cost, pth);
                    angleCosts.Add(angleCost, pth);
                    newP1.Add(p1New_current, pth);
                    topPlanes.AddRange(finalPlanes, pth);
                    meshes.Add(finalMesh, pth);
                    castingPositions.AddRange(castingIndex, pth);
                }




                else
                {
                    failed = true;
                    topAnglesFinal.AddRange(new List<double>(), pth);
                    sideAnglesFinal.AddRange(new List<double>(), pth);
                    sidePlanes.AddRange(new List<Plane>(), pth);
                    costs.Add(dist * 300, pth);
                    angleCosts.Add(double.PositiveInfinity, pth);
                    newP1.Add(new Point3d(), pth);
                    topPlanes.AddRange(new List<Plane>(), pth);
                    meshes.Add(new Mesh(), pth);

                }
                pLast = p1New_current;
                tl_current += 1;
                dl_current += 2;
                branch_current += 1;


            }

            DA.SetDataTree(0, originalLines);
            DA.SetDataTree(1, newP);
            DA.SetDataTree(2, newLines);
            DA.SetDataTree(3, topAnglesFinal);
            DA.SetDataTree(4, sideAnglesFinal);
            DA.SetDataTree(5, newP1);
            DA.SetDataTree(6, costs);
            DA.SetDataTree(7, angleCosts);
            DA.SetDataTree(8, sidePlanes);
            DA.SetDataTree(9, topPlanes);
            DA.SetDataTree(10, meshes);
            if (!failed)
                DA.SetData(11, dist);
            else
                DA.SetData(11, dist + 1 * 300);
            DA.SetDataTree(12, CastingsToTree(castings));
            DA.SetDataTree(13, CastingsToTree(machinedCastings));
            DA.SetDataTree(14, castingPositions);
            DA.SetDataList(15, moldInd);
        }

        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// You can add image files to your project resources and access them like this:
        /// return Resources.IconForThisComponent;
        /// </summary>
        protected override System.Drawing.Bitmap Icon => null;

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid => new Guid("28277592-c9f7-4799-95ff-c5881d439827");
    }
}