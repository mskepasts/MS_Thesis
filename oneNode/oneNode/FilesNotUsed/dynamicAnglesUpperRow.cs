using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Rhino;
using Rhino.Geometry;

namespace oneNode
{
    public class dynamicAnglesUpperRowComponent : GH_Component
    {
        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public dynamicAnglesUpperRowComponent()
          : base("dynamic angles upper Row", "DA upper Row",
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
            pManager.AddNumberParameter("Top Angles", "TAs", "All the top angles used to test", GH_ParamAccess.list);
            pManager.AddNumberParameter("Side Angles", "SAs", "All the side angles used to test", GH_ParamAccess.list);
            pManager.AddMeshParameter("Dome mesh", "DM", "Dome mesh", GH_ParamAccess.item);
            

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager) {
            pManager.Register_LineParam("Original Lines", "oLines", "The lines coming into the component");
            pManager.Register_PointParam("Changed Points", "pOut", "The new points changed based on the angles");
            pManager.Register_LineParam("Changed Lines", "newLines", "The lines after angle has been changed");
            pManager.Register_DoubleParam("AnglesTop", "aOutT", "The angles of the node going from top corner counterclockwise");
            pManager.Register_DoubleParam("AnglesSide", "aOutS", "The angles of the node going from top corner counterclockwise");
            pManager.Register_PointParam("Testing point", "pOutTest", "The new points changed based on the angles");
            pManager.Register_DoubleParam("Cost", "Cost", "Cost in order to change the node to fit the required angles");
            pManager.Register_DoubleParam("AngleCost", "ACost", "Cost in order to change the node to fit the required angles");
            pManager.Register_PlaneParam("side planes", "testSP", "sideplane test");
            pManager.Register_PlaneParam("top planes", "testTP", "topplane test");
            pManager.Register_MeshParam("New meshes", "testMesh", "mesh test");
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

            if(points.Contains(diagonalOne.From))
                points.Add(diagonalOne.To);
            else
                points.Add(diagonalOne.From);

            return points;
        }

        // Change out point from last angle being found
        public List<Point3d> GetAlteredPoints(Point3d pNew, List<Point3d> p) {
            p[2] = pNew;
            return p;
        }

        // Find side plane angle
        public double FindingSideAngle(Point3d p, Point3d pOther, Plane basePlane) {
            Point3d holderPt = p;
            var xformPt = Rhino.Geometry.Transform.Translation(basePlane.ZAxis);
            holderPt.Transform(xformPt);
            Plane planeSide = new Plane();
            Plane.FitPlaneToPoints(new List<Point3d> { p, basePlane.ClosestPoint(pOther), holderPt }, out planeSide);
            planeSide.Origin = p;
            double angle = Vector3d.VectorAngle(pOther - p, basePlane.ClosestPoint(pOther) - p, planeSide);
            if (angle > Math.PI) {angle = 2*Math.PI - angle;}
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
        public Point3d FindXandYandZ(double theta2, double thetaS2, double theta42, Plane plane2, Plane plane4, List<Point3d> p) {
            theta42 = RhinoMath.ToRadians(theta42);
            theta2 = RhinoMath.ToRadians(theta2);
            thetaS2 = RhinoMath.ToRadians(thetaS2);


            Line holderLn = new Line(p[1], plane2.ClosestPoint(p[3]));
            Plane planeS = new Plane();
            Plane.FitPlaneToPoints(new List<Point3d>{p[1], p[3], plane2.ClosestPoint(p[3])}, out planeS);
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
            Plane.FitPlaneToPoints(new List<Point3d> { p[3], holderLn4.To, holderPt}, out planeInt);
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
            double theta43 = Vector3d.VectorAngle(pts[3] - pts[0], pts[3] - pts[2],  plns[3]);

            // Node three
            double theta3 = Vector3d.VectorAngle(pts[2] - pts[3], pts[2] - pts[0], plns[2]);
            
            // Node two
            double theta2 = Vector3d.VectorAngle(pts[1] - pts[0], pts[1] - pts[3],  plns[1]);

            // Node one
            double theta12 = Vector3d.VectorAngle(pts[0] - pts[3], pts[0] - pts[1],  plns[0]);
            double theta13 = Vector3d.VectorAngle(pts[0] - pts[2], pts[0] - pts[3],  plns[0]);

            List<double> angles = new List<double>();
            angles.Add(Math.Round(RhinoMath.ToDegrees(theta3),2));
            angles.Add(Math.Round(RhinoMath.ToDegrees(theta12),2));
            angles.Add(Math.Round(RhinoMath.ToDegrees(theta2),2));
            angles.Add(Math.Round(RhinoMath.ToDegrees(theta13),2));
            angles.Add(Math.Round(RhinoMath.ToDegrees(theta42),2));
            angles.Add(Math.Round(RhinoMath.ToDegrees(theta43),2));




            return angles;

        }

        // Find the side planes (debugging)
        public List<Plane> FindPlanesSide(List<Point3d> pts, List<Plane> plns)
        {

            List<Plane> planes = new List<Plane>();



            //planes.Add((FindingSidePlane(pts[2], pts[3], plns[2])));
            planes.Add((FindingSidePlane(pts[2], pts[0], plns[2])));
            planes.Add((FindingSidePlane(pts[0], pts[2], plns[0])));
            planes.Add((FindingSidePlane(pts[0], pts[3], plns[0])));
            planes.Add((FindingSidePlane(pts[0], pts[1], plns[0])));
            planes.Add((FindingSidePlane(pts[1], pts[0], plns[1])));
            //planes.Add((FindingSidePlane(pts[1], pts[3], plns[1])));
            //planes.Add((FindingSidePlane(pts[3], pts[1], plns[3])));
            planes.Add((FindingSidePlane(pts[3], pts[0], plns[3])));
            //planes.Add((FindingSidePlane(pts[3], pts[2], plns[3])));

            return planes;

        }

        // Find angles side (sigma)
        public List<double> FindAnglesSide(List<Point3d> pts, List<Plane> plns)
        {
         
            List<double> angles = new List<double>();



            angles.Add(Math.Round(RhinoMath.ToDegrees(FindingSideAngle(pts[2], pts[3], plns[2])),2));
            angles.Add(Math.Round(RhinoMath.ToDegrees(FindingSideAngle(pts[2], pts[0], plns[2])),2));
            angles.Add(Math.Round(RhinoMath.ToDegrees(FindingSideAngle(pts[0], pts[2], plns[0])),2));
            angles.Add(Math.Round(RhinoMath.ToDegrees(FindingSideAngle(pts[0], pts[3], plns[0])),2));
            angles.Add(Math.Round(RhinoMath.ToDegrees(FindingSideAngle(pts[0], pts[1], plns[0])),2));
            angles.Add(Math.Round(RhinoMath.ToDegrees(FindingSideAngle(pts[1], pts[0], plns[1])),2));
            angles.Add(Math.Round(RhinoMath.ToDegrees(FindingSideAngle(pts[1], pts[3], plns[1])),2));
            angles.Add(Math.Round(RhinoMath.ToDegrees(FindingSideAngle(pts[3], pts[1], plns[3])),2));
            angles.Add(Math.Round(RhinoMath.ToDegrees(FindingSideAngle(pts[3], pts[0], plns[3])),2));
            angles.Add(Math.Round(RhinoMath.ToDegrees(FindingSideAngle(pts[3], pts[2], plns[3])),2));

            return angles;

        }

        // Find x and y based on the angles given
        public double FindXandY(double theta2_4, double theta42, List<Point3d> p,  ref double y, ref double theta43, ref double theta12)
        {
        double x2 = (p[1] - p[3])[0];
        double y2 = (p[1] - p[3])[1];

        double thetaRest = 360 - RhinoMath.ToDegrees(Vector3d.VectorAngle((p[2] - p[3]), (p[1] - p[3])));
        theta12 = 180 - RhinoMath.ToDegrees(theta42) - RhinoMath.ToDegrees(theta2_4);
        theta43 = 360 - RhinoMath.ToDegrees(theta42) - thetaRest;
        theta43 = RhinoMath.ToRadians(theta43);
        theta12 = RhinoMath.ToRadians(theta12);

        double l42 = Math.Sqrt(Math.Pow(x2, 2) + Math.Pow(y2, 2));



        double c1 = Math.Sin(theta2_4) * l42 / Math.Sin(theta12);
        double c2 = Math.Sin(theta42) * l42 / Math.Sin(theta12);
        double c3 = (Math.Pow(c1, 2) - Math.Pow(c2, 2) + Math.Pow(x2, 2) + Math.Pow(y2, 2)) / (2 * x2);
        double c4 = -y2 / x2;
        double c5 = Math.Pow(c4, 2) + 1;
        double c6 = Math.Pow(c3, 2) - Math.Pow(c1, 2);
        double c7 = c3 * c4 * 2;
        
        y = (-c7 + Math.Sqrt(Math.Pow(c7, 2) - (4 * c6 * c5))) / (2 * c5);
        double y_ = (-c7 - Math.Sqrt(Math.Pow(c7, 2) - (4 * c6 * c5))) / (2 * c5);

        double x = c3 - ((y * y2) / (x2));
        double x_ = c3 - ((y_ * y2) / (x2));

            if (x< 0)
            {
                x = x_;
                y = y_;
            }
            return x;

        }

        // Finding the max of the min
        public double FindMinMax(List<double> anglesNew, List<double> angles) {
            double minCost = 0.0;
            double holder = 0.0;
            int factor = 10;
            List<double> anglesT = new List<double>();
            int length = angles.Count;
            for (int j = 0; j < length; j++)
            {
                double angle = angles[j];
                anglesT.Add(angle);
                for (int i = 0; i < 5; i++)
                {
                    anglesT.Add(angle + ((i*1.0) / factor));
                    anglesT.Add(angle - ((i * 1.0) / factor));
                }
            }

            foreach (double angleNew in anglesNew)
            {
               
                holder = double.PositiveInfinity;
                foreach (double angle in anglesT)
                    if (Math.Abs(holder) > Math.Abs(angle - angleNew))
                        holder = (angle - angleNew);
                if (Math.Abs(holder) > Math.Abs(minCost))
                    minCost = holder;
            }
            
            return Math.Round(minCost, 2);
        }

        public double FindMinMaxBoth(List<double> anglesNewTop, List<double> anglesNewSide, List<double> anglesTop, List<double> anglesSide)
        {
            double minCost = 0.0;
            double holder = 0.0;
            int factor = 10;
            List<double> anglesT = new List<double>();
            List<double> anglesS = new List<double>();
            int length = anglesSide.Count;
            for (int j = 0; j < length; j++)
            {
                double angleTop = anglesTop[j];
                double angleSide = anglesSide[j];
                anglesT.Add(angleTop);
                anglesS.Add(angleSide);
                for (int i = 0; i < 10; i++)
                {
                    anglesT.Add(angleTop + ((i * 1.0) / factor));
                    anglesT.Add(angleTop - ((i * 1.0) / factor));
                    anglesS.Add(angleSide + ((i * 1.0) / factor));
                    anglesS.Add(angleSide - ((i * 1.0) / factor));
                }
            }

            for (int i = 0; i < anglesNewSide.Count; i++)
            {
                

                holder = double.PositiveInfinity;
                for (int j = 0; j < anglesT.Count; j++)
                {
                    double holderT = Math.Abs(anglesT[j] - anglesNewTop[i]);
                    double holderS = Math.Abs(anglesS[j] - anglesNewSide[i]);

                    if (holderT < 0.02)
                        holderT = 0;
                    if (holderS < 0.02)
                        holderS = 0;

                    if (Math.Abs(holder) > (holderT + holderS))
                        holder = holderT + holderS;
                }
            if (Math.Abs(holder) > Math.Abs(minCost))
                    minCost = holder;
            }

            return Math.Round(minCost, 2);
        }


        // Finding the cost of the function
        public double Cost(List<double> anglesTop, List<double> anglesSide, List<Plane> planes, List<Point3d> p, Point3d pointFound, ref double angleCost) {
            List<Point3d> newP = new List<Point3d>();
            newP.Add(pointFound);
            newP.Add(p[1]);
            newP.Add(p[2]);
            newP.Add(p[3]);
            double minCost = FindMinMax(FindAnglesTop(newP, planes), anglesTop) + FindMinMax(FindAnglesSide(newP, planes), anglesSide);
            

            double minCostDist = pointFound.DistanceTo(p[0]);
            double scaleAngles = 20.0;
            double scaleDist = 1.0;

            angleCost = minCost;
            
            double cost = scaleAngles * Math.Abs(angleCost) + scaleDist * minCostDist;
            return cost; 
        }

        // Finding the top planes
        public List<Plane> GetTopPlanes(List<Point3d> pts, Mesh domeMesh) {
            List<Plane> planes = new List<Plane>();
            foreach(Point3d pt in pts)
            {
                Vector3d normal = domeMesh.NormalAt(domeMesh.ClosestMeshPoint(pt,1));
                planes.Add(new Plane(pt, normal));
            }

            return planes;
        }

        public Mesh FindNewMesh(Point3d pNew, List<Point3d> p, Mesh grid)
        {
            Mesh newGrid = new Mesh();
            var nodes = new Grasshopper.Kernel.Geometry.Node2List();

            List<Point3d> pts = new List<Point3d>(); 
            foreach(Point3d point in grid.Vertices.ToPoint3dArray())
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
            

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<double> topAngles = new List<double>();
            List<double> sideAngles = new List<double>();

            
            
            List<Line> topLines = new List<Line> ();
            List<Line> diagonalLines = new List<Line>();

            int topLineIndexStart = 0;
            int diagonalLineIndexStart = 0;            
            
            Mesh domeMesh = new Mesh();

            if (!DA.GetDataList(0, topLines)) { return; }
            if (!DA.GetDataList(1, diagonalLines)) { return; }
            if (!DA.GetData(2, ref topLineIndexStart)) { return; }
            if (!DA.GetData(3, ref diagonalLineIndexStart)) { return; }
            if (!DA.GetDataList(4, topAngles)) { return; }
            if (!DA.GetDataList(5, sideAngles)) { return; }
            if (!DA.GetData(6, ref domeMesh)) { return; }


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
            DataTree<Mesh> meshes = new DataTree<Mesh>();
            Mesh finalMesh = domeMesh;
            DataTree<int> numOfAngles = new DataTree<int>();
            for (int i = topLineIndexStart; i < topLines.Count -1; i++) 
            {
                
                GH_Path pth = new GH_Path(branch_current);
                domeMesh = finalMesh;
                List<Point3d> p = GetEvalPoints(topLines[tl_current], diagonalLines[dl_current], diagonalLines[dl_current + 1]);
                List<Plane> planes = GetTopPlanes(p, domeMesh);
                List<Plane> finalPlanes = new List<Plane>();
                List<Plane> holderPlanes = new List<Plane>();
                double holderCost = 0.0;
                double angleCost = double.PositiveInfinity;
                double cost = double.PositiveInfinity;
                double otherCost = 0.0;
                double dist = double.PositiveInfinity;

                if (tl_current != topLineIndexStart)
                    p = GetAlteredPoints(pLast, p);
                
                Point3d holderPoint = new Point3d(); 
                Point3d p1New_current = p[0];
                Mesh holderMesh = new Mesh();
                
                
                
                foreach (double sideAngle in sideAngles)
                {
                    foreach (double topAngle1 in topAngles)
                    {
                        foreach (double topAngle2 in topAngles)
                        {

                            holderPoint = FindXandYandZ(topAngle1, sideAngle, topAngle2, planes[1], planes[3], p);
                            if(holderPoint.DistanceTo(p[0]) < dist)
                                dist = holderPoint.DistanceTo(p[0]);
                            if (holderPoint.DistanceTo(p[0]) < 0.3)
                            {
                                if (Cost(topAngles, sideAngles, planes, p, holderPoint, ref otherCost) < 4.0)
                                {

                                    holderMesh = FindNewMesh(holderPoint, p, domeMesh);
                                    //holderMesh = domeMesh;

                                    List<Point3d> holderP_current = new List<Point3d>();
                                    holderP_current.Add(holderPoint);
                                    holderP_current.Add(p[1]);
                                    holderP_current.Add(p[2]);
                                    holderP_current.Add(p[3]);

                                    holderPlanes = GetTopPlanes(holderP_current, holderMesh);
                                    holderCost = Cost(topAngles, sideAngles, planes, p, holderPoint, ref otherCost);
                                    if (holderCost < cost)
                                    {
                                        finalPlanes = holderPlanes;
                                        p1New_current = holderPoint;
                                        cost = holderCost;
                                        finalMesh = holderMesh;
                                    }
                                    if (otherCost < angleCost)
                                    {
                                        angleCost = otherCost;
                                    }
                                }
                                else
                                {
                                    holderMesh = domeMesh;
                                    holderPlanes = GetTopPlanes(p, holderMesh);
                                    holderCost = Cost(topAngles, sideAngles, planes, p, holderPoint, ref otherCost);
                                    if (holderCost < cost)
                                    {
                                        finalPlanes = holderPlanes;
                                        p1New_current = holderPoint;
                                        cost = holderCost;
                                        finalMesh = holderMesh;
                                    }
                                    if (otherCost < angleCost)
                                    {
                                        angleCost = otherCost;
                                    }


                                }
                            }
                        }


                    }
                }

                if (cost == double.PositiveInfinity)
                {
                    topAnglesFinal.AddRange(new List<double>(), pth);
                    sideAnglesFinal.AddRange(new List<double>(), pth);
                    sidePlanes.AddRange(new List<Plane>(), pth);
                    costs.Add(dist*30, pth);
                    angleCosts.Add(double.PositiveInfinity, pth);
                    newP1.Add(new Point3d(), pth);
                    topPlanes.AddRange(new List<Plane>(), pth) ;
                    meshes.Add(new Mesh(), pth);

                }
                else
                {
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
        public override Guid ComponentGuid => new Guid("1e371173-b523-4dd0-89b8-4e76e3107aaf");
    }
}