using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using Rhino;
using Rhino.Geometry;

namespace oneNode
{
    public class dynamicAnglesUpperSingleComponent : GH_Component
    {
        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public dynamicAnglesUpperSingleComponent()
          : base("dynamic angles upper Single", "DA upper Single",
            "Description",
            "Category", "Subcategory")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Points", "p", "Points to define geometry", GH_ParamAccess.list);
            pManager.AddNumberParameter("Top Angles", "TAs", "All the top angles used to test", GH_ParamAccess.list);
            pManager.AddNumberParameter("Side Angles", "SAs", "All the side angles used to test", GH_ParamAccess.list);
            pManager.AddPlaneParameter("Planes", "PLs", "Planes nodes from 1-4", GH_ParamAccess.list);
            

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
            angles.Add(RhinoMath.ToDegrees(theta3));
            angles.Add(RhinoMath.ToDegrees(theta13));
            angles.Add(RhinoMath.ToDegrees(theta12));
            angles.Add(RhinoMath.ToDegrees(theta2));
            angles.Add(RhinoMath.ToDegrees(theta42));
            angles.Add(RhinoMath.ToDegrees(theta43));




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



            angles.Add(RhinoMath.ToDegrees(FindingSideAngle(pts[2], pts[3], plns[2])));
            angles.Add(RhinoMath.ToDegrees(FindingSideAngle(pts[2], pts[0], plns[2])));
            angles.Add(RhinoMath.ToDegrees(FindingSideAngle(pts[0], pts[2], plns[0])));
            angles.Add(RhinoMath.ToDegrees(FindingSideAngle(pts[0], pts[3], plns[0])));
            angles.Add(RhinoMath.ToDegrees(FindingSideAngle(pts[0], pts[1], plns[0])));
            angles.Add(RhinoMath.ToDegrees(FindingSideAngle(pts[1], pts[0], plns[1])));
            angles.Add(RhinoMath.ToDegrees(FindingSideAngle(pts[1], pts[3], plns[1])));
            angles.Add(RhinoMath.ToDegrees(FindingSideAngle(pts[3], pts[1], plns[3])));
            angles.Add(RhinoMath.ToDegrees(FindingSideAngle(pts[3], pts[0], plns[3])));
            angles.Add(RhinoMath.ToDegrees(FindingSideAngle(pts[3], pts[2], plns[3])));

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
            foreach (double angleNew in anglesNew)
            {
               
                holder = double.PositiveInfinity;
                foreach (double angle in angles)
                    if (holder > Math.Abs(angle - angleNew))
                        holder = Math.Abs(angle - angleNew);
                if (holder > minCost)
                    minCost = holder;
            }
            return minCost;
        }

        // Finding the cost of the function
        public double Cost(List<double> anglesTop, List<double> anglesSide, List<Plane> planes, List<Point3d> p, Point3d pointFound, ref double angleCost) {
            List<Point3d> newP = new List<Point3d>();
            newP.Add(pointFound);
            newP.Add(p[1]);
            newP.Add(p[2]);
            newP.Add(p[3]);
            double minCostTop = FindMinMax(FindAnglesTop(newP, planes), anglesTop);
            double minCostSide = FindMinMax(FindAnglesSide(newP, planes), anglesSide);
            double minCostDist = pointFound.DistanceTo(p[0]);
            double scaleAngles = 4.0;
            double scaleDist = 5.0;
            angleCost = (minCostTop + minCostSide);
            
            double cost = scaleAngles * angleCost + scaleDist * minCostDist;
            return cost; 
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

            double cost = double.PositiveInfinity;
            List<Plane> planes = new List<Plane>();
            
            List<Point3d> p = new List<Point3d> ();
            double angleCost = 0.0;

            if (!DA.GetDataList(0, p)) { return; }
            if (!DA.GetDataList(1, topAngles)) { return; }
            if (!DA.GetDataList(2, sideAngles)) { return; }
            if (!DA.GetDataList(3, planes)) { return; }
            Point3d p1New = p[0];
            double holderCost = 0.0;
            double otherCost = 0.0;
            foreach (double sideAngle in sideAngles)
            {
                foreach (double topAngle1 in topAngles)
                {
                    foreach (double topAngle2 in topAngles)
                    {
                        holderCost = Cost(topAngles, sideAngles, planes, p, FindXandYandZ(topAngle1, sideAngle, topAngle2, planes[1], planes[3], p), ref otherCost);
                        if (holderCost < cost)
                        {
                            Cost(topAngles, sideAngles, planes, p, FindXandYandZ(topAngle1, sideAngle, topAngle2, planes[1], planes[3], p), ref angleCost);
                            p1New = FindXandYandZ(topAngle1, sideAngle, topAngle2, planes[1], planes[3], p);
                            cost = holderCost;
                        }
                    }


                }
            }



            List<Line> holderOriginalLines = new List<Line>();
            holderOriginalLines.Add(new Line(p[2], p[3]));
            holderOriginalLines.Add(new Line(p[1], p[3]));
            holderOriginalLines.Add(new Line(p[2], p[0]));
            holderOriginalLines.Add(new Line(p[0], p[3]));
            holderOriginalLines.Add(new Line(p[0], p[1]));


            List<Line> holderChangeLines = new List<Line>();
            holderChangeLines.Add(new Line(p[2], p[3]));
            holderChangeLines.Add(new Line(p[1], p[3]));
            holderChangeLines.Add(new Line(p[2], p1New));
            holderChangeLines.Add(new Line(p1New, p[3]));
            holderChangeLines.Add(new Line(p1New, p[1]));

            List<Point3d> newP = new List<Point3d>();
            newP.Add(p1New);
            newP.Add(p[1]);
            newP.Add(p[2]);
            newP.Add(p[3]);
           
           
            DA.SetDataList(0, holderOriginalLines);
            DA.SetDataList(1, newP);
            DA.SetDataList(2, holderChangeLines);
            DA.SetDataList(3, FindAnglesTop(newP, planes));
            DA.SetDataList(4, FindAnglesSide(newP, planes));
            DA.SetData(5, p1New);
            DA.SetData(6, cost);
            DA.SetData(7, angleCost);
            DA.SetDataList(8, FindPlanesSide(newP, planes));



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
        public override Guid ComponentGuid => new Guid("cb570e5c-c2f4-4685-8dfb-563b8fe75f96");
    }
}