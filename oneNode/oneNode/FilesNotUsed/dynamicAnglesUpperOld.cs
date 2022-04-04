using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using Rhino;
using Rhino.Geometry;

namespace oneNode
{
    public class dynamicAnglesUpperComponentOld : GH_Component
    {
        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public dynamicAnglesUpperComponentOld()
          : base("dynamic angles upper OLD", "DA upper OLD",
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
            pManager.AddNumberParameter("Angle1", "A1", "Angle 1 (known as theta2)", GH_ParamAccess.item);
            pManager.AddNumberParameter("Angle2", "A2", "Angle 2 (known as theta42)", GH_ParamAccess.item);
            pManager.AddNumberParameter("SigmaAngle2", "SA2", "Sigma Angle 2", GH_ParamAccess.item);
            pManager.AddPlaneParameter("Plane2", "PL2", "Plane for the 2 node", GH_ParamAccess.item);
            pManager.AddPlaneParameter("Plane4", "PL4", "Plane for the 4 node", GH_ParamAccess.item);

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager) {
            pManager.Register_LineParam("Original Lines", "oLines", "The lines coming into the component");
            pManager.Register_PointParam("Changed Points", "pOut", "The new points changed based on the angles");
            pManager.Register_LineParam("Changed Lines", "newLines", "The lines after angle has been changed");
            pManager.Register_DoubleParam("Angles", "aOut", "The angles of the node going from top corner counterclockwise");
           pManager.Register_PointParam("Testing point", "pOutTest", "The new points changed based on the angles");
            pManager.Register_PlaneParam("Testing plane parm", "testPL", "The new points changed based on the angles");
            pManager.Register_LineParam("Testing plane parm", "testLN", "The new points changed based on the angles");


        }

        // Find angles based on 4 joint
        public Point3d FindXandYandZ(double theta2, double thetaS2, double theta42, Plane plane4, Plane plane2, List<Point3d> p, ref Plane test, ref Line testln) {
            
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
            testln = holderLn;
            test = planeInt;



            return holderLn.PointAt(t);
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

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            double theta2 = 0.0;
            double theta42 = 0.0;
            double thetaS2 = 0.0;
            Plane plane2 = new Plane();
            Plane plane4 = new Plane();

            List<Point3d> p = new List<Point3d> ();


            if (!DA.GetDataList(0, p)) { return; }
            if (!DA.GetData(1, ref theta2)) { return; }
            if (!DA.GetData(2, ref theta42)) { return; }
            if (!DA.GetData(3, ref thetaS2)) { return; }
            if (!DA.GetData(4, ref plane2)) { return; }
            if (!DA.GetData(5, ref plane4)) { return; }

            theta42 = RhinoMath.ToRadians(theta42);
            theta2 = RhinoMath.ToRadians(theta2);
            thetaS2 = RhinoMath.ToRadians(thetaS2);

            List<Line> holderOriginalLines = new List<Line>();
            holderOriginalLines.Add(new Line(p[2], p[3]));
            holderOriginalLines.Add(new Line(p[1], p[3]));
            holderOriginalLines.Add(new Line(p[2], p[0]));
            holderOriginalLines.Add(new Line(p[0], p[3]));
            holderOriginalLines.Add(new Line(p[0], p[1]));

            double y = 0.0;
            double theta43 = 0.0;
            double theta12 = 0.0;
            
            double x = FindXandY(theta2, theta42, p, ref y, ref theta43, ref theta12);

            double x3 = (p[2] - p[3])[0];
            double y3 = (p[2] - p[3])[1];
            double l41 = Math.Sqrt(Math.Pow(x, 2) + Math.Pow(y, 2));
            double l31 = Math.Sqrt(Math.Pow((x - x3), 2) + Math.Pow((y - y3), 2));
            double theta3 = Math.Asin(l41 * Math.Sin(theta43) / l31);
            double theta13 = Math.PI - theta43 - theta3;


            

            Point3d p1New = new Point3d(p[3].X + x, p[3].Y + y, p[3].Z);

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
           


            List<double> angles = new List<double>();
            angles.Add(RhinoMath.ToDegrees(theta3));
            angles.Add(RhinoMath.ToDegrees(theta13));
            angles.Add(RhinoMath.ToDegrees(theta12));
            angles.Add(RhinoMath.ToDegrees(theta2));
            angles.Add(RhinoMath.ToDegrees(theta42));
            angles.Add(RhinoMath.ToDegrees(theta43));
            Plane test = new Plane();
            Line testln = new Line();   

            DA.SetDataList(0, holderOriginalLines);
            DA.SetDataList(1, newP);
            DA.SetDataList(2, holderChangeLines);
            DA.SetDataList(3, angles);
            DA.SetData(4, FindXandYandZ(theta2, thetaS2, theta42, plane4, plane2, p, ref test, ref testln));
            DA.SetData(5, test);
            DA.SetData(6, testln);



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
        public override Guid ComponentGuid => new Guid("9436da66-23de-476d-a475-6cac2aa44a96");
    }
}