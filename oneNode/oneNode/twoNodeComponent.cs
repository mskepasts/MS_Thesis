using System;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Rhino;
using Rhino.Geometry;

namespace oneNode
{
    public class twoNodeComponent : GH_Component
    {
        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public twoNodeComponent()
          : base("twoNode", "2n",
            "Description",
            "Category", "Subcategory")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("angle10","a10", "angle 10", GH_ParamAccess.list);
            pManager.AddNumberParameter("angle11", "a11", "angle 10", GH_ParamAccess.list);
            pManager.AddNumberParameter("angle20", "a20", "angle 10", GH_ParamAccess.list);
            pManager.AddNumberParameter("angle21", "a21", "angle 10", GH_ParamAccess.list);
            pManager.AddLineParameter("lines at node", "Ln", "lines at the node", GH_ParamAccess.list);
            pManager.AddPlaneParameter("topPlane","tp","top plane", GH_ParamAccess.item);
            pManager.AddPlaneParameter("sidePlanes", "sp", "side plane", GH_ParamAccess.list);


        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.Register_LineParam("lines out", "ln", "the new fittd lines");
            pManager.Register_DoubleParam("cost", "cost", "the cost of the new fittd lines");
        }
        // <Custom additional code> 

        // returns the angle of the projected lines
        private List<double> angle1Projected(List<Line> lns, Plane pl)
        {
            List<double> angles = new List<double>();
            if (RhinoMath.ToDegrees(Vector3d.VectorAngle(lns[0].Direction, lns[1].Direction, pl)) > 180)
            {
                for (int i = 0; i < lns.Count; i++)
                    if (i < lns.Count - 1)
                        angles.Add(360.0 - RhinoMath.ToDegrees(Vector3d.VectorAngle(lns[i].Direction, lns[i + 1].Direction, pl)));
                    else
                        angles.Add(360.0 - RhinoMath.ToDegrees(Vector3d.VectorAngle(lns[i].Direction, lns[0].Direction, pl)));
                Globals.flipAngle = 1;
            }
            else
            {

                for (int i = 0; i < lns.Count; i++)
                    if (i < lns.Count - 1)
                        angles.Add(RhinoMath.ToDegrees(Vector3d.VectorAngle(lns[i].Direction, lns[i + 1].Direction, pl)));
                    else
                        angles.Add(RhinoMath.ToDegrees(Vector3d.VectorAngle(lns[i].Direction, lns[0].Direction, pl)));
                Globals.flipAngle = 0;
            }
            return angles;
        }

        // returns the starting angle and then the relative angles
        private List<double> angle2Projected(List<Line> lns, Plane pl, List<Plane> plns, out double angleStart)
        {
            List<double> angles = new List<double>();
            List<double> anglesOriginal = new List<double>();
            angleStart = 0.0;
            for (int i = 0; i < lns.Count; i++)
            {
                Line holderLn = new Line(lns[i].From, pl.ClosestPoint(lns[i].To));
                Plane holderPl = plns[i];
                anglesOriginal.Add(RhinoMath.ToDegrees(Vector3d.VectorAngle(lns[i].Direction, holderLn.Direction, holderPl)));
            }

            angleStart = anglesOriginal[0];
            for (int i = 0; i < lns.Count; i++)
            {
                var holder = i + 1;
                if (holder == lns.Count)
                    holder = 0;

                angles.Add(anglesOriginal[holder] - anglesOriginal[i]);

            }


            return angles;
        }

        // cost of first operation angle 1
        private double cost1(List<double> angles1, double angle1Single)
        {
            return Math.Abs(angles1[0] - angle1Single);
        }

        //cost of second operation angle 1
        private double cost2(List<double> angles1, double angle1Single, double angle1Rel)
        {
            return Math.Abs(angles1[1] - angle1Rel) + Math.Abs(angles1[2] - angle1Single);
        }

        // cost of third operation angle 1
        private double cost3(List<double> angles1, double angle1Single, double angle1Rel, double angle1RelStart)
        {
            return Math.Abs(angles1[3] - angle1Rel) + Math.Abs(angles1[4] - angle1Single) + Math.Abs(angles1[5] - angle1RelStart);
        }
        // cost of first operation angle 2
        private double cost1_2(List<double> angles2, double angle2Start_act, double angle2Single, double angle2Start)
        {
            return Math.Abs(angles2[0] - angle2Single) + Math.Abs(angle2Start_act - angle2Start);
        }

        //cost of second operation angle 2
        private double cost2_2(List<double> angles2, double angle2Single, double angle2Rel)
        {
            return Math.Abs(angles2[1] - angle2Rel) + Math.Abs(angles2[2] - angle2Single);
        }

        // cost of third operation angle 2
        private double cost3_2(List<double> angles2, double angle2Single, double angle2Rel, double angle2RelStart)
        {
            return Math.Abs(angles2[3] - angle2Rel) + Math.Abs(angles2[4] - angle2Single) + Math.Abs(angles2[5] - angle2RelStart);
        }

        //updates lines to current state
        private List<Line> changeLines(List<Line> lns, List<int> path, List<double> angles1Single, List<List<double>> angles1Rel, List<double> angles2Single, List<List<double>> angles2Rel, List<double> angles2Start, Plane pl, List<Plane> plns2)
        {
            List<Line> lns_Updated = new List<Line>();
            var xform = Rhino.Geometry.Transform.Rotation(RhinoMath.ToRadians(360 - angles1Single[path[0]]), pl.ZAxis, pl.Origin); ;
            var xform2 = Rhino.Geometry.Transform.Rotation(RhinoMath.ToRadians(angles2Start[path[0]]), plns2[0].ZAxis, plns2[0].Origin);
            Line holderLn = new Line(lns[0].From, pl.ClosestPoint(lns[0].To));
            holderLn.Transform(xform);
            lns_Updated.Add(holderLn);

            for (int i = 0; i < lns.Count - 1; i++)
            {

                xform = Rhino.Geometry.Transform.Rotation(RhinoMath.ToRadians(360 - angles1Single[path[i]]), pl.ZAxis, pl.Origin);
                if (i % 2 == 1)
                {
                    xform = Rhino.Geometry.Transform.Rotation(RhinoMath.ToRadians(360 - angles1Rel[path[i]][path[i + 1]]), pl.ZAxis, pl.Origin);
                }
                if (Globals.flipAngle == 0)
                {
                    xform = Rhino.Geometry.Transform.Rotation(RhinoMath.ToRadians(angles1Single[path[i]]), pl.ZAxis, pl.Origin);
                    if (i % 2 == 1)
                    {
                        xform = Rhino.Geometry.Transform.Rotation(RhinoMath.ToRadians(angles1Rel[path[i]][path[i + 1]]), pl.ZAxis, pl.Origin);
                    }
                }


                Line holder = new Line(lns_Updated[lns_Updated.Count - 1].From, lns_Updated[lns_Updated.Count - 1].To);
                holder.Transform(xform);
                if (i % 2 == 1)
                    xform2 = Rhino.Geometry.Transform.Rotation(RhinoMath.ToRadians(angles2Rel[path[i]][path[i + 1]]), plns2[i + 1].ZAxis, plns2[i + 1].Origin);
                else
                    xform2 = Rhino.Geometry.Transform.Rotation(RhinoMath.ToRadians(angles2Single[path[i]]), plns2[i + 1].ZAxis, plns2[i + 1].Origin);
                holder.Transform(xform2);


                lns_Updated.Add(holder);

            }
            return lns_Updated;
        }
        public static class Globals
        {
            public static Int32 flipAngle = 0;
        }


        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<Line> linesAtNodes = new List<Line>();
            List<double> angle10 = new List<double>();
            List<double> angle11 = new List<double>();
            List<double> angle20 = new List<double>();
            List<double> angle21 = new List<double>();
            Plane plane = new Plane();
            List<Plane> plane2 = new List<Plane>();

            DA.GetDataList(0, angle10);
            DA.GetDataList(1, angle11);
            DA.GetDataList(2, angle20);
            DA.GetDataList(3, angle21);
            DA.GetDataList(4, linesAtNodes);
            DA.GetData(5, ref plane);
            DA.GetDataList(6, plane2);
            int numberOfNodes = angle10.Count;
            // Set up the angle arrays
            double angle = 360 / 3;
            List<double> angles10_total = new List<double>();
            for (int i = 0; i < numberOfNodes; i++)
            {
                angles10_total.Add(angle10[i]);
            }
            for (int i = 0; i < numberOfNodes; i++)
            {
                angles10_total.Add(angle11[i]);
            }
            List<double> angles11_total = new List<double>();
            for (int i = 0; i < numberOfNodes; i++)
            {
                angles11_total.Add(angle11[i]);
            }
            for (int i = 0; i < numberOfNodes; i++)
            {
                angles11_total.Add(angle10[i]);
            }

            List<double> angles1Single = new List<double>();
            List<List<double>> angles1Rel = new List<List<double>>();
            for (int i = 0; i < numberOfNodes * 2; i++)
            {
                List<double> sublist = new List<double>();
                angles1Single.Add(angle - angles10_total[i] - angles11_total[i]);
                for (int j = 0; j < numberOfNodes * 2; j++)
                {
                    sublist.Add(angles11_total[i] + angles10_total[j]);
                }
                angles1Rel.Add(sublist);
            }
           


            //Setting up angle2
            numberOfNodes = angle20.Count;
            List<double> angles20_total = new List<double>();
            for (int i = 0; i < numberOfNodes; i++)
            {
                angles20_total.Add(angle20[i]);
            }
            for (int i = 0; i < numberOfNodes; i++)
            {
                angles20_total.Add(-angle21[i]);
            }
            List<double> angles21_total = new List<double>();
            for (int i = 0; i < numberOfNodes; i++)
            {
                angles21_total.Add(angle21[i]);
            }
            for (int i = 0; i < numberOfNodes; i++)
            {
                angles21_total.Add(-angle20[i]);
            }

            List<double> angles2Single = new List<double>();
            List<List<double>> angles2Rel = new List<List<double>>();
            for (int i = 0; i < numberOfNodes * 2; i++)
            {
                List<double> sublist = new List<double>();
                angles2Single.Add(angles21_total[i] - angles20_total[i]);
                for (int j = 0; j < numberOfNodes * 2; j++)
                {
                    sublist.Add(angles20_total[j] - angles21_total[i]);
                }
                angles2Rel.Add(sublist);
            }
            





            List<double> angles1 = angle1Projected(linesAtNodes, plane);
            double angle2Start;
            List<double> angles2 = angle2Projected(linesAtNodes, plane, plane2, out angle2Start);

            List<int> path = new List<int>();
            int count = angles1Single.Count;

            double optCost = Double.PositiveInfinity;
            for (int i = 0; i < count; i++)
            {
                double cost = 0;
                double cost_1 = cost1(angles1, angles1Single[i]) + cost1_2(angles2, angle2Start, angles2Single[i], angles20_total[i]);
                cost += cost_1;
                if (cost < optCost)
                {
                    for (int j = 0; j < count; j++)
                    {
                        cost = cost_1;
                        double cost_2 = cost2(angles1, angles1Single[j], angles1Rel[i][j]) + cost2_2(angles2, angles2Single[j], angles2Rel[i][j]);
                        cost += cost_2;
                        if (cost < optCost)
                        {
                            for (int q = 0; q < count; q++)
                            {
                                cost = cost_1 + cost_2;
                                double cost_3 = cost3(angles1, angles1Single[q], angles1Rel[j][q], angles1Rel[i][q]) + cost3_2(angles2, angles2Single[q], angles2Rel[j][q], angles2Rel[i][q]);
                                cost += cost_3;
                                if (cost < optCost)
                                {
                                    path.Clear();
                                    path.Add(i);
                                    path.Add(i);
                                    path.Add(j);
                                    path.Add(j);
                                    path.Add(q);
                                    optCost = cost;
                                    
                                }
                            }
                        }
                    }
                }
            }

           // E = path;
            DA.SetDataList(0, changeLines(linesAtNodes, path, angles1Single, angles1Rel, angles2Single, angles2Rel, angles20_total, plane, plane2));
            DA.SetData(1, optCost); 


        }


        // </Custom additional code> 
    



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
public override Guid ComponentGuid => new Guid("622d4c1d-39a6-4ad4-81d7-6d9d3f8c7398");
    }
}