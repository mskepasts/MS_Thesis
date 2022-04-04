using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using Rhino;
using Rhino.Geometry;

namespace oneNode
{
    public class oneNodeComponent : GH_Component
    {
        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public oneNodeComponent()
          : base("oneNode", "1N",
            "Description",
            "Category", "Subcategory")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("angleTop", "angleTop", "Top Angles as defined in the paper", GH_ParamAccess.list);
            pManager.AddNumberParameter("angleSide", "angleTop", "Side anlges as defined in the paper", GH_ParamAccess.list);
            pManager.AddLineParameter("linesAtNodes", "linesAtNodes", "TOO FILL", GH_ParamAccess.list);
            pManager.AddPlaneParameter("plane", "plane", "TOO FILL", GH_ParamAccess.item);
            pManager.AddPlaneParameter("plane2", "plane2", "TOO FILL", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.Register_DoubleParam("Cost", "C", "Cost of the function");
            pManager.Register_LineParam("New Lines", "lns", "Lines found from the function");
        }

        public static class Globals
        {
            public static Int32 flipAngle = 0;
        }

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
        private List<double> angle2Projected(List<Line> lns, Plane pl, List<Plane> plns, out List<double> anglesOriginal)
        {
            List<double> angles = new List<double>();
            anglesOriginal = new List<double>();
            double angleStart = 0.0;
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


        //updates lines to current state
        private List<Line> changeLines(List<Line> lns, List<int> path, List<List<double>> angles1Rel, List<List<double>> angles2Rel, List<double> angles2Start, Plane pl, List<Plane> plns2)
        {
            List<Line> lns_Updated = new List<Line>();
            List<Line> lns_Updated_length = new List<Line>();
            var xform = Rhino.Geometry.Transform.Rotation(RhinoMath.ToRadians(360 - angles1Rel[path[0]][path[1]]), pl.ZAxis, pl.Origin); ;
            var xform2 = Rhino.Geometry.Transform.Rotation(-1 * RhinoMath.ToRadians(angles2Start[path[0]]), plns2[0].ZAxis, plns2[0].Origin);
            Line holderLn = new Line(lns[0].From, pl.ClosestPoint(lns[0].To));
            holderLn.Transform(xform2);
            lns_Updated.Add(holderLn);

            for (int i = 0; i < lns.Count - 1; i++)
            {
                xform = Rhino.Geometry.Transform.Rotation(RhinoMath.ToRadians(360 - angles1Rel[path[i]][path[i + 1]]), pl.ZAxis, pl.Origin);
                if (Globals.flipAngle == 0)
                    xform = Rhino.Geometry.Transform.Rotation(RhinoMath.ToRadians(angles1Rel[path[i]][path[i + 1]]), pl.ZAxis, pl.Origin);


                Line holder = new Line(lns_Updated[lns_Updated.Count - 1].From, lns_Updated[lns_Updated.Count - 1].To);
                holder.Transform(xform);
                xform2 = Rhino.Geometry.Transform.Rotation(-1 * RhinoMath.ToRadians(angles2Rel[path[i]][path[i + 1]]), plns2[i + 1].ZAxis, plns2[i + 1].Origin);
                holder.Transform(xform2);


                lns_Updated.Add(holder);

            }
            for (int i = 0; i < lns.Count; i++)
            {
                lns_Updated_length.Add(new Line(lns_Updated[i].From, lns_Updated[i].Direction, lns[i].Length));

            }
            return lns_Updated_length;
        }

        //Set up angle arrays

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
           
            List<double> angle10 = new List<double>();
            if (!DA.GetDataList(0, angle10)) { return; }
            List<double> angle20 = new List<double>();
            if (!DA.GetDataList(1, angle20)) { return; }
            List<Line> linesAtNodes = new List<Line>();
            if (!DA.GetDataList(2, linesAtNodes)) { return; }
            Plane plane = new Plane();
            if (!DA.GetData(3, ref plane)) { return; }
            List<Plane> plane2 = new List<Plane>();
            if (!DA.GetDataList(4, plane2)) { return; }



            int numberOfNodes = angle10.Count;


            // Set up the angle arrays
            double angle = 360 / 6;
            List<double> angles10_total = new List<double>();
            List<double> angle11 = new List<double>();

            for (int i = 0; i < numberOfNodes; i++)
            {
                angles10_total.Add(angle10[i]);
                angle11.Add(angle - angle10[i]);
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

            for (int i = 0; i < angles10_total.Count; i++)
            {
                List<double> sublist = new List<double>();
                angles1Single.Add(angle - angles10_total[i] - angles11_total[i]);
                for (int j = 0; j < angles10_total.Count; j++)
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
                angles20_total.Add(-angle20[i]);
            }


            List<double> angles2Single = new List<double>();
            List<List<double>> angles2Rel = new List<List<double>>();
            for (int i = 0; i < numberOfNodes * 2; i++)
            {
                List<double> sublist = new List<double>();
                for (int j = 0; j < numberOfNodes * 2; j++)
                {
                    sublist.Add(angles20_total[j] - angles20_total[i]);
                }
                angles2Rel.Add(sublist);
            }
            





            List<double> angles1 = angle1Projected(linesAtNodes, plane);
            List<double> angles2Start;
            List<double> angles2 = angle2Projected(linesAtNodes, plane, plane2, out angles2Start);

/*            B = angles2;
            C = angles20_total;*/
            List<int> path = new List<int>();
            int count = angles1Single.Count;
            double angleCheck = 0.0;
            double cost = 0.0;
            double optCost = Double.PositiveInfinity;
            for (int i = 0; i < count; i++)
            {
                cost = 0.0;
                angleCheck = angles20_total[i];
                //double cost1 = Math.Abs(angles20_total[i] - angle2Start[0]);
                double cost1 = Math.Abs(angles2Start[0] - angleCheck);
                cost += cost1;
                if (cost < optCost)
                    for (int j = 0; j < count; j++)
                    {
                        cost = cost1;
                        angleCheck = angles20_total[i];
                        angleCheck += angles2Rel[i][j];
                        double cost2 = Math.Abs(angleCheck - angles2Start[1]) + Math.Abs(angles1Rel[i][j] - angles1[0]);
                        //double cost2 = Math.Abs(angles2Rel[i][j] - angles2[0]) + Math.Abs(angles1Rel[i][j] - angles1[0]);
                        //double cost2 = Math.Abs(angles1Rel[i][j] - angles1[0]);
                        cost += cost2;
                        if (cost < optCost && cost < 10.00)
                            for (int k = 0; k < count; k++)
                            {
                                angleCheck = angles20_total[i] + angles2Rel[i][j];
                                angleCheck += angles2Rel[j][k];
                                cost = cost1 + cost2;
                                double cost3 = Math.Abs(angleCheck - angles2Start[2]) + Math.Abs(angles1Rel[j][k] - angles1[1]);

                                //double cost3 = Math.Abs(angles2Rel[j][k] - angles2[1]) + Math.Abs(angles1Rel[j][k] - angles1[1]);
                                //double cost3 = Math.Abs(angles1Rel[j][k] - angles1[1]);
                                cost += cost3;
                                if (cost < optCost && cost < 10.00)
                                    for (int l = 0; l < count; l++)
                                    {
                                        cost = cost1 + cost2 + cost3;
                                        angleCheck = angles20_total[i] + angles2Rel[i][j] + angles2Rel[j][k];
                                        angleCheck += angles2Rel[k][l];
                                        double cost4 = Math.Abs(angleCheck - angles2Start[3]) + Math.Abs(angles1Rel[k][l] - angles1[2]);
                                        //double cost4 = Math.Abs(angles1Rel[k][l] - angles1[2]);
                                        //double cost4 = Math.Abs(angles2Rel[k][l] - angles2[2]) + Math.Abs(angles1Rel[k][l] - angles1[2]);
                                        cost += cost4;
                                        if (cost < optCost && cost < 10.00)
                                            for (int m = 0; m < count; m++)
                                            {
                                                cost = cost1 + cost2 + cost3 + cost4;
                                                angleCheck = angles20_total[i] + angles2Rel[i][j] + angles2Rel[j][k] + angles2Rel[k][l];
                                                angleCheck += angles2Rel[l][m];
                                                double cost5 = Math.Abs(angleCheck - angles2Start[4]) + Math.Abs(angles1Rel[l][m] - angles1[3]);
                                                //double cost5 = Math.Abs(angles2Rel[l][m] - angles2[3]) + Math.Abs(angles1Rel[l][m] - angles1[3]);
                                                //double cost5 = Math.Abs(angles1Rel[l][m] - angles1[3]);
                                                cost += cost5;
                                                if (cost < optCost && cost < 10.00)
                                                    for (int n = 0; n < count; n++)
                                                    {
                                                        cost = cost1 + cost2 + cost3 + cost4 + cost5;
                                                        angleCheck = angles20_total[i] + angles2Rel[i][j] + angles2Rel[j][k] + angles2Rel[k][l] + angles2Rel[l][m];
                                                        angleCheck += angles2Rel[m][n];
                                                        double cost6 = Math.Abs(angleCheck - angles2Start[5]) + Math.Abs(angles1Rel[m][n] - angles1[4]) + Math.Abs(angles1Rel[n][i] - angles1[5]);
                                                        //double cost6 = Math.Abs(angles2Rel[m][n] - angles2[4]) + Math.Abs(angles1Rel[m][n] - angles1[4]) + Math.Abs(angles2Rel[n][i] - angles2[5]) + Math.Abs(angles1Rel[n][i] - angles1[5]);
                                                        //double cost6 = Math.Abs(angles1Rel[m][n] - angles1[4]) + Math.Abs(angles1Rel[n][i] - angles1[5]);

                                                        cost += cost6;

                                                        if (cost < optCost)
                                                        {
                                                            path.Clear();
                                                            path.Add(i);
                                                            path.Add(j);
                                                            path.Add(k);
                                                            path.Add(l);
                                                            path.Add(m);
                                                            path.Add(n);
                                                            optCost = cost;
                                                            
                                                        }

                                                    }
                                            }
                                    }
                            }
                    }
            }
            DA.SetData(0, optCost);
            //E = path;
            DA.SetDataList(1, changeLines(linesAtNodes, path, angles1Rel, angles2Rel, angles20_total, plane, plane2));
            List<double> check = new List<double>();
            for (int i = 0; i < path.Count - 1; i++)
            {
                check.Add(angles2Rel[path[i]][path[i + 1]]);
            }
            check.Add(angles2Rel[path[path.Count - 1]][path[0]]);
           /* F = check;*/
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
        public override Guid ComponentGuid => new Guid("21972360-E278-4D1B-B6DE-461BA81E0C26");
    }
}