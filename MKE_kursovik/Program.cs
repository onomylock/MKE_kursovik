using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;

namespace MKE_kursovik
{

    public class Program
    {
        IO io;
        Solve solve;


        int[] findKE(Vertex Point)
        {
            if (io.Elements == null) return null;
            foreach (var elem in io.Elements)
            {
                int[] elemArr = elem.VertexArr.ToArray();
                if (io.RZ[elemArr[0]].R <= Point.R && Point.R <= io.RZ[elemArr[1]].R)
                    if (io.RZ[elemArr[0]].Z <= Point.Z && Point.Z <= io.RZ[elemArr[3]].Z)
                        return elemArr;
            }

            return null;
        }

        // double BPoint(Vertex Point, int layerNum)
        // {

        // }
        double APoint(Vertex Point, int layerNum)
        {
            int[] s = findKE(Point);

            if (s == null)
            {
                Console.WriteLine("Error: Wrong point");
                return -1;
            }

            double f1, f2, f3, f4;
            double q1, q2, q3, q4;
            double x1, x2, y1, y2;

            x1 = (io.RZ[s[1]].R - Point.R) / (io.RZ[s[1]].R - io.RZ[s[0]].R);
            y1 = (io.RZ[s[3]].Z - Point.Z) / (io.RZ[s[3]].Z - io.RZ[s[0]].Z);
            x2 = (Point.R - io.RZ[s[0]].R) / (io.RZ[s[1]].R - io.RZ[s[0]].R);
            y2 = (Point.Z - io.RZ[s[0]].Z) / (io.RZ[s[3]].Z - io.RZ[s[0]].Z);

            f1 = x1 * y1;
            f2 = x2 * y1;
            f3 = x1 * y2;
            f4 = x2 * y2;

            q1 = solve.Q[layerNum][s[0]];
            q2 = solve.Q[layerNum][s[1]];
            q3 = solve.Q[layerNum][s[2]];
            q4 = solve.Q[layerNum][s[3]];

            double v = q1 * f1 + q2 * f2 + q3 * f3 + q4 * f4;
            return v;
        }

        static void Main(string[] args)
        {
            Program program = new Program();
            
            program.io = new IO(true);
            program.solve = new(program.io);
            Console.WriteLine("y");
            // List<Vertex> Points = new List<Vertex>();
            // int layer;

            // Points.Add(new Vertex(1.3, 1.2));
            // Points.Add(new Vertex(1.7, 1.6));
            // Points.Add(new Vertex(1.1, 1.4));
            // Points.Add(new Vertex(1.4, 1.2));
            // Points.Add(new Vertex(1.1, 1.9));
            // Points.Add(new Vertex(1.0, 1.0));

            // Console.WriteLine("Write layer:");
            // layer = Convert.ToInt32(Console.ReadLine());


            // using (StreamWriter writer = new StreamWriter("Point.txt", false))
            // {
            //     double? Afi;
            //     writer.WriteLine("Layer = " + layer.ToString());
            //     foreach (var Point in Points)
            //     {
            //         Afi = program.APoint(Point, layer);
            //         writer.WriteLine("R = " + Point.R.ToString() + "\tZ = " + Point.Z.ToString() + "\t Function = " + Afi.ToString());
            //     }
            //     writer.Close();
            // }

        }
    }
}
