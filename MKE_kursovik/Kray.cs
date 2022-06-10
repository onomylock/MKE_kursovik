using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MKE_kursovik
{
    public class S1
    {
        public int NumVertex { get; }
        public int Side { get; }
        public double S1Fun(Vertex Bound, double t, int side)
        {
            //switch (side)
            //{
            //    case 0: return 5 * r * r; //bottom
            //    case 1: return 5 * r * r; //up
            //    case 2: return 5 * r * r; //left
            //    case 3: return 5 * r * r; //right
            //    default: return 0;
            //}
            //return 5; //1 test
            //return Bound.R; //2 test
            //return Bound.R * Bound.R; // test 3
            //return Bound.R * Bound.R * Bound.R; //4 test
            //return t; //5 test        
            //return t * t; //6 test
            return 0;
        }

        public S1(int NumVertex, int Side)
        {
            this.NumVertex = NumVertex;
            this.Side = Side;
        }
    }

    public class S2
    {
        public S2(int NumVertex1, int NumVertex2, int Side)
        {
            this.NumVertex1 = NumVertex1;
            this.NumVertex2 = NumVertex2;
            this.Side = Side;
        }

        public int Side { get; }
        public int NumVertex1 { get; }
        public int NumVertex2 { get; }

        private double TettaFun(Vertex P)
        {
            return 0;
        }

        public double[] GenVectorBS2(Vertex P1, Vertex P2, double t, double Hz)
        {
            double[] B_S2 = new double[2];

            return B_S2;
        }
    }
}
