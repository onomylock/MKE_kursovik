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
            //return Bound.R * Bound.R + Bound.Z + t;
            return Bound.R;
            //return Bound.R * Bound.R * Bound.R * Bound.R;
            //return Bound.R + t * t;
            //return Math.Cos(Bound.R * Bound.R + 3 * Bound.Z) - Math.Sin(t);
            //return t * t;
            //return 5;
        }

        public S1(int NumVertex, int Side)
        {
            this.NumVertex = NumVertex;
            this.Side = Side;
        }
    }

    public class S2
    {

    }
}
