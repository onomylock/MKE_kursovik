using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MKE_kursovik
{
    public interface IElement
    {
        public IEnumerable<int> VertexArr { get; }
        public int NMat { get; }
        //double Sigma { get; set; }
        //double Mu { get; set; }
    }

    public class Element : IElement
    {
        public Element(int[] elem)
        {
            VertexArr = elem;
        }

        public int NMat { get; set; }
        public IEnumerable<int> VertexArr { get; }
    }

    public class Vertex
    {
        public double R { get; set; }
        public double Z { get; set; }
        public int NumOfFun { get; set; }

        public Vertex(double r, double z, int num)
        {
            R = r;
            Z = z;
            NumOfFun = num;
        }

        public Vertex(double r, double z)
        {
            R = r;
            Z = z;
            NumOfFun = 0;
        }
    }

    public interface IGlobalMatrix
    {
        public double[] ggl { get; set; }
        public double[] ggu { get; set; }
        public double[] di { get; set; }
        public double[,] GenLocal(double rp, Parameter param, double Hz, double Hr);
    }

    public class Parameter
    {
        public double Sigma { get; }
        public double Mu { get; }

        public Parameter(double sigma, double mu)
        {
            Sigma = sigma;
            Mu = 1.0 / (mu * 4 * Math.PI * 10e-7);
            //Mu = mu;
        }
    }

    public class Function
    {
        public double fun(Vertex rz, double t, Parameter param)
        {
			//return 5 / (rz.R * rz.R);   //1 test
			//return 0; //2 test
			//return -3; //3 test
			//return -8 * rz.R; //4 test
			//return 1; //5 test
			//return 2 * t + t * t / (rz.R * rz.R); //6 test

			switch (rz.NumOfFun)
			{
				case 1: return param.Mu / (rz.R * rz.R);
				//case 1: return param.Mu;
				case 2: return 0;
				default: return 0;
			}
		}

		public double AzTrue(Vertex rz, double t)
        {
			//return 5;   //1 test
			//return rz.R; //2 test
			//return rz.R * rz.R; //3 test
			//return rz.R * rz.R * rz.R; //4 test
			//return t; //5 test
			//return t * t; //6 test

			switch (rz.NumOfFun)
			{
				case 1: return 1;
				//case 1: return rz.R * rz.R + 4 * t;
				case 2: return 0;
				default: return 0;
			}
		}
    }

    public class GlobalVectorB
    {
        public double[] B { get; set; }

        public void GenGlobalB(IGlobalMatrix M, List<Vertex> RZ, List<Element> Elements, int[] ia, int[] ja, double[] QLast, double tNow, double th, List<Parameter> Params)
        {
            B = new double[RZ.Count];
            double[] Q_tmp = new double[RZ.Count];
            Function function = new Function();


            foreach (var elem in Elements)
            {
                int[] elemArr = elem.VertexArr.ToArray();
                double hz = RZ[elemArr[3]].Z - RZ[elemArr[0]].Z;
                double hr = RZ[elemArr[1]].R - RZ[elemArr[0]].R;
                double[] LocalF = new double[4];
                double[] LocalQ = new double[4];

                LocalF[0] = function.fun(RZ[elemArr[0]], tNow, Params[elem.NMat]);
                LocalF[1] = function.fun(RZ[elemArr[1]], tNow, Params[elem.NMat]);
                LocalF[2] = function.fun(RZ[elemArr[2]], tNow, Params[elem.NMat]);
                LocalF[3] = function.fun(RZ[elemArr[3]], tNow, Params[elem.NMat]);

                double[] LocalB = GenLocal(LocalF, RZ[elemArr[0]].R, hr, hz);
                for (int i = 0; i < elemArr.Length; i++)
                {
                    LocalQ[i] = QLast[elemArr[i]];
                }
                LocalQ = GenLocal(LocalQ, RZ[elemArr[0]].R, hr, hz);

                for (int i = 0; i < elemArr.Length; i++)
                {
                    B[elemArr[i]] += LocalB[i] + Params[elem.NMat].Sigma * LocalQ[i] / th;
                    //B[elemArr[i]] += LocalB[i];
                    //Q_tmp[elemArr[i]] += LocalQ[i];
                }
            }


            //for (int i = 0; i < B.Length; i++)
            //{
            //	B[i] += Q_tmp[i] / th;
            //}
            //double[] B2 = new double[B.Length];
            //for (int i = 0; i < B.Length; i++)
            //{
            //	B2[i] = M.di[i] * QLast[i];
            //	for (int j = ia[i]; j < ia[i + 1]; j++)
            //	{
            //		int k = ja[j];
            //		B2[i] += M.ggl[j] * QLast[k];
            //		B2[k] += M.ggu[j] * QLast[i];
            //	}
            //}

            //for (int i = 0; i < B.Length; i++)
            //{
            //	B[i] += B2[i] / th;
            //}
        }

        private double[] GenLocal(double[] LocalVec, double rp, double Hr, double Hz)
        {
            double[] LocalQ = new double[4];
            double b1 = Hz / 36.0,
                    b2 = Hr * rp,
                    //b2 = Hr,
                    b3 = Math.Pow(Hr, 2);

            LocalQ[0] = LocalVec[0] * b1 * (4 * b2 + b3) +
                        LocalVec[1] * b1 * (2 * b2 + b3) +
                        LocalVec[2] * b1 * (4 * b2 + b3) / 2.0 +
                        LocalVec[3] * b1 * (2 * b2 + b3) / 2.0;

            LocalQ[1] = LocalVec[0] * b1 * (2 * b2 + b3) +
                        LocalVec[1] * b1 * (4 * b2 + 3 * b3) +
                        LocalVec[2] * b1 * (2 * b2 + b3) / 2.0 +
                        LocalVec[3] * b1 * (4 * b2 + 3 * b3) / 2.0;

            LocalQ[2] = LocalVec[0] * b1 * (4 * b2 + b3) / 2.0 +
                        LocalVec[1] * b1 * (2 * b2 + b3) / 2.0 +
                        LocalVec[2] * b1 * (4 * b2 + b3) +
                        LocalVec[3] * b1 * (2 * b2 + b3);

            LocalQ[3] = LocalVec[0] * b1 * (2 * b2 + b3) / 2.0 +
                        LocalVec[1] * b1 * (4 * b2 + 3 * b3) / 2.0 +
                        LocalVec[2] * b1 * (2 * b2 + b3) +
                        LocalVec[3] * b1 * (4 * b2 + 3 * b3);

            return LocalQ;
        }
    }

    public class GenGlobalMatrix
    {
        public void AddToGlobal(IGlobalMatrix global, List<Vertex> RZ, List<Element> Elements, List<Parameter> Params, int[] ia, int[] ja)
        {
            foreach (var elem in Elements)
            {
                int[] elemArr = elem.VertexArr.ToArray();
                double hz = RZ[elemArr[2]].Z - RZ[elemArr[0]].Z;
                double hr = RZ[elemArr[1]].R - RZ[elemArr[0]].R;
                double[,] Local = global.GenLocal(RZ[elemArr[0]].R, Params[elem.NMat], hz, hr);

                for (int i = 0; i < elemArr.Length; i++)
                {
                    global.di[elemArr[i]] += Local[i, i];
                    int tmp = ia[elemArr[i]];
                    for (int j = 0; j < i; j++)
                    {
                        for (int k = tmp; k < ia[elemArr[i] + 1]; k++)
                        {
                            if (ja[k] == elemArr[j])
                            {
                                global.ggu[k] += Local[j, i];
                                global.ggl[k] += Local[i, j];
                                k++;
                                break;
                            }
                        }

                    }
                }
            }
        }
    }
    
    public class GlobalMatrixG : IGlobalMatrix
    {
        public double[] ggl { get; set; }
        public double[] ggu { get; set; }
        public double[] di { get; set; }

        public GlobalMatrixG(int iaLen, int jaLen)
        {
            ggl = new double[jaLen];
            ggu = new double[jaLen];
            di = new double[iaLen - 1];
        }

        public double[,] GenLocal(double rp, Parameter param, double Hz, double Hr)
        {
            double[,] LocalG = new double[4, 4];
            //double a1 = (Hz * rp) / ( 6 * Hr),
            //        a2 = Hz / (param.Mu * 12),
            //        a3 = (Hr * rp) / (param.Mu * 6 * Hz),
            //        a4 = (Hr * Hr) / (param.Mu * 12 * Hz);
            double a1 = (param.Mu * Hz * rp) / (6 * Hr),
                    a2 = (param.Mu * Hz) / (12),
                    a3 = (param.Mu * Hr * rp) / (6 * Hz),
                    a4 = (param.Mu * Hr) / (12 * Hz);
            LocalG[0, 0] = 2 * a1 + 2 * a2 + 2 * a3 + 1 * a4;
            LocalG[0, 1] = -2 * a1 - 2 * a2 + 1 * a3 + 1 * a4;
            LocalG[0, 2] = 1 * a1 + 1 * a2 - 2 * a3 - 1 * a4;
            LocalG[0, 3] = -1 * a1 - 1 * a2 - 1 * a3 - 1 * a4;

            LocalG[1, 0] = -2 * a1 - 2 * a2 + 1 * a3 + 1 * a4;
            LocalG[1, 1] = 2 * a1 + 2 * a2 + 2 * a3 + 3 * a4;
            LocalG[1, 2] = -1 * a1 - 1 * a2 - 1 * a3 - 1 * a4;
            LocalG[1, 3] = 1 * a1 + 1 * a2 - 2 * a3 - 3 * a4;

            LocalG[2, 0] = 1 * a1 + 1 * a2 - 2 * a3 - 1 * a4;
            LocalG[2, 1] = -1 * a1 - 1 * a2 - 1 * a3 - 1 * a4;
            LocalG[2, 2] = 2 * a1 + 2 * a2 + 2 * a3 + 1 * a4;
            LocalG[2, 3] = -2 * a1 - 2 * a2 + 1 * a3 + 1 * a4;

            LocalG[3, 0] = -1 * a1 - 1 * a2 - 1 * a3 - 1 * a4;
            LocalG[3, 1] = 1 * a1 + 1 * a2 - 2 * a3 - 3 * a4;
            LocalG[3, 2] = -2 * a1 - 2 * a2 + 1 * a3 + 1 * a4;
            LocalG[3, 3] = 2 * a1 + 2 * a2 + 2 * a3 + 3 * a4;

            return LocalG;
        }
    }

    public class GlobalMatrixM0 : IGlobalMatrix
    {
        public double[] ggl { get; set; }
        public double[] ggu { get; set; }
        public double[] di { get; set; }

        public GlobalMatrixM0(int iaLen, int jaLen)
        {
            ggl = new double[jaLen];
            ggu = new double[jaLen];
            di = new double[iaLen - 1];
        }

        public double[,] GenLocal(double rp, Parameter param, double Hz, double Hr)
        {
            double[,] LocalM = new double[4, 4];
            double[,] LocalMr = new double[2, 2];
            double[,] LocalMz = new double[2, 2];
            double dr = rp / Hr,
                    log = Math.Log(1 + 1.0 / dr);

            LocalMr[0, 0] = log * Math.Pow((1 + dr), 2) - dr - 3.0 / 2.0;
            LocalMr[0, 1] = -log * dr * (1 + dr) + dr + 1.0 / 2.0;

            LocalMr[1, 0] = -log * dr * (1 + dr) + dr + 1.0 / 2.0;
            LocalMr[1, 1] = log * Math.Pow(dr, 2) - dr + 1.0 / 2.0;

            LocalMz[0, 0] = Hz / 3.0;
            LocalMz[0, 1] = Hz / 6.0;

            LocalMz[1, 0] = Hz / 6.0;
            LocalMz[1, 1] = Hz / 3.0;

            LocalM[0, 0] = param.Mu * LocalMr[0, 0] * LocalMz[0, 0];
            LocalM[0, 1] = param.Mu * LocalMr[0, 1] * LocalMz[0, 0];
            LocalM[0, 2] = param.Mu * LocalMr[0, 0] * LocalMz[0, 1];
            LocalM[0, 3] = param.Mu * LocalMr[0, 1] * LocalMz[0, 1];

            LocalM[1, 0] = param.Mu * LocalMr[0, 1] * LocalMz[0, 0];
            LocalM[1, 1] = param.Mu * LocalMr[1, 1] * LocalMz[0, 0];
            LocalM[1, 2] = param.Mu * LocalMr[1, 0] * LocalMz[0, 1];
            LocalM[1, 3] = param.Mu * LocalMr[1, 1] * LocalMz[0, 1];

            LocalM[2, 0] = param.Mu * LocalMr[0, 0] * LocalMz[0, 1];
            LocalM[2, 1] = param.Mu * LocalMr[1, 0] * LocalMz[0, 1];
            LocalM[2, 2] = param.Mu * LocalMr[0, 0] * LocalMz[1, 1];
            LocalM[2, 3] = param.Mu * LocalMr[0, 1] * LocalMz[1, 1];

            LocalM[3, 0] = param.Mu * LocalMr[0, 1] * LocalMz[0, 1];
            LocalM[3, 1] = param.Mu * LocalMr[1, 1] * LocalMz[0, 1];
            LocalM[3, 2] = param.Mu * LocalMr[0, 1] * LocalMz[1, 1];
            LocalM[3, 3] = param.Mu * LocalMr[1, 1] * LocalMz[1, 1];

            return LocalM;
        }
    }
    //public double[,] GenLocal(double rp, Parameter param, double Hz, double Hr)
    //{
    //    double[,] LocalM = new double[4, 4];
    //    double rp1 = rp + Hr,
    //            log = Math.Log(Math.Abs(rp1 / rp));
    //    LocalM[0, 0] = param.Mu * Hz * (Math.Pow(rp1, 2) * log - 2 * (Math.Pow(rp1, 2) - rp * rp1) + (Math.Pow(rp1, 2) - Math.Pow(rp, 2)) / 2) / (3 * Math.Pow(Hr, 2));
    //    LocalM[0, 1] = param.Mu * Hz * ((Math.Pow(rp1, 2) - Math.Pow(rp, 2)) / 2 - 2 * rp1 * rp * log) / (3 * Math.Pow(Hr, 2));
    //    LocalM[0, 2] = param.Mu * Hz * (Math.Pow(rp1, 2) * log - 2 * (Math.Pow(rp1, 2) - rp * rp1) + (Math.Pow(rp1, 2) - Math.Pow(rp, 2)) / 2) / (6 * Math.Pow(Hr, 2));
    //    LocalM[0, 3] = param.Mu * Hz * ((Math.Pow(rp1, 2) - Math.Pow(rp, 2)) / 2 - 2 * rp1 * rp * log) / (6 * Math.Pow(Hr, 2));

    //    LocalM[1, 0] = LocalM[0, 1];
    //    LocalM[1, 1] = param.Mu * Hz * (Math.Pow(rp, 2) * log - 2 * (rp * rp1 - Math.Pow(rp, 2)) + (Math.Pow(rp1, 2) - Math.Pow(rp, 2)) / 2) / (3 * Math.Pow(Hr, 2));
    //    LocalM[1, 2] = LocalM[0, 3];
    //    LocalM[1, 3] = param.Mu * Hz * (Math.Pow(rp, 2) * log - 2 * (rp * rp1 - Math.Pow(rp, 2)) + (Math.Pow(rp1, 2) - Math.Pow(rp, 2)) / 2) / (6 * Math.Pow(Hr, 2));

    //    LocalM[2, 0] = LocalM[0, 2];
    //    LocalM[2, 1] = LocalM[1, 2];
    //    LocalM[2, 2] = param.Mu * Hz * (Math.Pow(rp1, 2) * log - 2 * (Math.Pow(rp1, 2) - rp * rp1) + (Math.Pow(rp1, 2) - Math.Pow(rp, 2)) / 2) / (3 * Math.Pow(Hr, 2));
    //    LocalM[2, 3] = param.Mu * Hz * ((Math.Pow(rp1, 2) - Math.Pow(rp, 2)) / 2 - 2 * rp1 * rp * log) / (3 * Math.Pow(Hr, 2));

    //    LocalM[3, 0] = LocalM[3, 0];
    //    LocalM[3, 1] = LocalM[1, 3];
    //    LocalM[3, 2] = LocalM[2, 3];
    //    LocalM[3, 3] = param.Mu * Hz * (Math.Pow(rp, 2) * log - 2 * (rp * rp1 - Math.Pow(rp, 2)) + (Math.Pow(rp1, 2) - Math.Pow(rp, 2)) / 2) / (3 * Math.Pow(Hr, 2));

    //    return LocalM;
    //}

    // public double[,] GenLocal(double rp, Parameter param, double Hz, double Hr)
    // {
    //     double[,] LocalM = new double[4, 4];
    //     double[,] LocalMr = new double[2, 2];
    //     double[,] LocalMz = new double[2, 2];
    //     double dr = rp / Hr,
    //             log = Math.Log(1.0 + 1.0 / dr);

    //     LocalMr[0, 0] = log * Math.Pow((1 + dr), 2) - dr - 3.0 / 2.0;
    //     LocalMr[0, 1] = -log * dr * (1 + dr) + dr + 1.0 / 2.0;

    //     LocalMr[1, 0] = -log * dr * (1 + dr) + dr + 1.0 / 2.0;
    //     LocalMr[1, 1] = log * Math.Pow(dr, 2) - dr + 1.0 / 2.0;

    //     LocalMz[0, 0] = Hz / 3.0;
    //     LocalMz[0, 1] = Hz / 6.0;

    //     LocalMz[1, 0] = Hz / 6.0;
    //     LocalMz[1, 1] = Hz / 3.0;

    //     LocalM[0, 0] = param.Mu * LocalMr[0, 0] * LocalMz[0, 0];
    //     LocalM[0, 1] = param.Mu * LocalMr[0, 1] * LocalMz[0, 0];
    //     LocalM[0, 2] = param.Mu * LocalMr[0, 0] * LocalMz[0, 1];
    //     LocalM[0, 3] = param.Mu * LocalMr[0, 1] * LocalMz[0, 1];

    //     LocalM[1, 0] = param.Mu * LocalMr[0, 1] * LocalMz[0, 0];
    //     LocalM[1, 1] = param.Mu * LocalMr[1, 1] * LocalMz[0, 0];
    //     LocalM[1, 2] = param.Mu * LocalMr[1, 0] * LocalMz[0, 1];
    //     LocalM[1, 3] = param.Mu * LocalMr[1, 1] * LocalMz[0, 1];

    //     LocalM[2, 0] = param.Mu * LocalMr[0, 0] * LocalMz[0, 1];
    //     LocalM[2, 1] = param.Mu * LocalMr[1, 0] * LocalMz[0, 1];
    //     LocalM[2, 2] = param.Mu * LocalMr[0, 0] * LocalMz[1, 1];
    //     LocalM[2, 3] = param.Mu * LocalMr[0, 1] * LocalMz[1, 1];

    //     LocalM[3, 0] = param.Mu * LocalMr[0, 1] * LocalMz[0, 1];
    //     LocalM[3, 1] = param.Mu * LocalMr[1, 1] * LocalMz[0, 1];
    //     LocalM[3, 2] = param.Mu * LocalMr[0, 1] * LocalMz[1, 1];
    //     LocalM[3, 3] = param.Mu * LocalMr[1, 1] * LocalMz[1, 1];

    //     return LocalM;
    // }


    public class GlobalMatrixM : IGlobalMatrix
    {
        public double[] ggl { get; set; }
        public double[] ggu { get; set; }
        public double[] di { get; set; }

        public GlobalMatrixM(int iaLen, int jaLen)
        {
            ggl = new double[jaLen];
            ggu = new double[jaLen];
            di = new double[iaLen - 1];
        }

        public double[,] GenLocal(double rp, Parameter param, double Hz, double Hr)
        {
            double[,] LocalM = new double[4, 4];
            double b1 = param.Sigma * Hz / 36.0,
                    b2 = Hr * rp,
                    //b2 = Hr,
                    b3 = Math.Pow(Hr, 2);
            LocalM[0, 0] = b1 * (4 * b2 + b3);
            LocalM[0, 1] = b1 * (2 * b2 + b3);
            LocalM[0, 2] = b1 * (4 * b2 + b3) / 2.0;
            LocalM[0, 3] = b1 * (2 * b2 + b3) / 2.0;

            LocalM[1, 0] = b1 * (2 * b2 + b3);
            LocalM[1, 1] = b1 * (4 * b2 + 3 * b3);
            LocalM[1, 2] = b1 * (2 * b2 + b3) / 2.0;
            LocalM[1, 3] = b1 * (4 * b2 + 3 * b3) / 2.0;

            LocalM[2, 0] = b1 * (4 * b2 + b3) / 2.0;
            LocalM[2, 1] = b1 * (2 * b2 + b3) / 2.0;
            LocalM[2, 2] = b1 * (4 * b2 + b3);
            LocalM[2, 3] = b1 * (2 * b2 + b3);

            LocalM[3, 0] = b1 * (2 * b2 + b3) / 2.0;
            LocalM[3, 1] = b1 * (4 * b2 + 3 * b3) / 2.0;
            LocalM[3, 2] = b1 * (2 * b2 + b3);
            LocalM[3, 3] = b1 * (4 * b2 + 3 * b3);

            return LocalM;
        }
    }

    public class GlobalMatrixA
    {
        public double[] ggl { get; set; }
        public double[] ggu { get; set; }
        public double[] di { get; set; }

        public GlobalMatrixA(int iaLen, int jaLen)
        {
            ggl = new double[jaLen];
            ggu = new double[jaLen];
            di = new double[iaLen - 1];
        }

        public void GenGolbalMatrixA(IGlobalMatrix G, IGlobalMatrix M, IGlobalMatrix M0, double Th)
        {
			for (int i = 0; i < ggl.Length; i++)
			{
				ggl[i] = G.ggl[i] + M.ggl[i] / Th + M0.ggl[i];
				ggu[i] = G.ggu[i] + M.ggu[i] / Th + M0.ggu[i];
			}

			for (int i = 0; i < di.Length; i++)
			{
				di[i] = G.di[i] + M.di[i] / Th + M0.di[i];
			}

			//for (int i = 0; i < ggl.Length; i++)
			//{
			//	ggl[i] = G.ggl[i] + M.ggl[i] / Th;
			//	ggu[i] = G.ggu[i] + M.ggu[i] / Th;
			//}

			//for (int i = 0; i < di.Length; i++)
			//{
			//	di[i] = G.di[i] + M.di[i] / Th;
			//}
		}
    }

    public class Mesh
    {
        readonly IEnumerable<IElement> Elements;
        readonly IList<Vertex> rz;

        public Mesh(List<Element> Elements, List<Vertex> rz)
        {
            this.Elements = Elements;
            this.rz = rz;
        }

        public (int[] ia, int[] ja) BuildPotrait()
        {
            var map = new SortedSet<int>[rz.Count];
            for (int i = 0; i < map.Length; i++)
            {
                map[i] = new SortedSet<int>();
            }

            foreach (var element in Elements)
            {
                foreach (var i in element.VertexArr)
                {
                    foreach (var j in element.VertexArr)
                        if (i > j) map[i].Add(j);
                }
            }

            var ia = new int[map.Length + 1];
            ia[0] = 0;
            for (int i = 0; i < map.Length; i++)
            {
                ia[i + 1] = ia[i] + map[i].Count;
            }

            var ja = new int[ia[ia.Length - 1]];
            for (int i = 0; i < map.Length; i++)
            {
                var jind = map[i].ToArray();
                for (int j = 0; j < jind.Length; j++)
                {
                    ja[ia[i] + j] = jind[j];
                }
            }
            return (ia, ja);
        }
    }
}
