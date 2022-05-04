using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace MKE_kursovik
{
    public class Solve
    {
        public List<double[]> Q { get; set; }
        GlobalMatrixM M;
        GlobalMatrixG G;
        GlobalMatrixA A;
        GlobalVectorB rightB;
        int[] ia, ja;
        int N;
        //IO io;

        public Solve(IO io)
        {
            //this.io = io;
            Mesh mesh = new Mesh(io.Elements, io.RZ);
            (ia, ja) = mesh.BuildPotrait();
            M = new GlobalMatrixM(ia.Length, ja.Length);
            G = new GlobalMatrixG(ia.Length, ja.Length);
            A = new GlobalMatrixA(ia.Length, ja.Length);
            GenGlobalMatrix gen = new GenGlobalMatrix();
            gen.AddToGlobal(G, io.RZ, io.Elements, io.Params, ia, ja);
            gen.AddToGlobal(M, io.RZ, io.Elements, io.Params, ia, ja);
            A.GenGolbalMatrixA(G, M, io.Time[1] - io.Time[0]);
            N = io.RZ.Count();
            Function function = new Function();

            Q = new List<double[]>();
            //A.GenGolbalMatrixA(G, M, io.Time[1] - io.Time[0]);
            //ResultVectors Q = new ResultVectors();
            double[] Q_tmp = new double[ia.Length - 1];
            for (int i = 0; i < Q_tmp.Length; i++)
            {
                Q_tmp[i] = function.AzTrue(io.RZ[i], io.Time[0]);
                //Q_tmp[i] = function.fun(io.RZ[i].R, io.RZ[i].Z, io.Time[0]);
            }
            Q.Add(Q_tmp);

            //foreach (var s1 in io.Bound1)
            //{
            //    A.di[s1.NumVertex] = 1.0e50;
            //}

            for (int i = 1; i < io.Time.Count; i++)
            {
                rightB = new GlobalVectorB();
                A.GenGolbalMatrixA(G, M, io.Time[i] - io.Time[i - 1]);
                rightB.GenGlobalB(M, io.RZ, io.Elements, ia, ja, Q[i - 1], io.Time[i], io.Time[i] - io.Time[i - 1], io.Params);
                Q_tmp = new double[N];

                for (int j = 0; j < N; j++)
                {
                    Q_tmp[j] += M.di[j] * Q[i - 1][j] / (io.Time[i] - io.Time[i - 1]);
                    for (int k = ia[j]; k < ia[j + 1]; k++)
                    {
                        int l = ja[k];
                        Q_tmp[j] += M.ggl[k] * Q[i - 1][l] / (io.Time[i] - io.Time[i - 1]);
                        Q_tmp[l] += M.ggu[k] * Q[i - 1][j] / (io.Time[i] - io.Time[i - 1]);
                    }
                }
                for (int j = 0; j < N; j++)
                {
                    rightB.B[j] += Q_tmp[j];
                }

                foreach (var s2 in io.Bound2)
                {
                    double HZtmp = io.RZ[s2.NumVertex2].Z - io.RZ[s2.NumVertex1].Z;
                    double[] S2Tmp = s2.GenVectorBS2(io.RZ[s2.NumVertex1], io.RZ[s2.NumVertex2], io.Time[i], HZtmp);
                    rightB.B[s2.NumVertex1] += S2Tmp[0];
                    rightB.B[s2.NumVertex2] += S2Tmp[1];
                }

                foreach (var s1 in io.Bound1)
                {
                    A.di[s1.NumVertex] = 1.0e50;
                    rightB.B[s1.NumVertex] = 1.0e50 * s1.S1Fun(io.RZ[s1.NumVertex], io.Time[i], s1.Side);
                }
                Q.Add(LOS());
            }

            using (StreamWriter writer = new StreamWriter("output.txt", false))
            {
                //string text = "Afi      Afi(true)      Afi(err)";
                //writer.WriteLine(text);

                //for (int i = 1; i < io.Time.Count; i++)
                //{
                //    //writer.WriteLine("Time = " + io.Time[i].ToString());
                //    //double[] Az_true = new double[Q[i].Length];
                //    //double sum = 0;
                //    //for (int k = 0; k < io.Elements.Count; k++)
                //    //{
                //    //    for (int j = 0; j < 4; j++)
                //    //    {
                //    //        int l = io.Elements[k].VertexArr.ToArray()[j];
                //    //        Az_true[l] = function.AzTrue(io.RZ[l], io.Time[i]);
                //    //    }
                //    //}

                //    for (int j = 0; j < Q[i].Length; j++)
                //    {
                //        //double sum_tmp = Math.Abs(Q[i][j] - Az_true[j]);
                //        //sum += Math.Pow(sum_tmp, 2);
                //        writer.WriteLine(io.RZ[j].R.ToString() + "\t\t\t\t\t" + Q[i][j].ToString());
                //    }
                //    //sum = Math.Sqrt(sum);
                //    //sum = Math.Sqrt(sum) / (double)Q[i].Length;
                //    //writer.WriteLine("square error = " + sum.ToString());
                //}
                int start = io.NumR * (io.NumZDown - 1);
                for (int j = start + io.NumRDown; j < start + io.NumR; j++)
                //for (int j = 0; j < Q[1].Length; j++)
                {
                    writer.WriteLine(io.RZ[j].R.ToString() + "\t\t\t\t\t" + Q[2][j].ToString());
                    //writer.WriteLine(io.RZ[j].R.ToString() + "\t\t\t\t\t" + io.RZ[j].Z.ToString() + "\t\t\t\t\t" + Q[2][j].ToString());
                }

                writer.Close();
            }

            //using (StreamWriter writer = new StreamWriter("output.txt", false))
            //{
            //    //string text = "Afi      Afi(true)      Afi(err)";
            //    //writer.WriteLine(text);

            //    for (int i = 0; i < io.Time.Count; i++)
            //    {
            //        writer.WriteLine("Time = " + io.Time[i].ToString());
            //        double[] Az_true = new double[Q[i].Length];
            //        double sum = 0;
            //        for (int k = 0; k < io.Elements.Count; k++)
            //        {
            //            for (int j = 0; j < 4; j++)
            //            {
            //                int l = io.Elements[k].VertexArr.ToArray()[j];
            //                Az_true[l] = function.AzTrue(io.RZ[l], io.Time[i]);
            //            }
            //        }
            //        for (int j = 0; j < Q[i].Length; j++)
            //        {
            //            double sum_tmp = Math.Abs(Q[i][j] - Az_true[j]);
            //            sum += Math.Pow(sum_tmp, 2);
            //            //writer.WriteLine(Q[i][j].ToString() + "\t\t\t\t\t" + Az_true[j].ToString() + "\t\t\t\t\t" + sum_tmp.ToString());
            //        }
            //        //sum = Math.Sqrt(sum);
            //        sum = Math.Sqrt(sum) / (double)Q[i].Length;
            //        writer.WriteLine("square error = " + sum.ToString());
            //    }
            //    writer.Close();
            //}
        }

        private void LUsq(double[] ggl, double[] ggu, double[] di)
        {
            for (int i = 0; i < N; i++)
            {
                double sumdi = 0;
                int i0 = ia[i];
                int i1 = ia[i + 1];
                for (int k = i0; k < i1; k++) // к показывает сколько элементов мы обработали в i строке 
                {
                    int j = ja[k];// номер столбца  к-го элемента i строки 
                    double sumal = 0;
                    double sumau = 0;

                    int j0 = ia[j];
                    int j1 = ia[j + 1];

                    int ki = i0;
                    int kj = j0;
                    for (; ki < k && kj < j1;) // пока есть элементы которые предшествуют к-му элементу
                    {
                        int j_kl = ja[ki];
                        int j_ku = ja[kj];
                        if (j_kl == j_ku) // чтобы рассматриваемые элементы не были нулевыми 
                        {
                            sumal += ggu[kj] * ggl[ki];
                            sumau += ggu[ki] * ggl[kj];
                            ki++; kj++;
                        }
                        else if (j_kl < j_ku) ki++;
                        else kj++;
                    }

                    ggu[k] = (ggu[k] - sumau) / di[j];
                    ggl[k] = (ggl[k] - sumal) / di[j];
                    sumdi += ggl[k] * ggu[k];
                }
                di[i] = Math.Sqrt(di[i] - sumdi);
            }
        }

        private double[] MultMatrixVector(double[] param)
        {
            double[] res = new double[param.Count()];
            for (int i = 0; i < N; i++)
            {
                int i0 = ia[i];
                int i1 = ia[i + 1];
                res[i] = A.di[i] * param[i];

                for (int j = i0; j < i1; j++)
                {
                    int st_j = ja[j];
                    res[i] += param[st_j] * A.ggl[j];
                    res[st_j] += param[i] * A.ggu[j];
                }
            }
            return res;
        }

        private double[] Forward(double[] ggl, double[] di, double[] F)
        {
            double[] X = new double[N];
            for (int i = 0; i < N; i++)
            {
                double sum = 0;

                for (int j = ia[i]; j < ia[i + 1]; j++)
                {
                    int st_j = ja[j];
                    sum += ggl[j] * X[st_j];
                }
                X[i] = (F[i] - sum) / di[i];
            }
            return X;
        }

        private double[] Reverse(double[] ggu, double[] di, double[] F)
        {
            double[] X = new double[N];
            Array.Copy(F, X, F.Count());
            for (int i = N - 1; i >= 0; i--)
            {
                X[i] /= di[i];

                for (int j = ia[i + 1] - 1; j >= ia[i]; j--)
                {
                    int st_j = ja[j];
                    X[st_j] -= ggu[j] * X[i];
                }
            }
            return X;
        }

        private double ScalarProduct(double[] A1, double[] A2)
        {
            double sum = 0;

            for (int i = 0; i < A1.Count(); i++)
            {
                sum += A1[i] * A2[i];
            }

            return sum;
        }

        private double[] LOS()
        {
            double eps = 1e-60, norm, alfa = 0, beta = 0;
            int k = 0, max_iter = 1000;

            double[] mult1, mult2, R, Z, P, ggl, ggu, di, Q_tmp;
            mult1 = new double[N];
            mult2 = new double[N];
            Q_tmp = new double[N];
            Array.Fill(Q_tmp, 1);
            R = new double[N];
            Z = new double[N];
            P = new double[N];
            di = new double[N];
            ggu = new double[A.ggu.Count()];
            ggl = new double[A.ggl.Count()];
            Array.Copy(A.ggu, ggu, A.ggu.Count());
            Array.Copy(A.di, di, A.di.Count());
            Array.Copy(A.ggl, ggl, A.ggl.Count());

            //foreach (var s1 in io.Bound1)
            //{
            //    di[s1.NumVertex] = 1.0e50;
            //}

            LUsq(ggl, ggu, di);
            mult1 = MultMatrixVector(Q_tmp);
            for (int i = 0; i < N; i++) mult2[i] = rightB.B[i] - mult1[i];

            R = Forward(ggl, di, mult2);
            Z = Reverse(ggu, di, R);

            mult1 = MultMatrixVector(Z);
            P = Forward(ggl, di, mult1);
            norm = Math.Sqrt(ScalarProduct(R, R));

            for (; k < max_iter && norm >= eps; k++)
            {
                double scal_a = ScalarProduct(P, R);
                double scal_b = ScalarProduct(P, P);
                alfa = scal_a / scal_b;

                for (int i = 0; i < N; i++)
                {
                    Q_tmp[i] += alfa * Z[i];
                    R[i] -= alfa * P[i];
                }
                mult1 = Reverse(ggu, di, R);
                mult2 = MultMatrixVector(mult1);
                mult1 = Forward(ggl, di, mult2);

                scal_a = ScalarProduct(P, mult1);
                beta = -scal_a / scal_b;

                mult2 = Reverse(ggu, di, R);

                for (int i = 0; i < N; i++)
                {
                    Z[i] = beta * Z[i] + mult2[i];
                    P[i] = mult1[i] + beta * P[i];
                }
                norm = Math.Sqrt(ScalarProduct(R, R));

                //Console.WriteLine("Iteration: " + (k + 1).ToString() + "\tError: " + norm.ToString());
            }

            return Q_tmp;
        }

    }
}
