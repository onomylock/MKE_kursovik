using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace MKE_kursovik
{
    public class IO
    {
        public int NumR;
        public int NumZ;
        public int NumElem;
        public int NumPointSource;

        public int NumRDown;
        public int NumRUp;
        public int NumZDown;
        public int NumZUp;

        public double DiscRUp;
        public double DiscZUp;
        public double DiscRDown;
        public double DiscZDown;

        //public double Sigma;
        //public double Mu;
        private double TimeStart;
        private double TimeEnd;
        private double R0;
        private double R1;
        private double Z0;
        private double Z1;
        private double DiscR;
        private double DiscZ;
        private double TimeH;
        private double HrDown;
        private double HrUp;
        private double HzDown;
        private double HzUp;
        private double Hr;
        private double Hz;
        private int[] NumOfParameters;

        public List<S1> Bound1 { get; set; }

        public List<S2> Bound2 { get; set; }
        public List<Parameter> Params { get; set; }
        public List<double> Time { get; set; }
        public List<Vertex> RZ { get; set; }
        public List<Element> Elements { get; set; }

        public Vertex PointSource { get; set; }

        public IO(bool PointSourseBool)
        {
            if (PointSourseBool)
            {
                InputGridPointSource();
                GetGridPointSource();
            }
            else
            {
                InputGrid();
                GetGrid();
            }

            InputTime();
            GetTime();
            InputNumofPram();
            GetElements();
            InputBound1();
            InputBound2();
            InputPrams();
        }

        private void GetTime()
        {
            Time = new List<double>();
            double TimeTmp = TimeStart;
            int TimeIter = ((int)((int)(TimeEnd - TimeStart) / TimeH));

            for (int i = 0; i <= TimeIter; i++, TimeTmp += TimeH)
            {
                Time.Add(TimeTmp);
            }
        }

        private void GetGrid()
        {
            double HrTmp = Hr;
            double HzTmp = Hz;
            double Ri, Zi;
            RZ = new List<Vertex>();
            RZ.Add(new Vertex(R0, Z0));
            try
            {
                for (int i = 0; i < NumZ; i++)
                {
                    if (i == 0) Zi = Z0;
                    else if (i == NumZ - 1) Zi = Z1;
                    else Zi = RZ[RZ.Count - 1].Z + HzTmp;

                    for (int j = 0; j < NumR; j++)
                    {
                        if (j == 0 && i == 0) continue;
                        else if (j == 0) Ri = R0;
                        else if (j == NumR - 1) Ri = R1;
                        else Ri = RZ[RZ.Count - 1].R + HrTmp;

                        RZ.Add(new Vertex(Ri, Zi));
                        HrTmp *= DiscR;
                    }
                    HzTmp *= DiscZ;
                }
            }
            catch (Exception)
            {

                throw;
            }
        }

        private void GetGridPointSource()
        {
            double HrTmpDown = HrDown;
            double HrTmpUp = HrUp;
            double HzTmp = HzDown;
            double Ri, Zi;
            RZ = new List<Vertex>();
            RZ.Add(new Vertex(R0, Z0));
            try
            {
                for (int i = 0; i < NumZDown; i++)
                {
                    if (i == 0) Zi = Z0;
                    else if (i == NumZDown - 1) Zi = PointSource.Z;
                    else Zi = RZ[RZ.Count - 1].Z + HzTmp;
                    HrTmpUp = HrUp;
                    HrTmpDown = HrDown;
                    for (int j = 0; j < NumRDown; j++)
                    {
                        if (j == 0 && i == 0) continue;
                        else if (j == 0)
                        {
                            RZ.Add(new Vertex(R0, Zi));
                            continue;
                        }
                        else if (j == NumRDown - 1) Ri = PointSource.R;
                        else Ri = RZ[RZ.Count - 1].R + HrTmpDown;

                        RZ.Add(new Vertex(Ri, Zi));
                        HrTmpDown *= DiscRDown;
                    }

                    for (int j = 0; j < NumRUp; j++)
                    {
                        if (j == 0 && i == 0) continue;
                        else if (j == 0)
                        {
                            //RZ.Add(new Vertex(PointSource.R, Zi));
                            continue;
                        }
                        else if (j == NumRUp - 1) Ri = R1;
                        else Ri = RZ[RZ.Count - 1].R + HrTmpUp;

                        RZ.Add(new Vertex(Ri, Zi));
                        HrTmpUp *= DiscRUp;
                    }
                    if (i == 0) continue;
                    HzTmp *= DiscZDown;
                }

                HzTmp = HzUp;
                for (int i = 0; i < NumZUp; i++)
                {
                    if (i == 0) continue;
                    //if (i == 0) Zi = PointSource.Z;
                    //if (i == 0)
                    //{
                    //    HzTmp /= DiscZ;
                    //    continue;
                    //}
                    if (i == NumZUp - 1) Zi = Z1;
                    else Zi = RZ[RZ.Count - 1].Z + HzTmp;
                    HrTmpUp = HrUp;
                    HrTmpDown = HrDown;
                    for (int j = 0; j < NumRDown; j++)
                    {
                        //if (j == 0 && i == 0) continue;
                        if (j == 0)
                        {
                            RZ.Add(new Vertex(R0, Zi));
                            continue;
                        }
                        else if (j == NumRDown - 1) Ri = PointSource.R;
                        else Ri = RZ[RZ.Count - 1].R + HrTmpDown;

                        RZ.Add(new Vertex(Ri, Zi));
                        HrTmpDown *= DiscRDown;
                    }

                    for (int j = 0; j < NumRUp; j++)
                    {
                        //if (j == 0 && i == 0) continue;
                        if (j == 0) continue;
                        else if (j == NumRUp - 1) Ri = R1;
                        else Ri = RZ[RZ.Count - 1].R + HrTmpUp;

                        RZ.Add(new Vertex(Ri, Zi));
                        HrTmpUp *= DiscRUp;
                    }
                    HzTmp *= DiscZUp;
                }
                NumPointSource = NumR * (NumZDown - 1) - 1 + NumRDown;
                RZ[NumPointSource].NumOfFun = 1;
            }


            catch (Exception)
            {

                throw;
            }
        }

        private void InputBound1()
        {
            Bound1 = new List<S1>();
            for (int i = 0; i < NumR - 1; i++)
            {
                Bound1.Add(new S1(i, 0));
            }

            int num = NumR * NumZ;

            for (int i = num - 2; i > num - 1 - NumZ; i--)
            {
                Bound1.Add(new S1(i, 1));
            }

            for (int i = NumR - 1; i < num; i += NumR)
            {
                Bound1.Add(new S1(i, 2));
            }


            //string path = "Bound1.txt";
            //try
            //{
            //    using (StreamReader sr = new StreamReader(path))
            //    {
            //        int NumPeaks = int.Parse(sr.ReadLine());
            //        int NumVertex, Side;
            //        Bound1 = new List<S1>();
            //        for (int i = 0; i < NumPeaks; i++)
            //        {
            //            var a = sr.ReadLine().Split();
            //            NumVertex = int.Parse(a[0]);
            //            Side = int.Parse(a[1]);
            //            Bound1.Add(new S1(NumVertex, Side));

            //        }
            //    }
            //}
            //catch (IOException e)
            //{
            //    Console.WriteLine("S1 input exeption");
            //    Console.WriteLine(e.Message);
            //}
        }

        private void InputBound2()
        {
            Bound2 = new List<S2>();

            int num = NumR * NumZ;

            for (int i = NumR; i < num; i += NumR)
            {
                Bound2.Add(new S2(i - NumR, i, 0));
            }

            //string path = "Bound2.txt";
            //try
            //{
            //    using (StreamReader sr = new StreamReader(path))
            //    {
            //        int NumBound = int.Parse(sr.ReadLine());
            //        int NumVertex1, NumVertex2, Side;
            //        Bound2 = new List<S2>();
            //        for (int i = 0; i < NumBound; i++)
            //        {
            //            var a = sr.ReadLine().Split();
            //            NumVertex1 = int.Parse(a[0]);
            //            NumVertex2 = int.Parse(a[1]);
            //            Side = int.Parse(a[2]);
            //            Bound2.Add(new S2(NumVertex1, NumVertex2, Side));
            //        }
            //    }
            //}
            //catch (IOException e)
            //{
            //    Console.WriteLine("S1 input exeption");
            //    Console.WriteLine(e.Message);
            //}


        }

        private void InputPrams()
        {
            string path = "Params.txt";
            try
            {
                using (StreamReader sr = new StreamReader(path, Encoding.Default))
                {
                    int num = int.Parse(sr.ReadLine());
                    Params = new List<Parameter>();
                    for (int i = 0; i < num; i++)
                    {
                        var tmp = sr.ReadLine().Split();
                        double sigma = double.Parse(tmp[0]);
                        double Mu = double.Parse(tmp[1]);
                        Params.Add(new Parameter(sigma, Mu));
                    }
                }
            }
            catch (IOException e)
            {
                Console.WriteLine("Param input exeption");
                Console.WriteLine(e.Message);
            }
        }

        private void InputNumofPram()
        {
            NumOfParameters = new int[(NumZ - 1) * (NumR - 1)];
            //string path = "Param_on_Element.txt";
            //try
            //{
            //    using (StreamReader sr = new StreamReader(path, Encoding.Default))
            //    {
            //        for (int i = 0; i < NumOfParameters.Length; i++)
            //        {
            //            NumOfParameters[i] = int.Parse(sr.ReadLine());
            //        }
            //    }
            //}
            //catch (IOException e)
            //{
            //    Console.WriteLine("Param input exeption");
            //    Console.WriteLine(e.Message);
            //}
        }

        private void InputGrid()
        {
            string path = "Grid.txt";
            try
            {
                using (StreamReader sr = new StreamReader(path))
                {
                    NumR = int.Parse(sr.ReadLine());
                    NumZ = int.Parse(sr.ReadLine());
                    NumElem = int.Parse(sr.ReadLine());
                    DiscR = double.Parse(sr.ReadLine());
                    DiscZ = double.Parse(sr.ReadLine());
                    R0 = double.Parse(sr.ReadLine());
                    R1 = double.Parse(sr.ReadLine());
                    Z0 = double.Parse(sr.ReadLine());
                    Z1 = double.Parse(sr.ReadLine());
                    Hr = double.Parse(sr.ReadLine());
                    Hz = double.Parse(sr.ReadLine());

                }
            }
            catch (IOException e)
            {
                Console.WriteLine("Grid input exeption");
                Console.WriteLine(e.Message);
            }
        }

        private void InputGridPointSource()
        {
            string path = "GridPoint.txt";
            try
            {
                using (StreamReader sr = new StreamReader(path))
                {
                    var tmp = sr.ReadLine().Split();

                    //NumPointSource = int.Parse(tmp[2]);
                    PointSource = new Vertex(double.Parse(tmp[0]), double.Parse(tmp[1]), 1);

                    tmp = sr.ReadLine().Split();
                    NumRDown = int.Parse(tmp[0]);
                    NumZDown = int.Parse(tmp[1]);
                    NumRUp = int.Parse(tmp[2]);
                    NumZUp = int.Parse(tmp[3]);
                    NumR = NumRDown + NumRUp - 1;
                    NumZ = NumZDown + NumZUp - 1;

                    tmp = sr.ReadLine().Split();
                    HrDown = double.Parse(tmp[0]);
                    HrUp = double.Parse(tmp[1]);
                    HzDown = double.Parse(tmp[2]);
                    HzUp = double.Parse(tmp[3]);

                    tmp = sr.ReadLine().Split();
                    DiscRDown = double.Parse(tmp[0]);
                    DiscRUp = double.Parse(tmp[1]);
                    DiscZDown = double.Parse(tmp[2]);
                    DiscZUp = double.Parse(tmp[3]);

                    tmp = sr.ReadLine().Split();
                    R0 = double.Parse(tmp[0]);
                    R1 = double.Parse(tmp[1]);

                    tmp = sr.ReadLine().Split();
                    Z0 = double.Parse(tmp[0]);
                    Z1 = double.Parse(tmp[1]);
                }
            }
            catch (IOException e)
            {
                Console.WriteLine("Grid input exeption");
                Console.WriteLine(e.Message);
            }
        }

        private void InputTime()
        {
            string path = "Time.txt";
            try
            {
                using (StreamReader sr = new StreamReader(path, Encoding.Default))
                {
                    TimeH = double.Parse(sr.ReadLine());
                    TimeStart = double.Parse(sr.ReadLine());
                    TimeEnd = double.Parse(sr.ReadLine());
                }
            }
            catch (IOException e)
            {
                Console.WriteLine("Time input exeption");
                Console.WriteLine(e.Message);
            }
        }

        private void GetElements()
        {

            Elements = new List<Element>();
            Element elementTmp;


            for (int i = 0; i < NumZ - 1; i++)
            {
                for (int j = i * NumR; j < i * NumR + NumR - 1; j++)
                {
                    int[] arrTmp = new int[4];
                    arrTmp[0] = j;
                    arrTmp[1] = j + 1;
                    arrTmp[2] = j + NumR;
                    arrTmp[3] = j + NumR + 1;
                    elementTmp = new Element(arrTmp);
                    Elements.Add(elementTmp);
                }
            }
            for (int i = 0; i < NumOfParameters.Count(); i++)
            {
                Elements[i].NMat = NumOfParameters[i];
            }

        }
    }
}
