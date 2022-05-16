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
        private double Hr;
        private double Hz;
        private int[] NumOfParameters;

        public List<S1> Bound1 { get; set; }
        public List<Parameter> Params { get; set; }
        public List<double> Time { get; set; }
        public List<Vertex> RZ { get; set; }
        public List<Element> Elements { get; set; }

        public IO()
        {
            InputGrid();
            InputTime();
            GetGrid();
            GetTime();
            InputNumofPram();
            GetElements();
            InputBound1();
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

        private void InputBound1()
        {
            string path = "Bound1.txt";
            try
            {
                using (StreamReader sr = new StreamReader(path))
                {
                    int NumPeaks = int.Parse(sr.ReadLine());
                    int NumVertex, Side;
                    Bound1 = new List<S1>();
                    for (int i = 0; i < NumPeaks; i++)
                    {
                        var a = sr.ReadLine().Split();
                        NumVertex = int.Parse(a[0]);
                        Side = int.Parse(a[1]);
                        Bound1.Add(new S1(NumVertex, Side));
                    }
                }
            }
            catch (IOException e)
            {
                Console.WriteLine("S1 input exeption");
                Console.WriteLine(e.Message);
            }
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
                        double Mu = double.Parse(tmp[0]);
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
            string path = "Param_on_Element.txt";
            try
            {
                using (StreamReader sr = new StreamReader(path, Encoding.Default))
                {
                    for (int i = 0; i < NumOfParameters.Length; i++)
                    {
                        NumOfParameters[i] = int.Parse(sr.ReadLine());
                    }
                }
            }
            catch (IOException e)
            {
                Console.WriteLine("Param input exeption");
                Console.WriteLine(e.Message);
            }
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
