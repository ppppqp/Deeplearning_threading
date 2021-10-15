using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using BioNet;

namespace CommonFunction
{
    public class BondParameter
    {
        public String ResidueType;
        public String BondAtom1;
        public String BondAtom2;
        public Double BondLengthAvg;
        public Double BondLengthStd;
        public Double ForceConstant;
        public BondParameter() { }

        public BondParameter(String residuetype, String bondatom1, String bondatom2, Double lengthavg, Double lengthstd, Double forceconstant)
        {
            this.ResidueType = residuetype;
            this.BondAtom1 = bondatom1;
            this.BondAtom2 = bondatom2;
            this.BondLengthAvg = lengthavg;
            this.BondLengthStd = lengthstd;
            this.ForceConstant = forceconstant;
        }
    }

    public class RamaParameter
    {
        public String ResidueType;
        public Char SecondStructure;
        public Double[,] H;
        public List<Double> PhiGrid;
        public List<Double> PsiGrid;

        public RamaParameter()
        {
            this.H = new Double[2, 2];
            this.PhiGrid = new List<Double>();
            this.PsiGrid = new List<Double>();
        }

        public RamaParameter(String residuetype, Char ss,Double[,] H,List<Double> PhiGrid, List<Double> PsiGrid)
        {
            this.ResidueType = residuetype;
            this.SecondStructure = ss;
            this.H = new Double[2, 2];
            this.H = H;
            this.PhiGrid = new List<Double>();
            this.PhiGrid = PhiGrid;
            this.PsiGrid = new List<Double>();
            this.PsiGrid = PsiGrid;
        }
    }

    public class CADistanceRestraintParameter
    {
        public int headResidueNum;
        public int tailResidueNum;
        public Double di;
        public Double dt;
        public Double TM;
        public Double sim;
        public Double lambda1;
        public Double mu1;
        public Double sigma1;
        public Double lambda2;
        public Double mu2;
        public Double sigma2;
        public Double lambda3;
        public Double mu3;
        public Double sigma3;
    }
    class Common
    {
        public const Double impossibleMaxDouble = 99999999;
        public const Double impossibleMinDouble = -99999999;
        public const Double BoltzmannConstant= 1.38e-23;
        public static String[] AAtype = { "ARG", "GLN", "PHE", "TYR", "TRP", "LYS", "GLY", "ALA", "HIS", "SER", "PRO", "GLU", "ASP", "THR", "CYS", "MET", "LEU", "ASN", "ILE", "VAL", "ASX", "GLX", "UNK" };
        public static Char[] AAtypeChar = { 'R', 'Q', 'F', 'Y', 'W', 'K', 'G', 'A', 'H', 'S', 'P', 'E', 'D', 'T', 'C', 'M', 'L', 'N', 'I', 'V', 'B', 'Z', 'X' };
        public static Char[] SStype = { 'H', 'C', 'E' };
        public const String BondLengthParametersDataPath = @"../../data/Bond/bondlength_pdb.dat";
        public static BondParameter[] BondLengthParameters;
        public static String RamaParametersDataPath = @"../../data/Rama/rama.dat";
        public static RamaParameter[] RamaParameters;
        public static CADistanceRestraintParameter[] CADistanceRestraintParameters;
        public static String CADistanceRestraintDirPath = @"../../data/CARestraint";
        public Common()
        {
            InitialBondLengthParameters();
            InitialRamaParameters();
        }

        public Common(StreamReader CADisResfile)
        {
            InitialBondLengthParameters();
            InitialRamaParameters();
            InitialCADistanceRestraintParameters(CADisResfile);
        }

        public void InitialBondLengthParameters()
        {
            StreamReader bondlengthsr = new StreamReader(BondLengthParametersDataPath);
            String title = bondlengthsr.ReadLine();
            List<String> templines = new List<String>();
            int length = 0;
            while (bondlengthsr.Peek() > -1)
            {
                String line = bondlengthsr.ReadLine().Trim();
                if (line.Length > 0)
                {
                    length++;
                    templines.Add(line);
                }
            }
            bondlengthsr.Close();
            BondLengthParameters = new BondParameter[length];
            for(int i=0;i<length;i++)
            {
                String[] temp = templines.ElementAt(i).Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                BondLengthParameters[i] = new BondParameter(temp[0], temp[1], temp[2], Convert.ToDouble(temp[3]), Convert.ToDouble(temp[4]), Convert.ToDouble(temp[5]));
            }
        }

        public void InitialRamaParameters()
        {
            StreamReader ramasr = new StreamReader(RamaParametersDataPath);
            List<RamaParameter> tmpRamaParamList = new List<RamaParameter>();

            while (ramasr.Peek()>-1)
            {
                String line = ramasr.ReadLine().Trim();
                if (line.Length > 0)
                {
                    if(line.StartsWith(">"))
                    {
                        RamaParameter rpNode = new RamaParameter();
                        rpNode.ResidueType = line.Substring(1, 3);
                        rpNode.SecondStructure = line.ElementAt(5);
                        tmpRamaParamList.Add(rpNode);
                    }
                    else if (line.StartsWith("#Hmatrix"))
                    {
                        String[] H1 = ramasr.ReadLine().Split(new Char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                        String[] H2 = ramasr.ReadLine().Split(new Char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                        tmpRamaParamList.Last().H[0, 0] = Convert.ToDouble(H1[0]);
                        tmpRamaParamList.Last().H[0, 1] = Convert.ToDouble(H1[1]);
                        tmpRamaParamList.Last().H[1, 0] = Convert.ToDouble(H2[0]);
                        tmpRamaParamList.Last().H[1, 1] = Convert.ToDouble(H2[1]);
                    }
                    else if(line.StartsWith("#PHI PSI"))
                    {
                        continue;
                    }
                    else
                    {
                        String[] tmpstr = line.Split(new Char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                        tmpRamaParamList.Last().PhiGrid.Add(Convert.ToDouble(tmpstr[0]));
                        tmpRamaParamList.Last().PsiGrid.Add(Convert.ToDouble(tmpstr[1]));
                    }

                }
            }
            ramasr.Close();
            RamaParameters = tmpRamaParamList.ToArray();
        }

        public void InitialCADistanceRestraintParameters(StreamReader CADisResfile)
        {
        }
        static public Double Variance(List<Double> datalist)
        {
            Double variance = 0;
            Double mean = datalist.Average();
            foreach (var data in datalist)
            {
                variance += Math.Pow((data - mean), 2);
            }
            variance /= datalist.Count;
            return variance;
        }

        static public Double Standard(List<Double> datalist)
        {
            return Math.Sqrt(Variance(datalist));
        }

        static public Double GaussianKernel(Double[] x, Double[,] H)
        {
            //2-d
            //x=(x1,x2) H=|a11 a12|
            //            |a21 a22|
            Double k = 0;
            //xT*H-1*x=(a11*x2^2+a22*x1^2-(a12+a21)x1*x2)/(a11*a22-a12*a21)
            Double powsum = (H[0, 0] * x[1] * x[1] + H[1, 1] * x[0] * x[0] - (H[0, 1] + H[1, 0]) * x[0] * x[1]) / (H[0, 0] * H[1, 1] - H[1, 0] * H[0, 1]);

            //Console.WriteLine(powsum);
            k = Math.Exp(-0.5 * powsum) / (2.0 * Math.PI * Math.Sqrt(H[0, 0] * H[1, 1] - H[1, 0] * H[0, 1]));
            return k;
        }

        static public Double KDE2calculate(Double[,] H, List<Double> PhiGrid,List<Double> PsiGrid, Double phi, Double psi)
        {
            Double p = 0;
            for (int i = 0; i < PhiGrid.Count; i++)
            {
                Double[] x = new Double[2];
                x[0] = phi - PhiGrid.ElementAt(i);
                x[1] = psi - PsiGrid.ElementAt(i);
                p += GaussianKernel(x, H);
            }
            int n = PhiGrid.Count;
            //Console.WriteLine(n);
            p /= n;
            return p;
        }

        static public Double GaussianDensity(Double d,Double mu,Double sigma)
        {
            Double p = 0;
            Double exp = Math.Exp(-0.5 * Math.Pow((d - mu) / sigma, 2));
            p = exp / (sigma * Math.Sqrt(2 * Math.PI));
            return p;
        }

        static public Double MixGaussianDensity(Double d, Double lambda1, Double mu1, Double sigma1, Double lambda2, Double mu2, Double sigma2, Double lambda3, Double mu3, Double sigma3)
        {
            Double p = 0;
            Double p1 = GaussianDensity(d, mu1, sigma1);
            Double p2 = GaussianDensity(d, mu2, sigma2);
            Double p3 = GaussianDensity(d, mu3, sigma3);
            p = lambda1 * p1 + lambda2 * p2 + lambda3 * p3;
            return p;
        }
       
        static public bool Kabsch(List<Point3D> xset,List<Point3D> yset,int mode,ref double rms,ref double[] t, ref double[,] u)
        {
            int n = xset.Count;
            double[,] x = new double[n, 3];
            double[,] y = new double[n, 3];
            for(int row=0;row<n;row++)
            {
                x[row, 0] = xset.ElementAt(row).X;
                x[row, 1] = xset.ElementAt(row).Y;
                x[row, 2] = xset.ElementAt(row).Z;
                y[row, 0] = yset.ElementAt(row).X;
                y[row, 1] = yset.ElementAt(row).Y;
                y[row, 2] = yset.ElementAt(row).Z;
            }
            int i, j, m, m1, l, k;
            double e0, rms1, d, h, g;
            double cth, sth, sqrth, p, det, sigma;
            double[] xc = new double[3];
            double[] yc = new double[3];
            double[,] a = new double[3, 3];
            double[,] b = new double[3, 3];
            double[,] r = new double[3, 3];
            double[] e = new double[3];
            double[] rr = new double[6];
            double[] ss = new double[6];
	        double sqrt3 = 1.73205080756888, tol = 0.01;
            int[] ip = { 0, 1, 3, 1, 2, 4, 3, 4, 5 };
            int[] ip2312 = { 1, 2, 0, 1 };
            int a_failed = 0, b_failed = 0;
            double epsilon = 0.00000001;
            //initializtation
            rms = 0;
            rms1 = 0;
            e0 = 0;
            double[] c1 = new double[3];
            double[] c2 = new double[3];
            double[] s1 = new double[3];
            double[] s2 = new double[3];
            double[] sx = new double[3];
            double[] sy = new double[3];
            double[] sz = new double[3];
            for (i = 0; i < 3; i++)
            {
                s1[i] = 0.0;
                s2[i] = 0.0;

                sx[i] = 0.0;
                sy[i] = 0.0;
                sz[i] = 0.0;
            }

            for (i = 0; i < 3; i++)
            {
                xc[i] = 0.0;
                yc[i] = 0.0;
                t[i] = 0.0;
                for (j = 0; j < 3; j++)
                {
                    u[i,j] = 0.0;
                    r[i,j] = 0.0;
                    a[i,j] = 0.0;
                    if (i == j)
                    {
                        u[i,j] = 1.0;
                        a[i,j] = 1.0;
                    }
                }
            }

            if (n < 1)
            {
                return false;
            }
            //compute centers for vector sets x, y
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < 3; j++)
                {
                    c1[j] = x[i,j];
                    c2[j] = y[i,j];

                    s1[j] += c1[j];
                    s2[j] += c2[j];
                }

                for (j = 0; j < 3; j++)
                {
                    sx[j] += c1[0] * c2[j];
                    sy[j] += c1[1] * c2[j];
                    sz[j] += c1[2] * c2[j];
                }
            }
            for (i = 0; i < 3; i++)
            {
                xc[i] = s1[i] / n;
                yc[i] = s2[i] / n;
            }
            if (mode == 2 || mode == 0)
            {
                for (int mm = 0; mm < n; mm++)
                {
                    for (int nn = 0; nn < 3; nn++)
                    {
                        e0 += (x[mm,nn] - xc[nn]) * (x[mm,nn] - xc[nn]) + (y[mm,nn] - yc[nn]) * (y[mm,nn] - yc[nn]);
                    }
                }
            }
            for (j = 0; j < 3; j++)
            {
                r[j,0] = sx[j] - s1[0] * s2[j] / n;
                r[j,1] = sy[j] - s1[1] * s2[j] / n;
                r[j,2] = sz[j] - s1[2] * s2[j] / n;
            }
            //compute determinat of matrix r
            det = r[0,0] * (r[1,1] * r[2,2] - r[1,2] * r[2,1])
		        - r[0,1] * (r[1,0] * r[2,2] - r[1,2] * r[2,0])
		        + r[0,2] * (r[1,0] * r[2,1] - r[1,1] * r[2,0]);
            sigma = det;
            //compute tras(r)*r
            m = 0;
            for (j = 0; j < 3; j++)
            {
                for (i = 0; i <= j; i++)
                {
                    rr[m] = r[0,i] * r[0,j] + r[1,i] * r[1,j] + r[2,i] * r[2,j];
                    m++;
                }
            }
            double spur = (rr[0] + rr[2] + rr[5]) / 3.0;
            double cof = (((((rr[2] * rr[5] - rr[4] * rr[4]) + rr[0] * rr[5])
                       - rr[3] * rr[3]) +rr[0] * rr[2]) -rr[1] * rr[1]) / 3.0;
            det = det * det;

            if (spur > 0)
            {
                d = spur * spur;
                h = d - cof;
                g = (spur * cof - det) / 2.0 - spur * h;

                if (h > 0)
                {
                    sqrth = Math.Sqrt(h);
                    d = h * h * h - g * g;
                    if (d < 0.0) d = 0.0;
                    d = Math.Atan2(Math.Sqrt(d), -g) / 3.0;
                    cth = sqrth * Math.Cos(d);
                    sth = sqrth * sqrt3 * Math.Sin(d);
                    e[0] = (spur + cth) + cth;
                    e[1] = (spur - cth) + sth;
                    e[2] = (spur - cth) - sth;

                    if (mode != 0)
                    {//compute a                
                        for (l = 0; l < 3; l = l + 2)
                        {
                            d = e[l];
                            ss[0] = (d - rr[2]) * (d - rr[5]) - rr[4] * rr[4];
                            ss[1] = (d - rr[5]) * rr[1] + rr[3] * rr[4];
                            ss[2] = (d - rr[0]) * (d - rr[5]) - rr[3] * rr[3];
                            ss[3] = (d - rr[2]) * rr[3] + rr[1] * rr[4];
                            ss[4] = (d - rr[0]) * rr[4] + rr[1] * rr[3];
                            ss[5] = (d - rr[0]) * (d - rr[2]) - rr[1] * rr[1];

                            if (Math.Abs(ss[0]) <= epsilon) ss[0] = 0.0;
                            if (Math.Abs(ss[1]) <= epsilon) ss[1] = 0.0;
                            if (Math.Abs(ss[2]) <= epsilon) ss[2] = 0.0;
                            if (Math.Abs(ss[3]) <= epsilon) ss[3] = 0.0;
                            if (Math.Abs(ss[4]) <= epsilon) ss[4] = 0.0;
                            if (Math.Abs(ss[5]) <= epsilon) ss[5] = 0.0;

                            if (Math.Abs(ss[0]) >= Math.Abs(ss[2]))
                            {
                                j = 0;
                                if (Math.Abs(ss[0]) < Math.Abs(ss[5]))
                                {
                                    j = 2;
                                }
                            }
                            else if (Math.Abs(ss[2]) >= Math.Abs(ss[5]))
                            {
                                j = 1;
                            }
                            else
                            {
                                j = 2;
                            }

                            d = 0.0;
                            j = 3 * j;
                            for (i = 0; i < 3; i++)
                            {
                                k = ip[i + j];
                                a[i,l] = ss[k];
                                d = d + ss[k] * ss[k];
                            }


                            //if( d > 0.0 ) d = 1.0 / sqrt(d);
                            if (d > epsilon) d = 1.0 / Math.Sqrt(d);
                            else d = 0.0;
                            for (i = 0; i < 3; i++)
                            {
                                a[i,l] = a[i,l] * d;
                            }
                        }//for l

                        d = a[0,0] * a[0,2] + a[1,0] * a[1,2] + a[2,0] * a[2,2];
                        if ((e[0] - e[1]) > (e[1] - e[2]))
                        {
                            m1 = 2;
                            m = 0;
                        }
                        else
                        {
                            m1 = 0;
                            m = 2;
                        }
                        p = 0;
                        for (i = 0; i < 3; i++)
                        {
                            a[i,m1] = a[i,m1] - d * a[i,m];
                            p = p + a[i,m1] * a[i,m1];
                        }
                        if (p <= tol)
                        {
                            p = 1.0;
                            for (i = 0; i < 3; i++)
                            {
                                if (p < Math.Abs(a[i,m]))
                                {
                                    continue;
                                }
                                p = Math.Abs(a[i,m]);
                                j = i;
                            }
                            k = ip2312[j];
                            l = ip2312[j + 1];
                            p = Math.Sqrt(a[k,m] * a[k,m] + a[l,m] * a[l,m]);
                            if (p > tol)
                            {
                                a[j,m1] = 0.0;
                                a[k,m1] = -a[l,m] / p;
                                a[l,m1] = a[k,m] / p;
                            }
                            else
                            {//goto 40
                                a_failed = 1;
                            }
                        }//if p<=tol
                        else
                        {
                            p = 1.0 / Math.Sqrt(p);
                            for (i = 0; i < 3; i++)
                            {
                                a[i,m1] = a[i,m1] * p;
                            }
                        }//else p<=tol  
                        if (a_failed != 1)
                        {
                            a[0,1] = a[1,2] * a[2,0] - a[1,0] * a[2,2];
                            a[1,1] = a[2,2] * a[0,0] - a[2,0] * a[0,2];
                            a[2,1] = a[0,2] * a[1,0] - a[0,0] * a[1,2];
                        }
                    }//if(mode!=0)       
                }//h>0

                //compute b anyway
                if (mode != 0 && a_failed != 1)//a is computed correctly
                {
                    //compute b
                    for (l = 0; l < 2; l++)
                    {
                        d = 0.0;
                        for (i = 0; i < 3; i++)
                        {
                            b[i,l] = r[i,0] * a[0,l] + r[i,1] * a[1,l] + r[i,2] * a[2,l];
                            d = d + b[i,l] * b[i,l];
                        }
                        //if( d > 0 ) d = 1.0 / sqrt(d);
                        if (d > epsilon) d = 1.0 / Math.Sqrt(d);
                        else d = 0.0;
                        for (i = 0; i < 3; i++)
                        {
                            b[i,l] = b[i,l] * d;
                        }
                    }
                    d = b[0,0] * b[0,1] + b[1,0] * b[1,1] + b[2,0] * b[2,1];
                    p = 0.0;

                    for (i = 0; i < 3; i++)
                    {
                        b[i,1] = b[i,1] - d * b[i,0];
                        p += b[i,1] * b[i,1];
                    }

                    if (p <= tol)
                    {
                        p = 1.0;
                        for (i = 0; i < 3; i++)
                        {
                            if (p < Math.Abs(b[i,0]))
                            {
                                continue;
                            }
                            p = Math.Abs(b[i,0]);
                            j = i;
                        }
                        k = ip2312[j];
                        l = ip2312[j + 1];
                        p = Math.Sqrt(b[k,0] * b[k,0] + b[l,0] * b[l,0]);
                        if (p > tol)
                        {
                            b[j,1] = 0.0;
                            b[k,1] = -b[l,0] / p;
                            b[l,1] = b[k,0] / p;
                        }
                        else
                        {
                            //goto 40
                            b_failed = 1;
                        }
                    }//if( p <= tol )
                    else
                    {
                        p = 1.0 / Math.Sqrt(p);
                        for (i = 0; i < 3; i++)
                        {
                            b[i,1] = b[i,1] * p;
                        }
                    }
                    if (b_failed != 1)
                    {
                        b[0,2] = b[1,0] * b[2,1] - b[1,1] * b[2,0];
                        b[1,2] = b[2,0] * b[0,1] - b[2,1] * b[0,0];
                        b[2,2] = b[0,0] * b[1,1] - b[0,1] * b[1,0];
                        //compute u
                        for (i = 0; i < 3; i++)
                        {
                            for (j = 0; j < 3; j++)
                            {
                                u[i,j] = b[i,0] * a[j,0] + b[i,1] * a[j,1]+b[i,2] * a[j,2];
                            }
                        }
                    }

                    //compute t
                    for (i = 0; i < 3; i++)
                    {
                        t[i] = ((yc[i] - u[i,0] * xc[0]) - u[i,1] * xc[1])-u[i,2] * xc[2];
                    }
                }//if(mode!=0 && a_failed!=1)
            }//spur>0
            else //just compute t and errors
            {
                //compute t
                for (i = 0; i < 3; i++)
                {
                    t[i] = ((yc[i] - u[i,0] * xc[0]) - u[i,1] * xc[1]) - u[i,2] * xc[2];
                }
            }//else spur>0 

            //compute rms
            for (i = 0; i < 3; i++)
            {
                if (e[i] < 0) e[i] = 0;
                e[i] = Math.Sqrt(e[i]);
            }
            d = e[2];
            if (sigma < 0.0)
            {
                d = -d;
            }
            d = (d + e[1]) + e[0];

            if (mode == 2 || mode == 0)
            {
                rms1 = (e0 - d) - d;
                if (rms1 < 0.0)
                    rms1 = 0.0;
            }

            rms = rms1;

            return true;
        }
    }
}
