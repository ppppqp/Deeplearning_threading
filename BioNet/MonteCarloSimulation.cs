using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.IO;
using BioNet;
using CommonFunction;

namespace MonteCarloSimulation
{
    public class MonteCarloResult
    {
        public Double T;
        public Double startE;
        public List<Double> decoyE;
        public List<Protein> decoys;
        public Double acceptRadio;
        public MonteCarloResult()
        {
            this.decoyE= new List<Double>(); ;
            this.decoys= new List<Protein>();
        }
    }

    public class MonteCarloSimulation
    {
        public MonteCarloResult MCR;
        //public const Double kb = 1.38e-23;
        
        public MonteCarloSimulation() { }
        public MonteCarloSimulation(Double Ti,int MaxStep,Protein initprotein, StreamWriter swinfo)
        {
            
            this.MCR = new MonteCarloResult();
            this.MCR.T = Ti;

            Protein pro = initprotein;
            Energy energy = new Energy();
            Double eng = energy.TestEnergy(initprotein);
            this.MCR.startE = eng;
            //int MaxStep = 3000;//100000;
            Double kb = Common.BoltzmannConstant;
            Double rate = 1.0;//0.9
            int decoy = 0;
            swinfo.WriteLine("step       decoy      energy     temperature accept_radio");
            for (int step = 1; step <= MaxStep; step++)
            {
                Double T;
                int w = Convert.ToInt32(Math.Floor(step / 10000.0));
                T = Ti * Math.Pow(rate, w);//模拟退火
                this.MCR.T = T;

                Move mvp = new Move();
                Protein newpro =mvp.Move3(pro);
                Energy newenergy = new Energy();
                Double neweng = newenergy.TestEnergy(newpro);
                if (neweng <= eng)
                {
                    pro = newpro;
                    eng = neweng;
                    this.MCR.decoyE.Add(eng);
                    this.MCR.decoys.Add(pro);
                    decoy++;
                    swinfo.WriteLine("{0,-10:f0} {1,-10:f0} {2,-10:f5} {3,-11:e3} {4,-12:f4}", step, decoy, eng, T, Convert.ToDouble(decoy) / step);
                }
                else
                {
                    Double E = Math.Exp(-(neweng - eng) / (kb * T));
                    RandomNumber rand = new RandomNumber();
                    Double r = rand.GetUnifromRandomNumber(10);
                    if (r < E)
                    {
                        this.MCR.decoyE.Add(neweng);
                        this.MCR.decoys.Add(newpro);
                        decoy++;
                        swinfo.WriteLine("{0,-10:f0} {1,-10:f0} {2,-10:f5} {3,-11:e3} {4,-12:f4}", step, decoy, eng, T, Convert.ToDouble(decoy) / step);
                        pro = newpro;
                        eng = neweng;
                    }
                }

            }
            swinfo.Close();
            this.MCR.acceptRadio = Convert.ToDouble(this.MCR.decoyE.Count) / MaxStep;
        }

    }
    public class ReplicaExchangeMonteCarloSimulation
    {
        //public const Double kb = 1.38e-23;
        public ReplicaExchangeMonteCarloSimulation() { }

        public ReplicaExchangeMonteCarloSimulation(Protein p, Double[] T,int loopTimes, int cpuNum, int MCMaxStep)
        {
            this.Run(p, T, loopTimes, cpuNum, MCMaxStep);
        }

        public void Run(Protein p, Double[] T, int loopTimes, int cpuNum,int MCMaxStep)
        {

            StreamWriter swlog = new StreamWriter(@"REMC.log");
            swlog.WriteLine("initial parameters");
            swlog.WriteLine("Max loop times:{0}", loopTimes);
            swlog.WriteLine("Max cpu numbers:{0}", cpuNum);
            swlog.WriteLine("Single Monte Carlo Simulation loop times:{0}", MCMaxStep);
            swlog.WriteLine("initial temperatures:");
            for(int t=0;t<T.Length;t++)
            {
                swlog.WriteLine("T" + (t + 1) + "={0,-15:e3}", T[t]);
            }
            swlog.WriteLine();
            swlog.WriteLine();

            //int loopTimes = 10;
            //int cpuNum = 4;
            //int MCMaxStep=3000;
            int taskNum = T.Length;
            List<Double> SwapRadios = new List<Double>();
            Protein[] pro = new Protein[taskNum];
            for (int i = 0; i < taskNum; i++)
            {
                pro[i] = p;
            }
            for (int i = 0; i < loopTimes; i++)
            {
                swlog.WriteLine("########################################## Replica Exchange Monte Carlo Simulation Round {0} ##########################################", i + 1);
                swlog.WriteLine("Replica Temperature Start_Energy Last_Energy Decoys_Numbers Accept_Radio");

                MonteCarloResult[] mcrs = new MonteCarloResult[taskNum];
                Double[] lasteng = new Double[taskNum];
                StreamWriter[] swinfo = new StreamWriter[taskNum];
                for (int k = 0; k < taskNum; k++)
                {
                    swinfo[k] = new StreamWriter(@"round" + (i + 1) + "_task" + (k + 1) + ".enginfo.tab");
                }

                //Run Monte Carlo Simulations
                Parallel.For(0, taskNum, new ParallelOptions() { MaxDegreeOfParallelism = cpuNum }, index => {

                    Console.WriteLine("round:" + (i + 1) + " task:" + (index + 1));
                    MonteCarloSimulation MC = new MonteCarloSimulation(T[index],MCMaxStep ,pro[index], swinfo[index]);
                    mcrs[index] = MC.MCR;

                });
                for (int k = 0; k < taskNum; k++)
                {
                    swinfo[k].Close();
                }
                for(int t=0;t<taskNum;t++)
                {
                    Double LowestE = mcrs[t].decoyE.Min();
                    int LowestEIndex = mcrs[t].decoyE.FindIndex(e => e == LowestE);
                    Protein LowestEdecoy = mcrs[t].decoys.ElementAt(LowestEIndex);
                    StreamWriter decoysw = new StreamWriter(@"round" + (i + 1) + "_task" + (t + 1) + ".lowestE.pdb");
                    //LowestEIndex+1 是为了显示下标从1开始，实际计算是从0开始
                    decoysw.WriteLine("MODEL {0}, Energy={1,-10:f4}", LowestEIndex+1, LowestE);
                    LowestEdecoy.WritePDB(decoysw);
                    decoysw.Close();
                }
                //Swap Temperatures
                for (int j = 0; j < taskNum; j++)
                {
                    T[j] = mcrs[j].T;
                    pro[j] = mcrs[j].decoys.Last();
                    lasteng[j] = mcrs[j].decoyE.Last();
                    //j+1 是为了显示下标从1开始，实际计算是从0开始
                    swlog.WriteLine("{0,7:f0} {1,11:e3} {2,12:f5} {3,11:f5} {4,14:f0} {5,12:f3}", j+1, T[j], mcrs[j].startE, mcrs[j].decoyE.Last(), mcrs[j].decoys.Count, mcrs[j].acceptRadio);
                }
                swlog.WriteLine();
                swlog.WriteLine();
                Console.WriteLine("swap T at times:" + (i + 1));
                //Double[] newT = SwapT(lasteng, T);
                //Double[] newT = ISwapT(lasteng, T, 4);
                Double[] newT = NSwapT(lasteng, T, i);
                //Swap Infomation
                if (i != (loopTimes - 1))
                {
                    swlog.WriteLine("Swap Replica Temperature");
                    int swaptimes = 0;
                    for (int t = 0; t < taskNum; t++)
                    {
                        int tt = T.ToList().FindIndex(value => value == newT[t]);
                        if (tt == t)
                        {    
                            //t+1 tt+1 是为了显示下标从1开始，实际计算是从0开始
                            swlog.WriteLine("keep {0,5:f0}<-->{1,-5:f0} {2,10:e3}<-->{3,-10:e3} {4,10:f3}<-->{5,-10:f3}", t+1, tt+1, mcrs[t].T, mcrs[tt].T, mcrs[t].decoyE.Last(), mcrs[tt].decoyE.Last());
                        }
                        else
                        {
                            swaptimes++;
                            //t+1 tt+1 是为了显示下标从1开始，实际计算是从0开始
                            swlog.WriteLine("swap {0,5:f0}<-->{1,-5:f0} {2,10:e3}<-->{3,-10:e3} {4,10:f3}<-->{5,-10:f3}", t+1, tt+1, mcrs[t].T, mcrs[tt].T, mcrs[t].decoyE.Last(), mcrs[tt].decoyE.Last());
                        }
                    }

                    SwapRadios.Add(Convert.ToDouble(swaptimes) / taskNum);
                    swlog.WriteLine();
                    swlog.WriteLine();
                }
                T = newT;
            }
            swlog.WriteLine();
            swlog.WriteLine("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Total Swap Information @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
            swlog.WriteLine("LoopTimes SwapRadio");
            for(int i=0;i<loopTimes-1;i++)
            {
                swlog.WriteLine("{0,9:f0} {1,9:f4}",i+1,SwapRadios.ElementAt(i));
            }
            swlog.Close();
        }

        public Double[] NSwapT(Double[] energy,Double[] T,int SwapFlag)
        {
            Double kb = Common.BoltzmannConstant;
            Double[] SwapT = new Double[T.Length];
            for(int t=0;t<T.Length;t++)
            {
                SwapT[t] = T[t];
            }
            int StartIndex = SwapFlag % 2;// 0 or 1   N length T, after N steps, complete inverse order. 
            //0 case 1 2 3 4 5 6-->2 1 4 3 6 5
            //1 case 1 2 3 4 5 6-->1 3 2 5 4 6 
            //Console.WriteLine(StartIndex);
            for (int i = StartIndex; i < T.Length - StartIndex; i += 2)
            {
                Double p0 = Math.Exp((energy[i] - energy[i+1]) * (1.0 / (kb * T[i]) - 1.0 / (kb * T[i+1])));
                Double p = Math.Min(1, p0);
                RandomNumber rand = new RandomNumber();
                Double r = rand.GetUnifromRandomNumber(16);
                if(r<p)
                {
                    SwapT[i] = T[i + 1];
                    SwapT[i + 1] = T[i];
                    //i+1 i+2 是为了显示下标从1开始，实际计算是从0开始
                    Console.WriteLine("swap {0} {1}", i+1, i+2);
                }
                else
                {
                    SwapT[i] = T[i];
                    SwapT[i + 1] = T[i + 1];
                    //i+1 i+2 是为了显示下标从1开始，实际计算是从0开始
                    Console.WriteLine("reject {0} {1}", i+1, i+2);
                }
            }
            //foreach (var node in T)
            //{
            //    Console.Write(" "+node);
            //}
            //Console.WriteLine();
            //foreach (var node in SwapT)
            //{
            //    Console.Write(" "+node);
            //}
            //Console.WriteLine();
            return SwapT;

        }

        public Double[] ISwapT(Double[] energy, Double[] T,int SwapTimes)
        {
            Double[] NewSwapT = new Double[T.Length];
            NewSwapT = SwapT(energy, T);
            if (SwapTimes==1)
            {                
                return NewSwapT;
            }
            else
            {
                for (int repeat = 2; repeat <= SwapTimes; repeat++)
                {
                    Double[] OldSwapT = NewSwapT;
                    List<Double> leftT = new List<Double>();
                    List<Double> leftE = new List<Double>();
                    int count = 0;
                    List<int> index = new List<int>();
                    for (int i = 0; i < T.Length; i++)
                    {
                        if (OldSwapT.ElementAt(i) == T.ElementAt(i))
                        {
                            count++;
                            index.Add(i);
                            leftT.Add(T.ElementAt(i));
                            leftE.Add(energy.ElementAt(i));
                        }
                    }
                    if (count <= 1)
                    {
                        return NewSwapT;
                    }
                    else
                    {
                        Double[] tempT = SwapT(leftE.ToArray(), leftT.ToArray());
                        for (int j = 0; j < tempT.Length; j++)
                        {
                            NewSwapT[index.ElementAt(j)] = tempT[j];
                        }
                    }
                }
            }
            return NewSwapT;
        }
        public Double[] SwapT(Double[] energy, Double[] T)
        {
            //T enegy should be 2n length arrary
            Double[] newT = new Double[T.Length];
            List<int> index = new List<int>();
            int memberNum = 2;
            int groupNum = T.Length / memberNum;
            int[,] groups = new int[groupNum, memberNum];
            Double kb = Common.BoltzmannConstant;
            for (int i = 0; i < T.Length; i++)
            {
                index.Add(i);
            }
            for (int j = 0; j < groupNum; j++)
            {
                for (int k = 0; k < memberNum; k++)
                {
                    int value = index[new Random((int)DateTime.Now.Ticks).Next(0, index.Count)];
                    index.Remove(value);
                    groups[j, k] = value;
                    Thread.Sleep(20);
                }
            }
            for (int i = 0; i < groupNum; i++)
            {
                Double p0 = Math.Exp((energy[groups[i, 0]] - energy[groups[i, 1]]) * (1.0 / (kb * T[groups[i, 0]]) - 1.0 / (kb * T[groups[i, 1]])));
                Double p = Math.Min(1, p0);
                RandomNumber rand = new RandomNumber();
                Double r = rand.GetUnifromRandomNumber(16);
                //Console.WriteLine(p);
                if (r < p)
                {
                    //swap
                    newT[groups[i, 0]] = T[groups[i, 1]];
                    newT[groups[i, 1]] = T[groups[i, 0]];
                    //groups[i, 0]+1 groups[i, 1]+1 是为了显示下标从1开始，实际计算是从0开始
                    Console.WriteLine("swap {0} {1}", groups[i, 0] + 1, groups[i, 1] + 1);
                }
                else
                {
                    newT[groups[i, 0]] = T[groups[i, 0]];
                    newT[groups[i, 1]] = T[groups[i, 1]];
                    //groups[i, 0]+1 groups[i, 1]+1 是为了显示下标从1开始，实际计算是从0开始
                    Console.WriteLine("reject {0} {1}", groups[i, 0] + 1, groups[i, 1] + 1);
                }

            }
            return newT;
        }

        static Protein XYZToProtein(List<Point3D> Points)
        {
            Protein initps = new Protein("");
            Chain chain = new Chain();
            chain.chainID = 'A';
            for (int i = 0; i < Points.Count; i++)
            {
                Point3D p = Points.ElementAt(i);
                Atom atom = new Atom();
                atom.atom = 'C';
                atom.atomtype = "CA";
                atom.atomserial = i + 1;
                atom.altLoc = ' ';
                atom.Xlaber = p.X;
                atom.Ylaber = p.Y;
                atom.Zlaber = p.Z;
                atom.occupancy = 1.00;
                atom.tempfactor = 1.00;
                atom.element = atom.atom.ToString();
                atom.charge = "";
                Residue residue = new Residue();
                residue.atoms.Add(atom);
                residue.iCode = ' ';
                residue.residuename = "ALA";
                residue.residueserial = atom.atomserial;
                chain.residues.Add(residue);

            }
            initps.chains.Add(chain);
            return initps;
        }

    }
}
