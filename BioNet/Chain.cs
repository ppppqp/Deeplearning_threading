using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CommonFunction;

namespace BioNet
{
    public class Chain
    {
        //member
        public Char chainID;//key
        public List<Residue> residues;
        //function

        /// <summary>
        /// 初始化蛋白链
        /// </summary>
        public Chain()
        {
            this.residues = new List<Residue>();
        }

        /// <summary>
        /// 初始化蛋白链，以chainID命名
        /// </summary>
        /// <param name="chainID">蛋白质链名</param>
        public Chain(Char chainID)
        {
            this.chainID = chainID;
            this.residues = new List<Residue>();
        }

        /// <summary>
        /// 初始化蛋白链，以chainID命名，设置residues
        /// </summary>
        /// <param name="chainID">蛋白质链名</param>
        /// <param name="residues">氨基酸链</param>
        public Chain(Char chainID,List<Residue> residues)
        {
            this.chainID = chainID;
            this.residues = new List<Residue>();
            foreach(Residue residue in residues.ToList())
            {
                this.residues.Add(residue);
            }
        }

        /// <summary>
        /// 返回链名字
        /// </summary>
        /// <returns>链名</returns>
        public Char GetChainID()
        {
            return this.chainID;
        }

        /// <summary>
        /// 返回氨基酸链
        /// </summary>
        /// <returns>氨基酸链</returns>
        public List<Residue> GetResidues()
        {
            return this.residues;
        }

        /// <summary>
        /// 返回氨基酸序号为residueserial，icode为iCode的氨基酸
        /// </summary>
        /// <param name="residueserial">氨基酸序号</param>
        /// <param name="iCode">氨基酸iCode</param>
        /// <returns>氨基酸</returns>
        public Residue GetResidue(int residueserial,Char iCode=' ')
        {
            Residue residue = this.residues.Find(node => node.residueserial==residueserial && node.iCode == iCode);
            return residue;
        }

        /// <summary>
        /// 返回氨基酸序号为residueserial的第一个氨基酸
        /// </summary>
        /// <param name="residueserial">氨基酸序号</param>
        /// <returns>氨基酸</returns>
        public Residue GetResidue(int residueserial)
        {
            Residue residue = this.residues.Find(node => node.residueserial == residueserial);
            return residue;
        }

        /// <summary>
        /// 返回链内occupancy
        /// </summary>
        /// <returns>List，occupancy</returns>
        public List<Double> GetOccupancy()
        {
            List<Double> occupancy = new List<Double>();
            foreach(Residue residue in this.residues)
            {
                foreach(Atom atom in residue.atoms)
                {
                    occupancy.Add(atom.occupancy);
                }
            }
            return occupancy;
        }

        /// <summary>
        /// 返回链内tempfactor（beta factor）链表
        /// </summary>
        /// <returns>List，beta factor</returns>
        public List<Double> GetTempfactor()
        {
            List<Double> tempfactor = new List<Double>();
            foreach (Residue residue in this.residues)
            {
                foreach (Atom atom in residue.atoms)
                {
                    tempfactor.Add(atom.occupancy);
                }
            }
            return tempfactor;
        }

        /// <summary>
        /// 返回CA原子链，氨基酸没有CA原子的被略过
        /// </summary>
        /// <returns>Chain</returns>
        public Chain GetCAChain()
        {
            Chain newChain = new Chain(this.chainID);
            foreach(Residue residue in this.residues)
            {
                Atom CA = residue.GetAtom("CA");
                if(CA!=null)
                {
                    Residue newResidue = new Residue(residue.residuename, residue.residueserial, residue.iCode);
                    newResidue.atoms.Add(CA);
                    newChain.residues.Add(newResidue);
                }
            }
            return newChain;
        }

        /// <summary>
        /// 返回链内主链原子组成的链，无主链原子氨基酸被略过
        /// </summary>
        /// <returns>Chain</returns>
        public Chain GetMainChain()
        {
            Chain newChain = new Chain(this.chainID);
            foreach(Residue residue in this.residues)
            {
                Residue newResidue = residue.GetMainChainAtoms();
                if(newResidue.atoms.Count!=0)
                {
                    newChain.residues.Add(newResidue);
                }
            }
            return newChain;
        }

        /// <summary>
        /// 自动选择所有residueserial相同的氨基酸中的第一个氨基酸组成的链（select icode）
        /// </summary>
        /// <returns>Chain</returns>
        public Chain AutoSelectiCode()
        {
            List<int> residueserials = new List<int>();
            Chain newChain = new Chain(this.chainID);
            foreach(Residue residue in this.residues)
            {
                if(residueserials.Count==0||!residueserials.Contains(residue.residueserial))
                {
                    newChain.residues.Add(residue);
                    residueserials.Add(residue.residueserial);
                }
            }
            return newChain;
        }

        /// <summary>
        /// 自动选择所有residueserial相同的氨基酸中的第一个氨基酸组成的链（select icode）
        /// </summary>
        /// <returns>Chain</returns>
        public Chain AutoSelectAltlociCode()
        {
            List<int> residueserials = new List<int>();
            Chain newChain = new Chain(this.chainID);
            foreach (Residue residue in this.residues)
            {
                if (residueserials.Count == 0 || !residueserials.Contains(residue.residueserial))
                {
                    Residue newresidue = residue.AutoSelectAltloc();
                    newChain.residues.Add(newresidue);
                    residueserials.Add(residue.residueserial);
                }
            }
            return newChain;
        }

        /// <summary>
        /// 设置链名为chainID
        /// </summary>
        /// <param name="chainID">链名</param>
        public void SetChainID(Char chainID)
        {
            this.chainID = chainID;
        }

        /// <summary>
        /// 清空链内氨基酸并设置链内氨基酸为residues
        /// </summary>
        /// <param name="residues">氨基酸链表</param>
        public void SetResidues(List<Residue> residues)
        {
            foreach(Residue residue in residues.ToList())
            {
                this.residues.Clear();
                this.residues.Add(residue);
            }
        }

        /// <summary>
        /// 复制链到新的链
        /// </summary>
        /// <returns>Chain</returns>
        public Chain Copy()
        {
            Chain newChain = new Chain(this.chainID);
            foreach(Residue residue in this.residues.ToList())
            {
                newChain.residues.Add(residue.Copy());
            }
            return newChain;
        }

        /// <summary>
        /// 按向量vector平移链
        /// </summary>
        /// <param name="vector">平移向量</param>
        public void Move(Vector3D vector)
        {
            foreach(Residue residue in this.residues)
            {
                residue.Move(vector);
            }
        }

        /// <summary>
        /// 按向量平移链到新的链
        /// </summary>
        /// <param name="vector">平移向量</param>
        /// <returns>平移后的链</returns>
        public Chain MoveToNewChain(Vector3D vector)
        {
            Chain newChain = new Chain(this.chainID);  
            foreach(Residue residue in this.residues)
            {
                newChain.residues.Add(residue.MoveToNewResidue(vector));
            }
            return newChain;
        }

        /// <summary>
        /// 按旋转矩阵旋转链
        /// </summary>
        /// <param name="RotationMatrix">旋转矩阵</param>
        public void Rotate(Double[,] RotationMatrix)
        {
            if (RotationMatrix.GetLength(0) != 3 && RotationMatrix.GetLength(1) != 3)
            {
                Console.WriteLine("rotation matrix size is not illegal");
                Environment.Exit(1);
            }
            foreach (var node in RotationMatrix)
            {
                if (node < -1 || node > 1)
                {
                    Console.WriteLine("rotation matrix value is not illegal");
                    Environment.Exit(1);
                }
            }
            foreach(Residue residue in this.residues)
            {
                residue.Rotate(RotationMatrix);
            }
        }

        /// <summary>
        /// 按旋转矩阵旋转链到新的链
        /// </summary>
        /// <param name="RotationMatrix">旋转矩阵</param>
        /// <returns>旋转后的链</returns>
        public Chain RotateToNewChain(Double[,] RotationMatrix)
        {
            if (RotationMatrix.GetLength(0) != 3 && RotationMatrix.GetLength(1) != 3)
            {
                Console.WriteLine("rotation matrix size is not illegal");
                Environment.Exit(1);
            }
            foreach (var node in RotationMatrix)
            {
                if (node < -1 || node > 1)
                {
                    Console.WriteLine("rotation matrix value is not illegal");
                    Environment.Exit(1);
                }
            }
            Chain newChain = new Chain(this.chainID);
            foreach (Residue residue in this.residues)
            {
                Residue newResidue = residue.RotateToNewResidue(RotationMatrix);
                newChain.residues.Add(newResidue);
            }
            return newChain;
        }

        /// <summary>
        /// 以静态方式计算蛋白质单链L1深度
        /// </summary>
        /// <param name="ChainData">输入的蛋白质链</param>
        /// <param name="Mode">模式，atom-atom，residue-atom，residue-residue</param>
        /// <param name="Pattern">类型，global，local</param>
        /// <param name="LocalR">local L1 depth 半径，当类型为local时</param>
        /// <returns></returns>
        public static Chain GetLoneDepth(Chain ChainData, String Mode, String Pattern, Double LocalR = 99999.0)
        {
            //Mode=atom-atom,residue-atom,residue-residue ;Pattern=global,local ;occupancy=atomdepth  tempfactor=residuedepth;
            Chain ChainResult = new Chain(ChainData.chainID);
            Chain ChainTempLast = new Chain();
            Chain ChainTempFirst = new Chain();
            int ModeNum = Mode.IndexOf('-');
            String FirstMode = Mode.Substring(0, ModeNum).ToLower();
            String LastMode = Mode.Substring(ModeNum + 1).ToLower();
            if (FirstMode == "atom")
            {
                ChainTempFirst = ChainData.Copy();
            }
            else if (FirstMode == "residue")
            {
                ChainTempFirst = ChainData.GetCAChain();
            }
            else
            {
                Console.WriteLine("mode first argument is wrong!");
                Environment.Exit(1);
            }
            if (LastMode == "atom")
            {
                ChainTempLast = ChainData.Copy();
            }
            else if (LastMode == "residue")
            {
                ChainTempLast = ChainData.GetCAChain();
            }
            else
            {
                Console.WriteLine("mode last argument is wrong!");
                Environment.Exit(1);
            }

            // Compute L1Depth
            #region head
            foreach (Residue ResFirst in ChainTempFirst.residues)
            {
                Residue ResidueResult = new Residue();
                ResidueResult = ResFirst.Copy();
                ResidueResult.atoms.Clear();
                foreach (Atom AtomFirst in ResFirst.atoms)
                {
                    Atom AtomResult = new Atom();
                    AtomResult = AtomFirst.Copy();
                    AtomResult.occupancy = 0;
                    AtomResult.tempfactor = 0;
                    int N = 0;
                    Vector3D First = new Vector3D();
                    First.X = AtomFirst.Xlaber;
                    First.Y = AtomFirst.Ylaber;
                    First.Z = AtomFirst.Zlaber;
                    Vector3D SubVector = new Vector3D();
                    SubVector.X = 0;
                    SubVector.Y = 0;
                    SubVector.Z = 0;
                    //tail node
                    #region tail
                    foreach (Residue ResLast in ChainTempLast.residues)
                    {
                        foreach (Atom AtomLast in ResLast.atoms)
                        {
                            Vector3D Last = new Vector3D();
                            Last.X = AtomLast.Xlaber;
                            Last.Y = AtomLast.Ylaber;
                            Last.Z = AtomLast.Zlaber;
                            Vector3D tempVector = Vector3D.VectorSubtract(Last, First);
                            Double norm = Vector3D.VectorNorm(tempVector);
                            if (norm == 0)
                            {
                                continue;
                            }
                            else
                            {
                                tempVector.X = tempVector.X / norm;
                                tempVector.Y = tempVector.Y / norm;
                                tempVector.Z = tempVector.Z / norm;
                            }
                            if (Pattern.ToLower() == "global")
                            {
                                SubVector = Vector3D.VectorPlus(SubVector, tempVector);
                                N++;
                            }
                            else if (Pattern.ToLower() == "local")
                            {
                                if (norm <= LocalR)
                                {
                                    SubVector = Vector3D.VectorPlus(SubVector, tempVector);
                                    N++;
                                }
                            }
                            else
                            {
                                Console.WriteLine("pattern argument is wrong!");
                                System.Environment.Exit(1);
                            }
                        }
                    }
                    #endregion tail
                    //没有自身节点，在加1
                    N++;
                    AtomResult.occupancy = 1 - Vector3D.VectorNorm(SubVector) / N;
                    if (AtomResult.atomtype == "CA")
                    {
                        AtomResult.tempfactor = AtomResult.occupancy;
                    }

                    ResidueResult.atoms.Add(AtomResult);
                }
                ChainResult.residues.Add(ResidueResult);
            }
            #endregion head
            Chain ChainOut = new Chain(ChainResult.chainID);
            foreach (Residue resnode in ChainResult.residues)
            {
                Residue newres = resnode.Copy();
                newres.atoms.Clear();
                Atom CA = resnode.GetAtom("CA");
                foreach (Atom atomnode in resnode.atoms)
                {
                    Atom newatom = atomnode.Copy();
                    newatom.tempfactor = CA.tempfactor;
                    newres.atoms.Add(newatom);
                }
                ChainOut.residues.Add(newres);
            }
            return ChainOut;
        }

        /// <summary>
        /// 计算蛋白质单链L1深度
        /// </summary>
        /// <param name="Mode">模式，atom-atom，residue-atom，residue-residue</param>
        /// <param name="Pattern">类型，global，local</param>
        /// <param name="LocalR">local L1 depth 半径，当类型为local时</param>
        /// <returns></returns>
        public Chain GetLoneDepth(String Mode, String Pattern, Double LocalR = 99999.0)
        {
            //Mode=atom-atom,residue-atom,residue-residue ;Pattern=global,local ;occupancy=atomdepth  tempfactor=residuedepth;
            Chain ChainResult = new Chain(this.chainID);
            Chain ChainTempLast = new Chain();
            Chain ChainTempFirst = new Chain();
            int ModeNum = Mode.IndexOf('-');
            String FirstMode = Mode.Substring(0, ModeNum).ToLower();
            String LastMode = Mode.Substring(ModeNum + 1).ToLower();
            if (FirstMode == "atom")
            {
                ChainTempFirst = this.Copy();
            }
            else if (FirstMode == "residue")
            {
                ChainTempFirst = this.GetCAChain();
            }
            else
            {
                Console.WriteLine("mode first argument is wrong!");
                Environment.Exit(1);
            }
            if (LastMode == "atom")
            {
                ChainTempLast = this.Copy();
            }
            else if (LastMode == "residue")
            {
                ChainTempLast = this.GetCAChain();
            }
            else
            {
                Console.WriteLine("mode last argument is wrong!");
                Environment.Exit(1);
            }

            // Compute L1Depth
            #region head
            foreach (Residue ResFirst in ChainTempFirst.residues)
            {
                Residue ResidueResult = new Residue();
                ResidueResult = ResFirst.Copy();
                ResidueResult.atoms.Clear();
                foreach (Atom AtomFirst in ResFirst.atoms)
                {
                    Atom AtomResult = new Atom();
                    AtomResult = AtomFirst.Copy();
                    AtomResult.occupancy = 0;
                    AtomResult.tempfactor = 0;
                    int N = 0;
                    Vector3D First = new Vector3D();
                    First.X = AtomFirst.Xlaber;
                    First.Y = AtomFirst.Ylaber;
                    First.Z = AtomFirst.Zlaber;
                    Vector3D SubVector = new Vector3D();
                    SubVector.X = 0;
                    SubVector.Y = 0;
                    SubVector.Z = 0;
                    //tail node
                    #region tail
                    foreach (Residue ResLast in ChainTempLast.residues)
                    {
                        foreach (Atom AtomLast in ResLast.atoms)
                        {
                            Vector3D Last = new Vector3D();
                            Last.X = AtomLast.Xlaber;
                            Last.Y = AtomLast.Ylaber;
                            Last.Z = AtomLast.Zlaber;
                            Vector3D tempVector = Vector3D.VectorSubtract(Last, First);
                            Double norm = Vector3D.VectorNorm(tempVector);
                            if (norm == 0)
                            {
                                continue;
                            }
                            else
                            {
                                tempVector.X = tempVector.X / norm;
                                tempVector.Y = tempVector.Y / norm;
                                tempVector.Z = tempVector.Z / norm;
                            }
                            if (Pattern.ToLower() == "global")
                            {
                                SubVector = Vector3D.VectorPlus(SubVector, tempVector);
                                N++;
                            }
                            else if (Pattern.ToLower() == "local")
                            {
                                if (norm <= LocalR)
                                {
                                    SubVector = Vector3D.VectorPlus(SubVector, tempVector);
                                    N++;
                                }
                            }
                            else
                            {
                                Console.WriteLine("pattern argument is wrong!");
                                System.Environment.Exit(1);
                            }
                        }
                    }
                    #endregion tail
                    //没有自身节点，在加1
                    N++;
                    AtomResult.occupancy = 1 - Vector3D.VectorNorm(SubVector) / N;
                    if (AtomResult.atomtype == "CA")
                    {
                        AtomResult.tempfactor = AtomResult.occupancy;
                    }
                    ResidueResult.atoms.Add(AtomResult);
                    Console.WriteLine(AtomResult.occupancy);
                }
                ChainResult.residues.Add(ResidueResult);
            }
            #endregion head
            Chain ChainOut = new Chain(ChainResult.chainID);
            foreach (Residue resnode in ChainResult.residues)
            {
                Residue newres = resnode.Copy();
                newres.atoms.Clear();
                Atom CA = resnode.GetAtom("CA");
                foreach (Atom atomnode in resnode.atoms)
                {
                    Atom newatom = atomnode.Copy();
                    newatom.tempfactor = CA.tempfactor;
                    newres.atoms.Add(newatom);
                }
                ChainOut.residues.Add(newres);
            }
            return ChainOut;
        }

        public List<Double> GetPsiList()
        {
            //n length chain, n-1 phi, n-1 psi
            //      1      2     3      4      5            i-1     i      i+1       n-1    n
            // H-N-CA-C-N-CA-C-N-CA-C-N-CA-C-N-CA-C-......-N-CA-C-N-CA-C-N-CA-C...-N-CA-C-N-CA-C=O    psi=s phi=h
            //      s1   h1 s2  h2 s3  h3 s4 h4  s5                                  sn-1  hn-1
            // psi_1 in residue_1    psi_n-1 in residue_n-1   phi_1 in residue_2  phi_n-1 in residue_n
            // only return n-2 psi ignore s1, from s2 (in r2) to s_n-1 (in r_n-1)
            List<Double> PsiList = new List<Double>();
            
            for (int i=1;i<this.residues.Count-1;i++)
            {
                Residue current = this.residues.ElementAt(i);
                Residue next = this.residues.ElementAt(i+1);
                Atom N = current.atoms.FirstOrDefault(a => a.atomtype == "N");
                Atom CA = current.atoms.FirstOrDefault(a => a.atomtype == "CA");
                Atom C = current.atoms.FirstOrDefault(a => a.atomtype == "C");
                Atom nextN = next.atoms.FirstOrDefault(a => a.atomtype == "N");
                if (N != null && CA != null && C != null && nextN != null)
                {
                    Double disCN = C.GetDistance(nextN);
                    //if (disCN < 1.5)// for broken chain
                    {
                        Double psi = Point3D.DiheDralAngle(N.ToPoint3D(), CA.ToPoint3D(), C.ToPoint3D(), nextN.ToPoint3D());
                        if(psi==Common.impossibleMaxDouble)
                        {
                            psi = 2.0*Math.PI;
                        }
                        if(psi==Common.impossibleMinDouble)
                        {
                            psi = -2.0 * Math.PI;
                        }
                        PsiList.Add(psi);
                    }
                }
                //else
                //{
                //    Console.WriteLine(current.residueserial);
                //}
            }

            return PsiList;
        }

        public List<Double> GetPhiList()
        {
            //n length chain, n-1 phi, n-1 psi
            //      1      2     3      4      5            i-1     i      i+1       n-1    n
            // H-N-CA-C-N-CA-C-N-CA-C-N-CA-C-N-CA-C-......-N-CA-C-N-CA-C-N-CA-C...-N-CA-C-N-CA-C=O    psi=s phi=h
            //      s1   h1 s2  h2 s3  h3 s4 h4  s5                                  sn-1  hn-1
            // psi_1 in residue_1    psi_n-1 in residue_n-1   phi_1 in residue_2  phi_n-1 in residue_n
            // return n-2 phi, ignore phi_n-1, from h1 (in r2) to h_n-2 (in r_n-1) 
            List<Double> PhiList = new List<Double>();

            for(int i=1;i<this.residues.Count-1;i++)
            {
                Residue previous = this.residues.ElementAt(i - 1);
                Residue current = this.residues.ElementAt(i);
                Atom previousC = previous.atoms.FirstOrDefault(a => a.atomtype == "C");
                Atom N = current.atoms.FirstOrDefault(a => a.atomtype == "N");
                Atom CA = current.atoms.FirstOrDefault(a => a.atomtype == "CA");
                Atom C = current.atoms.FirstOrDefault(a => a.atomtype == "C");
                if(previousC!=null&&N!=null&&CA!=null&&C!=null)
                {
                    Double disCN = previousC.GetDistance(N);
                    //if (disCN < 1.5)// for broken chain
                    {
                        Double phi = Point3D.DiheDralAngle(previousC.ToPoint3D(), N.ToPoint3D(), CA.ToPoint3D(), C.ToPoint3D());
                        if(phi==Common.impossibleMaxDouble)
                        {
                            phi = 2.0 * Math.PI;
                        }
                        if(phi==Common.impossibleMinDouble)
                        {
                            phi = -2.0 * Math.PI;
                        }
                        PhiList.Add(phi);
                    }
                }
            }

            return PhiList;

        }

        public void CalculateRamaPsiPhi()
        {
            //n length chain, n-1 phi, n-1 psi
            //      1      2     3      4      5            i-1     i      i+1       n-1    n
            // H-N-CA-C-N-CA-C-N-CA-C-N-CA-C-N-CA-C-......-N-CA-C-N-CA-C-N-CA-C...-N-CA-C-N-CA-C=O    psi=s phi=h
            //      s1   h1 s2  h2 s3  h3 s4 h4  s5                                  sn-1  hn-1
            // psi_1 in residue_1    psi_n-1 in residue_n-1   phi_1 in residue_2  phi_n-1 in residue_n
            // return n-2 phi, ignore phi_n-1, from h1 (in r2) to h_n-2 (in r_n-1) 
            for (int i=1;i<this.residues.Count-1;i++)
            {
                Residue previous = this.residues.ElementAt(i - 1);
                Residue current = this.residues.ElementAt(i);
                Residue next = this.residues.ElementAt(i + 1);

                Atom previousC = previous.atoms.FirstOrDefault(a => a.atomtype == "C");
                Atom N = current.atoms.FirstOrDefault(a => a.atomtype == "N");
                Atom CA = current.atoms.FirstOrDefault(a => a.atomtype == "CA");
                Atom C = current.atoms.FirstOrDefault(a => a.atomtype == "C");
                Atom nextN = next.atoms.FirstOrDefault(a => a.atomtype == "N");

                if (N != null && CA != null && C != null && nextN != null)
                {
                    Double disCN = C.GetDistance(nextN);
                    //if (disCN < 1.5)// for broken chain
                    {
                        Double psi = Point3D.DiheDralAngle(N.ToPoint3D(), CA.ToPoint3D(), C.ToPoint3D(), nextN.ToPoint3D());
                        if (psi == Common.impossibleMaxDouble)
                        {
                            psi = 2.0 * Math.PI;
                        }
                        if (psi == Common.impossibleMinDouble)
                        {
                            psi = -2.0 * Math.PI;
                        }
                        this.residues.ElementAt(i).psi = psi;
                    }
                }
                if (previousC != null && N != null && CA != null && C != null)
                {
                    Double disCN = previousC.GetDistance(N);
                    //if (disCN < 1.5)// for broken chain
                    {
                        Double phi = Point3D.DiheDralAngle(previousC.ToPoint3D(), N.ToPoint3D(), CA.ToPoint3D(), C.ToPoint3D());
                        if (phi == Common.impossibleMaxDouble)
                        {
                            phi = 2.0 * Math.PI;
                        }
                        if (phi == Common.impossibleMinDouble)
                        {
                            phi = -2.0 * Math.PI;
                        }
                        this.residues.ElementAt(i).phi = phi;
                    }
                }
            }
        }


        public List<Char> GetPdbFasta()
        {
            List<Char> fasta = new List<Char>();
            foreach (Residue node in this.residues)
            {
                fasta.Add(node.GetResidueType());
            }
            return fasta;
        }
         

    }
}
