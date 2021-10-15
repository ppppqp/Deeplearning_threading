using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace BioNet
{
    public class Protein
    {
        //member
        public String proteinname;
        public List<Chain> chains;
        //function

        /// <summary>
        /// 初始化Protein并初始化chains
        /// </summary>
        public Protein()
        {
            this.chains = new List<Chain>();
        }

        /// <summary>
        /// 初始化Protein并命名为proteinname出，初始化chains
        /// </summary>
        /// <param name="proteinname">蛋白质名字</param>
        public Protein(String proteinname)
        {
            this.proteinname = proteinname;
            this.chains = new List<Chain>();
        }

        /// <summary>
        /// 初始化Protein，命名为proteinname，设置蛋白链为chains
        /// </summary>
        /// <param name="proteinname">蛋白质名字</param>
        /// <param name="chains">蛋白质链</param>
        public Protein(String proteinname,List<Chain> chains)
        {
            this.proteinname = proteinname;
            this.chains = new List<Chain>();
            foreach(Chain chain in chains.ToList())
            {
                this.chains.Add(chain);
            }
        }

        /// <summary>
        /// 初始化蛋白，命名为proteinname
        /// </summary>
        /// <param name="pdb">stream reader，PDB format file</param>
        /// <param name="proteinname">蛋白质名字</param>
        public Protein(StreamReader pdb,String proteinname)
        {
            this.proteinname = proteinname;
            this.chains = new List<Chain>();
            while(pdb.Peek()>-1)
            {
                String line = pdb.ReadLine();

                #region Format ATOM
                if (line.StartsWith("ATOM"))
                {
                    //chain info
                    Char chainid = line.ElementAt(21);
                    //resisue info
                    String residuename= line.Substring(17, 3).Trim();
                    int residueserial= Convert.ToInt32(line.Substring(22, 4));
                    Char icode = line.ElementAt(26);
                    //atom info
                    Atom atomnode = new Atom();
                    atomnode.atomserial = Convert.ToInt32(line.Substring(6, 5));
                    atomnode.atomtype = line.Substring(12, 4).Trim();
                    atomnode.altLoc = line.ElementAt(16);
                    atomnode.atom = Convert.ToChar(atomnode.atomtype.ElementAt(0));
                    atomnode.Xlaber = Convert.ToDouble(line.Substring(30, 8));
                    atomnode.Ylaber = Convert.ToDouble(line.Substring(38, 8));
                    atomnode.Zlaber = Convert.ToDouble(line.Substring(46, 8));
                    if (line.Length > 60)
                    {
                        atomnode.occupancy = Convert.ToDouble(line.Substring(54, 6));
                    }

                    if (line.Length > 66)
                    {
                        atomnode.tempfactor = Convert.ToDouble(line.Substring(60, 6));
                    }

                    if (line.Count() >= 78 && line.Substring(76, 2).Trim().Length != 0)
                    {
                        atomnode.element = line.Substring(76, 2).Trim();
                    }
                    else if (!Char.IsNumber(atomnode.atomtype.ElementAt(0)))
                    {
                        atomnode.element = Convert.ToString(atomnode.atomtype.ElementAt(0));
                    }
                    else
                    {
                        atomnode.element = Convert.ToString(atomnode.atomtype.ElementAt(1));
                    }

                    if (line.Trim().Count() >= 80)
                    {
                        atomnode.charge = line.Substring(78, 2).Trim();
                    }
                    else
                    {
                        atomnode.charge = "";
                    }
                    //format data
                    if (this.chains.Count == 0 || this.chains.Last().chainID != chainid || this.chains.Last().residues.Last().residueserial > residueserial)
                    {
                        Chain chain = new Chain(chainid);
                        Residue residue = new Residue(residuename, residueserial, icode);
                        residue.atoms.Add(atomnode);
                        chain.residues.Add(residue);
                        this.chains.Add(chain);
                    }
                    else
                    {
                        if(this.chains.Last().residues.Last().residueserial==residueserial)
                        {
                            Residue residue = this.chains.Last().residues.Find(node => node.residueserial == residueserial && node.iCode == icode);
                            if(residue!=null)
                            {
                                residue.atoms.Add(atomnode);
                            }
                            else
                            {
                                Residue newresidue = new Residue(residuename, residueserial, icode);
                                newresidue.atoms.Add(atomnode);
                                this.chains.Last().residues.Add(newresidue);
                            }
                        }
                        else
                        {
                            Residue newresidue = new Residue(residuename, residueserial, icode);
                            newresidue.atoms.Add(atomnode);
                            this.chains.Last().residues.Add(newresidue);
                        }
                    }      
                }
                #endregion Format ATOM

            }
        }

        /// <summary>
        /// 返回蛋白质的名字
        /// </summary>
        /// <returns>String，蛋白质名字</returns>
        public String GetProteinName()
        {
            return this.proteinname;
        }

        /// <summary>
        /// 返回蛋白质中所有的链
        /// </summary>
        /// <returns>List，蛋白质链</returns>
        public List<Chain> GetChains()
        {
            return this.chains;
        }

        /// <summary>
        /// 返回所有以chainID为id的链
        /// </summary>
        /// <param name="chainID">链id</param>
        /// <returns>List，蛋白质链</returns>
        public List<Chain> GetChains(Char chainID)
        {
            return this.chains.FindAll(node => node.chainID == chainID).ToList();
        }

        /// <summary>
        /// 返回蛋白质中所有以chainID为id的第一条链
        /// </summary>
        /// <param name="chainID">链id</param>
        /// <returns>Chain，蛋白质链</returns>
        public Chain GetChain(Char chainID)
        {
            return this.chains.Find(node => node.chainID == chainID);
        }

        /// <summary>
        /// 自动选择氨基酸icode，所有residue serial相同中的第一个，返回到新的Protein
        /// </summary>
        /// <returns>蛋白质</returns>
        public Protein AutoSelectiCode()
        {
            Protein newProtein = new Protein(this.proteinname);
            foreach(Chain chain in this.chains)
            {
                Chain newChain = chain.AutoSelectiCode();
                newProtein.chains.Add(newChain);
            }
            return newProtein;
        }

        /// <summary>
        /// 设置蛋白质名字
        /// </summary>
        /// <param name="proteinname">蛋白质名字</param>
        public void SetProteinName(String proteinname)
        {
            this.proteinname = proteinname;
        }

        /// <summary>
        /// 清空蛋白全部链并设置蛋白链为chains
        /// </summary>
        /// <param name="chains">蛋白质链</param>
        public void SetChains(List<Chain> chains)
        {
            foreach(Chain chain in chains.ToList())
            {
                this.chains.Clear();
                this.chains.Add(chain);
            }
        }

        /// <summary>
        /// 复制蛋白到新的蛋白
        /// </summary>
        /// <returns>蛋白质</returns>
        public Protein Copy()
        {
            Protein newProtein = new Protein(this.proteinname);
            foreach(Chain chain in this.chains)
            {
                Chain newChain = chain.Copy();
                newProtein.chains.Add(newChain);
            }
            return newProtein;
        }

        /// <summary>
        /// 返回NMR结构中所有相同链中的第一条链
        /// </summary>
        /// <returns>去重后的蛋白质</returns>
        public Protein NMRTrimChain()
        {
            Protein newProtein = new Protein(this.proteinname);
            List<Char> chainIDs = new List<Char>();
            foreach(Chain chain in this.chains)
            {
                if(chainIDs.Count==0||!chainIDs.Contains(chain.chainID))
                {
                    newProtein.chains.Add(chain);
                    chainIDs.Add(chain.chainID);
                }
            }
            return newProtein;
        }

        /// <summary>
        /// 按向量平移蛋白
        /// </summary>
        /// <param name="vector">平移向量</param>
        public void Move(Vector3D vector)
        {
            foreach(Chain chain in this.chains)
            {
                chain.Move(vector);
            }
        }

        /// <summary>
        /// 按向量平移蛋白返回到新的蛋白
        /// </summary>
        /// <param name="vector">平移向量</param>
        /// <returns>平移后的蛋白质</returns>
        public Protein MoveToNewProtein(Vector3D vector)
        {
            Protein newProtein = new Protein(this.proteinname);
            foreach(Chain chain in this.chains)
            {
                Chain newChain = chain.MoveToNewChain(vector);
                newProtein.chains.Add(newChain);
            }
            return newProtein;
        }

        /// <summary>
        /// 按旋转矩阵旋转蛋白
        /// </summary>
        /// <param name="RotationMatrix">旋转矩阵</param>
        public void Rotate(Double [,] RotationMatrix)
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
            foreach(Chain chain in this.chains)
            {
                chain.Rotate(RotationMatrix);
            }
        }

        /// <summary>
        /// 按旋转矩阵旋转蛋白返回到新的蛋白
        /// </summary>
        /// <param name="RotationMatrix">旋转矩阵</param>
        /// <returns>旋转后的蛋白质</returns>
        public Protein RotateToNewProtein(Double [,] RotationMatrix)
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
            Protein newProtein = new Protein(this.proteinname);
            foreach(Chain chain in this.chains)
            {
                Chain newChain = chain.RotateToNewChain(RotationMatrix);
                newProtein.chains.Add(newChain);
            }
            return newProtein;
        }

        /// <summary>
        /// 将当前蛋白写入pdb文件
        /// </summary>
        /// <param name="pdb">stream writer，pdb</param>
        public void WritePDB(StreamWriter pdb)
        {
            foreach (Chain chainnode in this.chains)
            {
                foreach (Residue residuenode in chainnode.residues)
                {
                    foreach (Atom atomnode in residuenode.atoms)
                    {
                        pdb.Write("ATOM  ");
                        pdb.Write("{0,5}", atomnode.atomserial);
                        pdb.Write(" ");
                        if (atomnode.atomtype.Length == 4 || Char.IsNumber(atomnode.atomtype.ElementAt(0)))
                        {
                            pdb.Write(atomnode.atomtype.PadRight(4, ' '));
                        }
                        else
                        {
                            pdb.Write(" ");
                            pdb.Write(atomnode.atomtype.PadRight(3, ' '));
                        }
                        pdb.Write(atomnode.altLoc);
                        pdb.Write(residuenode.residuename);
                        pdb.Write(" ");
                        pdb.Write(chainnode.chainID);
                        pdb.Write("{0,4}", residuenode.residueserial);
                        pdb.Write(residuenode.iCode);
                        pdb.Write("   ");
                        pdb.Write("{0,8:f3}", atomnode.Xlaber);
                        pdb.Write("{0,8:f3}", atomnode.Ylaber);
                        pdb.Write("{0,8:f3}", atomnode.Zlaber);
                        pdb.Write("{0,6:f2}", atomnode.occupancy);
                        pdb.Write("{0,6:f2}", atomnode.tempfactor);
                        pdb.Write("          ");
                        pdb.Write(atomnode.element.PadLeft(2, ' '));
                        pdb.Write(atomnode.charge.PadLeft(2, ' '));
                        pdb.WriteLine();
                    }
                }
            }
        }

        /// <summary>
        /// 以静态方式将Protein pdb写入sw文件，pdb format file 
        /// </summary>
        /// <param name="pdb">Protein pdb</param>
        /// <param name="sw">StreamWriter sw</param>
        public static void WritePDB(Protein pdb,StreamWriter sw)
        {
            pdb.WritePDB(sw);
        }

        /// <summary>
        /// 将蛋白的坐标转化为point3D list
        /// </summary>
        /// <returns></returns>
        public List<Point3D> ToPoint3DList()
        {
            List<Point3D> ProteinPoints = new List<Point3D>();
            foreach(var chain in this.chains)
            {
                foreach(var residue in chain.residues)
                {
                    foreach(var atom in residue.atoms)
                    {
                        ProteinPoints.Add(atom.ToPoint3D());
                    }
                }
            }
            return ProteinPoints;
        }

    }
}
