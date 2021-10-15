using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CommonFunction;

namespace BioNet
{
    public class Residue
    {
        //class member
        public String residuename;
        public Char residuetype;
        public int residueserial;//key(can't decide all residue)
        public Char iCode;//second key
        public List<Atom> atoms;
        public Double phi;
        public Double psi;
        //function

        /// <summary>
        /// 初始化氨基酸
        /// </summary>
        public Residue()
        {
            this.atoms = new List<Atom>();
            this.iCode = ' ';
            this.phi = Common.impossibleMaxDouble;
            this.psi = Common.impossibleMaxDouble;
        }

        /// <summary>
        /// 初始化氨基酸，氨基酸名字为residuename，icode设置为iCode
        /// </summary>
        /// <param name="residuename">氨基酸名字</param>
        /// <param name="iCode">icode</param>
        public Residue(String residuename,Char iCode)
        {
            this.residuename = residuename;
            this.residuetype = InitResidueType(residuename);
            this.iCode = iCode;
            this.atoms = new List<Atom>();
            this.phi = Common.impossibleMaxDouble;
            this.psi = Common.impossibleMaxDouble;
        }

        /// <summary>
        /// 初始化氨基酸，氨基酸序号为residueserial，氨基酸icode设置为iCode
        /// </summary>
        /// <param name="residueserial">氨基酸序号</param>
        /// <param name="iCode">icode</param>
        public Residue(int residueserial,Char iCode=' ')
        {
            this.residueserial = residueserial;
            this.iCode = iCode;
            this.atoms = new List<Atom>();
            this.phi = Common.impossibleMaxDouble;
            this.psi = Common.impossibleMaxDouble;
        }

        /// <summary>
        /// 初始化氨基酸，氨基酸序号为residueserial，氨基酸icode设置为iCode，氨基酸名字为residuename
        /// </summary>
        /// <param name="residuename">氨基酸名字</param>
        /// <param name="residueserial">氨基酸序号</param>
        /// <param name="iCode">icode</param>
        public Residue(String residuename, int residueserial, Char iCode=' ')
        {
            this.residuename = residuename;
            this.residuetype = InitResidueType(residuename);
            this.residueserial = residueserial;
            this.iCode = iCode;
            this.atoms = new List<Atom>();
            this.phi = Common.impossibleMaxDouble;
            this.psi = Common.impossibleMaxDouble;
        }

        /// <summary>
        /// 初始化氨基酸，氨基酸序号为residueserial，氨基酸icode设置为iCode，氨基酸名字为residuename，设置原子为atoms
        /// </summary>
        /// <param name="residuename">氨基酸名字</param>
        /// <param name="residueserial">氨基酸序号</param>
        /// <param name="atoms">氨基酸原子</param>
        /// <param name="iCode">icode</param>
        public Residue(String residuename, int residueserial, List<Atom> atoms, Char iCode = ' ')
        {
            this.residuename = residuename;
            this.residuetype = InitResidueType(residuename);
            this.residueserial = residueserial;
            this.iCode = iCode;
            this.atoms = new List<Atom>();
            foreach (Atom atom in atoms.ToList())
            {
                this.atoms.Add(atom.Copy());
            }
            this.phi = Common.impossibleMaxDouble;
            this.psi = Common.impossibleMaxDouble;
        }

        /// <summary>
        /// 返回氨基酸名字
        /// </summary>
        /// <returns>氨基酸名字</returns>
        public String GetResidueName()
        {
            return this.residuename;
        }

        /// <summary>
        /// 返回氨基酸序号
        /// </summary>
        /// <returns>氨基酸序号</returns>
        public int GetResidueSerial()
        {
            return this.residueserial;
        }

        /// <summary>
        /// 返回氨基酸icode
        /// </summary>
        /// <returns>icode</returns>
        public Char GetiCode()
        {
            return this.iCode;
        }

        /// <summary>
        /// 返回氨基酸中所有的原子
        /// </summary>
        /// <returns>氨基酸原子</returns>
        public List<Atom> GetAtoms()
        {
            return this.atoms;
        }

        /// <summary>
        /// 返回原子序号为atomserial的第一个原子，没有返回null
        /// </summary>
        /// <param name="atomserial">原子序号</param>
        /// <returns>原子</returns>
        public Atom GetAtom(int atomserial)
        {
            Atom atom = this.atoms.Find(node => node.atomserial == atomserial);
            return atom;
        }

        /// <summary>
        /// 返回Type为atomtype的第一个原子，没有返回null
        /// </summary>
        /// <param name="atomtype"></param>
        /// <returns>原子</returns>
        public Atom GetAtom(String atomtype)
        {
            Atom atom = this.atoms.Find(node => node.atomtype == atomtype);
            return atom;
        }

        /// <summary>
        /// 返回氨基酸中的主链原子CA，C，N，O，返回到新氨基酸，并不一定完整
        /// </summary>
        /// <returns>主链原子组成的氨基酸</returns>
        public Residue GetMainChainAtoms()
        {
            Residue newResidue = new Residue(this.residuename, this.residueserial, this.iCode);
            Atom C = GetAtom("C");
            if (C != null)
            {
                newResidue.atoms.Add(C);
            }
            Atom N = GetAtom("N");
            if (N != null)
            {
                newResidue.atoms.Add(N);
            }
            Atom CA = GetAtom("CA");
            if (CA != null)
            {
                newResidue.atoms.Add(CA);
            }
            Atom O = new Atom();
			//Atom O = GetAtom("O");
            if (O != null)
            {
                newResidue.atoms.Add(O);
            }
            return newResidue;
        }

        /// <summary>
        /// 返回氨基酸侧链原子到新氨基酸
        /// </summary>
        /// <returns>侧链原子组成的氨基酸</returns>
        public Residue GetSideChainAtoms()
        {
            List<Atom> allatoms = this.atoms;
            List<Atom> mainchain = this.GetMainChainAtoms().atoms;
            List<Atom> sidechain = allatoms.Except(mainchain).ToList();
            Residue newResidue = new Residue(this.residuename, this.residueserial, sidechain, this.iCode);
            return newResidue;
        }

        /// <summary>
        /// 设置氨基酸名字的residuename
        /// </summary>
        /// <param name="residuename">氨基酸名字</param>
        public void SetResidueName(String residuename)
        {
            this.residuename = residuename;
			this.residuetype = InitResidueType(residuename);
        }

        /// <summary>
        /// 设置氨基酸序号为residueserial
        /// </summary>
        /// <param name="residueserial">氨基酸序号</param>
        public void SetResidueSerial(int residueserial)
        {
            this.residueserial = residueserial;
        }

        /// <summary>
        /// 设置氨基酸icode为iCode
        /// </summary>
        /// <param name="iCode">icode</param>
        public void SetiCode(Char iCode)
        {
            this.iCode = iCode;
        }

        /// <summary>
        /// 设置氨基酸原子为atoms
        /// </summary>
        /// <param name="atoms">氨基酸原子</param>
        public void SetAtoms(List<Atom> atoms)
        {
            this.atoms.Clear();
            foreach (Atom atom in atoms.ToList())
            {
                this.atoms.Add(atom.Copy());
            }
        }

        /// <summary>
        /// 判断当前氨基酸主链原子是否完整
        /// </summary>
        /// <returns>布尔值，完整true，不完整false</returns>
        public bool IsMainChainComplete()
        {
            bool complete = false;
            bool contain = true;
            if (this.atoms.Find(node => node.atomtype == "CA") == null)
            {
                contain = false;
            }
            if (this.atoms.Find(node => node.atomtype == "C") == null)
            {
                contain = false;
            }
            if (this.atoms.Find(node => node.atomtype == "O") == null)
            {
                contain = false;
            }
            if (this.atoms.Find(node => node.atomtype == "N") == null)
            {
                contain = false;
            }
            if (contain)
            {
                complete = true;
            }
            return complete;
        }

        /// <summary>
        /// 复制当前氨基酸为新氨基酸
        /// </summary>
        /// <returns>复制后的氨基酸</returns>
        public Residue Copy()
        {
            Residue residue = new Residue();
            residue.residuename = this.residuename;
            residue.residueserial = this.residueserial;
            residue.iCode = this.iCode;
            foreach (Atom atom in this.atoms.ToList())
            {
                residue.atoms.Add(atom.Copy());
            }
            return residue;
        }

        /// <summary>
        /// 返回当前氨基酸与另一氨基酸的CA距离，如果其中一个CA不存在，返回-9999.0
        /// </summary>
        /// <param name="residue">目标氨基酸</param>
        /// <returns>CA距离</returns>
        public Double GetCADistance(Residue residue)
        {
            Atom thisCA = this.GetAtom("CA");
            Atom residueCA = residue.GetAtom("CA");
            if(thisCA==null||residueCA==null)
            {
                return -9999.0;
            }
            return thisCA.GetDistance(residueCA);
        }

        /// <summary>
        /// 以静态方式返回两个氨基酸的CA距离，如果其中一个CA不存在，返回-9999.0
        /// </summary>
        /// <param name="residueOne">氨基酸1</param>
        /// <param name="residueTwo">氨基酸2</param>
        /// <returns>CA距离</returns>
        public static Double GetCADistance(Residue residueOne, Residue residueTwo)
        {
            return residueOne.GetCADistance(residueTwo);
        }

        /// <summary>
        /// 返回当前氨基酸和目标氨基酸间的最小距离
        /// </summary>
        /// <param name="residue">目标氨基酸</param>
        /// <returns>最小距离</returns>
        public Double GetLeastDistance(Residue residue)
        {
            Double leastDistance = 9999;
            foreach(Atom thisAtom in this.atoms)
            {
                foreach(Atom atom in residue.atoms)
                {
                    Double newDistance = thisAtom.GetDistance(atom);
                    if (newDistance<leastDistance)
                    {
                        leastDistance = newDistance;
                    }
                }
            }
            return leastDistance;
        }

        /// <summary>
        /// 以静态方式返回两个氨基酸的最小距离
        /// </summary>
        /// <param name="residueOne">氨基酸1</param>
        /// <param name="residueTwo">氨基酸2</param>
        /// <returns>最小距离</returns>
        public static Double GetLeastDistance(Residue residueOne,Residue residueTwo)
        {
            return residueOne.GetLeastDistance(residueTwo);
        }

        /// <summary>
        /// 按向量vector平移当前氨基酸到新的位置
        /// </summary>
        /// <param name="vector">平移向量</param>
        public void Move(Vector3D vector)
        {
            foreach(Atom atom in this.atoms)
            {
                atom.Move(vector);
            }
        }

        /// <summary>
        /// 按平移向量vector平移当前氨基酸到新的氨基酸
        /// </summary>
        /// <param name="vector">平移向量</param>
        /// <returns>平移后的氨基酸</returns>
        public Residue MoveToNewResidue(Vector3D vector)
        {
            Residue newResidue = this.Copy();
            newResidue.atoms.Clear();
            foreach(Atom atom in this.atoms)
            {
                newResidue.atoms.Add(atom.MoveToNewAtom(vector));
            }
            return newResidue;
        }

        /// <summary>
        /// 按旋转矩阵RotationMatrix旋转当前氨基酸到新的位置
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
            foreach(Atom atom in this.atoms)
            {
                atom.Rotate(RotationMatrix);
            }
        }

        /// <summary>
        /// 按旋转矩阵RotationMatrix旋转当前氨基酸到新的位置，并返回到新的氨基酸
        /// </summary>
        /// <param name="RotationMatrix">旋转矩阵</param>
        /// <returns>旋转后的氨基酸</returns>
        public Residue RotateToNewResidue(Double[,] RotationMatrix)
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
            Residue newResidue = this.Copy();
            newResidue.atoms.Clear();
            foreach (Atom atom in this.atoms)
            {
                Atom newAtom = atom.RotateToNewAtom(RotationMatrix);
                newResidue.atoms.Add(newAtom);
            }
            return newResidue;
        }

        /// <summary>
        /// 以静态方式返回三个原子形成的键角的余弦值[-1,1]
        /// </summary>
        /// <param name="one">边1上的原子</param>
        /// <param name="common">顶角原子</param>
        /// <param name="two">边2上的原子</param>
        /// <returns>键角的余弦值</returns>
        public static Double CosBondAngle(Atom one,Atom common, Atom two)
        {
            Vector3D vc1 = new Vector3D(one, common);
            Vector3D vc2 = new Vector3D(two, common);
            return Vector3D.VectorCosAngle(vc1, vc2);
        }

        /// <summary>
        /// 以静态方式返回三个原子形成的键角[0,pi]
        /// </summary>
        /// <param name="one">边1上的原子</param>
        /// <param name="common">顶角原子</param>
        /// <param name="two">边2上的原子</param>
        /// <returns>键角</returns>
        public static Double BondAngle(Atom one, Atom common, Atom two)
        {
            Vector3D vc1 = new Vector3D(one, common);
            Vector3D vc2 = new Vector3D(two, common);
            return Vector3D.VectorAngle(vc1, vc2);
        }

        /// <summary>
        /// 以静态方式返回四个原子形成的二面角的余弦值[-1,1]
        /// </summary>
        /// <param name="innerA1">平面1上的点</param>
        /// <param name="commonA1">棱上的靠近innerA1点</param>
        /// <param name="commonA2">棱上的靠近innerA2点</param>
        /// <param name="innerA2">平面2上的点</param>
        /// <returns>二面角余弦值</returns>
        public static Double CosDihedralAngle(Atom innerA1,Atom commonA1,Atom commonA2,Atom innerA2)
        {
            Point3D innerP1 = new Point3D(innerA1);
            Point3D commonP1 = new Point3D(commonA1);
            Point3D commonP2 = new Point3D(commonA2);
            Point3D innerP2 = new Point3D(innerA2);
            return Point3D.CosDiheDralAngle(innerP1, commonP1, commonP2, innerP2);
        }

        /// <summary>
        /// 以静态方式返回四个原子形成的二面角[0,pi]
        /// </summary>
        /// <param name="innerA1">平面1上的点</param>
        /// <param name="commonA1">棱上的靠近innerA1点</param>
        /// <param name="commonA2">棱上的靠近innerA2点</param>
        /// <param name="innerA2">平面2上的点</param>
        /// <returns>二面角</returns>
        public static Double DihedralAngle(Atom innerA1, Atom commonA1, Atom commonA2, Atom innerA2)
        {
            Point3D innerP1 = new Point3D(innerA1);
            Point3D commonP1 = new Point3D(commonA1);
            Point3D commonP2 = new Point3D(commonA2);
            Point3D innerP2 = new Point3D(innerA2);
            return Point3D.DiheDralAngle(innerP1, commonP1, commonP2, innerP2);
        }

        Char InitResidueType(String residuename)
        {
            Char restype = 'X';
            int index = Common.AAtype.ToList().IndexOf(residuename);
            if (index != -1)
            {
                restype = Common.AAtypeChar[index];
                return restype;
            }
            else
            {
                return restype;
            }
        }

        public Char GetResidueType()
        {
            return this.residuetype;
        }

        public Double GetPhi()
        {
            return this.phi;
        }
        public Double GetPsi()
        {
            return this.psi;
        }

        public bool IsContact(Residue residue)
        {
            Double distanceCutoff = 8.0;
            bool contact = false;
            Atom thisContactAtom = new Atom();
            if(this.residuename=="GLY")
            {
                thisContactAtom = this.GetAtom("CA");
            }
            else
            {
                thisContactAtom = this.GetAtom("CB");
            }
            Atom residueContactAtom = new Atom();
            if(residue.residuename=="GLY")
            {
                residueContactAtom = residue.GetAtom("CA");
            }
            else
            {
                residueContactAtom = residue.GetAtom("CB");
            }
            if(thisContactAtom==null||residueContactAtom==null)
            {
                return contact;
            }
            Double distance = thisContactAtom.GetDistance(residueContactAtom);
            if(distance<=distanceCutoff)
            {
                contact = true;
            }
            return contact;
        }

        public Residue AutoSelectAltloc()
        {
            Residue NewResidue = this.Copy();
            NewResidue.atoms.Clear();
            char pushcode = ' ';
            int breakpoint = this.atoms.Count+1;
            for(int i=0;i<this.atoms.Count;i++)
            {
                if(this.atoms.ElementAt(i).altLoc==' ')
                {
                    NewResidue.atoms.Add(this.atoms.ElementAt(i));
                }
                else
                {
                    pushcode = this.atoms.ElementAt(i).altLoc;
                    breakpoint = i;
                    break;
                }
            }
            for(int j=breakpoint;j<this.atoms.Count;j++)
            {
                if(this.atoms.ElementAt(j).altLoc==pushcode)
                {
                    NewResidue.atoms.Add(this.atoms.ElementAt(j));
                }
            }

            return NewResidue;
        }
    }
}
