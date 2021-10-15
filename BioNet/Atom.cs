using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BioNet
{
    public class Atom
    {
        //class member
        public Char atom;
        public String atomtype;
        public int atomserial;
        public Char altLoc;
        public Double Xlaber;
        public Double Ylaber;
        public Double Zlaber;
        public Double occupancy;
        public Double tempfactor;
        public String element;
        public String charge;
        //function

        /// <summary>
        /// 初始化原子
        /// </summary>
        public Atom() { }

        /// <summary>
        /// 初始化原子
        /// </summary>
        /// <param name="atom">原子</param>
        /// <param name="atomtype">原子类型</param>
        /// <param name="atomserial">原子序号</param>
        /// <param name="altLoc"></param>
        /// <param name="Xlaber">X坐标</param>
        /// <param name="Ylaber">Y坐标</param>
        /// <param name="Zlaber">Z坐标</param>
        public Atom(Char atom, String atomtype, int atomserial, Char altLoc, Double Xlaber, Double Ylaber, Double Zlaber)
        {
            this.atom = atom;
            this.atomtype = atomtype;
            this.atomserial = atomserial;
            this.altLoc = altLoc;
            this.Xlaber = Xlaber;
            this.Ylaber = Ylaber;
            this.Zlaber = Zlaber;
        }

        /// <summary>
        /// 初始化原子
        /// </summary>
        /// <param name="atom">原子</param>
        /// <param name="atomtype">原子类型</param>
        /// <param name="atomserial">原子序号</param>
        /// <param name="altLoc"></param>
        /// <param name="Xlaber">X坐标</param>
        /// <param name="Ylaber">Y坐标</param>
        /// <param name="Zlaber">Z坐标</param>
        /// <param name="occupancy">occupancy</param>
        /// <param name="tempfactor">beta factor</param>
        public Atom(Char atom, String atomtype, int atomserial, Char altLoc, Double Xlaber, Double Ylaber, Double Zlaber, Double occupancy, Double tempfactor)
        {
            this.atom = atom;
            this.atomtype = atomtype;
            this.atomserial = atomserial;
            this.altLoc = altLoc;
            this.Xlaber = Xlaber;
            this.Ylaber = Ylaber;
            this.Zlaber = Zlaber;
            this.occupancy = occupancy;
            this.tempfactor = tempfactor;
        }

        /// <summary>
        /// 初始化原子
        /// </summary>
        /// <param name="atom">原子</param>
        /// <param name="atomtype">原子类型</param>
        /// <param name="atomserial">原子序号</param>
        /// <param name="altLoc"></param>
        /// <param name="Xlaber">X坐标</param>
        /// <param name="Ylaber">Y坐标</param>
        /// <param name="Zlaber">Z坐标</param>
        /// <param name="occupancy">occupancy</param>
        /// <param name="tempfactor">beta factor</param>
        /// <param name="element">元素</param>
        /// <param name="charge">带电</param>
        public Atom(Char atom, String atomtype, int atomserial, Char altLoc, Double Xlaber, Double Ylaber, Double Zlaber, Double occupancy, Double tempfactor, String element, String charge)
        {
            this.atom = atom;
            this.atomtype = atomtype;
            this.atomserial = atomserial;
            this.altLoc = altLoc;
            this.Xlaber = Xlaber;
            this.Ylaber = Ylaber;
            this.Zlaber = Zlaber;
            this.occupancy = occupancy;
            this.tempfactor = tempfactor;
            this.element = element;
            this.charge = charge;
        }

        /// <summary>
        /// 返回原子
        /// </summary>
        /// <returns>原子</returns>
        public Char GetAtom()
        {
            return this.atom;
        }

        /// <summary>
        /// 返回原子类型
        /// </summary>
        /// <returns>原子类型</returns>
        public String GetAtomType()
        {
            return this.atomtype;
        }

        /// <summary>
        /// 返回原子序号
        /// </summary>
        /// <returns>原子序号</returns>
        public int GetAtomSerial()
        {
            return this.atomserial;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public Char GetAltloc()
        {
            return this.altLoc;
        }

        /// <summary>
        /// 返回X坐标
        /// </summary>
        /// <returns>X坐标</returns>
        public Double GetXlaber()
        {
            return this.Xlaber;
        }

        /// <summary>
        /// 返回Y坐标
        /// </summary>
        /// <returns>Y坐标</returns>
        public Double GetYLaber()
        {
            return this.Ylaber;
        }

        /// <summary>
        /// 返回Z坐标
        /// </summary>
        /// <returns>Z坐标</returns>
        public Double GetZlaber()
        {
            return this.Zlaber;
        }

        /// <summary>
        /// 返回occupancy
        /// </summary>
        /// <returns>occupancy</returns>
        public Double GetOccupancy()
        {
            return this.occupancy;
        }

        /// <summary>
        /// 返回beta factor
        /// </summary>
        /// <returns>beta factor</returns>
        public Double GetTempfactor()
        {
            return this.tempfactor;
        }

        /// <summary>
        /// 返回元素
        /// </summary>
        /// <returns>元素</returns>
        public String GetElement()
        {
            return this.element;
        }

        /// <summary>
        /// 返回带电
        /// </summary>
        /// <returns>charge</returns>
        public String GetCharge()
        {
            return this.charge;
        }

        /// <summary>
        /// 设置原子
        /// </summary>
        /// <param name="atom">原子</param>
        public void SetAtom(Char atom)
        {
            this.atom = atom;
        }

        /// <summary>
        /// 设置原子类型
        /// </summary>
        /// <param name="atomtype">原子类型</param>
        public void SetAtomType(String atomtype)
        {
            this.atomtype = atomtype;
        }

        /// <summary>
        /// 设置原子序号
        /// </summary>
        /// <param name="atomserial">原子序号</param>
        public void SetAtomSerial(int atomserial)
        {
            this.atomserial = atomserial;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="AltLoc"></param>
        public void SetAltloc(Char AltLoc)
        {
            this.altLoc = AltLoc;
        }

        /// <summary>
        /// 设置X坐标
        /// </summary>
        /// <param name="Xlaber">X坐标</param>
        public void SetXlaber(Double Xlaber)
        {
            this.Xlaber = Xlaber;
        }

        /// <summary>
        /// 设置Y坐标
        /// </summary>
        /// <param name="Ylaber">Y坐标</param>
        public void SetYlaber(Double Ylaber)
        {
            this.Ylaber = Ylaber;
        }

        /// <summary>
        /// 设置Z坐标
        /// </summary>
        /// <param name="Zlaber">Z坐标</param>
        public void SetZlaber(Double Zlaber)
        {
            this.Zlaber = Zlaber;
        }

        /// <summary>
        /// 设置occupancy
        /// </summary>
        /// <param name="Occupancy">occupancy</param>
        public void SetOccupancy(Double Occupancy)
        {
            this.occupancy = Occupancy;
        }

        /// <summary>
        /// 设置tempfactor，beta factor
        /// </summary>
        /// <param name="tempfactor">beta factor</param>
        public void SetTempfactor(Double tempfactor)
        {
            this.tempfactor = tempfactor;
        }

        /// <summary>
        /// 设置元素
        /// </summary>
        /// <param name="element">元素</param>
        public void SetElement(String element)
        {
            this.element = element;
        }

        /// <summary>
        /// 设置带电
        /// </summary>
        /// <param name="charge">带电</param>
        public void SetCharge(String charge)
        {
            this.charge = charge;
        }

        /// <summary>
        /// 将当前原子转化为3维向量
        /// </summary>
        /// <returns>3维向量</returns>
        public Vector3D ToVector3D()
        {
            return new Vector3D(this);
        }

        /// <summary>
        /// 将当前原子转化为3维点
        /// </summary>
        /// <returns>3维点</returns>
        public Point3D ToPoint3D()
        {
            return new Point3D(this);
        }

        /// <summary>
        /// 复制当前原子到新的原子
        /// </summary>
        /// <returns>新的原子</returns>
        public Atom Copy()
        {
            Atom atom = new Atom();
            atom.atom = this.atom;
            atom.atomserial = this.atomserial;
            atom.atomtype = this.atomtype;
            atom.altLoc = this.altLoc;
            atom.Xlaber = this.Xlaber;
            atom.Ylaber = this.Ylaber;
            atom.Zlaber = this.Zlaber;
            atom.occupancy = this.occupancy;
            atom.tempfactor = this.tempfactor;
            atom.element = this.element;
            atom.charge = this.charge;
            return atom;
        }

        /// <summary>
        /// 返回当前原子与目标原子间的距离
        /// </summary>
        /// <param name="atom">目标原子</param>
        /// <returns>距离</returns>
        public Double GetDistance(Atom atom)
        {
            Vector3D vc = new Vector3D(this, atom);
            return vc.Norm();
        }

        /// <summary>
        /// 以静态方式返回两个原子间的距离
        /// </summary>
        /// <param name="AtomOne">原子1</param>
        /// <param name="AtomTwo">原子2</param>
        /// <returns>距离</returns>
        public static Double AtomGetDistance(Atom AtomOne, Atom AtomTwo)
        {
            return AtomOne.GetDistance(AtomTwo);
        }

        /// <summary>
        /// 按平移向量vector平移当前原子到新的位置
        /// </summary>
        /// <param name="vector">平移向量</param>
        public void Move(Vector3D vector)
        {
            this.Xlaber += vector.X;
            this.Ylaber += vector.Y;
            this.Zlaber += vector.Z;
        }

        /// <summary>
        /// 按照平移向量平移原子到新的原子
        /// </summary>
        /// <param name="vector">平移向量</param>
        /// <returns>平移后的原子</returns>
        public Atom MoveToNewAtom(Vector3D vector)
        {
            Atom newAtom = this.Copy();
            newAtom.Xlaber += vector.X;
            newAtom.Ylaber += vector.Y;
            newAtom.Zlaber += vector.Z;
            return newAtom;
        }

        /// <summary>
        /// 按旋转矩阵RotationMatrix旋转当前原子到新的位置
        /// </summary>
        /// <param name="RotationMatrix">旋转矩阵</param>
        public void Rotate(Double[,] RotationMatrix)
        {
            if(RotationMatrix.GetLength(0)!=3&&RotationMatrix.GetLength(1)!=3)
            {
                Console.WriteLine("rotation matrix size is not illegal");
                Environment.Exit(1);
            }
            foreach(var node in RotationMatrix)
            {
                if(node<-1||node>1)
                {
                    Console.WriteLine("rotation matrix value is not illegal");
                    Environment.Exit(1);
                }
            }
            //        --------Rotation matrix to rotate  ------
            //m          u(m, 1)         u(m, 2)         u(m, 3)
            //1    - 0.4949315125   0.1654038329   0.8530441782
            //2    - 0.1440690798   0.9525078212 - 0.2682777492
            //3    - 0.8569054196 - 0.2556764020 - 0.4475965586
            //Code for rotating (x, y, z) to(X, Y, Z):
            //     X(i) = u(1, 1) * x(i) + u(1, 2) * y(i) + u(1, 3) * z(i)
            //     Y(i) = u(2, 1) * x(i) + u(2, 2) * y(i) + u(2, 3) * z(i)
            //     Z(i) = u(3, 1) * x(i) + u(3, 2) * y(i) + u(3, 3) * z(i)
            Double[,] u = RotationMatrix;
            Double newX = u[0, 0] * this.Xlaber + u[0, 1] * this.Ylaber + u[0, 2] * this.Zlaber;
            Double newY = u[1, 0] * this.Xlaber + u[1, 1] * this.Ylaber + u[1, 2] * this.Zlaber;
            Double newZ = u[2, 0] * this.Xlaber + u[2, 1] * this.Ylaber + u[2, 2] * this.Zlaber;
            this.Xlaber = newX;
            this.Ylaber = newY;
            this.Zlaber = newZ;
        }

        /// <summary>
        /// 按旋转矩阵RotationMatrix旋转当前原子到新的位置，返回到新的原子
        /// </summary>
        /// <param name="RotationMatrix">旋转矩阵</param>
        /// <returns>旋转后的原子</returns>
        public Atom RotateToNewAtom(Double[,] RotationMatrix)
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
            //        --------Rotation matrix to rotate  ------
            //m          u(m, 1)         u(m, 2)         u(m, 3)
            //1    - 0.4949315125   0.1654038329   0.8530441782
            //2    - 0.1440690798   0.9525078212 - 0.2682777492
            //3    - 0.8569054196 - 0.2556764020 - 0.4475965586
            //Code for rotating (x, y, z) to(X, Y, Z):
            //     X(i) = u(1, 1) * x(i) + u(1, 2) * y(i) + u(1, 3) * z(i)
            //     Y(i) = u(2, 1) * x(i) + u(2, 2) * y(i) + u(2, 3) * z(i)
            //     Z(i) = u(3, 1) * x(i) + u(3, 2) * y(i) + u(3, 3) * z(i)
            Double[,] u = RotationMatrix;
            Atom newAtom = this.Copy();
            newAtom.Xlaber = u[0, 0] * this.Xlaber + u[0, 1] * this.Ylaber + u[0, 2] * this.Zlaber;
            newAtom.Ylaber = u[1, 0] * this.Xlaber + u[1, 1] * this.Ylaber + u[1, 2] * this.Zlaber;
            newAtom.Zlaber = u[2, 0] * this.Xlaber + u[2, 1] * this.Ylaber + u[2, 2] * this.Zlaber;
            return newAtom;
        }

    }
}
