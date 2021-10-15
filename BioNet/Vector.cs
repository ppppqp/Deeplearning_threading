using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BioNet
{
    public class Vector3D
    {
        //member
        public Double X;
        public Double Y;
        public Double Z;
        //function

        /// <summary>
        /// 初始化3维向量
        /// </summary>
        public Vector3D() { }

        /// <summary>
        /// 初始化3维向量，坐标为x，y，z
        /// </summary>
        /// <param name="X"></param>
        /// <param name="Y"></param>
        /// <param name="Z"></param>
        public Vector3D(Double X,Double Y,Double Z)
        {
            this.X = X;
            this.Y = Y;
            this.Z = Z;
        }

        /// <summary>
        /// 以原子初始化3维向量，原子的x，y，z坐标为向量的3维坐标
        /// </summary>
        /// <param name="atom">初始化原子</param>
        public Vector3D(Atom atom)
        {
            this.X = atom.Xlaber;
            this.Y = atom.Ylaber;
            this.Z = atom.Zlaber;
        }

        /// <summary>
        /// 以3维点初始化3维向量
        /// </summary>
        /// <param name="point">Point3D点</param>
        public Vector3D(Point3D point)
        {
            this.X = point.X;
            this.Y = point.Y;
            this.Z = point.Z;
        }

        /// <summary>
        /// 以两点初始化3维向量，head-tail，tail->head
        /// </summary>
        /// <param name="head">向量头结点</param>
        /// <param name="tail">向量尾结点</param>
        public Vector3D(Point3D head,Point3D tail)
        {
            this.X = head.X - tail.X;
            this.Y = head.Y - tail.Y;
            this.Z = head.Z - tail.Z;
        }

        /// <summary>
        /// 以两个原子初始化3维向量，tail->head
        /// </summary>
        /// <param name="head">向量头结点</param>
        /// <param name="tail">向量尾结点</param>
        public Vector3D(Atom head, Atom tail)
        {
            this.X = head.Xlaber - tail.Xlaber;
            this.Y = head.Ylaber - tail.Ylaber;
            this.Z = head.Zlaber - tail.Zlaber;
        }

        /// <summary>
        /// 复制3维向量到新的向量
        /// </summary>
        /// <returns>复制后的向量</returns>
        public Vector3D Copy()
        {
            Vector3D vc = new Vector3D(this.X,this.Y,this.Z);
            return vc;
        }

        /// <summary>
        /// 返回三维向量的模
        /// </summary>
        /// <returns>模</returns>
        public Double Norm()
        {
            return Math.Sqrt(this.X * this.X + this.Y * this.Y + this.Z * this.Z);
        }

        /// <summary>
        /// 静态方法返回向量的模
        /// </summary>
        /// <param name="vector">输入向量</param>
        /// <returns>模</returns>
        public static Double VectorNorm(Vector3D vector)
        {
            return vector.Norm();
        }

        /// <summary>
        /// 返回两个向量的和，this+vector
        /// </summary>
        /// <param name="vector">被加向量</param>
        /// <returns>返回的和向量</returns>
        public Vector3D Plus(Vector3D vector)
        {
            Vector3D vc = new Vector3D();
            vc.X = this.X + vector.X;
            vc.Y = this.Y + vector.Y;
            vc.Z = this.Z + vector.Z;
            return vc;
        }

        /// <summary>
        /// 静态方式返回两个向量的和
        /// </summary>
        /// <param name="VectorOne">向量1</param>
        /// <param name="VectorTwo">向量2</param>
        /// <returns>返回的和向量</returns>
        public static Vector3D VectorPlus(Vector3D VectorOne, Vector3D VectorTwo)
        {
            //Vector3D VectorOut = new Vector3D();
            //VectorOut.X = VectorOne.X + VectorTwo.X;
            //VectorOut.Y = VectorOne.Y + VectorTwo.Y;
            //VectorOut.Z = VectorOne.Z + VectorTwo.Z;
            //return VectorOut;
            return VectorOne.Plus(VectorTwo);
        }

        /// <summary>
        /// 返回两个向量的差，this-vector
        /// </summary>
        /// <param name="vector">减向量</param>
        /// <returns>返回的差向量</returns>
        public Vector3D Subtract(Vector3D vector)
        {
            Vector3D vc = new Vector3D();
            vc.X = this.X - vector.X;
            vc.Y = this.Y - vector.Y;
            vc.Z = this.Z - vector.Z;
            return vc;
        }

        /// <summary>
        /// 以静态方式返回两个向量的差向量，vectorOne-vectorTwo
        /// </summary>
        /// <param name="VectorOne">被减向量</param>
        /// <param name="VectorTwo">减向量</param>
        /// <returns>返回的差向量</returns>
        public static Vector3D VectorSubtract(Vector3D VectorOne, Vector3D VectorTwo)
        {
            //Vector3D VectorOut = new Vector3D();
            //VectorOut.X = VectorOne.X - VectorTwo.X;
            //VectorOut.Y = VectorOne.Y - VectorTwo.Y;
            //VectorOut.Z = VectorOne.Z - VectorTwo.Z;
            //return VectorOut;
            return VectorOne.Subtract(VectorTwo);
        }

        /// <summary>
        /// 返回两个向量的内积，this*vector
        /// </summary>
        /// <param name="vector">向量</param>
        /// <returns>返回的内积</returns>
        public Double Multiply(Vector3D vector)
        {
            return this.X * vector.X + this.Y * vector.Y + this.Z * vector.Z;
        }

        /// <summary>
        /// 以静态方法返回两个向量的内积，vectorOne*vectorTwo
        /// </summary>
        /// <param name="VectorOne">向量1</param>
        /// <param name="VectorTwo">向量2</param>
        /// <returns>返回的内积</returns>
        public static Double VectorMultiply(Vector3D VectorOne, Vector3D VectorTwo)
        {
            //return VectorOne.X * VectorTwo.X + VectorOne.Y * VectorTwo.Y + VectorOne.Z * VectorTwo.Z;
            return VectorOne.Multiply(VectorTwo);
        }

        /// <summary>
        /// 返回两个向量的歪积（叉积），thisXvector，，返回向量方向服从右手系
        /// </summary>
        /// <param name="vector">乘积向量</param>
        /// <returns>返回的乘积向量</returns>
        public Vector3D CrossProduct(Vector3D vector)
        {
            //this X vector
            Vector3D vc = new Vector3D();
            vc.X = this.Y * vector.Z - this.Z * vector.Y;
            vc.Y = this.Z * vector.X - this.X * vector.Z;
            vc.Z = this.X * vector.Y - this.Y * vector.X;
            return vc;
        }

        /// <summary>
        /// 以静态方法返回两个向量的歪积（叉积），VectorOneXVectorTwo，返回向量方向服从右手系
        /// </summary>
        /// <param name="VectorOne">向量1</param>
        /// <param name="VectorTwo">向量2</param>
        /// <returns>返回的乘积向量</returns>
        public static Vector3D VectorCrossProduct(Vector3D VectorOne, Vector3D VectorTwo)
        {
            //vector1 X vector2
            //Vector3D VectorOut = new Vector3D();
            //VectorOut.X = VectorOne.Y * VectorTwo.Z - VectorOne.Z * VectorTwo.Y;
            //VectorOut.Y = VectorOne.Z * VectorTwo.X - VectorOne.X * VectorTwo.Z;
            //VectorOut.Z = VectorOne.X * VectorTwo.Y - VectorOne.Y * VectorTwo.X;
            //return VectorOut;
            return VectorOne.CrossProduct(VectorTwo);
        }

        /// <summary>
        /// 返回向量和向量vector的夹角的余弦值[-1,1]
        /// </summary>
        /// <param name="vector">向量</param>
        /// <returns>夹角的余弦值</returns>
        public Double CosAngle(Vector3D vector)
        {
            if(this.Norm()==0||vector.Norm()==0)
            {
                Console.WriteLine("one vector is zero!");
                Environment.Exit(0);
            }
            return this.Multiply(vector) / (this.Norm() * vector.Norm());
        }

        /// <summary>
        /// 以静态方式返回向量VectorOne和向量VectorTwo的夹角的余弦值[-1,1]
        /// </summary>
        /// <param name="VectorOne">向量1</param>
        /// <param name="VectorTwo">向量2</param>
        /// <returns>夹角的余弦值</returns>
        public static Double VectorCosAngle(Vector3D VectorOne, Vector3D VectorTwo)
        {
            //if (VectorOne.Norm() == 0 || VectorTwo.Norm() == 0)
            //{
            //    Console.WriteLine("one vector is zero!");
            //    Environment.Exit(0);
            //}
            //return Vector3D.VectorMultiply(VectorOne, VectorTwo) / (VectorOne.Norm() * VectorTwo.Norm());
            return VectorOne.CosAngle(VectorTwo);
        }

        /// <summary>
        /// 返回两个向量的夹角[0,pi]
        /// </summary>
        /// <param name="vector">向量</param>
        /// <returns>返回的夹角</returns>
        public Double Angle(Vector3D vector)
        {
            return Math.Acos(this.CosAngle(vector));
        }

        /// <summary>
        /// 以静态方式返回两个向量VectorOne，VectorTwo的夹角[0,pi]
        /// </summary>
        /// <param name="VectorOne">向量1</param>
        /// <param name="VectorTwo">向量2</param>
        /// <returns>返回的夹角</returns>
        public static Double VectorAngle(Vector3D VectorOne,Vector3D VectorTwo)
        {
            return VectorOne.Angle(VectorTwo);
        }

        /// <summary>
        /// 返回和当前向量同方向的单位向量
        /// </summary>
        /// <returns>返回的单位向量</returns>
        public Vector3D GetUnitVector()
        {
            Double Norm = this.Norm();
            if(Norm==0)
            {
                return this;
            }
            return new Vector3D(this.X / Norm, this.Y / Norm, this.Z / Norm);
        }

        public static Vector3D operator *(Double lambda, Vector3D vc)
        {
            Vector3D result = new Vector3D();
            result.X = vc.X * lambda;
            result.Y = vc.Y * lambda;
            result.Z = vc.Z * lambda;
            return result;
        }

        public static Vector3D operator *(Vector3D vc, Double lambda)
        {
            Vector3D result = new Vector3D();
            result.X = vc.X * lambda;
            result.Y = vc.Y * lambda;
            result.Z = vc.Z * lambda;
            return result;
        }

    }
}
