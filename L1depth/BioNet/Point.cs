using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BioNet
{
    public class Point3D
    {
        //member
        public Double X;
        public Double Y;
        public Double Z;
        //function

        /// <summary>
        /// 初始化3维点
        /// </summary>
        public Point3D() { }

        /// <summary>
        /// 初始化3维点，设置x，y，z
        /// </summary>
        /// <param name="X"></param>
        /// <param name="Y"></param>
        /// <param name="Z"></param>
        public Point3D(Double X,Double Y,Double Z)
        {
            this.X = X;
            this.Y = Y;
            this.Z = Z;
        }

        /// <summary>
        /// 以原子Atom坐标初始化3维点
        /// </summary>
        /// <param name="atom">原子</param>
        public Point3D(Atom atom)
        {
            this.X = atom.Xlaber;
            this.Y = atom.Ylaber;
            this.Z = atom.Zlaber;
        }

        /// <summary>
        /// 用3维向量初始化3维点，向量端点坐标为点坐标
        /// </summary>
        /// <param name="vector"></param>
        public Point3D(Vector3D vector)
        {
            this.X = vector.X;
            this.Y = vector.Y;
            this.Z = vector.Z;
        }

        /// <summary>
        /// 用球面坐标初始化
        /// </summary>
        /// <param name="point">球面坐标点</param>
        public Point3D(PointSphere point)
        {
            this.X = point.R * Math.Sin(point.Theta) * Math.Cos(point.Phi);
            this.Y = point.R * Math.Sin(point.Theta) * Math.Sin(point.Phi);
            this.Z = point.R * Math.Cos(point.Theta);
        }

        /// <summary>
        /// 将欧式坐标转化为球面坐标
        /// </summary>
        /// <returns>球面坐标</returns>
        public PointSphere ToPointSphere()
        {
            PointSphere point = new PointSphere(this);
            return point;
        }

        /// <summary>
        /// 返回当前点到原点的距离
        /// </summary>
        /// <returns>距离</returns>
        public Double Norm()
        {
            return Math.Sqrt(this.X * this.X + this.Y * this.Y + this.Z * this.Z);
        }

        /// <summary>
        /// 返回当前点到point点的距离
        /// </summary>
        /// <param name="point">目标点</param>
        /// <returns>距离</returns>
        public Double GetDistance(Point3D point)
        {
            Vector3D vc = new Vector3D(point, this);
            return vc.Norm();
        }

        /// <summary>
        /// 以静态方式返回两点间距离
        /// </summary>
        /// <param name="pointOne">点1</param>
        /// <param name="pointTwo">点2</param>
        /// <returns>距离</returns>
        public static Double GetDistance(Point3D pointOne, Point3D pointTwo)
        {
            return pointOne.GetDistance(pointTwo);
        }

        /// <summary>
        /// 返回三个点形成的键角的余弦值[-1,1]
        /// </summary>
        /// <param name="one">边1上的点</param>
        /// <param name="common">顶角点</param>
        /// <param name="two">边2上的点</param>
        /// <returns>键角余弦值</returns>
        public static Double CosAngle(Point3D one,Point3D common ,Point3D two)
        {
            Vector3D vc1 = new Vector3D(one, common);
            Vector3D vc2 = new Vector3D(two, common);
            return Vector3D.VectorCosAngle(vc1, vc2);
        }

        /// <summary>
        /// 返回三个点形成的键角[0,pi]
        /// </summary>
        /// <param name="one">边1上的点</param>
        /// <param name="common">顶角点</param>
        /// <param name="two">边2上的点</param>
        /// <returns>键角</returns>
        public static Double Angle(Point3D one, Point3D common, Point3D two)
        {
            return Math.Acos(Point3D.CosAngle(one, common, two));
        }

        /// <summary>
        /// 返回4个点形成的二面角的余弦值[-1,1]
        /// </summary>
        /// <param name="innerP1">平面1上的点</param>
        /// <param name="commonP1">棱线上的点，靠近P1</param>
        /// <param name="commonP2">棱线上的点，靠近P2</param>
        /// <param name="innerP2">平面2上的点</param>
        /// <returns>二面角余弦值</returns>
        public static Double CosDiheDralAngle(Point3D innerP1,Point3D commonP1,Point3D commonP2, Point3D innerP2)
        {
            Vector3D innerVc1 = new Vector3D(innerP1, commonP1);
            Vector3D commonVc1 = new Vector3D(commonP2, commonP1);
            Vector3D commonVc2 = new Vector3D(commonP1, commonP2);
            Vector3D innerVc2 = new Vector3D(innerP2, commonP2);
            Vector3D directVc1 = Vector3D.VectorCrossProduct(commonVc1, innerVc1);
            Vector3D directVc2 = Vector3D.VectorCrossProduct(commonVc2, innerVc2);
            return -Vector3D.VectorCosAngle(directVc1, directVc2);
        }

        /// <summary>
        /// 返回4个点形成的二面角[0,pi]
        /// </summary>
        /// <param name="innerP1">平面1上的点</param>
        /// <param name="commonP1">棱线上的点，靠近P1</param>
        /// <param name="commonP2">棱线上的点，靠近P2</param>
        /// <param name="innerP2">平面2上的点</param>
        /// <returns>二面角</returns>
        public static Double DiheDralAngle(Point3D innerP1, Point3D commonP1, Point3D commonP2, Point3D innerP2)
        {
            Vector3D innerVc1 = new Vector3D(innerP1, commonP1);
            Vector3D commonVc1 = new Vector3D(commonP2, commonP1);
            Vector3D commonVc2 = new Vector3D(commonP1, commonP2);
            Vector3D innerVc2 = new Vector3D(innerP2, commonP2);
            Vector3D directVc1 = Vector3D.VectorCrossProduct(commonVc1, innerVc1);
            Vector3D directVc2 = Vector3D.VectorCrossProduct(commonVc2, innerVc2);
            Double lamda2 = innerVc2.Multiply(commonVc2) / Math.Pow(commonVc2.Norm(), 2);
            Vector3D VcinPlane2 = innerVc2.Subtract(lamda2 * commonVc2);
            if (VcinPlane2.Multiply(directVc1) >= 0)
            {
                return Math.Acos(Point3D.CosDiheDralAngle(innerP1, commonP1, commonP2, innerP2));
            }
            else
            {
                return -Math.Acos(Point3D.CosDiheDralAngle(innerP1, commonP1, commonP2, innerP2));
            }
        }

    }

    public class PointSphere
    {
        //member
        public Double R;// (0, inf);
        public Double Theta;// angle to z [0,pi];
        public Double Phi;// angle in XOY plane [0,2*pi);
        //function

        /// <summary>
        /// 初始化球坐标
        /// </summary>
        public PointSphere() { }

        /// <summary>
        /// 初始化球坐标，半径r，与z正半轴夹角theta，xoy平面角phi
        /// </summary>
        /// <param name="r">半径，r>0</param>
        /// <param name="theta">与z轴正半轴夹角，范围[0,pi]</param>
        /// <param name="phi">xoy平面角，范围[0,2*pi)</param>
        public PointSphere(Double r,Double theta,Double phi)
        {
            if(r<=0)
            {
                Console.WriteLine("r is minus!");
                Environment.Exit(1);
            }
            this.R = r;
            if (theta >= 2 * Math.PI || theta < 0)
            {
                int k = Convert.ToInt32(Math.Floor((theta / (2 * Math.PI))));
                theta -= 2 * k * Math.PI;
            }
            if(theta>Math.PI)
            {
                theta = 2 * Math.PI - theta;
            }
            this.Theta = theta;
            if (phi >= 2 * Math.PI || phi < 0)
            {
                int k = Convert.ToInt32(Math.Floor((phi / (2 * Math.PI))));
                phi -= 2 * k * Math.PI;
            }
            this.Phi = phi;
        }

        /// <summary>
        /// 以欧式坐标初始化为球面坐标，当x=0，y=0时设置phi=0
        /// </summary>
        /// <param name="point">欧式坐标点</param>
        public PointSphere(Point3D point)
        {
            this.R = point.Norm();
            this.Theta = Math.Acos(point.Z / this.R);
            if(point.X==0&&point.Y==0)
            {
                this.Phi = 0;
            }
            else if(point.X==0&&point.Y>0)
            {
                this.Phi = Math.PI / 2;
            }
            else if(point.X==0&&point.Y<0)
            {
                this.Phi = 3 * Math.PI / 2;
            }
            else if(point.X<0&&point.Y==0)
            {
                this.Phi = Math.PI;
            }
            else if(point.X>0&&point.Y==0)
            {
                this.Phi = 0;
            }
            else
            {
                Double phi = Math.Atan(Math.Abs(point.Y / point.X));
                if(point.X>0&&point.Y>0)
                {
                    this.Phi = phi;
                }
                if (point.X < 0 && point.Y > 0)
                {
                    this.Phi = Math.PI - phi;
                }
                if (point.X < 0 && point.Y < 0)
                {
                    this.Phi = Math.PI + phi;
                }
                if (point.X > 0 && point.Y < 0)
                {
                    this.Phi = 2 * Math.PI - phi;
                }
            }

        }

        /// <summary>
        /// 将球面坐标转化为3维欧式坐标
        /// </summary>
        /// <returns>返回坐标点</returns>
        public Point3D ToPoint3D()
        {
            Point3D point = new Point3D();
            point.X = this.R * Math.Sin(this.Theta) * Math.Cos(this.Phi);
            point.Y = this.R * Math.Sin(this.Theta) * Math.Sin(this.Phi);
            point.Z = this.R * Math.Cos(this.Theta);
            return point;
        }

        /// <summary>
        /// 按球坐标移动当前点到新的位置，R增加addR，theta增加addthetha，
        /// phi增加addphi，增加量可正可负，theta在0或pi附近有跳跃，越界后返回
        /// </summary>
        /// <param name="addR">R的增加量</param>
        /// <param name="addtheta">theta的增加量</param>
        /// <param name="addphi">phi的增加量</param>
        public void Move(Double addR,Double addtheta, Double addphi)
        {
            this.R += addR;
            this.Theta += addtheta;
            this.Phi += addphi;
            PointSphere newpoint = new PointSphere(this.R, this.Theta, this.Phi);
            this.R = newpoint.R;
            this.Theta = newpoint.Theta;
            this.Phi = newpoint.Phi;
        }

        /// <summary>
        /// 按球坐标移动当前点到新的位置，返回新的点，R增加addR，theta增加addthetha，
        /// phi增加addphi，增加量可正可负，theta在0或pi附近有跳跃，越界后返回
        /// </summary>
        /// <param name="addR">R的增加量</param>
        /// <param name="addtheta">theta的增加量</param>
        /// <param name="addphi">phi的增加量</param>
        /// <returns>移动后的点</returns>
        public PointSphere MoveToNewPoint(Double addR, Double addtheta, Double addphi)
        {
            this.R += addR;
            this.Theta += addtheta;
            this.Phi += addphi;
            PointSphere newpoint = new PointSphere(this.R, this.Theta, this.Phi);
            return newpoint;
        }

        /// <summary>
        /// 返回当前点与目标点间球面距离
        /// </summary>
        /// <param name="point">目标点</param>
        /// <returns>球面距离</returns>
        public Double GetSphereDistance(PointSphere point)
        {
            if (Math.Abs(this.R - point.R) > 1e-1)
            {
                Console.WriteLine("two points are not in the same sphere!");
                Console.ReadKey();
                Environment.Exit(1);
                
            }
            Double longitude1 = this.Phi;
            Double latitude1 = Math.PI / 2 - this.Theta;
            Double longitude2 = point.Phi;
            Double latitude2 = Math.PI / 2 - point.Theta;
            Double value = Math.Cos(longitude1 - longitude2) * Math.Cos(latitude1) * Math.Cos(latitude2) + Math.Sin(latitude1) * Math.Sin(latitude2);
            if(value<=-1)
            {
                value = -1;
            }
            if(value>=1)
            {
                value = 1;
            }
            Double angle = Math.Acos(value);
            Double spheredistance = this.R * angle;
            return spheredistance;
        }

        /// <summary>
        /// 以静态方式返回球面两点间球面距离
        /// </summary>
        /// <param name="point1">球面坐标点1</param>
        /// <param name="point2">球面坐标点2</param>
        /// <returns>球面距离</returns>
        public static Double GetSphereDistance(PointSphere point1,PointSphere point2)
        {
            return point1.GetSphereDistance(point2);
        }


    }
}
