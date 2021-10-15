using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MonteCarloSimulation
{
    public class RandomNumber
    {
        public RandomNumber() { }

        /// <summary>
        /// 返回0到1之间的均匀分布的随机数
        /// </summary>
        /// <param name="len">长度</param>
        /// <returns></returns>
        public Double GetUnifromRandomNumber(int len)
        {
            string k = "0.";
            Random rand = new Random(GetRandomSeed());
            for (int i = 0; i < len; i++)
            {
                k += rand.Next(0, 10).ToString();
            }
            return Convert.ToDouble(k);
        }

        /// <summary>
        /// 取随机数种子
        /// </summary>
        /// <returns></returns>
        public int GetRandomSeed()
        {
            byte[] bytes = new byte[4];
            System.Security.Cryptography.RNGCryptoServiceProvider rng = new System.Security.Cryptography.RNGCryptoServiceProvider();
            rng.GetBytes(bytes);
            return BitConverter.ToInt32(bytes, 0);
        }
    }
}
