using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CommonFunction;
using System.IO;
using BioNet;

namespace MonteCarloSimulation
{
    
    public class Energy
    {
        
        public Double ProteinTotalEnergy;
        public const Double ProteinBondWeight = 1.0;
        public Double ProteinBondEnergy;
        public const Double ProteinRamaWeight = 1.0;
        public Double ProteinRamaEnergy;

        public Energy() { }

        public Energy(Protein protein,List<Char> SecondStructure)
        {
            BondEnergy(protein);
            RamaEnergy(protein, SecondStructure);
            this.ProteinTotalEnergy = ProteinBondWeight * this.ProteinBondEnergy + ProteinRamaWeight * this.ProteinRamaEnergy;
        }

        /// <summary>
        /// 只计算第一条链的成键能量
        /// </summary>
        /// <param name="protein"></param>
        /// <returns></returns>
        public Double BondEnergy(Protein protein)
        {
            StreamWriter sw = new StreamWriter(@protein.proteinname + ".bond.eng");
            sw.WriteLine("resNum_i resNum_j resType BondAtom1 BondAtom2 distance bondEnergy totalEnergy");
            Double bondenergy = 0;
            for(int i=0;i<protein.chains.ElementAt(0).residues.Count;i++)
            {
                Residue resNode = protein.chains.ElementAt(0).residues.ElementAt(i);
                //bond in the residue; 
                for(int j=0;j<resNode.atoms.Count-1;j++)
                {
                    for(int k=j+1;k<resNode.atoms.Count;k++)
                    {
                        Atom a1 = resNode.atoms.ElementAt(j);
                        Atom a2 = resNode.atoms.ElementAt(k);
                        BondParameter bp = Common.BondLengthParameters.FirstOrDefault(blp => blp.ResidueType == resNode.residuename && (a1.atomtype == blp.BondAtom1 && a2.atomtype == blp.BondAtom2 || a1.atomtype == blp.BondAtom2 && a2.atomtype == blp.BondAtom1));
                        if(bp!=null)
                        {
                            Double dis = a1.GetDistance(a2);
                            Double eng=  bp.ForceConstant * Math.Pow((dis - bp.BondLengthAvg), 2);
                            bondenergy += eng;
                            sw.WriteLine("{0,-8} {1,-8} {2,-7} {3,-9} {4,-9} {5,-8:f5} {6,-10:f6} {7,-11:f6}", resNode.residueserial, resNode.residueserial, resNode.residuename, bp.BondAtom1, bp.BondAtom2, dis, eng, bondenergy);
                        }
                    }
                }
                //i i+1 C-N bond
                if(i< protein.chains.ElementAt(0).residues.Count-1)
                {
                    Residue nextResNode = protein.chains.ElementAt(0).residues.ElementAt(i + 1);
                    if (nextResNode.residueserial - resNode.residueserial == 1)
                    {
                        Atom Catom = resNode.atoms.FirstOrDefault(a => a.atomtype == "C");
                        Atom nextNatom = nextResNode.atoms.FirstOrDefault(a => a.atomtype == "N");
                        if (Catom != null && nextNatom != null)
                        {
                            BondParameter bp = Common.BondLengthParameters.FirstOrDefault(blp => blp.ResidueType == "ALL" && blp.BondAtom1 == "C" && blp.BondAtom2 == "N");
                            Double dis = Catom.GetDistance(nextNatom);
                            Double eng = 0.5 * bp.ForceConstant * Math.Pow((dis - bp.BondLengthAvg), 2);
                            bondenergy += eng;
                            sw.WriteLine("{0,-8} {1,-8} {2,-7} {3,-9} {4,-9} {5,-8:f5} {6,-10:f6} {7,-11:f6}", resNode.residueserial, nextResNode.residueserial, resNode.residuename, bp.BondAtom1, bp.BondAtom2, dis, eng, bondenergy);
                        }
                    }
                }
                //C terminal C-OXT bond
                if (i == protein.chains.ElementAt(0).residues.Count - 1)
                {
                    Atom Catom = resNode.atoms.FirstOrDefault(a => a.atomtype == "C");
                    Atom OXTatom = resNode.atoms.FirstOrDefault(a => a.atomtype == "OXT");
                    if(Catom!=null&&OXTatom!=null)
                    {
                        BondParameter bp = Common.BondLengthParameters.FirstOrDefault(blp => blp.ResidueType == "ALL" && blp.BondAtom1 == "C" && blp.BondAtom2 == "OXT");
                        Double dis = Catom.GetDistance(OXTatom);
                        Double eng = 0.5 * bp.ForceConstant * Math.Pow((dis - bp.BondLengthAvg), 2);
                        bondenergy += eng;
                        sw.WriteLine("{0,-8} {1,-8} {2,-7} {3,-9} {4,-9} {5,-8:f5} {6,-10:f6} {7,-11:f6}", resNode.residueserial, resNode.residueserial, resNode.residuename, bp.BondAtom1, bp.BondAtom2, dis, eng, bondenergy);
                    }
                }
            }
            this.ProteinBondEnergy = bondenergy;
            sw.Close();
            return bondenergy;
        }

        /// <summary>
        /// 只计算第一条链的Rama能量
        /// </summary>
        /// <param name="protein"></param>
        /// <param name="SecondStructure"></param>
        /// <returns></returns>
        public Double RamaEnergy(Protein protein,List<Char> SecondStructure)
        {
            Double infinityPsedoCount = 1e-40;
            Double ramaenergy = 0;
            Chain ChainOne = protein.chains.ElementAt(0);
            ChainOne.CalculateRamaPsiPhi();
            // residue 2 to n-1
            for(int i=1;i<ChainOne.residues.Count-1;i++)
            {
                Residue residuenode = ChainOne.residues.ElementAt(i);
                RamaParameter rp = Common.RamaParameters.FirstOrDefault(node => node.ResidueType == residuenode.residuename && node.SecondStructure == SecondStructure.ElementAt(i));
                Double phi = residuenode.phi * 180 / Math.PI;
                Double psi = residuenode.psi * 180 / Math.PI;
                //if (residuenode.phi >= -Math.PI && residuenode.phi <= Math.PI && residuenode.psi >= -Math.PI && residuenode.psi <= Math.PI)
                {
                    Double P = Common.KDE2calculate(rp.H, rp.PhiGrid, rp.PsiGrid, phi, psi);
                    if(P==0)
                    {
                        P += infinityPsedoCount;
                    }
                    ramaenergy += -Math.Log(P);
                }
            }
            this.ProteinRamaEnergy = ramaenergy;
            return ramaenergy;
        }

        public Double TestEnergy(Protein protein)
        {
            Double testenergy = 0;
            int count = protein.chains.ElementAt(0).residues.Count;// one chain 
            for (int i = 0; i < count - 1; i++)
            {
                for (int j = i + 1; j < count; j++)
                {
                    Atom a1 = protein.chains.ElementAt(0).residues.ElementAt(i).atoms.ElementAt(0);// one CA atom
                    Atom a2 = protein.chains.ElementAt(0).residues.ElementAt(j).atoms.ElementAt(0);//one CA atom
                    PointSphere p1 = new PointSphere(new Point3D(a1));
                    PointSphere p2 = new PointSphere(new Point3D(a2));
                    Double d = PointSphere.GetSphereDistance(p1, p2);
                    //Console.WriteLine(a1.atomserial + " " + a2.atomserial + " " + d);
                    Console.WriteLine(Math.Acos(Vector3D.VectorCosAngle(a1.ToVector3D(), a2.ToVector3D())) * 180 / Math.PI);
                    Double eng=1 / Math.Pow(d, 1);
                    testenergy += eng;
                }
            }

            return testenergy;
        }
    }
}
