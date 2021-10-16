using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using CommonFunction;
using BioNet;

namespace CalForceFieldParameters
{
    class BondKey
    {
        public String ResType;
        public String bondAtom1;
        public String bondAtom2;
        public Double forceConstant;
    }
    class CalculateForceFieldParameters
    {
        #region bond key list
        public static BondKey[] keylist = {
                                            //ARG 10 bond 9 type, CZ-NH1=CZ-NH2
                                            new BondKey { ResType="ARG",bondAtom1="N",bondAtom2="CA",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="ARG",bondAtom1="CA",bondAtom2="C",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="ARG",bondAtom1="C",bondAtom2="O",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="ARG",bondAtom1="CA",bondAtom2="CB",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="ARG",bondAtom1="CB",bondAtom2="CG",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="ARG",bondAtom1="CG",bondAtom2="CD",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="ARG",bondAtom1="CD",bondAtom2="NE",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="ARG",bondAtom1="NE",bondAtom2="CZ",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="ARG",bondAtom1="CZ",bondAtom2="NH",forceConstant=1000.00},//side chain NH1, NH2
                                            //new BondKey { ResType="ARG",bondAtom1="C",bondAtom2="OXT"},//main chain C terminal
                                            //GLN 8 bond 8 type
                                            new BondKey { ResType="GLN",bondAtom1="N",bondAtom2="CA",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="GLN",bondAtom1="CA",bondAtom2="C",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="GLN",bondAtom1="C",bondAtom2="O",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="GLN",bondAtom1="CA",bondAtom2="CB",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="GLN",bondAtom1="CB",bondAtom2="CG",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="GLN",bondAtom1="CG",bondAtom2="CD",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="GLN",bondAtom1="CD",bondAtom2="OE1",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="GLN",bondAtom1="CD",bondAtom2="NE2",forceConstant=1000.00},//side chain
                                            //new BondKey { ResType="GLN",bondAtom1="C",bondAtom2="OXT"},//main chain C terminal
                                            //PHE 11 bond 8 type
                                            new BondKey { ResType="PHE",bondAtom1="N",bondAtom2="CA",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="PHE",bondAtom1="CA",bondAtom2="C",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="PHE",bondAtom1="C",bondAtom2="O",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="PHE",bondAtom1="CA",bondAtom2="CB",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="PHE",bondAtom1="CB",bondAtom2="CG",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="PHE",bondAtom1="CG",bondAtom2="CD",forceConstant=1000.00},//side chain CG-CD1 CG-CD2
                                            new BondKey { ResType="PHE",bondAtom1="CD",bondAtom2="CE",forceConstant=1000.00},//side chain CD1-CE1 CD2-CE2
                                            new BondKey { ResType="PHE",bondAtom1="CE",bondAtom2="CZ",forceConstant=1000.00},//side chain CE1-CZ CE2-CZ
                                            //new BondKey { ResType="PHE",bondAtom1="C",bondAtom2="OXT"},//main chain C terminal
                                            //TYR 12 bond 9 type
                                            new BondKey { ResType="TYR",bondAtom1="N",bondAtom2="CA",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="TYR",bondAtom1="CA",bondAtom2="C",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="TYR",bondAtom1="C",bondAtom2="O",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="TYR",bondAtom1="CA",bondAtom2="CB",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="TYR",bondAtom1="CB",bondAtom2="CG",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="TYR",bondAtom1="CG",bondAtom2="CD",forceConstant=1000.00},//side chain CG-CD1 CG-CD2
                                            new BondKey { ResType="TYR",bondAtom1="CD",bondAtom2="CE",forceConstant=1000.00},//side chain CD1-CE1 CD2-CE2
                                            new BondKey { ResType="TYR",bondAtom1="CE",bondAtom2="CZ",forceConstant=1000.00},//side chain CE1-CZ CE2-CZ
                                            new BondKey { ResType="TYR",bondAtom1="CZ",bondAtom2="OH",forceConstant=1000.00},//side chain
                                            //new BondKey { ResType="TYR",bondAtom1="C",bondAtom2="OXT"},//main chain C terminal
                                            //TRY 15 bond 15 type
                                            new BondKey { ResType="TRP",bondAtom1="N",bondAtom2="CA",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="TRP",bondAtom1="CA",bondAtom2="C",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="TRP",bondAtom1="C",bondAtom2="O",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="TRP",bondAtom1="CA",bondAtom2="CB",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="TRP",bondAtom1="CB",bondAtom2="CG",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="TRP",bondAtom1="CG",bondAtom2="CD1",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="TRP",bondAtom1="CG",bondAtom2="CD2",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="TRP",bondAtom1="CD1",bondAtom2="NE1",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="TRP",bondAtom1="NE1",bondAtom2="CE2",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="TRP",bondAtom1="CD2",bondAtom2="CE2",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="TRP",bondAtom1="CD2",bondAtom2="CE3",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="TRP",bondAtom1="CE2",bondAtom2="CZ2",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="TRP",bondAtom1="CE3",bondAtom2="CZ3",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="TRP",bondAtom1="CZ2",bondAtom2="CH2",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="TRP",bondAtom1="CZ3",bondAtom2="CH2",forceConstant=1000.00},//side chain
                                            //new BondKey { ResType="TRP",bondAtom1="C",bondAtom2="OXT"},//main chain C terminal
                                            //LYS 8 bond 8 type
                                            new BondKey { ResType="LYS",bondAtom1="N",bondAtom2="CA",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="LYS",bondAtom1="CA",bondAtom2="C",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="LYS",bondAtom1="C",bondAtom2="O",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="LYS",bondAtom1="CA",bondAtom2="CB",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="LYS",bondAtom1="CB",bondAtom2="CG",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="LYS",bondAtom1="CG",bondAtom2="CD",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="LYS",bondAtom1="CD",bondAtom2="CE",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="LYS",bondAtom1="CE",bondAtom2="NZ",forceConstant=1000.00},//side chain
                                            //new BondKey { ResType="LYS",bondAtom1="C",bondAtom2="OXT"},//main chain C terminal
                                            //GLY 3 bond 3 type
                                            new BondKey { ResType="GLY",bondAtom1="N",bondAtom2="CA",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="GLY",bondAtom1="CA",bondAtom2="C",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="GLY",bondAtom1="C",bondAtom2="O",forceConstant=1000.00},//main chain
                                            //new BondKey { ResType="GLY",bondAtom1="C",bondAtom2="OXT"},//main chain C terminal
                                            //ALA 4 bond 4 type
                                            new BondKey { ResType="ALA",bondAtom1="N",bondAtom2="CA",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="ALA",bondAtom1="CA",bondAtom2="C",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="ALA",bondAtom1="C",bondAtom2="O",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="ALA",bondAtom1="CA",bondAtom2="CB",forceConstant=1000.00},//side chain
                                            //new BondKey { ResType="ALA",bondAtom1="C",bondAtom2="OXT"},//main chain C terminal
                                            //HIS 10 bond 10 type
                                            new BondKey { ResType="HIS",bondAtom1="N",bondAtom2="CA",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="HIS",bondAtom1="CA",bondAtom2="C",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="HIS",bondAtom1="C",bondAtom2="O",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="HIS",bondAtom1="CA",bondAtom2="CB",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="HIS",bondAtom1="CB",bondAtom2="CG",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="HIS",bondAtom1="CG",bondAtom2="ND1",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="HIS",bondAtom1="CG",bondAtom2="CD2",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="HIS",bondAtom1="ND1",bondAtom2="CE1",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="HIS",bondAtom1="CD2",bondAtom2="NE2",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="HIS",bondAtom1="CE1",bondAtom2="NE2",forceConstant=1000.00},//side chain
                                            //new BondKey { ResType="HIS",bondAtom1="C",bondAtom2="OXT"},//main chain C terminal
                                            //SER 5 bond 5 type
                                            new BondKey { ResType="SER",bondAtom1="N",bondAtom2="CA",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="SER",bondAtom1="CA",bondAtom2="C",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="SER",bondAtom1="C",bondAtom2="O",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="SER",bondAtom1="CA",bondAtom2="CB",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="SER",bondAtom1="CB",bondAtom2="OG",forceConstant=1000.00},//side chain
                                            //new BondKey { ResType="SER",bondAtom1="C",bondAtom2="OXT"},//main chain C terminal
                                            //PRO 7 bond 7 type
                                            new BondKey { ResType="PRO",bondAtom1="N",bondAtom2="CA",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="PRO",bondAtom1="CA",bondAtom2="C",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="PRO",bondAtom1="C",bondAtom2="O",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="PRO",bondAtom1="CA",bondAtom2="CB",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="PRO",bondAtom1="CB",bondAtom2="CG",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="PRO",bondAtom1="CG",bondAtom2="CD",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="PRO",bondAtom1="CD",bondAtom2="N",forceConstant=1000.00},//side chain
                                            //new BondKey { ResType="PRO",bondAtom1="C",bondAtom2="OXT"},//main chain C terminal
                                            //GLU 8 bond 7 type
                                            new BondKey { ResType="GLU",bondAtom1="N",bondAtom2="CA",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="GLU",bondAtom1="CA",bondAtom2="C",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="GLU",bondAtom1="C",bondAtom2="O",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="GLU",bondAtom1="CA",bondAtom2="CB",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="GLU",bondAtom1="CB",bondAtom2="CG",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="GLU",bondAtom1="CG",bondAtom2="CD",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="GLU",bondAtom1="CD",bondAtom2="OE",forceConstant=1000.00},//side chain CD-OE1 CD-OE2
                                            //new BondKey { ResType="GLU",bondAtom1="C",bondAtom2="OXT"},//main chain C terminal
                                            //ASP 7 bond 6 type
                                            new BondKey { ResType="ASP",bondAtom1="N",bondAtom2="CA",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="ASP",bondAtom1="CA",bondAtom2="C",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="ASP",bondAtom1="C",bondAtom2="O",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="ASP",bondAtom1="CA",bondAtom2="CB",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="ASP",bondAtom1="CB",bondAtom2="CG",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="ASP",bondAtom1="CG",bondAtom2="OD",forceConstant=1000.00},//side chain CG-OD1 CG-OD2
                                            //new BondKey { ResType="ASP",bondAtom1="C",bondAtom2="OXT"},//main chain C terminal
                                            //THR 6 bond 6 type
                                            new BondKey { ResType="THR",bondAtom1="N",bondAtom2="CA",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="THR",bondAtom1="CA",bondAtom2="C",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="THR",bondAtom1="C",bondAtom2="O",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="THR",bondAtom1="CA",bondAtom2="CB",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="THR",bondAtom1="CB",bondAtom2="OG1",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="THR",bondAtom1="CB",bondAtom2="CG2",forceConstant=1000.00},//side chain
                                            //new BondKey { ResType="THR",bondAtom1="C",bondAtom2="OXT"},//main chain C terminal
                                            //CYS 5 bond 5 type
                                            new BondKey { ResType="CYS",bondAtom1="N",bondAtom2="CA",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="CYS",bondAtom1="CA",bondAtom2="C",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="CYS",bondAtom1="C",bondAtom2="O",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="CYS",bondAtom1="CA",bondAtom2="CB",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="CYS",bondAtom1="CB",bondAtom2="SG",forceConstant=1000.00},//side chain Disulfide bond will be calculated by other function
                                            //new BondKey { ResType="CYS",bondAtom1="C",bondAtom2="OXT"},//main chain C terminal
                                            //MET 7 bond 7 type
                                            new BondKey { ResType="MET",bondAtom1="N",bondAtom2="CA",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="MET",bondAtom1="CA",bondAtom2="C",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="MET",bondAtom1="C",bondAtom2="O",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="MET",bondAtom1="CA",bondAtom2="CB",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="MET",bondAtom1="CB",bondAtom2="CG",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="MET",bondAtom1="CG",bondAtom2="SD",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="MET",bondAtom1="SD",bondAtom2="CE",forceConstant=1000.00},//side chain
                                            //new BondKey { ResType="MET",bondAtom1="C",bondAtom2="OXT"},//main chain C terminal
                                            //LEU 7 bond 6 type
                                            new BondKey { ResType="LEU",bondAtom1="N",bondAtom2="CA",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="LEU",bondAtom1="CA",bondAtom2="C",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="LEU",bondAtom1="C",bondAtom2="O",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="LEU",bondAtom1="CA",bondAtom2="CB",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="LEU",bondAtom1="CB",bondAtom2="CG",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="LEU",bondAtom1="CG",bondAtom2="CD",forceConstant=1000.00},//side chain CG-CD1 CG-CD2
                                            //new BondKey { ResType="LEU",bondAtom1="C",bondAtom2="OXT"},//main chain C terminal
                                             //ASN 7 bond 7 type
                                            new BondKey { ResType="ASN",bondAtom1="N",bondAtom2="CA",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="ASN",bondAtom1="CA",bondAtom2="C",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="ASN",bondAtom1="C",bondAtom2="O",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="ASN",bondAtom1="CA",bondAtom2="CB",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="ASN",bondAtom1="CB",bondAtom2="CG",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="ASN",bondAtom1="CG",bondAtom2="OD1",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="ASN",bondAtom1="CG",bondAtom2="ND2",forceConstant=1000.00},//side chain
                                            //new BondKey { ResType="ASN",bondAtom1="C",bondAtom2="OXT"},//main chain C terminal
                                             //ILE 7 bond 7 type
                                            new BondKey { ResType="ILE",bondAtom1="N",bondAtom2="CA",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="ILE",bondAtom1="CA",bondAtom2="C",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="ILE",bondAtom1="C",bondAtom2="O",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="ILE",bondAtom1="CA",bondAtom2="CB",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="ILE",bondAtom1="CB",bondAtom2="CG2",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="ILE",bondAtom1="CB",bondAtom2="CG1",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="ILE",bondAtom1="CG1",bondAtom2="CD1",forceConstant=1000.00},//side chain
                                            //new BondKey { ResType="ILE",bondAtom1="C",bondAtom2="OXT"},//main chain C terminal
                                            //VAL 6 bond 5 type
                                            new BondKey { ResType="VAL",bondAtom1="N",bondAtom2="CA",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="VAL",bondAtom1="CA",bondAtom2="C",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="VAL",bondAtom1="C",bondAtom2="O",forceConstant=1000.00},//main chain
                                            new BondKey { ResType="VAL",bondAtom1="CA",bondAtom2="CB",forceConstant=1000.00},//side chain
                                            new BondKey { ResType="VAL",bondAtom1="CB",bondAtom2="CG",forceConstant=1000.00},//side chain CB-CG1 CB-CG2
                                            //new BondKey { ResType="VAL",bondAtom1="C",bondAtom2="OXT"},//main chain C terminal
                                          };
        #endregion
        public CalculateForceFieldParameters() { }

        public void CalculateBondLengthParameters(String datadir,String outdir)
        {
            //ARG
            Dictionary<BondKey, List<Double>> bonddatadic = new Dictionary<BondKey, List<Double>>();
            foreach(var keynode in keylist)
            {
                List<Double> bondlengthlist = new List<Double>();
                bonddatadic.Add(keynode, bondlengthlist);
            }
            String[] pdbfiles = Directory.GetFiles(datadir);
            List<Double> CNlength = new List<Double>();
            List<Double> COXTlength = new List<Double>();
            foreach(var pdbfile in pdbfiles)
            {
                StreamReader sr = new StreamReader(pdbfile);
                Console.WriteLine(pdbfile);
                Protein pro = new Protein(sr, "");
                sr.Close();
                for(int i=0;i<pro.chains.ElementAt(0).residues.Count;i++)
                {
                    Residue res = pro.chains.ElementAt(0).residues.ElementAt(i);

                    if(i< pro.chains.ElementAt(0).residues.Count-1)
                    {
                        Residue resC = res;
                        Residue resN = pro.chains.ElementAt(0).residues.ElementAt(i + 1);
                        if(resN.residueserial-resC.residueserial==1)
                        {
                            Atom CinResC = resC.atoms.FirstOrDefault(a => a.atomtype == "C");
                            Atom NinResN = resN.atoms.FirstOrDefault(a => a.atomtype == "N");
                            if (CinResC != null && NinResN != null)
                            {
                                Double dis = CinResC.GetDistance(NinResN);
                                //if (dis > 1.5)
                                //{
                                //    Console.WriteLine(pdbfile + " " + resC.residueserial + " " + dis);
                                //}
                                if (dis < 1.5)//else chain break?
                                {
                                    CNlength.Add(dis);
                                }
                                
                            }
                        }
                    }
                    if (i == pro.chains.ElementAt(0).residues.Count - 1)
                    {
                        Atom C = res.atoms.FirstOrDefault(a => a.atomtype == "C");
                        Atom OXT = res.atoms.FirstOrDefault(a => a.atomtype == "OXT");
                        if (C != null && OXT != null)
                        {
                            Double dis = C.GetDistance(OXT);
                            if (dis > 1.0 && dis < 1.5)
                            {
                                COXTlength.Add(dis);
                            }
                        }
                    }

                    #region add data to datadic
                    #region ARG
                    if (res.residuename=="ARG")
                    {
                        for (int j = 0; j < res.atoms.Count - 1; j++)
                        {
                            for (int k = j + 1; k < res.atoms.Count; k++)
                            {
                                Atom a1 = res.atoms.ElementAt(j);
                                Atom a2 = res.atoms.ElementAt(k);
                                if ((a1.atomtype=="CZ"&&a2.atomtype.StartsWith("NH"))|| (a2.atomtype == "CZ" && a1.atomtype.StartsWith("NH")))
                                {
                                    BondKey bk = keylist.FirstOrDefault(key => key.ResType == res.residuename && key.bondAtom1 == "CZ" && key.bondAtom2 == "NH");
                                    if (bk != null)
                                    {
                                        Double dis = a1.GetDistance(a2);
                                        if (dis < 2.0)
                                        {
                                            bonddatadic[bk].Add(dis);
                                        }
                                    }
                                }
                                else
                                {
                                    BondKey bk = keylist.FirstOrDefault(key => key.ResType == res.residuename && ((key.bondAtom1 == a1.atomtype && key.bondAtom2 == a2.atomtype) || (key.bondAtom1 == a2.atomtype && key.bondAtom2 == a1.atomtype)));
                                    if (bk != null)
                                    {
                                        Double dis = a1.GetDistance(a2);
                                        //if (a1.atomtype == "OXT" || a2.atomtype == "OXT")
                                        //{
                                        //    if (dis < 1.5 && dis > 1.0)
                                        //    {
                                        //        bonddatadic[bk].Add(dis);
                                        //    }
                                        //}
                                        //else
                                        {
                                            if (dis < 2.0)
                                            {
                                                bonddatadic[bk].Add(dis);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    #endregion
                    #region PHE TYR
                    else if (res.residuename=="PHE"|| res.residuename == "TYR")
                    {
                        for (int j = 0; j < res.atoms.Count - 1; j++)
                        {
                            for (int k = j + 1; k < res.atoms.Count; k++)
                            {
                                Atom a1 = res.atoms.ElementAt(j);
                                Atom a2 = res.atoms.ElementAt(k);
                                if ((a1.atomtype == "CG" && a2.atomtype.StartsWith("CD")) || (a2.atomtype == "CG" && a1.atomtype.StartsWith("CD")))
                                {
                                    BondKey bk = keylist.FirstOrDefault(key => key.ResType == res.residuename && key.bondAtom1 == "CG" && key.bondAtom2 == "CD");
                                    if (bk != null)
                                    {
                                        Double dis = a1.GetDistance(a2);
                                        if (dis < 2.0)
                                        {
                                            bonddatadic[bk].Add(dis);
                                        }
                                    }
                                }
                                else if((a1.atomtype == "CZ" && a2.atomtype.StartsWith("CE")) || (a2.atomtype == "CZ" && a1.atomtype.StartsWith("CE")))
                                {
                                    BondKey bk = keylist.FirstOrDefault(key => key.ResType == res.residuename && key.bondAtom1 == "CE" && key.bondAtom2 == "CZ");
                                    if (bk != null)
                                    {
                                        Double dis = a1.GetDistance(a2);
                                        if (dis < 2.0)
                                        {
                                            bonddatadic[bk].Add(dis);
                                        }
                                    }
                                }
                                else if((a1.atomtype=="CD1"&&a2.atomtype=="CE1")|| (a1.atomtype == "CD2" && a2.atomtype == "CE2"))
                                {
                                    BondKey bk = keylist.FirstOrDefault(key => key.ResType == res.residuename && key.bondAtom1 == "CD" && key.bondAtom2 == "CE");
                                    if (bk != null)
                                    {
                                        Double dis = a1.GetDistance(a2);
                                        if (dis < 2.0)
                                        {
                                            bonddatadic[bk].Add(dis);
                                        }
                                    }
                                }
                                else
                                {
                                    BondKey bk = keylist.FirstOrDefault(key => key.ResType == res.residuename && ((key.bondAtom1 == a1.atomtype && key.bondAtom2 == a2.atomtype) || (key.bondAtom1 == a2.atomtype && key.bondAtom2 == a1.atomtype)));
                                    if (bk != null)
                                    {
                                        Double dis = a1.GetDistance(a2);
                                        //if (a1.atomtype == "OXT" || a2.atomtype == "OXT")
                                        //{
                                        //    if (dis < 1.5 && dis > 1.0)
                                        //    {
                                        //        bonddatadic[bk].Add(dis);
                                        //    }
                                        //}
                                        //else
                                        {
                                            if (dis < 2.0)
                                            {
                                                bonddatadic[bk].Add(dis);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    #endregion
                    #region GLU
                    else if (res.residuename == "GLU")
                    {
                        for (int j = 0; j < res.atoms.Count - 1; j++)
                        {
                            for (int k = j + 1; k < res.atoms.Count; k++)
                            {
                                Atom a1 = res.atoms.ElementAt(j);
                                Atom a2 = res.atoms.ElementAt(k);
                                if ((a1.atomtype == "CD" && a2.atomtype.StartsWith("OE")) || (a2.atomtype == "CD" && a1.atomtype.StartsWith("OE")))
                                {
                                    BondKey bk = keylist.FirstOrDefault(key => key.ResType == res.residuename && key.bondAtom1 == "CD" && key.bondAtom2 == "OE");
                                    if (bk != null)
                                    {
                                        Double dis = a1.GetDistance(a2);
                                        if (dis < 2.0)
                                        {
                                            bonddatadic[bk].Add(dis);
                                        }
                                    }
                                }
                                else
                                {
                                    BondKey bk = keylist.FirstOrDefault(key => key.ResType == res.residuename && ((key.bondAtom1 == a1.atomtype && key.bondAtom2 == a2.atomtype) || (key.bondAtom1 == a2.atomtype && key.bondAtom2 == a1.atomtype)));
                                    if (bk != null)
                                    {
                                        Double dis = a1.GetDistance(a2);
                                        //if (a1.atomtype == "OXT" || a2.atomtype == "OXT")
                                        //{
                                        //    if (dis < 1.5 && dis > 1.0)
                                        //    {
                                        //        bonddatadic[bk].Add(dis);
                                        //    }
                                        //}
                                        //else
                                        {
                                            if (dis < 2.0)
                                            {
                                                bonddatadic[bk].Add(dis);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    #endregion
                    #region ASP
                    else if (res.residuename == "ASP")
                    {
                        for (int j = 0; j < res.atoms.Count - 1; j++)
                        {
                            for (int k = j + 1; k < res.atoms.Count; k++)
                            {
                                Atom a1 = res.atoms.ElementAt(j);
                                Atom a2 = res.atoms.ElementAt(k);
                                if ((a1.atomtype == "CG" && a2.atomtype.StartsWith("OD")) || (a2.atomtype == "CG" && a1.atomtype.StartsWith("OD")))
                                {
                                    BondKey bk = keylist.FirstOrDefault(key => key.ResType == res.residuename && key.bondAtom1 == "CG" && key.bondAtom2 == "OD");
                                    if (bk != null)
                                    {
                                        Double dis = a1.GetDistance(a2);
                                        if (dis < 2.0)
                                        {
                                            bonddatadic[bk].Add(dis);
                                        }
                                    }
                                }
                                else
                                {
                                    BondKey bk = keylist.FirstOrDefault(key => key.ResType == res.residuename && ((key.bondAtom1 == a1.atomtype && key.bondAtom2 == a2.atomtype) || (key.bondAtom1 == a2.atomtype && key.bondAtom2 == a1.atomtype)));
                                    if (bk != null)
                                    {
                                        Double dis = a1.GetDistance(a2);
                                        //if (a1.atomtype == "OXT" || a2.atomtype == "OXT")
                                        //{
                                        //    if (dis < 1.5 && dis > 1.0)
                                        //    {
                                        //        bonddatadic[bk].Add(dis);
                                        //    }
                                        //}
                                        //else
                                        {
                                            if (dis < 2.0)
                                            {
                                                bonddatadic[bk].Add(dis);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    #endregion
                    #region LEU
                    else if (res.residuename == "LEU")
                    {
                        for (int j = 0; j < res.atoms.Count - 1; j++)
                        {
                            for (int k = j + 1; k < res.atoms.Count; k++)
                            {
                                Atom a1 = res.atoms.ElementAt(j);
                                Atom a2 = res.atoms.ElementAt(k);
                                if ((a1.atomtype == "CG" && a2.atomtype.StartsWith("CD")) || (a2.atomtype == "CG" && a1.atomtype.StartsWith("CD")))
                                {
                                    BondKey bk = keylist.FirstOrDefault(key => key.ResType == res.residuename && key.bondAtom1 == "CG" && key.bondAtom2 == "CD");
                                    if (bk != null)
                                    {
                                        Double dis = a1.GetDistance(a2);
                                        if (dis < 2.0)
                                        {
                                            bonddatadic[bk].Add(dis);
                                        }
                                    }
                                }
                                else
                                {
                                    BondKey bk = keylist.FirstOrDefault(key => key.ResType == res.residuename && ((key.bondAtom1 == a1.atomtype && key.bondAtom2 == a2.atomtype) || (key.bondAtom1 == a2.atomtype && key.bondAtom2 == a1.atomtype)));
                                    if (bk != null)
                                    {
                                        Double dis = a1.GetDistance(a2);
                                        //if (a1.atomtype == "OXT" || a2.atomtype == "OXT")
                                        //{
                                        //    if (dis < 1.5 && dis > 1.0)
                                        //    {
                                        //        bonddatadic[bk].Add(dis);
                                        //    }
                                        //}
                                        //else
                                        {
                                            if (dis < 2.0)
                                            {
                                                bonddatadic[bk].Add(dis);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    #endregion
                    #region VAL
                    else if (res.residuename == "VAL")
                    {
                        for (int j = 0; j < res.atoms.Count - 1; j++)
                        {
                            for (int k = j + 1; k < res.atoms.Count; k++)
                            {
                                Atom a1 = res.atoms.ElementAt(j);
                                Atom a2 = res.atoms.ElementAt(k);
                                if ((a1.atomtype == "CB" && a2.atomtype.StartsWith("CG")) || (a2.atomtype == "CB" && a1.atomtype.StartsWith("CG")))
                                {
                                    BondKey bk = keylist.FirstOrDefault(key => key.ResType == res.residuename && key.bondAtom1 == "CB" && key.bondAtom2 == "CG");
                                    if (bk != null)
                                    {
                                        Double dis = a1.GetDistance(a2);                                   
                                        if (dis < 2.0)
                                        {
                                            bonddatadic[bk].Add(dis);
                                        }

                                    }
                                }
                                else
                                {
                                    BondKey bk = keylist.FirstOrDefault(key => key.ResType == res.residuename && ((key.bondAtom1 == a1.atomtype && key.bondAtom2 == a2.atomtype) || (key.bondAtom1 == a2.atomtype && key.bondAtom2 == a1.atomtype)));
                                    if (bk != null)
                                    {
                                        Double dis = a1.GetDistance(a2);
                                        //if (a1.atomtype == "OXT" || a2.atomtype == "OXT")
                                        //{
                                        //    if (dis < 1.5 && dis > 1.0)
                                        //    {
                                        //        bonddatadic[bk].Add(dis);
                                        //    }
                                        //}
                                        //else
                                        {
                                            if (dis < 2.0)
                                            {
                                                bonddatadic[bk].Add(dis);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    #endregion
                    #region ELSE
                    else
                    {
                        for(int j=0;j<res.atoms.Count-1;j++)
                        {
                            for(int k=j+1;k<res.atoms.Count;k++)
                            {
                                Atom a1 = res.atoms.ElementAt(j);
                                Atom a2 = res.atoms.ElementAt(k);
                                BondKey bk = keylist.FirstOrDefault(key => key.ResType == res.residuename && ((key.bondAtom1 == a1.atomtype && key.bondAtom2 == a2.atomtype) || (key.bondAtom1 == a2.atomtype && key.bondAtom2 == a1.atomtype)));
                                if (bk != null)
                                {
                                    Double dis = a1.GetDistance(a2);
                                    //if (a1.atomtype == "OXT" || a2.atomtype == "OXT")
                                    //{
                                    //    if (dis < 1.5 && dis > 1.0)
                                    //    {
                                    //        bonddatadic[bk].Add(dis);
                                    //    }
                                    //}
                                    //else
                                    {
                                        if (dis < 2.0)
                                        {
                                            bonddatadic[bk].Add(dis);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    #endregion
                    #endregion
                }
            }
            Console.WriteLine("finish collect data");
            StreamWriter bondlengthsw = new StreamWriter(@outdir+"/bondlength.dat");
            bondlengthsw.WriteLine("ResidueType BondAtom1 BondAtom2 AvgLength StdLength ForceConstant");
            //bondlengthsw.WriteLine("{0,-11} {1,-9} {2,-9} {3,-9:f5} {4,-9:f5} {5,-13:f2}", "ALL", "C", "N", Math.Round(CNlength.Average(), 5), Math.Round(Common.Standard(CNlength), 5),1000.00);
            //bondlengthsw.WriteLine("{0,-11} {1,-9} {2,-9} {3,-9:f5} {4,-9:f5} {5,-13:f2}", "ALL", "C", "OXT", Math.Round(COXTlength.Average(), 5), Math.Round(Common.Standard(COXTlength), 5),1000.00);
            //StreamWriter CNsw = new StreamWriter(@outdir + "/ALL_C-N.bondlength.tab");
            //Console.WriteLine("write file ALL C-N");
            //foreach (var node in CNlength)
            //{
            //    CNsw.WriteLine(node);
            //}
            //CNsw.Close();
            //StreamWriter COXTsw = new StreamWriter(@outdir+"/ALL_C-OXT.bondlength.tab");
            //Console.WriteLine("write file ALL C-OXT");
            //foreach (var node in COXTlength)
            //{
            //    COXTsw.WriteLine(node);
            //}
            //COXTsw.Close();
            BondKey allcn = new BondKey();
            allcn.ResType = "ALL";
            allcn.bondAtom1 = "C";
            allcn.bondAtom2 = "N";
            bonddatadic.Add(allcn, CNlength);
            BondKey allcoxt = new BondKey();
            allcoxt.ResType = "ALL";
            allcoxt.bondAtom1 = "C";
            allcoxt.bondAtom2 = "OXT";
            bonddatadic.Add(allcoxt, COXTlength);
            
            //String currentResType = "";
            //remove error data out range from [mu-10*sigma,mu+10*sigma];
            Dictionary<BondKey, List<Double>> bonddatadic2 = new Dictionary<BondKey, List<Double>>();
            foreach(var bondnode in bonddatadic)
            {
                Double mean = bondnode.Value.Average();
                Double std = Common.Standard(bondnode.Value);
                List<Double> data = new List<Double>();
                foreach(Double rawdata in bondnode.Value)
                {
                    if(rawdata<=mean+10*std&&rawdata>=mean-10*std)
                    {
                        data.Add(rawdata);
                    }
                }
                bondnode.Key.forceConstant = 1.0 / (2.0 * Common.Variance(data));
                bonddatadic2.Add(bondnode.Key, data);
            }
            foreach (var bonddata in bonddatadic2)
            {
                Console.WriteLine("write file " + bonddata.Key.ResType+" "+bonddata.Key.bondAtom1+"-"+bonddata.Key.bondAtom2);
                StreamWriter sw = new StreamWriter(@outdir+"/"+bonddata.Key.ResType+"_"+bonddata.Key.bondAtom1+"-"+bonddata.Key.bondAtom2+".bondlength.tab");
                foreach(var node in bonddata.Value)
                {
                    sw.WriteLine(node);
                }
                sw.Close();
                if(bonddata.Key.ResType=="ARG")
                {
                    if(bonddata.Key.bondAtom1=="CZ"&&bonddata.Key.bondAtom2=="NH")
                    {
                        bondlengthsw.WriteLine("{0,-11} {1,-9} {2,-9} {3,-9:f5} {4,-9:f5} {5,-13:f2}", bonddata.Key.ResType, bonddata.Key.bondAtom1, "NH1", Math.Round(bonddata.Value.Average(), 5), Math.Round(Common.Standard(bonddata.Value), 5), bonddata.Key.forceConstant);
                        bondlengthsw.WriteLine("{0,-11} {1,-9} {2,-9} {3,-9:f5} {4,-9:f5} {5,-13:f2}", bonddata.Key.ResType, bonddata.Key.bondAtom1, "NH2", Math.Round(bonddata.Value.Average(), 5), Math.Round(Common.Standard(bonddata.Value), 5), bonddata.Key.forceConstant);

                    }
                    else
                    {
                        bondlengthsw.WriteLine("{0,-11} {1,-9} {2,-9} {3,-9:f5} {4,-9:f5} {5,-13:f2}", bonddata.Key.ResType, bonddata.Key.bondAtom1, bonddata.Key.bondAtom2, Math.Round(bonddata.Value.Average(), 5), Math.Round(Common.Standard(bonddata.Value), 5), bonddata.Key.forceConstant);

                    }
                }
                else if(bonddata.Key.ResType=="PHE"|| bonddata.Key.ResType == "TYR")
                {
                    if (bonddata.Key.bondAtom1 == "CG" && bonddata.Key.bondAtom2 == "CD")
                    {
                        bondlengthsw.WriteLine("{0,-11} {1,-9} {2,-9} {3,-9:f5} {4,-9:f5} {5,-13:f2}", bonddata.Key.ResType, bonddata.Key.bondAtom1, "CD1", Math.Round(bonddata.Value.Average(), 5), Math.Round(Common.Standard(bonddata.Value), 5), bonddata.Key.forceConstant);
                        bondlengthsw.WriteLine("{0,-11} {1,-9} {2,-9} {3,-9:f5} {4,-9:f5} {5,-13:f2}", bonddata.Key.ResType, bonddata.Key.bondAtom1, "CD2", Math.Round(bonddata.Value.Average(), 5), Math.Round(Common.Standard(bonddata.Value), 5), bonddata.Key.forceConstant);

                    }
                    else if (bonddata.Key.bondAtom1 == "CD" && bonddata.Key.bondAtom2 == "CE")
                    {
                        bondlengthsw.WriteLine("{0,-11} {1,-9} {2,-9} {3,-9:f5} {4,-9:f5} {5,-13:f2}", bonddata.Key.ResType, "CD1", "CE1", Math.Round(bonddata.Value.Average(), 5), Math.Round(Common.Standard(bonddata.Value), 5), bonddata.Key.forceConstant);
                        bondlengthsw.WriteLine("{0,-11} {1,-9} {2,-9} {3,-9:f5} {4,-9:f5} {5,-13:f2}", bonddata.Key.ResType, "CD2", "CE2", Math.Round(bonddata.Value.Average(), 5), Math.Round(Common.Standard(bonddata.Value), 5), bonddata.Key.forceConstant);
                    }
                    else if (bonddata.Key.bondAtom1 == "CE" && bonddata.Key.bondAtom2 == "CZ")
                    {
                        bondlengthsw.WriteLine("{0,-11} {1,-9} {2,-9} {3,-9:f5} {4,-9:f5} {5,-13:f2}", bonddata.Key.ResType, "CE1", bonddata.Key.bondAtom2, Math.Round(bonddata.Value.Average(), 5), Math.Round(Common.Standard(bonddata.Value), 5), bonddata.Key.forceConstant);
                        bondlengthsw.WriteLine("{0,-11} {1,-9} {2,-9} {3,-9:f5} {4,-9:f5} {5,-13:f2}", bonddata.Key.ResType, "CE2", bonddata.Key.bondAtom2, Math.Round(bonddata.Value.Average(), 5), Math.Round(Common.Standard(bonddata.Value), 5), bonddata.Key.forceConstant);
                    }
                    else
                    {
                        bondlengthsw.WriteLine("{0,-11} {1,-9} {2,-9} {3,-9:f5} {4,-9:f5} {5,-13:f2}", bonddata.Key.ResType, bonddata.Key.bondAtom1, bonddata.Key.bondAtom2, Math.Round(bonddata.Value.Average(), 5), Math.Round(Common.Standard(bonddata.Value), 5), bonddata.Key.forceConstant);

                    }
                }
                else if(bonddata.Key.ResType == "GLU")
                {
                    if (bonddata.Key.bondAtom1 == "CD" && bonddata.Key.bondAtom2 == "OE")
                    {
                        bondlengthsw.WriteLine("{0,-11} {1,-9} {2,-9} {3,-9:f5} {4,-9:f5} {5,-13:f2}", bonddata.Key.ResType, bonddata.Key.bondAtom1, "OE1", Math.Round(bonddata.Value.Average(), 5), Math.Round(Common.Standard(bonddata.Value), 5), bonddata.Key.forceConstant);
                        bondlengthsw.WriteLine("{0,-11} {1,-9} {2,-9} {3,-9:f5} {4,-9:f5} {5,-13:f2}", bonddata.Key.ResType, bonddata.Key.bondAtom1, "OE2", Math.Round(bonddata.Value.Average(), 5), Math.Round(Common.Standard(bonddata.Value), 5), bonddata.Key.forceConstant);

                    }
                    else
                    {
                        bondlengthsw.WriteLine("{0,-11} {1,-9} {2,-9} {3,-9:f5} {4,-9:f5} {5,-13:f2}", bonddata.Key.ResType, bonddata.Key.bondAtom1, bonddata.Key.bondAtom2, Math.Round(bonddata.Value.Average(), 5), Math.Round(Common.Standard(bonddata.Value), 5), bonddata.Key.forceConstant);

                    }
                }
                else if(bonddata.Key.ResType == "ASP")
                {
                    if (bonddata.Key.bondAtom1 == "CG" && bonddata.Key.bondAtom2 == "OD")
                    {
                        bondlengthsw.WriteLine("{0,-11} {1,-9} {2,-9} {3,-9:f5} {4,-9:f5} {5,-13:f2}", bonddata.Key.ResType, bonddata.Key.bondAtom1, "OD1", Math.Round(bonddata.Value.Average(), 5), Math.Round(Common.Standard(bonddata.Value), 5), bonddata.Key.forceConstant);
                        bondlengthsw.WriteLine("{0,-11} {1,-9} {2,-9} {3,-9:f5} {4,-9:f5} {5,-13:f2}", bonddata.Key.ResType, bonddata.Key.bondAtom1, "OD2", Math.Round(bonddata.Value.Average(), 5), Math.Round(Common.Standard(bonddata.Value), 5), bonddata.Key.forceConstant);

                    }
                    else
                    {
                        bondlengthsw.WriteLine("{0,-11} {1,-9} {2,-9} {3,-9:f5} {4,-9:f5} {5,-13:f2}", bonddata.Key.ResType, bonddata.Key.bondAtom1, bonddata.Key.bondAtom2, Math.Round(bonddata.Value.Average(), 5), Math.Round(Common.Standard(bonddata.Value), 5), bonddata.Key.forceConstant);

                    }
                }
                else if(bonddata.Key.ResType == "LEU")
                {
                    if (bonddata.Key.bondAtom1 == "CG" && bonddata.Key.bondAtom2 == "CD")
                    {
                        bondlengthsw.WriteLine("{0,-11} {1,-9} {2,-9} {3,-9:f5} {4,-9:f5} {5,-13:f2}", bonddata.Key.ResType, bonddata.Key.bondAtom1, "CD1", Math.Round(bonddata.Value.Average(), 5), Math.Round(Common.Standard(bonddata.Value), 5), bonddata.Key.forceConstant);
                        bondlengthsw.WriteLine("{0,-11} {1,-9} {2,-9} {3,-9:f5} {4,-9:f5} {5,-13:f2}", bonddata.Key.ResType, bonddata.Key.bondAtom1, "CD2", Math.Round(bonddata.Value.Average(), 5), Math.Round(Common.Standard(bonddata.Value), 5), bonddata.Key.forceConstant);

                    }
                    else
                    {
                        bondlengthsw.WriteLine("{0,-11} {1,-9} {2,-9} {3,-9:f5} {4,-9:f5} {5,-13:f2}", bonddata.Key.ResType, bonddata.Key.bondAtom1, bonddata.Key.bondAtom2, Math.Round(bonddata.Value.Average(), 5), Math.Round(Common.Standard(bonddata.Value), 5), bonddata.Key.forceConstant);

                    }
                }
                else if(bonddata.Key.ResType == "VAL")
                {
                    if (bonddata.Key.bondAtom1 == "CB" && bonddata.Key.bondAtom2 == "CG")
                    {
                        bondlengthsw.WriteLine("{0,-11} {1,-9} {2,-9} {3,-9:f5} {4,-9:f5} {5,-13:f2}", bonddata.Key.ResType, bonddata.Key.bondAtom1, "CG1", Math.Round(bonddata.Value.Average(), 5), Math.Round(Common.Standard(bonddata.Value), 5), bonddata.Key.forceConstant);
                        bondlengthsw.WriteLine("{0,-11} {1,-9} {2,-9} {3,-9:f5} {4,-9:f5} {5,-13:f2}", bonddata.Key.ResType, bonddata.Key.bondAtom1, "CG2", Math.Round(bonddata.Value.Average(), 5), Math.Round(Common.Standard(bonddata.Value), 5), bonddata.Key.forceConstant);

                    }
                    else
                    {
                        bondlengthsw.WriteLine("{0,-11} {1,-9} {2,-9} {3,-9:f5} {4,-9:f5} {5,-13:f2}", bonddata.Key.ResType, bonddata.Key.bondAtom1, bonddata.Key.bondAtom2, Math.Round(bonddata.Value.Average(), 5), Math.Round(Common.Standard(bonddata.Value), 5), bonddata.Key.forceConstant);

                    }
                }
                else
                {
                    bondlengthsw.WriteLine("{0,-11} {1,-9} {2,-9} {3,-9:f5} {4,-9:f5} {5,-13:f2}", bonddata.Key.ResType, bonddata.Key.bondAtom1, bonddata.Key.bondAtom2, Math.Round(bonddata.Value.Average(), 5), Math.Round(Common.Standard(bonddata.Value), 5), bonddata.Key.forceConstant);
                }
                //if (currentResType != bonddata.Key.ResType)
                //{
                //    bondlengthsw.WriteLine("{0,-11} {1,-9} {2,-9} {3,-9:f5} {4,-9:f5} {5,-13:f2}", bonddata.Key.ResType, "C", "N", Math.Round(CNlength.Average(), 5), Math.Round(Common.Standard(CNlength), 5),1000.00);
                //    bondlengthsw.WriteLine("{0,-11} {1,-9} {2,-9} {3,-9:f5} {4,-9:f5} {5,-13:f2}", bonddata.Key.ResType, "C", "OXT", Math.Round(COXTlength.Average(), 5), Math.Round(Common.Standard(COXTlength), 5),1000.00);
                //    currentResType = bonddata.Key.ResType;
                //}

            }
            bondlengthsw.Close();

        }

        public void CalculateRamaParameters(String datadir, String dsspdir, String outdir)
        {
            String[] paths = Directory.GetFiles(@datadir);
            StreamWriter sw = new StreamWriter(@outdir+"/rama.dat");
            sw.WriteLine("phi psi");
            
            StreamWriter[,] AAtypeSSsw = new StreamWriter[20, 3];
            for(int i=0;i<20;i++)
            {
                for(int j=0;j<3;j++)
                {
                    String swpath = @outdir + "/" + Common.AAtype[i] + "_" + Common.SStype[j] + ".rama.dat";
                    AAtypeSSsw[i, j] = new StreamWriter(swpath);
                    AAtypeSSsw[i, j].WriteLine("phi psi");
                }
            }
            
            foreach (var path in paths)
            {
                String pname = Path.GetFileNameWithoutExtension(path);
                Console.WriteLine(pname);
                StreamReader sr = new StreamReader(path);
                Protein p = new Protein(sr, pname);
                sr.Close();
                StreamReader dsspsr = new StreamReader(@dsspdir + "/" + pname + ".dssp");
                DSSP dssp = new DSSP(dsspsr, pname);

                p.chains.ElementAt(0).CalculateRamaPsiPhi();

                for (int i = 1; i < p.chains.ElementAt(0).residues.Count - 1; i++)
                {
                    Residue r = p.chains.ElementAt(0).residues.ElementAt(i);
                    DSSPResidue dsspr = dssp.dsspchains.ElementAt(0).dsspresidues.FirstOrDefault(dr => dr.residueserial == r.residueserial);
                    if (r.psi != Common.impossibleMaxDouble && r.phi != Common.impossibleMaxDouble && dsspr != null)
                    {
                        sw.WriteLine("{0} {1}", r.phi * 180 / Math.PI, r.psi * 180 / Math.PI);
                        int aaindex = Common.AAtype.ToList().IndexOf(r.residuename);
                        int ssindex = Common.SStype.ToList().IndexOf(dsspr.ss);
                        if(aaindex>-1&&aaindex<20&&ssindex>-1)
                        {
                            if (r.phi <= Math.PI && r.phi >= -Math.PI && r.psi <= Math.PI && r.psi >= -Math.PI)
                            {
                                AAtypeSSsw[aaindex, ssindex].WriteLine("{0} {1}", r.phi * 180 / Math.PI, r.psi * 180 / Math.PI);
                            }
                        }
                    }
                }
            }
            sw.Close();
            for (int i = 0; i < 20; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    AAtypeSSsw[i, j].Close();
                }
            }
        }



    }
}
