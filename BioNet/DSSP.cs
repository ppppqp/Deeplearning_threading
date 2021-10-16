using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CommonFunction;
using System.IO;
//read dssp file and build a DSSP class, dssp version=2.0.4
namespace BioNet
{
    class DSSPResidue
    {
        public int residueserial;
        public Char residuetype;
        public Char dssp;
        public Char ss;
        public Double psi;
        public Double phi;

        public DSSPResidue() { }

        public DSSPResidue(int residueserial, Char residuetype, Char dssp, Double psi,Double phi)
        {
            this.residueserial = residueserial;
            this.residuetype = residuetype;
            this.dssp = dssp;
            this.ss = DSSPtoSS(dssp);
            this.psi = psi;
            this.phi = phi;
        }

        public Char DSSPtoSS(Char dssp)
        {
            Char structure = ' ';
            switch (dssp)
            {
                case 'H': { structure = 'H'; }; break;
                case 'G': { structure = 'H'; }; break;
                case 'I': { structure = 'C'; }; break;
                case 'B': { structure = 'E'; }; break;
                case 'E': { structure = 'E'; }; break;
                case 'T': { structure = 'C'; }; break;
                case 'S': { structure = 'C'; }; break;
                case ' ': { structure = 'C'; }; break;
                default: break;
            }
            return structure;
        }

        public DSSPResidue Copy()
        {
            DSSPResidue dsspresidue = new DSSPResidue();
            dsspresidue.residueserial = this.residueserial;
            dsspresidue.residuetype = this.residuetype;
            dsspresidue.dssp = this.dssp;
            dsspresidue.ss = this.ss;
            dsspresidue.phi = this.phi;
            dsspresidue.psi = this.psi;
            return dsspresidue;
        }
        
    }
    class DSSPChain
    {
        public Char chainID;
        public List<DSSPResidue> dsspresidues;

        public DSSPChain()
        {
            this.dsspresidues = new List<DSSPResidue>();
        }
        public DSSPChain(Char chainID)
        {
            this.chainID = chainID;
            this.dsspresidues = new List<DSSPResidue>();
        }

        public DSSPChain(char chainID,List<DSSPResidue> dsspresidues)
        {
            this.dsspresidues = new List<DSSPResidue>();
            this.chainID = chainID;
            foreach(DSSPResidue dsspresidue in dsspresidues.ToList())
            {
                this.dsspresidues.Add(dsspresidue.Copy());
            }
        }

        public DSSPChain Copy()
        {
            DSSPChain dsspchain = new DSSPChain();
            dsspchain.chainID = this.chainID;
            foreach(DSSPResidue dsspresidue in this.dsspresidues.ToList())
            {
                dsspchain.dsspresidues.Add(dsspresidue.Copy());
            }
            return dsspchain;
        }


    }
    class DSSP
    {
        public String proteinname;
        public List<DSSPChain> dsspchains;

        public DSSP()
        {
            this.dsspchains = new List<DSSPChain>();
        }

        public DSSP(String proteinname)
        {
            this.proteinname = proteinname;
            this.dsspchains = new List<DSSPChain>();
        }

        public DSSP(String proteinname, List<DSSPChain> dsspchains)
        {
            this.proteinname = proteinname;
            this.dsspchains = new List<DSSPChain>();
            foreach(DSSPChain dsspchain in dsspchains.ToList())
            {
                this.dsspchains.Add(dsspchain.Copy());
            }
        }

        public DSSP(StreamReader sr, String proteinname)
        {
            this.proteinname = proteinname;
            this.dsspchains = new List<DSSPChain>();
            //read header
            while(sr.Peek()>-1)
            {
                String line = sr.ReadLine();
                if(line.StartsWith("  #  RESIDUE"))
                {
                    break;
                }
            }
            //read dssp
            while(sr.Peek()>-1)
            {
                String line = sr.ReadLine().TrimEnd();
                if (line.Length > 0)
                {
                    String residueserial_str = line.Substring(5, 5).Trim();
                    if(residueserial_str.Length==0)
                    {
                        continue;
                    }
                    int residueserial = Convert.ToInt16(residueserial_str);
                    Char chainID = line.ElementAt(11);
                    Char residuetype = line.ElementAt(13);
                    Char dssp = line.ElementAt(16);
                    Double phi = Convert.ToDouble(line.Substring(103, 6));
                    Double psi = Convert.ToDouble(line.Substring(109, 6));
                    DSSPResidue dsspresidue = new DSSPResidue(residueserial, residuetype, dssp, psi, phi);
                    if (this.dsspchains.Count == 0 || this.dsspchains.Last().chainID != chainID)
                    {
                        DSSPChain dsspchain = new DSSPChain();
                        dsspchain.chainID = chainID;
                        dsspchain.dsspresidues.Add(dsspresidue);
                        this.dsspchains.Add(dsspchain);
                    }
                    else
                    {
                        this.dsspchains.Last().dsspresidues.Add(dsspresidue);
                    }
                }

            }
        }

        public DSSP Copy()
        {
            DSSP dssp = new DSSP();
            dssp.proteinname = this.proteinname;
            foreach(DSSPChain dsspchain in this.dsspchains.ToList())
            {
                dssp.dsspchains.Add(dsspchain.Copy());
            }
            return dssp;
        }
    }
}
