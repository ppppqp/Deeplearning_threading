using System;
using System.IO;

namespace BioNet
{
    class Program
    {
        static void Main(string[] args)
        {
            StreamReader sr = new StreamReader("7n3oA.pdb");
            Protein protein = new Protein(sr, "7n3oA");
            Chain a = protein.GetChain(' ').GetLoneDepth("residue-residue", "Global");
        }
    }
}
