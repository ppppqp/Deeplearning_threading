using System;
using System.IO;
namespace BioNet
{
    class Program
    {
        static void Main(string[] args)
        {
            StreamReader sr = new StreamReader("../../protein_stru/testFiles/1a4z.pdb");
            String name = "name";
            Protein protein = new Protein(sr, name);
            Chain chainA = protein.GetChain('A');
            Chain result = chainA.GetLoneDepth("residue-residue", "global");
        }
    }
}
