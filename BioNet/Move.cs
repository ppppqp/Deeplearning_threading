using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BioNet;

namespace MonteCarloSimulation
{
    public class Move
    {
        public Move() { }

        public Protein Move1(Protein p)
        {
            Protein newp = new Protein(p.proteinname);
            Chain newchain = new Chain(p.chains.ElementAt(0).chainID);
            foreach (Residue residue in p.chains.ElementAt(0).residues)
            {
                Residue newresidue = new Residue(residue.residuename, residue.residueserial, residue.iCode);
                Atom newatom = residue.atoms.ElementAt(0).Copy();
                PointSphere ps = new PointSphere(new Point3D(residue.atoms.ElementAt(0)));
                RandomNumber rand = new RandomNumber();
                Double r1 = rand.GetUnifromRandomNumber(10);
                Double r2 = rand.GetUnifromRandomNumber(10);
                Double addtheta = (r1 * 2 - 1) * Math.PI / 50;
                Double addphi = (r2 * 2 - 1) * 2 * Math.PI / 50;
                Point3D newpoint = ps.MoveToNewPoint(0, addtheta, addphi).ToPoint3D();
                newatom.Xlaber = newpoint.X;
                newatom.Ylaber = newpoint.Y;
                newatom.Zlaber = newpoint.Z;
                newresidue.atoms.Add(newatom);
                newchain.residues.Add(newresidue);
            }
            newp.chains.Add(newchain);
            return newp;
        }
        public Protein Move2(Protein p)
        {
            Protein newp = new Protein(p.proteinname);
            Chain newchain = new Chain(p.chains.ElementAt(0).chainID);
            foreach (Residue residue in p.chains.ElementAt(0).residues)
            {
                Residue newresidue = new Residue(residue.residuename, residue.residueserial, residue.iCode);
                Atom newatom = residue.atoms.ElementAt(0).Copy();
                PointSphere ps = new PointSphere(new Point3D(residue.atoms.ElementAt(0)));
                RandomNumber rand = new RandomNumber();
                Double r1 = rand.GetUnifromRandomNumber(10);
                Double r2 = rand.GetUnifromRandomNumber(10);
                Double addtheta = (r1 * 2 - 1) * Math.PI / 360;
                Double addphi = (r2 * 2 - 1) * 2 * Math.PI / 360;
                Point3D newpoint = ps.MoveToNewPoint(0, addtheta, addphi).ToPoint3D();
                newatom.Xlaber = newpoint.X;
                newatom.Ylaber = newpoint.Y;
                newatom.Zlaber = newpoint.Z;
                newresidue.atoms.Add(newatom);
                newchain.residues.Add(newresidue);
            }
            newp.chains.Add(newchain);
            return newp;
        }

        public Protein Move3(Protein p)
        {
            Protein newp = new Protein(p.proteinname);
            Chain newchain = new Chain(p.chains.ElementAt(0).chainID);
            foreach (Residue residue in p.chains.ElementAt(0).residues)
            {
                Residue newresidue = new Residue(residue.residuename, residue.residueserial, residue.iCode);
                Atom newatom = residue.atoms.ElementAt(0).Copy();
                PointSphere ps = new PointSphere(new Point3D(residue.atoms.ElementAt(0)));
                RandomNumber rand = new RandomNumber();
                Double r1 = rand.GetUnifromRandomNumber(10);
                Double r2 = rand.GetUnifromRandomNumber(10);
                Double addtheta = (r1 * 2 - 1) * Math.PI / 720;
                Double addphi = (r2 * 2 - 1) * 2 * Math.PI / 720;
                Point3D newpoint = ps.MoveToNewPoint(0, addtheta, addphi).ToPoint3D();
                newatom.Xlaber = newpoint.X;
                newatom.Ylaber = newpoint.Y;
                newatom.Zlaber = newpoint.Z;
                newresidue.atoms.Add(newatom);
                newchain.residues.Add(newresidue);
            }
            newp.chains.Add(newchain);
            return newp;
        }
        public Protein Move4(Protein p)
        {
            Protein newp = new Protein(p.proteinname);
            Chain newchain = new Chain(p.chains.ElementAt(0).chainID);
            foreach (Residue residue in p.chains.ElementAt(0).residues)
            {
                Residue newresidue = new Residue(residue.residuename, residue.residueserial, residue.iCode);
                Atom newatom = residue.atoms.ElementAt(0).Copy();
                PointSphere ps = new PointSphere(new Point3D(residue.atoms.ElementAt(0)));
                RandomNumber rand = new RandomNumber();
                Double r1 = rand.GetUnifromRandomNumber(10);
                Double r2 = rand.GetUnifromRandomNumber(10);
                Double addtheta = (r1 * 2 - 1) * Math.PI / 900;
                Double addphi = (r2 * 2 - 1) * 2 * Math.PI / 1800;
                Point3D newpoint = ps.MoveToNewPoint(0, addtheta, addphi).ToPoint3D();
                newatom.Xlaber = newpoint.X;
                newatom.Ylaber = newpoint.Y;
                newatom.Zlaber = newpoint.Z;
                newresidue.atoms.Add(newatom);
                newchain.residues.Add(newresidue);
            }
            newp.chains.Add(newchain);
            return newp;
        }

        public Protein Move5(Protein p)
        {
            Protein newp = new Protein(p.proteinname);
            Chain newchain = new Chain(p.chains.ElementAt(0).chainID);
            foreach (Residue residue in p.chains.ElementAt(0).residues)
            {
                Residue newresidue = new Residue(residue.residuename, residue.residueserial, residue.iCode);
                Atom newatom = residue.atoms.ElementAt(0).Copy();
                PointSphere ps = new PointSphere(new Point3D(residue.atoms.ElementAt(0)));
                RandomNumber rand = new RandomNumber();
                Double r1 = rand.GetUnifromRandomNumber(10);
                Double r2 = rand.GetUnifromRandomNumber(10);
                Double addtheta = (r1 * 2 - 1) * Math.PI / 2000;
                Double addphi = (r2 * 2 - 1) * 2 * Math.PI / 4000;
                Point3D newpoint = ps.MoveToNewPoint(0, addtheta, addphi).ToPoint3D();
                newatom.Xlaber = newpoint.X;
                newatom.Ylaber = newpoint.Y;
                newatom.Zlaber = newpoint.Z;
                newresidue.atoms.Add(newatom);
                newchain.residues.Add(newresidue);
            }
            newp.chains.Add(newchain);
            return newp;
        }
    }
}
