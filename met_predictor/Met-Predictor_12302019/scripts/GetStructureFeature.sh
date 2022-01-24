#!/bin/tcsh -f

############################################# set your own folder ##################################################################
# set Met_predictor home directory to your own path 
setenv METHOME /nfs/amino-home/panqp/protein_stru/repo/met_predictor/Met-Predictor_12302019
set mypython = /nfs/amino-library/anaconda/bin/python
##################### soft already contained in package, if you want to use your own version, then change them ##################### 
# Where the HHsuite/HHpred programs have been installed
setenv HHLIB /nfs/amino-library/HHpred/HHsuite_lometsD

# The name of the PfamScan dir
set pfamscandir = /nfs/amino-home/zhengwei/Pfam/PfamScan

# pfamscan database dir
set pfamdatadir = /nfs/amino-home/zhengwei/Pfam/PfamScan/data/
set mono = /nfs/amino-home/zhengwei/bin/Mono/mono-4.2.1/bin/mono
# Where the hmmer pragrams have been installed
set hmmerdir = $METHOME/lib/hmmer
setenv PATH $hmmerdir/bin:$PATH
#setenv PERL5LIB $METHOME/lib/PfamScan:$METHOME/lib/perl5lib/lib/perl5
#setenv PERL5LIB $METHOME/lib/PfamScan:$METHOME/lib/perl5lib/lib64/perl5:$PERL5LIB
#setenv PERL5LIB $METHOME/lib/PfamScan:$METHOME/lib/perl5lib/share/perl5:$PERL5LIB

setenv PERL5LIB $METHOME/lib/PfamScan:$PERL5LIB

# Where the structure bin dir
set structurebindir = $METHOME/lib/Structure
###################################################### don't change ##################################################################


# run pfam
set fastabasename = $1:t:r
echo $fastabasename
set hhoutdir = ./
\cp $1 $hhoutdir$fastabasename.fasta
echo "run pfamscan..."

$pfamscandir/pfam_scan.pl -fasta $fastabasename.fasta -dir $pfamdatadir -as -outfile $fastabasename.psout -cpu 1

$mono $METHOME/bin/PfamScanSplit.exe ./$fastabasename.fasta ./$fastabasename.psout ./

mkdir -p $fastabasename\_fasta
mkdir -p $fastabasename\_pdb
mkdir -p $fastabasename\_features
foreach i (5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30)
mkdir -p $fastabasename\_features/$i
end

cp $fastabasename*.fasta $fastabasename\_fasta

################### run hhpred ##################
set fastafiles=`ls $fastabasename\_fasta`
foreach fastafile ( $fastafiles )

set fastaname=$fastafile:t:r
echo $fastaname
$HHLIB/scripts/hhpred/hhpmodeller_align.pl -i $fastabasename\_fasta/$fastaname.fasta -o $fastaname.pdb   >/dev/null
mv $fastaname.pdb $fastabasename\_pdb
cp $fastabasename\_fasta/$fastaname.fasta $fastabasename\_features/
cp $fastabasename\_pdb/$fastaname.pdb $fastabasename\_features/
echo "run dssp..."
$METHOME/lib/Structure/dssp-2.0.4-linux-amd64 -i $fastabasename\_pdb/$fastaname.pdb -o $fastabasename\_features/$fastaname.dssp
echo "run Naccess and Chops..."
cp $METHOME/lib/Structure/accall ./
cp $METHOME/lib/Structure/accall.input ./
cp $METHOME/lib/Structure/standard.data ./
cp $METHOME/lib/Structure/vdw.radii ./
cp $fastabasename\_pdb/$fastaname.pdb ./input.pdb
./accall <accall.input
mv input.asa $fastabasename\_features/$fastaname.asa
mv input.rsa $fastabasename\_features/$fastaname.rsa
rm input.pdb
$mono $METHOME/lib/Structure/chops.exe $fastabasename\_features/$fastaname.asa $fastabasename\_features/ single
$mono $METHOME/lib/Structure/kthCH.exe $fastabasename\_pdb/$fastaname.pdb $fastabasename\_features/
echo "run residue depth..."
echo $fastabasename
$METHOME/lib/depth-1.0/DEPTH -i $fastabasename\_pdb/$fastaname.pdb -o $fastabasename\_features/$fastaname -n 10 -survive 3 -keep $fastabasename\_features/$fastaname-sol

echo "run L1depth and HSE..."
   foreach R (5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30)
   $mono $METHOME/lib/Structure/L1depth.exe $fastabasename\_pdb/  $fastabasename\_features/$R -getldepth atom-atom local $R  >/dev/null
      foreach hsetype (HSEAU HSEBU HSEAD HSEBD)
      $mypython $METHOME/lib/Structure/calculateHSE.py -t $hsetype -r $R -o $fastabasename\_features/$R/$fastaname\_$hsetype$R.pdb $fastabasename\_pdb/$fastaname.pdb  >/dev/null
      end
   end

end
rm accall.input  standard.data  vdw.radii *.log accall


