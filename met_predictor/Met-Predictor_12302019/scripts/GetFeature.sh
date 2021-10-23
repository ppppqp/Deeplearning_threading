#!/bin/tcsh -f

############################################# set your own folder ##################################################################
# set Met_predictor home directory to your own path 
setenv METHOME /nfs/amino-home/panqp/protein_stru/repo/met_predictor/Met-Predictor_12302019
set mypython = /mnt/home/wuyunqiq/bin/python/anaconda2/bin/python
##################### soft already contained in package, if you want to use your own version, then change them ##################### 
# Where the HHsuite programs have been installed
setenv HHLIB $METHOME/lib/hhsuite-2.0.16-linux-x86_64/

# The name of hhblits database uniprot20
set uniprot20 = /mnt/home/wuyunqiq/scratch/jlspzw/database/DeepMSADB/uniclust30_2017_04/uniclust30_2017_04 #uniprot20_2015_06

# The name of the BLAST data bank
set dbname = /mnt/home/wuyunqiq/scratch/jlspzw/database/BlastNR/nr

# Where the NCBI programs have been installed
set ncbidir = $METHOME/lib/blast-2.2.26/bin

# Where the PSIPRED V2 programs have been installed
set Psi_execdir = $METHOME/lib/psipred3.3/bin

# Where the PSIPRED V2 data files have been installed
set Psi_datadir = $METHOME/lib/psipred3.3/data

# Where the DISOPRED V2 programs have been installed
set Diso_execdir = $METHOME/lib/disopred/bin

# Where the DISOPRED V2 data files have been installed
set Diso_datadir = $METHOME/lib/disopred/data

# Where the spineX script files have been installed
set SpineXhome = $METHOME/lib/spineX

# Where the spider-hse script files have been installed
set spiderdir = $METHOME/lib/SPIDER2_local/misc

set hsedir = $METHOME/lib/HSE

###################################################### don't change ##################################################################
# run hhblits
set hhtarget = $1:t:r
set hhoutdir = ./
\cp $1 $hhoutdir$hhtarget.fasta
echo "run hhblits..."
$HHLIB/bin/hhblits -i $hhoutdir$hhtarget.fasta -d $uniprot20 -oa3m $hhoutdir$hhtarget.a3m -n 3 -maxfilt 500000 -diff inf -id 99 -cov 60 -cpu 1 > $hhoutdir$hhtarget.hhblog
egrep -v "^>" $hhoutdir$hhtarget.a3m | sed 's/[a-z]//g' > $hhoutdir$hhtarget.aln

set basename = $1:r
set rootname = $basename:t

# Generate a "unique" temporary filename root
set hostid = `hostid`
set tmproot = temp$$$hostid

\cp -f $1 $tmproot.fasta

#$1 /home/jlspzw/Psipred/test.fasta
#echo $tmproot
#temp31897007f0101 hostid
#echo $basename
#/home/jlspzw/Psipred/test
#echo $rootname
#test

echo "Running PSI-BLAST with sequence" $1 "..."

$ncbidir/blastpgp -b 0 -j 3 -h 0.001 -d $dbname -i $tmproot.fasta -a 4 -C $tmproot.chk -Q $rootname.mat >& $rootname.blast

if ($status != 0) then
    tail $rootname.blast
    echo "FATAL: Error whilst running blastpgp - script terminated!"
    exit $status
endif

echo "Predicting secondary structure..."

echo $tmproot.chk > $tmproot.pn
echo $tmproot.fasta > $tmproot.sn

$ncbidir/makemat -P $tmproot

if ($status != 0) then
    echo "FATAL: Error whilst running makemat - script terminated!"
    exit $status
endif

echo Pass1 ...

$Psi_execdir/psipred $tmproot.mtx $Psi_datadir/weights.dat $Psi_datadir/weights.dat2 $Psi_datadir/weights.dat3 > $rootname.ss

if ($status != 0) then
    echo "FATAL: Error whilst running psipred - script terminated!"
    exit $status
endif

echo Pass2 ...

$Psi_execdir/psipass2 $Psi_datadir/weights_p2.dat 1 1.0 1.0 $rootname.ss2 $rootname.ss > $rootname.horiz_p

if ($status != 0) then
    echo "FATAL: Error whilst running psipass2 - script terminated!"
    exit $status
endif

echo "Predicting Protein Disorder..."

$Diso_execdir/disopred $rootname $tmproot.mtx $Diso_datadir/ ./

echo "Get ASA etc from spineX..."

echo $rootname > $tmproot.id
mkdir $tmproot
\cp -f $rootname.mat $tmproot/$rootname.mat
$SpineXhome/spX.pl $tmproot.id $tmproot
\mv ./spXout/$rootname.spXout ./

echo "Get HSE from spider-hse..."
\cp $rootname.mat $rootname.pssm
$mypython $spiderdir/pred_pssm.py $rootname.pssm
$hsedir/run1.sh ./ $rootname.spd3 $mypython

# Remove temporary files

echo Cleaning up ...
\rm -f $tmproot.* error.log
\rm -f $tmproot/$rootname.*
\rmdir $tmproot


echo "Final output files:" $rootname.ss2 $rootname.horiz_p $rootname.diso $rootname.horiz_d $rootname.mat $rootname.spineXout $rootname.hsa2 $rootname.hsb2 $rootname.a3m
echo "Finished."
