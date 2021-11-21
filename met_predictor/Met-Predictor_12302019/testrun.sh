#!/bin/bash
met_HOME=/nfs/amino-home/panqp/protein_stru/repo/met_predictor/Met-Predictor_12302019
PtmGetFeatures=$met_HOME/bin/PtmGetFeatures.exe;
PtmGetFeatures2=$met_HOME/bin/addstructurefeature.py;
mono=$met_HOME/lib/mono/bin/mono;
prefix=.
fasta=$prefix/example/P0CX53.fasta
aa=$prefix/data/AAindex.dat
feature=$prefix/example/example_features/P0CX53
sample=$prefix/data/sample/K.sample
window=17
residue=K
outfile=K.libsvm.old
# fasta=$prefix/test_fasta.fasta
# aa=$prefix/data/AAindex.dat
# feature=$prefix
# sample=$prefix/data/sample/K.sample
# window=5
# residue=K
# outfile=K.libsvm.old
$mono $PtmGetFeatures $fasta $aa $feature $sample $window $residue $outfile