#!/bin/bash
prefix=../met_predictor/Met-Predictor_12302019
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
echo $fasta\n$aa\n$feature\n$sample\n$window\n$residue\n$outfile\n
python3 PtmGetFeatures.py $fasta $aa $feature $sample $window $residue $outfile