#!/bin/bash
mono=/nfs/amino-home/zhengwei/bin/Mono/mono-4.2.1/bin/mono
METHOME=/nfs/amino-home/panqp/protein_stru/repo/met_predictor/Met-Predictor_12302019
fastabasename=../test_fasta
fastaname=test_fasta
$METHOME/lib/depth-1.0/DEPTH -i $fastabasename\_pdb/$fastaname.pdb -o $fastabasename\_features/$fastaname -n 10 -survive 3 -keep $fastabasename\_features/$fastaname-sol