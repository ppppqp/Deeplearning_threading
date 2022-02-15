#!/bin/bash

mode=$1

if [[ $mode = "train" ]]
then
  echo cleaning train
  input="../datasets/tr3672.txt"
  src="./train/"
  dest="./train_clean/"
elif [[ $mode = "validation" ]]
then
  echo cleaning validation
  input="../datasets/val918.txt"
  src="./validation/"
  dest="./validation_clean/"
else
  echo cleaning test
  input="../datasets/ts1199.txt"
  src="./test/"
  dest="./test_clean/"
fi



# /nfs/amino-home/panqp/protein_stru/repo/protein_stru/test_pdb ../datasets/tr3672.txt ./train/ ./train_clean/
/nfs/amino-home/panqp/protein_stru/repo/protein_stru/test_pdb $input $src $dest
echo end