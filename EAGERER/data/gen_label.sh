#!/bin/bash

input="../datasets/tr3672.txt"
mode="train"
count=0
input="../datasets/tr3672.txt"
mode=$1
if [[ $mode = "train" ]]
then
  echo generating train label
  input="../datasets/tr3672.txt"
  src="./train_clean/"
  dest="./train_label_raw/"
elif [[ $mode = "validation" ]]
then
  echo generating validation label
  input="../datasets/val918.txt"
  src="./validation_clean/"
  dest="./validation_label_raw/"
else
  echo generating test label
  input="../datasets/ts1199.txt"
  src="./test_clean/"
  dest="./test_label_raw/"
fi

while read -r line
do
  echo $src${line}
  /nfs/amino-home/zhengwei/dssp/dssp-2.0.4-linux-amd64 -i $src${line} -o $dest${line}
  if [ $? == 0 ]; then 
    ((count++));
    echo $count
  fi
done < "$input"
echo total pdb files: $count
