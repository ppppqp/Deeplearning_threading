#!/bin/bash

# mode = "train" | "validation" | "test"
input="../datasets/tr3672.txt"
mode=$1
if [[ $mode = "train" ]]
then
  echo downloading train
  input="../datasets/tr3672.txt"
elif [[ $mode = "validation" ]]
then
  echo downloading validation
  input="../datasets/val918.txt"
else
  echo downloading test
  input="../datasets/ts1199.txt"
fi

count=0

while read -r line
do
  pdb=${line%?}
  # echo ${pdb%?}
  wget -P ./$mode https://files.rcsb.org/download/${pdb^^}.pdb
  if [ $? == 0 ]
  then 
    ((count++));
  fi
done < "$input"
echo total pdb files: $count
# 9021