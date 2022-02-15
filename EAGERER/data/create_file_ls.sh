#!/bin/bash
input="../datasets/tr3672.txt"

while read -r line
do
  echo data/train_clean/$line 
  
done <"$input" > ../eg.list