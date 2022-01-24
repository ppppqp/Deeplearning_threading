#!/bin/bash

mono=/nfs/amino-home/zhengwei/bin/Mono/mono-4.2.1/bin/mono
mcs=/nfs/amino-home/zhengwei/bin/Mono/mono-4.2.1/bin/mcs

$mcs -out:output.exe /nfs/amino-home/panqp/protein_stru/repo/met_predictor/PtmGetFeatures/PtmGetFeatures/Properties/Program.cs
echo compile completed
# $mono output.exe