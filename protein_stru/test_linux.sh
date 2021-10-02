#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
#SBATCH --job-name=test-job
#SBATCH --output=test_linux.out
module purge
echo "JOB START"
pwd
cd protein_stru/repo/protein_stru
make test_linux
./test_linux
echo "JOB END"