#!/bin/bash

projectdir=$1 # project directory name, full path
marginal=$2  # hygeSSI: marginal=1; mhygeSSI: marginal=0
N1=$3 # sample permutation starting number
N2=$4 # sample permutation ending number
cd ${projectdir}

# compute SNP-SNP interaction networks
# BridGE/scripts/run_hygessi.sh
for R in `seq $N1 $N2`
do
     run_hygessi.sh --randRun=$R --marginal=$marginal --plinkCluster2=gwas_data_final.cluster2 --nWorker=8
done
