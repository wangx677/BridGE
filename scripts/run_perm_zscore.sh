#!/bin/bash
n1=$1
n2=$2
m=$3
k=$5

for i in `seq $n1 $n2`
do
     ssmFile=ssM_hygeSSI_alpha10.05_alpha20.05_${m}_R${i}
     nice matlab -nodisplay -nodesktop -nosplash -r "genstats_zscore('$ssmFile','BPMind.mat',0,${k});exit" </dev/null> /dev/null
done
