#!/bin/bash
n1=$1
n2=$2     

for i in `seq ${n1} ${n2}`
do
     nice matlab -nodisplay -nodesktop -nosplash -r "genstats('ssM_hygeSSI_alpha10.05_alpha20.05_RR_R${i}','','BPMind.mat',10,10);exit" </dev/null> /dev/null
done
