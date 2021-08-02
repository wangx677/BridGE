#!/bin/bash

projectdir=$1 # project directory name, full path
gi=$2  # genetic interaction type: hygeSSI, mhygeSSI
model=$3 # disease model: RR, DD, RD, combined

cd ${projectdir}


for i in `seq 0 1000`
do
     matlab -nodisplay -nodesktop -nosplash -r "getdensity('ssM_${gi}_alpha10.05_alpha20.05_${model}',$i);exit" </dev/null> /dev/null
done

matlab -nodisplay -nodesktop -nosplash -r "density2pvalue('ssM_${gi}_alpha10.05_alpha20.05_${model}',1000);exit" </dev/null> /dev/null
matlab -nodisplay -nodesktop -nosplash -r "computeFDR('ssM_${gi}_alpha10.05_alpha20.05_${model}',1000);exit" </dev/null> /dev/null

