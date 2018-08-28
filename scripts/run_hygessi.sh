#!/bin/bash
n1=$1
n2=$2

for R in `seq $n1 $n2`
do
          nice matlab -nodisplay -nodesktop -nosplash -r "computessi(4,0.05,0.05,'plinkFile.cluster2',10,${R});exit" </dev/null> /dev/null
done
