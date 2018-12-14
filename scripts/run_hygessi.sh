#!/bin/bash
n1=$1
n2=$2
marginal=$3
cluster2file=$4

for R in `seq $n1 $n2`
do
          if [ -z "$cluster2file" ]; then
               nice matlab -nodisplay -nodesktop -nosplash -r "computessi('combined',$marginal,0.05,0.05,'',10,${R});exit" </dev/null> /dev/null
          else
               nice matlab -nodisplay -nodesktop -nosplash -r "computessi('combined',$marginal,0.05,0.05,'$cluster2file',10,${R});exit" </dev/null> /dev/null
          fi
done
