#!/bin/bash

# this script is used to generate random pheno for plink
# without matching cases and controls

plinkfam=$1
phenofile=$2  #name for output phenofile

n=`wc -l ${plinkfam}|awk '{print $1}'`
cat ${plinkfam} | awk '{print $1"\t"$2}' > tmptmp1
shuf -n ${n} ${plinkfam} | awk '{print $6}' >tmptmp2

pr -m -t -s tmptmp1 tmptmp2 >${phenofile}

rm tmptmp1 tmptmp2
