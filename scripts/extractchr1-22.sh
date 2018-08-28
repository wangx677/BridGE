#!/bin/bash

# this script extracts SNPs on chromosome 1 to 22

inputfile=$1
outputfile=$2

chr=`cat ${inputfile}.bim | awk '{print $1}' | sort | uniq`
if [ "${chr}" == "0" ]; then
	plink1.9 --bfile ${inputfile} --allow-no-sex --make-bed --out ${outputfile} > /dev/null
else
	plink1.9 --bfile ${inputfile} --chr 1-22 --allow-no-sex --make-bed --out ${outputfile} > /dev/null
fi
