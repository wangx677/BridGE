#!/bin/bash

# this script is used to handle duplicated SNPs in plink file

plinkfile=$1 # original plink file
output=$2 # new plink file after removing duplicated SNPs

cat ${plinkfile}.bim | awk '{print $2}'|sort|uniq -c | awk '$1>1'| awk '{print $2}' > ${output}
