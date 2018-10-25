#!/bin/bash

snpdata=$1

# To perform principal components or MDS analysis, it is also important that we use SNPs that are not too correlated. 
plink1.9 --bfile ${snpdata} --indep 50 5 2 --out ${snpdata}_tmp > /dev/null
plink1.9 --bfile ${snpdata} --extract ${snpdata}_tmp.prune.in --make-bed --out ${snpdata}_tmp > /dev/null

# Calculate genome-wide estimates of IBD sharing
plink1.9 --bfile ${snpdata}_tmp --freq --allow-no-sex --out ${snpdata}_tmp > /dev/null       
plink1.9 --bfile ${snpdata}_tmp --read-freq ${snpdata}_tmp.frq --genome --allow-no-sex --out ${snpdata}_tmp > /dev/null

# Compute MDS scores
plink1.9 --noweb --bfile ${snpdata}_tmp --read-genome ${snpdata}_tmp.genome --cluster --mds-plot 2 --out ${snpdata}_tmp > /dev/null

cp ${snpdata}_tmp.mds ${PlinkFile}_HapMap.mds
plink_plot_mds_HapMap.r ${snpdata}_tmp allpopid.txt MDSplot_${PlinkFile}_HapMap ${pos} ${inputpop} ${hapmappop}

