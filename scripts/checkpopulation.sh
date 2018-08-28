#!/bin/bash

# this script is used to make MDS plot based on given plink data

PlinkFile=$1 # plink file
pos=$2  # legend position

# To perform principal components or MDS analysis, it is also important that we use SNPs that are not too correlated. 
plink1.9 --bfile ${PlinkFile} --indep 50 5 2 --out prunedsnps_tmp > /dev/null
plink1.9 --bfile ${PlinkFile} --extract prunedsnps_tmp.prune.in --make-bed --out prunedsnps_tmp > /dev/null

if [ ! -f ${PlinkFile}.genome ];
then
	plink1.9 --bfile prunedsnps_tmp --freq --allow-no-sex --out prunedsnps_tmp > /dev/null       
	plink1.9 --bfile prunedsnps_tmp --read-freq prunedsnps_tmp.frq --genome --allow-no-sex --out prunedsnps_tmp > /dev/null
fi


plink1.9 --noweb --bfile prunedsnps_tmp --read-genome prunedsnps_tmp.genome --cluster --mds-plot 2 --out prunedsnps_tmp > /dev/null

cp prunedsnps_tmp.mds ${PlinkFile}.mds
plink_plot_mds.r prunedsnps_tmp MDSplot_${PlinkFile} ${pos}

rm prunedsnps_tmp*
