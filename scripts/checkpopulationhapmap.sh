#!/bin/bash

# this script combines study data with hapmap data and check if the population 
# in the study data is consistent with the population in the hapmap data. 

PlinkFile=$1 # study data plink file 
inputpop=$2 # inputpop: study data's race, CEU-European ASW-African CHB-Chinese YRI-Youruba ${hapmappop}-Indian
pos=$3  #legend position on the MDS plot
hapmappop=$4  # HapMap population 
HapmapFile=$5  # HapMap plink file
popID=$6  # HapMap file population ids

# PlinkFile: plink file
# inputpop: study data's race, CEU-European ASW-African CHB-Chinese YRI-Youruba ${hapmappop}-Indian
# pos: legend position
# HapmapFile: HapMap plink file
# popID: HapMap file population ids


if [ -z ${pos} ]; then
	pos="topleft"
fi

if [ -z ${hapmappop} ]; then
	hapmappop=GIH
fi

if [ -z ${HapmapFile} ]; then
	HapmapFile=/project/csbio/HapMap_PhaseIII_2010/hapmap3_r3_b36_fwd.consensus.qc.poly
fi

if [ -z ${popID} ]; then
	popID=/project/csbio/HapMap_PhaseIII_2010/allpopid.txt
fi

# 1. merge with HapMap data
# only keep snps from study
cat ${PlinkFile}.bim | awk '{print $2}' > snplist.tmp
plink1.9 --bfile ${HapmapFile} --extract snplist.tmp --make-bed --out Hapmap_tmp > /dev/null

# update Hapmap data SNP information to build 37
cat ${PlinkFile}.bim | awk '{print $2 "\t" $1 "\t" $4}' > strand.tmp
plink_update_build_illumina.sh Hapmap_tmp strand.tmp Hapmap_tmp0
# only keep CEU, CHB, JPI and ASW populations

grep -w CEU ${popID} | awk '{print $1 "\t" $2}' > sublist.tmp
grep -w CHB ${popID} | awk '{print $1 "\t" $2}' >> sublist.tmp
grep -w ASW ${popID} | awk '{print $1 "\t" $2}' >> sublist.tmp
grep -w ${hapmappop} ${popID} | awk '{print $1 "\t" $2}' >> sublist.tmp
grep -w YRI ${popID} | awk '{print $1 "\t" $2}' >> sublist.tmp

# filter out nonfounders
plink1.9 --bfile Hapmap_tmp0 -filter-founders --make-bed --out Hapmap_tmp1 > /dev/null
plink1.9 --bfile Hapmap_tmp1 --keep sublist.tmp --make-bed --out Hapmap_tmp2 > /dev/null
plink1.9 --bfile Hapmap_tmp2 --allow-no-sex --bmerge ${PlinkFile} --make-bed --out allpop_tmp > /dev/null
while [ ! -f allpop_tmp.bim ] && [ -f allpop_tmp-merge.missnp ];
do
	 plink1.9 --bfile Hapmap_tmp2 --exclude allpop_tmp-merge.missnp --make-bed --out plink_tmp > /dev/null
         plink1.9 --bfile ${PlinkFile} --exclude allpop_tmp-merge.missnp --make-bed --out ${PlinkFile}  > /dev/null   	 
         rm allpop_tmp-merge.missnp
         plink1.9 --bfile plink_tmp --allow-no-sex --bmerge ${PlinkFile} --make-bed --out allpop_tmp > /dev/null

	#if [ -f allpop_tmp-merge.missnp ];then
	#	plink1.9 --bfile plink_tmp --exclude allpop_tmp-merge.missnp --make-bed --out plink_tmp > /dev/null
	#	plink1.9 --bfile ${PlinkFile} --exclude allpop_tmp-merge.missnp --make-bed --out ${PlinkFile}  > /dev/null
	#	rm allpop_tmp-merge.missnp
	#	plink1.9 --bfile plink_tmp --allow-no-sex --bmerge ${PlinkFile} --make-bed --out allpop_tmp > /dev/null
	#fi
done

plink1.9 --bfile allpop_tmp --geno 0.05 --make-bed --out allpop_tmp1 > /dev/null
qc.sh --PlinkFile=allpop_tmp1 --OutputFile=allpop_qc_tmp 

# ASW CEU CHB CHD ${hapmappop} JPT LWK MEX MKK TSI ASW
# genereate population idex
z=0
for i in `seq 1 5` ; do z=`echo ${z} 0`; done

echo FID SID CEU CHB ASW ${hapmappop} YRI StudyPop > allpopid.txt

while read line
do
        pattern=`echo ${line} | awk '{print $1 " " $2}'`
        pop=`grep -w "${pattern}" ${popID} | awk '{print $3}'`
        if [ "${pop}" == "CEU" ]; then
                echo ${pattern} `echo $z | awk '{ $1=1; print }'` >> allpopid.txt
        elif [ "${pop}" == "CHB" ]; then
                echo ${pattern} `echo $z | awk '{ $2=1; print }'` >> allpopid.txt
        elif [ "${pop}" == "ASW" ]; then
                echo ${pattern} `echo $z | awk '{ $3=1; print }'` >> allpopid.txt
        elif [ "${pop}" == "${hapmappop}" ]; then
                echo ${pattern} `echo $z | awk '{ $4=1;print }'` >> allpopid.txt
        elif [ "${pop}" == "YRI" ]; then
                echo ${pattern} `echo $z | awk '{ $5=1;print }'` >> allpopid.txt
        else
		echo ${pattern} `echo $z | awk -v x=6, '{ $x=1; print }'` >> allpopid.txt
        fi
done < allpop_qc_tmp.fam

# To perform principal components or MDS analysis, it is also important that we use SNPs that are not too correlated. 
plink1.9 --bfile allpop_qc_tmp --indep 50 5 2 --out prunedsnps_tmp > /dev/null
plink1.9 --bfile allpop_qc_tmp --extract prunedsnps_tmp.prune.in --make-bed --out prunedsnps_tmp > /dev/null

# Calculate genome-wide estimates of IBD sharing
plink1.9 --bfile prunedsnps_tmp --freq --allow-no-sex --out prunedsnps_tmp > /dev/null       
plink1.9 --bfile prunedsnps_tmp --read-freq prunedsnps_tmp.frq --genome --allow-no-sex --out prunedsnps_tmp > /dev/null

# Compute MDS scores
plink1.9 --noweb --bfile prunedsnps_tmp --read-genome prunedsnps_tmp.genome --cluster --mds-plot 2 --out prunedsnps_tmp > /dev/null

cp prunedsnps_tmp.mds ${PlinkFile}_HapMap.mds
plink_plot_mds_HapMap.r prunedsnps_tmp allpopid.txt MDSplot_${PlinkFile}_HapMap ${pos} ${inputpop} ${hapmappop}

# rm prunedsnps_tmp*
