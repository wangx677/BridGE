#!/bin/bash
# This script is used to download gene annotation from Sanger's ftp site
# Only genes on chromosome 1-22 are included

# get release list from ftp site
curl  ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/ > tmplist

# find the newest release number
n=`cat tmplist | awk '{print $9}' | grep -v README | awk -F_ '{print $2}' | sort -nrk1,1 | head -1`

# download and extract file
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_${n}/gencode.v${n}.annotation.gtf.gz
gunzip gencode.v${n}.annotation.gtf.gz

# extract gene level information and remove genes not on chromosome 1-22
grep -w gene gencode.v${n}.annotation.gtf |grep -v chrX | grep -v chrY |grep -v chrM > tmp

# extract gene name, chromosome id, start/end position
cat tmp | awk '{print $1 " " $4 " " $5}' > tmp1
cat tmp|awk -F "gene_name" '{print $2}'|awk -F"\"" '{print $2}' >tmp2
paste -d " " tmp1 tmp2 > gencode.v${n}.annotation.csv
sed -i 's/chr//' gencode.v${n}.annotation.csv

# remove tmp files
rm gencode.v${n}.annotation.gtf
rm tmp tmp1 tmp2 tmplist
