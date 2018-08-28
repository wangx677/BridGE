#!/bin/sh

#Required parameters:
#1. The original plink bed file bed (not including file extension suffix)
#2. The strand file to apply
#3. The new plink bed file name (not including file extension suffix)

#Unpack the parameters into labelled variables
stem=$1
strand_file=$2
outstem=$3

#Cut the strand file into a series of Plink slices
chr_file=tmp.chr
pos_file=tmp.pos
rs_file=tmp.rs
flip_file=tmp.flip
snp_file=tmp.snp

cat ${strand_file} | cut -f 1,3 |grep -v "\-\-\-" > ${chr_file}
cat ${strand_file} | cut -f 1,4 |grep -v "\-\-\-" > ${pos_file}
cat ${strand_file} | cut -f 1,2 |grep -v "\-\-\-" > ${rs_file}
cat ${strand_file} | cut -f 2 |grep -v "\-\-\-" > ${snp_file}
cat ${strand_file} | awk '{if ($5=="-") print $1}' | cut -f 1 > ${flip_file}

#Because Plink only allows you to update one attribute at a time, we need lots of temp
#Plink files
temp_prefix=TEMP_FILE_
temp1=${temp_prefix}"1"
temp2=${temp_prefix}"2"
temp3=${temp_prefix}"3"
temp4=${temp_prefix}"4"

#1. Apply the chr
plink1.9 --allow-no-sex --bfile ${stem} --update-chr ${chr_file} --make-bed --out ${temp1} > /dev/null
#2. Apply the pos
plink1.9 --allow-no-sex --bfile ${temp1} --update-map ${pos_file} --make-bed --out ${temp2} > /dev/null
#3. Apply the flip
plink1.9 --allow-no-sex --bfile ${temp2} --flip ${flip_file} --make-bed --out ${temp3} > /dev/null
#4. Apply the rs number
 plink1.9 --allow-no-sex --bfile ${temp3} --update-name ${rs_file} --make-bed --out ${temp4} > /dev/null
#5. Extract the SNPs in the snp file, we don't want SNPs that aren't in the strand file
plink1.9 --allow-no-sex --bfile ${temp4} --extract ${snp_file} --make-bed --out ${outstem} > /dev/null

#Now delete any temporary artefacts produced
rm -f ${temp_prefix}* tmp*
