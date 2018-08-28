#!/bin/bash

# this script is used to remove related individuals from the data

PlinkFile=$1 # plink file
OutputFile=$2 # new plink file
pi_hat=$3 # threshold for relatedness, default=0.2

if [ -z "${PlinkFile}" ]; then echo "Plink file is not provided"; exit; fi
if [ -z "${OutputFile}" ]; then echo "Output Plink file is not provided"; exit; fi	
if [ -z "${pi_hat}" ]; then pi_hat=0.2; fi

plink1.9 --bfile ${PlinkFile} --freq --allow-no-sex --out ${PlinkFile} > /dev/null       

plink1.9 --bfile ${PlinkFile} --read-freq ${PlinkFile}.frq --genome --allow-no-sex --out ${PlinkFile} > /dev/null

# find close related individuals and remove them

awk -v x=${pi_hat} '$10>x' ${PlinkFile}.genome |tail -n +2|awk '{print $1 "\t" $2}' |grep -v FID1 |sort |uniq > related_1.tmp
awk -v x=${pi_hat} '$10>x' ${PlinkFile}.genome |tail -n +2|awk '{print $3 "\t" $4}' |grep -v FID1 |sort |uniq > related_2.tmp


cat related_1.tmp > related_subject2remove.tmp
comm -12 related_1.tmp related_2.tmp > related_common.tmp

tmpfile2=plink_tmp

while [ -s related_common.tmp ]
do 
        plink1.9 --bfile ${PlinkFile} --keep related_2.tmp --make-bed --out ${tmpfile2} > /dev/null
        plink1.9 --bfile ${tmpfile2} --read-freq ${PlinkFile}.frq --genome --allow-no-sex --out ${tmpfile2} > /dev/null
        awk -v x=${pi_hat} '$10>x' ${tmpfile2}.genome |tail -n +2|awk '{print $1 "\t" $2}' |grep -v FID1 |sort |uniq > related_1.tmp
        awk -v x=${pi_hat} '$10>x' ${tmpfile2}.genome |tail -n +2|awk '{print $3 "\t" $4}' |grep -v FID1 |sort |uniq > related_2.tmp
        cat related_1.tmp >> related_subject2remove.tmp
        comm -12 related_1.tmp related_2.tmp > related_common.tmp
done
 
cat related_subject2remove.tmp |sort |uniq > related_subject2remove
plink1.9 --bfile ${PlinkFile} --remove related_subject2remove --allow-no-sex --make-bed --out ${OutputFile} > /dev/null

rm related_1.tmp related_2.tmp related_common.tmp   
