#!/bin/bash

# this script is used to randomly split data to two subsets

file=$1 # input plink data without suffix 

mkdir splitdata1
mkdir splitdata2

n1=`cat  ${file}.fam|awk '{print $6}'|grep -n 1|wc -l`
n2=`cat  ${file}.fam|awk '{print $6}'|grep -n 2|wc -l`

n1=`echo $(($n1/2))`
n2=`echo $(($n2/2))`

cat ${file}.fam | awk '$6 == 1 { print $0 }' > tmplist1
cat ${file}.fam | awk '$6 == 2 { print $0 }' > tmplist2

shuf -n ${n1} tmplist1 > tmpfam1
shuf -n ${n2} tmplist2 > tmpfam2

rm tmplist1 tmplist2

cat tmpfam2 >> tmpfam1
mv tmpfam1 tmplist
rm tmpfam2

plink --bfile ${file} --keep tmplist --noweb --make-bed --out splitdata1/SNPdata
plink --bfile ${file} --remove tmplist --noweb --make-bed --out splitdata2/SNPdata

rm tmpfam1



