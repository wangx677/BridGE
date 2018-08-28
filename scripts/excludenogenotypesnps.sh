#!/bin/bash
# remove SNPs without assigned genotype

PlinkFile=$1
OutputFile=$2

awk '$5 != "0" && $6 != "0" { print $0 }' ${PlinkFile}.bim > tmp.snp
awk '$5 == "0" || $6 == "0" { print $0 }' ${PlinkFile}.bim > nogenotype_snp
plink1.9 --bfile ${PlinkFile} --extract tmp.snp --make-bed --out ${OutputFile} > /dev/null

rm tmp.snp nogenotype_snp
