#!/bin/bash

# this script is used to remove related individuals from the data

PlinkFile=$1 # plink file
OutputFile=$2 # new plink file
pi_hat=$3 # threshold for relatedness, default=0.2


if [ -z "${PlinkFile}" ]; then echo "Plink file is not provided"; exit; fi
if [ -z "${OutputFile}" ]; then echo "Output Plink file is not provided"; exit; fi	
if [ -z "${pi_hat}" ]; then pi_hat=0.2; fi

# LD prune before calculate IBD
plink1.9 --bfile ${PlinkFile} --allow-no-sex --indep-pairwise 50 5 0.2
plink1.9 --bfile ${PlinkFile} --extract plink.prune.in --allow-no-sex --make-bed --out ${PlinkFile}_pruned_tmp

# calculate IBD
if [ ! -f ${PlinkFile}.genome ]; then
     plink1.9 --bfile ${PlinkFile}_pruned_tmp --freq --allow-no-sex --out ${PlinkFile} > /dev/null       
     plink1.9 --bfile ${PlinkFile}_pruned_tmp --read-freq ${PlinkFile}.frq --genome --min ${pi_hat} --allow-no-sex --out ${PlinkFile} > /dev/null
fi

# if *.genome file has more than 2 lines
if [[ $(wc -l <${PlinkFile}.genome) -ge 2 ]];then
     # find close related individuals and remove them
     if [ ! -f related_subject2remove.txt ]; then
          rm related_subject2remove.txt
     fi
     nice matlab -nodisplay -nodesktop -nosplash -r "removerelatedindividual('${PlinkFile}.genome',${pi_hat});exit" </dev/null> /dev/null
     plink1.9 --bfile ${PlinkFile} --remove related_subject2remove.txt --allow-no-sex --make-bed --out ${OutputFile} > /dev/null
else
     # no need to remove any individuals
     plink1.9 --bfile ${PlinkFile} --allow-no-sex --make-bed --out ${OutputFile} > /dev/null
fi
rm ${PlinkFile}_pruned_tmp*
