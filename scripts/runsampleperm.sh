#!/bin/bash

_usage() {
###### USAGE ######
cat <<EOF
$*
 Usage: runsampleperm.sh <[options]>
 DESCRIPTION
 This script is used to summarize results and list signficant BPM/WPM/PATH
 interactions with corresponding statsitics.

 INPUTS:
  --projectDir=PROJECTDIR
     Projectdir is the data directory where results will also be stored at.
     Default is the current directory.

  --model=MODEL
     Prefered interaction model which can be chosen from:
          RR recessive-recessive interaction
          DD dominant-dominant interaction
          RD recessive-dominant interaction
          combined combined all three interactions
          AA additive-additive interaction
     Default is combined.

  --minPath=MINPATH
     Minimum number of genes in the gene set.
     Default is 10.

  --snpPerms=SNPPERMS
     Number of SNP shuffling randomization.
     Default is 10000.

  --binaryNetwork=BINARYNETWORK
    binaryNetwork=1: run BridGE based on binarized network
    binaryNetwork=0: run BridGE based on weighted network
    Default is binaryNetwork=0

  --randRun=RANDRUN
     For random permutation run identification.
     randRun=0 for real run
     randRun=1,2,... for random runs
     This input is required.

  --alpha1=ALPHA1
     Marginal significance level for joint effect.
     Defalut is 0.05.

  --alpha2=ALPHA2
     Marginal significance level for individual effects.
     Default is 0.05.
    
EOF
}

ARGS=$@

# check inputs
options=$@
arguments=${options}

if [ -z "${arguments}" ]; then _usage; exit; fi

for argument in $options
do
     case $argument in
     --projectDir=*) projectDir=${argument/*=/""} ;;
     --model=*) model=${argument/*=/""} ;;
     --minPath=*) minPath=${argument/*=/""} ;;
     --snpPerms=*) snpPerms=${argument/*=/""} ;;
     --binaryNetwork=*) binaryNetwork=${argument/*=/""} ;;
     --randRun=*) randRun=${argument/*=/""} ;;
     --alpha1=*) alpha1=${argument/*=/""} ;;
     --alpha2=*) alpha2=${argument/*=/""} ;;
     esac
done

# set defaults
if [ -z "${projectDir}" ]; then projectDir=`pwd`; fi
if [ -z "${model}" ]; then model=combined; fi
if [ -z "${minPath}" ]; then minPath=10; fi
if [ -z "${snpPerms}" ]; then snpPerms=10000; fi
if [ -z "${bnaryNetwork}" ]; then binaryNetwork=0; fi
if [ -z "${randRun}" ]; then echo "Please specify randRun to proceed."; exit; fi
if [ -z "${alpha1}" ]; then alpha1=0.05; fi
if [ -z "${alpha2}" ]; then alpha2=0.05; fi

cd ${projectDir}

if [ -z "${ssmFile}" ]; then
     if [ "${model}" = "RR" ]; then
          ssmFile=ssM_hygeSSI_alpha1${alpha1}_alpha2${alpha2}_RR
     elif [ "${model}" = "DD" ]; then
          ssmFile=ssM_hygeSSI_alpha1${alpha1}_alpha2${alpha2}_DD
     elif [ "${model}" = "RD" ]; then
          ssmFile=ssM_hygeSSI_alpha1${alpha1}_alpha2${alpha2}_RD
     elif [ "${model}" = "combined" ]; then
          ssmFile=ssM_hygeSSI_alpha1${alpha1}_alpha2${alpha2}_combined
     elif [ "${model}" = "AA" ]; then
          ssmFile=ssM_lr_cassi_pv0.05
     fi
fi

nice matlab -nodisplay -nodesktop -nosplash -r "genstats('${ssmFile}_R${randRun}','BPMind.mat',${binaryNetwork},${snpPerms},${minPath});exit" </dev/null> /dev/null
