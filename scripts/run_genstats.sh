#!/bin/bash

_usage() {
###### USAGE ######
cat <<EOF
$*
 Usage: run_genstats.sh <[options]>
 DESCRIPTION
 This script is used to evaluate pathway-level interactions for genetic networks

REQUIRED INPUTS
  --ssmFile=SSMFILE
    ssmFile is the SNP-SNP interaction file.

  --randRun=RANDRUN
     For random permutation run identification.
     randRun=0 for real run
     randRun=1,2,... for random runs
     This input is required.

  
OPTIONAL INPUTS
  --projectDir=PROJECTDIR
     projectDir is the data directory where results will also be stored at.
     Default is the current directory.

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
     --ssmFile=*) ssmFile=${argument/*=/""} ;;
     --minPath=*) minPath=${argument/*=/""} ;;
     --snpPerms=*) snpPerms=${argument/*=/""} ;;
     --binaryNetwork=*) binaryNetwork=${argument/*=/""} ;;
     --randRun=*) randRun=${argument/*=/""} ;;
     esac
done

# set defaults
if [ -z "${projectDir}" ]; then projectDir=`pwd`; fi
if [ -z "${minPath}" ]; then minPath=10; fi
if [ -z "${snpPerms}" ]; then snpPerms=10000; fi
if [ -z "${bnaryNetwork}" ]; then binaryNetwork=0; fi

# check required inuts
if [ -z "${ssmFile}" ]; then echo "Please specify ssmFile to proceed."; exit; fi
if [ -z "${randRun}" ]; then echo "Please specify randRun to proceed."; exit; fi

cd ${projectDir}

nice matlab -nodisplay -nodesktop -nosplash -r "genstats('${ssmFile}_R${randRun}','BPMind.mat',${binaryNetwork},${snpPerms},${minPath});exit" </dev/null> /dev/null
