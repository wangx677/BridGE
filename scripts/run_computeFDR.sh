#!/bin/bash

_usage() {
###### USAGE ######
cat <<EOF
$*
 Usage: run_computeFDR.sh <[options]>
 DESCRIPTION
 This script is used to compute FDR for pathway-level interactions

REQUIRED INPUTS
  --ssmFile=SSMFILE
    ssmFile is the SNP-SNP interaction file.

  --samplePerms=SAMPLEPERMS
     number of sample permutations.
  
OPTIONAL INPUTS
  --projectDir=PROJECTDIR
     projectDir is the data directory where results will also be stored at.
     Default is the current directory.

  --minPath=MINPATH
     Minimum number of genes in the gene set.
     Default is 10.

  --pvalueCutoff=PVALUECUTOFF
    Compute FDRs for BPMs with permutation p-value less than pvalueCutoff
    Default is 0.005

  --bpmindFile=BPMINDFILE
    index file for pathway-level interactions
    Default is BPMind.mat

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
     --pvalueCutoff=*) pvalueCutoff=${argument/*=/""} ;;
     --samplePerms=*) samplePerms=${argument/*=/""} ;;
     --bpmindFile=*) bpmindFile=${argument/*=/""} ;;
     esac
done

# set defaults
if [ -z "${projectDir}" ]; then projectDir=`pwd`; fi
if [ -z "${minPath}" ]; then minPath=10; fi
if [ -z "${pvalueCutoff}" ]; then pvalueCutoff=0.005; fi
if [ -z "${bpmindFile}" ]; then bpmindFile=BPMind.mat; fi
if [ -z "${samplePerms}" ]; then samplePerms=10; fi

# check required inuts
if [ -z "${ssmFile}" ]; then echo "Please specify ssmFile to proceed."; exit; fi
if [ -z "${samplePerms}" ]; then echo "Please specify samplePerms to proceed."; exit; fi

cd ${projectDir}
nice matlab -nodisplay -nodesktop -nosplash -r "fdrsampleperm('${ssmFile}','${bpmindFile}',${pvalueCutoff},${minPath},${samplePerms});exit" </dev/null> /dev/null
