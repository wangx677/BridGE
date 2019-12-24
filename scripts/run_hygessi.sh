#!/bin/bash

ARGS=$@

_usage() {
###### USAGE ######
cat <<EOF
$*
Usage: run_hygessi.sh <[options]>
DESCRIPTION
This script provides an alternative way (bridge.sh also can be used) to compute hygeSSI based 
genetic interaction networks. 

REQUIRED OPTIONS
  --randRun=RANDRUN
    This parameter is used to define if the interaction is calculated based on real data or 
    randomized data. 0 for real data, 1,2,3,...N are for randomized data

OPTIONAL OPTIONS
  --marginal=MARGINAL
     marginal=0: impact of joint mutation doesn't need to be greater than single SNPs
     marginal=1: impact of joint mutation need to be greater than single SNPs
     Default is 1.

  --plinkCluster2=PLINKCLUSTERCLUSTEI2
     Output from plink size-2 clustering after removing individuals that
     are not paired with others.
     Default is PlinkFile.cluster2.

  --nWorker=NWORKER
     Number of workers for parallel computing in matlab.
     Default is 10 

EOF
}


# check inputs
options=$@
arguments=${options}

# if [ -z "${arguments}" ]; then _usage; exit; fi

for argument in $options
do
        case $argument in
        --randRun=*) randRun=${argument/*=/""} ;;
        --marginal=*) marginal=${argument/*=/""} ;;
        --plinkCluster2=*) plinkCluster2=${argument/*=/""} ;;
        --nWorker=*) nWorker=${argument/*=/""} ;;
        --help) _usage; exit;;
        esac
done

# set defaults
if [ -z "${marginal}" ]; then marginal=1; fi
if [ -z "${plinkCluster2}" ]; then plinkCluster2=PlinkFile.cluster2; fi
if [ -z "${nWorker}" ]; then nWorker=10; fi

# compute SNP-SNP interaction networks use 10 matlab workers
nice matlab -nodisplay -nodesktop -nosplash -r "computessi('combined',${marginal},0.05,0.05,'${plinkCluster2}',${nWorker},${randRun});exit" </dev/null> /dev/null
