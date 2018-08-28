#!/bin/bash

_usage() {
###### USAGE ######
cat <<EOF
$*
 Usage: runsampleperm_serial_pre.sh <[options]>
 DESCRIPTION  
 This script is used to run specific sample permutation runs.
 INPUTS:
  --projectDir=PROJECTDIR
	Projectdir is the data directory where results will also be stored at.
	Default is the current directory.

  --nWorker=NWORKER
     Number of workers for parallel computing in matlab.
     Default is 10 or (number of available CPUs - 2).

  --minPath=MINPATH
	Minimum number of genes in the gene set. 
	Default is 10.

  --interaction=INTERACTION
     hygeSSI: hypergeometric based interaction
     lrSSI: logistic regression based interaction
     Default is hygeSSI.

  --model=MODEL
	Prefered interaction model which can be chosen from: 
		1 recessive-recessive interaction
		2 dominant-dominant interaction
		3 recessive-dominant interaction
		4 combined all three interactions
		5 additive-additive interaction
	This input is required.

  --alpha1=ALPHA1
     Marginal significance level for joint effect.
     Defalut is 0.05.

  --alpha2=ALPHA2
     Marginal significance level for individual effects.
     Default is 0.05.

  --snpPerms=SNPPERMS
     Number of SNP shuffling randomization.
     Default is 100000.

  --netDensity=NETDENSITY
	Network binarization density based on considering both protecitve/risk interactions.
	This input is required.

  --plinkCluster2=PLINKCLUSTERCLUSTEI2
	Output from plink size-2 clustering after removing individuals that 
	are not paired with others.
	Default is PlinkFile.cluster2.

  --randRun=RANDRUN
     For random permutation run identification.
     R=0 for real run
     R=1,2,... for random runs
     This input is required.
     
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
     --interaction=*) interaction=${argument/*=/""} ;;
     --minPath=*) minPath=${argument/*=/""} ;;
     --snpPerms=*) snpPerms=${argument/*=/""} ;;
     --plinkCluster2=*) plinkCluster2=${argument/*=/""} ;;
     --model=*) model=${argument/*=/""} ;;
     --alpha1=*) alpha1=${argument/*=/""} ;;
     --alpha2=*) alpha2=${argument/*=/""} ;;
     --netDensity=*) netDensity=${argument/*=/""} ;;
     --nWorker=*) nWorker=${argument/*=/""} ;;
     --randRun=*) randRun=${argument/*=/""} ;;
     --help) _usage; exit;;
     esac
done

# set defaults
if [ -z "${projectDir}" ]; then projectDir=`pwd`; fi
if [ -z "${interaction}" ]; then interaction=hygeSSI; fi
if [ -z "${minPath}" ]; then minPath=10; fi
if [ -z "${snpPerms}" ]; then snpPerms=100000; fi
if [ -z "${plinkCluster2}" ]; then plinkCluster2=plinkFile.cluster2; fi
if [ -z "${alpha1}" ]; then alpha1=0.05; fi
if [ -z "${alpha2}" ]; then alpha2=0.05; fi
if [ -z "${nWorker}" ]; then nWorker=10; fi


cd ${projectDir}/

# check inputs
if [ -z "${netDensity}" ]; then echo 'netDensity is not defined'; exit; fi
if [ -z "${randRun}" ]; then echo 'randRun is not defined'; exit; fi
if [ -z "${model}" ]; then echo 'model is not defined'; exit; fi

# define ssmFile based on given inputs
if [ "$model" -eq "1" ]; then
     ssmFile=ssM_hygeSSI_alpha1${alpha1}_alpha2${alpha2}_RR
elif [ "$model" -eq "2" ]; then
     ssmFile=ssM_hygeSSI_alpha1${alpha1}_alpha2${alpha2}_DD
elif [ "$model" -eq "3" ]; then
     ssmFile=ssM_hygeSSI_alpha1${alpha1}_alpha2${alpha2}_RD
elif [ "$model" -eq "4" ]; then
     ssmFile=ssM_hygeSSI_alpha1${alpha1}_alpha2${alpha2}_combined
elif [ "$model" -eq "5" ]; then
     ssmFile=ssM_lr_cassi_pv0.05
fi


# start permutation run
bpmindFile=BPMind.mat
permBlockSize=10000

printf "BridGE is running sample permutaion with the followng parameters...\n"
printf "runsampleperm_serial.sh --interaction=${interaction} --model=${model} --netDensity=${netDensity} \n"
printf "          --snpPerms=${snpPerms} --alpha1=${alpha1} --alpha2=${alpha2}\n"
printf "          --plinkCluster2=${plinkCluster2} --nWorker=${nWorker} --randRun=${randRun}\n\n"


printf "computing SNP-SNP interaction for run #${randRun} ...\n\n"
nice matlab -nodisplay -nodesktop -nosplash -r "computessi(${model},${alpha1},${alpha2},'${plinkCluster2}',${nWorker},${randRun});exit" </dev/null> /dev/null

printf "Preparing SNP serial permutaton for sample permutation run #${randRun} ...\n\n"
nice matlab -nodisplay -nodesktop -nosplash -r "sampleperm_serial_pre('${ssmFile}','${bpmindFile}',${netDensity},${snpPerms},${permBlockSize},${minPath},${nWorker},${randRun});exit" </dev/null> /dev/null

