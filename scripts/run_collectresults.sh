#!/bin/bash

_usage() {
###### USAGE ######
cat <<EOF
$*
 Usage: run_collectresults.sh <[options]>
 DESCRIPTION
 This script is used to summarize results and list signficant BPM/WPM/PATH
 interactions with corresponding statsitics.

 INPUTS:
  --projectDir=PROJECTDIR
     Projectdir is the data directory where results will also be stored at.
     Default is the current directory.

  --interaction=INTERACTION
     hygeSSI: hypergeometric based interaction
     lrSSI: logistic regression based interaction
     Default is hygeSSI.

  --model=MODEL
     Prefered interaction model which can be chosen from:
          RR recessive-recessive interaction
          DD dominant-dominant interaction
          RD recessive-dominant interaction
          combined combined all three interactions
          AA additive-additive interaction
     Default is combined.
    
  --validationDir=VALIDATIONDIR
     Validation is the data directory where results from validation dataset
     will also be stored at.
     Optional.

  --fdrCutoff=FDRCUT
    Fdrcutoff is the FDR cutoff used to determine the signficant interactions.
    Default is 0.4

  --bpmind=BPMIND
    BPMind is a mat file that indicates SNP to BPM, SNP to WPM membership.
    Default is BPMind.mat

  --snpPathwayFile=SNP2PATHWAY
    Snp2pathway is a mat file that indicates SNP to pathway membership.
    Default is snp_pathway_min10_max300.mat

  --snpGeneMappingFile=SNP2GENEMAPPING
    snpGeneMappingFile is a mat file that indicate SNP to gene membership.
    Default is snpgenemapping_50kb.mat

EOF
}

# check inputs
options=$@
arguments=${options}

if [ -z "${arguments}" ]; then _usage; exit; fi

for argument in $options
do
     case $argument in
     --projectDir=*) projectDir=${argument/*=/""} ;;
     --interaction=*) interaction=${argument/*=/""} ;;
     --model=*) model=${argument/*=/""} ;;
     --validationDir=*) validationDir=${argument/*=/""} ;;
     --fdrCutoff=*) fdrCutoff=${argument/*=/""} ;;
     --bpmind=*) bpmind=${argument/*=/""} ;;
     --snpPathwayFile=*) snpPathwayFile=${argument/*=/""} ;;
     --snpGeneMappingFile=*) snpGeneMappingFile=${argument/*=/""} ;;
     --help) _usage; exit;;
     esac
done

# set defaults
if [ -z "${projectDir}" ]; then printf "Please specify projectDir. \n"; exit; fi
if [ -z "${interaction}" ]; then interaction=hygeSSI; fi
if [ -z "${model}" ]; then printf "Please specify model. \n"; exit; fi
if [ -z "${fdrCutoff}" ]; then fdrCutoff=0.4; fi
if [ -z "${bpmind}" ]; then bpmind=BPMind.mat; fi
if [ -z "${snpPathwayFile}" ]; then snpPathwayFile=snp_pathway_min10_max300.mat; fi
if [ -z "${snpGeneMappingFile}" ]; then snpGeneMappingFile=snpgenemapping_50kb.mat; fi


if [ "$interaction" == "hygeSSI" ]; then
    resultfile=results_ssM_hygeSSI_alpha10.05_alpha20.05_${model}_R0.mat
    ssmFile=ssM_hygeSSI_alpha10.05_alpha20.05_${model}_R0.mat
elif [ "$interaction" == "regSSI" ]; then
    resultfile=results_ssM_regSSI_${model}_R0.mat
    ssmFile=ssM_regSSI_${model}_R0.mat
elif [ "$interaction" == "lrSSI" ]; then
     resultfile=results_ssM_lr_cassi_pv0.05_R0.mat
     ssmFile=ssM_lr_cassi_pv0.05_R0.mat
fi

if [ -z "${validationDir}" ]; then
     matbg "collectresults('${resultfile}',${fdrCutoff},'${ssmFile}','${bpmind}','${snpPathwayFile}','${snpGeneMappingFile}')" matbg_$model.log
else
     if [ "$interaction" == "hygeSSI" ]; then
         validationfile=${BRIDGEPATH}/${validationDir}/results_ssM_hygeSSI_alpha10.05_alpha20.05_${model}_R0.mat
     elif [ "$interaction" == "regSSI" ]; then
         validationfile=${BRIDGEPATH}/${validationDir}/results_ssM_regSSI_${model}_R0.mat
     elif [ "$interaction" == "lrSSI" ]; then
          validationfile=${BRIDGEPATH}/${validationDir}/results_ssM_lr_cassi_pv0.05_R0.mat
     fi
     matbg "collectresults('${resultfile}',${fdrCutoff},'${ssmFile}','${bpmind}','${snpPathwayFile}','${snpGeneMappingFile}','${validationfile}')" matbg_$model.log
fi

