#!/bin/bash

_usage() {
###### USAGE ######
cat <<EOF
$*
 Usage: run_collectresults.sh <[options]>
 DESCRIPTION
 This script is used to summarize results and list signficant BPM/WPM/PATH
 interactions with corresponding statsitics.

 REQUIRED INPUTS:
 --ssmFile=SSMFILE
   SNP-SNP interaction file name


 OPTIONAL INPUTS:
  --projectDir=PROJECTDIR
     Projectdir is the data directory where results will also be stored at.
     Default is the current directory.
    
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

  --snpPathwayFile=snpPathwayFile
    Snp2pathway is a mat file that indicates SNP to pathway membership.
    Default is snp_pathway_min10_max300.mat

  --snpGeneMappingFile=SNP2GENEMAPPING
    snpGeneMappingFile is a mat file that indicate SNP to gene membership.
    Default is snpgenemapping_50kb.mat

  --densityCutoff=DENSITYCUTOFF
    densityCutoff is a density cutoff used to binarize SNP-SNP interaction network for 
    identification of driver SNPs/genes. For example, densityCutoff=0.01 means using top 1%
    of the interactions to derive driver SNPs and genes.
    if densityCutoff is not setting, then a SNP-SNP interaction score cutoff of 0.2 will be 
    used to binarize the network. This is default for mhygeSSI and hygeSSI networks.
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
     --ssmFile=*) ssmFile=${argument/*=/""} ;;
     --validationDir=*) validationDir=${argument/*=/""} ;;
     --fdrCutoff=*) fdrCutoff=${argument/*=/""} ;;
     --bpmind=*) bpmind=${argument/*=/""} ;;
     --snpPathwayFile=*) snpPathwayFile=${argument/*=/""} ;;
     --snpGeneMappingFile=*) snpGeneMappingFile=${argument/*=/""} ;;
     --densityCutoff=*) densityCutoff=${argument/*=/""} ;;
     --help) _usage; exit;;
     esac
done

# set defaults
if [ -z "${projectDir}" ]; then projectDir=`pwd`; fi
if [ -z "${fdrCutoff}" ]; then fdrCutoff=0.4; fi
if [ -z "${bpmind}" ]; then bpmind=BPMind.mat; fi
if [ -z "${snpPathwayFile}" ]; then snpPathwayFile=snp_pathway_min10_max300.mat; fi
if [ -z "${snpGeneMappingFile}" ]; then snpGeneMappingFile=snpgenemapping_50kb.mat; fi

# check required inputs
if [ -z "${ssmFile}" ]; then echo "Please specify ssmFile to proceed."; exit; fi


if [ -z "${validationDir}" ]; then
    if [ -z "${densityCutoff}" ]; then
     nice matlab -nodisplay -nodesktop -nosplash -r "collectresults('results_${ssmFile}',${fdrCutoff},'${ssmFile}','${bpmind}','${snpPathwayFile}','${snpGeneMappingFile}')" </dev/null> /dev/null
     else
         nice matlab -nodisplay -nodesktop -nosplash -r "collectresults('results_${ssmFile}',${fdrCutoff},'${ssmFile}','${bpmind}','${snpPathwayFile}','${snpGeneMappingFile}',${densityCutoff})" </dev/null> /dev/null

     fi
else
     validationfile=${BRIDGEPATH}/${validationDir}/results_${ssmFile}
     if [ -z "${densityCutoff}" ]; then
          nice matlab -nodisplay -nodesktop -nosplash -r "collectresults('results_${ssmFile}',${fdrCutoff},'${ssmFile}','${bpmind}','${snpPathwayFile}','${snpGeneMappingFile}','${validationfile}')" </dev/null> /dev/null
      else
          nice matlab -nodisplay -nodesktop -nosplash -r "collectresults('results_${ssmFile}',${fdrCutoff},'${ssmFile}','${bpmind}','${snpPathwayFile}','${snpGeneMappingFile}','${validationfile}',${densityCutoff})" </dev/null> /dev/null

      fi
fi

