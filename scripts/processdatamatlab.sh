#!/bin/bash

ARGS=$@

_usage() {
###### USAGE ######
cat <<EOF
$*
Usage: processdatamatlab.sh <[options]>
DESCRIPTION  
This script is used to process data (from matlab format) and prepare for 
running BridGE algorithm.

REQUIRED OPTIONS
  --matlabFile=MATLABFILE
	Processed SNP data in MATLAB mat file. The .mat file consists a structure 
     array "SNPdata" with the following fields:
	   rsid:snp names
	   data:genotype data
	   chr: chromosome id
	   loc: physical location
	   pheno: sample's phenotype
	   fid: sample's family id
	   pid: sample id
	   gender: smaple's gender (optional)
	
OPTIONAL OPTIONS
  --genesets=GENESETS
	Gene sets in MATLAB mat file. The mat file has to include a 
	structure array 'GeneSet' with the following fields:
		gpMatrix: a binary matrix with rows as genes and columns as pathways
		PathwayNames: a cell array includes pathway names with exact same 
			order as in gp_matrix
		GeneNames: a cell array includes gene names with exact same order 
			as in gp_matrix
		EntrezIDs: a cell array includes gene entrez ids with exact same 
			order as in gp_matrix (optional)
	Default is the gene sets in <Your BridGE PATH>/refdata/c2.cp.v6.0.mat 
	(MsigDB CP: Canonical pathways version 6.0).
	
  --minPath=MINPATH
	Minimum number of genes in the gene set. 
	Default is 10.

  --maxPath=MAXPATH
	Maximum number of genes in the gene set. 
	Default is 300.

  --geneAnnotation=GENEANNOTATION
	Gene anotation file should be a text file in the following format: 
	one row per gene, each row has chromosome, start and stop positions 
	(base-pair) and then gene name:
		7 20140803 20223538 7A5
		19 63549983 63556677 A1BG
		10 52236330 52315441 A1CF
		8 43266741 43337485 A26A1
	Default is the <Your BridGE PATH>/refdata/gencode.v27.annotation.csv  
	downloaded from ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/.
	
  --mappingDistance=MAPPINGDISTANCE
	SNP to gene mapping distance in base-pair. This parameter is used to 
	add additional distance to the gene bounds. 
	Default is 50000 base-pair.

EOF
}


# check inputs
options=$@
arguments=${options}

if [ -z "${arguments}" ]; then _usage; exit; fi

for argument in $options
do
        case $argument in
        --matlabFile=*) matlabFile=${argument/*=/""} ;;
        --genesets=*) genesets=${argument/*=/""} ;;
        --minPath=*) minPath=${argument/*=/""} ;;
        --maxPath=*) maxPath=${argument/*=/""} ;;
        --geneAnnotation=*) geneAnnotation=${argument/*=/""} ;;
        --mappingDistance=*) mappingDistance=${argument/*=/""} ;;
        --help) _usage; exit;;
        esac
done

# set defaults
if [ -z "${genesets}" ]; then genesets=${BRIDGEPATH}/refdata/c2.cp.v6.0.mat; fi
if [ -z "${minPath}" ]; then minPath=10; fi
if [ -z "${maxPath}" ]; then maxPath=300; fi

if [ -z "${geneAnnotation}" ]; then 
	geneAnnotation=${BRIDGEPATH}/refdata/gencode.v27.annotation.csv; 
fi

if [ -z "${mappingDistance}" ]; then mappingDistance=50000; fi


# matlabFile is required
if [ -z "${matlabFile}" ]; then echo "matlabFile is not provided"; exit; fi

# conver SNP data from 0,1,2 format to 0,1 format

nice matlab -nodisplay -nodesktop -nosplash -r "bindataar('${matlabFile}'); \
	exit" </dev/null> /dev/null

nice matlab -nodisplay -nodesktop -nosplash -r "bindataad('${matlabFile}'); \
	exit" </dev/null> /dev/null

# get SNP information from matlabFile
nice matlab -nodisplay -nodesktop -nosplash -r "matsnpinfo('${matlabFile}', \
	'SNPdata.txt');exit" </dev/null> /dev/null

# prepare gene-set data
if [ "${mappingDistance}" -ge "1000" ]; then
     snpGeneMappingFile=snpgenemapping_$((${mappingDistance}/1000))kb.mat
else
     snpGeneMappingFile=snpgenemapping_${mappingDistance}bp.mat
fi

nice matlab -nodisplay -nodesktop -nojvm -r "mapsnp2gene('SNPdata.txt', \
	'${geneAnnotation}',$mappingDistance,'matrix','${snpGeneMappingFile}'); \
	exit"  </dev/null> /dev/null
	
nice matlab -nodisplay -nodesktop -r "snppathway('${matlabFile}', \
	'${snpGeneMappingFile}','${genesets}',${minPath},${maxPath});exit" \
	</dev/null> /dev/null

snpPathwayFile=snp_pathway_min${minPath}_max${maxPath}.mat

nice matlab -nodisplay -nodesktop -r "bpmind('${snpPathwayFile}'); \
	exit" </dev/null> /dev/null
