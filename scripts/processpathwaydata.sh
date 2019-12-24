#!/bin/bash

ARGS=$@

_usage() {
###### USAGE ######
cat <<EOF
$*
Usage: processpathwaydata.sh <[options]>
DESCRIPTION  
This script is used to process pathway data and prepare for 
running BridGE algorithm.

OPTIONAL OPTIONS
  --plinkFile=PLINKFILE
     PLINK filename without the extension.
     For example, your plink files are plinkexample.bim, plinkexample.bed,
     plinkexample.fam. Please use "--plinkFile=plinkexample".

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

# if [ -z "${arguments}" ]; then _usage; exit; fi

for argument in $options
do
        case $argument in
        --plinkFile=*) plinkFile=${argument/*=/""} ;;
        --genesets=*) genesets=${argument/*=/""} ;;
        --minPath=*) minPath=${argument/*=/""} ;;
        --maxPath=*) maxPath=${argument/*=/""} ;;
        --geneAnnotation=*) geneAnnotation=${argument/*=/""} ;;
        --mappingDistance=*) mappingDistance=${argument/*=/""} ;;
        --help) _usage; exit;;
        esac
done

# set defaults
if [ -z "${plinkFile}" ]; then plinkFile=gwas_data_final; fi
if [ -z "${genesets}" ]; then genesets=${BRIDGEPATH}/refdata/c2.cp.v6.0.mat; fi
if [ -z "${minPath}" ]; then minPath=10; fi
if [ -z "${maxPath}" ]; then maxPath=300; fi

if [ -z "${geneAnnotation}" ]; then 
	geneAnnotation=${BRIDGEPATH}/refdata/gencode.v27.annotation.csv; 
fi

if [ -z "${mappingDistance}" ]; then mappingDistance=50000; fi


# get the list of Genes included in the gene sets
nice matlab -nodisplay -nodesktop -nojvm -r "genesetglist('${genesets}'); \
	exit"  </dev/null> /dev/null
nice matlab -nodisplay -nodesktop -nojvm -r "genesetsantn('geneList_from_genesets', \
	'${geneAnnotation}','geneList_from_genesets_annotation'); \
	exit" </dev/null> /dev/null

# map SNPs to Genes
nice matlab -nodisplay -nodesktop -nojvm -r "mapsnp2gene('${plinkFile}.bim', \
	'geneList_from_genesets_annotation',$mappingDistance,'snplist', \
	'snp2keep_geneset');exit"  </dev/null> /dev/null  

# prepare gene-set data
# define output file name
if [ "${mappingDistance}" -ge "1000" ]; then
     snpGeneMappingFile=snpgenemapping_$((${mappingDistance}/1000))kb.mat
else
     snpGeneMappingFile=snpgenemapping_${mappingDistance}bp.mat
fi

# map snp to genes
nice matlab -nodisplay -nodesktop -nojvm -r "mapsnp2gene('gwas_data_final.bim', \
'${geneAnnotation}',$mappingDistance,'matrix','${snpGeneMappingFile}'); \
exit"  </dev/null> /dev/null
	
# get snp to pathway relationship
nice matlab -nodisplay -nodesktop -r "snppathway('gwas_data_final.mat', \
'${snpGeneMappingFile}','${genesets}',${minPath},${maxPath});exit"  \
</dev/null> /dev/null

# define output file name
snpPathwayFile=snp_pathway_min${minPath}_max${maxPath}.mat

# get SNP indexes for between/within pathway interactions
nice matlab -nodisplay -nodesktop -r "bpmind('${snpPathwayFile}');exit" \
	</dev/null> /dev/null
