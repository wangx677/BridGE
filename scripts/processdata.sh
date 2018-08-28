#!/bin/bash

ARGS=$@

_usage() {
###### USAGE ######
cat <<EOF
$*
Usage: processdata.sh <[options]>
DESCRIPTION  
This script is used to process data (from plink format) and prepare for 
running BridGE algorithm.

REQUIRED OPTIONS
  --plinkFile=PLINKFILE
	PLINK filename without the extension. 
	For example, your plink files are plinkexample.bim, plinkexample.bed, 
	plinkexample.fam. Please use "--plinkFile=plinkexample".

OPTIONAL OPTIONS
  --mind=MIND
	PLINK parameter used to control missing genotypes per individual. 
	Default is 0.02.

  --geno=GENO
	PLINK parameter used to control missing genotypes per SNP. 
	Default is 0.02. 

  --maf=MAF 
	PLINK parameter used to control allele frequency. 
	Default is 0.05.

  --hwe=HWE
	PLINK parameter used to control Hardy-Weinberg equilibrium. 
	Default is 0.000001.

  --pihat=PIHAT
	PIHAT: Identical twins, and duplicates, are 100% identical by descent (Pihat 1.0); 
		First-degree relatives are 50% IBD (Pihat 0.5); 
		Second-degree relatives are 25% IBD (Pihat 0.25); 
		Third-degree relatives are 12.5% equal IBD (Pihat 0.125); 
		4th degree Pihat = 0.0625; 5th degree Pihat=0.03125 level.        
	Default is 0.05.
	
  --matchCC=MATCHCC
	This parameter tells if cases and controls need to be matched.
	1: matched; 0: unmatched. 
	Default is 1.

  --ldWindow=LDWINDOW
	PLINK parameter used for linkage disequilibrium based SNP pruning. 
	Default is 50.

  --ldShift=LDSHIFT
	PLINK parameter used for linkage disequilibrium based SNP pruning. 
	Default is 5.

  --ldR2=LDR2
	PLINK parameter used for linkage disequilibrium based SNP pruning. 
	Default is 0.1. 

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
        --plinkFile=*) plinkFile=${argument/*=/""} ;;
        --geno=*) geno=${argument/*=/""} ;;
        --mind=*) mind=${argument/*=/""} ;;
        --maf=*) maf=${argument/*=/""} ;;
        --hwe=*) hwe=${argument/*=/""} ;;
        --pihat=*) pihat=${argument/*=/""} ;;
        --matchCC=*) matchCC=${argument/*=/""} ;;
        --ldWindow=*) ldWindow=${argument/*=/""} ;;
        --ldShift=*) ldShift=${argument/*=/""} ;;
        --ldR2=*) ldR2=${argument/*=/""} ;;
        --genesets=*) genesets=${argument/*=/""} ;;
        --minPath=*) minPath=${argument/*=/""} ;;
        --maxPath=*) maxPath=${argument/*=/""} ;;
        --geneAnnotation=*) geneAnnotation=${argument/*=/""} ;;
        --mappingDistance=*) mappingDistance=${argument/*=/""} ;;
        --help) _usage; exit;;
        esac
done

# set defaults
if [ -z "${plinkFile}" ]; then echo "Plink file is not provided"; exit; fi
if [ -z "${mind}" ]; then mind=0.02; fi
if [ -z "${geno}" ]; then geno=0.02; fi
if [ -z "${maf}" ]; then maf=0.05; fi
if [ -z "${hwe}" ]; then hwe=0.000001; fi
if [ -z "${pihat}" ]; then pihat=0.2; fi
if [ -z "${matchCC}" ]; then matchCC=1; fi
if [ -z "${ldWindow}" ]; then ldWindow=50; fi
if [ -z "${ldShift}" ]; then ldShift=5; fi
if [ -z "${ldR2}" ]; then ldR2=0.1; fi
if [ -z "${genesets}" ]; then genesets=${BRIDGEPATH}/refdata/c2.cp.v6.0.mat; fi
if [ -z "${minPath}" ]; then minPath=10; fi
if [ -z "${maxPath}" ]; then maxPath=300; fi

if [ -z "${geneAnnotation}" ]; then 
	geneAnnotation=${BRIDGEPATH}/refdata/gencode.v27.annotation.csv; 
fi

if [ -z "${mappingDistance}" ]; then mappingDistance=50000; fi

# qc based on mind geno maf and hwe
plink1.9 --bfile ${plinkFile} --noweb --allow-no-sex --mind ${mind} \
	--geno ${geno} --maf ${maf} --hwe ${hwe} --make-bed --out ${plinkFile}_tmp0 > /dev/null

# only keep autosomal chromosomes
extractchr1-22.sh ${plinkFile}_tmp0 ${plinkFile}_tmp1

# remove snps without correct genotype assigned
excludenogenotypesnps.sh ${plinkFile}_tmp1 ${plinkFile}_tmp2

# remove related individuals
removerelatedindividual.sh ${plinkFile}_tmp2 ${plinkFile}_tmp3 ${pihat}

# use size 2 clustering to match case and control
if [ "${matchCC}" -eq 1 ] 
then
     matchcasecontrol.sh ${plinkFile}_tmp3 gwas_data_all 
     mv ${plinkFile}_tmp3.cluster1 plinkFile.cluster1
     mv ${plinkFile}_tmp3.cluster2 plinkFile.cluster2
     mv ${plinkFile}_tmp3.cluster2.orig plinkFile.cluster2.orig
else
     plink1.9 --bfile ${plinkFile}_tmp3 --make-bed --out gwas_data_all
fi

# get the list of Genes included in the gene sets
nice matlab -nodisplay -nodesktop -nojvm -r "genesetglist('${genesets}'); \
	exit"  </dev/null> /dev/null
nice matlab -nodisplay -nodesktop -nojvm -r "genesetsantn('geneList_from_genesets', \
	'${geneAnnotation}','geneList_from_genesets_annotation'); \
	exit" </dev/null> /dev/null

# get the list of snps which can be mapped to the Genes
nice matlab -nodisplay -nodesktop -nojvm -r "mapsnp2gene('gwas_data_all.bim', \
	'geneList_from_genesets_annotation',$mappingDistance,'snplist', \
	'snp2keep_geneset');exit"  </dev/null> /dev/null  

# get less redundant SNP set
plink1.9 --bfile gwas_data_all --extract snp2keep_geneset --allow-no-sex \
	--indep-pairwise ${ldWindow} ${ldShift} ${ldR2} --noweb \
	--out gwas_data_all > /dev/null

# generate new plink data
plink1.9 --bfile gwas_data_all --allow-no-sex --extract gwas_data_all.prune.in \
	--make-bed --out gwas_data_final  > /dev/null

# conver plink to matlab
plink1.9 --bfile gwas_data_final --allow-no-sex --noweb --recodeA \
	--out RecodeA_file  > /dev/null

nice matlab -nodisplay -nodesktop -r "plink2mat('RecodeA_file.raw', \
	'gwas_data_final.bim','gwas_data_final.fam','gwas_data_final.mat');exit" \
	</dev/null> /dev/null

rm *tmp*

# conver SNP data from 0,1,2 format to 0,1 format
nice matlab -nodisplay -nodesktop -nosplash -r \
	"bindataar('gwas_data_final.mat');exit" </dev/null> /dev/null
nice matlab -nodisplay -nodesktop -nosplash -r \
	"bindataad('gwas_data_final.mat');exit" </dev/null> /dev/null

# prepare gene-set data
geneAnnotationFile=geneList_from_genesets_annotation

if [ "${mappingDistance}" -ge "1000" ]; then
     snpGeneMappingFile=snpgenemapping_$((${mappingDistance}/1000))kb.mat
else
     snpGeneMappingFile=snpgenemapping_${mappingDistance}bp.mat
fi

nice matlab -nodisplay -nodesktop -nojvm -r "mapsnp2gene('${plinkFile}.bim', \
	'${geneAnnotationFile}',$mappingDistance,'matrix','${snpGeneMappingFile}'); \
	exit"  </dev/null> /dev/null
	
nice matlab -nodisplay -nodesktop -r "snppathway('gwas_data_final.mat', \
	'${snpGeneMappingFile}','${genesets}',${minPath},${maxPath});exit"  \
	</dev/null> /dev/null

snpPathwayFile=snp_pathway_min${minPath}_max${maxPath}.mat

nice matlab -nodisplay -nodesktop -r "bpmind('${snpPathwayFile}');exit" \
	</dev/null> /dev/null
