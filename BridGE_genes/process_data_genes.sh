#!/bin/bash
projectdir=$1 # project directory
genelist=$2 # genelist is a text file with each line is a gene
gwasfile=$3 # gwas plink file without file extension 
mappingDistance=$4 # SNP to gene mapping distance, if 50k, use input 50000

if [ ! -d "${projectdir}" ]; then 
     mkdir ${projectdir}
fi

cd ${projectdir}

# get gene information
if [ -f "gene_annotation" ]; then
     rm gene_annotation
fi

for genename in `cat ${genelist}`
do
     grep -w "${genename}$" /project/csbio/wwang/BridGE/refdata/glist-hg38 >> gene_annotation
done

# get SNP list with predefined mapping distance
geneAnnotation=gene_annotation
option=snplist

snpAnnotation=${gwasfile}.bim
z=$((${mappingDistance}/1000))kb
outputFile=snplist_$z

if [ -f "${outputFile}" ]; then
     rm ${outputFile}
fi

nice matlab -nodisplay -nodesktop -nosplash -r " mapsnp2gene('${snpAnnotation}','${geneAnnotation}',${mappingDistance},'${option}','${outputFile}');exit" </dev/null> /dev/null

# prepare plink data associated with genes of interest
plink1.9 --bfile ${gwasfile} --extract snplist_$z --allow-no-sex --make-bed --out gwas_data_final

# convert plink to matlab
plink1.9 --bfile gwas_data_final --allow-no-sex --recodeA --out RecodeA_file

file=gwas_data_final
nice matlab -nodisplay -nodesktop -r "plink2mat('RecodeA_file.raw','${file}.bim','${file}.fam','${file}.mat');exit" </dev/null> /dev/null

cp -p ${gwasfile}.cluster2 gwas_data_final.cluster2

nice matlab -nodisplay -nodesktop -nosplash -r "bindataar('gwas_data_final.mat');exit" </dev/null> /dev/null
nice matlab -nodisplay -nodesktop -nosplash -r "bindataad('gwas_data_final.mat');exit" </dev/null> /dev/null

# compute LD using plink
plink1.9 --bfile gwas_data_final --allow-no-sex --ld-window-kb 300 --list-all --r2

