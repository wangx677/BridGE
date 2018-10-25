#!/bin/bash

_usage() {
###### USAGE ######
cat <<EOF
$*
 Usage: bridge.sh <[options]>
 DESCRIPTION  
 This script is used to process data and run BridGE algorithm to identify 
 between-pathway model(BPM), within-pathway model (WPM), and interaction 
 hub pathways (PATH).

 GENERAL INPUTS:
  --job=JOB
	DataProcess: prepare data to run BridGE
     ComputeInteraction: compute SNP-SNP interaction
	SamplePermutation: run sample shuffling permutation
	Analysis: run post analysis 
	This input is always required.
	
  --projectDir=PROJECTDIR
	Projectdir is the data directory where results will also be stored at.
	Default is the current directory.

  --nWorker=NWORKER
     Number of workers for parallel computing in matlab.
     Default is 10 or (number of available CPUs - 2).

 DATA PROCESS INPUTS:
  --plinkFile=PLINKFILE
	PLINK filename without the extension. 
	For example, your plink files are plinkexample.bim, plinkexample.bed, 
	plinkexample.fam. Please use "--plinkFile=plinkexample".
	This input is required when "--job=DataProcess" and "matlabFile" is not set.

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
     pihat=1 means no filtering based on sample relatedness
	Default is 0.2.

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
	structure array 'geneset' with the following fields:
		gpatrix: a binary matrix with rows as genes and columns as pathways
		pathwaynames: a cell array includes pathway names with exact same 
			order as in gp_matrix
		genenames: a cell array includes gene names with exact same order 
			as in gp_matrix
		entrezids: a cell array includes gene entrez ids with exact same 
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
     This is an optional input. If a matlabFile is given, when "job=DataProcess",
     bridge.sh will directly use the given data to map SNPs to pathways.

COMPUTE INTERACTION INPUTS:
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
	This input is not required when "--job=DataProcess"

  --marginal==MARGINAL
     marginal=0: impact of joint mutation doesn't need to be greater than single SNPs
     marginal=1: impact of joint mutation need to be greater than single SNPs
     Default is 0.

  --alpha1=ALPHA1
     Marginal significance level for joint effect.
     Defalut is 0.05.

  --alpha2=ALPHA2
     Marginal significance level for individual effects.
     Default is 0.05.

  --samplePerms=SAMPLEPERMS
     Number of sample shuffling randomizations.
     Default is 10.

SAMPLE PERMUTATION INPUTS:
  --samplePerms=SAMPLEPERMS
	Number of sample shuffling randomizations. 
	Default is 10.

  --snpPerms=SNPPERMS
     Number of SNP shuffling randomization.
     Default is 10000.

  --binaryNetwork=BINARYNETWORK
    binaryNetwork=1: run BridGE based on binarized network
    binaryNetwork=0: run BridGE based on weighted network
    Default is binaryNetwork=0

  --netDensity=NETDENSITY
	Network binarization density based on considering both protecitve/risk interactions.
     It is required when binaryNetwork=1

  --plinkCluster2=PLINKCLUSTERCLUSTEI2
	Output from plink size-2 clustering after removing individuals that 
	are not paired with others.
	Default is PlinkFile.cluster2.

ANALYSIS INPUTS:
  --ssmFile=SSMFILE
	SNP-SNP interaction file name without suffix '_R<R>.mat'

  --fdrCutoff=FDRCUTOFF
     Significant FDR cutoff.
     Information for BPM, WPM, and PATH that pass this significant
     threshold will be collected and write into an excel file
     Default is 0.4.

  --pvalueCutoff=PVALUECUTOFF
    Compute FDRs for BPMs with permutation p-value less than pvalueCutoff
    Default is 0.005

  --validationDir=VALIDATIONDIR
     Validation is the data directory where results from validation dataset
     will also be stored at.
     Optional.
     
EOF
}

# check inputs
options=$@
arguments=${options}

if [ -z "${arguments}" ]; then _usage; exit; fi

for argument in $options
do
     case $argument in
     --job=*) job=${argument/*=/""} ;;
     --projectDir=*) projectDir=${argument/*=/""} ;;
     --plinkFile=*) plinkFile=${argument/*=/""} ;;
     --interaction=*) interaction=${argument/*=/""} ;;
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
     --geneAnnotation=*) geneAnnotation=${argument/*=/""} ;;
     --mappingDistance=*) mappingDistance=${argument/*=/""} ;;
     --minPath=*) minPath=${argument/*=/""} ;;
     --maxPath=*) maxPath=${argument/*=/""} ;;
     --snpPerms=*) snpPerms=${argument/*=/""} ;;
     --samplePerms=*) samplePerms=${argument/*=/""} ;;
     --plinkCluster2=*) plinkCluster2=${argument/*=/""} ;;
     --model=*) model=${argument/*=/""} ;;
     --marginal=*) marginal=${argument/*=/""} ;;
     --alpha1=*) alpha1=${argument/*=/""} ;;
     --alpha2=*) alpha2=${argument/*=/""} ;;
     --binaryNetwork=*) binaryNetwork=${argument/*=/""} ;;
     --netDensity=*) netDensity=${argument/*=/""} ;;
     --matlabFile=*) matlabFile=${argument/*=/""} ;;
     --ssmFile=*) ssmFile=${argument/*=/""} ;;
     --nWorker=*) nWorker=${argument/*=/""} ;;
     --fdrCutoff=*) fdrCutoff=${argument/*=/""} ;;
     --pvalueCutoff=*) pvalueCutoff=${argument/*=/""} ;;
     --help) _usage; exit;;
     esac
done

# set defaults
if [ -z "${job}" ]; then printf "Please specify which job to proceed."; exit; fi
if [ -z "${projectDir}" ]; then projectDir=`pwd`; fi
if [ -z "${interaction}" ]; then interaction=hygeSSI; fi
if [ -z "${mind}" ]; then mind=0.02; fi
if [ -z "${geno}" ]; then geno=0.02; fi
if [ -z "${maf}" ]; then maf=0.05; fi
if [ -z "${hwe}" ]; then hwe=0.000001; fi
if [ -z "${pihat}" ]; then pihat=0.2; fi
if [ -z "${matchCC}" ]; then matchCC=1; fi
if [ -z "${ldWindow}" ]; then ldWindow=50; fi
if [ -z "${ldShift}" ]; then ldShift=5; fi
if [ -z "${ldR2}" ]; then ldR2=0.1; fi

if [ -z "${genesets}" ]; then 
	genesets=${BRIDGEPATH}/refdata/c2.cp.v6.0.mat; 
fi

if [ -z "${geneAnnotation}" ]; then 
	geneAnnotation=${BRIDGEPATH}/refdata/gencode.v27.annotation.csv; 
fi

if [ -z "${mappingDistance}" ]; then mappingDistance=50000; fi
if [ -z "${minPath}" ]; then minPath=10; fi
if [ -z "${maxPath}" ]; then maxPath=300; fi
if [ -z "${snpPerms}" ]; then snpPerms=10000; fi
if [ -z "${samplePerms}" ]; then samplePerms=10; fi
if [ -z "${plinkCluster2}" ]; then plinkCluster2=plinkFile.cluster2; fi
if [ -z "${marginal}" ]; then marginal=0; fi
if [ -z "${alpha1}" ]; then alpha1=0.05; fi
if [ -z "${alpha2}" ]; then alpha2=0.05; fi
if [ -z "${fdrCutoff}" ]; then fdrCutoff=0.4; fi
if [ -z "${pvalueCutoff}" ]; then pvalueCutoff=0.005; fi
if [ -z "${nWorker}" ]; then nWorker=1; fi
if [ -z "${binaryNetwork}" ]; then binaryNetwork=0; fi
if [ "${interaction}" = "hygeSSI" ] && [ -z "${model}" ]; then model=combined; fi
if [ "${interaction}" = "lrSSI" ] && [ -z "${model}" ]; then model=AA; fi
if [ "${mdel}" = "AA" ]; then interaction=lrSSI; fi

# define file names
snpPathwayFile=snp_pathway_min${minPath}_max${maxPath}.mat
bpmindFile=BPMind.mat

if [ "${mappingDistance}" -ge "1000" ]; then
     snpGeneMappingFile=snpgenemapping_$((${mappingDistance}/1000))kb.mat
else
     snpGeneMappingFile=snpgenemapping_${mappingDistance}bp.mat
fi

case "${job}" in

# Plink data processing and  get snp-pathway matrix
DataProcess)
     if [ -z "${plinkFile}" ] && [ -z "${matlabFile}" ]; then 
		printf "Data file (plinkFile or matlabFile) is not provided"; exit; 
	fi

     cd ${projectDir}
     echo ${pihat} > pihat

     if [ -n "${plinkFile}" ]; then
          printf "BridGE is processing data ${plinkFile} by calling ... \n"
          printf "processdata.sh --job=DataProcess --projectDir=${projectDir} --plinkFile=${plinkFile} \n"
          printf "               --mind=${mind} --geno=${geno} --maf=${maf} --hwe=${hwe} --pihat=${pihat} \n"
          printf "               --matchCC=${matchCC} --ldWindow=${ldWindow} --ldShift=${ldShift} --ldR2=${ldR2} \n"
          printf "               --genesets=${genesets} --minPath=${minPath} --maxPath=${maxPath} \n"
          printf "               --geneAnnotation=${geneAnnotation} --mappingDistance=${mappingDistance} \n\n"

          processdata.sh --plinkFile=${plinkFile} --mind=${mind} --geno=${geno} --maf=${maf} --hwe=${hwe} \
               --pihat=${pihat} --matchCC=${matchCC} --ldWindow=${ldWindow} --ldShift=${ldShift} \
               --ldR2=${ldR2}  --genesets=${genesets} --minPath=${minPath} --maxPath=${maxPath} \
               --mappingDistance=${mappingDistance} --geneAnnotation=${geneAnnotation}

     elif [ -n "${matlabFile}" ]; then
          printf "BridGE is processing data ${matlabFile} by calling ...\n"
          printf "processdatamatlab.sh --matlabFile=${matlabFile} \n"
          printf "                     --genesets=${genesets} --minPath=${minPath} --maxPath=${maxPath} \n"
          printf "                     --geneAnnotation=${geneAnnotation} --mappingDistance=${mappingDistance} \n\n"

          processdatamatlab.sh --matlabFile=${matlabFile} --genesets=${genesets} --minPath=${minPath} \
                --maxPath=${maxPath} --geneAnnotation=${geneAnnotation} --mappingDistance=${mappingDistance}
     fi

     printf "Data processing is complete!\n"

     ;;

# compute SNP-SNP interaction
ComputeInteraction)
     cd ${projectDir}/

     printf "BridGE is computing SNP-SNP interaction with the followng parameters...\n"
     printf "bridge.sh --job=ComputeInteraction --projectDir=${projectDir} --interaction=${interaction} \n"
     printf "          --model=${model} --margnal=${marginal} --alpha1=${alpha1} --alpha2=${alpha1} --plinkCluster2=${plinkCluster2} \n"
     printf "          --samplePerms=${samplePerms} --nWorker=${nWorker}\n\n"
     
     for R in `seq 0 ${samplePerms}`
     do
          nice matlab -nodisplay -nodesktop -nosplash -r "computessi('${model}',${marginal},${alpha1},${alpha1},'plinkFile.cluster2',${nWorker},${R});exit" </dev/null> /dev/null
     done

     printf "Computing SNP-SNP interaction is complete!\n"

     ;;

# run sample permutation
SamplePermutation)
     cd ${projectDir}/

     if [ "${model}" = "RR" ]; then
          ssmFile=ssM_hygeSSI_alpha1${alpha1}_alpha2${alpha2}_RR
     elif [ "${model}" = "DD" ]; then
          ssmFile=ssM_hygeSSI_alpha1${alpha1}_alpha2${alpha2}_DD
     elif [ "${model}" = "RD" ]; then
          ssmFile=ssM_hygeSSI_alpha1${alpha1}_alpha2${alpha2}_RD
     elif [ "${model}" = "combined" ]; then
          ssmFile=ssM_hygeSSI_alpha1${alpha1}_alpha2${alpha2}_combined
     elif [ "${model}" = "AA" ]; then
          ssmFile=ssM_lr_cassi_pv0.05
     fi

     if [ -z "${bpmindFile}" ]; then bpmindFile='BPMind.mat'; fi

     printf "BridGE is running sample permutaion with the followng parameters...\n"
     printf "bridge.sh --job=SamplePermutation --projectDir=${projectDir} --interaction=${interaction} --model=${model} --marginal=${marginal} \n"
     printf "          --binaryNetwork=${binaryNetwork} --alpha1=${alpha1} --alpha2=${alpha1} --plinkCluster2=${plinkCluster2} \n"
     printf "          --samplePerms=${samplePerms} --snpPerms=${snpPerms} --minPath=${minPath} --nWorker=${nWorker}\n\n"

     for R in `seq 0 ${samplePerms}`
     do
          printf "computing SNP-SNP interaction for run #${R} ...\n\n"
          nice matlab -nodisplay -nodesktop -nosplash -r "computessi('${model}',${marginal},${alpha1},${alpha2},'${plinkCluster2}',${nWorker},${R});exit" </dev/null> /dev/null

          printf "calling SNP permutation for sample permutation run #${R} ...\n\n"
          nice matlab -nodisplay -nodesktop -nosplash -r "genstats('${ssmFile}_R${R}','${bpmindFile}',${binaryNetwork},${snpPerms},${minPath});exit" </dev/null> /dev/null
     done
     
     printf "Sample permutation is complete!\n"
     ;;

# post analysis
Analysis)
     cd ${projectDir}/
     if [ -z "${ssmFile}" ];then
          if [ "${model}" = "RR" ]; then
               ssmFile=ssM_hygeSSI_alpha1${alpha1}_alpha2${alpha2}_RR
          elif [ "${model}" = "DD" ]; then
               ssmFile=ssM_hygeSSI_alpha1${alpha1}_alpha2${alpha2}_DD
          elif [ "${model}" = "RD" ]; then
               ssmFile=ssM_hygeSSI_alpha1${alpha1}_alpha2${alpha2}_RD
          elif [ "${model}" = "combined" ]; then
               ssmFile=ssM_hygeSSI_alpha1${alpha1}_alpha2${alpha2}_combined
          elif [ "${model}" = "AA" ]; then
               ssmFile=ssM_lr_cassi_pv0.05
          fi
     fi

     printf "BridGE is running post analysis with the followng parameters...\n"
     printf "bridge.sh --job=Analysis --projectDir=${projectDir} --ssmFile=${ssmFile} --interaction=${interaction} \n"
     printf "          --model=${model}  --samplePerms=${samplePerms} --snpPerms=${snpPerms} --validationDir=${validationDir} \n"
     printf "          --fdrCutoff=${fdrCutoff} --pvalueCutoff=${pvalueCutoff} --minPath=${minPath} --snpGeneMappingFile=${snpGeneMappingFile} \n"     
     printf "          --snpPathwayFile=${snpPathwayFile} \n"

     if [ -z "${bpmindFile}" ]; then bpmindFile='BPMind.mat'; fi
 
     # compute false discovery rate based on sample permutation
     # output file is results_<ssmFile>_R0.mat
     printf "computing false discovery rate ...\n\n"
     nice matlab -nodisplay -nodesktop -nosplash -r "fdrsampleperm('${ssmFile}','${bpmindFile}',${pvalueCutoff},${minPath},${samplePerms});exit" </dev/null> /dev/null

     # write signficant results to excel file
     # excel file is named "output_results_<ssmFile>_R0.xls"
     # there is also a corresponding mat-file.
     printf "getting detailed information for significant discoveries and writing to excel ...\n\n"
     run_collectresults.sh --projectDir=${projectDir} --interaction=${interaction} --model=${model} --validationDir=${validationDir} \
               --fdrCutoff=${fdrCutoff} --snpPathwayFile=${snpPathwayFile} --snpGeneMappingFile=${snpGeneMappingFile}
     
     printf "Post analysis is complete!\n"
     ;;

*)
     printf "Unknown job request"
esac
