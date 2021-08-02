#!/bin/bash

ARGS=$@

_usage() {
###### USAGE ######
cat <<EOF
$*
        Usage: qc.sh <[options]>
        DESCRIPTION  
        This script is used to process plink data based on standard procedure.

        REQUIRED OPTIONS
              --PlinkFile=PLINKFILE
                      Your pink filename without the extension. 
                      For example, your plink files are plinkexample.bim, plinkexample.bed, plinkexample.fam. Please use "--PlinkFile=plinkexample"
              --OutputFIle=OUTPUTFILE
		      Your output pink filename without the extension.		
        OPTIONAL OPTIONS
              --mind=MIND
                      PLINK parameter used to control missing genotypes per individual. Default is 0.02.
              --geno=GENO
                      PLINK parameter used to control missing genotypes per SNP. Default is 0.02. 
              --maf=MAF 
                      PLINK parameter used to control allele frequency. Default is 0.05.
              --hwe=HWE
                      PLINK parameter used to control Hardy-Weinberg equilibrium. Default is 0.000001.
EOF
}


# check inputs
options=$@
arguments=${options}

if [ -z "${arguments}" ]; then _usage; exit; fi

for argument in $options
do
        case $argument in
        --PlinkFile=*) PlinkFile=${argument/*=/""} ;;
        --OutputFile=*) OutputFile=${argument/*=/""} ;;
        --geno=*) geno=${argument/*=/""} ;;
        --mind=*) mind=${argument/*=/""} ;;
        --maf=*) maf=${argument/*=/""} ;;
        --hwe=*) hwe=${argument/*=/""} ;;
        --help) _usage; exit;;
        esac
done

# set defaults
if [ -z "${PlinkFile}" ]; then echo "Plink file is not provided"; exit; fi
if [ -z "${OutputFile}" ]; then echo "Output file is not provided"; exit; fi
if [ -z "${mind}" ]; then mind=0.02; fi
if [ -z "${geno}" ]; then geno=0.02; fi
if [ -z "${maf}" ]; then maf=0.05; fi
if [ -z "${hwe}" ]; then hwe=0.000001; fi

# qc based on mind geno maf and hwe
plink1.9 --bfile ${PlinkFile} --noweb --allow-no-sex --mind ${mind} --geno ${geno} --maf ${maf} --hwe ${hwe} --snps-only 'just-acgt' --make-bed --out ${OutputFile} > /dev/null

