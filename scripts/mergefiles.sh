#!/bin/bash

# this script merges plink data files based on the given file list

FileList=$1 # plink file list, each line is one plink data
OutputFile=$2 # plink file name of merged data

if [ -z "${FileList}" ]; then echo "Merge list file is not provided"; exit; fi
if [ -z "${OutputFile}" ]; then echo "Output Plink file is not provided"; exit; fi	

PlinkFile=`head -1 ${FileList}`

tail -n +2 ${FileList} > FileList.tmp

if [ -f MergeList.tmp ]
then
	rm MergeList.tmp
fi

for file in `cat FileList.tmp`
do
	echo ${file}.bed ${file}.bim ${file}.fam >> MergeList.tmp
done

plink1.9 --bfile ${PlinkFile} --allow-no-sex --merge-list MergeList.tmp --mind 1 --geno 1 --make-bed  --out ${OutputFile} > /dev/null

if [ -f missnp.tmp ]
then
	rm missnp.tmp
fi

while [ ! -f ${OutputFile}.bim ];
do
        if [ -f ${OutputFile}-merge.missnp ];then
		cat ${OutputFile}-merge.missnp >> missnp.tmp
		grep Multiple ${OutputFile}.log |awk -F\' '{print $2}' >> missnp.tmp

                plink1.9 --bfile ${PlinkFile} --exclude missnp.tmp --geno 1 --mind 1 --make-bed --out plink_tmp > /dev/null
		mm=1
		while read line
		do
			plink1.9 --bfile ${line} --exclude missnp.tmp --geno 1 --mind 1 --make-bed --out plink_tmp_$mm > /dev/null
			echo plink_tmp_${mm}.bed plink_tmp_${mm}.bim plink_tmp_${mm}.fam >> MergeList_new.tmp
			mm=`expr $mm + 1`
		done<FileList.tmp	

		rm ${OutputFile}-merge.missnp

                plink1.9 --bfile plink_tmp --allow-no-sex  --merge-list MergeList_new.tmp --geno 1 --mind 1 --make-bed  --out ${OutputFile} > /dev/null
		rm MergeList_new.tmp
        fi
done


rm MergeList.tmp FileList.tmp 
