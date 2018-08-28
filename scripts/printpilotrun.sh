#!/bin/bash

# this script collects all results from pilotrun and summarize them

nperms=$1
outputfile=pretmp.txt

# if file exist, remove 
if [ -f pre_result_snpP${nperms}.txt ];then
     rm pre_result_snpP${nperms}.txt
     rm pre_result_snpP${nperms}_summary.txt
fi

# collect results from different disease model and density combinations
echo -e filename "\t"BPM_minpv_no"\t"BPM_minfdr"\t"WPM_minpv_no"\t"WPM_minfdr"\t"PATH_minpv_no"\t"PATH_minfdr > pre_result_snpP${nperms}.txt
for file in `ls results*R0*snpP*run1*total${nperms}.mat`
do
     nice matlab -nodisplay -nodesktop -r "printpilotrun('${file}',${nperms},'${outputfile}');exit" </dev/null> /dev/null
     while read p; do
          echo -e $file "\t" $p >> pre_result_snpP${nperms}.txt
     done<${outputfile}
     rm ${outputfile}
done

echo -e interaction"\t"density"\t"BPM_minpv_no"\t"BPM_minfdr"\t"WPM_minpv_no"\t"WPM_minfdr"\t"PATH_minpv_no"\t"PATH_minfdr > pre_result_snpP${nperms}_summary.txt
  
for interaction in RR DD combined lr
do
     echo `grep ${interaction} pre_result_snpP${nperms}.txt | awk -Fdensity '{print $2}' | awk -F_ '{print $1}'| sort |uniq ` >> ${outputfile}
     for density in `cat ${outputfile}`
     do
          grep $interaction pre_result_snpP${nperms}.txt | grep density${density}_> tmp.txt
          while read line
             do
                     tmp1=`echo $line |awk '{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7}'`
             done<tmp.txt
          echo ${interaction} ${density} $tmp1 >>pre_result_snpP${nperms}_summary.txt
          rm tmp.txt
     done
     rm ${outputfile}
done

rm pre_result_snpP${nperms}.txt
