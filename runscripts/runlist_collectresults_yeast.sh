#!/bin/bash

pheno=$1

## yeast
interaction=hygeSSI
pvalueCutoff=0.05
minPath=5
samplePerms=10
snpPathwayFile=snp_pathway_min5_max300.mat
snpGeneMappingFile=snpgenemapping_500bp.mat

#for  pheno in SC4NQO01ugml_38h SCCHX05ugml_38h SCpH3_38h SCpH8_38h YPD42_40h YPDCHX05_40h YPDSDS_40h YPGLYCEROL_40h
#do

cd /project/csbio/wwang/BridGE/project_yeast_${pheno}_complex_t25_b50_mhygeSSI
for model in RR DD RD combined
do
#     nice matlab -nodisplay -nodesktop -nosplash -r "fdrsampleperm('${model}','${interaction}',${pvalueCutoff},${minPath},${samplePerms});exit" </dev/null> /dev/null &&
     run_collectresults.sh --projectDir=project_yeast_${pheno}_complex_t25_b50_mhygeSSI --model=${model} --fdrcut=0.4 --snpPathwayFile=${snpPathwayFile} --snpGeneMappingFile=${snpGeneMappingFile} &
done

cd /project/csbio/wwang/BridGE/project_yeast_${pheno}_GI_t25_b50_mhygeSSI 
for model in RR DD RD combined
do
#     nice matlab -nodisplay -nodesktop -nosplash -r "fdrsampleperm('${model}','${interaction}',${pvalueCutoff},${minPath},${samplePerms});exit" </dev/null> /dev/null &&
     run_collectresults.sh --projectDir=project_yeast_${pheno}_GI_t25_b50_mhygeSSI --model=${model} --fdrcut=0.4 --snpPathwayFile=${snpPathwayFile} --snpGeneMappingFile=${snpGeneMappingFile} &
done

cd /project/csbio/wwang/BridGE/project_yeast_${pheno}_complex_t50_b25_mhygeSSI 
for model in RR DD RD combined
do
#     nice matlab -nodisplay -nodesktop -nosplash -r "fdrsampleperm('${model}','${interaction}',${pvalueCutoff},${minPath},${samplePerms});exit" </dev/null> /dev/null &&
     run_collectresults.sh --projectDir=project_yeast_${pheno}_complex_t50_b25_mhygeSSI --model=${model} --fdrcut=0.4 --snpPathwayFile=${snpPathwayFile} --snpGeneMappingFile=${snpGeneMappingFile} &
done

cd /project/csbio/wwang/BridGE/project_yeast_${pheno}_GI_t50_b25_mhygeSSI 
for model in RR DD RD combined
do
#     nice matlab -nodisplay -nodesktop -nosplash -r "fdrsampleperm('${model}','${interaction}',${pvalueCutoff},${minPath},${samplePerms});exit" </dev/null> /dev/null &&
     run_collectresults.sh --projectDir=project_yeast_${pheno}_GI_t50_b25_mhygeSSI --model=${model} --fdrcut=0.4 --snpPathwayFile=${snpPathwayFile} --snpGeneMappingFile=${snpGeneMappingFile} &
done 

cd /project/csbio/wwang/BridGE
#done
