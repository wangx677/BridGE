#!/bin/bash

interaction=hygeSSI
pvalueCutoff=0.005
minPath=10
samplePerms=10

# Prostate Cancer
cd /project/csbio/wwang/BridGE/project_ProstateC_812_hygeSSI
for model in RR DD RD combined
do
     nice matlab -nodisplay -nodesktop -nosplash -r "fdrsampleperm('${model}','${interaction}',${pvalueCutoff},${minPath},${samplePerms});exit" </dev/null> /dev/null && 
     run_collectresults.sh --projectDir=project_ProstateC_812_hygeSSI --model=${model} --fdrcut=0.4 --validationDir=project_ProstateC_207_all_hygeSSI &
done

cd /project/csbio/wwang/BridGE/project_ProstateC_812_mhygeSSI
for model in RR DD RD combined
do
     nice matlab -nodisplay -nodesktop -nosplash -r "fdrsampleperm('${model}','${interaction}',${pvalueCutoff},${minPath},${samplePerms});exit" </dev/null> /dev/null &&
     run_collectresults.sh --projectDir=project_ProstateC_812_mhygeSSI --model=${model} --fdrcut=0.4 --validationDir=project_ProstateC_207_all_mhygeSSI &
done

cd /project/csbio/wwang/BridGE/project_ProstateC_207_all_hygeSSI
for model in RR DD RD combined
do  
     nice matlab -nodisplay -nodesktop -nosplash -r "fdrsampleperm('${model}','${interaction}',${pvalueCutoff},${minPath},${samplePerms});exit" </dev/null> /dev/null &&
     run_collectresults.sh --projectDir=project_ProstateC_207_all_hygeSSI --model=${model} --fdrcut=0.4 --validationDir=project_ProstateC_812_hygeSSI &
done

cd /project/csbio/wwang/BridGE/project_ProstateC_207_all_mhygeSSI
for model in RR DD RD combined
do
     nice matlab -nodisplay -nodesktop -nosplash -r "fdrsampleperm('${model}','${interaction}',${pvalueCutoff},${minPath},${samplePerms});exit" </dev/null> /dev/null &&
     run_collectresults.sh --projectDir=project_ProstateC_207_all_mhygeSSI --model=${model} --fdrcut=0.4 --validationDir=project_ProstateC_812_mhygeSSI &
done

# Parkinson's Disease
cd /project/csbio/wwang/BridGE/project_PD_Simon_hygeSSI
for model in RR DD RD combined
do
     nice matlab -nodisplay -nodesktop -nosplash -r "fdrsampleperm('${model}','${interaction}',${pvalueCutoff},${minPath},${samplePerms});exit" </dev/null> /dev/null &&     
     run_collectresults.sh --projectDir=project_PD_Simon_hygeSSI --model=${model} --fdrcut=0.4 --validationDir=project_PD_NGRC_hygeSSI &
done

cd /project/csbio/wwang/BridGE/project_PD_Simon_mhygeSSI
for model in RR DD RD combined
do
     nice matlab -nodisplay -nodesktop -nosplash -r "fdrsampleperm('${model}','${interaction}',${pvalueCutoff},${minPath},${samplePerms});exit" </dev/null> /dev/null &&
     run_collectresults.sh --projectDir=project_PD_Simon_mhygeSSI --model=${model} --fdrcut=0.4 --validationDir=project_PD_NGRC_mhygeSSI &
done

cd /project/csbio/wwang/BridGE/project_PD_NGRC_hygeSSI
for model in RR DD RD combined
do
     nice matlab -nodisplay -nodesktop -nosplash -r "fdrsampleperm('${model}','${interaction}',${pvalueCutoff},${minPath},${samplePerms});exit" </dev/null> /dev/null && 
     run_collectresults.sh --projectDir=project_PD_NGRC_hygeSSI --model=${model} --fdrcut=0.4 --validationDir=project_PD_Simon_hygeSSI &
done

cd /project/csbio/wwang/BridGE/project_PD_NGRC_mhygeSSI
for model in RR DD RD combined
do
     nice matlab -nodisplay -nodesktop -nosplash -r "fdrsampleperm('${model}','${interaction}',${pvalueCutoff},${minPath},${samplePerms});exit" </dev/null> /dev/null &&
     run_collectresults.sh --projectDir=project_PD_NGRC_mhygeSSI --model=${model} --fdrcut=0.4 --validationDir=project_PD_Simon_mhygeSSI &
done

interaction=regSSI
cd /project/csbio/wwang/BridGE/project_PD_Simon_LR
for model in RR DD RD combined
do
     nice matlab -nodisplay -nodesktop -nosplash -r "fdrsampleperm('${model}','${interaction}',${pvalueCutoff},${minPath},${samplePerms});exit" </dev/null> /dev/null &&
     run_collectresults.sh --projectDir=project_PD_Simon_LR --model=${model} --fdrcut=0.4 --validationDir=project_PD_NGRC_LR --interaction=regSSI &
done

cd /project/csbio/wwang/BridGE/project_PD_NGRC_LR
for model in RR DD RD combined
do
     nice matlab -nodisplay -nodesktop -nosplash -r "fdrsampleperm('${model}','${interaction}',${pvalueCutoff},${minPath},${samplePerms});exit" </dev/null> /dev/null &&
     run_collectresults.sh --projectDir=project_PD_NGRC_LR --model=${model} --fdrcut=0.4 --validationDir=project_PD_Simon_LR --interaction=regSSI &
done

# Autism
interaction=hygeSSI
cd /project/csbio/wwang/BridGE/project_MSSNG_hiseqx_CP_mhygeSSI
for model in RR DD RD combined
do
     nice matlab -nodisplay -nodesktop -nosplash -r "fdrsampleperm('${model}','${interaction}',${pvalueCutoff},${minPath},${samplePerms});exit" </dev/null> /dev/null &&
     run_collectresults.sh --projectDir=project_MSSNG_hiseqx_CP_mhygeSSI --model=${model} --fdrcut=0.4 &
done

cd /project/csbio/wwang/BridGE/project_MSSNG_hiseqx_hallmark_mhygeSSI
for model in RR DD RD combined
do
     nice matlab -nodisplay -nodesktop -nosplash -r "fdrsampleperm('${model}','${interaction}',${pvalueCutoff},${minPath},${samplePerms});exit" </dev/null> /dev/null &&
     run_collectresults.sh --projectDir=project_MSSNG_hiseqx_hallmark_mhygeSSI --model=${model} --fdrcut=0.4 &
done

# Diabetes
cd  /project/csbio/wwang/BridGE/project_MedGenome_Diabetes_MF_hygeSSI
for model in RR DD RD combined
do
     nice matlab -nodisplay -nodesktop -nosplash -r "fdrsampleperm('${model}','${interaction}',${pvalueCutoff},${minPath},${samplePerms});exit" </dev/null> /dev/null &&
     run_collectresults.sh --projectDir=project_MedGenome_Diabetes_MF_hygeSSI --model=${model} --fdrcut=0.4 &
done

cd  /project/csbio/wwang/BridGE/project_MedGenome_Diabetes_MF_mhygeSSI
for model in RR DD RD combined
do
     nice matlab -nodisplay -nodesktop -nosplash -r "fdrsampleperm('${model}','${interaction}',${pvalueCutoff},${minPath},${samplePerms});exit" </dev/null> /dev/null &&
     run_collectresults.sh --projectDir=project_MedGenome_Diabetes_MF_mhygeSSI --model=${model} --fdrcut=0.4 &
done

cd  /project/csbio/wwang/BridGE/project_MedGenome_Diabetes_hygeSSI
for model in RR DD RD combined
do
     nice matlab -nodisplay -nodesktop -nosplash -r "fdrsampleperm('${model}','${interaction}',${pvalueCutoff},${minPath},${samplePerms});exit" </dev/null> /dev/null &&
     run_collectresults.sh --projectDir=project_MedGenome_Diabetes_hygeSSI --model=${model} --fdrcut=0.4 &
done

cd  /project/csbio/wwang/BridGE/project_MedGenome_Diabetes_mhygeSSI
for model in RR DD RD combined
do
     nice matlab -nodisplay -nodesktop -nosplash -r "fdrsampleperm('${model}','${interaction}',${pvalueCutoff},${minPath},${samplePerms});exit" </dev/null> /dev/null &&
     run_collectresults.sh --projectDir=project_MedGenome_Diabetes_mhygeSSI --model=${model} --fdrcut=0.4 &
done

