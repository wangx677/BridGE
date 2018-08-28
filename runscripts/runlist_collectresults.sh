#!/bin/bash
# Prostate Cancer
cd /project/csbio/wwang/BridGE/project_ProstateC_812_hygeSSI
run_collectresults.sh --projectDir=project_ProstateC_812_hygeSSI --fdrcut=0.4 --validationDir=project_ProstateC_207_all_hygeSSI

cd /project/csbio/wwang/BridGE/project_ProstateC_812_mhygeSSI
run_collectresults.sh --projectDir=project_ProstateC_812_mhygeSSI --fdrcut=0.4 --validationDir=project_ProstateC_207_all_mhygeSSI

cd /project/csbio/wwang/BridGE/project_ProstateC_207_all_hygeSSI
run_collectresults.sh --projectDir=project_ProstateC_207_all_hygeSSI --fdrcut=0.4 --validationDir=project_ProstateC_812_hygeSSI

cd /project/csbio/wwang/BridGE/project_ProstateC_207_all_mhygeSSI
run_collectresults.sh --projectDir=project_ProstateC_207_all_mhygeSSI --fdrcut=0.4 --validationDir=project_ProstateC_812_mhygeSSI

# Parkinson's Disease
cd /project/csbio/wwang/BridGE/project_PD_Simon_hygeSSI
run_collectresults.sh --projectDir=project_PD_Simon_hygeSSI --fdrcut=0.4 --validationDir=project_PD_NGRC_hygeSSI

cd /project/csbio/wwang/BridGE/project_PD_Simon_mhygeSSI
run_collectresults.sh --projectDir=project_PD_Simon_mhygeSSI --fdrcut=0.4 --validationDir=project_PD_NGRC_mhygeSSI

cd /project/csbio/wwang/BridGE/project_PD_Simon_LR
run_collectresults.sh --projectDir=project_PD_Simon_LR --fdrcut=0.4 --validationDir=project_PD_NGRC_LR --interaction=regSSI

cd /project/csbio/wwang/BridGE/project_PD_NGRC_hygeSSI
run_collectresults.sh --projectDir=project_PD_NGRC_hygeSSI --fdrcut=0.4 --validationDir=project_PD_Simon_hygeSSI

cd /project/csbio/wwang/BridGE/project_PD_NGRC_mhygeSSI
run_collectresults.sh --projectDir=project_PD_NGRC_mhygeSSI --fdrcut=0.4 --validationDir=project_PD_Simon_mhygeSSI

cd /project/csbio/wwang/BridGE/project_PD_NGRC_LR
run_collectresults.sh --projectDir=project_PD_NGRC_LR --fdrcut=0.4 --validationDir=project_PD_Simon_LR --interaction=regSSI

# Autism
cd /project/csbio/wwang/BridGE/project_MSSNG_hiseqx_CP_mhygeSSI
run_collectresults.sh --projectDir=project_MSSNG_hiseqx_CP_mhygeSSI --fdrcut=0.4 

cd /project/csbio/wwang/BridGE/project_MSSNG_hiseqx_hallmark_mhygeSSI
run_collectresults.sh --projectDir=project_MSSNG_hiseqx_hallmark_mhygeSSI --fdrcut=0.4

# Diabetes
cd  /project/csbio/wwang/BridGE/project_MedGenome_Diabetes_MF_hygeSSI
run_collectresults.sh --projectDir=project_MedGenome_Diabetes_MF_hygeSSI --fdrcut=0.4

cd  /project/csbio/wwang/BridGE/project_MedGenome_Diabetes_MF_mhygeSSI
run_collectresults.sh --projectDir=project_MedGenome_Diabetes_MF_mhygeSSI --fdrcut=0.4

cd  /project/csbio/wwang/BridGE/project_MedGenome_Diabetes_hygeSSI
run_collectresults.sh --projectDir=project_MedGenome_Diabetes_hygeSSI --fdrcut=0.4

cd  /project/csbio/wwang/BridGE/project_MedGenome_Diabetes_mhygeSSI
run_collectresults.sh --projectDir=project_MedGenome_Diabetes_mhygeSSI --fdrcut=0.4

# yeast


