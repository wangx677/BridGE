#!/bin/bash

for  pheno in SC4NQO01ugml_38h SCCHX05ugml_38h SCpH3_38h SCpH8_38h YPD42_40h YPDCHX05_40h YPDSDS_40h YPGLYCEROL_40h; 
do 
     for binary in t25_b50 t50_b25 t25_b25; 
     do 
     cd /project/csbio/wwang/BridGE/project_yeast_${pheno}_complex_${binary}_mhygeSSI; 
          for model in RR DD RD combined; 
          do  
               nice matlab -nodisplay -nodesktop -nosplash -r "add_yeast_GI_enrich('output_results_ssM_hygeSSI_alpha10.05_alpha20.05_${model}_R0.mat.xls');exit" </dev/null> /dev/null; 
          done; 
     done; 
done 
