export CURRENTDIR=`pwd`
export BRIDGEPATH=${CURRENTDIR} 

export PATH=${BRIDGEPATH}/testfuns:${BRIDGEPATH}/testfuns/bin:${BRIDGEPATH}:${BRIDGEPATH}/scripts_local:${BRIDGEPATH}/runscripts_local:${BRIDGEPATH}/datatools_local:${BRIDGEPATH}/cassi:${BRIDGEPATH}/PPIExtend:${BRIDGEPATH}/scripts_data_processing:/project/csbio/wwang/BridGE/scripts_ALS:${BRIDGEPATH}/scripts_yeast:${BRIDGEPATH}/scripts_NDDs:${BRIDGEPATH}/BridGE_genes:${PATH}

export MATLABPATH=${BRIDGEPATH}/testfuns:${BRIDGEPATH}/datatools_local:${BRIDGEPATH}/corefuns_local:${BRIDGEPATH}/analysis_tool_local:${BRIDGEPATH}/miscs_local:${BRIDGEPATH}/runscripts_local:${BRIDGEPATH}/scripts_yeast:${BRIDGEPATH}/scripts_evaluation:${BRIDGEPATH}/scripts_pbody:${BRIDGEPATH}/BridGE_genes:${BRIDGEPATH}/scripts_NDDs

alias plink=${BRIDGEPATH}/scripts/plink1.9
