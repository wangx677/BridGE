setenv CURRENTDIR `pwd`
setenv BRIDGEPATH ${CURRENTDIR} 
setenv PATH ${BRIDGEPATH}/testfuns:${BRIDGEPATH}/testfuns/bin:${BRIDGEPATH}:${BRIDGEPATH}/scripts:${BRIDGEPATH}/runscripts:${BRIDGEPATH}/datatools:${BRIDGEPATH}/cassi:${BRIDGEPATH}/PPIExtend:${PATH}:${BRIDGEPATH}/scripts_yeast
setenv MATLABPATH ${BRIDGEPATH}/datatools:${BRIDGEPATH}/corefuns:${BRIDGEPATH}/analysis_tool:${BRIDGEPATH}/miscs:${BRIDGEPATH}/runscripts:${BRIDGEPATH}/scripts_yeast:${BRIDGEPATH}/scripts_evaluation:${BRIDGEPATH}/scripts_PD_pbody
alias plink ${BRIDGEPATH}/scripts/plink1.9
