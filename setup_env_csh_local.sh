setenv CURRENTDIR `pwd`
setenv BRIDGEPATH ${CURRENTDIR} 
setenv PATH ${BRIDGEPATH}/testfuns:${BRIDGEPATH}/testfuns/bin:${BRIDGEPATH}:${BRIDGEPATH}/scripts:${BRIDGEPATH}/runscripts:${BRIDGEPATH}/datatools:${BRIDGEPATH}/cassi:${BRIDGEPATH}/PPIExtend:${PATH}:${BRIDGEPATH}/script_yeast
setenv MATLABPATH ${BRIDGEPATH}/datatools:${BRIDGEPATH}/corefuns:${BRIDGEPATH}/analysis_tool:${BRIDGEPATH}/miscs:${BRIDGEPATH}/runscripts:${BRIDGEPATH}/script_yeast
alias plink ${BRIDGEPATH}/scripts/plink1.9
