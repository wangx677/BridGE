setenv CURRENTDIR `pwd`
setenv BRIDGEPATH ${CURRENTDIR} 
setenv PATH ${BRIDGEPATH}:${BRIDGEPATH}/scripts:${BRIDGEPATH}/runscripts:${BRIDGEPATH}/datatools:${BRIDGEPATH}/cassi:${BRIDGEPATH}/PPIExtend:${BRIDGEPATH}/BridGE_genes:${PATH}
setenv MATLABPATH ${BRIDGEPATH}/datatools:${BRIDGEPATH}/corefuns:${BRIDGEPATH}/analysis_tool:${BRIDGEPATH}/miscs:${BRIDGEPATH}/runscripts:${BRIDGEPATH}/BridGE_genes
alias plink ${BRIDGEPATH}/scripts/plink1.9
