setenv CURRENTDIR `pwd`
setenv BRIDGEPATH ${CURRENTDIR} 
setenv PATH ${BRIDGEPATH}/testfuns:${BRIDGEPATH}/testfuns/bin:${BRIDGEPATH}:${BRIDGEPATH}/scripts:${BRIDGEPATH}/runscripts:${BRIDGEPATH}/datatools:${BRIDGEPATH}/cassi:${PATH}
setenv MATLABPATH ${BRIDGEPATH}/datatools:${BRIDGEPATH}/corefuns:${BRIDGEPATH}/miscs:${BRIDGEPATH}/runscripts
alias plink ${BRIDGEPATH}/scripts/plink1.9
