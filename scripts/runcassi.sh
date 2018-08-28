#!/bin/bash

# this script is used to compute additive-additive interaction using logistic regression

data=$1 # data in plink format without suffix
interaction=$2 # define interaction type for CASSI (http://www.staff.ncl.ac.uk/richard.howey/cassi)
pvcut=$3 # CASSI only print out interactions whose p-value is less than pvcut
R=$4 # R=0 for real phenotype label, otherwise for random labels

if [ "${interaction}" = "LR" ]; then
	# Logistic Regression
	if [ ! -f cassi_pv${pvcut}_R${R}.lr ]; then
	     cassi -i ${data}.bed -mem2 -lr -lr-th ${pvcut} -max 0 -o cassi_pv${pvcut}_R${R}.lr
	fi
	matlab -nodisplay -nodesktop -r "addpath ${BRIDGEPATH}/corefuns;cassissm('cassi_pv${pvcut}_R${R}.lr');exit"  </dev/null> /dev/null
fi
