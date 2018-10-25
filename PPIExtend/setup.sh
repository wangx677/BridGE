#!/bin/sh

# Setup PPIExtend.sh 


if test '!' -f `pwd`/PPIE.jar ; then
  echo 'You must run this setup script in the directory containing PPIE.jar'
  exit 1
fi

sed <PPIExtend.in >PPIExtend -e 's;PPIE_JAR_DIR=.;PPIE_JAR_DIR='`pwd`';'

chmod +x PPIExtend

echo 'The shell script PPIExtend is configured.' 
echo 'Please move the shell script to a /bin directory and leave' 
echo 'all other files in this directory.'
echo
echo 'To run PPIExtend, use'
echo 'PPIExtend <input file name> <delta> <tau> <alpha> <result file name>'
echo
