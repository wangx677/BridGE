#!/bin/bash
for dir in `ls -d project*yeast*complex*`; 
do 
     if [ "$(ls -A /project/csbio/wwang/BridGE/$dir/BPM_WPM_info)" ];then 
          cd /project/csbio/wwang/BridGE/$dir/BPM_WPM_info 
          for file in `ls *.xls`; 
          do 
               nice matlab -nodisplay -nodesktop -nosplash -r "yeast_flip_positive_negative_info_xls('$file');exit" </dev/null> /dev/null; 
          done; 
     fi & 
done
