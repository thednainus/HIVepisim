#!/bin/bash

#R_wrapper_script=wrapper_script.R
R_script=test_exit_code.R
#mkdir -p $dirname
#cd $dirname

Rscript $R_script $line_number

#to capture error insite Rscript
error=$?
echo $error

if [[ $error -eq 85 ]]
then
 echo $error
 exit 85
elif [[ $error -eq 0 ]]
then 
 mkdir -p results
 echo $error
 rm *.rdata
 
fi
