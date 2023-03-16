#!/bin/bash
# Testrun notebooks
for f in ./vignettes/*.Rmd
do if [[ "$f" != *"landing"* ]] && [[ "$f" != *"install"* ]]
then
arrf=(${f//./ })
# Try to execute notebook
jupyter nbconvert --to html --ExecutePreprocessor.kernel_name=ir --execute .${arrf[0]}.ipynb
rm .${arrf[0]}.html
fi
done
