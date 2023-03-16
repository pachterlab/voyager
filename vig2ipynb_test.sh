#!/bin/bash
# Testrun notebooks
for f in ./vignettes/*.Rmd
do if [[ "$f" != *"landing"* ]] && [[ "$f" != *"install"* ]]
then
arrf=(${f//./ })
# Try to execute notebook
jupyter nbconvert --execute .${arrf[0]}.ipynb
fi
done
