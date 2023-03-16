#!/bin/bash
# Testrun notebooks
for f in ./vignettes/*.Rmd
do if [[ "$f" != *"landing"* ]] && [[ "$f" != *"install"* ]]
then
arrf=(${f//./ })
# Try to execute notebook
jupyter nbconvert --to 'html' --execute ./${arrf[0]}.ipynb
# Clean up
rm ${arrf[0]}.html
fi
done
