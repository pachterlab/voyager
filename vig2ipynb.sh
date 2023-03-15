#!/bin/bash
for f in ./vignettes/*.Rmd
do if [[ "$f" != *"landing"* ]] && [[ "$f" != *"install"* ]]
then
jupytext --to notebook $f
fi
done
