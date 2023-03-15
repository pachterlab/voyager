#!/bin/bash
for f in ./vignettes/*.Rmd
do 
jupytext --to notebook $f
done
