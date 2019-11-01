#!/bin/csh -f

set filename = `echo $1 | sed 's/.tex//'`

pdflatex $filename
bibtex $filename
~/Dropbox/Python/useful_scripts/ads_importer.py $filename
bibtex $filename
pdflatex $filename
pdflatex $filename

open $filename.pdf
