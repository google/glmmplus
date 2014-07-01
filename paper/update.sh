#!/bin/bash

R CMD Sweave main.Rnw  # creates main.tex
# The BibTeX sandwhich
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex

rm main.tex  # So I don't open it instead of main.Rnw

# Copy pdf to Drive
cp main.pdf ~/Google\ Drive/Data_Science/glmmplus
