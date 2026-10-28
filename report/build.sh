#!/bin/bash
set -e
pandoc body.md -o report.tex --standalone \
  --include-in-header=preamble.tex --include-before-body=titlepage.tex \
  -V documentclass=article -V papersize=a4 -V geometry:margin=2.5cm \
  -V fontsize=11pt -V fontfamily=mathptmx --top-level-division=section
sed -i 's/\\\^{}/\\textasciicircum{}/g' report.tex
latexmk -pdf -interaction=nonstopmode report.tex
