#!/bin/bash
# Compile TeX document to PDF
cd /mnt/akita/nagasaku/WORK/docs
pdflatex von_mises_block_newton_algorithm.tex
pdflatex von_mises_block_newton_algorithm.tex  # Run twice for references
echo "PDF generated: von_mises_block_newton_algorithm.pdf"