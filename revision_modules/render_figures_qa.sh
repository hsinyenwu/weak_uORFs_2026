#!/bin/bash
# Render every final figure PDF (and the standalone ladder plot/table) at 110 dpi for inspection into output_revised_figs_20261028/_renders_qa/.
F=/Users/wu/Desktop/Weak_uORF_submission/output_revised_figs_20261028
Q=/Users/wu/Desktop/Weak_uORF_submission/output_revised_figs_20261028/_renders_qa
rm -rf "$Q"; mkdir -p "$Q"
for n in 1 2 3 4; do pdftoppm -r 110 -png -singlefile $F/Figure${n}_compact_composite.pdf $Q/Figure$n; done
pdftoppm -r 110 -png -singlefile $F/Figure5_compact_panels_A-F.pdf $Q/Figure5
pdftoppm -r 110 -png -singlefile $F/Figure6_compact_panels_A-D.pdf $Q/Figure6
pdftoppm -r 110 -png -singlefile $F/Figure3_compact_panels_A-B.pdf $Q/Figure3_AB
for s in 1 2 3 4 5 6 7 8 9 11 13 14 15 16 17 18 19; do f=$(ls $F/FigureS${s}_revised_*.pdf | grep -v "_plot\|_table" | head -1); pdftoppm -r 110 -png -singlefile "$f" $Q/S$s; done
ls $Q | wc -l
