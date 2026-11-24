#!/bin/bash
# V48 (2026-11-17): 21 supplemental figures and 12 supplemental tables, numbered in the order of first citation
# (the number is printed in the title inside each figure and in the README or legend of each table; revision_docx/renumber_v48_modules.py changed the strings).
# This script redraws these figures and rewrites these tables. No analysis changed. Then run: FIGS_ONLY=1 bash revision_modules/finalize_figures.sh.
# The figure driver is not rerun into its usual folder: the two pages drawn from saved objects (S9, S10) are redrawn by fig_S10_S12_compact.R, and the
# logo page (S11) is taken from a driver run in a separate folder and copied into the driver folder.
# Order: A31 reads the per-gene tables of A2, A11 and A12; A25 (full run) reads the table of A2 and writes Supplemental Table S5; A20 writes the legends
# of Supplemental Tables S1, S3, S5 and S6 last.
cd /Users/wu/Desktop/Weak_uORF_submission || exit 1
echo "start $(date)"
for m in A2_leader_splicing_audit A3_uIF_abundance_and_independence A6_nulls_and_uncertainty A7_ATG_composition_controls A8_CPuORF_matched_controls A11_uncapped_vs_uIF A12_PARE_uncapped_vs_uIF A31_decay_by_inclusion_and_splicing A18_mRNA_weighted_TSS_splicing A24_GO_by_max_pUP_figures_table FigS2_TSSseq_panels Fig_schematic_leader_splicing_audit; do
  echo "== $m $(date)"; Rscript revision_modules/$m.R > output/${m%%_*}_run_v48.log 2>&1; echo "exit $?"
done
echo "== ladder TSSseq $(date)"; LADDER_SOURCE=TSSseq Rscript revision_modules/Fig3C_prevalence_ladder_panel.R > output/ladder_TSSseq_v48.log 2>&1; echo "exit $?"
echo "== A25 (full run) $(date)"; Rscript revision_modules/A25_ouORF_leader_splicing_audit.R > output/A25_run_v48.log 2>&1; echo "exit $?"
echo "== A20 $(date)"; Rscript revision_modules/A20_flag_neighbouring_gene_CDS_and_table_legends.R > output/A20_run_v48.log 2>&1; echo "exit $?"
D=_to_delete/driver_compact_20261028; mkdir -p $D/_before_V48_20261117
cp $D/FigureS11_revised_multi_uORF_genes.pdf $D/FigureS12_revised_undetected_uORFs.pdf $D/FigureS13_compact_logos_50pct.pdf $D/_before_V48_20261117/
echo "== S9, S10 from the saved driver objects $(date)"; Rscript revision_modules/fig_S10_S12_compact.R $D > output/fig_S10_S12_v48.log 2>&1; echo "exit $?"
D2=_to_delete/driver_V48_logos_20261117; mkdir -p $D2
echo "== driver in a separate folder (logo page only) $(date)"; python3 revision_modules/build_revised_figure_driver.py 20260827_weak_uORFs.R $D2/revised_figure_driver.R $D2 > /dev/null && Rscript $D2/revised_figure_driver.R > $D2/driver_run.log 2>&1; echo "driver exit $?"
cp $D2/FigureS13_compact_logos_50pct.pdf $D/FigureS13_compact_logos_50pct.pdf && echo "logo page copied into $D"
echo "end $(date)"
