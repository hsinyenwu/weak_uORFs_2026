#!/bin/bash
# Reruns every revision figure module and the compact driver (A3, A2, A5, A6, A7, A8, A9, A11, A12, A31, FigS2, the ladder in both modes, then the driver); ~11 min. Then run finalize_figures.sh.
# V48 (2026-11-17): the figure numbers printed in the titles follow figure_renumbering_v48.json. This script does not run A18, A24, A25_figure or the schematic;
# rerun_for_renumbering_v48.sh lists every script that draws a supplemental figure and redraws them without touching the panels of the main figures.
cd /Users/wu/Desktop/Weak_uORF_submission || exit 1
echo "start $(date)"
for m in A3_uIF_abundance_and_independence A2_leader_splicing_audit A5_riboseq_QC_figure A6_nulls_and_uncertainty A7_ATG_composition_controls A8_CPuORF_matched_controls A9_GO_figures_FDR A11_uncapped_vs_uIF A12_PARE_uncapped_vs_uIF A31_decay_by_inclusion_and_splicing FigS2_TSSseq_panels; do   # V38 (2026-11-11): A31 reads the per-gene tables of A2, A11 and A12, so it runs after them
  echo "== $m $(date)"; Rscript revision_modules/$m.R > output/${m%%_*}_run_compact.log 2>&1; echo "exit $?"
done
echo "== ladder CAGE"; Rscript revision_modules/Fig3C_prevalence_ladder_panel.R > output/ladder_CAGE_compact.log 2>&1; echo "exit $?"
echo "== ladder TSSseq"; LADDER_SOURCE=TSSseq Rscript revision_modules/Fig3C_prevalence_ladder_panel.R > output/ladder_TSSseq_compact.log 2>&1; echo "exit $?"
D=_to_delete/driver_compact_20261028; python3 revision_modules/build_revised_figure_driver.py 20260827_weak_uORFs.R $D/revised_figure_driver.R $D > /dev/null && Rscript $D/revised_figure_driver.R > $D/driver_run.log 2>&1; echo "driver exit $? $(date)"
echo "end $(date)"
