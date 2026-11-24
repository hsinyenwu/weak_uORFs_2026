#!/bin/bash
# After the second-pass reruns: final figure files, PNG renders, docx with embedded figures, copies, package.
# V39 (2026-11-11, Larry: 'Renumber'): the decay figures are cited first in the last Results section, so they are numbered last: S4 -> S21 (A11), S5 -> S22 (A12), S6 ... S22 -> S4 ... S20 (S1-S3 and S23 unchanged).
#   Only the final names on the left of each pair below changed. The two pages drawn by fig_S10_S12_compact.R and the logo page keep their old numbers in the driver folder (FigureS11_..., FigureS12_..., FigureS13_compact_logos...)
#   and are copied to FigureS9_, FigureS10_ and FigureS11_. S13 (the TSS-seq ladder) is written by Fig3C_prevalence_ladder_panel.R itself. The V38 files are in output_revised_figs_20261028/_renumbered_20261111/.
# V38 (2026-11-11, Larry: the decay analyses as a Results section of their own): S23 = mRNA half-life and uncapped mRNA by uORF inclusion and leader splicing (A31); S4 and S5 redrawn with the colours of Figures S2B and 2F (values unchanged).
# V21 (2026-11-07, Larry): S8 = schematic of the leader-splicing audit (Fig_schematic_leader_splicing_audit.R); the former S8-S21 are S9-S22 (old files in _renumbered_20261107). FIGS_ONLY=1 stops after the PNG map (the docx chain then runs from finalize_v21.sh).
# V21 (2026-11-06): S16 = overlapping uORFs (A25), GO S18-S20, CPuORF S21.
# V19 (2026-11-05): S16-S18 = the gene-level Gene Ontology grouping (A24).
# V48 (2026-11-17): 21 supplemental figures (final names below; the names in the driver folder $D do not change).
#   S13 (the TSS-seq ladder) is written by Fig3C_prevalence_ladder_panel.R itself. The V47 files are in output_revised_figs_20261028/_replaced_by_V48_renumbering_20261117/.
# V46 (2026-11-14; Larry: "group all ORFs in the first section"): S16 -> S3 (overlapping uORFs, now in the first Results section), S3..S15 -> S4..S16 (final names below).
#   S15 (the TSS-seq ladder) is written by Fig3C_prevalence_ladder_panel.R itself. The V45 files of S3-S16 are in output_revised_figs_20261028/_replaced_by_V46_renumbering_20261114/.
# V16 (2026-11-04): Figure 1B = the uIF distribution (Figure1_compact_B_dist.pdf), S2 = the uIF-TE panels (FigureS2_revised_uIF_vs_TE.pdf).
# V15 (2026-11-04): S3 = the inclusion-adjusted uIF only (A3_Figure_uIF_inclusion_adjusted.pdf).
# V14 (2026-11-03): supplemental numbering S1-S19 (the mORF figure removed; S12-S20 -> S11-S19); the docx steps below build the version named in N.
# Usage: bash revision_modules/finalize_figures.sh  (after rerun_figure_modules.sh). Copies the module figures to their final names, builds the Figure 1-4 composites,
# renders the 200-dpi PNGs for the docx, rebuilds the docx (apply -> embed -> zip -> validate), the marked/clean copies and the submission package.
set -e
cd /Users/wu/Desktop/Weak_uORF_submission
F=output_revised_figs_20261028; D=_to_delete/driver_compact_20261028; OLD=$F/_replaced_by_compact_20261028; mkdir -p $OLD $F/compact_panels $F/png_for_docx
# 1) compact module figures replace the previous full-size files (same final names)
PAIRS="FigureS1_revised_riboseq_QC.pdf|output/A5_Figure_riboseq_QC.pdf FigureS4_revised_uIF_inclusion_adjusted.pdf|output/A3_Figure_uIF_inclusion_adjusted.pdf
FigureS19_revised_uncapped_halflife.pdf|output/A11_Figure_uncapped_vs_uIF.pdf
FigureS20_revised_PARE_uncapped.pdf|output/A12_Figure_PARE_uncapped_vs_uIF.pdf FigureS3_revised_nulls_and_uncertainty.pdf|output/A6_Figure_nulls_and_uncertainty.pdf
FigureS6_revised_leader_splicing_schematic.pdf|output/Fig_schematic_leader_splicing_audit.pdf FigureS7_revised_leader_splicing_audit.pdf|output/A2_Figure_leader_splicing_audit.pdf FigureS8_revised_panels_A-E_TSSseq.pdf|output/FigureS2_TSSseq_panels_A-E.pdf
FigureS12_revised_ATG_composition_controls.pdf|output/A7_Figure_ATG_composition_controls.pdf FigureS15_revised_GO_MF_FDR.pdf|output/A24_GO_MF_by_max_pUP_FDR.pdf
FigureS16_revised_GO_BP_FDR.pdf|output/A24_GO_BP_by_max_pUP_FDR.pdf FigureS17_revised_GO_CC_FDR.pdf|output/A24_GO_CC_by_max_pUP_FDR.pdf
FigureS18_revised_CPuORF_matched_controls.pdf|output/A8_Figure_CPuORF_matched_controls.pdf FigureS5_revised_ouORF_TSS_splicing.pdf|output/A25_Figure_ouORF_vs_TuORF.pdf FigureS14_revised_mRNA_weighted_TSS_splicing.pdf|output/A18_Figure_mRNA_weighted_TSS_splicing.pdf FigureS21_revised_decay_inclusion_splicing.pdf|output/A31_Figure_decay_by_inclusion_and_splicing.pdf
FigureS2_revised_uIF_vs_TE.pdf|$D/FigureS2_revised_uIF_vs_TE.pdf FigureS9_revised_multi_uORF_genes.pdf|$D/FigureS11_revised_multi_uORF_genes.pdf FigureS10_revised_undetected_uORFs.pdf|$D/FigureS12_revised_undetected_uORFs.pdf
FigureS11_revised_TSS_logos_50pct.pdf|$D/FigureS13_compact_logos_50pct.pdf Figure3_compact_panels_A-B.pdf|$D/Figure3_compact_panels_A-B.pdf
Figure5_compact_panels_A-F.pdf|$D/Figure5_compact_panels_A-F.pdf Figure6_compact_panels_A-D.pdf|$D/Figure6_compact_panels_A-D.pdf"
for pr in $PAIRS; do dst=${pr%%|*}; src=${pr##*|}; [ -f "$F/$dst" ] && mv "$F/$dst" "$OLD/"; cp "$src" "$F/$dst"; done
cp $D/Figure1_compact_B_dist.pdf $D/Figure2_compact_panels_B-F.pdf $D/Figure4_compact_panels_D-F.pdf $F/compact_panels/   # V30 (2026-11-09): Figures 2 B-F and 4 D-F are aligned grids (fig_main_layout.R)
# 2) composites (Figures 1-4) from the compact panels + submitted crops
Rscript revision_modules/compose_main_figures.R $D 2>&1 | grep "saved\|bottom\|row at"
python3 revision_modules/compose_main_figures_vector.py $D 2>&1 | grep "saved"   # V21 (2026-11-07): vector composites for the package (Illustrator-editable)
# 3) PNG renders for the docx (200 dpi)
P=$F/png_for_docx; rm -f $P/*.png
for n in 1 2 3 4; do pdftoppm -r 200 -png -singlefile $F/Figure${n}_compact_composite.pdf $P/Figure$n; done
pdftoppm -r 200 -png -singlefile $F/Figure5_compact_panels_A-F.pdf $P/Figure5; pdftoppm -r 200 -png -singlefile $F/Figure6_compact_panels_A-D.pdf $P/Figure6; pdftoppm -r 200 -png -singlefile $F/Figure7_revised_model.pdf $P/Figure7
for s in $(seq 1 22); do f=$(ls $F/FigureS${s}_revised_*.pdf | grep -v "_plot\|_table" | head -1); pdftoppm -r 200 -png -singlefile "$f" $P/S$s; done   # S13 = the TSS-seq ladder written by Fig3C_prevalence_ladder_panel.R itself
python3 - <<'PY'
import json, os
P = '/Users/wu/Desktop/Weak_uORF_submission/output_revised_figs_20261028/png_for_docx'
m = {'main': {str(n): f'{P}/Figure{n}.png' for n in range(1, 8)}, 'supp': {str(s): f'{P}/S{s}.png' for s in range(1, 23)}}   # V14: S1-S19 (the mORF figure removed), Figure 7 from Figure7_revised_model.pdf
for k in m.values():
    for f in k.values(): assert os.path.exists(f), f
json.dump(m, open('/Users/wu/Desktop/weak_uORF_revision/revision_docx/figure_images.json', 'w'), indent=1); print('mapping written:', len(m['main']), 'main,', len(m['supp']), 'supplemental')
PY
[ "${FIGS_ONLY:-0}" = 1 ] && { echo "FIGS_ONLY=1: figures, composites and PNG map done; run finalize_v21.sh for the docx"; exit 0; }
# 4) docx: build, embed, zip, validate, copies, package
cd /Users/wu/Desktop/weak_uORF_revision/revision_docx; N=${N:-Weak_uORF_2026_11_03_revision_v14}
rm -rf unpacked_work && cp -r unpacked unpacked_work && python3 apply_revision_edits.py unpacked_work . > apply.log 2>&1 && tail -1 apply.log
python3 embed_figure_images.py unpacked_work figure_images.json 200 8.5
python3 recolor_minimal.py unpacked_work unpacked
(cd unpacked_work && rm -f ../${N}_blue.docx && zip -Xr ../${N}_blue.docx . > /dev/null)
V="/Users/wu/Library/Application Support/Claude/local-agent-mode-sessions/skills-plugin/1bcf2be6-180f-44a5-8d3f-5c567f226d97/5b3ab5c7-ab57-4641-b1e3-530ff268e340/skills/docx/scripts/office/validate.py"
python3 "$V" ${N}_blue.docx --original base_for_validation.docx 2>&1 | tail -1
python3 make_submission_copies.py unpacked_work ${N}_marked.docx ${N}_clean.docx | cut -c1-80
for f in ${N}_marked.docx ${N}_clean.docx; do printf "%s: " $f; python3 "$V" $f --original base_for_validation.docx 2>&1 | tail -1; done
cp ${N}_blue.docx ../${N}_blue.docx
python3 build_submission_package.py unpacked_work /Users/wu/Desktop/weak_uORF_revision/submission_package_20261028 | tail -2
ls -la ../${N}_blue.docx ${N}_clean.docx | awk '{print $5, $9}'
