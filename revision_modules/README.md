# Revision analyses (`revision_modules/`)

This folder holds the analyses and figure scripts added for the revision of the manuscript "Transcription Start Site Heterogeneity Confounds the Landscape and Functional Interpretation of uORFs". The main pipeline, `20260827_weak_uORFs.R` in the repository root, is the code of the submitted manuscript, updated early in the revision (2026-10-27 to 2026-11-01: the latest CiPS calls, the %uORF+ window of the 5′ leader, the names of the uIF classes); it has not changed since. The modules here read its output tables (`output/Table_S1_all_TuORFs.xlsx`, `output/Table_S2_single_TuORF_genes.xlsx` and others), add the new analyses, and draw the revised figures.

Figure and table numbers below are those of the revised manuscript (version V56, 2026-11-24). The submitted manuscript and the main script use other numbers.

## How to run

Run every command from the repository root (the folder that contains `20260827_weak_uORFs.R`), with the inputs in `Files/`.

1. If the CAGE-seq file was replaced, run `Rscript revision_modules/correct_CAGE_positions.R` first. It moves every CAGE-seq position 1 nucleotide upstream (Supplemental Methods, "Positions of the CAGE-seq signal"). It reads the file as received, so it cannot shift the positions twice.
2. `python3 revision_modules/rerun_all_modules.py --fresh` runs the main script, the figure driver and every analysis module below. The steps run in the order that their inputs need, up to six at the same time. The logs are in `output/rerun_logs/`. Use this command after an input file or a script has changed. The script records how each step ended (`output/rerun_logs/<step>.json`). `python3 revision_modules/rerun_all_modules.py`, without `--fresh`, runs only the steps that are not recorded as ended without error. It does not check whether an input file or a script has changed. Use it to continue a run that stopped. `python3 revision_modules/rerun_all_modules.py A31 A42` runs only the named steps. Each step that a named step needs must have ended without error or be named in the same command. If not, the script stops with a message. The steps that read the outputs of the named steps are not run again. Do not give step names together with `--fresh`. The schematic of Supplemental Fig. S6 (`Fig_schematic_leader_splicing_audit.R`) and the relabelling of the submitted drawings (`relabel_submitted_figure_labels.py`) are run separately.
3. The final figure files are made by `revision_modules/finalize_figures.sh`. This script also builds the manuscript files on the authors' computer, with paths outside this repository.

What each step needs before it can run (from `rerun_all_modules.py`). The names in the column "Step" are the names that the command accepts:

| Step | Needs |
|---|---|
| A5, A6, A7, A8, A9, A15, A16, A17, A21, A22, A23, A30, A33, A35, A36, FigS2 (`FigS2_TSSseq_panels.R`) | main |
| main (the main script, `20260827_weak_uORFs.R`) | nothing |
| driver_build (`build_revised_figure_driver.py` builds the figure driver from the main script) | nothing |
| driver (the figure driver) | driver_build and main |
| A2 | main |
| A3, A14, A19, A25 | A2 |
| A20 | A3 |
| A11, A29 | A20 |
| A12 | A11 |
| A13, A39, A40, A41 | A12 |
| A31 | A12 and A2 |
| A38 | A12 and A31 |
| A42 | A31 |
| A25_tables | A25 and A20 |
| A25_figure, A26 | A25 |
| Fig3C_CAGE (`Fig3C_prevalence_ladder_panel.R` with CAGE-seq) | A26 |
| Fig3C_TSSseq (the same script with TSS-seq) | Fig3C_CAGE |
| A34 | A33 and A2 |
| A18 | A34 |
| A24 | A23 |
| A37 | A35 and A6 |

Requirements: the revision modules were last run with R 4.5.1. They use the packages listed in the main README, plus the Bioconductor packages `rtracklayer`, `GenomicAlignments` and `SummarizedExperiment` (and `txdbmaker` when it is installed; otherwise `GenomicFeatures`), and the CRAN packages `patchwork`, `gridExtra`, `gtable`, `cowplot`, `readxl` and `png`. Python 3 runs `rerun_all_modules.py` and the figure driver builder; `pypdf` and `PyMuPDF` are needed for `compose_main_figures_vector.py` and `relabel_submitted_figure_labels.py`; `pdftoppm` (Poppler) renders the figures.

Some inputs are not in `Files/`, and the scripts read them through absolute paths set at their top: the Ribo-seq and RNA-seq BAM files (A2, A5, A12, A14, A25, A34), the per-replicate RSEM tables (A3), the Ribo-seQC results (A5), the PARE files (A12), the half-lives of Sorenson et al. 2018 (A41, A42), the expression atlas (A35), the submitted figure files (`compose_main_figures.R`, `compose_main_figures_vector.py`, `fig_S10_S12_compact.R`, `relabel_submitted_figure_labels.py`) and the supplemental table files of the manuscript folder (A11, A12, A20, A24, A25_tables, A29, A42). Change these paths before running the modules elsewhere. A few helper scripts (for example `finalize_figures.sh`) also contain paths of the authors' computer for building the manuscript files.

## Where each figure and table comes from

### Main figures

| Figure | Made by |
|---|---|
| Figure 1 | A (schematic) and C–F (genome browser panels): the submitted figure with the new uIF-class wording (`relabel_submitted_figure_labels.py`) and the row letters C–F (file `composite_work/Fig1_submitted_relabelled_v16.pdf`; the script that changed the row letters is no longer on disk); B: figure driver. Assembled by `compose_main_figures.R` (raster) and `compose_main_figures_vector.py` (vector). |
| Figure 2 | A: submitted schematic; B, C, E, F and G: figure driver, laid out by `fig_main_layout.R`; D (table of TuORFs affected by TSS selection and by leader splicing): `A19_TuORF_TSS_vs_splicing_counts.R`. Assembled as Figure 1. |
| Figure 3 | A–B: figure driver; C (prevalence ladder): `Fig3C_prevalence_ladder_panel.R`, with the values of `A26_uORF_pool_fraction_TSS_splicing.R`. |
| Figure 4 | A–C: submitted genome browser panels; D–F: figure driver. |
| Figure 5 | figure driver |
| Figure 6 | figure driver (dual-luciferase assays and RT-qPCR) |
| Figure 7 | submitted model drawing, relabelled by `relabel_submitted_figure_labels.py` |

The figure driver is a copy of the main script with the revised panel letters and page layout. `build_revised_figure_driver.py` builds it, and it uses `fig_main_layout.R` and `fig_S10_S12_compact.R`.

### Supplemental figures

| Figure | Content | Made by |
|---|---|---|
| S1 | Ribo-seq quality metrics | `A5_riboseq_QC_figure.R` |
| S2 | uIF and mORF TE | figure driver (`fig_main_layout.R`) |
| S3 | Null models and uncertainty | `A6_nulls_and_uncertainty.R` |
| S4 | Inclusion-adjusted uIF | `A3_uIF_abundance_and_independence.R` |
| S5 | Translated overlapping uORFs | `A25_ouORF_leader_splicing_audit.R` and `A25_figure.R` |
| S6 | Schematic of the leader-splicing audit | `Fig_schematic_leader_splicing_audit.R` |
| S7 | Leader-splicing audit | `A2_leader_splicing_audit.R` |
| S8 | Figure 2 panels with TSS-seq | `FigS2_TSSseq_panels.R` |
| S9 | Genes with two, three or four TuORFs | figure driver (`fig_S10_S12_compact.R`) |
| S10 | Undetected uORFs | figure driver (`fig_S10_S12_compact.R`) |
| S11 | Sequence logos at the TSS | figure driver |
| S12 | ATG frequency with composition controls | `A7_ATG_composition_controls.R` |
| S13 | Prevalence ladder with TSS-seq | `Fig3C_prevalence_ladder_panel.R` (`LADDER_SOURCE=TSSseq`) |
| S14 | Composition of the mRNA of genes with a translated uORF | `A18_mRNA_weighted_TSS_splicing.R` |
| S15–S17 | Gene Ontology enrichment (MF, BP, CC) | `A24_GO_by_max_pUP_figures_table.R` (input from `A23_GO_pooled_mid_high.R`) |
| S18 | CPuORFs and matched controls | `A8_CPuORF_matched_controls.R` |
| S19 | Uncapped mRNA (GMUCT) and half-life | `A11_uncapped_vs_uIF.R` |
| S20 | Uncapped mRNA (PARE) | `A12_PARE_uncapped_vs_uIF.R` |
| S21 | Decay by uORF inclusion and leader splicing | `A31_decay_by_inclusion_and_splicing.R` |
| S22 | Second half-life dataset; uncapped mRNA versus half-life and mORF TE | `A42_FigureS22_Sorenson_uncapped.R` |

Most figure files in `output/` keep the name of their module (for example `A6_Figure_nulls_and_uncertainty.pdf`), and `finalize_figures.sh` copies them to their final names in `output_revised_figs_20261028/` (for example `FigureS3_revised_nulls_and_uncertainty.pdf`). Supplemental Figs. S13 and S22 are written under their final names by their scripts, and the pages of the figure driver are taken from its output folder.

### Supplemental tables

| Table | Content | Made by |
|---|---|---|
| S1 | All TuORFs | main script (`Table_S1_all_TuORFs.xlsx`), with the flag columns and legend of `A20` |
| S2 | Ribo-seq quality metrics | `A5` (`A5_Table_riboseq_QC.xlsx`) |
| S3 | Single-TuORF genes | main script (`Table_S2_single_TuORF_genes.xlsx`), with columns of `A3` and the flags and legend of `A20` |
| S4 | Null models and uncertainty | `A6` (`A6_Table_nulls_and_uncertainty.xlsx`) |
| S5 | Leader-splicing audit | `A2` (TuORFs) and `A25_tables.R` (overlapping uORFs) |
| S6 | Undetected-uORF genes | main script (`Table_S3_UuORF_genes.xlsx`) |
| S7 | ATG composition controls | `A7` (`A7_Table_ATG_composition_controls.xlsx`) |
| S8 | Gene Ontology enrichment | `A24` |
| S9 | CPuORF matched controls | `A8` (`A8_Table_CPuORF_matched_controls.xlsx`) |
| S10 | Uncapped mRNA (GMUCT) and half-life | `A11` (`A11_Table_uncapped_halflife.xlsx`) |
| S11 | Uncapped mRNA (PARE) | `A12` (`A12_Table_PARE_uncapped.xlsx`) |
| S12 | Decay by uORF inclusion and leader splicing; second half-life dataset | `A31` (`A31_Table_decay_inclusion_splicing.xlsx`) and `A42` |

## The modules

### Shared code (sourced by other scripts)

| File | What it does |
|---|---|
| `fig_compact.R` | Figure style: pages 6.5 in wide and at most 9 in tall, Arial, no text below 8 pt. |
| `fig_main_layout.R` | Layout of the data panels of Figures 2–5 and Supplemental Fig. S2 (legends inside the plots, aligned panels). |
| `splicing_rule.R` | The splicing factor s_k of a uORF (fraction of RNA-seq reads that keep its start codon), used by A18, A26, A33, A34 and A37. |
| `psite_extract.R`, `psite_extract.py` | P-site counts per position from a Ribo-seq BAM file. Rule A (reads of 24–28 nt, offsets 8–12 nt) is the rule of this study. The two versions give the same output. |
| `correct_CAGE_positions.R` | Moves every CAGE-seq position 1 nucleotide upstream (see "How to run"). |
| `A3_roi_TuORF_and_CDS_blocks.bed` | Data file: TuORF and CDS regions, made for the per-replicate P-site counts of A3 with `psite_extract.R --roi`. No script reads it in the current run. |

### Leader splicing

| Module | Question | Used in |
|---|---|---|
| `A2_leader_splicing_audit.R` | How often does leader splicing, rather than TSS selection, decide whether a TuORF is in the mRNA? Annotation classes and RNA-seq junction reads at the start and stop codons. | Supplemental Fig. S7, Table S5 |
| `A14_junction_supported_splicing.R` | How does the share of TuORFs affected by splicing change if an annotated splicing event counts only when junction reads support it? | check |
| `A19_TuORF_TSS_vs_splicing_counts.R` | How many TuORFs are excluded by TSS selection and how many by leader splicing? | Figure 2D |
| `A25_ouORF_leader_splicing_audit.R`, `A25_tables.R`, `A25_figure.R` | The audit of A2 for the translated overlapping uORFs (start codon only). | Supplemental Fig. S5, Table S5 |

### TSS selection, uIF and the null models

| Module | Question | Used in |
|---|---|---|
| `A3_uIF_abundance_and_independence.R` | How much of uIF comes from uORF inclusion? uIF and TE from independent Ribo-seq libraries; inclusion-adjusted uIF. | Supplemental Fig. S4, Table S3 |
| `A6_nulls_and_uncertainty.R` | Null models for multi-TuORF genes and for uORF length; uncertainty and sequencing-depth dependence of %uORF+. | Supplemental Fig. S3, Table S4 |
| `A13_union_leader_window_sensitivity.R` | Does %uORF+ change with a wider window (the union of the leaders of all isoforms)? | check |
| `A15_intronic_TSS_in_leaders.R` | How much TSS signal lies inside leader introns? | check |
| `A20_flag_neighbouring_gene_CDS_and_table_legends.R` | Which TuORFs lie inside the CDS of a neighboring gene? Writes the flag columns and legends of Supplemental Tables S1, S3, S5 and S6. | Tables S1, S3, S5, S6 |
| `A21_mORF_window_check.R` | Are %uORF+ and the mORF percentage calculated over comparable windows? How much CAGE-seq signal lies outside the %uORF+ window? | Supplemental Methods (signal inside the coding sequence or 3′ UTR) |
| `A22_uIF_TE_abundance_null.R` | Is part of the uIF–TE correlation expected from mRNA abundance alone? | check |
| `A29_pUP_TE_abundance_adjusted.R` | Is the %uORF+–TE relation explained by RNA abundance? Partial correlation and thirds of RNA abundance. | Supplemental Methods ("%uORF+ and mORF TE within strata of RNA abundance") |
| `A30_same_stop_TuORF_merge_sensitivity.R` | Do the results change when TuORFs that share a start or a stop codon are merged? | Supplemental Methods ("TuORFs that share a start or a stop codon") |
| `FigS2_TSSseq_panels.R` | The Figure 2 panels with TSS-seq instead of CAGE-seq. | Supplemental Fig. S8 |

### Prevalence in the transcriptome

| Module | Question | Used in |
|---|---|---|
| `A17_Figure3C_same_gene_mapping.R` | Figure 3C with each uORF mapped only to the transcripts of its own gene (the rule now used by A26). | check |
| `A26_uORF_pool_fraction_TSS_splicing.R` | The fraction of the nuclear mRNA that contains a uORF, before and after TSS weighting and leader splicing. | Figure 3C, Supplemental Fig. S13 |
| `A33_predicted_uORFs_with_overlapping.R` | The same fractions for sequence-predicted uORFs, with and without predicted overlapping uORFs. | Results (prevalence) |
| `A34_predicted_uORFs_TSS_splicing.R` | Leader splicing for the sequence-predicted uORFs. | input of A18 (Supplemental Fig. S14D) |
| `A18_mRNA_weighted_TSS_splicing.R` | Of the mRNA of genes with a translated uORF, how much contains it, lacks it because of TSS selection, or lacks it because of splicing? | Supplemental Fig. S14 |
| `A16_predicted_only_uORF_TE.R` | Are the predicted but undetected uORFs that remain after TSS weighting associated with lower TE? | check |
| `Fig3C_prevalence_ladder_panel.R` | Draws the prevalence ladder. | Figure 3C, Supplemental Fig. S13 |

### Sequence composition, conserved peptide uORFs and gene function

| Module | Question | Used in |
|---|---|---|
| `A7_ATG_composition_controls.R` | Is the drop of ATG frequency at the TSS larger than that of trinucleotides with the same base composition? | Supplemental Fig. S12, Table S7 |
| `A8_CPuORF_matched_controls.R` | Do CPuORFs have higher inclusion than non-conserved TuORFs from matched gene classes? | Supplemental Fig. S18, Table S9 |
| `A9_GO_figures_FDR.R` | Gene Ontology bubble plots with the FDR marked (grouping of the submitted manuscript). | earlier version of Supplemental Figs. S15–S17 |
| `A23_GO_pooled_mid_high.R` | Gene Ontology enrichment with all TuORF genes grouped by their highest %uORF+. | input of A24 |
| `A24_GO_by_max_pUP_figures_table.R` | Gene Ontology figures and table for that grouping. | Supplemental Figs. S15–S17, Table S8 |
| `A35_housekeeping_tau.R`, `A36_GO_umbrella_terms_low_group.R` | Are the genes of the low %uORF+ group housekeeping genes, and which broad terms describe them? | not part of the manuscript |
| `A37_low_group_expression_level.R` | Are the genes of the low %uORF+ group highly expressed? | Results (gene function), Supplemental Methods ("Expression level of the %uORF+ groups") |

### Ribo-seq quality

| Module | Question | Used in |
|---|---|---|
| `A5_riboseq_QC_figure.R` | Read lengths, reading frame per read length, metagenes at start and stop codons. | Supplemental Fig. S1, Table S2 |

### Uncapped mRNA and mRNA half-life

| Module | Question | Used in |
|---|---|---|
| `A11_uncapped_vs_uIF.R` | Uncapped mRNA (GMUCT, Yu et al. 2021) and half-life (Szabo et al. 2020) versus uIF, %uORF+ and mORF TE. | Supplemental Fig. S19, Table S10 |
| `A12_PARE_uncapped_vs_uIF.R` | The same with PARE (Nagarajan et al. 2019). | Supplemental Fig. S20, Table S11 |
| `A31_decay_by_inclusion_and_splicing.R` | Half-life and uncapped mRNA by uORF inclusion and by leader splicing at the uORF start codon. | Supplemental Fig. S21, Table S12 |
| `A38_matched_decay_comparison.R` | The comparison of A31 with matched genes (closest in exon length and RNA abundance). | Supplemental Methods |
| `A39_uncapped_TE_halflife_partial.R` | Does the uncapped ratio follow mORF TE or half-life? Partial correlations with the 5-EU half-lives. | check (the values in the Supplemental Methods come from A42) |
| `A40_uncapped_halflife_TE_figure.R` | First version of Supplemental Fig. S22C–H. | replaced by A42 |
| `A41_sorenson_halflife_check.R` | The decay results with a second half-life dataset (Sorenson et al. 2018). | check (A42 makes the figure) |
| `A42_FigureS22_Sorenson_uncapped.R` | Supplemental Fig. S22 and its sheets in Table S12. | Supplemental Fig. S22, Table S12 |

### Figure assembly and checks (no analysis)

| File | What it does |
|---|---|
| `build_revised_figure_driver.py` | Builds the figure driver from the main script. |
| `fig_S10_S12_compact.R` | Pages of Supplemental Figs. S9 and S10, sourced by the figure driver (the file name keeps old numbers). |
| `Fig_schematic_leader_splicing_audit.R` | Schematic of Supplemental Fig. S6. |
| `compose_main_figures.R`, `compose_main_figures_vector.py` | Composites of Figures 1–4 (raster and vector). |
| `relabel_submitted_figure_labels.py` | Replaces the old uORF class wording inside the submitted Figure 1 and Figure 7 drawings. |
| `finalize_figures.sh` | Copies the figure files to their final names and renders PNG images. |
| `rerun_all_modules.py` | Runs the steps in order (see "How to run"; `--fresh` alone runs every step again). |
| `rerun_figure_modules.sh`, `rerun_for_renumbering_v48.sh` | Older helper scripts that redraw groups of figures. |
| `render_figures_qa.sh`, `check_figs.py` | Render the figures for inspection; check page size, fonts and text size. |
| `compare_output_runs.R` | Compares two output folders of the main script, file by file. |

Folders `_superseded/` and `_pre_github_backup/` hold earlier copies of these scripts and are not part of the analysis.
