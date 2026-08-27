# TSS heterogeneity and the uORF landscape in Arabidopsis

Analysis code for the manuscript **"Transcription Start Site Heterogeneity Confounds the Landscape and Functional Interpretation of uORFs"**.

Upstream open reading frames (uORFs) are usually annotated on the longest transcript isoform, but transcription start site (TSS) selection determines whether a uORF is actually present in an mRNA. This repository contains a single R script, `20260827_weak_uORFs.R`, that computes, for every translated uORF (TuORF) identified by Ribo-seq (RiboTaper + CiPS), the fraction of transcripts that include the uORF (**%uORF+**, called `pUP` in the code) from TSS profiling data (CAGE-seq by default, or TSS-seq from Nielsen et al. 2019), and combines it with RNA-seq/Ribo-seq quantification, reporter assays, and public datasets to ask how TSS selection shapes uORF inclusion and function. The script reproduces the figures, supplementary tables, and headline numbers of the manuscript.

## What the script does, section by section

The script runs top to bottom as one pipeline; later modules reuse objects built in sections 1–4.

1. **Setup (section 0)** — input/output locations (`DATA_DIR = "Files"`, `OUT_DIR_MAIN = "output"`), output helpers, libraries, and existence checks for every input file.
2. **Section 1: Load TuORFs and annotation** — reads RiboTaper ORF calls and CiPS TuORF calls, merges them into the TuORF set, and builds the Araport11 TxDb.
3. **Section 2: mORF translation metrics** — Ribo-seq P-site counts per main ORF (mORF) CDS, giving mORF translation levels.
4. **Section 3: %uORF+ (pUP) for TuORFs** — for each TuORF, the fraction of TSS signal (CAGE-seq or TSS-seq, chosen by `TSS_SOURCE`) upstream of the uORF start within the annotated 5′ UTR + 100 nt window.
5. **Section 4: Translation efficiency and uIF** — RSEM RNA/Ribo TPM per gene, the single-TuORF gene set, and the uORF impact factor (uIF = uORF translation level / mORF translation level).
6. **Table S2 export** — the single-TuORF gene analysis set with uIF groups (Strong / Medium / Weak).
7. **Table S1 export** — all TuORFs with %uORF+ computed from both CAGE-seq and TSS-seq.
8. **Figure 1 C–E** — uIF distribution and uIF vs mORF translation efficiency (TE).
9. **Figure 3 B–H (script numbering)** — %uORF+ distributions, %uORF+ vs uIF and TE, and the mORF-start pUP control; these panels feed manuscript Figure 2 and Figure S4.
10. **Figure 1&3 add-on** — RNA half-life analyses joined to the uIF/%uORF+ sets.
11. **Figure 4: UuORF module** — undetected uORFs (UuORFs): AUG-containing 5′ UTRs without detected translation, their %uORF+, TE, and RNA levels (manuscript Figure S5).
12. **Table S3 export** — UuORF genes.
13. **Figure S5D** — %uORF+ distribution of UuORF genes, styled to match Figure 2B.
14. **Figure 5: CPuORF module** — conserved peptide uORFs vs other TuORFs vs no-TuORF genes (manuscript Figure 4).
15. **Figure 6: ATG frequency around TSSs** — ATG frequency upstream/downstream of the 50% cumulative CAGE position in TuORF-less vs TuORF-containing genes (manuscript Figure 3).
16. **Figure 7: reporter assays** — dual-luciferase (DLA) and qPCR measurements for weak-uORF/UuORF constructs, including perfect-Kozak variants (manuscript Figure 6).
17. **Figure 8: uORF length analyses** — length distributions by detection set and by TSS filtering, including panel H (uORF length before vs after removing uORFs with %uORF+ ≤ 30) (manuscript Figure 5).
18. **GO SLIM enrichment** — BP/MF/CC bubble plots for two gene groupings (by uIF and by %uORF+), each with a shared fourth reference group (manuscript Figures S7–S9).
19. **Figure 9: mRNA half-life groups** — half-life by uIF groups and by %uORF+ groups.
20. **Figure 11: uORF order** — log2 uIF by uORF position in multi-uORF genes, plus per-rank %uORF+ summaries (manuscript Figure S3).
21. **Transcriptome-wide TuORF impact** — multi-isoform-aware estimate of the transcriptome fraction containing TuORFs before and after accounting for TSS selection (the manuscript's 21.6% vs 9.6%, Figure 3C); runs for both CAGE-seq and TSS-seq and writes one table and figure per source.
22. **uORF length-bin plots** — RiboTaper vs CiPS length-bin comparison.
23. **Kozak module** — Kozak context scores, ECDFs, and sequence logos for uORFs (by uIF group) and mORFs (five groups).
24. **Correlation module** — distances (uORF stop → CDS start; 50% CAGE position → CDS start or uORF start) and Kozak score vs TE, before and after %uORF+/RNA filters.
25. **Reconciliation module** — recomputes every headline count in the manuscript (TuORF totals, UuORF groups, uORFs per gene, length bins) as a traceable consistency check.

Script figure numbers are internal; where they differ from the manuscript: script Figure 3 → manuscript Figure 2/S4, script Figure 4 → S5, script Figure 5 → 4, script Figure 6 → 3, script Figure 7 → 6, script Figure 8 → 5, script Figure 11 → S3.

## Requirements

R ≥ 4.3 with Bioconductor (tested with R 4.3.3 / Bioconductor 3.18).

- Bioconductor: `GenomicFeatures`, `GenomicRanges`, `GenomeInfoDb`, `IRanges`, `S4Vectors`, `BiocGenerics`, `Biostrings`, `Rsamtools`, `ORFik`, `GO.db`, `AnnotationDbi`
- CRAN: `dplyr`, `tidyr`, `data.table`, `readr`, `ggplot2`, `ggpubr`, `rstatix`, `scales`, `openxlsx`, `knitr`, `rcompanion`, `multcompView`, `ggseqlogo`, `forcats`, `stringr`, `RColorBrewer`

## Input files

All inputs go in `Files/` (see `DATA_DIR` in section 0 of the script).

| File | Description |
|---|---|
| `ORFs_max_filt_expressed` | RiboTaper ORF calls (expressed set) |
| `CiPS_uORF_Final.xlsx` | CiPS TuORF calls |
| `Araport11+CTRL_20181206.gtf` | Araport11 annotation used throughout |
| `Araport11_20181206_max2_expressed_isoform_TPM0.5_rm_diff_CDSstart.gtf` | Expressed-isoform annotation (Figure 4 TxDb) |
| `CTRL_expressed_P_sites_sort_count` | Ribo-seq P-site positions/counts |
| `wt_R123_genome_RiboPlotR.txt` | CAGE-seq TSS signal (default `TSS_SOURCE`) |
| `wt_TSS_Nielsen_2019.txt` | TSS-seq signal (Nielsen et al. 2019) |
| `RSEM_Araport11_1st+2nd_Riboseq_RNA.genes.results_CDS_only` | RNA-seq RSEM, gene level |
| `RSEM_Araport11_1st+2nd_Riboseq_Ribo.genes.results_CDS_only` | Ribo-seq RSEM, gene level |
| `genes_w_5UTR.tsv` | Genes with an annotated 5′ UTR |
| `CPuORF_collection_3_manully_fix.xlsx` | Curated CPuORF collection |
| `TAIR10_chr_all_2.fas` | TAIR10 genome sequence |
| `RNA-half-lives-5EU.tsv` / `RNA-half-lives-5EU.txt` | RNA half-life dataset (5-EU labeling) |
| `ATH_GO_GOSLIM.tsv` | TAIR GO / GOslim annotations |
| `RSEM_isoforms_TPM_6samples.tsv` | Per-isoform RSEM TPM (6 samples in one table) |
| `DLA_weak_uORF_UuORFs.xlsx` | Dual-luciferase assay measurements |
| `qPCR_weak_uORF_uuORF.xlsx`, `qPCR_weak_uORF_uuORF_perfect_Kozak.xlsx` | qPCR measurements |

Note: `CTRL_expressed_P_sites_sort_count` (~220 MB) and `TAIR10_chr_all_2.fas` (~116 MB) exceed GitHub's 100 MB per-file limit, so this repository ships them as `CTRL_expressed_P_sites_sort_count.zip` and `TAIR10_chr_all_2.fas.zip` instead — **unzip both into `Files/` before running the script**. The TAIR10 genome can also be downloaded directly from [TAIR](https://www.arabidopsis.org).

## How to run

```bash
git clone https://github.com/hsinyenwu/weak_uORFs_2026.git
cd weak_uORFs_2026               # the folder containing the script, with inputs in Files/
Rscript 20260827_weak_uORFs.R
```

The script runs as a single pass and writes everything to `output/` (created automatically). To recompute %uORF+ from TSS-seq instead of CAGE-seq, set `TSS_SOURCE <- "TSS"` in section 0 (Table S1 and the transcriptome-impact module always report both sources in one run). A full run takes roughly 25–35 minutes and up to ~10 GB of RAM, most of it in the TxDb/FASTA and CAGE overlap steps.

## Outputs

All results are written to `output/`: figure PDFs (e.g., `Figure1_CDE.pdf`, `Figure3_BCDEFGH.pdf`, `Figure8_ABCDEFGH_Full_Violin.pdf`, GO bubble plots), supplementary tables (`Table_S1_all_TuORFs.xlsx`, `Table_S2_single_TuORF_genes.xlsx`, `Table_S3_UuORF_genes.xlsx`), the transcriptome-impact tables (`Table_Transcriptome_Impact_MultiIsoform_CAGE.csv` / `_TSSseq.csv`), per-rank %uORF+ summaries, and the reconciliation files (`Reconciliation_A_numbers.txt`, `Reconcile_A*.csv`).

## Citation

If you use this code, please cite the manuscript above (in preparation; citation will be updated upon publication).
