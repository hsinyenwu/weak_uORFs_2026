# A25_tables.R — the shipped tables of A25: output/A25_Table_ouORF_leader_splicing.xlsx and the second sheet of Supplemental Table S5 with its legend.
# Sourced by A25_ouORF_leader_splicing_audit.R after it has written output/A25_Table_ouORF_leader_splicing.csv, or run alone (V41, 2026-11-12):
#   Rscript revision_modules/A25_tables.R        (from Weak_uORF_submission/; reads the saved csv, so the BAM scan of A25 is not repeated)
if (!exists("out")) {
  suppressPackageStartupMessages({ library(data.table); library(openxlsx) })
  OUT_DIR <- "output"; .log <- character(); .add <- function(m) { cat(m, "\n"); .log <<- c(.log, m) }
  out <- as.data.frame(fread(file.path(OUT_DIR, "A25_Table_ouORF_leader_splicing.csv"), na.strings = c("NA", ""), colClasses = c(chr = "character")))   # chr is text in the table, as in A25
}
legend <- data.frame(Column = c("gene_id, tx_id, ORF_id2", "n_isoform_calls", "chr, Strand", "start_codon_genomic", "overlapping_uORF_genome_range_min / _max", "AA",
  "pct_uORF_plus_CAGE / _TSSseq", "gene_has_TuORF", "start_shared_with_TuORF", "n_isoforms", "gene_has_leader_intron", "src_has_leader_intron",
  "other_TIS_in_intron (_ex)", "src_intron_in_body", "other_body_overlap", "annot_class", "annot_class_expressed_isoforms", "annot_TIS_affected",
  "R_TIS, E_TIS, AS_TIS, cov_TIS", "S_TIS", "TIS_spliced", "junction_status", "pUP_eff_CAGE / _TSSseq", "proxy_status"),
  Meaning = c("Gene, source isoform and identifier of the overlapping uORF (CiPS list of overlapping uORFs, Files/CiPS_ouORF_all_Final.xlsx; one row per unique start and stop codon).", "Number of isoforms on which the same overlapping uORF was called (2 = listed twice in the CiPS list).", "Chromosome and strand.",
  "Genomic position of the first nucleotide of the start codon of the overlapping uORF.", "Genomic span of the overlapping uORF (start codon to stop codon inside the main ORF).", "Length of the overlapping uORF in amino acids.",
  "%uORF+ of the overlapping uORF (percentage of the CAGE-seq or TSS-seq signal at or upstream of its start codon; same window and formula as for the TuORFs of Supplemental Table S1; the higher value where it was called on two isoforms).",
  "The gene also has a translated uORF in Supplemental Table S1.", "The same start codon is listed as a TuORF start codon on another isoform (the stop codon depends on the isoform).",
  "Number of coding isoforms of the gene (Araport11+CTRL).", "At least one isoform of the gene has an intron in its 5' leader.", "The source isoform (tx_id) of the overlapping uORF has an intron in its 5' leader.",
  "The start codon lies inside an intron of another isoform of the gene (_ex: of another isoform in the expressed two-isoform set).",
  "The body of the overlapping uORF spans an intron of its own isoform (start codon kept; includes coding-sequence introns).", "An intron of another isoform overlaps the body of the overlapping uORF without covering the start codon.",
  "Annotation class: D = start codon inside an intron of another isoform; C1/C2 = introns inside the body, start codon kept; B = leader intron(s) elsewhere in the gene; A = no leader intron in any isoform.",
  "The same, counting only alternative isoforms in the expressed two-isoform set.", "TRUE for class D (the only class in which annotated splicing removes the start codon).",
  "RNA-seq reads (unique, sense) whose aligned block covers the start codon (R), whose splice gap spans it (E; gaps that are annotated introns of another gene only are not counted), antisense reads discarded (AS); cov = R + E.",
  "Fraction of reads keeping the start codon, R / (R + E); NA when fewer than 10 reads.", "S_TIS < 0.90.",
  "retained / spliced over / unassessable (fewer than 10 reads).", "Effective inclusion = %uORF+ x S_TIS (NA when S_TIS is NA).",
  "Start codon affected by splicing (class D or S_TIS < 0.90) / annotation-clean but junctions unassessable / clean."))
# V41 (2026-11-12): the abbreviation ouORF is written out in the shipped sheets (sheet name, two column names, two class labels); the csv file and the figure code keep the short names
out_x <- as.data.frame(out); names(out_x) <- sub("^ouORF_", "overlapping_uORF_", names(out_x)); out_x$annot_class <- gsub("ouORF body", "body of the overlapping uORF", out_x$annot_class, fixed = TRUE)
SHEET_OU <- "A25 overlapping uORF splicing"
wb <- createWorkbook(); addWorksheet(wb, SHEET_OU); writeData(wb, 1, out_x); addWorksheet(wb, "Column legend"); writeData(wb, 2, legend)
setColWidths(wb, 2, cols = 1:2, widths = c(34, 120)); saveWorkbook(wb, file.path(OUT_DIR, "A25_Table_ouORF_leader_splicing.xlsx"), overwrite = TRUE)
################################################################################
# 6) Supplemental Table S5: add the ouORF sheet and its legend rows next to the TuORF audit (revision_tables copy; backup kept)
################################################################################
TABS <- "/Users/wu/Desktop/weak_uORF_revision/revision_tables"; f8 <- file.path(TABS, "TableS5_TuORF_leader_splicing_audit.xlsx")   # Supplemental Table S5
if (file.exists(f8)) {
  bk <- file.path(TABS, "_replaced_by_V23_20261108"); dir.create(bk, showWarnings = FALSE)
  if (!file.exists(file.path(bk, basename(f8)))) file.copy(f8, file.path(bk, basename(f8)))
  a2x <- file.path(OUT_DIR, "A2_Table_TuORF_leader_splicing.xlsx")   # V23 (2026-11-08): TuORF sheet and legend from the current A2 output (before: from the previous copy of S8, so an A2 rerun never reached the table)
  s8_tu <- read.xlsx(a2x, sheet = "A2 TuORF leader splicing"); l8 <- read.xlsx(a2x, sheet = "Column legend")
  l8 <- l8[!grepl("^Sheet 'A25", l8$Column) & !l8$Column %in% legend$Column | seq_len(nrow(l8)) <= which(l8$Column == "proxy_status")[1], ]   # keep the TuORF legend only (idempotent)
  l8 <- l8[seq_len(which(l8$Column == "proxy_status")[1]), ]
  add <- data.frame(Column = c("", sprintf("Sheet '%s'", SHEET_OU), legend$Column),
                    Description = c("", "The same audit for the start codons of the 958 unique translated overlapping uORFs (CiPS; 965 calls, 7 duplicated across isoforms). Only the start codon is audited; the stop codon lies inside the main ORF. Columns:", legend$Meaning))
  l8 <- rbind(setNames(l8[, 1:2], c("Column", "Description")), add)
  wb8 <- createWorkbook(); addWorksheet(wb8, "A2 TuORF leader splicing"); writeData(wb8, 1, s8_tu)
  addWorksheet(wb8, SHEET_OU); writeData(wb8, 2, out_x)
  addWorksheet(wb8, "Column legend"); writeData(wb8, 3, l8); setColWidths(wb8, 3, cols = 1:2, widths = c(34, 120))
  saveWorkbook(wb8, f8, overwrite = TRUE)
  .add(sprintf("  Supplemental Table S5 rewritten with the sheet 'A25 overlapping uORF splicing' (%d rows) and its legend (%d legend rows in all); backup in %s", nrow(out), nrow(l8), bk))
}
