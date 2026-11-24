################################################################################
# A20 (2026-11-01, Larry) — TuORFs that coincide with the CDS of a neighbouring gene, and the shipped table files            #
#                                                                                                                              #
# (1) For every TuORF of Table S1A, the share of its spliced length that lies inside the CDS exons of another gene on the same  #
#     strand (expressed-isoform annotation). 14 TuORFs lie (>= 70%) inside a neighbouring gene's CDS: the Ribo-seq signal is  #
#     that gene's main ORF, read through an overlapping transcript model. Larry's decision (2026-11-01): keep the numbers,    #
#     flag the entries in the tables (new column 'neighbouring_gene_CDS') and say so in the Methods and table legends.        #
# (2) Writes the revised Supplemental Tables S1 and S3 to weak_uORF_revision/revision_tables/ with the flag column and the     #
#     corrected legend rows (ORF_type, TL_region_nt, uIF = uORF intensity factor), and patches the legends / README rows of   #
#     the shipped Tables S5 and S6 (group label, table pointers). Values are unchanged.                                      #
# Outputs: output/A20_Table_TuORFs_in_neighbouring_gene_CDS.csv, output/A20_summary.txt, revision_tables/TableS1_*, S3_*,    #
#          S5_*, S6_* (rewritten).                                                                                            #
# Run from Weak_uORF_submission/:  Rscript revision_modules/A20_flag_neighbouring_gene_CDS_and_table_legends.R               #
################################################################################
suppressPackageStartupMessages({ library(data.table); library(readxl); library(openxlsx); library(GenomicRanges); library(rtracklayer) })
.log <- character(); .add <- function(...) { m <- sprintf(...); cat(m, "\n"); .log <<- c(.log, m) }
OUT <- "output"; TABS <- "/Users/wu/Desktop/weak_uORF_revision/revision_tables"
PATH_GTF2 <- "Files/Araport11_20181206_max2_expressed_isoform_TPM0.5_rm_diff_CDSstart.gtf"
MIN_FRAC <- 0.70
S1 <- "output/Table_S1_all_TuORFs.xlsx"; S2 <- "output/Table_S2_single_TuORF_genes.xlsx"
s1a <- as.data.table(read_excel(S1, sheet = "Table S1A (CAGE)")); s1b <- as.data.table(read_excel(S1, sheet = "Table S1B (TSS-seq)")); leg1 <- as.data.frame(read_excel(S1, sheet = "Column legend"))
s3 <- as.data.table(read_excel(S2, sheet = 1)); leg3 <- as.data.frame(read_excel(S2, sheet = "Column legend"))

## (1) overlap with another gene's CDS, same strand; fraction of the spliced uORF length
gtf <- import(PATH_GTF2); cds <- gtf[gtf$type == "CDS"]; cds$gene <- cds$gene_id
gr <- GRanges(as.character(s1a$chr), IRanges(pmin(s1a$uORF_start_genomic, s1a$uORF_stop_genomic), pmax(s1a$uORF_start_genomic, s1a$uORF_stop_genomic)), strand = s1a$Strand)
ov <- findOverlaps(gr, cds); ov <- ov[s1a$gene_id[queryHits(ov)] != cds$gene[subjectHits(ov)]]
other <- character(nrow(s1a)); frac <- numeric(nrow(s1a))
for (i in unique(queryHits(ov))) {
  gs <- unique(cds$gene[subjectHits(ov)[queryHits(ov) == i]]); o <- reduce(cds[cds$gene %in% gs])
  frac[i] <- min(1, sum(width(GenomicRanges::intersect(gr[i], o))) / s1a$uORF_length_nt[i]); other[i] <- paste(gs, collapse = ";")
}
s1a[, `:=`(neighbouring_gene_CDS = fifelse(frac >= MIN_FRAC, other, ""), frac_in_neighbouring_CDS = round(frac, 2))]
flag <- s1a[neighbouring_gene_CDS != "", .(gene_id, tx_id, uORF_genome_range, neighbouring_gene_CDS, frac_in_neighbouring_CDS, uORF_length_aa, n_TuORFs_in_gene, uORF_P_sites, uIF, pct_uORF_plus, in_single_TuORF_cohort = gene_id %in% s3$gene_id)]
fwrite(flag, file.path(OUT, "A20_Table_TuORFs_in_neighbouring_gene_CDS.csv"))
.add("A20: %d of %d TuORFs lie (>= %.0f%% of their spliced length) inside the CDS of a neighbouring gene on the same strand; %d of them are in the single-TuORF cohort (Table S3); any overlap: %d",
     nrow(flag), nrow(s1a), 100 * MIN_FRAC, sum(flag$in_single_TuORF_cohort), sum(frac > 0))
for (i in seq_len(nrow(flag))) .add("  %s %s = CDS of %s (%.0f%%), %d aa, uIF %.3f, %%uORF+ %.1f%s", flag$gene_id[i], flag$uORF_genome_range[i], flag$neighbouring_gene_CDS[i], 100 * flag$frac_in_neighbouring_CDS[i], flag$uORF_length_aa[i], flag$uIF[i], flag$pct_uORF_plus[i], ifelse(flag$in_single_TuORF_cohort[i], " (cohort)", ""))

## (1b) V23 (2026-11-08): TuORFs whose listed isoform (Araport11+CTRL) does not contain the start and stop codon in its exons, and groups of TuORFs that share
##      gene, strand and stop codon (calls of the two callers that start at different AUGs; listed as separate TuORFs, Methods)
suppressPackageStartupMessages({ library(GenomicFeatures) })
txdbA <- if (requireNamespace("txdbmaker", quietly = TRUE)) txdbmaker::makeTxDbFromGFF("Files/Araport11+CTRL_20181206.gtf", format = "gtf") else GenomicFeatures::makeTxDbFromGFF("Files/Araport11+CTRL_20181206.gtf", format = "gtf")
exA <- exonsBy(txdbA, by = "tx", use.names = TRUE)
inex <- function(t, pos) { e <- exA[[t]]; !is.null(e) && any(start(e) <= pos & end(e) >= pos) }
s1a[, not_in_isoform_exons := !mapply(function(t, a, b) t %in% names(exA) && inex(t, a) && inex(t, b), tx_id, uORF_start_genomic, uORF_stop_genomic)]
s1a[, stop_key := paste(gene_id, Strand, uORF_stop_genomic)]; s1a[, n_same_stop := .N, by = stop_key]
grp_ids <- s1a[n_same_stop > 1, unique(stop_key)]; s1a[, same_stop_group := fifelse(n_same_stop > 1, paste0("S", match(stop_key, grp_ids)), "")]; s1a[, c("stop_key", "n_same_stop") := NULL]
.add("V23: %d TuORFs lie outside the exons of their listed isoform; %d TuORFs in %d groups share gene, strand and stop codon (different start codons)",
     sum(s1a$not_in_isoform_exons), sum(s1a$same_stop_group != ""), length(grp_ids))
## (1c) V23: inclusion-adjusted uIF for Table S3 (module A3; genes with at least 1% uORF+)
a3 <- fread(file.path(OUT, "A3_single_TuORF_cohort_uIF_variants.csv"))[, .(gene_id, uORF_genome_range, uIF_incl)]; a3[, uORF_genome_range := sub("^[^|]*\\|", "", uORF_genome_range)]   # A3 writes gene|range

## (2) tables: flag column + legend rows
key <- function(d) paste(d$gene_id, d$uORF_genome_range)
s1b[, neighbouring_gene_CDS := s1a$neighbouring_gene_CDS[match(key(s1b), key(s1a))]]
s3[, neighbouring_gene_CDS := s1a$neighbouring_gene_CDS[match(key(s3), key(s1a))]]
s3[, uIF_inclusion_adjusted := round(a3$uIF_incl[match(key(s3), key(a3))], 3)]; s3[, not_in_isoform_exons := s1a$not_in_isoform_exons[match(key(s3), key(s1a))]]
s1b[, `:=`(not_in_isoform_exons = s1a$not_in_isoform_exons[match(key(s1b), key(s1a))], same_stop_group = s1a$same_stop_group[match(key(s1b), key(s1a))])]
.add("  Table S3: inclusion-adjusted uIF for %d of %d genes (the others have < 1%% uORF+)", sum(is.finite(s3$uIF_inclusion_adjusted)), nrow(s3))
s1a[, frac_in_neighbouring_CDS := NULL]
LEG_NEW <- sprintf("Gene whose annotated CDS (same strand) contains at least %.0f%% of the TuORF's spliced length; the Ribo-seq signal of such an ORF is the main ORF of that neighbouring gene, read through an overlapping transcript model (%d TuORFs; the entries are kept and flagged, Methods)", 100 * MIN_FRAC, nrow(flag))
fix_legend <- function(leg) {
  stopifnot(ncol(leg) >= 2); cn <- names(leg); col <- leg[[1]]; des <- leg[[2]]
  des[col == "ORF_type"] <- "uORF = non-overlapping uORF (stop codon within the 5' leader); overlapping uORFs (start codon in the leader, stop codon inside the coding sequence (CDS) of the mORF) are not treated as TuORFs and are absent from this table (Methods)"
  des[col == "TL_region_nt"] <- "Length (nt) of the region used for uORF_TL: the coding length of the uORF (equal to uORF_length_nt)"
  des[col == "uIF"] <- "uORF intensity factor (relative ribosome occupancy): uORF_TL / mORF_TL (rounded to 3 decimals)"
  des[col == "pct_uORF_plus"] <- "%uORF+: percentage of the CAGE-seq (Table S1A) or TSS-seq (Table S1B) signal in the annotated 5' leader of the listed isoform plus 100 nt upstream that lies at or upstream of the first nucleotide of the start codon"
  des[col == "RNA_TPM"] <- sub("5'UTR-containing, RNA TPM > 1", "genes with an annotated 5' leader and RNA TPM > 1", des[col == "RNA_TPM"], fixed = TRUE)   # V41 (2026-11-12): the term of the manuscript
  leg[[2]] <- des; add <- function(l, c1, c2) rbind(l, setNames(as.data.frame(as.list(c(c1, c2, rep("", ncol(leg) - 2)))), cn))
  out <- add(leg, "neighbouring_gene_CDS", LEG_NEW)
  out <- add(out, "not_in_isoform_exons", "TRUE where the exons of the listed isoform do not contain the start and stop codon (15 TuORFs; their %uORF+ is computed on that isoform's window, which does not contain the uORF)")
  out }
leg_s1 <- function(leg) add_row2(fix_legend(leg), "same_stop_group", "TuORFs with the same label share gene, strand and stop codon but start at different AUGs (calls of the two callers or of different isoforms; listed as separate TuORFs, Methods)")
add_row2 <- function(l, c1, c2) rbind(l, setNames(as.data.frame(as.list(c(c1, c2, rep("", ncol(l) - 2)))), names(l)))
leg_s3 <- function(leg) {   # V40 (2026-11-11): four rows are written for Table S3 (they had the text of Table S1, and the uIF classes were said to be in Figure 1)
  out <- add_row2(fix_legend(leg), "uIF_inclusion_adjusted", "Inclusion-adjusted uIF: uIF divided by %uORF+/100 (ribosome occupancy per uORF-containing transcript), for genes with at least 1% uORF+ (Supplemental Fig. S4)")
  set <- function(col, txt) { stopifnot(sum(out[[1]] == col) == 1); out[[2]][out[[1]] == col] <<- txt }
  set("uIF_Group", "uIF class as in Supplemental Fig. S2B: High (uIF > 1), Medium (0.2 < uIF <= 1), Low (uIF <= 0.2)")
  s1l <- fix_legend(leg1); for (cc in c("RNA_TPM", "Ribo_TPM", "mORF_TE")) set(cc, s1l[[2]][s1l[[1]] == cc])   # V41 (2026-11-12): the three rows say when a value is NA, with the text of Table S1 (77 genes of Table S3 are not in the TE analysis set)
  set("pct_uORF_plus", "%uORF+: percentage of the CAGE-seq signal in the annotated 5' leader of the listed isoform plus 100 nt upstream that lies at or upstream of the first nucleotide of the start codon (the values of Table S1A)")
  set("neighbouring_gene_CDS", sprintf("Gene whose annotated CDS (same strand) contains at least %.0f%% of the TuORF's spliced length; the Ribo-seq signal of such an ORF is the main ORF of that neighbouring gene, read through an overlapping transcript model (%d TuORFs in Table S1A, %d of them in this table; the entries are kept and flagged, Methods)", 100 * MIN_FRAC, nrow(flag), sum(s3$neighbouring_gene_CDS != "", na.rm = TRUE)))
  set("not_in_isoform_exons", sprintf("TRUE where the exons of the listed isoform do not contain the start and stop codon (%d TuORFs in Table S1A, %d of them in this table; their %%uORF+ is computed on that isoform's window, which does not contain the uORF)", sum(s1a$not_in_isoform_exons), sum(s3$not_in_isoform_exons, na.rm = TRUE)))
  out }
wb <- createWorkbook(); addWorksheet(wb, "Table S1A (CAGE)"); writeData(wb, "Table S1A (CAGE)", s1a); addWorksheet(wb, "Table S1B (TSS-seq)"); writeData(wb, "Table S1B (TSS-seq)", s1b)
addWorksheet(wb, "Column legend"); writeData(wb, "Column legend", leg_s1(leg1)); saveWorkbook(wb, file.path(TABS, "TableS1_all_TuORFs_revised.xlsx"), overwrite = TRUE)
wb <- createWorkbook(); addWorksheet(wb, "Table S3"); writeData(wb, "Table S3", s3); addWorksheet(wb, "Column legend"); writeData(wb, "Column legend", leg_s3(leg3)); saveWorkbook(wb, file.path(TABS, "TableS3_single_TuORF_genes_revised.xlsx"), overwrite = TRUE)
.add("  Tables S1 and S3 rewritten with the column 'neighbouring_gene_CDS' (%d / %d / %d flagged rows in S1A / S1B / S3)", sum(s1a$neighbouring_gene_CDS != ""), sum(s1b$neighbouring_gene_CDS != "", na.rm = TRUE), sum(s3$neighbouring_gene_CDS != "", na.rm = TRUE))

# Table S5 (leader-splicing audit): uIF legend row
f8 <- file.path(TABS, "TableS5_TuORF_leader_splicing_audit.xlsx"); w8 <- loadWorkbook(f8); l8 <- read.xlsx(f8, sheet = "Column legend")
l8$Description[l8$Column == "uIF"] <- "uORF intensity factor (relative ribosome occupancy; Table S1A)"; writeData(w8, "Column legend", l8); saveWorkbook(w8, f8, overwrite = TRUE)
# Table S6 (undetected-uORF genes): the CAGE signal column is TPM-scaled signal, not a read count
f9 <- file.path(TABS, "TableS6_undetected_uORF_genes_revised.xlsx"); w9 <- loadWorkbook(f9); l9 <- read.xlsx(f9, sheet = "Column legend")   # V41 (2026-11-12): file TableS7_UuORF_genes_revised.xlsx renamed
# V41 (2026-11-12): the abbreviation UuORF is written out in the legend and in three column names of the data sheet (only the header row is written)
hd9 <- names(read.xlsx(f9, sheet = 1, rows = 1:2)); hd9n <- sub("^uuORF_", "undetected_uORF_", hd9); stopifnot(sum(grepl("^undetected_uORF_", hd9n)) == 3)
writeData(w9, sheet = 1, x = as.data.frame(t(hd9n), stringsAsFactors = FALSE), startRow = 1, colNames = FALSE)
l9[[1]] <- sub("^uuORF_", "undetected_uORF_", l9[[1]])
set9 <- function(col, txt) { stopifnot(sum(l9[[1]] == col) == 1); l9[[2]][l9[[1]] == col] <<- txt }
set9("gene_id", "Arabidopsis gene identifier (AGI). The genes have predicted but undetected uORFs and are not identified as TuORF genes by RiboTaper or CiPS")
set9("undetected_uORF_id", "Identifier of the undetected uORF (transcript_rank from ORF finding)")
set9("undetected_uORF_start_genomic", "Genomic coordinate of the start codon of the undetected uORF (biological orientation)")
set9("undetected_uORF_length_aa", "Peptide length of the undetected uORF in amino acids (stop codon excluded)")
set9("peptide_seq", "Predicted peptide sequence of the undetected uORF (trailing * = stop codon)")
l9[[2]][l9[[1]] == "RNA_TPM"] <- sub("5'UTR-containing, RNA TPM > 1", "genes with an annotated 5' leader and RNA TPM > 1", l9[[2]][l9[[1]] == "RNA_TPM"], fixed = TRUE)
l9[[2]][l9[[1]] == "TSS_total_count"] <- "Total TPM-scaled CAGE-seq signal within the annotated 5' leader + 100 nt upstream (0 = no detectable signal; such genes are scored as 0% uORF+ and fall in the <30% group; 58 genes)"
stopifnot(sum(l9[[1]] == "mORF_TE") == 1)
l9[[2]][l9[[1]] == "mORF_TE"] <- "mORF translation efficiency: Ribo_TPM / RNA_TPM; NA as above. Figure S10E uses a lower cutoff (RNA TPM > 0.1), so it includes genes that have NA here"   # V40 (2026-11-11): the key of S10E counts 2,309 / 117 / 136 genes with TE; this table has a TE value for 2,266 / 109 / 124
writeData(w9, "Column legend", l9); saveWorkbook(w9, f9, overwrite = TRUE)
# V23 (2026-11-08): the Gene Ontology table (Supplemental Table S8) is written by A24 only (the block that patched it here would restore the pre-V19 group legend and turn the text odds ratios into Inf).
.add("  Tables S5 and S6 legends patched")
writeLines(.log, file.path(OUT, "A20_summary.txt"))
