################################################################################
# A14 — Junction-supported leader splicing (Larry, 2026-10-29)                   #
# Question: the manuscript says TSS usage alone determines inclusion for 80% of #
#   TuORFs, counting a TuORF as affected by splicing when ANY annotated isoform  #
#   (Araport11+CTRL) carries an intron that removes its start codon or changes   #
#   its body/stop (A2 classes C1, C2, D), or when RNA-seq junction reads show    #
#   splicing over its start or stop codon (A2: S < 0.90). How does that change   #
#   if an annotated event counts only when STAR junction reads support it?      #
#   (Discussion of 2026-10-29: unsupported isoforms are set aside; no isoform    #
#   re-quantification.)                                                         #
# Junction reads: STAR alignments of the six merged control RNA-seq libraries   #
#   (the BAM of A2), unique mappers (MAPQ 255), sense strand (dUTP: read 2 is    #
#   sense), counted as fragments (read names) per exact intron coordinates.     #
# Support: K = 1, 3 (main) or 10 fragments, at two levels                       #
#   intron  — the intron that causes the event has >= K fragments               #
#   isoform — every intron of the isoform carrying that intron has >= K          #
# Classes and the junction view (S at the start / stop codon) are those of A2;  #
#   with all introns counted the classes must reproduce A2 exactly.             #
# Outputs (output/): A14_junction_supported_splicing.xlsx, A14_summary.txt.     #
# NOT part of the manuscript (results first).                                   #
# Run from Weak_uORF_submission/:  Rscript revision_modules/A14_junction_supported_splicing.R
################################################################################
suppressPackageStartupMessages({
  library(GenomicFeatures); library(GenomicRanges); library(GenomicAlignments); library(Rsamtools)
  library(data.table); library(openxlsx)
})
options(stringsAsFactors = FALSE, warn = 1)
t0 <- Sys.time(); .log <- character(0); .add <- function(...) { s <- paste0(...); message(s); .log <<- c(.log, s) }
DATA_DIR <- "Files"; OUT_DIR <- "output"
PATH_GTF <- file.path(DATA_DIR, "Araport11+CTRL_20181206.gtf")
PATH_BAM <- "/Users/wu/Desktop/CTRL_v1/RNA_CTRL_merged.bam"          # STAR 2.6.0c, D1/D2/D3/D5/D6/D7 merged, paired, dUTP (read 2 = sense)
PATH_S1  <- file.path(OUT_DIR, "Table_S1_all_TuORFs.xlsx")
PATH_S2  <- file.path(OUT_DIR, "Table_S2_single_TuORF_genes.xlsx")
PATH_A2  <- file.path(OUT_DIR, "A2_Table_TuORF_leader_splicing.csv")
K_SET <- c(1L, 3L, 10L); K_MAIN <- 3L
SPLICED <- "spliced over (S < 0.90 at TIS or stop)"; UNASS <- "unassessable (<10 reads at both sites)"
for (f in c(PATH_GTF, PATH_BAM, PATH_S1, PATH_S2, PATH_A2)) stopifnot(file.exists(f))
pc <- function(x) sprintf("%.1f%%", 100 * mean(x))
big <- function(x) trimws(format(x, big.mark = ","))

################################################################################
# 1) Annotation, as A2: introns of coding isoforms
################################################################################
.add("A14: TxDb from ", PATH_GTF)
txdb <- if (requireNamespace("txdbmaker", quietly = TRUE)) txdbmaker::makeTxDbFromGFF(PATH_GTF, format = "gtf") else GenomicFeatures::makeTxDbFromGFF(PATH_GTF, format = "gtf")
tx_all  <- transcripts(txdb, columns = c("tx_name", "gene_id"))
tx2gene <- setNames(as.character(unlist(mcols(tx_all)$gene_id)), mcols(tx_all)$tx_name)
cds_by  <- cdsBy(txdb, by = "tx", use.names = TRUE)
intr_by <- intronsByTranscript(txdb, use.names = TRUE); intr_by <- intr_by[names(intr_by) %in% names(cds_by)]
intr <- unlist(intr_by, use.names = TRUE); intr$tx_id <- names(intr); names(intr) <- NULL; intr$gene_id <- unname(tx2gene[intr$tx_id])

################################################################################
# 2) TuORFs (Table S1A) with A2's per-TuORF results
################################################################################
S1A <- as.data.table(read.xlsx(PATH_S1, sheet = "Table S1A (CAGE)"))
U <- S1A[, .(gene_id, tx_id, chr = as.character(chr), Strand, uORF_genome_range,
             g_min = as.integer(sub("^[^_]+_([0-9]+)_([0-9]+)$", "\\1", uORF_genome_range)),
             g_max = as.integer(sub("^[^_]+_([0-9]+)_([0-9]+)$", "\\2", uORF_genome_range)),
             TIS = as.integer(uORF_start_genomic))]
U[, uid := paste(gene_id, uORF_genome_range, sep = "|")]
A2 <- fread(PATH_A2)
A2 <- A2[, .(uid = paste(gene_id, uORF_genome_range, sep = "|"), pct_uORF_plus, uIF, mORF_TE, gene_has_leader_intron, annot_class_A2 = annot_class, junction_status, S_TIS)]
U <- merge(U, A2, by = "uid", all.x = TRUE, sort = FALSE)
stopifnot(!anyNA(U$junction_status), !any(duplicated(U$uid)))
.add(sprintf("  TuORFs: %s in %s genes; A2 results merged for all", big(nrow(U)), big(uniqueN(U$gene_id))))

################################################################################
# 3) TuORF-intron pairs in the same gene (A2 rules)
################################################################################
U_gr <- GRanges(U$chr, IRanges(U$g_min, U$g_max), strand = U$Strand)
ov <- findOverlaps(U_gr, intr, ignore.strand = FALSE)
pair <- data.table(i = queryHits(ov), j = subjectHits(ov))
pair[, `:=`(same_gene = U$gene_id[i] == intr$gene_id[j], same_tx = U$tx_id[i] == intr$tx_id[j],
            tis_inside = U$TIS[i] >= start(intr)[j] & U$TIS[i] <= end(intr)[j],
            within_body = start(intr)[j] > U$g_min[i] & end(intr)[j] < U$g_max[i])]
pair <- pair[same_gene == TRUE]
pair[, ev := fifelse(same_tx & within_body, "C1", fifelse(!same_tx & tis_inside, "D", fifelse(!same_tx & !tis_inside, "C2", NA_character_)))]
pair <- pair[!is.na(ev)]

################################################################################
# 4) Junction fragments per intron in the TuORF genes (STAR alignments)
################################################################################
gn <- genes(txdb); reg <- reduce(gn[names(gn) %in% unique(U$gene_id)], ignore.strand = TRUE)
bam_sl <- names(scanBamHeader(PATH_BAM)[[1]]$targets); reg <- reg[as.character(seqnames(reg)) %in% bam_sl]
.add(sprintf("A14: junction reads in %s TuORF-gene regions (%.1f Mb) from %s", big(length(reg)), sum(width(reg)) / 1e6, basename(PATH_BAM)))
J <- list()
for (ch in unique(as.character(seqnames(reg)))) {
  r <- reg[seqnames(reg) == ch]; seqlevels(r) <- bam_sl
  param <- ScanBamParam(which = r, what = c("flag", "qname"), mapqFilter = 255,
                        flag = scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE, isNotPassingQualityControls = FALSE))
  gal <- readGAlignments(PATH_BAM, param = param); n_all <- length(gal)
  gal <- gal[njunc(gal) > 0L]
  fl <- mcols(gal)$flag; s <- as.character(strand(gal))
  sense <- ifelse(bitwAnd(fl, 64L) > 0L, ifelse(s == "+", "-", "+"), s)          # dUTP: read 1 antisense, read 2 sense
  jn <- junctions(gal); nj <- elementNROWS(jn); ul <- unlist(jn, use.names = FALSE)
  d <- unique(data.table(chr = ch, start = start(ul), end = end(ul), strand = rep(sense, nj), qname = rep(mcols(gal)$qname, nj)))   # one count per fragment and junction
  J[[ch]] <- d[, .(frag = .N), by = .(chr, start, end, strand)]
  .add(sprintf("  chr %s: %s alignments, %s spliced, %s distinct junctions", ch, big(n_all), big(length(gal)), big(nrow(J[[ch]]))))
  rm(gal, jn, ul, d); invisible(gc())
}
J <- rbindlist(J)[, .(frag = sum(frag)), by = .(chr, start, end, strand)]
idt <- data.table(idx = seq_along(intr), chr = as.character(seqnames(intr)), start = start(intr), end = end(intr), strand = as.character(strand(intr)), tx_id = intr$tx_id, gene_id = intr$gene_id)
idt <- merge(idt, J, by = c("chr", "start", "end", "strand"), all.x = TRUE); idt[is.na(frag), frag := 0L]; setorder(idt, idx)
in_genes <- idt$gene_id %in% unique(U$gene_id)
.add(sprintf("  annotated introns of coding isoforms in TuORF genes: %s; with >= 1 / >= 3 / >= 10 fragments: %s / %s / %s",
             big(sum(in_genes)), pc(idt$frag[in_genes] >= 1), pc(idt$frag[in_genes] >= 3), pc(idt$frag[in_genes] >= 10)))
iso_min <- idt[, .(min_frag = min(frag)), by = tx_id]
pair[, `:=`(frag = idt$frag[j], iso_min = iso_min$min_frag[match(intr$tx_id[j], iso_min$tx_id)])]
for (e in c("D", "C2", "C1")) { z <- unique(pair[ev == e, .(j, frag)])
  .add(sprintf("  introns behind class %s events: %s; with 0 fragments %s, >= 3 %s, median %s fragments", e, big(nrow(z)), pc(z$frag == 0), pc(z$frag >= 3), big(median(z$frag)))) }

################################################################################
# 5) Classes and the TSS-proxy scope under each support rule
################################################################################
classify <- function(keep) {
  f <- pair[keep, .(C1 = any(ev == "C1"), D = any(ev == "D"), C2 = any(ev == "C2")), by = i]
  D <- C1 <- C2 <- rep(FALSE, nrow(U)); D[f$i] <- f$D; C1[f$i] <- f$C1; C2[f$i] <- f$C2
  fcase(D, "D", C1, "C1", C2, "C2", U$gene_has_leader_intron, "B", default = "A") }
status <- function(cls) fcase(cls %in% c("C1", "C2", "D") | U$junction_status == SPLICED, "affected by splicing",
                              U$junction_status == UNASS, "annotation-clean, junctions unassessable", default = "clean (annotation + junctions)")
gene_rollup <- function(st) { d <- data.table(gene_id = U$gene_id, st = st)
  d[, .(st = if (any(st == "affected by splicing")) "affected by splicing" else if (all(st == "clean (annotation + junctions)")) "clean (annotation + junctions)" else "annotation-clean, junctions unassessable"), by = gene_id] }
rules <- c(list(list(name = "all annotated introns (A2, manuscript)", keep = rep(TRUE, nrow(pair)))),
           lapply(K_SET, function(k) list(name = sprintf("intron with >= %d fragments", k), keep = pair$frag >= k)),
           lapply(K_SET, function(k) list(name = sprintf("isoform with all introns >= %d fragments", k), keep = pair$iso_min >= k)))
S2 <- as.data.table(read.xlsx(PATH_S2, sheet = 1)); coh <- U$gene_id %in% S2$gene_id
summ <- list(); per <- data.table(uid = U$uid, gene_id = U$gene_id, single_TuORF_cohort = coh, junction_status = U$junction_status)
grp <- function(p) cut(p, c(-Inf, 30, 80, Inf), labels = c("0~30", "30~80", "80~100"), right = FALSE)
for (r in rules) {
  cls <- classify(r$keep); st <- status(cls); g <- gene_rollup(st)
  if (identical(r$name, "all annotated introns (A2, manuscript)")) {
    a2c <- substr(U$annot_class_A2, 1, regexpr("\\.", U$annot_class_A2) - 1L)
    .add(sprintf("\n  check: classes with all annotated introns reproduce A2 for %d of %d TuORFs", sum(cls == a2c), nrow(U))); stopifnot(all(cls == a2c)) }
  tab <- table(factor(cls, levels = c("A", "B", "C1", "C2", "D")))
  cz <- U[coh & st != "affected by splicing"]; cz <- cz[!duplicated(gene_id)]
  ok1 <- is.finite(cz$uIF) & cz$uIF > 0; ok2 <- is.finite(cz$mORF_TE) & cz$mORF_TE > 0; gz <- grp(cz$pct_uORF_plus)
  row <- data.table(rule = r$name,
    TuORFs_affected = sum(st == "affected by splicing"), TuORFs_affected_pct = round(100 * mean(st == "affected by splicing"), 1),
    TuORFs_TSS_alone_pct = round(100 * mean(st != "affected by splicing"), 1), TuORFs_strictly_clean_pct = round(100 * mean(st == "clean (annotation + junctions)"), 1),
    genes_affected = sum(g$st == "affected by splicing"), genes_TSS_alone_pct = round(100 * mean(g$st != "affected by splicing"), 1),
    genes_strictly_clean_pct = round(100 * mean(g$st == "clean (annotation + junctions)"), 1),
    class_A = tab[["A"]], class_B = tab[["B"]], class_C1 = tab[["C1"]], class_C2 = tab[["C2"]], class_D = tab[["D"]],
    cohort_not_affected_n = nrow(cz), cohort_at_100_pct = round(100 * mean(cz$pct_uORF_plus >= 99.95), 1), cohort_below_90_pct = round(100 * mean(cz$pct_uORF_plus < 90), 1),
    cohort_groups = paste(sprintf("%.1f", 100 * as.numeric(table(gz)) / length(gz)), collapse = " / "),
    cohort_r_log2uIF = round(cor(cz$pct_uORF_plus[ok1], log2(cz$uIF[ok1])), 3), cohort_r_log2TE = round(cor(cz$pct_uORF_plus[ok2], log2(cz$mORF_TE[ok2])), 3))
  summ[[r$name]] <- row; per[, (paste0("class | ", r$name)) := cls]; per[, (paste0("status | ", r$name)) := st]
  .add(sprintf("  %-42s TuORFs affected %5d (%4.1f%%) -> TSS alone %4.1f%% (strictly clean %4.1f%%) | genes TSS alone %4.1f%% (strictly clean %4.1f%%) | classes A/B/C1/C2/D %d/%d/%d/%d/%d",
               r$name, row$TuORFs_affected, row$TuORFs_affected_pct, row$TuORFs_TSS_alone_pct, row$TuORFs_strictly_clean_pct, row$genes_TSS_alone_pct, row$genes_strictly_clean_pct,
               row$class_A, row$class_B, row$class_C1, row$class_C2, row$class_D))
  .add(sprintf("  %-42s single-TuORF genes not affected n = %d: 100%% %.1f%%, <90%% %.1f%%, groups %s, r(%%uORF+, log2 uIF) = %.3f, r(%%uORF+, log2 TE) = %.3f",
               "", row$cohort_not_affected_n, row$cohort_at_100_pct, row$cohort_below_90_pct, row$cohort_groups, row$cohort_r_log2uIF, row$cohort_r_log2TE))
}
summ <- rbindlist(summ)
# junction-only view for reference (A2): spliced over at the start or stop codon
.add(sprintf("\n  reference, junction reads only (A2): spliced over %s of TuORFs; unassessable %s", pc(U$junction_status == SPLICED), pc(U$junction_status == UNASS)))

################################################################################
# 6) Outputs
################################################################################
wb <- createWorkbook(); addWorksheet(wb, "README"); writeData(wb, "README", data.frame(x = c(
  "A14 — junction-supported leader splicing (2026-10-29; not part of the manuscript).",
  "An annotated splicing event (A2 classes: D = start codon inside an intron of another isoform; C1 = TuORF spans an intron of its own isoform; C2 = another isoform's intron alters the TuORF body/stop) counts only when the intron (or every intron of the isoform carrying it) has >= K junction fragments in the merged control RNA-seq (STAR, unique mappers, sense strand).",
  "affected by splicing = a counted annotated event, or junction reads splicing over the start or stop codon (A2: S < 0.90 with >= 10 reads); TSS alone = not affected (clean + junction-unassessable); genes: affected if any TuORF is affected.",
  "Sheets: summary (one row per rule), per_TuORF (class and status under each rule), intron_support (annotated introns of coding isoforms in TuORF genes with fragment counts), log."), stringsAsFactors = FALSE), colNames = FALSE)
addWorksheet(wb, "summary"); writeData(wb, "summary", summ)
addWorksheet(wb, "per_TuORF"); writeData(wb, "per_TuORF", per)
addWorksheet(wb, "intron_support"); writeData(wb, "intron_support", idt[gene_id %in% unique(U$gene_id), .(gene_id, tx_id, chr, start, end, strand, fragments = frag)])
addWorksheet(wb, "log"); writeData(wb, "log", data.frame(line = .log), colNames = FALSE)
saveWorkbook(wb, file.path(OUT_DIR, "A14_junction_supported_splicing.xlsx"), overwrite = TRUE)
.add(sprintf("\nA14 done in %.1f min; outputs: output/A14_junction_supported_splicing.xlsx, A14_summary.txt", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
writeLines(.log, file.path(OUT_DIR, "A14_summary.txt"))
