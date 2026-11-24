################################################################################
# A13 — Sensitivity of %uORF+ to the TSS window (Larry, 2026-10-28)              #
# Current definition (manuscript): the annotated 5' leader exons of the TuORF's  #
#   isoform + 100 nt upstream of that isoform's annotated 5' end.                #
# Alternative tested here: the UNION of the 5' leader exons of ALL isoforms of   #
#   the gene whose leader contains the TuORF start codon, extended 100 nt        #
#   upstream of the most upstream 5' end; it ends at the most downstream mORF   #
#   start among those isoforms (no CDS / 3' UTR signal).                          #
# Same formula as the pipeline: %uORF+ = 100 x signal at or upstream of the      #
#   uORF start / all signal in the window (0 when the window has no signal).     #
# Inputs: Files/Araport11+CTRL_20181206.gtf (TxDb), Files/wt_R123_genome_        #
#   RiboPlotR.txt (CAGE), Files/wt_TSS_Nielsen_2019.txt (TSS-seq), output/        #
#   Table_S1_all_TuORFs.xlsx (S1A), output/Table_S2_single_TuORF_genes.xlsx,     #
#   output/A12_Table_genes.csv (No-TuORF TE for the ECDF panel).                 #
# Outputs (output/): A13_union_window_sensitivity.xlsx, A13_Figure_union_window  #
#   _sensitivity.pdf, A13_summary.txt, A13_run.log. NOT part of the manuscript.  #
# Run from Weak_uORF_submission/:  Rscript revision_modules/A13_union_leader_window_sensitivity.R
################################################################################
suppressPackageStartupMessages({ library(GenomicFeatures); library(GenomicRanges); library(GenomeInfoDb); library(data.table); library(openxlsx); library(ggplot2); library(gridExtra); library(grid) })
options(warn = 1); t0 <- Sys.time(); .log <- character(0); .add <- function(...) { m <- paste0(...); cat(m, "\n"); .log <<- c(.log, m) }
DATA_DIR <- "Files"; OUT_DIR <- "output"; UP <- 100L
PATH_GTF <- file.path(DATA_DIR, "Araport11+CTRL_20181206.gtf"); PATH_CAGE <- file.path(DATA_DIR, "wt_R123_genome_RiboPlotR.txt"); PATH_TSS <- file.path(DATA_DIR, "wt_TSS_Nielsen_2019.txt")
stopifnot(file.exists(PATH_GTF), file.exists(PATH_CAGE), file.exists(PATH_TSS))
source("revision_modules/fig_compact.R")

# 1) transcript models
txdb <- if (requireNamespace("txdbmaker", quietly = TRUE)) txdbmaker::makeTxDbFromGFF(PATH_GTF, format = "gtf") else makeTxDbFromGFF(PATH_GTF, format = "gtf")
utr5 <- fiveUTRsByTranscript(txdb, use.names = TRUE)
txinfo <- as.data.table(transcripts(txdb, columns = c("tx_name", "gene_id")))[, .(tx_id = tx_name, gene_id = as.character(gene_id))]
.add(sprintf("A13: TxDb %d transcripts, %d with a 5' leader; %.1f min", nrow(txinfo), length(utr5), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# 2) TuORFs
S1A <- as.data.table(read.xlsx(file.path(OUT_DIR, "Table_S1_all_TuORFs.xlsx"), sheet = "Table S1A (CAGE)"))
S1B <- as.data.table(read.xlsx(file.path(OUT_DIR, "Table_S1_all_TuORFs.xlsx"), sheet = "Table S1B (TSS-seq)"))
S2  <- as.data.table(read.xlsx(file.path(OUT_DIR, "Table_S2_single_TuORF_genes.xlsx"), sheet = 1))
U <- unique(S1A[, .(gene_id, tx_id, chr = as.character(chr), Strand, uORF_start_genomic = as.integer(uORF_start_genomic), uORF_genome_range, pUP_iso_S1A = pct_uORF_plus)])
U[, key := paste(gene_id, uORF_genome_range, sep = "|")]
.add(sprintf("  TuORFs: %d (unique gene|range), %d genes; listed isoforms in the TxDb: %d; MSTRG isoforms among them: %d", nrow(U), uniqueN(U$gene_id), sum(U$tx_id %in% names(utr5)), sum(grepl("^MSTRG", U$tx_id))))

# 3) isoforms whose 5' leader contains the TuORF start codon (first nucleotide), among all transcripts of the gene (+ the listed isoform)
starts <- GRanges(U$chr, IRanges(U$uORF_start_genomic, width = 1), strand = U$Strand); names(starts) <- U$key
utr_ul <- unlist(utr5, use.names = TRUE); utr_ul$tx <- names(utr_ul); utr_ul$gene <- txinfo$gene_id[match(utr_ul$tx, txinfo$tx_id)]
ov <- findOverlaps(starts, utr_ul)
cand <- data.table(kk = names(starts)[queryHits(ov)], tx = utr_ul$tx[subjectHits(ov)], gene = utr_ul$gene[subjectHits(ov)]); setnames(cand, "kk", "key")   # "key" is a reserved argument of data.table()
cand <- merge(cand, U[, .(key, gene_id, tx_id)], by = "key")
cand <- cand[gene == gene_id | tx == tx_id]                                   # same gene (or the listed isoform itself)
iso <- rbind(cand[, .(key, tx)], U[, .(key, tx = tx_id)])[, unique(.SD)]      # always include the listed isoform
iso <- iso[tx %in% names(utr5)]
n_iso <- iso[, .(n_iso = .N), by = key]; U <- merge(U, n_iso, by = "key", all.x = TRUE); U[is.na(n_iso), n_iso := 0L]
gene_ntx <- txinfo[, .(n_tx_gene = .N), by = gene_id]; U <- merge(U, gene_ntx, by = "gene_id", all.x = TRUE)
.add(sprintf("  isoforms containing the start codon: 1 for %d TuORFs (%.1f%%), 2 for %d, 3+ for %d; genes with > 1 annotated transcript: %d of %d TuORF genes",
             sum(U$n_iso == 1), 100 * mean(U$n_iso == 1), sum(U$n_iso == 2), sum(U$n_iso >= 3), uniqueN(U[n_tx_gene > 1, gene_id]), uniqueN(U$gene_id)))

# 4) windows: isoform (current) and union; each a GRanges with names = key (possibly several ranges per key)
extend5 <- function(gr) {                                  # extend the 5'-most range of each key by UP nt upstream (strand-aware)
  d <- as.data.table(gr)[, i := .I]; d[, k := names(gr)]
  d[, pos5 := ifelse(strand == "+", start, -end)]; first <- d[order(k, pos5), .SD[1], by = k]$i
  gr[first] <- resize(gr[first], width = width(gr[first]) + UP, fix = "end"); restrict(gr, start = 1L) }
win_iso <- utr5[U$tx_id[U$tx_id %in% names(utr5)]]; names(win_iso) <- U$key[U$tx_id %in% names(utr5)]
win_iso <- unlist(win_iso, use.names = TRUE); win_iso <- extend5(win_iso)
ul <- utr5[iso$tx]; names(ul) <- iso$key; ul <- unlist(ul, use.names = TRUE)
win_union <- unlist(reduce(split(ul, names(ul))), use.names = TRUE)        # union of leader exons per key
win_union <- extend5(win_union)
wid <- function(w) as.data.table(w)[, .(width = sum(width)), by = .(key = names(w))]
U <- merge(U, setnames(wid(win_iso), "width", "window_iso_nt"), by = "key", all.x = TRUE); U <- merge(U, setnames(wid(win_union), "width", "window_union_nt"), by = "key", all.x = TRUE)
U[, window_differs := !is.na(window_union_nt) & !is.na(window_iso_nt) & window_union_nt > window_iso_nt]
.add(sprintf("  union window wider than the isoform window for %d TuORFs (%.1f%%); median extra width %d nt", sum(U$window_differs), 100 * mean(U$window_differs), as.integer(median(U[window_differs == TRUE, window_union_nt - window_iso_nt]))))

# 5) %uORF+ against a TSS file for a set of windows
read_tss <- function(path) { d <- fread(path, header = FALSE, select = 1:4, col.names = c("count", "chr", "pos", "strand")); d[, chr := sub("^[Cc]hr", "", as.character(chr))]
  gr <- GRanges(d$chr, IRanges(d$pos, width = 1), strand = d$strand, count = as.numeric(d$count)); gr[seqnames(gr) %in% seqlevels(txdb)] }
pup <- function(gr, win, U) {
  ov <- findOverlaps(gr, win); dt <- data.table(kk = names(win)[subjectHits(ov)], pos = start(gr)[queryHits(ov)], count = gr$count[queryHits(ov)]); setnames(dt, "kk", "key")
  dt <- merge(dt, U[, .(key, Strand, s0 = uORF_start_genomic)], by = "key")
  r <- dt[, .(total = sum(count), up = sum(count[if (Strand[1] == "+") pos <= s0[1] else pos >= s0[1]])), by = key]   # ifelse() would return a length-1 mask
  r[, pUP := ifelse(total > 0, 100 * up / total, 0)]; r }
grp <- function(p) cut(p, c(-Inf, 30, 80, Inf), labels = c("0~30", "30~80", "80~100"), right = FALSE)
results <- list()
for (src in c("CAGE", "TSSseq")) {
  gr <- read_tss(if (src == "CAGE") PATH_CAGE else PATH_TSS)
  ri <- pup(gr, win_iso, U); ru <- pup(gr, win_union, U)
  R <- merge(U, setnames(ri, c("key", "signal_iso", "up_iso", "pUP_iso")), by = "key", all.x = TRUE)
  R <- merge(R, setnames(ru, c("key", "signal_union", "up_union", "pUP_union")), by = "key", all.x = TRUE)
  for (cn in c("signal_iso", "up_iso", "pUP_iso", "signal_union", "up_union", "pUP_union")) R[is.na(get(cn)), (cn) := 0]
  ref <- if (src == "CAGE") U$pUP_iso_S1A else S1B$pct_uORF_plus[match(U$key, paste(S1B$gene_id, S1B$uORF_genome_range, sep = "|"))]
  R[, pUP_table := ref]
  R[, `:=`(diff = pUP_union - pUP_iso, group_iso = grp(pUP_iso), group_union = grp(pUP_union))]
  ok <- is.finite(R$pUP_table)
  .add(sprintf("\n== %s: check of the re-implementation against Table S1%s: max |pUP_iso - table| = %.3f points over %d TuORFs (%d differ by > 0.5)", src, if (src == "CAGE") "A" else "B", max(abs(R$pUP_iso - R$pUP_table)[ok]), sum(ok), sum(abs(R$pUP_iso - R$pUP_table)[ok] > 0.5)))
  ch <- R[window_differs == TRUE]
  .add(sprintf("  %s, all %d TuORFs: Spearman(pUP_iso, pUP_union) = %.4f; identical value for %d (%.1f%%); |change| > 5 points for %d (%.1f%%), > 10 for %d (%.1f%%), > 20 for %d (%.1f%%); mean change %+.2f points",
               src, nrow(R), cor(R$pUP_iso, R$pUP_union, method = "spearman"), sum(abs(R$diff) < 1e-6), 100 * mean(abs(R$diff) < 1e-6), sum(abs(R$diff) > 5), 100 * mean(abs(R$diff) > 5), sum(abs(R$diff) > 10), 100 * mean(abs(R$diff) > 10), sum(abs(R$diff) > 20), 100 * mean(abs(R$diff) > 20), mean(R$diff)))
  .add(sprintf("  %s, the %d TuORFs whose window widens: median change %+.1f points, mean %+.1f; rises > 10 points %d, drops > 10 points %d", src, nrow(ch), median(ch$diff), mean(ch$diff), sum(ch$diff > 10), sum(ch$diff < -10)))
  tab <- table(iso = R$group_iso, union = R$group_union); .add(paste0("  group transitions (rows isoform window, columns union window):\n", paste(capture.output(print(tab)), collapse = "\n")))
  # single-TuORF cohort (Table S3 genes): the manuscript's numbers under both windows
  C <- R[gene_id %in% S2$gene_id]; C <- merge(C, S2[, .(gene_id, uIF = as.numeric(uIF), TE = as.numeric(mORF_TE))], by = "gene_id", all.x = TRUE)
  if (!"mORF_TE" %in% names(S2)) C[, TE := as.numeric(S2$TE[match(C$gene_id, S2$gene_id)])]
  C <- C[!duplicated(gene_id)]
  for (w in c("iso", "union")) { p <- C[[paste0("pUP_", w)]]
    .add(sprintf("  %s single-TuORF cohort (n = %d), %s window: at 100%%: %d (%.1f%%); below 100%%: %.1f%%; below 90%%: %.1f%%; groups 0~30 / 30~80 / 80~100: %.1f / %.1f / %.1f%%", src, nrow(C), w, sum(p >= 99.95), 100 * mean(p >= 99.95), 100 * mean(p < 99.95), 100 * mean(p < 90), 100 * mean(p < 30), 100 * mean(p >= 30 & p < 80), 100 * mean(p >= 80)))
    ok1 <- is.finite(C$uIF) & C$uIF > 0; ok2 <- is.finite(C$TE) & C$TE > 0
    .add(sprintf("     r(log2 uIF, %%uORF+) = %.3f (n = %d); r(log2 TE, %%uORF+) = %.3f (n = %d); Spearman %.3f / %.3f", cor(log2(C$uIF[ok1]), p[ok1]), sum(ok1), cor(log2(C$TE[ok2]), p[ok2]), sum(ok2), cor(C$uIF[ok1], p[ok1], method = "spearman"), cor(C$TE[ok2], p[ok2], method = "spearman")))
    g <- grp(p); med <- tapply(C$TE[ok2], g[ok2], median); ks <- function(a, b) tryCatch(suppressWarnings(ks.test(C$TE[ok2 & g == a], C$TE[ok2 & g == b])$p.value), error = function(e) NA_real_)
    .add(sprintf("     median mORF TE by group: 0~30 %.2f (n = %d), 30~80 %.2f (n = %d), 80~100 %.2f (n = %d); KS p: 0~30 vs 30~80 %.2g, 30~80 vs 80~100 %.2g, 0~30 vs 80~100 %.2g",
                 med["0~30"], sum(g[ok2] == "0~30"), med["30~80"], sum(g[ok2] == "30~80"), med["80~100"], sum(g[ok2] == "80~100"), ks("0~30", "30~80"), ks("30~80", "80~100"), ks("0~30", "80~100"))) }
  .add(sprintf("  %s cohort: Spearman(pUP_iso, pUP_union) = %.4f; genes changing %%uORF+ group: %d (%.1f%%); genes whose window widens: %d", src, cor(C$pUP_iso, C$pUP_union, method = "spearman"), sum(C$group_iso != C$group_union), 100 * mean(C$group_iso != C$group_union), sum(C$window_differs)))
  results[[src]] <- list(all = R, cohort = C)
}

# 6) figure (CAGE, single-TuORF cohort): A scatter of the two %uORF+ values, B %uORF+ (union) distribution, C TE ECDF by union groups, D uIF vs %uORF+ (union)
C <- results$CAGE$cohort; A12 <- tryCatch(fread(file.path(OUT_DIR, "A12_Table_genes.csv")), error = function(e) NULL)
noTu <- if (!is.null(A12) && all(c("cohort", "TE") %in% names(A12))) A12[cohort == "No TuORF" & is.finite(TE) & TE > 0, TE] else numeric(0)
pA <- ggplot(C, aes(pUP_iso, pUP_union)) + geom_point(alpha = 0.25, size = 0.7, colour = "#2166AC") + geom_abline(linetype = 2) + theme_classic() +
  labs(title = "A  %uORF+: union of isoform leaders vs the TuORF's isoform", subtitle = sprintf("single-TuORF genes, CAGE-seq; Spearman %.3f; %d genes change value, %d change group", cor(C$pUP_iso, C$pUP_union, method = "spearman"), sum(abs(C$diff) > 1e-6), sum(C$group_iso != C$group_union)), x = "%uORF+ (isoform leader + 100 nt; manuscript)", y = "%uORF+ (union of leaders + 100 nt)")
gshare <- C[, .(n = .N), by = .(g = grp(pUP_union))][, pct := 100 * n / sum(n)]
pB <- ggplot(C, aes(pUP_union)) + geom_histogram(binwidth = 2, fill = "#3498db", colour = "white", linewidth = 0.1) + geom_vline(xintercept = c(30, 80), colour = "black", linewidth = 0.4) + theme_classic() +
  annotate("text", x = c(15, 55, 90), y = Inf, vjust = 1.3, size = 8 / .pt, label = sprintf("%.1f%%\n(%d)", gshare$pct[match(c("0~30", "30~80", "80~100"), gshare$g)], gshare$n[match(c("0~30", "30~80", "80~100"), gshare$g)])) +
  labs(title = "B  %uORF+ with the union window", subtitle = sprintf("%d genes; at 100%%: %.1f%%; below 90%%: %.1f%%", nrow(C), 100 * mean(C$pUP_union >= 99.95), 100 * mean(C$pUP_union < 90)), x = "%uORF+ (CAGE-seq, union window)", y = "Number of genes")
dE <- rbind(C[is.finite(TE) & TE > 0, .(TE, g = as.character(grp(pUP_union)))], if (length(noTu)) data.table(TE = noTu, g = "No TuORF") else NULL)
dE[, g := factor(g, c("80~100", "30~80", "0~30", "No TuORF"))]; ng <- dE[, .N, by = g]
pC <- ggplot(dE[TE <= 3.5], aes(TE, colour = g)) + stat_ecdf(linewidth = 0.8) + scale_colour_manual(values = c("80~100" = "#F8766D", "30~80" = "#00BA38", "0~30" = "#619CFF", "No TuORF" = "grey40"), labels = sprintf("%s (n = %d)", ng$g, ng$N)[match(c("80~100", "30~80", "0~30", "No TuORF"), ng$g)], name = "%uORF+ group (union window)") +
  theme_classic() + theme(legend.position = c(0.68, 0.35)) + labs(title = "C  mORF TE by %uORF+ group, union window", subtitle = "as Figure 2G", x = "mORF TE", y = "Fraction")
ok1 <- is.finite(C$uIF) & C$uIF > 0
pD <- ggplot(C[ok1], aes(pUP_union, log2(uIF))) + geom_point(alpha = 0.25, size = 0.7, colour = "#2166AC") + geom_smooth(method = "lm", se = FALSE, colour = "red", linewidth = 0.8) + theme_classic() +
  labs(title = "D  uIF vs %uORF+, union window", subtitle = sprintf("as Figure 2E; r = %.2f (isoform window r = %.2f), n = %d", cor(C$pUP_union[ok1], log2(C$uIF[ok1])), cor(C$pUP_iso[ok1], log2(C$uIF[ok1])), sum(ok1)), x = "%uORF+ (CAGE-seq, union window)", y = bquote(log[2]~"(uIF)"))
g <- arrangeGrob(grobs = lapply(list(pA, pB, pC, pD), compact), ncol = 2, top = fig_title("A13. Sensitivity of %uORF+ to the TSS window (not part of the manuscript)"))
fig_save(file.path(OUT_DIR, "A13_Figure_union_window_sensitivity.pdf"), g, h = 6.5)

# 7) tables
wb <- createWorkbook(); addWorksheet(wb, "README"); writeData(wb, "README", data.frame(x = c(
  "A13 — sensitivity of %uORF+ to the TSS window (2026-10-28; not part of the manuscript).",
  "isoform window (manuscript): annotated 5' leader exons of the TuORF's isoform + 100 nt upstream of its 5' end; union window: union of the 5' leader exons of all isoforms of the gene whose leader contains the TuORF start codon, + 100 nt upstream of the most upstream 5' end (ends at the most downstream mORF start among them).",
  "%uORF+ = 100 x signal at or upstream of the uORF start / all signal in the window (0 if no signal); CAGE-seq = wt_R123_genome_RiboPlotR.txt, TSS-seq = wt_TSS_Nielsen_2019.txt.",
  "Sheets: TuORFs_CAGE / TuORFs_TSSseq (all 7,423 TuORFs: isoform counts, window widths, signal, %uORF+ under both windows, groups), cohort_CAGE / cohort_TSSseq (single-TuORF genes with uIF and TE), summary (log)."), stringsAsFactors = FALSE), colNames = FALSE)
for (src in names(results)) { addWorksheet(wb, paste0("TuORFs_", src)); writeData(wb, paste0("TuORFs_", src), results[[src]]$all[order(gene_id)]); addWorksheet(wb, paste0("cohort_", src)); writeData(wb, paste0("cohort_", src), results[[src]]$cohort[order(gene_id)]) }
addWorksheet(wb, "summary"); writeData(wb, "summary", data.frame(line = .log), colNames = FALSE)
saveWorkbook(wb, file.path(OUT_DIR, "A13_union_window_sensitivity.xlsx"), overwrite = TRUE)
.add(sprintf("A13 done in %.1f min; outputs: output/A13_union_window_sensitivity.xlsx, A13_Figure_union_window_sensitivity.pdf, A13_summary.txt", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
writeLines(.log, file.path(OUT_DIR, "A13_summary.txt"))
