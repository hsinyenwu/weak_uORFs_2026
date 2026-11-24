################################################################################
# A15 — TSS signal inside 5' leader introns (Larry, 2026-10-29)                 #
# The %uORF+ window is the 5' leader exons of the TuORF's isoform + 100 nt       #
#   upstream (A13 reproduces it); TSS signal inside the leader introns of that   #
#   isoform is not counted. A TSS inside a leader intron DOWNSTREAM of the uORF  #
#   start (an alternative promoter) makes transcripts without the uORF that      #
#   count neither as TSS loss nor as splicing loss; one UPSTREAM of the uORF     #
#   start makes transcripts that keep it.                                        #
# Measured here, for CAGE-seq and TSS-seq: how much signal lies in those         #
#   introns, upstream or downstream of the uORF start, how much of the           #
#   downstream part sits at annotated 5' ends of other isoforms of the gene      #
#   (+/- 50 nt), and how %uORF+ changes when the intron signal is counted:      #
#   %uORF+ (with introns) = (window + intron signal at or upstream of the uORF   #
#   start) / (window + intron signal); 0 when there is no signal.              #
# Inputs as A13: Files/Araport11+CTRL_20181206.gtf, Files/wt_R123_genome_        #
#   RiboPlotR.txt (CAGE), Files/wt_TSS_Nielsen_2019.txt (TSS-seq), output/        #
#   Table_S1_all_TuORFs.xlsx, output/Table_S2_single_TuORF_genes.xlsx.           #
# Outputs (output/): A15_intronic_TSS_in_leaders.xlsx, A15_summary.txt.          #
# NOT part of the manuscript (results first).                                   #
# Run from Weak_uORF_submission/:  Rscript revision_modules/A15_intronic_TSS_in_leaders.R
################################################################################
suppressPackageStartupMessages({ library(GenomicFeatures); library(GenomicRanges); library(data.table); library(openxlsx) })
options(warn = 1); t0 <- Sys.time(); .log <- character(0); .add <- function(...) { m <- paste0(...); cat(m, "\n"); .log <<- c(.log, m) }
DATA_DIR <- "Files"; OUT_DIR <- "output"; UP <- 100L; NEAR <- 50L
PATH_GTF <- file.path(DATA_DIR, "Araport11+CTRL_20181206.gtf"); PATH_CAGE <- file.path(DATA_DIR, "wt_R123_genome_RiboPlotR.txt"); PATH_TSS <- file.path(DATA_DIR, "wt_TSS_Nielsen_2019.txt")
stopifnot(file.exists(PATH_GTF), file.exists(PATH_CAGE), file.exists(PATH_TSS))

# 1) transcript models
txdb <- if (requireNamespace("txdbmaker", quietly = TRUE)) txdbmaker::makeTxDbFromGFF(PATH_GTF, format = "gtf") else makeTxDbFromGFF(PATH_GTF, format = "gtf")
utr5 <- fiveUTRsByTranscript(txdb, use.names = TRUE)
tx_all <- transcripts(txdb, columns = c("tx_name", "gene_id"))
txinfo <- data.table(tx_id = mcols(tx_all)$tx_name, gene_id = as.character(unlist(mcols(tx_all)$gene_id)))
cr <- unlist(range(cdsBy(txdb, by = "tx", use.names = TRUE)))
cds_start <- setNames(ifelse(as.character(strand(cr)) == "+", start(cr), end(cr)), names(cr))
intr_by <- intronsByTranscript(txdb, use.names = TRUE)

# 2) TuORFs (as A13)
S1A <- as.data.table(read.xlsx(file.path(OUT_DIR, "Table_S1_all_TuORFs.xlsx"), sheet = "Table S1A (CAGE)"))
S1B <- as.data.table(read.xlsx(file.path(OUT_DIR, "Table_S1_all_TuORFs.xlsx"), sheet = "Table S1B (TSS-seq)"))
S2  <- as.data.table(read.xlsx(file.path(OUT_DIR, "Table_S2_single_TuORF_genes.xlsx"), sheet = 1))
U <- unique(S1A[, .(gene_id, tx_id, chr = as.character(chr), Strand, uORF_start_genomic = as.integer(uORF_start_genomic), uORF_genome_range, pUP_S1A = pct_uORF_plus)])
U[, key := paste(gene_id, uORF_genome_range, sep = "|")]

# 3) windows: the manuscript window (A13 code) and the leader introns of the same isoform
extend5 <- function(gr) {
  d <- as.data.table(gr)[, i := .I]; d[, k := names(gr)]
  d[, pos5 := ifelse(strand == "+", start, -end)]; first <- d[order(k, pos5), .SD[1], by = k]$i
  gr[first] <- resize(gr[first], width = width(gr[first]) + UP, fix = "end"); restrict(gr, start = 1L) }
has_utr <- U$tx_id %in% names(utr5)
win_iso <- utr5[U$tx_id[has_utr]]; names(win_iso) <- U$key[has_utr]; win_iso <- extend5(unlist(win_iso, use.names = TRUE))
src <- unique(U$tx_id[U$tx_id %in% names(intr_by) & U$tx_id %in% names(cds_start)])
li <- unlist(intr_by[src], use.names = TRUE); li$tx <- names(li); names(li) <- NULL
li <- li[ifelse(as.character(strand(li)) == "+", end(li) < cds_start[li$tx], start(li) > cds_start[li$tx])]            # leader introns of the source isoform
lid <- data.table(tx = li$tx, chr = as.character(seqnames(li)), start = start(li), end = end(li), strand = as.character(strand(li)))
wi <- merge(U[, .(key, tx = tx_id, s0 = uORF_start_genomic)], lid, by = "tx", allow.cartesian = TRUE)
wi[, intron_side := ifelse((strand == "+" & end < s0) | (strand == "-" & start > s0), "upstream", "downstream")]   # relative to the uORF start
win_intr <- GRanges(wi$chr, IRanges(wi$start, wi$end), strand = wi$strand); names(win_intr) <- wi$key; win_intr$side <- wi$intron_side
U[, `:=`(n_leader_introns = wi[, .N, by = key][match(U$key, key), N], n_introns_downstream = wi[intron_side == "downstream", .N, by = key][match(U$key, key), N])]
U[is.na(n_leader_introns), n_leader_introns := 0L]; U[is.na(n_introns_downstream), n_introns_downstream := 0L]
.add(sprintf("A15: %d TuORFs in %d genes; source isoform with >= 1 leader intron: %d TuORFs (%.1f%%); with a leader intron downstream of the uORF start: %d (%.1f%%)",
             nrow(U), uniqueN(U$gene_id), sum(U$n_leader_introns > 0), 100 * mean(U$n_leader_introns > 0), sum(U$n_introns_downstream > 0), 100 * mean(U$n_introns_downstream > 0)))
# annotated 5' ends of the other transcripts of the gene, inside a downstream leader intron of the source isoform
t5 <- resize(tx_all, width = 1, fix = "start"); t5d <- data.table(tx5 = mcols(tx_all)$tx_name, gene = txinfo$gene_id, chr = as.character(seqnames(t5)), pos = start(t5), strand = as.character(strand(t5)))
alt <- merge(wi[intron_side == "downstream", .(key, tx, chr, start, end, strand)], merge(U[, .(key, gene_id)], t5d, by.x = "gene_id", by.y = "gene", allow.cartesian = TRUE), by = c("key", "chr", "strand"), allow.cartesian = TRUE)
alt <- alt[tx5 != tx & pos >= start & pos <= end]
U[, alt_5end_in_downstream_intron := key %in% alt$key]
.add(sprintf("  TuORFs with an annotated 5' end of another isoform inside a downstream leader intron: %d (%.1f%%)", sum(U$alt_5end_in_downstream_intron), 100 * mean(U$alt_5end_in_downstream_intron)))

# 4) signal
read_tss <- function(path) { d <- fread(path, header = FALSE, select = 1:4, col.names = c("count", "chr", "pos", "strand")); d[, chr := sub("^[Cc]hr", "", as.character(chr))]
  gr <- GRanges(d$chr, IRanges(d$pos, width = 1), strand = d$strand, count = as.numeric(d$count)); gr[seqnames(gr) %in% seqlevels(txdb)] }
hits <- function(gr, win) { ov <- findOverlaps(gr, win); d <- data.table(kk = names(win)[subjectHits(ov)], pos = start(gr)[queryHits(ov)], count = gr$count[queryHits(ov)], strand = as.character(strand(gr))[queryHits(ov)]); setnames(d, "kk", "key") }   # "key" is a reserved argument of data.table()
grp <- function(p) cut(p, c(-Inf, 30, 80, Inf), labels = c("0~30", "30~80", "80~100"), right = FALSE)
results <- list()
for (src_name in c("CAGE", "TSSseq")) {
  gr <- read_tss(if (src_name == "CAGE") PATH_CAGE else PATH_TSS)
  hw <- merge(hits(gr, win_iso), U[, .(key, Strand, s0 = uORF_start_genomic)], by = "key")
  hi <- merge(unique(hits(gr, win_intr)), U[, .(key, Strand, s0 = uORF_start_genomic, gene_id, tx_id)], by = "key")
  up_flag <- function(d) ifelse(d$Strand == "+", d$pos <= d$s0, d$pos >= d$s0)
  hw[, up := up_flag(hw)]; hi[, up := up_flag(hi)]
  # downstream intron signal at annotated 5' ends of other isoforms of the gene (+/- NEAR nt)
  hd <- hi[up == FALSE]
  if (nrow(hd)) { nb <- merge(hd[, .(key, gene_id, tx_id, pos, count, strand, i = .I)], t5d, by.x = c("gene_id", "strand"), by.y = c("gene", "strand"), allow.cartesian = TRUE)
    near_i <- unique(nb[tx5 != tx_id & abs(pos.x - pos.y) <= NEAR, i]); hd[, near_alt5 := seq_len(.N) %in% near_i] } else hd[, near_alt5 := logical(0)]
  sw <- hw[, .(win_total = sum(count), win_up = sum(count[up])), by = key]
  si <- hi[, .(intr_total = sum(count), intr_up = sum(count[up]), intr_down = sum(count[!up])), by = key]
  sa <- hd[, .(intr_down_at_alt5 = sum(count[near_alt5])), by = key]
  R <- Reduce(function(a, b) merge(a, b, by = "key", all.x = TRUE), list(U, sw, si, sa))
  for (cn in c("win_total", "win_up", "intr_total", "intr_up", "intr_down", "intr_down_at_alt5")) R[is.na(get(cn)), (cn) := 0]
  R[, `:=`(pUP_window = ifelse(win_total > 0, 100 * win_up / win_total, 0),
           pUP_with_introns = ifelse(win_total + intr_total > 0, 100 * (win_up + intr_up) / (win_total + intr_total), 0))]
  R[, `:=`(diff = pUP_with_introns - pUP_window, intron_share = ifelse(win_total + intr_total > 0, intr_total / (win_total + intr_total), 0))]
  ref <- if (src_name == "CAGE") R$pUP_S1A else S1B$pct_uORF_plus[match(R$key, paste(S1B$gene_id, S1B$uORF_genome_range, sep = "|"))]
  ok <- is.finite(ref)
  .add(sprintf("\n== %s: check of the window re-implementation against Table S1%s: max |difference| = %.3f points over %d TuORFs", src_name, if (src_name == "CAGE") "A" else "B", max(abs(R$pUP_window - ref)[ok]), sum(ok)))
  wl <- R[n_leader_introns > 0]
  .add(sprintf("  %s signal in leader introns of the source isoform, all TuORFs: %.2f%% of window + intron signal summed over TuORFs; %.1f%% of it upstream and %.1f%% downstream of the uORF start; downstream part at annotated 5' ends of other isoforms (+/- %d nt): %.1f%%",
               src_name, 100 * sum(R$intr_total) / sum(R$win_total + R$intr_total), 100 * sum(R$intr_up) / max(sum(R$intr_total), 1), 100 * sum(R$intr_down) / max(sum(R$intr_total), 1), NEAR, 100 * sum(R$intr_down_at_alt5) / max(sum(R$intr_down), 1)))
  .add(sprintf("  %s TuORFs whose isoform has a leader intron (n = %d): any intron signal %.1f%%; intron share of the signal median %.1f%%, > 10%% in %.1f%%, > 25%% in %.1f%%",
               src_name, nrow(wl), 100 * mean(wl$intr_total > 0), 100 * median(wl$intron_share), 100 * mean(wl$intron_share > 0.10), 100 * mean(wl$intron_share > 0.25)))
  .add(sprintf("  %s all TuORFs: Spearman(window, with introns) = %.4f; |change| > 5 points %d (%.1f%%), > 10 %d (%.1f%%), > 20 %d (%.1f%%); drops > 10 points %d, rises > 10 points %d; mean change %+.2f points",
               src_name, cor(R$pUP_window, R$pUP_with_introns, method = "spearman"), sum(abs(R$diff) > 5), 100 * mean(abs(R$diff) > 5), sum(abs(R$diff) > 10), 100 * mean(abs(R$diff) > 10),
               sum(abs(R$diff) > 20), 100 * mean(abs(R$diff) > 20), sum(R$diff < -10), sum(R$diff > 10), mean(R$diff)))
  tab <- table(window = grp(R$pUP_window), with_introns = grp(R$pUP_with_introns)); .add(paste0("  group transitions, all TuORFs (rows window, columns with introns):\n", paste(capture.output(print(tab)), collapse = "\n")))
  C <- R[gene_id %in% S2$gene_id]; C <- merge(C, S2[, .(gene_id, uIF = as.numeric(uIF), TE = as.numeric(mORF_TE))], by = "gene_id", all.x = TRUE); C <- C[!duplicated(gene_id)]
  for (w in c("window", "with_introns")) { p <- C[[paste0("pUP_", w)]]; ok1 <- is.finite(C$uIF) & C$uIF > 0; ok2 <- is.finite(C$TE) & C$TE > 0; g <- grp(p)
    med <- tapply(C$TE[ok2], g[ok2], median); ks <- function(a, b) tryCatch(suppressWarnings(ks.test(C$TE[ok2 & g == a], C$TE[ok2 & g == b])$p.value), error = function(e) NA_real_)
    .add(sprintf("  %s single-TuORF cohort (n = %d), %s: at 100%% %.1f%%; below 90%% %.1f%%; groups 0~30 / 30~80 / 80~100 = %.1f / %.1f / %.1f%%; r(%%uORF+, log2 uIF) = %.3f; r(%%uORF+, log2 TE) = %.3f; median TE %.2f / %.2f / %.2f; KS p 0~30 vs 80~100 %.2g",
                 src_name, nrow(C), w, 100 * mean(p >= 99.95), 100 * mean(p < 90), 100 * mean(g == "0~30"), 100 * mean(g == "30~80"), 100 * mean(g == "80~100"),
                 cor(p[ok1], log2(C$uIF[ok1])), cor(p[ok2], log2(C$TE[ok2])), med["0~30"], med["30~80"], med["80~100"], ks("0~30", "80~100"))) }
  .add(sprintf("  %s cohort: genes changing group %d (%.1f%%); Spearman(window, with introns) = %.4f", src_name, sum(grp(C$pUP_window) != grp(C$pUP_with_introns)), 100 * mean(grp(C$pUP_window) != grp(C$pUP_with_introns)), cor(C$pUP_window, C$pUP_with_introns, method = "spearman")))
  results[[src_name]] <- list(all = R, cohort = C)
}

# 5) outputs
wb <- createWorkbook(); addWorksheet(wb, "README"); writeData(wb, "README", data.frame(x = c(
  "A15 — TSS signal inside the 5' leader introns of the TuORF's isoform (2026-10-29; not part of the manuscript).",
  "window = manuscript %uORF+ window (5' leader exons of the TuORF's isoform + 100 nt upstream); introns = leader introns of the same isoform (upstream of the CDS start), split into upstream / downstream of the uORF start codon.",
  sprintf("pUP_with_introns = 100 x (window + intron signal at or upstream of the uORF start) / (window + intron signal); intr_down_at_alt5 = downstream intron signal within %d nt of an annotated 5' end of another isoform of the gene.", NEAR),
  "Sheets: TuORFs_CAGE / TuORFs_TSSseq (all TuORFs), cohort_CAGE / cohort_TSSseq (single-TuORF genes with uIF and TE), log."), stringsAsFactors = FALSE), colNames = FALSE)
for (s in names(results)) { addWorksheet(wb, paste0("TuORFs_", s)); writeData(wb, paste0("TuORFs_", s), results[[s]]$all[order(gene_id)]); addWorksheet(wb, paste0("cohort_", s)); writeData(wb, paste0("cohort_", s), results[[s]]$cohort[order(gene_id)]) }
addWorksheet(wb, "log"); writeData(wb, "log", data.frame(line = .log), colNames = FALSE)
saveWorkbook(wb, file.path(OUT_DIR, "A15_intronic_TSS_in_leaders.xlsx"), overwrite = TRUE)
.add(sprintf("\nA15 done in %.1f min; outputs: output/A15_intronic_TSS_in_leaders.xlsx, A15_summary.txt", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
writeLines(.log, file.path(OUT_DIR, "A15_summary.txt"))
