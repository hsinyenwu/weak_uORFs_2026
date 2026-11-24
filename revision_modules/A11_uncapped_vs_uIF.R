################################################################################
# A11 — Uncapped mRNA (GMUCT) and mRNA half-life versus uIF and mORF TE          #
# (Reviewer 3: "compare percent uncapped from the degradome data of Yu et al.    #
#  2021 with uIF as in Figure 1B, and across the uIF groups as in Figure 1C")    #
#                                                                                #
# GMUCT and total RNA-seq of Col-0 12-day-old seedlings without ABA from Yu et   #
# al. 2021 Dev Cell 56:125-140 (GEO GSE142388, part of GSE142390): the deposited  #
# HTSeq gene counts (STAR, TAIR10) of GMUCT rep1/2 (GSM4226794/GSM4226795) and   #
# total RNA-seq rep1/2 (GSM4226810/GSM4226811; GEO titles Col0_NT_total_rep1/2). #
# Uncapped ratio U = GMUCT RPM / RNA-seq RPM per gene (= 2^"proportion uncapped" #
# of Yu et al. 2021, who report log2[RPM_GMUCT / RPM_RNA-seq]).                  #
# mRNA half-lives: 5-EU labelling data set (Files/RNA-half-lives-5EU.tsv; mean   #
# of three replicates), as used in the canonical script.                         #
#                                                                                #
# Cohorts: single-TuORF genes (Table S3; uIF classes as Figure 1B) and the       #
# No-TuORF control genes of Figure 1C / A5 (annotated 5' leader, RNA TPM >       #
# 1, TE defined, no TuORF and no ouORF).                                         #
# Outputs (output/): A11_Figure_uncapped_vs_uIF.pdf (Figure S19),               #
#   A11_Table_uncapped_halflife.xlsx (Table S10), A11_Table_genes.csv,                     #
#   A11_summary.txt, A11_run.log                                                 #
# Run from Weak_uORF_submission/:  Rscript revision_modules/A11_uncapped_vs_uIF.R#
################################################################################
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(ggplot2); library(openxlsx); library(gridExtra); library(grid); library(multcompView)
})
options(stringsAsFactors = FALSE, warn = 1, dplyr.summarise.inform = FALSE)
set.seed(11)
t0 <- Sys.time()
DATA_DIR <- "Files"; OUT_DIR <- "output"
GEO_DIR  <- file.path(DATA_DIR, "GSE142388_seedling_GMUCT", "raw")
GEO <- c(GMUCT_rep1 = "GSM4226794_Col_NT_GMUCT_rep1_gene_Count.txt.gz", GMUCT_rep2 = "GSM4226795_Col_NT_GMUCT_rep2_gene_Count.txt.gz",
         RNA_rep1   = "GSM4226810_16054X9geneCount.txt.gz",             RNA_rep2   = "GSM4226811_16054X10geneCount.txt.gz")
PATH_S1    <- file.path(OUT_DIR, "Table_S1_all_TuORFs.xlsx")
S2_CANDIDATES <- c(Sys.getenv("A11_S2", ""), "/Users/wu/Desktop/weak_uORF_revision/revision_tables/TableS3_single_TuORF_genes_revised.xlsx", "/Users/wu/Desktop/weak_uORF_revision/revision_tables/TableS2_single_TuORF_genes_revised.xlsx")
PATH_S2    <- S2_CANDIDATES[nzchar(S2_CANDIDATES) & file.exists(S2_CANDIDATES)][1]
PATH_OUORF <- file.path(DATA_DIR, "CiPS_ouORF_all_Final.xlsx")
PATH_RSEM_RNA  <- file.path(DATA_DIR, "RSEM_Araport11_1st+2nd_Riboseq_RNA.genes.results_CDS_only")
PATH_RSEM_RIBO <- file.path(DATA_DIR, "RSEM_Araport11_1st+2nd_Riboseq_Ribo.genes.results_CDS_only")
PATH_5UTR  <- file.path(DATA_DIR, "genes_w_5UTR.tsv")
PATH_HL    <- file.path(DATA_DIR, "RNA-half-lives-5EU.tsv")
for (f in c(file.path(GEO_DIR, GEO), PATH_S1, PATH_S2, PATH_OUORF, PATH_RSEM_RNA, PATH_RSEM_RIBO, PATH_5UTR, PATH_HL)) stopifnot(file.exists(f))
MIN_RNA <- 10L            # RNA-seq reads (both replicates summed) required for a gene's uncapped ratio
CUTOFFS <- c(5L, 10L, 20L); NBOOT <- 1000L
CLASSES <- c("High uIF (> 1)", "Medium uIF (0.2 to 1)", "Low uIF (0.2 or less)", "No TuORF")
PUP_CLASSES <- c("0 to 30% uORF+", "30 to 80% uORF+", "80 to 100% uORF+", "No TuORF")
.log <- character(0); .add <- function(...) { s <- paste0(...); message(s); .log <<- c(.log, s) }
big <- function(x) trimws(format(x, big.mark = ",")); f2 <- function(x, d = 2) formatC(x, format = "f", digits = d)

################################################################################
# 1) GEO counts -> RPM, uncapped ratio U
################################################################################
.add("A11: GEO counts (", GEO_DIR, ")")
rd <- function(f) { d <- fread(cmd = paste("gzip -dc", shQuote(file.path(GEO_DIR, f))), header = FALSE, col.names = c("gene_id", "n")); d[!grepl("^__", gene_id)] }
cnt <- Reduce(function(a, b) merge(a, b, by = "gene_id", all = TRUE), Map(function(nm, f) setnames(rd(f), "n", nm), names(GEO), GEO))
for (j in names(GEO)) set(cnt, which(is.na(cnt[[j]])), j, 0L)
nuclear <- grepl("^AT[1-5]G", cnt$gene_id)
libsize <- sapply(names(GEO), function(j) sum(cnt[[j]][nuclear]))       # reads assigned to nuclear genes (chromosomes 1-5)
for (j in names(GEO)) set(cnt, j = paste0("rpm_", j), value = cnt[[j]] / libsize[[j]] * 1e6)
cnt[, `:=`(GMUCT_reads = GMUCT_rep1 + GMUCT_rep2, RNA_reads = RNA_rep1 + RNA_rep2,
           GMUCT_rpm = (rpm_GMUCT_rep1 + rpm_GMUCT_rep2) / 2, RNA_rpm = (rpm_RNA_rep1 + rpm_RNA_rep2) / 2)]
cnt[, U := GMUCT_rpm / RNA_rpm]
cnt[, U_rep1 := ifelse(RNA_rep1 > 0 & GMUCT_rep1 > 0, rpm_GMUCT_rep1 / rpm_RNA_rep1, NA_real_)]
cnt[, U_rep2 := ifelse(RNA_rep2 > 0 & GMUCT_rep2 > 0, rpm_GMUCT_rep2 / rpm_RNA_rep2, NA_real_)]
.add(sprintf("  genes in the count tables: %s; nuclear-gene reads: GMUCT %s / %s, RNA-seq %s / %s (the RNA-seq libraries are total RNA; %.0f%% of their assigned reads are organellar)",
             big(nrow(cnt)), big(libsize[1]), big(libsize[2]), big(libsize[3]), big(libsize[4]),
             100 * (1 - (libsize[3] + libsize[4]) / (sum(cnt$RNA_rep1) + sum(cnt$RNA_rep2)))))
nuc <- cnt[nuclear]
.add(sprintf("  nuclear genes with >= %d RNA-seq reads: %s of %s; replicate agreement among them (Spearman): GMUCT RPM %.3f, RNA-seq RPM %.3f, uncapped ratio %.3f",
             MIN_RNA, big(nuc[RNA_reads >= MIN_RNA, .N]), big(nrow(nuc)),
             nuc[RNA_reads >= MIN_RNA, cor(rpm_GMUCT_rep1, rpm_GMUCT_rep2, method = "spearman")],
             nuc[RNA_reads >= MIN_RNA, cor(rpm_RNA_rep1, rpm_RNA_rep2, method = "spearman")],
             nuc[RNA_reads >= MIN_RNA & is.finite(U_rep1) & is.finite(U_rep2), cor(U_rep1, U_rep2, method = "spearman")]))

################################################################################
# 2) Cohorts: single-TuORF genes (Table S3) and the No-TuORF controls of Figure 1C
################################################################################
S2 <- as.data.table(read.xlsx(PATH_S2, 1))
S2 <- S2[, .(gene_id, uIF = ifelse(is.finite(mORF_TL) & mORF_TL > 0, uORF_TL / mORF_TL, uIF), uIF_Group, pUP = pct_uORF_plus, pUP_Group, RNA_TPM, Ribo_TPM, TE = mORF_TE)]   # unrounded uIF (2026-11-01)
S2[, class := fcase(grepl("^(High|Strong)", uIF_Group), CLASSES[1], grepl("^Medium", uIF_Group), CLASSES[2], grepl("^(Low|Weak)", uIF_Group), CLASSES[3])]   # Table S2 files carry either the old (Strong/Weak) or the renamed (High/Low) labels
S2[, pclass := fcase(pUP_Group == "0~30", PUP_CLASSES[1], pUP_Group == "30~80", PUP_CLASSES[2], pUP_Group == "80~100", PUP_CLASSES[3])]
stopifnot(!anyNA(S2$class), !anyNA(S2$pclass))
S1A <- read.xlsx(PATH_S1, sheet = "Table S1A (CAGE)"); ou <- read.xlsx(PATH_OUORF)
rna <- fread(PATH_RSEM_RNA, select = c("gene_id", "TPM")); ribo <- fread(PATH_RSEM_RIBO, select = c("gene_id", "TPM"))
ctrl <- data.table(gene_id = rna$gene_id, RNA_TPM = rna$TPM, Ribo_TPM = ribo$TPM[match(rna$gene_id, ribo$gene_id)])
ctrl[, TE := Ribo_TPM / RNA_TPM]
ctrl <- ctrl[gene_id %in% read.delim(PATH_5UTR)$gene_id & RNA_TPM > 1 & is.finite(TE) & !gene_id %in% unique(c(S1A$gene_id, as.character(ou$gene_id)))]
ctrl[, `:=`(uIF = NA_real_, uIF_Group = NA_character_, pUP = NA_real_, pUP_Group = NA_character_, class = CLASSES[4], pclass = PUP_CLASSES[4])]
G <- rbind(S2[, cohort := "single-TuORF"], ctrl[, cohort := "No TuORF"], fill = TRUE)
hl <- fread(PATH_HL); setnames(hl, 1, "gene_id"); hl <- hl[, .(gene_id, half_life = Mean)]
G <- merge(G, hl, by = "gene_id", all.x = TRUE)
G <- merge(G, cnt[, .(gene_id, GMUCT_rep1, GMUCT_rep2, RNA_rep1, RNA_rep2, GMUCT_reads, RNA_reads, GMUCT_rpm, RNA_rpm, U, U_rep1, U_rep2)], by = "gene_id", all.x = TRUE)
G[, U_labRNA := GMUCT_rpm / RNA_TPM]                                     # sensitivity: our own (deep) RNA-seq as denominator
G[, class := factor(class, CLASSES)]; G[, pclass := factor(pclass, PUP_CLASSES)]
.add(sprintf("  cohorts: %s single-TuORF genes (Table S2; %d / %d / %d high / medium / low uIF) and %s No-TuORF control genes; %s have GEO counts, %s have a half-life",
             big(sum(G$cohort == "single-TuORF")), sum(G$class == CLASSES[1]), sum(G$class == CLASSES[2]), sum(G$class == CLASSES[3]),
             big(sum(G$cohort == "No TuORF")), big(sum(!is.na(G$U))), big(sum(!is.na(G$half_life)))))
D <- G[RNA_reads >= MIN_RNA & GMUCT_reads >= 1]
.add(sprintf("  genes with >= %d RNA-seq reads (uncapped-ratio set): %s single-TuORF (%d / %d / %d) + %s No TuORF; median RNA-seq reads %s (single-TuORF) / %s (No TuORF), median GMUCT reads %s / %s",
             MIN_RNA, big(sum(D$cohort == "single-TuORF")), sum(D$class == CLASSES[1]), sum(D$class == CLASSES[2]), sum(D$class == CLASSES[3]), big(sum(D$cohort == "No TuORF")),
             big(median(D[cohort == "single-TuORF", RNA_reads])), big(median(D[cohort == "No TuORF", RNA_reads])), big(median(D[cohort == "single-TuORF", GMUCT_reads])), big(median(D[cohort == "No TuORF", GMUCT_reads]))))
for (cl in CLASSES) .add(sprintf("    %s: %d of %d genes kept (%.0f%%); median RNA TPM (our data) %.1f", cl, sum(D$class == cl), sum(G$class == cl), 100 * mean(G$class[G$RNA_reads >= MIN_RNA & !is.na(G$RNA_reads)] == cl) * sum(G$RNA_reads >= MIN_RNA, na.rm = TRUE) / sum(G$class == cl), median(D[class == cl, RNA_TPM], na.rm = TRUE)))

################################################################################
# 3) Statistics helpers
################################################################################
sp_ci <- function(x, y, B = NBOOT) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < 10) return(c(rho = NA, lo = NA, hi = NA, n = n))
  r <- cor(x, y, method = "spearman"); rx <- rank(x); ry <- rank(y)
  b <- replicate(B, { i <- sample.int(n, n, TRUE); cor(rx[i], ry[i]) })      # ranks of the resampled values are not re-ranked (fast, close to exact for n > 200)
  c(rho = r, lo = unname(quantile(b, 0.025)), hi = unname(quantile(b, 0.975)), n = n)
}
psp <- function(x, y, z) {                                                   # partial Spearman correlation of x and y given z
  ok <- is.finite(x) & is.finite(y) & is.finite(z); rx <- rank(x[ok]); ry <- rank(y[ok]); rz <- rank(z[ok])
  rxy <- cor(rx, ry); rxz <- cor(rx, rz); ryz <- cor(ry, rz); (rxy - rxz * ryz) / sqrt((1 - rxz^2) * (1 - ryz^2))
}
fold_ci <- function(a, b, B = NBOOT) {                                       # median(a) / median(b) with a bootstrap interval
  a <- a[is.finite(a)]; b <- b[is.finite(b)]; f <- median(a) / median(b)
  bb <- replicate(B, median(sample(a, length(a), TRUE)) / median(sample(b, length(b), TRUE)))
  c(fold = f, lo = unname(quantile(bb, 0.025)), hi = unname(quantile(bb, 0.975)))
}
ks_letters <- function(d, grp, val) {                                         # pairwise KS tests, letters at p < 0.05, as for Figure 1C
  d <- d[is.finite(d[[val]]) & !is.na(d[[grp]]), ]; lv <- levels(droplevels(d[[grp]])); d[[grp]] <- factor(d[[grp]], lv)
  cmb <- combn(lv, 2, simplify = FALSE)
  pw <- rbindlist(lapply(cmb, function(g) { a <- d[[val]][d[[grp]] == g[1]]; b <- d[[val]][d[[grp]] == g[2]]
    data.table(group1 = g[1], group2 = g[2], n1 = length(a), n2 = length(b), KS_p = suppressWarnings(ks.test(a, b)$p.value),
               Wilcoxon_p = wilcox.test(a, b)$p.value, median_fold_1_over_2 = median(a) / median(b)) }))
  ids <- setNames(paste0("G", seq_along(lv)), lv); p <- pw$KS_p; names(p) <- paste(ids[pw$group1], ids[pw$group2], sep = "-")
  L <- multcompLetters(p, compare = "<", threshold = 0.05, Letters = letters)$Letters
  L <- L[ids]; names(L) <- lv
  ch <- unique(unlist(strsplit(unname(L), "")))                              # re-letter by order of appearance (a for the first level), keeping shared letters such as "ab"
  mp <- setNames(letters[seq_along(ch)], ch); L[] <- vapply(unname(L), function(s) paste(mp[strsplit(s, "")[[1]]], collapse = ""), "")
  n <- table(d[[grp]])[lv]
  list(pairs = pw, letters = L, labels = setNames(sprintf("%s (n = %s) %s", lv, trimws(big(as.integer(n))), L), lv), medians = tapply(d[[val]], d[[grp]], median)[lv])
}

################################################################################
# 4) Correlations (single-TuORF genes with an uncapped ratio)
################################################################################
S <- D[cohort == "single-TuORF"]
cors <- list()
addc <- function(label, x, y, set = S, extra = "") { v <- sp_ci(x, y); cors[[length(cors) + 1]] <<- data.table(comparison = label, note = extra, n = v[["n"]], rho = v[["rho"]], ci_lo = v[["lo"]], ci_hi = v[["hi"]]); v }
v1 <- addc("uncapped ratio vs uIF", S$U, S$uIF, extra = "GEO RNA-seq denominator; Figure S19A")
v2 <- addc("uncapped ratio vs mORF TE", S$U, S$TE, extra = "GEO RNA-seq denominator; Figure S19C")
v3 <- addc("uncapped ratio vs %uORF+", S$U, S$pUP, extra = "GEO RNA-seq denominator")
v4 <- addc("uncapped ratio vs RNA-seq RPM (Yu et al.)", S$U, S$RNA_rpm, extra = "expression dependence of the ratio")
v5 <- addc("uncapped ratio vs RNA TPM (this study)", S$U, S$RNA_TPM, extra = "expression dependence of the ratio")
v6 <- addc("uncapped ratio (our RNA-seq as denominator) vs uIF", S$U_labRNA, S$uIF, extra = "sensitivity; TE shares this denominator, so TE is not tested with it")
v7 <- addc("uncapped ratio vs mORF TE, No-TuORF genes", D[cohort == "No TuORF", U], D[cohort == "No TuORF", TE], extra = "does decapping track TE without a uORF")
v8 <- addc("uncapped ratio vs half-life, all genes with both", D$U, D$half_life, extra = "validation of the ratio as a decay readout")
v9 <- addc("half-life vs uIF", G[cohort == "single-TuORF", half_life], G[cohort == "single-TuORF", uIF], extra = "all single-TuORF genes with a half-life")
v10 <- addc("half-life vs mORF TE", G[cohort == "single-TuORF", half_life], G[cohort == "single-TuORF", TE], extra = "Figure S19F")
v11 <- addc("half-life vs %uORF+", G[cohort == "single-TuORF", half_life], G[cohort == "single-TuORF", pUP], extra = "")
v12 <- addc("half-life vs mORF TE, No-TuORF genes", G[cohort == "No TuORF", half_life], G[cohort == "No TuORF", TE], extra = "")
v13 <- addc("uIF vs mORF TE, same genes", S$uIF, S$TE, extra = "the Supplemental Fig. S2A relationship within the uncapped-ratio set")
cor_tab <- rbindlist(cors)
pr_log <- function(x, y) { ok <- is.finite(x) & is.finite(y) & x > 0 & y > 0; cor(log10(x[ok]), log10(y[ok])) }   # V41 (2026-11-12): Pearson r on log10 values without the genes that have a value of 0 (log10(0) is not finite; one cell of the table was NaN and showed #NUM! in Excel)
cor_tab[, pearson_r_log := NA_real_]
cor_tab[comparison == "uncapped ratio vs uIF", pearson_r_log := pr_log(S$U, S$uIF)]
cor_tab[comparison == "uncapped ratio vs mORF TE", pearson_r_log := pr_log(S$U, S$TE)]
cor_tab[comparison == "half-life vs mORF TE", pearson_r_log := G[cohort == "single-TuORF", pr_log(half_life, TE)]]
pr_n0 <- function(x, y) { ok <- is.finite(x) & is.finite(y); c(all = sum(ok & !(x > 0 & y > 0)), y0 = sum(ok & y <= 0)) }   # genes that the Pearson r leaves out (all) and those of them whose second value is 0 (y0); the note of the row says how many (column n is the n of rho)
for (cmp in list(list("uncapped ratio vs uIF", S$U, S$uIF, "uIF"), list("uncapped ratio vs mORF TE", S$U, S$TE, "mORF TE"), list("half-life vs mORF TE", G[cohort == "single-TuORF", half_life], G[cohort == "single-TuORF", TE], "mORF TE"))) { k <- pr_n0(cmp[[2]], cmp[[3]])
  if (k[["all"]] > 0) cor_tab[comparison == cmp[[1]], note := paste0(note, ifelse(nzchar(note), "; ", ""), sprintf("Pearson r without %d gene%s with a %s of 0", k[["all"]], ifelse(k[["all"]] == 1, "", "s"), ifelse(k[["y0"]] == k[["all"]], cmp[[4]], "value")))] }
# partial correlations given expression
part <- data.table(comparison = c("uncapped ratio vs uIF | RNA-seq RPM (Yu et al.)", "uncapped ratio vs uIF | RNA TPM (this study)",
                                  "uncapped ratio vs mORF TE | RNA-seq RPM (Yu et al.)", "uncapped ratio vs mORF TE | RNA TPM (this study)",
                                  "half-life vs mORF TE | RNA TPM (this study)", "half-life vs uIF | RNA TPM (this study)"),
                   partial_rho = c(psp(S$U, S$uIF, S$RNA_rpm), psp(S$U, S$uIF, S$RNA_TPM), psp(S$U, S$TE, S$RNA_rpm), psp(S$U, S$TE, S$RNA_TPM),
                                   G[cohort == "single-TuORF", psp(half_life, TE, RNA_TPM)], G[cohort == "single-TuORF", psp(half_life, uIF, RNA_TPM)]))
.add("  Spearman correlations (single-TuORF genes unless stated):")
for (i in seq_len(nrow(cor_tab))) .add(sprintf("    %-58s rho = %6.3f (95%% CI %.3f to %.3f), n = %s", cor_tab$comparison[i], cor_tab$rho[i], cor_tab$ci_lo[i], cor_tab$ci_hi[i], big(cor_tab$n[i])))
for (i in seq_len(nrow(part))) .add(sprintf("    partial %-50s rho = %6.3f", part$comparison[i], part$partial_rho[i]))
# cutoff and replicate sensitivity
sens <- rbindlist(lapply(CUTOFFS, function(k) { s <- G[cohort == "single-TuORF" & RNA_reads >= k & GMUCT_reads >= 1]
  rbind(data.table(RNA_reads_min = k, denominator = "both replicates", n = nrow(s), rho_uIF = cor(s$U, s$uIF, method = "spearman"), rho_TE = cor(s$U, s$TE, method = "spearman", use = "complete.obs"), rho_pUP = cor(s$U, s$pUP, method = "spearman")),
        { s1 <- s[is.finite(U_rep1)]; data.table(RNA_reads_min = k, denominator = "replicate 1 only", n = nrow(s1), rho_uIF = cor(s1$U_rep1, s1$uIF, method = "spearman"), rho_TE = cor(s1$U_rep1, s1$TE, method = "spearman", use = "complete.obs"), rho_pUP = cor(s1$U_rep1, s1$pUP, method = "spearman")) },
        { s2 <- s[is.finite(U_rep2)]; data.table(RNA_reads_min = k, denominator = "replicate 2 only", n = nrow(s2), rho_uIF = cor(s2$U_rep2, s2$uIF, method = "spearman"), rho_TE = cor(s2$U_rep2, s2$TE, method = "spearman", use = "complete.obs"), rho_pUP = cor(s2$U_rep2, s2$pUP, method = "spearman")) }) }))
.add("  sensitivity to the RNA-seq read cutoff and to the replicate:")
for (i in seq_len(nrow(sens))) .add(sprintf("    >= %2d reads, %-17s n = %5d: rho(U, uIF) = %6.3f, rho(U, TE) = %6.3f, rho(U, %%uORF+) = %6.3f", sens$RNA_reads_min[i], sens$denominator[i], sens$n[i], sens$rho_uIF[i], sens$rho_TE[i], sens$rho_pUP[i]))

################################################################################
# 5) Group comparisons (Figure 1C style)
################################################################################
grp_summary <- function(d, grp, val, lv_ref) {
  k <- ks_letters(as.data.frame(d), grp, val); lv <- names(k$letters)
  tab <- rbindlist(lapply(lv, function(g) { x <- d[[val]][d[[grp]] == g & is.finite(d[[val]])]
    fc <- if (g != lv_ref) fold_ci(x, d[[val]][d[[grp]] == lv_ref & is.finite(d[[val]])]) else c(fold = 1, lo = NA, hi = NA)
    data.table(group = g, n = length(x), median = median(x), q25 = unname(quantile(x, .25)), q75 = unname(quantile(x, .75)),
               fold_vs_ref = fc[["fold"]], fold_lo = fc[["lo"]], fold_hi = fc[["hi"]], KS_letter = k$letters[[g]]) }))
  list(table = tab, pairs = k$pairs, labels = k$labels)
}
gU  <- grp_summary(D, "class", "U", CLASSES[4])          # uncapped ratio by uIF class
gUp <- grp_summary(D, "pclass", "U", PUP_CLASSES[4])     # uncapped ratio by %uORF+ group
gH  <- grp_summary(G[is.finite(half_life)], "class", "half_life", CLASSES[4])    # half-life by uIF class (all genes with a half-life)
gHp <- grp_summary(G[is.finite(half_life)], "pclass", "half_life", PUP_CLASSES[4])
gT  <- grp_summary(D, "class", "TE", CLASSES[4])         # TE in the same (uncapped-ratio) set, to show the Figure 1C pattern holds there
gR  <- grp_summary(D, "class", "RNA_rpm", CLASSES[4])    # expression by class in the set
show <- function(nm, g, unit = "") { .add("  ", nm, ":"); for (i in seq_len(nrow(g$table))) .add(sprintf("    %-24s n = %5d, median %s%s, %.2fx vs No TuORF (%.2f to %.2f), letter %s", g$table$group[i], g$table$n[i], f2(g$table$median[i], 3), unit, g$table$fold_vs_ref[i], g$table$fold_lo[i], g$table$fold_hi[i], g$table$KS_letter[i]))
  for (i in seq_len(nrow(g$pairs))) .add(sprintf("      %s vs %s: KS p = %.2g, Wilcoxon p = %.2g, median fold %.2f", g$pairs$group1[i], g$pairs$group2[i], g$pairs$KS_p[i], g$pairs$Wilcoxon_p[i], g$pairs$median_fold_1_over_2[i])) }
show("uncapped ratio by uIF class", gU); show("uncapped ratio by %uORF+ group", gUp); show("half-life (h) by uIF class", gH, " h"); show("half-life (h) by %uORF+ group", gHp, " h")
show("mORF TE by uIF class, uncapped-ratio set", gT); show("RNA-seq RPM (Yu et al.) by uIF class, uncapped-ratio set", gR)
# expression-stratified check: median U by class within RNA-seq RPM tertiles of the whole set
D[, rna_tertile := cut(RNA_rpm, quantile(RNA_rpm, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("low RNA", "mid RNA", "high RNA"))]
strat <- D[, .(n = .N, median_U = median(U)), by = .(rna_tertile, class)][order(rna_tertile, class)]
strat_w <- dcast(strat, class ~ rna_tertile, value.var = "median_U")
.add("  median uncapped ratio by class within RNA-seq abundance tertiles:"); for (i in seq_len(nrow(strat_w))) .add(sprintf("    %-24s low %.3f, mid %.3f, high %.3f", strat_w$class[i], strat_w$`low RNA`[i], strat_w$`mid RNA`[i], strat_w$`high RNA`[i]))

################################################################################
# 6) Figure S19
################################################################################
th <- theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), legend.position = "bottom", legend.direction = "vertical", legend.key.size = unit(0.3, "cm"),
                                      legend.text = element_text(size = 7), plot.title = element_text(size = 9.5, face = "bold"), plot.subtitle = element_text(size = 7.3), legend.margin = margin(0, 0, 0, 0))
pal4 <- c("#F8766D", "#00BA38", "#619CFF", "grey40")      # V38 (2026-11-11): the colours of Supplemental Fig. S2B (High, Medium, Low uIF; No TuORF)
pal_pup <- c("#619CFF", "#00BA38", "#F8766D", "grey40")   # V38: the colours of Figure 2F (0-30%, 30-80%, 80-100% uORF+; No TuORF). Before, the 0-30% group had the red of the 80-100% group of Figure 2F.
lab_plain <- function(x) format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)
scat <- function(d, x, y, xl, yl, title, sub, v) {
  d <- d[is.finite(d[[x]]) & is.finite(d[[y]]) & d[[x]] > 0 & d[[y]] > 0, ]                 # log axes: zeros are kept in the statistics, not drawn
  ggplot(d, aes(x = .data[[x]], y = .data[[y]])) + geom_point(alpha = 0.25, size = 0.7, colour = "grey30") +
    geom_smooth(method = "lm", se = FALSE, colour = "#D7301F", linewidth = 0.7, formula = y ~ x) + scale_x_log10(labels = lab_plain) + scale_y_log10(labels = lab_plain) +
    labs(title = title, subtitle = sprintf("%s; Spearman rho = %.2f (%.2f to %.2f), n = %s", sub, v[["rho"]], v[["lo"]], v[["hi"]], trimws(big(v[["n"]]))), x = xl, y = yl) + th   # statistics in the subtitle, never over the points
}
ecdfp <- function(d, grp, val, g, xl, title, sub, pal = pal4, xlim = NULL) {
  d <- d[is.finite(d[[val]]) & !is.na(d[[grp]]) & d[[val]] > 0, ]; d$lab <- factor(g$labels[as.character(d[[grp]])], g$labels)
  p <- ggplot(d, aes(x = .data[[val]], colour = lab)) + stat_ecdf(geom = "step", linewidth = 0.8) + scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
    labs(title = title, subtitle = sub, x = xl, y = "Fraction of genes") + th + guides(colour = guide_legend(ncol = 1)) + scale_x_log10(labels = lab_plain)
  if (!is.null(xlim)) p <- p + coord_cartesian(xlim = xlim)                                # ECDFs use all genes; the view is clipped to the central range
  p
}
U_lab <- "Uncapped ratio (GMUCT / RNA-seq, reads per million)"
pA <- scat(S, "uIF", "U", "uIF", "Uncapped ratio", "A", sprintf("Single-TuORF genes, >= %d RNA-seq reads (Yu et al. 2021, Col-0, no ABA); %d with uIF = 0 not drawn", MIN_RNA, sum(S$uIF == 0)), v1)
pB <- ecdfp(D, "class", "U", gU, U_lab, "B", "Letters: pairwise Kolmogorov–Smirnov tests, p < 0.05 (as Fig. S2B); x axis clipped to 0.1 to 50", xlim = c(0.1, 50))
pC <- scat(S, "TE", "U", "mORF TE", "Uncapped ratio", "C", "Same genes as A", v2)
pD <- ecdfp(D, "pclass", "U", gUp, U_lab, "D", "By %uORF+ group (single-TuORF genes) and No-TuORF genes; x axis clipped to 0.1 to 50", pal = pal_pup, xlim = c(0.1, 50))
pE <- ecdfp(G[is.finite(half_life)], "class", "half_life", gH, "mRNA half-life (h; 5-EU labeling)", "E", "All single-TuORF and No-TuORF genes with a half-life")
pF <- scat(G[cohort == "single-TuORF" & is.finite(half_life)], "TE", "half_life", "mORF TE", "mRNA half-life (h)", "F", "Single-TuORF genes with a half-life", v10)
source("revision_modules/fig_compact.R")   # 6.5 x <= 9 in, Arial, >= 8 pt (2026-10-28)
g <- arrangeGrob(patchwork::patchworkGrob(patchwork::wrap_plots(lapply(list(pA, pB, pC, pD, pE, pF), function(p) compact(p) + ggplot2::theme(panel.heights = grid::unit(1.25, "in"), panel.widths = grid::unit(2.0, "in"))), ncol = 2)), top = fig_title("Figure S19. Uncapped mRNA (GMUCT) and mRNA half-life versus uIF and mORF TE"))   # V21 (2026-11-07, Larry): six plotting areas of equal size, panels aligned
fig_save(file.path(OUT_DIR, "A11_Figure_uncapped_vs_uIF.pdf"), g, h = 9)

################################################################################
# 7) Tables and summary
################################################################################
genes_out <- G[order(cohort != "single-TuORF", class, gene_id), .(gene_id, cohort, uIF_class = class, pct_uORF_plus_group = pclass, uIF, pct_uORF_plus = pUP, mORF_TE = TE, RNA_TPM_this_study = RNA_TPM,
  GMUCT_reads_rep1 = GMUCT_rep1, GMUCT_reads_rep2 = GMUCT_rep2, RNA_reads_rep1 = RNA_rep1, RNA_reads_rep2 = RNA_rep2, GMUCT_RPM = GMUCT_rpm, RNA_RPM = RNA_rpm,
  uncapped_ratio = U, uncapped_ratio_rep1 = U_rep1, uncapped_ratio_rep2 = U_rep2, uncapped_ratio_ourRNA = U_labRNA, in_uncapped_ratio_set = RNA_reads >= MIN_RNA & GMUCT_reads >= 1, half_life_h = half_life)]
fwrite(genes_out, file.path(OUT_DIR, "A11_Table_genes.csv"))
wb <- createWorkbook()
addWorksheet(wb, "README"); writeData(wb, "README", data.frame(x = c(
  "Table S10. Uncapped mRNA (GMUCT) and mRNA half-life versus uIF and mORF TE (Figure S19).",
  "GMUCT and total RNA-seq of Col-0 12-day-old seedlings without ABA from Yu et al. 2021 Dev Cell 56:125-140 (GEO GSE142388; GMUCT GSM4226794/GSM4226795, total RNA-seq GSM4226810/GSM4226811), deposited HTSeq gene counts (STAR, TAIR10).",
  "Reads per million (RPM) use the reads assigned to nuclear genes of each library; the uncapped ratio of a gene is the mean GMUCT RPM of the two replicates divided by the mean RNA-seq RPM (= 2 to the power of the 'proportion uncapped' of Yu et al. 2021, who use log2[RPM_GMUCT/RPM_RNA-seq]).",
  sprintf("Genes enter the uncapped-ratio analyses with at least %d RNA-seq reads (both replicates summed) and at least one GMUCT read; cutoffs of 5 and 20 reads and single replicates are in 'cutoff_sensitivity' (rho_uIF, rho_TE and rho_pUP = Spearman rho of the uncapped ratio with uIF, mORF TE and %%uORF+). In the sheet 'genes', uncapped_ratio is NA (empty cell) for genes without RNA-seq reads and for the %d genes that are not in the deposited count tables; uncapped_ratio_rep1 and uncapped_ratio_rep2 are NA when that replicate has no RNA-seq read or no GMUCT read; uncapped_ratio_ourRNA = GMUCT RPM divided by the RNA TPM of this study.", MIN_RNA, sum(is.na(G$RNA_reads))),
  "Half-lives: 5-EU pulse-chase labelling data of Szabo et al. 2020 (mean of three replicates, hours).",
  "Cohorts: single-TuORF genes of Table S3 with their uIF classes (Supplemental Fig. S2B) and %uORF+ groups; No-TuORF genes as in Supplemental Fig. S2B (annotated 5' leader, RNA TPM > 1, TE defined, no TuORF and no overlapping uORF).",
  "Sheets: correlations (Spearman rho with bootstrap 95% intervals over genes, 1,000 resamples; the column pearson_r_log gives the Pearson r of the log10 values for some rows, and genes with a value of 0 are not used for it), partial_correlations (rank-based partial correlations given RNA abundance), cutoff_sensitivity, group_summaries (median, quartiles, median fold versus No TuORF with bootstrap 95% interval (1,000 resamples; fold_lo and fold_hi), letters from pairwise Kolmogorov–Smirnov (KS) tests at p < 0.05 as in Supplemental Fig. S2B), pairwise_tests, expression_strata (median uncapped ratio by class within RNA-seq abundance tertiles), genes (per-gene values).",
  "Code: revision_modules/A11_uncapped_vs_uIF.R in the Supplemental Code.")), colNames = FALSE)
addWorksheet(wb, "correlations"); writeData(wb, "correlations", cor_tab)
addWorksheet(wb, "partial_correlations"); writeData(wb, "partial_correlations", part)
addWorksheet(wb, "cutoff_sensitivity"); writeData(wb, "cutoff_sensitivity", sens)
gs <- rbind(cbind(measure = "uncapped ratio", grouping = "uIF class", gU$table), cbind(measure = "uncapped ratio", grouping = "%uORF+ group", gUp$table),
            cbind(measure = "half-life (h)", grouping = "uIF class", gH$table), cbind(measure = "half-life (h)", grouping = "%uORF+ group", gHp$table),
            cbind(measure = "mORF TE (uncapped-ratio set)", grouping = "uIF class", gT$table), cbind(measure = "RNA-seq RPM, Yu et al. (uncapped-ratio set)", grouping = "uIF class", gR$table))
addWorksheet(wb, "group_summaries"); writeData(wb, "group_summaries", gs)
pt <- rbind(cbind(measure = "uncapped ratio", grouping = "uIF class", gU$pairs), cbind(measure = "uncapped ratio", grouping = "%uORF+ group", gUp$pairs),
            cbind(measure = "half-life (h)", grouping = "uIF class", gH$pairs), cbind(measure = "half-life (h)", grouping = "%uORF+ group", gHp$pairs),
            cbind(measure = "mORF TE (uncapped-ratio set)", grouping = "uIF class", gT$pairs), cbind(measure = "RNA-seq RPM, Yu et al. (uncapped-ratio set)", grouping = "uIF class", gR$pairs))
addWorksheet(wb, "pairwise_tests"); writeData(wb, "pairwise_tests", pt)
addWorksheet(wb, "expression_strata"); writeData(wb, "expression_strata", strat)
genes_x <- copy(genes_out); genes_x[!is.finite(uncapped_ratio), uncapped_ratio := NA_real_]   # V41 (2026-11-12): the ratio is not defined without RNA-seq reads (Inf or NaN, shown as #NUM! in Excel); the sheet gets an empty cell there, the csv file is unchanged
addWorksheet(wb, "genes"); writeData(wb, "genes", genes_x)
saveWorkbook(wb, file.path(OUT_DIR, "A11_Table_uncapped_halflife.xlsx"), overwrite = TRUE)
.add(sprintf("A11 done in %.1f min; outputs: output/A11_Figure_uncapped_vs_uIF.pdf, A11_Table_uncapped_halflife.xlsx, A11_Table_genes.csv, A11_summary.txt", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
writeLines(.log, file.path(OUT_DIR, "A11_summary.txt"))
