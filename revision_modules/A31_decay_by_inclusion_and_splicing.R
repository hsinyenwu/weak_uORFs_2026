################################################################################
# A31 — mRNA half-life and uncapped mRNA by uORF inclusion and by leader splicing #
# (Larry, 2026-11-11: "For the degradome data, should we include %uORF+,          #
#  splicing vs degradome and also make an independent section at the end?"        #
#  -> "Make a separate section")                                                  #
#                                                                                #
# Question: A11 and A12 compare the decay measures with uIF and with %uORF+       #
# (uORF inclusion by transcription start site). Here the same measures are        #
# compared with leader splicing at the uORF start codon:                          #
#   1) the %uORF+ groups again, in the genes unaffected by leader splicing;       #
#   2) among genes with 80-100% uORF+, genes whose start codon is retained        #
#      (S >= 0.90) against genes in which leader splicing removes the start       #
#      codon from at least half of the transcripts (%uORF+/100 x (1 - S) >= 0.5,  #
#      the rule of Figure 2G);                                                    #
#   3) Spearman correlations with %uORF+, with S and with the effective           #
#      inclusion (%uORF+ x S);                                                    #
#   4) (after the independent read of 2026-11-11) the %uORF+ groups compared with #
#      genes without TuORFs after adjustment for transcript length (exon length   #
#      of the gene) and RNA abundance, and within thirds of RNA abundance. The    #
#      uncapped ratio divides 5'-end reads by RNA-seq reads, and RNA-seq reads    #
#      increase with transcript length; the groups differ in length.              #
# Measures: mRNA half-life (5-EU pulse-chase, Szabo et al. 2020), uncapped ratio  #
# from PARE (Nagarajan et al. 2019; A12) and from GMUCT (Yu et al. 2021; A11),    #
# and mORF TE for reference.                                                      #
#                                                                                #
# Inputs (per-gene tables written by other modules; run A2, A11 and A12 first):   #
#   output/A11_Table_genes.csv, output/A12_Table_genes.csv,                       #
#   output/A2_single_TuORF_cohort_splicing_flags.csv,                             #
#   output/A11_Table_uncapped_halflife.xlsx (to check the half-life groups),      #
#   output/A12_gene_counts.csv (exon length of each gene).                        #
# Outputs (output/): A31_Figure_decay_by_inclusion_and_splicing.pdf (Figure S21), #
#   A31_Table_decay_inclusion_splicing.xlsx (Table S12), A31_Table_genes.csv,     #
#   A31_summary.txt, A31_csv_*.csv (the sheets of the table as CSV files)         #
# Run from Weak_uORF_submission/:                                                 #
#   Rscript revision_modules/A31_decay_by_inclusion_and_splicing.R                #
################################################################################
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(openxlsx); library(gridExtra); library(grid); library(multcompView)
})
options(stringsAsFactors = FALSE, warn = 1)
set.seed(31)
t0 <- Sys.time()
OUT_DIR <- "output"
P11 <- file.path(OUT_DIR, "A11_Table_genes.csv"); P12 <- file.path(OUT_DIR, "A12_Table_genes.csv")
PFL <- file.path(OUT_DIR, "A2_single_TuORF_cohort_splicing_flags.csv"); P11X <- file.path(OUT_DIR, "A11_Table_uncapped_halflife.xlsx"); PLEN <- file.path(OUT_DIR, "A12_gene_counts.csv")
for (f in c(P11, P12, PFL, P11X, PLEN)) stopifnot(file.exists(f))
NBOOT <- 1000L
MIN_COV <- 10L        # RNA-seq reads across the start codon required for S (as in A2)
S_RETAINED <- 0.90    # start codon retained (as in A2: "retained (S >= 0.90)")
LOSS_SPLICED <- 0.5   # start codon removed from at least half of the transcripts: %uORF+/100 x (1 - S) >= 0.5 (as in Figure 2G)
PUP <- c("0 to 30% uORF+", "30 to 80% uORF+", "80 to 100% uORF+", "No TuORF")
PUP_KEY <- c("0–30% uORF+", "30–80% uORF+", "80–100% uORF+", "No TuORF")
SC <- c("Start codon retained", "Removed by splicing", "No TuORF")
.log <- character(0); .add <- function(...) { s <- paste0(...); message(s); .log <<- c(.log, s) }
big <- function(x) trimws(format(x, big.mark = ","))

################################################################################
# 1) Per-gene table
################################################################################
a11 <- fread(P11); a12 <- fread(P12); fl <- fread(PFL)
stopifnot(identical(sort(a11$gene_id), sort(a12$gene_id)), !anyDuplicated(a11$gene_id), !anyDuplicated(fl$gene_id))
a12 <- a12[match(a11$gene_id, gene_id)]
stopifnot(identical(a11$cohort, a12$cohort), isTRUE(all.equal(a11$half_life_h, a12$half_life_h)), isTRUE(all.equal(a11$pct_uORF_plus, a12$pct_uORF_plus)))
G <- data.table(gene_id = a11$gene_id, cohort = a11$cohort, uIF_class = a11$uIF_class, pclass = factor(a11$pct_uORF_plus_group, PUP), uIF = a11$uIF, pUP = a11$pct_uORF_plus,
                TE = a11$mORF_TE, RNA_TPM = a11$RNA_TPM_this_study, HL = a11$half_life_h,
                U_gmuct = ifelse(a11$in_uncapped_ratio_set %in% TRUE, a11$uncapped_ratio, NA_real_), U_pare = ifelse(a12$in_uncapped_ratio_set %in% TRUE, a12$uncapped_ratio, NA_real_),
                U_gmuct_tpm = ifelse(a11$in_uncapped_ratio_set %in% TRUE, a11$uncapped_ratio_ourRNA, NA_real_), U_pare_tpm = ifelse(a12$in_uncapped_ratio_set %in% TRUE, a12$uncapped_ratio_ourRNA, NA_real_))
G[!is.finite(TE), TE := NA_real_]                                              # genes with TE = 0 are kept for the groups and the plot (linear axis), as in Figure 2F; they are left out of the models on log values
G[!is.finite(U_gmuct_tpm), U_gmuct_tpm := NA_real_]; G[!is.finite(U_pare_tpm), U_pare_tpm := NA_real_]
len <- fread(PLEN, select = c("gene_id", "exon_length")); G[, exon_length := len$exon_length[match(gene_id, len$gene_id)]]; stopifnot(!anyNA(G$exon_length))
single <- G$cohort == "single-TuORF"
stopifnot(sum(single) == nrow(fl), setequal(G$gene_id[single], fl$gene_id), !anyNA(G$pclass))
G <- merge(G, fl[, .(gene_id, annot_class, junction_status, splicing_status = proxy_status, S = S_TIS, cov = cov_TIS, eff = pUP_eff_CAGE, pUP_flags = pct_uORF_plus)], by = "gene_id", all.x = TRUE)
stopifnot(G[cohort == "single-TuORF", max(abs(pUP - pUP_flags))] < 1e-6); G[, pUP_flags := NULL]
G[, clean := grepl("^clean", splicing_status)]; G[, affected := grepl("^affected", splicing_status)]
G[, S_ok := cohort == "single-TuORF" & is.finite(S) & !is.na(cov) & cov >= MIN_COV]
G[, loss_spliced := ifelse(S_ok, pUP / 100 * (1 - S), NA_real_)]                 # fraction of the transcripts that lose the start codon by leader splicing
G[, eff_incl := ifelse(S_ok, pUP * S, NA_real_)]                                 # effective inclusion (%), genes with S only
stopifnot(G[S_ok == TRUE, max(abs(eff_incl - eff))] < 0.01)                      # the same value as in the A2 table
G[, sc := fcase(cohort != "single-TuORF", SC[3],
                S_ok & pclass == PUP[3] & S >= S_RETAINED, SC[1],
                S_ok & pclass == PUP[3] & loss_spliced >= LOSS_SPLICED, SC[2], default = NA_character_)]
G[, sc := factor(sc, SC)]
.add("A31: ", big(sum(single)), " single-TuORF genes and ", big(sum(!single)), " No-TuORF genes (tables of A11, A12 and A2)")
.add(sprintf("  leader-splicing status of the single-TuORF genes: %s unaffected (clean), %s affected, %s not assessable; S available (>= %d RNA-seq reads across the start codon) for %s",
             big(sum(G$clean)), big(sum(G$affected)), big(sum(single) - sum(G$clean) - sum(G$affected)), MIN_COV, big(sum(G$S_ok))))
h <- G[cohort == "single-TuORF" & pclass == PUP[3]]
.add(sprintf("  genes with 80-100%% uORF+: %s; with S: %s; start codon retained (S >= %.2f): %s; removed from at least half of the transcripts by leader splicing: %s; in between: %s",
             big(nrow(h)), big(sum(h$S_ok)), S_RETAINED, big(sum(h$sc == SC[1], na.rm = TRUE)), big(sum(h$sc == SC[2], na.rm = TRUE)), big(sum(h$S_ok & is.na(h$sc)))))
.add(sprintf("    removed group: median S %.2f (range %.2f to %.2f), median %%uORF+ %.1f; splicing status: %s", h[sc == SC[2], median(S)], h[sc == SC[2], min(S)], h[sc == SC[2], max(S)], h[sc == SC[2], median(pUP)],
             paste(names(table(h[sc == SC[2], splicing_status])), table(h[sc == SC[2], splicing_status]), collapse = "; ")))
.add(sprintf("    removed group, annotation classes: %s", paste(substr(names(table(h[sc == SC[2], annot_class])), 1, 2), table(h[sc == SC[2], annot_class]), collapse = "; ")))

################################################################################
# 2) Statistics helpers (as in A11)
################################################################################
sp_ci <- function(x, y, B = NBOOT) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < 10) return(c(rho = NA, lo = NA, hi = NA, n = n, p = NA))
  r <- cor(x, y, method = "spearman"); rx <- rank(x); ry <- rank(y)
  b <- replicate(B, { i <- sample.int(n, n, TRUE); cor(rx[i], ry[i]) })
  c(rho = r, lo = unname(quantile(b, 0.025)), hi = unname(quantile(b, 0.975)), n = n, p = suppressWarnings(cor.test(x, y, method = "spearman")$p.value))
}
fold_ci <- function(a, b, B = NBOOT) {                                           # median(a) / median(b) with a bootstrap interval
  a <- a[is.finite(a)]; b <- b[is.finite(b)]; f <- median(a) / median(b)
  bb <- replicate(B, median(sample(a, length(a), TRUE)) / median(sample(b, length(b), TRUE)))
  c(fold = f, lo = unname(quantile(bb, 0.025)), hi = unname(quantile(bb, 0.975)))
}
ks_letters <- function(d, grp, val, key = NULL) {                                # pairwise KS tests, letters at p < 0.05
  d <- as.data.frame(d); d <- d[is.finite(d[[val]]) & !is.na(d[[grp]]), ]; lv <- levels(droplevels(d[[grp]])); d[[grp]] <- factor(d[[grp]], lv)
  pw <- rbindlist(lapply(combn(lv, 2, simplify = FALSE), function(g) { a <- d[[val]][d[[grp]] == g[1]]; b <- d[[val]][d[[grp]] == g[2]]
    data.table(group1 = g[1], group2 = g[2], n1 = length(a), n2 = length(b), KS_p = suppressWarnings(ks.test(a, b)$p.value),
               Wilcoxon_p = wilcox.test(a, b)$p.value, median_fold_1_over_2 = median(a) / median(b)) }))
  ids <- setNames(paste0("G", seq_along(lv)), lv); p <- pw$KS_p; names(p) <- paste(ids[pw$group1], ids[pw$group2], sep = "-")
  L <- multcompLetters(p, compare = "<", threshold = 0.05, Letters = letters)$Letters
  L <- L[ids]; names(L) <- lv
  ch <- unique(unlist(strsplit(unname(L), ""))); mp <- setNames(letters[seq_along(ch)], ch)
  L[] <- vapply(unname(L), function(s) paste(mp[strsplit(s, "")[[1]]], collapse = ""), "")
  n <- table(d[[grp]])[lv]; shown <- if (is.null(key)) lv else key[lv]
  list(pairs = pw, letters = L, labels = setNames(sprintf("%s (n = %s) %s", shown, trimws(big(as.integer(n))), L), lv))
}
grp_summary <- function(d, grp, val, lv_ref, key = NULL) {
  k <- ks_letters(d, grp, val, key); lv <- names(k$letters)
  tab <- rbindlist(lapply(lv, function(g) { x <- d[[val]][d[[grp]] == g & is.finite(d[[val]]) & !is.na(d[[grp]])]
    fc <- if (g != lv_ref) fold_ci(x, d[[val]][d[[grp]] == lv_ref & is.finite(d[[val]]) & !is.na(d[[grp]])]) else c(fold = 1, lo = NA, hi = NA)
    data.table(group = g, n = length(x), median = median(x), q25 = unname(quantile(x, .25)), q75 = unname(quantile(x, .75)),
               fold_vs_No_TuORF = fc[["fold"]], fold_lo = fc[["lo"]], fold_hi = fc[["hi"]], KS_letter = k$letters[[g]]) }))
  list(table = tab, pairs = k$pairs, labels = k$labels)
}
MEAS <- c(HL = "mRNA half-life (h)", U_pare = "uncapped ratio (PARE)", U_gmuct = "uncapped ratio (GMUCT)", TE = "mORF TE")
show <- function(nm, g) { .add("  ", nm, ":")
  for (i in seq_len(nrow(g$table))) .add(sprintf("    %-22s n = %5d, median %.3f, %.2fx vs No TuORF (%.2f to %.2f), letter %s", g$table$group[i], g$table$n[i], g$table$median[i], g$table$fold_vs_No_TuORF[i], g$table$fold_lo[i], g$table$fold_hi[i], g$table$KS_letter[i]))
  for (i in seq_len(nrow(g$pairs))) .add(sprintf("      %s vs %s: KS p = %.2g, Wilcoxon p = %.2g, median fold %.2f", g$pairs$group1[i], g$pairs$group2[i], g$pairs$KS_p[i], g$pairs$Wilcoxon_p[i], g$pairs$median_fold_1_over_2[i])) }

################################################################################
# 3) %uORF+ groups: all single-TuORF genes (as A11 / A12) and the genes unaffected by leader splicing
################################################################################
Gall   <- G                                                                     # single-TuORF genes + No-TuORF genes
Gclean <- G[cohort != "single-TuORF" | clean]
gs_all <- lapply(c("HL", "U_pare", "U_gmuct"), function(v) grp_summary(Gall, "pclass", v, PUP[4], setNames(PUP_KEY, PUP)));   names(gs_all) <- c("HL", "U_pare", "U_gmuct")
gs_cln <- lapply(c("HL", "U_pare", "U_gmuct"), function(v) grp_summary(Gclean, "pclass", v, PUP[4], setNames(PUP_KEY, PUP))); names(gs_cln) <- c("HL", "U_pare", "U_gmuct")
# the half-life groups of all single-TuORF genes are those of A11 (Table S10): same n, medians and letters
x11 <- as.data.table(read.xlsx(P11X, sheet = "group_summaries")); x11 <- x11[measure == "half-life (h)" & grouping == "%uORF+ group"]
stopifnot(identical(x11$group, gs_all$HL$table$group), identical(as.integer(x11$n), gs_all$HL$table$n), isTRUE(all.equal(x11$median, gs_all$HL$table$median)), identical(x11$KS_letter, gs_all$HL$table$KS_letter))
for (v in names(gs_all)) show(paste0(MEAS[[v]], " by %uORF+ group, all single-TuORF genes"), gs_all[[v]])
for (v in names(gs_cln)) show(paste0(MEAS[[v]], " by %uORF+ group, genes unaffected by leader splicing (n = ", big(sum(G$clean)), ")"), gs_cln[[v]])

################################################################################
# 4) Start codon retained or removed by leader splicing, genes with 80-100% uORF+
################################################################################
Gsc <- G[!is.na(sc)]
gs_sc <- lapply(names(MEAS), function(v) grp_summary(Gsc, "sc", v, SC[3])); names(gs_sc) <- names(MEAS)
for (v in names(gs_sc)) show(paste0(MEAS[[v]], ": start codon retained or removed by leader splicing (genes with 80-100% uORF+)"), gs_sc[[v]])
two <- Gsc[sc %in% SC[1:2]]; two[, removed := as.integer(sc == SC[2])]
fold_model <- function(dd, covars) {                                             # fold of the geometric mean, removed over retained, from a linear model of the log10 value
  m <- lm(as.formula(paste("y ~ removed", if (length(covars)) paste("+", paste(covars, collapse = " + ")) else "")), data = dd); s <- summary(m)$coefficients; ci <- confint(m)["removed", ]
  c(fold = 10^s["removed", "Estimate"], lo = 10^ci[[1]], hi = 10^ci[[2]], p = s["removed", "Pr(>|t|)"]) }
direct <- rbindlist(lapply(names(MEAS), function(v) { d <- two[is.finite(get(v))]
  fc <- fold_ci(d[removed == 1][[v]], d[removed == 0][[v]])
  dd <- data.frame(y = log10(d[[v]]), removed = d$removed, lrna = log10(d$RNA_TPM), llen = log10(d$exon_length)); dd <- dd[is.finite(dd$y) & is.finite(dd$lrna) & is.finite(dd$llen), ]   # log values: genes with a value of 0 are left out of the models
  m0 <- fold_model(dd, character(0)); mA <- fold_model(dd, "lrna"); mL <- fold_model(dd, "llen"); mLA <- fold_model(dd, c("llen", "lrna"))
  data.table(measure = MEAS[[v]], n_retained = sum(d$removed == 0), n_removed = sum(d$removed == 1), median_retained = median(d[removed == 0][[v]]), median_removed = median(d[removed == 1][[v]]),
             median_fold_removed_over_retained = fc[["fold"]], fold_lo = fc[["lo"]], fold_hi = fc[["hi"]], Wilcoxon_p = wilcox.test(d[removed == 1][[v]], d[removed == 0][[v]])$p.value,
             KS_p = suppressWarnings(ks.test(d[removed == 1][[v]], d[removed == 0][[v]])$p.value), n_in_models = nrow(dd),
             geometric_mean_fold_unadjusted = m0[["fold"]], unadjusted_lo = m0[["lo"]], unadjusted_hi = m0[["hi"]], unadjusted_p = m0[["p"]],
             fold_adjusted_for_RNA_abundance = mA[["fold"]], adjA_lo = mA[["lo"]], adjA_hi = mA[["hi"]], adjA_p = mA[["p"]],
             fold_adjusted_for_transcript_length = mL[["fold"]], adjL_lo = mL[["lo"]], adjL_hi = mL[["hi"]], adjL_p = mL[["p"]],
             fold_adjusted_for_length_and_abundance = mLA[["fold"]], adjLA_lo = mLA[["lo"]], adjLA_hi = mLA[["hi"]], adjLA_p = mLA[["p"]]) }))
direct[, KS_p_Holm_three_decay_measures := NA_real_]; direct[measure != MEAS[["TE"]], KS_p_Holm_three_decay_measures := p.adjust(KS_p, "holm")]
.add("  removed against retained (genes with 80-100% uORF+): median fold with bootstrap interval; folds of the geometric mean from linear models of the log10 value (unadjusted; adjusted for log10 RNA TPM; for log10 exon length; for both):")
for (i in seq_len(nrow(direct))) .add(sprintf("    %-24s n = %4d / %3d, medians %.3f / %.3f, %.2fx (%.2f to %.2f), Wilcoxon p = %.2g, KS p = %.2g; adjusted %.2fx (%.2f to %.2f), p = %.2g; unadjusted model %.2fx, p = %.2g; adjusted for length %.2fx, p = %.2g; for length and abundance %.2fx (%.2f to %.2f), p = %.2g; Holm KS p (three decay measures) %s",
  direct$measure[i], direct$n_retained[i], direct$n_removed[i], direct$median_retained[i], direct$median_removed[i], direct$median_fold_removed_over_retained[i], direct$fold_lo[i], direct$fold_hi[i],
  direct$Wilcoxon_p[i], direct$KS_p[i], direct$fold_adjusted_for_RNA_abundance[i], direct$adjA_lo[i], direct$adjA_hi[i], direct$adjA_p[i], direct$geometric_mean_fold_unadjusted[i], direct$unadjusted_p[i],
  direct$fold_adjusted_for_transcript_length[i], direct$adjL_p[i], direct$fold_adjusted_for_length_and_abundance[i], direct$adjLA_lo[i], direct$adjLA_hi[i], direct$adjLA_p[i],
  if (is.na(direct$KS_p_Holm_three_decay_measures[i])) "not applicable" else sprintf("%.2g", direct$KS_p_Holm_three_decay_measures[i])))
.add(sprintf("    of the %s genes with a retained start codon, %s are unaffected by leader splicing as defined in A2 (Supplemental Fig. S7); the others are affected at another site (classes: %s)",
             big(sum(two$removed == 0)), big(sum(two$removed == 0 & two$clean)), paste(substr(names(table(two[removed == 0 & clean == FALSE, annot_class])), 1, 2), table(two[removed == 0 & clean == FALSE, annot_class]), collapse = "; ")))
desc <- two[, .(n = .N, median_S = median(S), median_pct_uORF_plus = median(pUP), median_effective_inclusion = median(eff_incl), median_uIF = median(uIF), median_RNA_TPM = median(RNA_TPM, na.rm = TRUE), median_exon_length = median(exon_length)), by = .(group = as.character(sc))]
for (i in seq_len(nrow(desc))) .add(sprintf("    %-22s n = %4d: median S %.2f, median %%uORF+ %.1f, median effective inclusion %.1f%%, median uIF %.3f, median RNA TPM %.1f, median exon length %.0f nt", desc$group[i], desc$n[i], desc$median_S[i], desc$median_pct_uORF_plus[i], desc$median_effective_inclusion[i], desc$median_uIF[i], desc$median_RNA_TPM[i], desc$median_exon_length[i]))

################################################################################
# 5) Correlations with %uORF+, with S and with the effective inclusion
################################################################################
cors <- list()
addc <- function(set, d, v, x, xn) { r <- sp_ci(d[[v]], d[[x]]); cors[[length(cors) + 1]] <<- data.table(gene_set = set, measure = MEAS[[v]], versus = xn, n = r[["n"]], rho = r[["rho"]], ci_lo = r[["lo"]], ci_hi = r[["hi"]], p = r[["p"]]) }
S1 <- G[cohort == "single-TuORF"]
for (v in c("HL", "U_pare", "U_gmuct")) {
  addc("all single-TuORF genes", S1, v, "pUP", "%uORF+")
  addc("unaffected by leader splicing", S1[clean == TRUE], v, "pUP", "%uORF+")
  addc("affected by leader splicing", S1[affected == TRUE], v, "pUP", "%uORF+")
  addc(sprintf("at least %d RNA-seq reads across the start codon", MIN_COV), S1[S_ok == TRUE], v, "pUP", "%uORF+")
  addc(sprintf("at least %d RNA-seq reads across the start codon", MIN_COV), S1[S_ok == TRUE], v, "eff_incl", "effective inclusion (%uORF+ x S)")
  addc(sprintf("at least %d RNA-seq reads across the start codon", MIN_COV), S1[S_ok == TRUE], v, "S", "S (fraction of reads that retain the start codon)") }
cor_tab <- rbindlist(cors)
psp <- function(x, y, Z) {                                                       # rank-based partial correlation of x and y given the covariates in the list Z
  ok <- is.finite(x) & is.finite(y) & Reduce(`&`, lapply(Z, is.finite)); rx <- rank(x[ok]); ry <- rank(y[ok]); RZ <- as.data.frame(lapply(Z, function(z) rank(z[ok])))
  c(rho = cor(resid(lm(rx ~ ., data = RZ)), resid(lm(ry ~ ., data = RZ))), n = sum(ok)) }
part <- rbindlist(lapply(c("HL", "U_pare", "U_gmuct"), function(v) rbind(
  { r <- psp(S1[[v]], S1$pUP, list(S1$RNA_TPM)); data.table(gene_set = "all single-TuORF genes", measure = MEAS[[v]], versus = "%uORF+", given = "RNA TPM (this study)", n = r[["n"]], partial_rho = r[["rho"]]) },
  { r <- psp(S1[[v]], S1$pUP, list(S1$exon_length)); data.table(gene_set = "all single-TuORF genes", measure = MEAS[[v]], versus = "%uORF+", given = "exon length", n = r[["n"]], partial_rho = r[["rho"]]) },
  { r <- psp(S1[[v]], S1$pUP, list(S1$exon_length, S1$RNA_TPM)); data.table(gene_set = "all single-TuORF genes", measure = MEAS[[v]], versus = "%uORF+", given = "exon length and RNA TPM (this study)", n = r[["n"]], partial_rho = r[["rho"]]) },
  { k <- S1[S_ok == TRUE]; r <- psp(k[[v]], k$eff_incl, list(k$RNA_TPM)); data.table(gene_set = sprintf("at least %d RNA-seq reads across the start codon", MIN_COV), measure = MEAS[[v]], versus = "effective inclusion (%uORF+ x S)", given = "RNA TPM (this study)", n = r[["n"]], partial_rho = r[["rho"]]) })))
.add("  partial Spearman correlations (rank-based):")
for (i in seq_len(nrow(part))) .add(sprintf("    %-24s vs %-34s given %-38s %-48s partial rho = %6.3f, n = %s", part$measure[i], part$versus[i], part$given[i], part$gene_set[i], part$partial_rho[i], big(part$n[i])))
# correlations across both gene sets (single-TuORF and No-TuORF genes): the uncapped ratios with half-life and with TE; every measure with exon length
pool <- rbindlist(lapply(list(c("U_gmuct", "HL"), c("U_pare", "HL"), c("U_gmuct", "TE"), c("U_pare", "TE"), c("HL", "exon_length"), c("U_gmuct", "exon_length"), c("U_pare", "exon_length"), c("U_gmuct_tpm", "exon_length"), c("U_pare_tpm", "exon_length"), c("HL", "RNA_TPM")), function(pr) {
  r <- sp_ci(G[[pr[1]]], G[[pr[2]]]); data.table(x = pr[1], y = pr[2], n = r[["n"]], rho = r[["rho"]], ci_lo = r[["lo"]], ci_hi = r[["hi"]]) }))
NM <- c(HL = "mRNA half-life (h)", U_pare = "uncapped ratio (PARE)", U_gmuct = "uncapped ratio (GMUCT)", TE = "mORF TE", exon_length = "exon length", RNA_TPM = "RNA TPM (this study)",
        U_pare_tpm = "PARE reads per million / RNA TPM of this study", U_gmuct_tpm = "GMUCT reads per million / RNA TPM of this study")
pool[, `:=`(x = NM[x], y = NM[y])]
.add("  Spearman correlations across the single-TuORF and No-TuORF genes together:")
for (i in seq_len(nrow(pool))) .add(sprintf("    %-48s vs %-22s rho = %6.3f (%.3f to %.3f), n = %s", pool$x[i], pool$y[i], pool$rho[i], pool$ci_lo[i], pool$ci_hi[i], big(pool$n[i])))
gdesc <- G[, .(n = .N, median_exon_length = median(exon_length), median_RNA_TPM = median(RNA_TPM, na.rm = TRUE)), by = .(group = as.character(pclass))][match(PUP, group)]
.add("  transcript length (exon length of the gene) and RNA abundance of the groups:")
for (i in seq_len(nrow(gdesc))) .add(sprintf("    %-18s n = %5d, median exon length %.0f nt, median RNA TPM %.1f", gdesc$group[i], gdesc$n[i], gdesc$median_exon_length[i], gdesc$median_RNA_TPM[i]))
# %uORF+ groups against genes without TuORFs: linear models of the log10 value, adjusted for transcript length and RNA abundance
ADJ <- c(HL = "mRNA half-life (h)", U_pare = "uncapped ratio (PARE)", U_gmuct = "uncapped ratio (GMUCT)", U_pare_tpm = "PARE reads per million / RNA TPM of this study", U_gmuct_tpm = "GMUCT reads per million / RNA TPM of this study")
grp_model <- function(v, covars, label) {
  d <- G[is.finite(get(v)) & get(v) > 0 & is.finite(RNA_TPM) & RNA_TPM > 0]
  dd <- data.frame(y = log10(d[[v]]), g = relevel(factor(as.character(d$pclass), PUP), ref = PUP[4]), llen = log10(d$exon_length), lrna = log10(d$RNA_TPM))
  m <- lm(as.formula(paste("y ~ g", if (length(covars)) paste("+", paste(covars, collapse = " + ")) else "")), data = dd); s <- summary(m)$coefficients; ci <- confint(m)
  rbindlist(lapply(PUP[1:3], function(g) { r <- paste0("g", g); data.table(measure = ADJ[[v]], model = label, group = g, n_group = sum(dd$g == g), n_No_TuORF = sum(dd$g == PUP[4]),
    fold_vs_No_TuORF = 10^s[r, "Estimate"], ci_lo = 10^ci[r, 1], ci_hi = 10^ci[r, 2], p = s[r, "Pr(>|t|)"]) })) }
MODELS <- list(c("unadjusted (ratio of geometric means)"), c("adjusted for transcript length", "llen"), c("adjusted for RNA abundance", "lrna"), c("adjusted for transcript length and RNA abundance", "llen", "lrna"))
adj <- rbindlist(lapply(names(ADJ), function(v) rbindlist(lapply(MODELS, function(m) grp_model(v, m[-1], m[1])))))
.add("  %uORF+ groups against genes without TuORFs, linear models of the log10 value (fold of the geometric mean):")
for (v in ADJ) for (g in PUP[1:3]) { a <- adj[measure == v & group == g]
  .add(sprintf("    %-48s %-18s unadjusted %.2fx (p = %.2g); length %.2fx (p = %.2g); abundance %.2fx (p = %.2g); length and abundance %.2fx (%.2f to %.2f; p = %.2g)", v, g, a$fold_vs_No_TuORF[1], a$p[1], a$fold_vs_No_TuORF[2], a$p[2], a$fold_vs_No_TuORF[3], a$p[3], a$fold_vs_No_TuORF[4], a$ci_lo[4], a$ci_hi[4], a$p[4])) }
# %uORF+ groups against genes without TuORFs within thirds of RNA abundance
strata <- rbindlist(lapply(c("HL", "U_pare", "U_gmuct"), function(v) { d <- G[is.finite(get(v)) & is.finite(RNA_TPM)]
  d[, third := cut(RNA_TPM, quantile(RNA_TPM, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("low RNA abundance", "middle RNA abundance", "high RNA abundance"))]
  rbindlist(lapply(levels(d$third), function(t3) { ref <- d[third == t3 & pclass == PUP[4]][[v]]
    rbindlist(lapply(PUP[1:3], function(g) { x <- d[third == t3 & pclass == g][[v]]
      data.table(measure = MEAS[[v]], RNA_abundance_third = t3, group = g, n_group = length(x), n_No_TuORF = length(ref), median_group = median(x), median_No_TuORF = median(ref), median_fold_vs_No_TuORF = median(x) / median(ref),
                 KS_p = suppressWarnings(ks.test(x, ref)$p.value), Wilcoxon_p = wilcox.test(x, ref)$p.value) })) })) }))
.add("  %uORF+ groups against genes without TuORFs within thirds of RNA abundance (thirds of RNA TPM over all genes with the measure): median fold (KS p)")
for (v in MEAS[c("HL", "U_pare", "U_gmuct")]) for (g in PUP[1:3]) { a <- strata[measure == v & group == g]
  .add(sprintf("    %-24s %-18s low %.2fx (p = %.2g, n = %d); middle %.2fx (p = %.2g, n = %d); high %.2fx (p = %.2g, n = %d)", v, g, a$median_fold_vs_No_TuORF[1], a$KS_p[1], a$n_group[1], a$median_fold_vs_No_TuORF[2], a$KS_p[2], a$n_group[2], a$median_fold_vs_No_TuORF[3], a$KS_p[3], a$n_group[3])) }
.add("  Spearman correlations (single-TuORF genes):")
for (i in seq_len(nrow(cor_tab))) .add(sprintf("    %-24s vs %-50s %-48s rho = %6.3f (%.3f to %.3f), p = %.2g, n = %s", cor_tab$measure[i], cor_tab$versus[i], cor_tab$gene_set[i], cor_tab$rho[i], cor_tab$ci_lo[i], cor_tab$ci_hi[i], cor_tab$p[i], big(cor_tab$n[i])))

################################################################################
# 6) Figure S21
################################################################################
th <- theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), legend.position = "bottom", legend.direction = "vertical", legend.key.size = unit(0.3, "cm"),
                                      legend.text = element_text(size = 7), plot.title = element_text(size = 9.5, face = "bold"), plot.subtitle = element_text(size = 7.3), legend.margin = margin(0, 0, 0, 0))
pal_pup <- c("#619CFF", "#00BA38", "#F8766D", "grey40")       # the colours of Figure 2F: 0-30% blue, 30-80% green, 80-100% red
pal_sc  <- c("#F8766D", "#7B3294", "grey40")                  # retained = the colour of the 80-100% group; removed by splicing; No TuORF
lab_plain <- function(x) format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)
ecdfp <- function(d, grp, val, g, xl, title, sub, pal, xlim = NULL, logx = TRUE) {
  d <- as.data.frame(d); d <- d[is.finite(d[[val]]) & !is.na(d[[grp]]) & (!logx | d[[val]] > 0), ]; d$lab <- factor(g$labels[as.character(d[[grp]])], g$labels)
  p <- ggplot(d, aes(x = .data[[val]], colour = lab)) + stat_ecdf(geom = "step", linewidth = 0.8) + scale_colour_manual(values = setNames(pal, g$labels), name = NULL, drop = FALSE) +
    labs(title = title, subtitle = sub, x = xl, y = "Fraction of genes") + th + guides(colour = guide_legend(ncol = 1))
  if (logx) p <- p + scale_x_log10(labels = lab_plain)
  if (!is.null(xlim)) p <- p + coord_cartesian(xlim = xlim)                    # the curves use all genes; only the view is limited
  p
}
HL_lab <- "mRNA half-life (h)"
pA <- ecdfp(Gall, "pclass", "HL", gs_all$HL, HL_lab, "A", "Single-TuORF genes by %uORF+ group", pal_pup, xlim = c(0.4, 10))
pB <- ecdfp(Gclean, "pclass", "HL", gs_cln$HL, HL_lab, "B", "Genes unaffected by leader splicing", pal_pup, xlim = c(0.4, 10))
pC <- ecdfp(Gsc, "sc", "HL", gs_sc$HL, HL_lab, "C", "Genes with 80–100% uORF+", pal_sc, xlim = c(0.4, 10))
pD <- ecdfp(Gsc, "sc", "U_pare", gs_sc$U_pare, "Uncapped ratio (PARE)", "D", "Genes with 80–100% uORF+", pal_sc, xlim = c(0.16, 7.9))
pE <- ecdfp(Gsc, "sc", "U_gmuct", gs_sc$U_gmuct, "Uncapped ratio (GMUCT)", "E", "Genes with 80–100% uORF+", pal_sc, xlim = c(0.1, 50))
pF <- ecdfp(Gsc, "sc", "TE", gs_sc$TE, "mORF TE", "F", "Genes with 80–100% uORF+", pal_sc, xlim = c(0, 3.5), logx = FALSE)
source("revision_modules/fig_compact.R")   # 6.5 x <= 9 in, Arial, >= 8 pt; legends inside the plots (Larry, 2026-11-09)
wrap_count <- function(x) sub(" (n = ", "\n(n = ", x, fixed = TRUE)            # the count and the letter on a second line, so that the key fits in the free corner
cc <- lapply(list(pA, pB, pC, pD, pE, pF), function(p) legend_in(compact(p), 0.99, 0.03, wrap = wrap_count) + theme(legend.key.spacing.y = unit(3.5, "pt")))   # space between the entries, so that each count stays with its group
g <- fig_grid(cc, ncol = 2, top = fig_title("Figure S21. mRNA half-life, uncapped mRNA and mORF TE by uORF inclusion and leader splicing", wrap = 64), h = 8.2)
fig_save(file.path(OUT_DIR, "A31_Figure_decay_by_inclusion_and_splicing.pdf"), g, h = 8.2)

################################################################################
# 7) Tables and summary
################################################################################
genes_out <- G[cohort == "single-TuORF"][order(gene_id), .(gene_id, uIF_class, pct_uORF_plus_group = as.character(pclass), uIF, pct_uORF_plus = pUP, mORF_TE = TE, RNA_TPM_this_study = RNA_TPM,
  half_life_h = HL, uncapped_ratio_PARE = U_pare, uncapped_ratio_GMUCT = U_gmuct, annotation_class = annot_class, leader_splicing_status = splicing_status,
  RNA_reads_across_start_codon = cov, S_fraction_retaining_start_codon = ifelse(S_ok, S, NA_real_), effective_inclusion_pct = eff_incl, fraction_of_transcripts_losing_start_codon_by_splicing = loss_spliced,
  start_codon_group_Fig_S21C_F = as.character(sc))]
fwrite(genes_out, file.path(OUT_DIR, "A31_Table_genes.csv"))
lab_set <- function(x, set) cbind(gene_set = set, x)
gsum <- rbind(
  rbindlist(lapply(names(gs_cln), function(v) cbind(measure = MEAS[[v]], grouping = "%uORF+ group, genes unaffected by leader splicing (Figure S21B for half-life)", gs_cln[[v]]$table))),
  rbindlist(lapply(names(gs_sc),  function(v) cbind(measure = MEAS[[v]], grouping = "start codon retained or removed by leader splicing, genes with 80-100% uORF+ (Figure S21C-F)", gs_sc[[v]]$table))))
pair <- rbind(
  rbindlist(lapply(names(gs_cln), function(v) cbind(measure = MEAS[[v]], grouping = "%uORF+ group, genes unaffected by leader splicing", gs_cln[[v]]$pairs))),
  rbindlist(lapply(names(gs_sc),  function(v) cbind(measure = MEAS[[v]], grouping = "start codon retained or removed by leader splicing, genes with 80-100% uORF+", gs_sc[[v]]$pairs))))
wb <- createWorkbook()
addWorksheet(wb, "README"); writeData(wb, "README", data.frame(x = c(
  "Table S12. mRNA half-life, uncapped mRNA and mORF TE by uORF inclusion and leader splicing (Figure S21).",
  "Measures: mRNA half-life (5-EU pulse-chase labeling, Szabo et al. 2020; hours), uncapped ratio from PARE (Nagarajan et al. 2019; Table S11) and from GMUCT (Yu et al. 2021; Table S10), and mORF TE. The decay measures by %uORF+ group for all single-TuORF genes are in Tables S10 and S11.",
  sprintf("Leader splicing: S is the fraction of RNA-seq reads across the uORF start codon that retain it, used for genes with at least %d such reads (Table S5). A gene is unaffected by leader splicing when neither the annotation nor the junction reads indicate splicing at the TuORF (Supplemental Fig. S7).", MIN_COV),
  sprintf("Start-codon groups (genes with 80-100%% uORF+ and S): 'Start codon retained' = S >= %.2f; 'Removed by splicing' = leader splicing removes the start codon from at least half of the transcripts, %%uORF+/100 x (1 - S) >= %.1f (the rule of Figure 2D). Genes between the two groups are not used.", S_RETAINED, LOSS_SPLICED),
  "Sheets: group_summaries (n, median, quartiles, median fold relative to genes without TuORFs with a bootstrap 95% interval from 1,000 resamples, letters from pairwise Kolmogorov–Smirnov (KS) tests at p < 0.05); pairwise_tests (Kolmogorov–Smirnov and Wilcoxon rank-sum tests);",
  "removed_vs_retained (the two start-codon groups compared directly: medians, median fold with bootstrap interval, tests, the Kolmogorov–Smirnov p values of the three decay measures after Holm correction, and the folds of the geometric mean from linear models of the log10 value: unadjusted, adjusted for log10 RNA TPM, for log10 exon length, and for both; columns ending in _lo and _hi = limits of the 95% intervals; adjA_, adjL_ and adjLA_ = adjusted for RNA abundance, for transcript length, and for both); start_codon_groups (description of the two groups);",
  "adjusted_group_comparisons (each %uORF+ group against genes without TuORFs: fold of the geometric mean from linear models of the log10 value, unadjusted and adjusted for transcript length, for RNA abundance, and for both; ci_lo and ci_hi = limits of the 95% confidence interval of the fold; transcript length is the exon length of the gene, that is, the length of the union of its annotated exons); RNA_abundance_thirds (the same groups within thirds of RNA TPM); group_length_abundance (median exon length and RNA TPM of the groups); correlations_all_genes (single-TuORF and No-TuORF genes together).",
  { .r12 <- function(m) adj[measure == m & model == "adjusted for transcript length" & group == PUP[3]]; .p12 <- .r12("PARE reads per million / RNA TPM of this study"); .g12 <- .r12("GMUCT reads per million / RNA TPM of this study"); stopifnot(nrow(.p12) == 1, nrow(.g12) == 1)   # V53 (2026-11-22): the values of this note are taken from the table (they were typed: 0.93, 0.005, 0.22)
    sprintf("Note on transcript length: the uncapped ratio divides 5'-end reads per million by RNA-seq reads per million, and RNA-seq reads per gene increase with transcript length. With RNA TPM of this study as the denominator (rows 'reads per million / RNA TPM'), the ratio increases with transcript length instead. After adjustment for transcript length, this second ratio is %.2f times (PARE, p = %s) and %.2f times (GMUCT, p = %s) that of genes without TuORFs in the genes with 80-100%% uORF+ (sheet adjusted_group_comparisons), which is the direction of the uncapped ratio.", .p12$fold_vs_No_TuORF, format(signif(.p12$p, 2)), .g12$fold_vs_No_TuORF, format(signif(.g12$p, 2))) },
  "correlations (Spearman rho with bootstrap 95% intervals: the decay measures against %uORF+, against S and against the effective inclusion, %uORF+ x S); partial_correlations (rank-based partial correlations given RNA abundance, exon length, or both); genes (per-gene values of the single-TuORF genes).",
  "Note on the adjustment for RNA abundance: the steady-state abundance of an mRNA depends on its half-life, so an adjustment for abundance can remove part of a true difference in half-life.",
  "Code: revision_modules/A31_decay_by_inclusion_and_splicing.R in the Supplemental Code (reads the per-gene tables of A2, A11 and A12).")), colNames = FALSE)
addWorksheet(wb, "group_summaries"); writeData(wb, "group_summaries", gsum)
addWorksheet(wb, "pairwise_tests"); writeData(wb, "pairwise_tests", pair)
addWorksheet(wb, "removed_vs_retained"); writeData(wb, "removed_vs_retained", direct)
addWorksheet(wb, "start_codon_groups"); writeData(wb, "start_codon_groups", desc)
addWorksheet(wb, "correlations"); writeData(wb, "correlations", cor_tab)
addWorksheet(wb, "partial_correlations"); writeData(wb, "partial_correlations", part)
addWorksheet(wb, "correlations_all_genes"); writeData(wb, "correlations_all_genes", pool)
addWorksheet(wb, "adjusted_group_comparisons"); writeData(wb, "adjusted_group_comparisons", adj)
addWorksheet(wb, "RNA_abundance_thirds"); writeData(wb, "RNA_abundance_thirds", strata)
addWorksheet(wb, "group_length_abundance"); writeData(wb, "group_length_abundance", gdesc)
addWorksheet(wb, "genes"); writeData(wb, "genes", genes_out)
saveWorkbook(wb, file.path(OUT_DIR, "A31_Table_decay_inclusion_splicing.xlsx"), overwrite = TRUE)
# the same sheets as CSV files (read by the checks of the manuscript text: weak_uORF_revision/revision_docx/verify_v38.py)
for (nm in c("gsum", "pair", "direct", "desc", "cor_tab", "part", "pool", "adj", "strata", "gdesc")) fwrite(get(nm), file.path(OUT_DIR, sprintf("A31_csv_%s.csv", nm)))
.add(sprintf("A31 done in %.1f min; outputs: output/A31_Figure_decay_by_inclusion_and_splicing.pdf, A31_Table_decay_inclusion_splicing.xlsx, A31_Table_genes.csv, A31_summary.txt", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
writeLines(.log, file.path(OUT_DIR, "A31_summary.txt"))
