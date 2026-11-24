################################################################################
# A37 — Are the genes of the only-low %uORF+ group highly expressed? (Larry,    #
# 2026-11-20: "I want some highly expressed set of genes", as a replacement for #
# "housekeeping genes" in the Abstract). Not part of the manuscript.            #
# Expression = our seedling RNA-seq: gene TPM = sum of the expressed isoforms   #
# (RSEM TPM, mean of the six libraries; at most two isoforms per gene; the pool #
# of Figure 3C, output/Table_Transcriptome_Impact_PerTranscript_CAGE.csv).      #
# "Highly expressed" = top 10% (also top 5% and top 20%) of all expressed       #
# nuclear genes with a 5' leader. Groups as in Supplemental Figs. S15-S17       #
# (4,886 TuORF genes with a GO term; highest %uORF+ among a gene's TuORFs).     #
# Also: the same in the Klepikova 2016 atlas (median TPM over 54 organs and     #
# stages; output/A35_gene_tau.csv), and the GO plastid / cytosol enrichment of  #
# the only-low group after adjustment for expression level (logistic models).   #
# Output: output/A37_low_group_expression_summary.txt                           #
# Run from Weak_uORF_submission/:  Rscript revision_modules/A37_low_group_expression_level.R
################################################################################
suppressPackageStartupMessages({ library(openxlsx); library(dplyr); library(tidyr); library(GO.db); library(AnnotationDbi) })
options(stringsAsFactors = FALSE, warn = 1)
OUT_DIR <- "output"; PATH_GOSLIM <- "Files/ATH_GO_GOSLIM.tsv"
.log <- character(0); .add <- function(...) { s <- paste0(...); message(s); .log <<- c(.log, s) }
.add(sprintf("A37 expression level of the %%uORF+ groups — %s", format(Sys.time(), "%Y-%m-%d %H:%M")))
# ---- annotation, universe and groups exactly as A23 ------------------------------------------------------------------------
raw_go <- read.delim(PATH_GOSLIM, header = FALSE, sep = "\t", quote = "")
gene_col_idx <- which.max(apply(raw_go[1:100, ], 2, function(x) mean(grepl("^AT[1-5MC]G", x, ignore.case = TRUE))))
go_col_idx   <- which.max(apply(raw_go[1:100, ], 2, function(x) mean(grepl("GO:", x))))
go_long <- data.frame(gene_raw = trimws(raw_go[, gene_col_idx]), go_string = trimws(raw_go[, go_col_idx])) %>%
  separate_rows(go_string, sep = ",") %>% mutate(gene = toupper(gsub("\\.\\d+$", "", gene_raw)), go_id = trimws(go_string)) %>%
  dplyr::select(gene, go_id) %>% filter(grepl("^GO:", go_id)) %>% distinct()
go_info <- AnnotationDbi::select(GO.db, keys = unique(go_long$go_id), columns = c("TERM", "ONTOLOGY"), keytype = "GOID")
go_full <- left_join(go_long, go_info, by = c("go_id" = "GOID")) %>% filter(!is.na(ONTOLOGY))
S1A <- read.xlsx(file.path(OUT_DIR, "Table_S1_all_TuORFs.xlsx"), sheet = "Table S1A (CAGE)")
universe_genes <- intersect(unique(toupper(S1A$gene_id[is.finite(S1A$uIF)])), unique(go_full$gene))
g_all <- S1A %>% filter(is.finite(uIF)) %>% transmute(gene = toupper(gene_id), pUP = as.numeric(pct_uORF_plus)) %>%
  group_by(gene) %>% summarise(max_pUP = max(pUP, na.rm = TRUE), .groups = "drop") %>% filter(is.finite(max_pUP), gene %in% universe_genes) %>%
  mutate(group = case_when(max_pUP < 30 ~ "only low (<30%)", max_pUP < 80 ~ "at least one medium (30-80%)", TRUE ~ "at least one high (>=80%)"))
lev <- c("only low (<30%)", "at least one medium (30-80%)", "at least one high (>=80%)")
stopifnot(length(universe_genes) == 4886L, all(as.integer(table(factor(g_all$group, lev))) == c(1046L, 864L, 2976L)))   # V53 (2026-11-22): group sizes with the corrected CAGE positions (V52: 1,083 / 885 / 2,918)
# ---- seedling gene TPM and the highly expressed sets -----------------------------------------------------------------------
# 2026-11-20 (number check against the Methods): the pool of Figure 3C as module A26 defines it (nuclear transcripts with a leader from
# splicing_rule.R::leader_lengths(); 20,934 transcripts, 18,359 genes), not the Category column of the table (20,899 transcripts, 18,275 genes)
suppressPackageStartupMessages(library(data.table)); source("revision_modules/splicing_rule.R")
LL <- leader_lengths(); leader_tx <- LL[leader_len > 0, transcript_id]
tx <- read.csv(file.path(OUT_DIR, "Table_Transcriptome_Impact_PerTranscript_CAGE.csv"))
gt <- tx %>% filter(grepl("^Nuclear", Category), transcript_id %in% leader_tx) %>% group_by(gene = gene_id) %>% summarise(TPM = sum(TPM), .groups = "drop")
stopifnot(nrow(gt) == 18359L, abs(sum(gt$TPM) - 635103) < 1)
cut <- quantile(gt$TPM, c(0.80, 0.90, 0.95))
.add(sprintf("expressed nuclear genes with a 5' leader: %d; top 20%% from %.1f TPM, top 10%% from %.1f TPM, top 5%% from %.1f TPM", nrow(gt), cut[1], cut[2], cut[3]))
d <- left_join(g_all, gt, by = "gene"); .add(sprintf("TuORF GO genes with seedling TPM: %d of %d", sum(!is.na(d$TPM)), nrow(d))); d <- d[!is.na(d$TPM), ]
d$top20 <- d$TPM >= cut[1]; d$top10 <- d$TPM >= cut[2]; d$top5 <- d$TPM >= cut[3]
gs <- function(z, lab) sprintf("%-38s n = %4d | median %.1f TPM | top 5%%: %4.1f%% | top 10%%: %4.1f%% | top 20%%: %4.1f%% | share of the group's mRNA in top-10%% genes: %4.1f%%",
                               lab, nrow(z), median(z$TPM), 100 * mean(z$top5), 100 * mean(z$top10), 100 * mean(z$top20), 100 * sum(z$TPM[z$top10]) / sum(z$TPM))
.add(""); .add(gs(gt %>% mutate(top5 = TPM >= cut[3], top10 = TPM >= cut[2], top20 = TPM >= cut[1]), "all expressed genes with a leader"))
.add(gs(d, "all TuORF GO genes")); for (g in lev) .add(gs(d[d$group == g, ], g))
lo <- d$group == lev[1]; hi <- d$group == lev[3]
for (k in c("top5", "top10", "top20")) {
  f1 <- fisher.test(table(factor(lo, c(TRUE, FALSE)), factor(d[[k]], c(TRUE, FALSE))))
  f2 <- fisher.test(table(factor(d$group[lo | hi] == lev[1], c(TRUE, FALSE)), factor(d[[k]][lo | hi], c(TRUE, FALSE))))
  .add(sprintf("%s: only low vs the other TuORF genes odds ratio %.2f (95%% CI %.2f-%.2f, p = %.2g); vs the high group %.2f (%.2f-%.2f, p = %.2g)",
               k, f1$estimate, f1$conf.int[1], f1$conf.int[2], f1$p.value, f2$estimate, f2$conf.int[1], f2$conf.int[2], f2$p.value))
}
.add(sprintf("TPM, only low vs the other TuORF genes: Wilcoxon p = %.2g; median ratio %.2f", wilcox.test(d$TPM[lo], d$TPM[!lo])$p.value, median(d$TPM[lo]) / median(d$TPM[!lo])))
.add(sprintf("Spearman rho (highest %%uORF+ of the gene, seedling TPM) over the TuORF GO genes: %.3f", cor(d$max_pUP, d$TPM, method = "spearman")))
# what share of the top-10% TuORF genes is in each group
t10 <- d[d$top10, ]; .add(sprintf("TuORF GO genes in the top 10%%: %d; in only low %.1f%%, medium %.1f%%, high %.1f%% (all TuORF GO genes: %.1f%% / %.1f%% / %.1f%%)", nrow(t10),
                                  100 * mean(t10$group == lev[1]), 100 * mean(t10$group == lev[2]), 100 * mean(t10$group == lev[3]), 100 * mean(lo), 100 * mean(d$group == lev[2]), 100 * mean(hi)))
# ---- independent atlas (Klepikova 2016): median TPM over 54 organs and stages ---------------------------------------------
at <- read.csv(file.path(OUT_DIR, "A35_gene_tau.csv"))[, c("gene", "median_tpm", "max_tpm")]
d2 <- inner_join(d, at, by = "gene"); c10 <- quantile(d2$median_tpm, 0.9)
f3 <- fisher.test(table(factor(d2$group == lev[1], c(TRUE, FALSE)), factor(d2$median_tpm >= c10, c(TRUE, FALSE))))
.add(""); .add(sprintf("atlas (n = %d TuORF GO genes): median of the 54-group medians: only low %.1f, medium %.1f, high %.1f TPM; top 10%% of TuORF genes by atlas median: only low %.1f%% vs others %.1f%%, odds ratio %.2f (%.2f-%.2f, p = %.2g)",
                   nrow(d2), median(d2$median_tpm[d2$group == lev[1]]), median(d2$median_tpm[d2$group == lev[2]]), median(d2$median_tpm[d2$group == lev[3]]),
                   100 * mean(d2$median_tpm[d2$group == lev[1]] >= c10), 100 * mean(d2$median_tpm[d2$group != lev[1]] >= c10), f3$estimate, f3$conf.int[1], f3$conf.int[2], f3$p.value))
# ---- does expression level explain the plastid and cytosol enrichment? ----------------------------------------------------
gl <- go_full %>% filter(gene %in% universe_genes)
offs <- function(id, tab) c(id, unlist(AnnotationDbi::mget(id, tab), use.names = FALSE))
d$plastid <- d$gene %in% gl$gene[gl$go_id %in% offs("GO:0009536", GOCCOFFSPRING)]
d$cytosol <- d$gene %in% gl$gene[gl$go_id %in% offs("GO:0005829", GOCCOFFSPRING)]
d$lo <- lo; d$lTPM <- log10(d$TPM + 1)
.add("")
for (k in c("plastid", "cytosol")) {
  m0 <- glm(reformulate("lo", k), family = binomial, data = d); m1 <- glm(reformulate(c("lo", "lTPM"), k), family = binomial, data = d)
  .add(sprintf("%s term, only low vs others: odds ratio %.2f unadjusted, %.2f adjusted for seedling TPM (95%% CI %.2f-%.2f, p = %.2g)", k, exp(coef(m0)[["loTRUE"]]), exp(coef(m1)[["loTRUE"]]),
               exp(confint.default(m1)["loTRUE", 1]), exp(confint.default(m1)["loTRUE", 2]), summary(m1)$coefficients["loTRUE", 4]))
}
# ---- are the plastid genes of the only-low group highly expressed? (Larry's wording "plastid and other highly expressed genes") --
.add("")
ps <- function(z, lab) .add(sprintf("%-44s n = %4d | median %.1f TPM | top 5%%: %4.1f%% | top 10%%: %4.1f%% | top 20%%: %4.1f%% | above the median of expressed genes (%.1f TPM): %4.1f%%",
                                    lab, nrow(z), median(z$TPM), 100 * mean(z$top5), 100 * mean(z$top10), 100 * mean(z$top20), median(gt$TPM), 100 * mean(z$TPM > median(gt$TPM))))
ps(d[d$lo & d$plastid, ], "only low, plastid term"); ps(d[d$lo & !d$plastid, ], "only low, no plastid term")
ps(d[!d$lo & d$plastid, ], "medium + high, plastid term"); ps(d[!d$lo & !d$plastid, ], "medium + high, no plastid term")
.add(sprintf("only-low genes in the top 10%%: %d, of which %d (%.1f%%) have a plastid term", sum(d$lo & d$top10), sum(d$lo & d$top10 & d$plastid), 100 * mean(d$plastid[d$lo & d$top10])))
.add(sprintf("only-low genes with a plastid term or in the top 10%%: %d of %d (%.1f%%); other TuORF genes: %.1f%%", sum(d$lo & (d$plastid | d$top10)), sum(d$lo), 100 * mean(d$plastid[d$lo] | d$top10[d$lo]), 100 * mean(d$plastid[!d$lo] | d$top10[!d$lo])))
# ---- second TSS method: the same groups from TSS-seq %uORF+ (Table S1B) -----------------------------------------------------
S1B <- read.xlsx(file.path(OUT_DIR, "Table_S1_all_TuORFs.xlsx"), sheet = "Table S1B (TSS-seq)")
gB <- S1B %>% filter(is.finite(uIF)) %>% transmute(gene = toupper(gene_id), pUP = as.numeric(pct_uORF_plus)) %>% filter(is.finite(pUP)) %>%
  group_by(gene) %>% summarise(max_pUP = max(pUP), .groups = "drop") %>% filter(gene %in% universe_genes) %>%
  mutate(group = case_when(max_pUP < 30 ~ lev[1], max_pUP < 80 ~ lev[2], TRUE ~ lev[3])) %>% inner_join(gt, by = "gene")
gB$top10 <- gB$TPM >= cut[2]; loB <- gB$group == lev[1]
fB <- fisher.test(table(factor(loB, c(TRUE, FALSE)), factor(gB$top10, c(TRUE, FALSE))))
.add(""); .add(sprintf("TSS-seq groups (n = %d; only low %d, medium %d, high %d): top 10%% in only low %.1f%%, medium %.1f%%, high %.1f%%; only low vs others odds ratio %.2f (%.2f-%.2f, p = %.2g); median TPM only low %.1f vs others %.1f",
                   nrow(gB), sum(loB), sum(gB$group == lev[2]), sum(gB$group == lev[3]), 100 * mean(gB$top10[loB]), 100 * mean(gB$top10[gB$group == lev[2]]), 100 * mean(gB$top10[gB$group == lev[3]]),
                   fB$estimate, fB$conf.int[1], fB$conf.int[2], fB$p.value, median(gB$TPM[loB]), median(gB$TPM[!loB])))
# (V53, 2026-11-23) the values of this note are read from output/A6_summary.txt (run module A6 first); until V52 they were typed here and went out of date
a6 <- readLines(file.path(OUT_DIR, "A6_summary.txt"))
.a6v <- function(pattern) { l <- grep(pattern, a6, value = TRUE); stopifnot(length(l) == 1); as.numeric(sub(".*<= 30% ([0-9.]+)%.*", "\\1", l)) }
a6_full <- .a6v("^  all TuORFs \\(n = "); a6_ds <- vapply(c(20, 50, 100, 200), function(k) .a6v(sprintf("all TuORFs downsampled to %d tags", k)), numeric(1))
.add(sprintf("Read depth: module A6 downsampled every TuORF to a fixed number of CAGE tags; the share at <= 30%% uORF+ is %.1f%% at full depth and %.1f%%, %.1f%%, %.1f%% and %.1f%% at 20, 50, 100 and 200 tags (output/A6_summary.txt), so the low group is not a depth effect.",
             a6_full, a6_ds[1], a6_ds[2], a6_ds[3], a6_ds[4]))
writeLines(.log, file.path(OUT_DIR, "A37_low_group_expression_summary.txt")); .add("written: output/A37_low_group_expression_summary.txt")
