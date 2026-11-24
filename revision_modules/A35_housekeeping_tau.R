################################################################################
# A35 — Are the genes of the only-low %uORF+ group housekeeping genes?          #
# (Larry, 2026-11-20: "Can we say the genes in the low group are housekeeping   #
#  genes?")  Not part of the manuscript.                                          #
# Housekeeping = expressed in all organs and stages at a similar level. Gene    #
# Ontology terms do not test this, so this module measures expression breadth   #
# in an independent atlas: Klepikova et al. 2016 (Expression Atlas E-MTAB-7978, #
# 54 organ / developmental-stage groups, median TPM per group).                 #
#   tau   = Yanai tissue-specificity index over the 54 groups on log2(TPM + 1)  #
#           (0 = same level everywhere, 1 = one group only)                      #
#   broad = TPM >= 1 in all 54 groups                                            #
#   housekeeping-like = broad and tau < 0.2                                      #
# Groups: the Gene Ontology grouping of Supplemental Figs. S15-S17 (A23/A24):   #
# all TuORF genes with a GO term (n = 4,886) by the highest %uORF+ among their   #
# TuORFs: only low (< 30), at least one medium (30-80), at least one high (>=80) #
# Within the only-low group, genes are also split by the significant GO terms   #
# (plastid/photosynthesis, glycolysis/folding/cytosol, cell wall/secreted,      #
# stress responses).                                                             #
# Outputs: output/A35_housekeeping_tau_summary.txt, output/A35_gene_tau.csv,     #
#          output/A35_tau_by_group.pdf                                           #
# Run from Weak_uORF_submission/:  Rscript revision_modules/A35_housekeeping_tau.R
################################################################################
suppressPackageStartupMessages({ library(openxlsx); library(dplyr); library(tidyr); library(GO.db); library(AnnotationDbi); library(ggplot2) })
options(stringsAsFactors = FALSE, warn = 1)
OUT_DIR <- "output"; PATH_GOSLIM <- "Files/ATH_GO_GOSLIM.tsv"
ATLAS <- "/Users/wu/Desktop/Desktop_2026_July/DOC/E-MTAB-7978-query-results.tpms.tsv"
.log <- character(0); .add <- function(...) { s <- paste0(...); message(s); .log <<- c(.log, s) }
.add(sprintf("A35 housekeeping check (tau in the Klepikova 2016 atlas) — %s", format(Sys.time(), "%Y-%m-%d %H:%M")))

# ---- GO annotation and universe, exactly as A23 ---------------------------------------------------------------------------
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
n_grp <- table(g_all$group)
.add("groups: ", paste(sprintf("%s n=%d", names(n_grp), as.integer(n_grp)), collapse = "; "))
stopifnot(length(universe_genes) == 4886L, n_grp[["only low (<30%)"]] == 1046L, n_grp[["at least one medium (30-80%)"]] == 864L, n_grp[["at least one high (>=80%)"]] == 2976L)   # V53 (2026-11-22): group sizes with the corrected CAGE positions (V52: 1,083 / 885 / 2,918)

# ---- atlas: tau and breadth ----------------------------------------------------------------------------------------------
a <- read.delim(ATLAS, comment.char = "#", check.names = FALSE)
stopifnot(names(a)[1] == "Gene ID", ncol(a) == 56L)
a <- a[grepl("^AT[1-5CM]G\\d{5}$", a[["Gene ID"]]), ]; stopifnot(!anyDuplicated(a[["Gene ID"]]))
M <- as.matrix(a[, 3:56]); storage.mode(M) <- "double"; M[is.na(M)] <- 0   # blank = not detected
X <- log2(M + 1); xmax <- apply(X, 1, max)
tau <- ifelse(xmax > 0, rowSums(1 - X / xmax) / (ncol(X) - 1), NA_real_)
atl <- data.frame(gene = a[["Gene ID"]], tau = tau, n_ge1 = rowSums(M >= 1), max_tpm = apply(M, 1, max), median_tpm = apply(M, 1, median))
atl$broad <- atl$n_ge1 == 54L; atl$hk_like <- atl$broad & atl$tau < 0.2
.add(sprintf("atlas: %d genes, 54 groups; all atlas genes with TPM >= 1 somewhere: median tau %.2f, broad %.1f%%, housekeeping-like %.1f%%",
             nrow(atl), median(atl$tau[atl$max_tpm >= 1]), 100 * mean(atl$broad[atl$max_tpm >= 1]), 100 * mean(atl$hk_like[atl$max_tpm >= 1])))
d <- left_join(g_all, atl, by = "gene"); .add(sprintf("TuORF GO genes found in the atlas: %d of %d", sum(!is.na(d$tau)), nrow(d)))
d <- d[!is.na(d$tau), ]

# ---- low-group genes by the significant GO term categories ----------------------------------------------------------------
cat_terms <- list(
  "plastid / photosynthesis" = c("GO:0009507", "GO:0009570", "GO:0009579", "GO:0009941", "GO:0009534", "GO:0010319", "GO:0009535", "GO:0009536", "GO:0016168", "GO:0015995"),
  "cytosol, glycolysis, folding" = c("GO:0005829", "GO:0006096", "GO:0006002", "GO:0003872", "GO:0006457", "GO:0051082", "GO:0044183", "GO:0051287", "GO:0030170"),
  "cell wall / extracellular" = c("GO:0048046", "GO:0005576", "GO:0009505", "GO:0099503"),
  "stress (cold, wounding, hypoxia)" = c("GO:0009409", "GO:0009611", "GO:0071456"))
gl <- go_full %>% filter(gene %in% universe_genes) %>% distinct(gene, go_id)

summ <- function(z, label) data.frame(set = label, n = nrow(z), median_tau = median(z$tau), q25 = quantile(z$tau, 0.25), q75 = quantile(z$tau, 0.75),
                                     pct_broad = 100 * mean(z$broad), pct_hk_like = 100 * mean(z$hk_like), median_max_tpm = median(z$max_tpm), row.names = NULL)
lev <- c("only low (<30%)", "at least one medium (30-80%)", "at least one high (>=80%)")
res <- rbind(summ(d, "all TuORF GO genes (background)"), do.call(rbind, lapply(lev, function(g) summ(d[d$group == g, ], g))))
low <- d[d$group == lev[1], ]
for (k in names(cat_terms)) {
  gk <- unique(gl$gene[gl$go_id %in% cat_terms[[k]]])
  res <- rbind(res, summ(low[low$gene %in% gk, ], paste0("  only low, ", k)), summ(d[d$gene %in% gk & d$group != lev[1], ], paste0("  medium+high, ", k)))
}
.add(""); .add(sprintf("%-46s %5s %8s %13s %8s %8s %10s", "set", "n", "tau", "tau IQR", "broad%", "hk-like%", "max TPM"))
for (i in seq_len(nrow(res))) .add(sprintf("%-46s %5d %8.2f   %.2f-%.2f %8.1f %8.1f %10.0f", res$set[i], res$n[i], res$median_tau[i], res$q25[i], res$q75[i], res$pct_broad[i], res$pct_hk_like[i], res$median_max_tpm[i]))

# ---- tests: only low versus the other TuORF GO genes, and versus the high group -----------------------------------------
lo <- d$group == lev[1]; hi <- d$group == lev[3]
w1 <- wilcox.test(d$tau[lo], d$tau[!lo]); w2 <- wilcox.test(d$tau[lo], d$tau[hi])
f1 <- fisher.test(table(factor(lo, c(TRUE, FALSE)), factor(d$hk_like, c(TRUE, FALSE))))
f2 <- fisher.test(table(factor(d$group[lo | hi] == lev[1], c(TRUE, FALSE)), factor(d$hk_like[lo | hi], c(TRUE, FALSE))))
f3 <- fisher.test(table(factor(lo, c(TRUE, FALSE)), factor(d$broad, c(TRUE, FALSE))))
.add(""); .add(sprintf("tau, only low vs the other TuORF GO genes: Wilcoxon p = %.2g; vs the high group: p = %.2g", w1$p.value, w2$p.value))
.add(sprintf("housekeeping-like, only low vs the other TuORF GO genes: odds ratio %.2f (95%% CI %.2f-%.2f), p = %.2g", f1$estimate, f1$conf.int[1], f1$conf.int[2], f1$p.value))
.add(sprintf("housekeeping-like, only low vs the high group: odds ratio %.2f (95%% CI %.2f-%.2f), p = %.2g", f2$estimate, f2$conf.int[1], f2$conf.int[2], f2$p.value))
.add(sprintf("broad (TPM >= 1 in all 54 groups), only low vs the other TuORF GO genes: odds ratio %.2f (95%% CI %.2f-%.2f), p = %.2g", f3$estimate, f3$conf.int[1], f3$conf.int[2], f3$p.value))
# expression level as a possible confounder: logistic model of housekeeping-like on group with max TPM
m <- glm(hk_like ~ lo + log10(max_tpm + 1), family = binomial, data = transform(d, lo = group == lev[1]))
ci <- suppressMessages(exp(confint.default(m)["loTRUE", ]))
.add(sprintf("housekeeping-like, only low vs the others, adjusted for the highest atlas TPM: odds ratio %.2f (95%% CI %.2f-%.2f), p = %.2g",
             exp(coef(m)[["loTRUE"]]), ci[1], ci[2], summary(m)$coefficients["loTRUE", 4]))

write.csv(d, file.path(OUT_DIR, "A35_gene_tau.csv"), row.names = FALSE)
pd <- d; pd$group <- factor(pd$group, levels = lev)
p <- ggplot(pd, aes(tau, colour = group)) + stat_ecdf(linewidth = 0.6) + geom_vline(xintercept = 0.2, linetype = 3) +
  scale_colour_manual(values = c("#D55E00", "#999999", "#0072B2"), name = NULL) +
  labs(x = "tissue-specificity index tau (Klepikova 2016 atlas, 54 organs and stages)", y = "cumulative fraction of genes",
       title = "Expression breadth of TuORF genes by %uORF+ group", subtitle = "groups as in Supplemental Figs. S15-S17; dotted line: tau = 0.2") +
  theme_bw(base_size = 9) + theme(legend.position = c(0.72, 0.25), legend.background = element_blank(), panel.grid.minor = element_blank())
ggsave(file.path(OUT_DIR, "A35_tau_by_group.pdf"), p, width = 5.5, height = 4)
writeLines(.log, file.path(OUT_DIR, "A35_housekeeping_tau_summary.txt"))
.add("written: output/A35_housekeeping_tau_summary.txt, output/A35_gene_tau.csv, output/A35_tau_by_group.pdf")
