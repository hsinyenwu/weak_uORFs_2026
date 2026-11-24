################################################################################
# A23 — Gene Ontology enrichment for the pooled medium + high %uORF+ group       #
# (Larry's question, 2026-11-05: what if the medium and high %uORF+ groups are   #
# put together?)                                                                #
#                                                                               #
# Same test as the canonical script (20260827_weak_uORFs.R, GO section):        #
# TAIR GO annotations (Files/ATH_GO_GOSLIM.tsv), universe = genes with a        #
# detected TuORF and a GO term, terms with >= 3 background genes, one-sided     #
# Fisher's exact test, Benjamini-Hochberg within each group. The three          #
# %uORF+ groups of the single-TuORF cohort are recomputed as a check against    #
# output/GO_enrichment_results_by_pct_uORF_plus.csv, and the pooled group       #
# ">= 30% uORF+" (medium + high) is added.                                      #
# Inputs: output/Table_S1_all_TuORFs.xlsx (universe), output/Table_S2_single_   #
# (single-TuORF cohort with pct_uORF_plus). Run from Weak_uORF_submission/:     #
#   Rscript revision_modules/A23_GO_pooled_mid_high.R                           #
# Not part of the manuscript unless Larry decides otherwise.                    #
################################################################################
suppressPackageStartupMessages({ library(openxlsx); library(dplyr); library(tidyr); library(GO.db); library(AnnotationDbi) })
options(stringsAsFactors = FALSE, warn = 1)
OUT_DIR <- "output"; PATH_GOSLIM <- "Files/ATH_GO_GOSLIM.tsv"; FDR_CUT <- 0.05
ROOT <- c("GO:0008150", "GO:0003674", "GO:0005575")
.log <- character(0); .add <- function(...) { s <- paste0(...); message(s); .log <<- c(.log, s) }
.add(sprintf("A23 GO, pooled medium + high %%uORF+ — %s", format(Sys.time(), "%Y-%m-%d %H:%M")))

# ---- GO annotation, as in the canonical script ---------------------------------------------------------------------------
raw_go <- read.delim(PATH_GOSLIM, header = FALSE, sep = "\t", quote = "")
gene_col_idx <- which.max(apply(raw_go[1:100, ], 2, function(x) mean(grepl("^AT[1-5MC]G", x, ignore.case = TRUE))))
go_col_idx   <- which.max(apply(raw_go[1:100, ], 2, function(x) mean(grepl("GO:", x))))
go_long <- data.frame(gene_raw = trimws(raw_go[, gene_col_idx]), go_string = trimws(raw_go[, go_col_idx])) %>%
  separate_rows(go_string, sep = ",") %>% mutate(gene = toupper(gsub("\\.\\d+$", "", gene_raw)), go_id = trimws(go_string)) %>%
  dplyr::select(gene, go_id) %>% filter(grepl("^GO:", go_id)) %>% distinct()
go_info <- AnnotationDbi::select(GO.db, keys = unique(go_long$go_id), columns = c("TERM", "ONTOLOGY"), keytype = "GOID")
go_full <- left_join(go_long, go_info, by = c("go_id" = "GOID")) %>% filter(!is.na(ONTOLOGY))

# ---- universe and groups ---------------------------------------------------------------------------------------------------
S1A <- read.xlsx(file.path(OUT_DIR, "Table_S1_all_TuORFs.xlsx"), sheet = "Table S1A (CAGE)")
S3  <- read.xlsx(file.path(OUT_DIR, "Table_S2_single_TuORF_genes.xlsx"), sheet = 1)   # the single-TuORF cohort table of the pipeline output (= Supplemental Table S3)
all_uorf_genes <- unique(toupper(S1A$gene_id[is.finite(S1A$uIF)]))
universe_genes <- intersect(all_uorf_genes, unique(go_full$gene)); universe_go <- go_full %>% filter(gene %in% universe_genes)
.add(sprintf("universe (genes with a detected TuORF and a GO term): %d (canonical run: %s)", length(universe_genes),
             if (file.exists(file.path(OUT_DIR, "GO_universe_size.txt"))) readLines(file.path(OUT_DIR, "GO_universe_size.txt"))[1] else "4886"))
coh <- S3 %>% transmute(gene = toupper(gene_id), pUP = as.numeric(pct_uORF_plus), uIF = as.numeric(uIF)) %>% filter(is.finite(uIF)) %>% distinct(gene, .keep_all = TRUE)
grp <- list("1-uORF 0~30% uORF+"   = coh$gene[coh$pUP < 30],
            "1-uORF 30~80% uORF+"  = coh$gene[coh$pUP >= 30 & coh$pUP < 80],
            "1-uORF 80~100% uORF+" = coh$gene[coh$pUP >= 80],
            "1-uORF >=30% uORF+ (medium + high pooled)" = coh$gene[coh$pUP >= 30])
grp <- lapply(grp, function(g) intersect(unique(g), universe_genes))
.add(paste("group sizes:", paste(sprintf("%s n=%d", names(grp), sapply(grp, length)), collapse = "; ")))

# ---- Fisher test, as in the canonical script -------------------------------------------------------------------------------
run_fisher <- function(group_genes) {
  N <- length(universe_genes); n <- length(group_genes)
  bg <- universe_go %>% distinct(gene, go_id, TERM, ONTOLOGY) %>% count(go_id, TERM, ONTOLOGY, name = "bg_count") %>% filter(bg_count >= 3)
  fg <- universe_go %>% filter(gene %in% group_genes) %>% distinct(gene, go_id) %>% count(go_id, name = "fg_count")
  bg2 <- left_join(bg, fg, by = "go_id") %>% mutate(fg_count = coalesce(fg_count, 0L)) %>% filter(fg_count > 0)
  res <- lapply(seq_len(nrow(bg2)), function(i) { a <- bg2$fg_count[i]; b <- bg2$bg_count[i] - a; c <- n - a; d <- (N - n) - b
    ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
    data.frame(GO_ID = bg2$go_id[i], Term = bg2$TERM[i], Ontology = bg2$ONTOLOGY[i], FG_Count = a, BG_Count = bg2$bg_count[i], p_value = ft$p.value, Odds_Ratio = unname(ft$estimate)) })
  bind_rows(res) %>% mutate(p_adj_BH = p.adjust(p_value, method = "BH"), Group_size = n)
}
res <- bind_rows(lapply(names(grp), function(g) run_fisher(grp[[g]]) %>% mutate(Group = g))) %>% dplyr::select(Group, Group_size, Ontology, GO_ID, Term, FG_Count, BG_Count, Odds_Ratio, p_value, p_adj_BH)
write.csv(res, file.path(OUT_DIR, "A23_GO_pooled_mid_high.csv"), row.names = FALSE)

# ---- check against the canonical results, then the pooled group ----------------------------------------------------------
can <- read.csv(file.path(OUT_DIR, "GO_enrichment_results_by_pct_uORF_plus.csv"))
for (g in names(grp)[1:3]) {
  a <- res %>% filter(Group == g, !(GO_ID %in% ROOT), p_adj_BH < FDR_CUT); b <- can %>% filter(Group == g, !(GO_ID %in% ROOT), p_adj_BH < FDR_CUT)
  .add(sprintf("  %-24s FDR<0.05 terms: here %d, canonical %d, shared %d", g, nrow(a), nrow(b), length(intersect(a$GO_ID, b$GO_ID))))
}
pooled <- res %>% filter(Group == names(grp)[4], !(GO_ID %in% ROOT)) %>% arrange(p_value)
sigp <- pooled %>% filter(p_adj_BH < FDR_CUT)
.add(sprintf("pooled medium + high (>= 30%% uORF+, n = %d): %d terms at FDR < 0.05 of %d tested", grp[[4]] %>% length(), nrow(sigp), nrow(pooled)))
for (i in seq_len(nrow(sigp))) .add(sprintf("  %s %s %s | genes %d of %d in the universe | odds ratio %.2f | p %.2g | FDR %.3f", sigp$Ontology[i], sigp$GO_ID[i], sigp$Term[i], sigp$FG_Count[i], sigp$BG_Count[i], sigp$Odds_Ratio[i], sigp$p_value[i], sigp$p_adj_BH[i]))
.add("  ten smallest p in the pooled group (whatever the FDR):")
top <- head(pooled, 10)
for (i in seq_len(nrow(top))) .add(sprintf("  %s %s %s | %d of %d | odds ratio %.2f | p %.2g | FDR %.3f", top$Ontology[i], top$GO_ID[i], top$Term[i], top$FG_Count[i], top$BG_Count[i], top$Odds_Ratio[i], top$p_value[i], top$p_adj_BH[i]))
# terms of interest (nucleus, kinase, transcription factor) in the pooled group
keys <- c("GO:0005634", "GO:0016301", "GO:0003700", "GO:0006468", "GO:0006355", "GO:0043565")
ti <- pooled %>% filter(GO_ID %in% keys)
.add("  nucleus / kinase / transcription-factor terms in the pooled group:")
for (i in seq_len(nrow(ti))) .add(sprintf("  %s %s %s | %d of %d | odds ratio %.2f | p %.2g | FDR %.3f", ti$Ontology[i], ti$GO_ID[i], ti$Term[i], ti$FG_Count[i], ti$BG_Count[i], ti$Odds_Ratio[i], ti$p_value[i], ti$p_adj_BH[i]))
writeLines(.log, file.path(OUT_DIR, "A23_summary.txt")); .add("written: output/A23_GO_pooled_mid_high.csv, output/A23_summary.txt")

# ---- grouping over ALL TuORF genes by the highest %uORF+ among a gene's TuORFs (Larry, 2026-11-05) --------------------------
# only low: every TuORF of the gene below 30% uORF+; at least one medium: highest TuORF in 30-80%; at least one high: highest >= 80%
g_all <- S1A %>% filter(is.finite(uIF)) %>% transmute(gene = toupper(gene_id), pUP = as.numeric(pct_uORF_plus), n_tu = as.integer(n_TuORFs_in_gene)) %>%
  group_by(gene) %>% summarise(max_pUP = max(pUP, na.rm = TRUE), n_tu = max(n_tu), .groups = "drop") %>% filter(is.finite(max_pUP), gene %in% universe_genes)
grpC <- list("all TuORF genes: only low %uORF+ (<30)"      = g_all$gene[g_all$max_pUP < 30],
             "all TuORF genes: at least one medium (30-80)" = g_all$gene[g_all$max_pUP >= 30 & g_all$max_pUP < 80],
             "all TuORF genes: at least one high (>=80)"   = g_all$gene[g_all$max_pUP >= 80])
.add(""); .add(sprintf("grouping over all TuORF genes in the universe (%d genes; %d with one TuORF, %d with several): %s", nrow(g_all), sum(g_all$n_tu == 1), sum(g_all$n_tu > 1),
                      paste(sprintf("%s n=%d", names(grpC), sapply(grpC, length)), collapse = "; ")))
resC <- bind_rows(lapply(names(grpC), function(g) run_fisher(grpC[[g]]) %>% mutate(Group = g))) %>% dplyr::select(Group, Group_size, Ontology, GO_ID, Term, FG_Count, BG_Count, Odds_Ratio, p_value, p_adj_BH)
write.csv(resC, file.path(OUT_DIR, "A23_GO_all_TuORF_genes_by_max_pUP.csv"), row.names = FALSE)
for (g in names(grpC)) {
  sg <- resC %>% filter(Group == g, !(GO_ID %in% ROOT)) %>% arrange(p_value); sig_g <- sg %>% filter(p_adj_BH < FDR_CUT)
  .add(sprintf("  %s: %d terms at FDR < 0.05 of %d tested", g, nrow(sig_g), nrow(sg)))
  show <- if (nrow(sig_g)) sig_g else head(sg, 5)
  for (i in seq_len(nrow(show))) .add(sprintf("     %s %s %s | %d of %d | odds ratio %.2f | p %.2g | FDR %.3f", show$Ontology[i], show$GO_ID[i], show$Term[i], show$FG_Count[i], show$BG_Count[i], show$Odds_Ratio[i], show$p_value[i], show$p_adj_BH[i]))
  tk <- sg %>% filter(GO_ID %in% keys)
  for (i in seq_len(nrow(tk))) .add(sprintf("     [key term] %s %s | %d of %d | odds ratio %.2f | p %.2g | FDR %.3f", tk$GO_ID[i], tk$Term[i], tk$FG_Count[i], tk$BG_Count[i], tk$Odds_Ratio[i], tk$p_value[i], tk$p_adj_BH[i]))
}
writeLines(.log, file.path(OUT_DIR, "A23_summary.txt")); .add("written: output/A23_GO_all_TuORF_genes_by_max_pUP.csv")
