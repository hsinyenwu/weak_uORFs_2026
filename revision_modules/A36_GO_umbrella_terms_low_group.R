################################################################################
# A36 — One-word replacements for "housekeeping" in the Abstract (Larry,        #
# 2026-11-20: "Other one word terms instead of housekeeping?"). Not part of the  #
# manuscript.                                                                   #
# For each candidate word, the matching broad Gene Ontology term is propagated  #
# (a gene has the term if it is annotated with the term or any of its           #
# descendants; GO.db offspring tables) and tested in the only-low %uORF+ group  #
# against the other TuORF genes of the same universe (n = 4,886; groups as in   #
# Supplemental Figs. S15-S17; two-sided Fisher's exact test).                   #
# Output: output/A36_umbrella_terms_low_group.txt                               #
# Run from Weak_uORF_submission/:  Rscript revision_modules/A36_GO_umbrella_terms_low_group.R
################################################################################
suppressPackageStartupMessages({ library(openxlsx); library(dplyr); library(tidyr); library(GO.db); library(AnnotationDbi) })
options(stringsAsFactors = FALSE, warn = 1)
OUT_DIR <- "output"; PATH_GOSLIM <- "Files/ATH_GO_GOSLIM.tsv"
.log <- character(0); .add <- function(...) { s <- paste0(...); message(s); .log <<- c(.log, s) }
.add(sprintf("A36 broad Gene Ontology terms in the only-low %%uORF+ group — %s", format(Sys.time(), "%Y-%m-%d %H:%M")))
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
  group_by(gene) %>% summarise(max_pUP = max(pUP, na.rm = TRUE), .groups = "drop") %>% filter(is.finite(max_pUP), gene %in% universe_genes)
low <- g_all$gene[g_all$max_pUP < 30]; high <- g_all$gene[g_all$max_pUP >= 80]
stopifnot(length(universe_genes) == 4886L, length(low) == 1046L, length(high) == 2976L)   # V53 (2026-11-22): group sizes with the corrected CAGE positions (V52: 1,083 / 2,918)
gl <- go_full %>% filter(gene %in% universe_genes)
# ---- candidate words and their broad terms ---------------------------------------------------------------------------------
cand <- data.frame(word = c("metabolic", "enzyme", "plastid", "chloroplast", "photosynthesis", "cytoplasmic", "cytosolic", "secreted", "cell wall", "stress", "chaperone", "transcription regulator (contrast)", "kinase (contrast)"),
                   go = c("GO:0008152", "GO:0003824", "GO:0009536", "GO:0009507", "GO:0015979", "GO:0005737", "GO:0005829", "GO:0005576", "GO:0005618", "GO:0006950", "GO:0006457", "GO:0140110", "GO:0016301"))
offs <- function(id) {
  ont <- AnnotationDbi::select(GO.db, keys = id, columns = "ONTOLOGY", keytype = "GOID")$ONTOLOGY
  tab <- switch(ont, BP = GOBPOFFSPRING, MF = GOMFOFFSPRING, CC = GOCCOFFSPRING)
  c(id, unlist(AnnotationDbi::mget(id, tab), use.names = FALSE))
}
nl <- length(low); no <- length(universe_genes) - nl; nh <- length(high)
.add(sprintf("only low n = %d; other TuORF genes n = %d; high n = %d (universe %d)", nl, no, nh, length(universe_genes)))
.add(sprintf("%-36s %-11s %-36s %8s %8s %8s %7s %9s", "word", "GO", "term", "low %", "other %", "high %", "OR", "p"))
for (i in seq_len(nrow(cand))) {
  ids <- offs(cand$go[i]); hit <- unique(gl$gene[gl$go_id %in% ids])
  a <- sum(low %in% hit); b <- sum(setdiff(universe_genes, low) %in% hit); h <- sum(high %in% hit)
  ft <- fisher.test(matrix(c(a, nl - a, b, no - b), 2))
  term <- AnnotationDbi::select(GO.db, keys = cand$go[i], columns = "TERM", keytype = "GOID")$TERM
  .add(sprintf("%-36s %-11s %-36s %8.1f %8.1f %8.1f %7.2f %9.2g", cand$word[i], cand$go[i], substr(term, 1, 36), 100 * a / nl, 100 * b / no, 100 * h / nh, ft$estimate, ft$p.value))
}
writeLines(.log, file.path(OUT_DIR, "A36_umbrella_terms_low_group.txt")); .add("written: output/A36_umbrella_terms_low_group.txt")
