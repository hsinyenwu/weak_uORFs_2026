################################################################################
# A8 — CPuORFs versus matched non-conserved TuORFs (Reviewer 2, §2.5)           #
#                                                                               #
# Fig 4D-F compare the 76 CPuORF genes with all other TuORF genes. The reviewer  #
# notes that CPuORFs sit in regulatory genes with their own expression          #
# properties, so high inclusion is consistent with constraint but a control set #
# of non-conserved uORFs from matched gene classes is needed. Here:             #
#  * gene classes from TAIR GO annotations propagated with GO.db                #
#    (transcription regulation, kinase/phosphorylation, signalling, transport); #
#  * control set 1 = non-CPuORF TuORF genes with the same class profile,        #
#    sampled at random (k per CPuORF gene; 1,000 draws);                         #
#  * control set 2 = the same, additionally nearest in RNA abundance, leader    #
#    length and uORF length (k nearest neighbours; 200 draws of the greedy      #
#    order);                                                                     #
#  * regression over all TuORF genes: %uORF+ >= 90, uIF and TE on CPuORF status  #
#    adjusted for class and covariates.                                          #
# One representative TuORF per gene (CPuORF first, then highest uORF density),  #
# as in Fig 4. Outputs (output/A8_*): tables, summary txt, figure pdf.           #
# Run from Weak_uORF_submission/:  Rscript revision_modules/A8_CPuORF_matched_controls.R
################################################################################
suppressPackageStartupMessages({
  library(GenomicFeatures); library(GenomicRanges); library(data.table); library(dplyr); library(ggplot2); library(openxlsx); library(gridExtra); library(grid); library(GO.db); library(AnnotationDbi)
})
options(stringsAsFactors = FALSE, warn = 1, dplyr.summarise.inform = FALSE)
set.seed(20261028)
t0 <- Sys.time()
DATA_DIR <- "Files"; OUT_DIR <- "output"
PATH_GTF <- file.path(DATA_DIR, "Araport11+CTRL_20181206.gtf"); PATH_S1 <- file.path(OUT_DIR, "Table_S1_all_TuORFs.xlsx")
PATH_CP  <- file.path(DATA_DIR, "CPuORF_collection_3_manully_fix.xlsx"); PATH_GO <- file.path(DATA_DIR, "ATH_GO_GOSLIM.tsv")
K <- 5L; N_RANDOM <- 1000L; N_NN <- 200L; PUP_CUT <- 90
for (f in c(PATH_GTF, PATH_S1, PATH_CP, PATH_GO)) stopifnot(file.exists(f))
.log <- character(0); .add <- function(...) { s <- paste0(...); message(s); .log <<- c(.log, s) }
pct <- function(x, d = 1) sprintf(paste0("%.", d, "f%%"), 100 * x); big <- function(x) trimws(format(x, big.mark = ",")); med <- function(x) median(x[is.finite(x)])
clean_peptide <- function(x) { x <- toupper(as.character(x)); x <- gsub("\\s+", "", x); x <- gsub("[\\*X]", "", x); x[nchar(x) == 0] <- NA_character_; x }

################################################################################
# 1) TuORF genes (one representative TuORF each), CPuORF status, covariates
################################################################################
.add("A8: tables")
S1A <- read.xlsx(PATH_S1, sheet = "Table S1A (CAGE)")
CP <- read.xlsx(PATH_CP, sheet = 1) %>% mutate(gene_id = trimws(gene_id), pep = clean_peptide(peptide), len = nchar(pep)) %>% filter(!is.na(gene_id))
keys_strict <- paste0(CP$gene_id, "_", CP$pep); keys_len <- paste0(CP$gene_id, "_", CP$len)
S1A <- S1A %>% mutate(pep = clean_peptide(peptide_seq), pep_m1 = ifelse(substr(pep, 1, 2) == "MM", substring(pep, 2), NA_character_),   # V23 (2026-11-08): a TuORF called from an AUG one codon upstream (extra N-terminal Met) is the same CPuORF (SAMDC AT5G15950, ATHB1 AT3G01470)
                       is_CP = paste0(trimws(gene_id), "_", pep) %in% keys_strict | paste0(trimws(gene_id), "_", nchar(pep)) %in% keys_len | (!is.na(pep_m1) & paste0(trimws(gene_id), "_", pep_m1) %in% keys_strict))
R <- S1A %>% group_by(gene_id) %>% arrange(desc(is_CP), desc(uORF_TL), desc(uORF_P_sites), .by_group = TRUE) %>% slice(1) %>% ungroup() %>%
  transmute(gene_id, tx_id, is_CP, pUP = pct_uORF_plus, uIF = ifelse(is.finite(mORF_TL) & mORF_TL > 0, uORF_TL / mORF_TL, uIF), TE = mORF_TE, RNA_TPM, uORF_len_aa = uORF_length_aa, n_TuORFs = n_TuORFs_in_gene) %>% as.data.table()
txdb <- if (requireNamespace("txdbmaker", quietly = TRUE)) txdbmaker::makeTxDbFromGFF(PATH_GTF, format = "gtf") else GenomicFeatures::makeTxDbFromGFF(PATH_GTF, format = "gtf")
utr <- fiveUTRsByTranscript(txdb, use.names = TRUE); ll <- sum(width(utr)); R[, leader_len := as.integer(ll[tx_id])]
n_all <- nrow(R); n_cp_all <- sum(R$is_CP); R <- R[is.finite(TE)]                      # the Fig 4D-F universe: TuORF genes with TE data
.add(sprintf("  TuORF genes with any detected TuORF: %d (%d CPuORF); restricted to genes with TE data as in Fig 4D-F: %d (%d CPuORF)", n_all, n_cp_all, nrow(R), sum(R$is_CP)))
.add(sprintf("  TuORF genes: %d; CPuORF genes with a detected TuORF: %d (CPuORF collection: %d entries, %d genes); genes with TE: %d; with leader length: %d", nrow(R), sum(R$is_CP), nrow(CP), length(unique(CP$gene_id)), sum(is.finite(R$TE)), sum(!is.na(R$leader_len))))
.add(sprintf("  unmatched comparison (Fig 4D-F): >= %d%% uORF+ in %s of CPuORF genes vs %s of other TuORF genes; median %%uORF+ %.1f vs %.1f; median uIF %.2f vs %.2f (%.1fx); median TE %.2f vs %.2f (%.2fx)", PUP_CUT, pct(mean(R[is_CP == TRUE, pUP] >= PUP_CUT)), pct(mean(R[is_CP == FALSE, pUP] >= PUP_CUT)), med(R[is_CP == TRUE, pUP]), med(R[is_CP == FALSE, pUP]), med(R[is_CP == TRUE, uIF]), med(R[is_CP == FALSE, uIF]), med(R[is_CP == TRUE, uIF]) / med(R[is_CP == FALSE, uIF]), med(R[is_CP == TRUE, TE]), med(R[is_CP == FALSE, TE]), med(R[is_CP == TRUE, TE]) / med(R[is_CP == FALSE, TE])))

################################################################################
# 2) Gene classes from GO (TAIR annotations, propagated with GO.db)
################################################################################
.add("A8: GO classes")
raw <- fread(PATH_GO, header = FALSE, sep = "\t", quote = ""); setnames(raw, c("gene", "go"))
go_long <- raw[, .(go_id = trimws(unlist(strsplit(go, ",")))), by = gene][grepl("^GO:", go_id)]; go_long[, gene := toupper(sub("\\.[0-9]+$", "", gene))]; go_long <- unique(go_long)
offspring <- function(id, onto) { o <- switch(onto, MF = as.list(GOMFOFFSPRING), BP = as.list(GOBPOFFSPRING))[[id]]; unique(c(id, o[!is.na(o)])) }
CLASSES <- list(
  transcription = c(offspring("GO:0003700", "MF"), offspring("GO:0006355", "BP")),          # DNA-binding TF activity; regulation of DNA-templated transcription
  kinase        = c(offspring("GO:0016301", "MF"), offspring("GO:0006468", "BP")),          # kinase activity; protein phosphorylation
  signalling    = offspring("GO:0007165", "BP"),                                            # signal transduction
  transport     = c(offspring("GO:0005215", "MF"), offspring("GO:0055085", "BP")))          # transporter activity; transmembrane transport
for (cl in names(CLASSES)) { g <- unique(go_long[go_id %in% CLASSES[[cl]], gene]); R[, (cl) := gene_id %in% g] }
R[, annotated := gene_id %in% go_long$gene]
R[, profile := paste0(ifelse(transcription, "T", ""), ifelse(kinase, "K", ""), ifelse(signalling, "S", ""), ifelse(transport, "P", ""))]; R[profile == "", profile := ifelse(annotated, "none", "unannotated")]
R[, n_reg := transcription + kinase + signalling + transport]
cls_tab <- rbindlist(lapply(c("transcription", "kinase", "signalling", "transport"), function(cl) data.table(class = cl, CPuORF = 100 * mean(R[is_CP == TRUE][[cl]]), other_TuORF = 100 * mean(R[is_CP == FALSE][[cl]]))))
cls_tab <- rbind(cls_tab, data.table(class = "any of the four", CPuORF = 100 * mean(R[is_CP == TRUE, n_reg > 0]), other_TuORF = 100 * mean(R[is_CP == FALSE, n_reg > 0])), data.table(class = "no Gene Ontology annotation", CPuORF = 100 * mean(!R[is_CP == TRUE, annotated]), other_TuORF = 100 * mean(!R[is_CP == FALSE, annotated])))
for (i in seq_len(nrow(cls_tab))) with(cls_tab[i], .add(sprintf("  class %-18s CPuORF genes %.1f%% | other TuORF genes %.1f%%", class, CPuORF, other_TuORF)))
.add("  CPuORF class profiles: ", paste(sprintf("%s %d", names(table(R[is_CP == TRUE, profile])), as.integer(table(R[is_CP == TRUE, profile]))), collapse = ", "))
# covariates (z-scores over all TuORF genes)
R[, `:=`(z_tpm = as.numeric(scale(log2(RNA_TPM + 0.1))), z_lead = as.numeric(scale(log2(pmax(leader_len, 1)))), z_len = as.numeric(scale(log2(uORF_len_aa + 1))))]
R[is.na(z_lead), z_lead := 0]

################################################################################
# 3) Control sets
################################################################################
.add("A8: control sets")
ou_genes <- unique(trimws(read.xlsx(file.path(DATA_DIR, "CiPS_ouORF_all_Final.xlsx"))$gene_id))   # V23 (2026-11-08): genes with a translated overlapping uORF are left out of the control pool, as from every control group of the study
cp <- R[is_CP == TRUE]; others_all <- R[is_CP == FALSE]; others <- others_all[!gene_id %in% ou_genes]   # V25 (2026-11-09): `others` = the control pool only; the "other TuORF genes" layer of the figure (A subtitle, B-E) is `others_all`, as the legend says (n = 4,686); in V23/V24 it was the 4,357-gene pool
.add(sprintf("  control pool: %d non-CPuORF TuORF genes with TE data after removing %d genes with a translated overlapping uORF", nrow(others), sum(R$is_CP == FALSE & R$gene_id %in% ou_genes)))
candidates <- function(g) { same <- others[profile == g$profile]
  if (nrow(same) >= 2 * K) return(list(idx = same$gene_id, level = "same class profile"))
  if (g$n_reg > 0) { share <- others[(transcription & g$transcription) | (kinase & g$kinase) | (signalling & g$signalling) | (transport & g$transport)]; if (nrow(share) >= 2 * K) return(list(idx = share$gene_id, level = "shares a regulatory class")) }
  list(idx = others$gene_id, level = "the whole control pool") }   # V26 (2026-11-09): the pool without overlapping-uORF genes, not all other TuORF genes (never reached: 74 + 2 = 76)
cand <- lapply(seq_len(nrow(cp)), function(i) candidates(cp[i])); lev <- sapply(cand, `[[`, "level")
.add("  candidate pools: ", paste(sprintf("%s %d", names(table(lev)), as.integer(table(lev))), collapse = ", "))
summ <- function(d) c(ge90 = mean(d$pUP >= PUP_CUT), med_pUP = med(d$pUP), med_uIF = med(d$uIF), med_TE = med(d$TE), med_tpm = med(d$RNA_TPM), med_lead = med(d$leader_len), med_len = med(d$uORF_len_aa))
obs <- summ(cp)
# set 1: class-matched random controls
setkey(others, gene_id)
draw_random <- function() { ids <- unlist(lapply(cand, function(cc) sample(cc$idx, K))); summ(others[J(ids)]) }
S1 <- t(replicate(N_RANDOM, draw_random()))
# set 2: class + covariate nearest neighbours, greedy without replacement in random CPuORF order
draw_nn <- function() { used <- character(0); ids <- character(0)
  for (i in sample(nrow(cp))) { pool <- setdiff(cand[[i]]$idx, used); if (length(pool) < K) pool <- setdiff(others$gene_id, used)
    o <- others[J(pool)]; d <- (o$z_tpm - cp$z_tpm[i])^2 + (o$z_lead - cp$z_lead[i])^2 + (o$z_len - cp$z_len[i])^2
    pick <- o$gene_id[order(d)][1:K]; used <- c(used, pick); ids <- c(ids, pick) }
  list(summary = summ(others[J(ids)]), ids = ids) }
nn <- replicate(N_NN, draw_nn(), simplify = FALSE); S2 <- t(sapply(nn, `[[`, "summary")); nn_ids <- nn[[1]]$ids
report <- function(S, label) { q <- apply(S, 2, quantile, c(0.025, 0.5, 0.975))
  .add(sprintf("  %s (k = %d per CPuORF gene): >= %d%% uORF+ %.1f%% (95%% range %.1f-%.1f) vs CPuORF %.1f%% [controls >= CPuORF in %.1f%% of draws]; median %%uORF+ %.1f (%.1f-%.1f) vs %.1f; median uIF %.2f (%.2f-%.2f) vs %.2f; median TE %.2f (%.2f-%.2f) vs %.2f",
    label, K, PUP_CUT, 100 * q[2, "ge90"], 100 * q[1, "ge90"], 100 * q[3, "ge90"], 100 * obs["ge90"], 100 * mean(S[, "ge90"] >= obs["ge90"]), q[2, "med_pUP"], q[1, "med_pUP"], q[3, "med_pUP"], obs["med_pUP"], q[2, "med_uIF"], q[1, "med_uIF"], q[3, "med_uIF"], obs["med_uIF"], q[2, "med_TE"], q[1, "med_TE"], q[3, "med_TE"], obs["med_TE"]))
  .add(sprintf("    covariates, controls vs CPuORF: median RNA TPM %.1f vs %.1f; leader %d vs %d nt; uORF length %d vs %d aa; draws with median uIF >= CPuORF %.1f%%, median TE <= CPuORF %.1f%%", q[2, "med_tpm"], obs["med_tpm"], as.integer(q[2, "med_lead"]), as.integer(obs["med_lead"]), as.integer(q[2, "med_len"]), as.integer(obs["med_len"]), 100 * mean(S[, "med_uIF"] >= obs["med_uIF"]), 100 * mean(S[, "med_TE"] <= obs["med_TE"])))
  q }
q1 <- report(S1, "control set 1: same gene class, random"); q2 <- report(S2, "control set 2: same gene class + nearest in RNA abundance, leader length, uORF length")
# tests on one representative nearest-neighbour set (the first draw)
ctl <- others[J(nn_ids)]
ks <- function(a, b) suppressWarnings(ks.test(a, b)$p.value)
.add(sprintf("  first nearest-neighbour set (n = %d): KS p for %%uORF+ %.1e, uIF %.1e, TE %.1e; Fisher p for >= %d%% uORF+ %.1e", nrow(ctl), ks(cp$pUP, ctl$pUP), ks(cp$uIF, ctl$uIF), ks(cp$TE[is.finite(cp$TE)], ctl$TE[is.finite(ctl$TE)]), PUP_CUT, fisher.test(table(c(rep("CP", nrow(cp)), rep("ctl", nrow(ctl))), c(cp$pUP >= PUP_CUT, ctl$pUP >= PUP_CUT)))$p.value))

################################################################################
# 4) Regression over all TuORF genes
################################################################################
.add("A8: adjusted models (all TuORF genes)")
M <- R[is.finite(uIF) & uIF > 0]
m_pup <- glm(I(pUP >= PUP_CUT) ~ is_CP + transcription + kinase + signalling + transport + z_tpm + z_lead + z_len, data = M, family = binomial)
m_uif <- lm(log2(uIF) ~ is_CP + transcription + kinase + signalling + transport + z_tpm + z_lead + z_len, data = M)
m_te  <- lm(log2(TE) ~ is_CP + transcription + kinase + signalling + transport + z_tpm + z_lead + z_len, data = M[is.finite(TE) & TE > 0])
m_pup0 <- glm(I(pUP >= PUP_CUT) ~ is_CP, data = M, family = binomial); m_uif0 <- lm(log2(uIF) ~ is_CP, data = M); m_te0 <- lm(log2(TE) ~ is_CP, data = M[is.finite(TE) & TE > 0])
ex <- function(m, f = exp) { cf <- coef(summary(m))["is_CPTRUE", ]; ci <- suppressMessages(confint.default(m))["is_CPTRUE", ]; c(est = unname(f(cf[1])), lo = unname(f(ci[1])), hi = unname(f(ci[2])), p = unname(cf[4])) }
two <- function(x) 2^x
mod_tab <- rbind(data.frame(outcome = ">= 90% uORF+ (odds ratio)", model = "unadjusted", t(ex(m_pup0))), data.frame(outcome = ">= 90% uORF+ (odds ratio)", model = "adjusted for class and covariates", t(ex(m_pup))),
                 data.frame(outcome = "uIF (fold)", model = "unadjusted", t(ex(m_uif0, two))), data.frame(outcome = "uIF (fold)", model = "adjusted for class and covariates", t(ex(m_uif, two))),
                 data.frame(outcome = "mORF TE (fold)", model = "unadjusted", t(ex(m_te0, two))), data.frame(outcome = "mORF TE (fold)", model = "adjusted for class and covariates", t(ex(m_te, two))))
for (i in seq_len(nrow(mod_tab))) with(mod_tab[i, ], .add(sprintf("  %s, %s: CPuORF effect %.2f (95%% CI %.2f-%.2f), p = %.1e", outcome, model, est, lo, hi, p)))
.add(sprintf("  (n = %d TuORF genes with uIF; %d with TE; covariates: four GO classes, RNA TPM, leader length, uORF length)", nrow(M), sum(is.finite(M$TE) & M$TE > 0)))

################################################################################
# 5) Tables
################################################################################
fwrite(R[order(-is_CP, gene_id), .(gene_id, tx_id, is_CPuORF = is_CP, pUP_CAGE = pUP, uIF, mORF_TE = TE, RNA_TPM, leader_len, uORF_len_aa, n_TuORFs, transcription, kinase, signalling, transport, GO_annotated = annotated, class_profile = profile, in_first_NN_control_set = gene_id %in% nn_ids)], file.path(OUT_DIR, "A8_Table_TuORF_genes_classes.csv"))
ctl_tab <- rbind(data.table(set = "CPuORF", metric = names(obs), value = obs, lo = NA_real_, hi = NA_real_), data.table(set = "control 1 (class, random)", metric = colnames(S1), value = q1[2, ], lo = q1[1, ], hi = q1[3, ]), data.table(set = "control 2 (class + covariates, nearest)", metric = colnames(S2), value = q2[2, ], lo = q2[1, ], hi = q2[3, ]))
fwrite(ctl_tab, file.path(OUT_DIR, "A8_Table_control_summaries.csv")); fwrite(mod_tab, file.path(OUT_DIR, "A8_Table_adjusted_models.csv")); fwrite(cls_tab, file.path(OUT_DIR, "A8_Table_class_composition.csv"))
wb <- createWorkbook()
addWorksheet(wb, "README"); writeData(wb, "README", data.frame(x = c(
  "Table S9 — CPuORFs versus matched non-conserved TuORFs (Figure S18).",
  "genes: one representative TuORF per gene (CPuORF first, then highest uORF P-site density), %uORF+ (CAGE), uIF, mORF TE, RNA TPM, annotated leader length, uORF length, classes derived from Gene Ontology terms (transcription = DNA-binding transcription factor activity or regulation of DNA-templated transcription; kinase = kinase activity or protein phosphorylation; signalling = signal transduction; transport = transporter activity or transmembrane transport; each with all GO.db descendants), class profile, membership in the first nearest-neighbour control set. GO_annotated = the gene has a Gene Ontology annotation.",
  sprintf("control_summaries: CPuORF values and the median (2.5-97.5%% range over draws) of control-set summaries; control 1 = %d random non-CPuORF TuORF genes with the same class profile per CPuORF gene (%d draws); control 2 = the %d nearest in RNA abundance, leader length and uORF length within the class pool (%d draws of the greedy order). Metrics: ge90 = fraction of genes (0 to 1) with at least 90%% uORF+; med_pUP = median %%uORF+; med_uIF, med_TE, med_tpm, med_lead and med_len = medians of uIF, mORF TE, RNA TPM, leader length (nt) and uORF length (aa).", K, N_RANDOM, K, N_NN),
  "adjusted_models: CPuORF effect on >= 90% uORF+ (odds ratio), uIF (fold) and mORF TE (fold) over all TuORF genes, unadjusted and adjusted for the four classes, RNA TPM, leader length and uORF length; est = estimate, lo and hi = limits of its 95% confidence interval.",
  "class_composition: share of CPuORF and other TuORF genes in each class. class_profile codes: T = transcription, K = kinase, S = signalling, P = transport; none = annotated but in none of the four; unannotated = no Gene Ontology annotation.",
  "CPuORFs: TuORFs whose peptide matches an entry of the CPuORF collection of the same gene (identical sequence, identical length, or identical after removing one extra N-terminal Met); control genes exclude genes with a translated overlapping uORF.")), colNames = FALSE)
addWorksheet(wb, "control_summaries"); writeData(wb, "control_summaries", ctl_tab); addWorksheet(wb, "adjusted_models"); writeData(wb, "adjusted_models", mod_tab); addWorksheet(wb, "class_composition"); writeData(wb, "class_composition", cls_tab)
addWorksheet(wb, "genes"); writeData(wb, "genes", R[order(-is_CP, gene_id), .(gene_id, tx_id, is_CPuORF = is_CP, pUP_CAGE = pUP, uIF, mORF_TE = TE, RNA_TPM, leader_len, uORF_len_aa, n_TuORFs, transcription, kinase, signalling, transport, GO_annotated = annotated, class_profile = profile, in_first_NN_control_set = gene_id %in% nn_ids)])
saveWorkbook(wb, file.path(OUT_DIR, "A8_Table_CPuORF_matched_controls.xlsx"), overwrite = TRUE)

################################################################################
# 6) Figure
################################################################################
th <- theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), legend.position = "top", legend.key.size = unit(0.35, "cm"), plot.title = element_text(size = 9.5, face = "bold"), plot.subtitle = element_text(size = 7.3))
cols <- c(CPuORF = "#E76F69", `other TuORF genes` = "#48B5C4", `matched controls (class + covariates)` = "#4D4D4D")
# A class composition
da <- melt(cls_tab, id.vars = "class", variable.name = "set", value.name = "pct"); da[, set := factor(ifelse(set == "CPuORF", "CPuORF", "other TuORF genes"), levels = c("CPuORF", "other TuORF genes"))]; da[, class := factor(class, levels = cls_tab$class)]
pA <- ggplot(da, aes(class, pct, fill = set)) + geom_col(position = position_dodge(0.8), width = 0.75) + scale_fill_manual(values = cols[1:2], name = NULL) + labs(x = NULL, y = "% of genes", title = "A  Gene classes (Gene Ontology) of CPuORF and other TuORF genes", subtitle = sprintf("CPuORF genes n = %d; other TuORF genes n = %s", nrow(cp), big(nrow(others_all)))) + th + theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 7))
# B-D ECDFs: CPuORF, all other TuORF genes, pooled nearest-neighbour controls (all draws)
pool2 <- others[J(unlist(lapply(nn, `[[`, "ids")))]
ecdf_df <- function(col) rbind(data.table(set = "CPuORF", v = cp[[col]]), data.table(set = "other TuORF genes", v = others_all[[col]]), data.table(set = "matched controls (class + covariates)", v = pool2[[col]]))[is.finite(v)][, set := factor(set, levels = names(cols))]
pE_ <- function(col, xlab, ttl, sub, logx = FALSE) { d <- ecdf_df(col); p <- ggplot(d, aes(v, colour = set)) + stat_ecdf(linewidth = 0.7) + scale_colour_manual(values = cols, name = NULL) + guides(colour = guide_legend(nrow = 3)) + labs(x = xlab, y = "Cumulative\nfraction\nof genes", title = ttl, subtitle = sub) + th + theme(legend.direction = "vertical"); if (logx) p <- p + scale_x_log10(); p }
pB <- pE_("pUP", "%uORF+ (CAGE)", "B  %uORF+", sprintf(">= %d%%: CPuORF %.0f%%, other TuORF genes %.0f%%, matched controls %.0f%% (range %.0f-%.0f)", PUP_CUT, 100 * obs["ge90"], 100 * mean(others_all$pUP >= PUP_CUT), 100 * q2[2, "ge90"], 100 * q2[1, "ge90"], 100 * q2[3, "ge90"]))
pC <- pE_("uIF", "uIF (log scale)", "C  uIF", sprintf("medians: CPuORF %.2f, other TuORF genes %.2f, matched controls %.2f (range %.2f-%.2f)", obs["med_uIF"], med(others_all$uIF), q2[2, "med_uIF"], q2[1, "med_uIF"], q2[3, "med_uIF"]), logx = TRUE)
pD <- pE_("TE", "mORF TE (log scale)", "D  mORF translation efficiency", sprintf("medians: CPuORF %.2f, other TuORF genes %.2f, matched controls %.2f (range %.2f-%.2f)", obs["med_TE"], med(others_all$TE), q2[2, "med_TE"], q2[1, "med_TE"], q2[3, "med_TE"]), logx = TRUE)
# E matching quality
dq <- rbind(cp[, .(set = "CPuORF", `RNA (TPM)` = RNA_TPM, `leader length (nt)` = as.numeric(leader_len), `uORF length (aa)` = as.numeric(uORF_len_aa))], ctl[, .(set = "matched controls (class + covariates)", `RNA (TPM)` = RNA_TPM, `leader length (nt)` = as.numeric(leader_len), `uORF length (aa)` = as.numeric(uORF_len_aa))], others_all[, .(set = "other TuORF genes", `RNA (TPM)` = RNA_TPM, `leader length (nt)` = as.numeric(leader_len), `uORF length (aa)` = as.numeric(uORF_len_aa))])
dq <- melt(dq, id.vars = "set", variable.name = "covariate", value.name = "v")[is.finite(v) & v > 0]; dq[, set := factor(set, levels = names(cols))]
pE <- ggplot(dq, aes(set, v, fill = set)) + geom_boxplot(outlier.size = 0.4, linewidth = 0.3) + facet_wrap(~covariate, scales = "free_y") + scale_y_log10() + scale_fill_manual(values = cols, guide = "none") + labs(x = NULL, y = NULL, title = "E  Matching covariates", subtitle = sprintf("first nearest-neighbour control set (n = %d) vs CPuORF genes and all other TuORF genes", nrow(ctl))) + th + theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 6.5))
# F adjusted effects
mf <- as.data.table(mod_tab); mf[, lab1 := paste0(outcome, ", ", ifelse(model == "unadjusted", "unadjusted", "adjusted"))]; mf[, label := factor(lab1, levels = rev(lab1))]   # one-line labels; "adjusted" is defined in the subtitle
pF <- ggplot(mf, aes(est, label)) + geom_vline(xintercept = 1, linetype = 3) + geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.2, orientation = "y") + geom_point(size = 2, colour = "#E76F69") + scale_x_log10() + labs(x = "CPuORF effect (odds ratio or fold; 95% confidence interval, log scale)", y = NULL, title = "F  CPuORF effect over all TuORF genes", subtitle = "adjusted = four Gene Ontology classes + RNA abundance + leader length + uORF length") + th + theme(axis.text.y = element_text(size = 6.5))
source("revision_modules/fig_compact.R")   # 6.5 x <= 9 in, Arial, >= 8 pt (2026-10-28)
# V30 (2026-11-09, Larry): the panels of a row share their plotting area (fig_grid)
w22 <- function(x) sub(" (class", "\n(class", x, fixed = TRUE)
cS22 <- c(lapply(list(pA, pB, pC, pD), compact), lapply(list(pE + scale_x_discrete(labels = w22), pF), compact_wide))   # E: the long label on two lines, so that the rotated labels take less height
cS22[[2]] <- free_xlab(cS22[[2]])   # the legends stay above the panels: the curves of B and C leave no free room inside for the three long labels
g <- fig_grid(cS22, design = "AB\nCD\nEE\nFF", heights = c(1.3, 1.3, 0.9, 0.8), top = fig_title("Figure S18. CPuORFs versus non-conserved TuORFs from matched gene classes"))   # E and F need the full width
fig_save(file.path(OUT_DIR, "A8_Figure_CPuORF_matched_controls.pdf"), g, h = 9)
.add(sprintf("Done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
writeLines(.log, file.path(OUT_DIR, "A8_summary.txt"))
