################################################################################
# A3 — uIF: abundance vs repression, and the built-in uIF–TE dependence         #
# (Reviewer 2, §2.2)                                                            #
#                                                                               #
# (1) INDEPENDENT-DATA TEST of the uIF–TE anticorrelation.  mORF TE has mORF    #
#     RPF in its numerator and uIF has it in its denominator, so part of the    #
#     inverse correlation (Fig 1D) could be shared sampling noise.  Here uIF    #
#     and TE are computed from DIFFERENT Ribo-seq libraries: replicate D1 (its  #
#     own BAM) vs the other five (merged BAM minus D1), with TE from the         #
#     per-replicate RSEM tables.  Independent halves share no reads, so any     #
#     remaining correlation is biological.  A permutation null (uORF densities  #
#     shuffled across genes) isolates the purely definitional component.        #
# (2) INCLUSION-ADJUSTED uIF: uIF divided by the fraction of transcripts that  #
#     carry the TuORF (%uORF+/100, or effective inclusion = %uORF+ x S_TIS     #
#     from A2), i.e. RPF occupancy per uORF-containing transcript.              #
# (3) Stratified checks: uIF–TE within %uORF+ strata; numerator/denominator    #
#     alone.                                                                    #
# Inputs: P-site counts per position from revision_modules/psite_extract.py    #
#   (24–28 nt reads, offsets 8–12, all primary alignments; reproduces           #
#   Files/CTRL_expressed_P_sites_sort_count with r = 0.9999), per-replicate    #
#   RSEM (CTRL_v1/data/{Ribo,RNA}/D*/), Tables S1/S2, A2 flags.                 #
# Run from Weak_uORF_submission/:  Rscript revision_modules/A3_uIF_abundance_and_independence.R
################################################################################
suppressPackageStartupMessages({
  library(GenomicFeatures); library(GenomicRanges); library(data.table); library(dplyr); library(ggplot2); library(openxlsx)
})
options(stringsAsFactors = FALSE, warn = 1)
set.seed(20261027)
t0 <- Sys.time()
DATA_DIR <- "Files"; OUT_DIR <- "output"
PATH_GTF <- file.path(DATA_DIR, "Araport11+CTRL_20181206.gtf")
PATH_PS_D1  <- file.path(DATA_DIR, "A3_psites_D1_24-28nt.tsv.gz")      # from psite_extract.py on star_D1_Aligned.sortedByCoord.out.bam (ROI)
PATH_PS_ALL <- file.path(DATA_DIR, "A3_psites_merged_24-28nt.tsv.gz")  # same on ribo_CTRL_merged.bam (ROI)
RSEM_DIR <- "/Users/wu/Desktop/CTRL_v1/data"                            # per-replicate RSEM: Ribo/D*/D*.genes.results, RNA/D*/D*.genes.results (CDS-only)
PATH_S1 <- file.path(OUT_DIR, "Table_S1_all_TuORFs.xlsx"); PATH_S2 <- file.path(OUT_DIR, "Table_S2_single_TuORF_genes.xlsx")
PATH_A2 <- file.path(OUT_DIR, "A2_single_TuORF_cohort_splicing_flags.csv")
PATH_RSEM_RNA <- file.path(DATA_DIR, "RSEM_Araport11_1st+2nd_Riboseq_RNA.genes.results_CDS_only")
PATH_RSEM_RIBO <- file.path(DATA_DIR, "RSEM_Araport11_1st+2nd_Riboseq_Ribo.genes.results_CDS_only")
PATH_5UTR <- file.path(DATA_DIR, "genes_w_5UTR.tsv"); PATH_OUORF <- file.path(DATA_DIR, "CiPS_ouORF_all_Final.xlsx")
REPS <- c("D1", "D2", "D3", "D5", "D6", "D7"); REST <- setdiff(REPS, "D1")
N_PERM <- 2000; N_BOOT <- 1000
for (f in c(PATH_GTF, PATH_PS_D1, PATH_PS_ALL, PATH_S1, PATH_S2, PATH_A2, PATH_RSEM_RNA, PATH_RSEM_RIBO, PATH_5UTR, PATH_OUORF)) stopifnot(file.exists(f))
.log <- character(0); .add <- function(...) { s <- paste0(...); message(s); .log <<- c(.log, s) }
r2 <- function(x, y) { ok <- is.finite(x) & is.finite(y); c(n = sum(ok), r = unname(cor(x[ok], y[ok])), rho = unname(cor(x[ok], y[ok], method = "spearman"))) }
fmt <- function(v) sprintf("n=%4d  Pearson r=%6.3f  Spearman rho=%6.3f", v["n"], v["r"], v["rho"])

################################################################################
# 1) Regions: TuORF blocks on the source isoform, CDS blocks of the source isoform
################################################################################
.add("A3: TxDb + regions")
txdb <- if (requireNamespace("txdbmaker", quietly = TRUE)) txdbmaker::makeTxDbFromGFF(PATH_GTF, format = "gtf") else GenomicFeatures::makeTxDbFromGFF(PATH_GTF, format = "gtf")
ex <- exonsBy(txdb, by = "tx", use.names = TRUE); cds <- cdsBy(txdb, by = "tx", use.names = TRUE)
S1A <- read.xlsx(PATH_S1, sheet = "Table S1A (CAGE)") %>% mutate(uid = paste(gene_id, uORF_genome_range, sep = "|"), chr = as.character(chr),
  g_min = as.integer(sub("^[^_]+_([0-9]+)_([0-9]+)$", "\\1", uORF_genome_range)), g_max = as.integer(sub("^[^_]+_([0-9]+)_([0-9]+)$", "\\2", uORF_genome_range)))
S2 <- read.xlsx(PATH_S2, sheet = 1) %>% mutate(uid = paste(gene_id, uORF_genome_range, sep = "|"))
span <- GRanges(S1A$chr, IRanges(S1A$g_min, S1A$g_max), strand = S1A$Strand)
ublk <- pintersect(ex[S1A$tx_id], span); ublk <- ublk[width(ublk) > 0]; names(ublk) <- S1A$uid
u_len <- sum(width(ublk))                                   # includes the stop codon (uORF_length_nt + 3)
cblk <- cds[unique(S1A$tx_id)]; c_len <- sum(width(cblk))
ub <- unlist(ublk, use.names = FALSE); ub$uid <- rep(names(ublk), elementNROWS(ublk))
cb <- unlist(cblk, use.names = FALSE); cb$tx <- rep(names(cblk), elementNROWS(cblk))

################################################################################
# 2) P-site counts per TuORF and per CDS: D1, merged (all six), rest = all - D1
################################################################################
read_ps <- function(f) { d <- fread(f, header = FALSE, col.names = c("chr", "pos", "strand", "count")); d[, chr := as.character(chr)]; GRanges(d$chr, IRanges(d$pos, d$pos), strand = d$strand, count = d$count) }
ps_D1 <- read_ps(PATH_PS_D1); ps_all <- read_ps(PATH_PS_ALL)
.add(sprintf("  P-sites in the regions of interest: D1 %d (%d positions) | all six %d (%d positions)", sum(ps_D1$count), length(ps_D1), sum(ps_all$count), length(ps_all)))
count_in <- function(ps, feat, key) { h <- findOverlaps(ps, feat, ignore.strand = FALSE); tapply(ps$count[queryHits(h)], mcols(feat)[[key]][subjectHits(h)], sum) }
u_D1 <- count_in(ps_D1, ub, "uid"); u_all <- count_in(ps_all, ub, "uid"); c_D1 <- count_in(ps_D1, cb, "tx"); c_all <- count_in(ps_all, cb, "tx")
getc <- function(v, k) { x <- as.numeric(unname(v[k])); x[is.na(x)] <- 0; x }   # tapply returns 1-D arrays; dplyr/case_when need plain vectors
U <- S1A %>% transmute(uid, gene_id, tx_id, pct_uORF_plus, uORF_P_sites_paper = uORF_P_sites, mORF_P_sites_paper = mORF_P_sites, uORF_TL_paper = uORF_TL, mORF_TL_paper = mORF_TL,
                       u_len = as.numeric(u_len[uid]), c_len = as.numeric(c_len[tx_id]),
                       u_D1 = getc(u_D1, uid), u_all = getc(u_all, uid), c_D1 = getc(c_D1, tx_id), c_all = getc(c_all, tx_id)) %>%
  mutate(u_rest = u_all - u_D1, c_rest = c_all - c_D1)
stopifnot(all(U$u_rest >= 0), all(U$c_rest >= 0))
.add(sprintf("  check vs the paper's P-site file: mORF P-sites per isoform r = %.4f (merged counts vs Table S1A mORF_P_sites); uORF P-sites r = %.3f (CiPS/RiboTaper counting rules differ slightly)",
             cor(U$c_all, U$mORF_P_sites_paper), cor(U$u_all, U$uORF_P_sites_paper)))
.add(sprintf("  share of P-sites in D1: TuORFs %.1f%% | CDSs %.1f%%", 100 * sum(U$u_D1) / sum(U$u_all), 100 * sum(U$c_D1) / sum(U$c_all)))
dens <- function(cnt, len) ifelse(cnt > 0 & len > 0, cnt / len, NA_real_)
U <- U %>% mutate(uIF_paper = uORF_TL_paper / mORF_TL_paper,
                  uIF_all  = dens(u_all, u_len) / dens(c_all, c_len),
                  uIF_D1   = dens(u_D1, u_len) / dens(c_D1, c_len),
                  uIF_rest = dens(u_rest, u_len) / dens(c_rest, c_len))

################################################################################
# 3) TE per replicate (RSEM, CDS-only) and the paper's TE
################################################################################
rd <- function(kind, r) { d <- fread(sprintf("%s/%s/%s/%s.genes.results", RSEM_DIR, kind, r, r), select = c("gene_id", "TPM")); setnames(d, c("gene_id", paste0(kind, "_", r))); d }
rsem <- Reduce(function(a, b) merge(a, b, by = "gene_id", all = TRUE), c(lapply(REPS, function(r) rd("Ribo", r)), lapply(REPS, function(r) rd("RNA", r))))
paperRNA <- fread(PATH_RSEM_RNA, select = c("gene_id", "TPM")); paperRibo <- fread(PATH_RSEM_RIBO, select = c("gene_id", "TPM"))
utr5_genes <- read.delim(PATH_5UTR)$gene_id
TE_tbl <- data.frame(gene_id = paperRNA$gene_id, RNA_paper = paperRNA$TPM, Ribo_paper = paperRibo$TPM[match(paperRNA$gene_id, paperRibo$gene_id)]) %>%
  mutate(TE_paper = Ribo_paper / RNA_paper) %>% filter(gene_id %in% utr5_genes, RNA_paper > 1, is.finite(TE_paper))      # same universe as the main script
mean_cols <- function(d, cols) rowMeans(as.matrix(d[, cols, with = FALSE]), na.rm = FALSE)
rsem[, TE_D1   := ifelse(Ribo_D1 > 0 & RNA_D1 > 0, Ribo_D1 / RNA_D1, NA_real_)]
rsem[, TE_rest := { rb <- mean_cols(rsem, paste0("Ribo_", REST)); rn <- mean_cols(rsem, paste0("RNA_", REST)); ifelse(rb > 0 & rn > 0, rb / rn, NA_real_) }]
rsem[, TE_all6 := { rb <- mean_cols(rsem, paste0("Ribo_", REPS)); rn <- mean_cols(rsem, paste0("RNA_", REPS)); ifelse(rb > 0 & rn > 0, rb / rn, NA_real_) }]
rsem[, TE_D123 := { rb <- mean_cols(rsem, paste0("Ribo_", c("D1","D2","D3"))); rn <- mean_cols(rsem, paste0("RNA_", c("D1","D2","D3"))); ifelse(rb > 0 & rn > 0, rb / rn, NA_real_) }]
rsem[, TE_D567 := { rb <- mean_cols(rsem, paste0("Ribo_", c("D5","D6","D7"))); rn <- mean_cols(rsem, paste0("RNA_", c("D5","D6","D7"))); ifelse(rb > 0 & rn > 0, rb / rn, NA_real_) }]
TE_tbl <- TE_tbl %>% left_join(as.data.frame(rsem)[, c("gene_id", "TE_D1", "TE_rest", "TE_all6", "TE_D123", "TE_D567")], by = "gene_id")
.add(sprintf("  TE universe (5'UTR genes, RNA TPM > 1): %d genes; per-replicate RSEM available for %d; r(log2 TE_all6, log2 TE_paper) = %.4f; replicate reliability r(log2 TE D1, log2 TE rest) = %.3f, r(D1-3, D5-7) = %.3f",
             nrow(TE_tbl), sum(!is.na(TE_tbl$TE_all6)), r2(log2(TE_tbl$TE_all6), log2(TE_tbl$TE_paper))["r"], r2(log2(TE_tbl$TE_D1), log2(TE_tbl$TE_rest))["r"], r2(log2(TE_tbl$TE_D123), log2(TE_tbl$TE_D567))["r"]))

ou <- read.xlsx(PATH_OUORF); ou_genes <- unique(if ("gene_id" %in% names(ou)) as.character(ou$gene_id) else substr(as.character(ou$tx_id), 1, 9))
upstream_translated <- unique(c(S1A$gene_id, ou_genes))
noTu <- TE_tbl %>% filter(!(gene_id %in% upstream_translated))
.add(sprintf("  No-TuORF control genes (no TuORF, no ouORF): %d", nrow(noTu)))

################################################################################
# 4) Single-TuORF cohort table
################################################################################
A2 <- fread(PATH_A2) %>% mutate(uid = paste(gene_id, uORF_genome_range, sep = "|")) %>% select(uid, pUP_eff_CAGE, proxy_status, S_TIS)
C <- S2 %>% select(uid, gene_id, pct_uORF_plus, uIF_Group, uORF_length_aa) %>%
  left_join(U %>% select(-gene_id, -pct_uORF_plus), by = "uid") %>% left_join(TE_tbl, by = "gene_id") %>% left_join(A2, by = "uid") %>%
  mutate(uIF_incl = ifelse(pct_uORF_plus >= 1, uIF_paper / (pct_uORF_plus / 100), NA_real_), uIF_incl_eff = ifelse(!is.na(pUP_eff_CAGE) & pUP_eff_CAGE >= 1, uIF_paper / (pUP_eff_CAGE / 100), NA_real_),
         pUP_grp = cut(pct_uORF_plus, c(-Inf, 30, 80, Inf), right = FALSE, labels = c("<30", "30-80", "80-100")))
stopifnot(nrow(C) == nrow(S2))
.add(sprintf("  cohort: %d genes; uIF_paper reproduced from the merged BAM: r(log2 uIF_all, log2 uIF_paper) = %.3f; uIF available from D1 alone for %d genes (uORF or CDS without a D1 P-site otherwise), from the other five for %d",
             nrow(C), r2(log2(C$uIF_all), log2(C$uIF_paper))["r"], sum(is.finite(C$uIF_D1)), sum(is.finite(C$uIF_rest))))

################################################################################
# 5) Correlations: paper, same-data halves, independent halves, reliability
################################################################################
.add(""); .add("=== 5. uIF vs mORF TE in the single-TuORF cohort (log2 both; Pearson r as in Fig 1D, Spearman on raw values) ===")
rows <- list(
  c("Fig 1D as published: uIF (all six libraries) vs TE (all six; RSEM)", "uIF_paper", "TE_paper"),
  c("reconstructed from the merged BAM: uIF_all vs TE_all6", "uIF_all", "TE_all6"),
  c("SAME libraries, D1 only: uIF_D1 vs TE_D1", "uIF_D1", "TE_D1"),
  c("SAME libraries, other five: uIF_rest vs TE_rest", "uIF_rest", "TE_rest"),
  c("INDEPENDENT: uIF from D1 vs TE from the other five", "uIF_D1", "TE_rest"),
  c("INDEPENDENT: uIF from the other five vs TE from D1", "uIF_rest", "TE_D1"),
  c("reliability: uIF_D1 vs uIF_rest", "uIF_D1", "uIF_rest"),
  c("reliability: TE_D1 vs TE_rest", "TE_D1", "TE_rest"),
  c("numerator only: uORF P-site density (all) vs TE_paper", "uORF_TL_paper", "TE_paper"),
  c("denominator only: mORF P-site density (all) vs TE_paper (shared term)", "mORF_TL_paper", "TE_paper"))
cor_tab <- do.call(rbind, lapply(rows, function(r) { v <- r2(log2(C[[r[2]]]), log2(C[[r[3]]])); data.frame(comparison = r[1], x = r[2], y = r[3], n = v["n"], pearson_r_log2 = round(v["r"], 3), spearman_rho = round(v["rho"], 3)) }))
for (k in seq_len(nrow(cor_tab))) .add(sprintf("  %-72s %s", cor_tab$comparison[k], fmt(c(n = cor_tab$n[k], r = cor_tab$pearson_r_log2[k], rho = cor_tab$spearman_rho[k]))))
# bootstrap CIs for the two independent correlations and the published one
boot_ci <- function(x, y) { ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; b <- replicate(N_BOOT, { i <- sample.int(length(x), replace = TRUE); cor(x[i], y[i]) }); quantile(b, c(0.025, 0.975)) }
for (p in list(c("uIF_paper", "TE_paper"), c("uIF_D1", "TE_rest"), c("uIF_rest", "TE_D1"))) { ci <- boot_ci(log2(C[[p[1]]]), log2(C[[p[2]]])); .add(sprintf("  95%% bootstrap CI, %s vs %s: [%.3f, %.3f]", p[1], p[2], ci[1], ci[2])) }
# same genes for a fair comparison (those with all four values)
cc <- C %>% filter(is.finite(uIF_D1), is.finite(uIF_rest), is.finite(TE_D1), is.finite(TE_rest))
.add(sprintf("  on the %d genes with all four values: same-data r = %.3f (D1) / %.3f (rest); independent r = %.3f (uIF_D1~TE_rest) / %.3f (uIF_rest~TE_D1); published pairs on the same genes r = %.3f",
             nrow(cc), r2(log2(cc$uIF_D1), log2(cc$TE_D1))["r"], r2(log2(cc$uIF_rest), log2(cc$TE_rest))["r"], r2(log2(cc$uIF_D1), log2(cc$TE_rest))["r"], r2(log2(cc$uIF_rest), log2(cc$TE_D1))["r"], r2(log2(cc$uIF_paper), log2(cc$TE_paper))["r"]))

################################################################################
# 6) Permutation null: how much anticorrelation does the shared mORF term alone produce?
################################################################################
.add(""); .add(sprintf("=== 6. Permutation null (%d shuffles of the uORF density across genes; mORF density and TE kept per gene) ===", N_PERM))
perm_r <- function(uTL, mTL, TE) { ok <- is.finite(uTL) & is.finite(mTL) & is.finite(TE) & uTL > 0 & mTL > 0 & TE > 0; uTL <- uTL[ok]; mTL <- mTL[ok]; TE <- TE[ok]
  obs <- cor(log2(uTL / mTL), log2(TE)); nul <- replicate(N_PERM, cor(log2(sample(uTL) / mTL), log2(TE))); c(n = sum(ok), observed = obs, null_mean = mean(nul), null_lo = quantile(nul, 0.025), null_hi = quantile(nul, 0.975), p_two_sided = mean(abs(nul - mean(nul)) >= abs(obs - mean(nul)))) }
pn1 <- perm_r(C$uORF_TL_paper, C$mORF_TL_paper, C$TE_paper)
.add(sprintf("  published design (all libraries): observed r = %.3f | null from the shared mORF term alone: mean %.3f [%.3f, %.3f] | n = %d | p(perm) = %.4f", pn1["observed"], pn1["null_mean"], pn1["null_lo.2.5%"], pn1["null_hi.97.5%"], pn1["n"], pn1["p_two_sided"]))
.add(sprintf("  => biological component = observed - null = %.3f of the %.3f", pn1["observed"] - pn1["null_mean"], pn1["observed"]))
pn2 <- perm_r(dens(C$u_D1, C$u_len), dens(C$c_D1, C$c_len), C$TE_rest)
.add(sprintf("  independent design (uIF from D1, TE from the other five): observed r = %.3f | null mean %.3f [%.3f, %.3f] | n = %d", pn2["observed"], pn2["null_mean"], pn2["null_lo.2.5%"], pn2["null_hi.97.5%"], pn2["n"]))
pn3 <- perm_r(dens(C$u_rest, C$u_len), dens(C$c_rest, C$c_len), C$TE_D1)
.add(sprintf("  independent design (uIF from the other five, TE from D1): observed r = %.3f | null mean %.3f [%.3f, %.3f] | n = %d", pn3["observed"], pn3["null_mean"], pn3["null_lo.2.5%"], pn3["null_hi.97.5%"], pn3["n"]))

# Partial correlations: association between uIF (equivalently uORF density) and TE at FIXED mORF density = the uORF-specific component
pc <- function(x, y, z) { ok <- is.finite(x) & is.finite(y) & is.finite(z); rx <- resid(lm(x[ok] ~ z[ok])); ry <- resid(lm(y[ok] ~ z[ok])); c(n = sum(ok), r = cor(rx, ry)) }
p1 <- pc(log2(C$uIF_paper), log2(C$TE_paper), log2(C$mORF_TL_paper))
p2 <- pc(log2(C$uIF_D1), log2(C$TE_rest), log2(dens(C$c_D1, C$c_len)))
p3 <- pc(log2(C$uIF_rest), log2(C$TE_D1), log2(dens(C$c_rest, C$c_len)))
.add(sprintf("  partial correlation of uIF with TE given the gene's mORF P-site density (shared term held fixed): published %.3f (n=%d) | independent %.3f (n=%d) / %.3f (n=%d)", p1["r"], p1["n"], p2["r"], p2["n"], p3["r"], p3["n"]))
.add(sprintf("  coupling of the two densities: r(log2 uORF density, log2 mORF density) = %.3f; r(log2 mORF density, log2 RNA TPM) = %.3f", r2(log2(C$uORF_TL_paper), log2(C$mORF_TL_paper))["r"], r2(log2(C$mORF_TL_paper), log2(C$RNA_paper))["r"]))
uTE <- C$uORF_TL_paper / C$RNA_paper   # uORF P-site density per RNA (shares the RNA denominator with mORF TE -> positive bias)
.add(sprintf("  uORF translation per RNA (uORF density / RNA TPM) vs mORF TE (shared RNA denominator, bias towards +): %s", fmt(r2(log2(uTE), log2(C$TE_paper)))))

################################################################################
# 7) uIF groups -> TE of the other libraries (Fig 1E with independent data)
################################################################################
grp3 <- function(u) factor(case_when(u > 1 ~ "High (uIF > 1)", u > 0.2 ~ "Medium (0.2 < uIF <= 1)", TRUE ~ "Low (uIF <= 0.2)"), levels = c("High (uIF > 1)", "Medium (0.2 < uIF <= 1)", "Low (uIF <= 0.2)"))
ecdf_df <- function(ucol, tecol, label) bind_rows(C %>% filter(is.finite(.data[[ucol]]), is.finite(.data[[tecol]])) %>% transmute(design = label, Group = as.character(grp3(.data[[ucol]])), TE = .data[[tecol]]),
                                                    noTu %>% filter(is.finite(.data[[tecol]])) %>% transmute(design = label, Group = "No TuORF", TE = .data[[tecol]]))
E1 <- ecdf_df("uIF_paper", "TE_paper", "published: uIF and TE from all six libraries")
E2 <- ecdf_df("uIF_D1", "TE_rest", "independent: uIF from D1, TE from the other five")
E3 <- ecdf_df("uIF_rest", "TE_D1", "independent: uIF from the other five, TE from D1")
.add(""); .add("=== 7. Median mORF TE by uIF group (High / Medium / Low / No TuORF) ===")
for (E in list(E1, E2, E3)) { m <- E %>% group_by(Group) %>% summarise(n = n(), med = median(TE), .groups = "drop"); m <- m[match(c("High (uIF > 1)", "Medium (0.2 < uIF <= 1)", "Low (uIF <= 0.2)", "No TuORF"), m$Group), ]
  pw <- pairwise.wilcox.test(E$TE, E$Group, p.adjust.method = "BH", exact = FALSE)$p.value
  .add(sprintf("  %-55s %s | BH-adjusted Wilcoxon High vs Low p = %.2g, Low vs No TuORF p = %.2g", E$design[1], paste(sprintf("%s: %.2f (n=%d)", c("High", "Medium", "Low", "No TuORF"), m$med, m$n), collapse = " | "),
               pw["Low (uIF <= 0.2)", "High (uIF > 1)"], tryCatch(pw["Low (uIF <= 0.2)", "No TuORF"], error = function(e) pw["No TuORF", "Low (uIF <= 0.2)"]))) }

################################################################################
# 8) Inclusion-adjusted uIF (occupancy per uORF-containing transcript)
################################################################################
.add(""); .add("=== 8. Inclusion-adjusted uIF: uIF / (%uORF+/100) [and / (effective inclusion/100) from A2] ===")
d8 <- C %>% filter(is.finite(TE_paper), TE_paper > 0, is.finite(uIF_paper), uIF_paper > 0)
.add(sprintf("  uIF vs TE          : %s", fmt(r2(log2(d8$uIF_paper), log2(d8$TE_paper)))))
.add(sprintf("  uIF_incl vs TE     : %s", fmt(r2(log2(d8$uIF_incl), log2(d8$TE_paper)))))
d8e <- d8 %>% filter(is.finite(uIF_incl_eff))
.add(sprintf("  on genes with effective inclusion (>= 10 RNA-seq reads at the TIS): uIF vs TE %s | uIF_incl vs TE %s | uIF_incl_eff vs TE %s",
             fmt(r2(log2(d8e$uIF_paper), log2(d8e$TE_paper))), fmt(r2(log2(d8e$uIF_incl), log2(d8e$TE_paper))), fmt(r2(log2(d8e$uIF_incl_eff), log2(d8e$TE_paper)))))
.add(sprintf("  %%uORF+ vs uIF: r(pUP, log2 uIF) = %.3f | r(pUP, log2 uIF_incl) = %.3f (division by pUP removes the abundance share; a residual positive value means uORFs on fully-included transcripts are also more occupied per transcript)",
             r2(C$pct_uORF_plus, log2(C$uIF_paper))["r"], r2(C$pct_uORF_plus, log2(C$uIF_incl))["r"]))
okx <- is.finite(C$uIF_incl); xt <- table(uIF = grp3(C$uIF_paper[okx]), uIF_incl = grp3(C$uIF_incl[okx]))   # V23: genes without a defined value (%uORF+ < 1) left out of the table (they were counted as Low before)
.add("  group membership, rows = uIF groups, columns = uIF_incl groups (High / Medium / Low):")
for (r in rownames(xt)) .add(sprintf("    %-26s %s", r, paste(xt[r, ], collapse = " / ")))
.add(sprintf("  genes changing group: %d of %d with a defined value (%.1f%%); all moves are upward since uIF_incl >= uIF; Low -> High: %d (%.1f%%)", sum(xt) - sum(diag(xt)), sum(xt), 100 * (sum(xt) - sum(diag(xt))) / sum(xt), xt[3, 1], 100 * xt[3, 1] / sum(xt)))
d7 <- C %>% filter(is.finite(uIF_incl), is.finite(uIF_paper), uIF_paper > 0)
.add(sprintf("  on the same %d genes (%%uORF+ >= 1): r(%%uORF+, log2 uIF) = %.3f | r(%%uORF+, log2 uIF_incl) = %.3f (Supplemental Fig. S4A-B)", nrow(d7), r2(d7$pct_uORF_plus, log2(d7$uIF_paper))["r"], r2(d7$pct_uORF_plus, log2(d7$uIF_incl))["r"]))
E4 <- bind_rows(C %>% filter(is.finite(uIF_incl), is.finite(TE_paper), TE_paper > 0) %>% transmute(design = "inclusion-adjusted uIF groups, TE from all six", Group = as.character(grp3(uIF_incl)), TE = TE_paper),
                noTu %>% transmute(design = "inclusion-adjusted uIF groups, TE from all six", Group = "No TuORF", TE = TE_paper))
m <- E4 %>% group_by(Group) %>% summarise(n = n(), med = median(TE), .groups = "drop"); m <- m[match(c("High (uIF > 1)", "Medium (0.2 < uIF <= 1)", "Low (uIF <= 0.2)", "No TuORF"), m$Group), ]
.add(sprintf("  median TE by uIF_incl group: %s", paste(sprintf("%s: %.2f (n=%d)", c("High", "Medium", "Low", "No TuORF"), m$med, m$n), collapse = " | ")))
# V31 (2026-11-10, Larry: 'add the table to Figure S7'; the genes with at least 90% uORF+): how much of the differences in uIF among genes %uORF+ accounts for.
# d7 = the 2,856 genes of Supplemental Fig. S5A-B (%uORF+ >= 1, uIF > 0). (1) r^2 of panel A; (2) medians by %uORF+ group (panel D); (3) the range from the 5th to the 95th percentile of genes (panel E).
r_AB <- r2(d7$pct_uORF_plus, log2(d7$uIF_paper))["r"]
.add(sprintf("  V31: share of the variance of log uIF accounted for by %%uORF+ (r^2 of Supplemental Fig. S4A): %.1f%% (r = %.3f); variance of the log inclusion-adjusted uIF = %.0f%% of the variance of log uIF",
             100 * r_AB^2, r_AB, 100 * var(log2(d7$uIF_incl)) / var(log2(d7$uIF_paper))))
BIN_BRK <- c(1, 10, 30, 50, 80, 90, Inf); BIN_LAB <- c("1-10%", "10-30%", "30-50%", "50-80%", "80-90%", "90-100%")   # each group includes its lower bound; the last one includes 100%
d7$pUP_bin <- cut(d7$pct_uORF_plus, BIN_BRK, right = FALSE, labels = BIN_LAB)
tabD <- d7 %>% group_by(pUP_bin) %>% summarise(n = n(), median_uIF = median(uIF_paper), median_uIF_incl = median(uIF_incl), .groups = "drop") %>% as.data.frame()
spread <- function(x, label) data.frame(measure = label, n = length(x), p05 = unname(quantile(x, 0.05)), median = median(x), p95 = unname(quantile(x, 0.95)), ratio_p95_p05 = unname(quantile(x, 0.95) / quantile(x, 0.05)))
tabE <- rbind(spread(d7$uIF_paper, "uIF, all genes"), spread(d7$uIF_incl, "inclusion-adjusted uIF, all genes"), spread(d7$uIF_paper[d7$pct_uORF_plus >= 90], "uIF, genes with at least 90% uORF+"))
.add("  V31: median uIF and median inclusion-adjusted uIF by %uORF+ group (Supplemental Fig. S4D):")
for (k in seq_len(nrow(tabD))) .add(sprintf("    %-8s n=%4d  median uIF %.2f  median inclusion-adjusted uIF %.2f", tabD$pUP_bin[k], tabD$n[k], tabD$median_uIF[k], tabD$median_uIF_incl[k]))
.add("  V31: 5th and 95th percentile of genes (Supplemental Fig. S4E):")
for (k in seq_len(nrow(tabE))) .add(sprintf("    %-36s n=%4d  5th %.3g  median %.3g  95th %.3g  95th/5th = %.0f-fold", tabE$measure[k], tabE$n[k], tabE$p05[k], tabE$median[k], tabE$p95[k], tabE$ratio_p95_p05[k]))
d7b <- d7[d7$pct_uORF_plus >= 10, ]   # without the 1-10% group, whose median inclusion-adjusted uIF is the highest
.add(sprintf("  V31: ratio of the median uIF, 90-100%% group / 1-10%% group = %.1f; r(%%uORF+, log uIF_incl) without the 1-10%% group = %.3f (n = %d; all genes %.3f)",
             tabD$median_uIF[6] / tabD$median_uIF[1], r2(d7b$pct_uORF_plus, log2(d7b$uIF_incl))["r"], nrow(d7b), r2(d7$pct_uORF_plus, log2(d7$uIF_incl))["r"]))
fwrite(tabD, file.path(OUT_DIR, "A3_uIF_by_inclusion_group.csv")); fwrite(tabE, file.path(OUT_DIR, "A3_uIF_range_among_genes.csv"))

################################################################################
# 9) Within %uORF+ strata: does uIF still track TE when inclusion is held (nearly) constant?
################################################################################
.add(""); .add("=== 9. uIF vs TE within %uORF+ strata (published values; and independent halves) ===")
for (g in levels(C$pUP_grp)) { d <- C %>% filter(pUP_grp == g)
  .add(sprintf("  %%uORF+ %-6s published %s | independent (uIF_D1 ~ TE_rest) %s | (uIF_rest ~ TE_D1) %s", g, fmt(r2(log2(d$uIF_paper), log2(d$TE_paper))), fmt(r2(log2(d$uIF_D1), log2(d$TE_rest))), fmt(r2(log2(d$uIF_rest), log2(d$TE_D1))))) }
d90 <- C %>% filter(pct_uORF_plus >= 90)
.add(sprintf("  %%uORF+ >= 90 (inclusion essentially complete, n=%d): published r = %.3f | independent r = %.3f / %.3f", nrow(d90), r2(log2(d90$uIF_paper), log2(d90$TE_paper))["r"], r2(log2(d90$uIF_D1), log2(d90$TE_rest))["r"], r2(log2(d90$uIF_rest), log2(d90$TE_D1))["r"]))

################################################################################
# 10) Outputs
################################################################################
fwrite(cor_tab, file.path(OUT_DIR, "A3_correlation_table.csv"))
fwrite(C %>% select(gene_id, uORF_genome_range = uid, pct_uORF_plus, pUP_eff_CAGE, proxy_status, uIF_paper, uIF_all, uIF_D1, uIF_rest, uIF_incl, uIF_incl_eff, u_D1, u_rest, c_D1, c_rest, u_len, c_len,
                    TE_paper, TE_D1, TE_rest, TE_all6), file.path(OUT_DIR, "A3_single_TuORF_cohort_uIF_variants.csv"))
perm_tab <- rbind(published = pn1, independent_uIF_D1_TE_rest = pn2, independent_uIF_rest_TE_D1 = pn3)
write.csv(perm_tab, file.path(OUT_DIR, "A3_permutation_null.csv"))

pal4 <- c("High (uIF > 1)" = "#F8766D", "Medium (0.2 < uIF <= 1)" = "#00BA38", "Low (uIF <= 0.2)" = "#619CFF", "No TuORF" = "#999999")
sc <- function(d, xcol, ycol, title, xl, yl) { d <- d %>% filter(is.finite(.data[[xcol]]), is.finite(.data[[ycol]]), .data[[xcol]] > 0, .data[[ycol]] > 0); v <- r2(log2(d[[xcol]]), log2(d[[ycol]]))
  ggplot(d, aes(x = log2(.data[[xcol]]), y = log2(.data[[ycol]]))) + geom_point(colour = "#3498db", alpha = 0.5, size = 1.2) + geom_smooth(method = "lm", colour = "red", se = FALSE, linewidth = 1) +
    theme_classic() + labs(title = title, subtitle = sprintf("r = %.2f, rho = %.2f, n = %d", v["r"], v["rho"], v["n"]), x = xl, y = yl) + theme(plot.title = element_text(face = "bold", size = 14)) }
ec <- function(E, title) { lv <- c("High (uIF > 1)", "Medium (0.2 < uIF <= 1)", "Low (uIF <= 0.2)", "No TuORF"); m <- E %>% group_by(Group) %>% summarise(n = n(), .groups = "drop"); E$Group <- factor(E$Group, levels = lv)
  labs_n <- setNames(sprintf("%s (n=%d)", lv, m$n[match(lv, m$Group)]), lv)
  ggplot(E %>% filter(TE <= 3.5), aes(x = TE, colour = Group)) + stat_ecdf(geom = "step", linewidth = 1.1) + coord_cartesian(xlim = c(0, 3.5)) + scale_colour_manual(values = pal4, labels = labs_n, name = NULL) +
    theme_classic() + labs(title = title, x = "mORF TE", y = "Fraction") + theme(plot.title = element_text(face = "bold", size = 14), legend.position = c(0.62, 0.35), legend.background = element_blank(), legend.text = element_text(size = 8)) }
pA <- sc(C, "uIF_paper", "TE_paper", "A  published (shared libraries)", bquote(log[2]~"uIF (all six libraries)"), bquote(log[2]~"mORF TE (all six)"))
pB <- sc(C, "uIF_D1", "TE_rest", "B  independent libraries", bquote(log[2]~"uIF (D1 only)"), bquote(atop(log[2]~"mORF TE", "(D2, D3, D5, D6, D7)")))
pC <- sc(C, "uIF_rest", "TE_D1", "C  independent libraries, reversed", bquote(log[2]~"uIF (D2, D3, D5, D6, D7)"), bquote(log[2]~"mORF TE (D1 only)"))
nul <- { ok <- is.finite(C$uORF_TL_paper) & is.finite(C$mORF_TL_paper) & is.finite(C$TE_paper) & C$uORF_TL_paper > 0 & C$TE_paper > 0; u <- C$uORF_TL_paper[ok]; m0 <- C$mORF_TL_paper[ok]; te <- C$TE_paper[ok]; replicate(N_PERM, cor(log2(sample(u) / m0), log2(te))) }
pD <- ggplot(data.frame(r = nul), aes(x = r)) + geom_histogram(bins = 60, fill = "#999999", colour = "white") + geom_vline(xintercept = pn1["observed"], colour = "red", linewidth = 1) +
  theme_classic() +
  labs(title = "D  built-in component (uORF density shuffled)", subtitle = sprintf("observed r = %.2f (red line); shared-term null: mean r = %.2f", pn1["observed"], pn1["null_mean"]), x = bquote("r("*log[2]~"uIF,"~log[2]~"TE) under the null"), y = "Permutations") + theme(plot.title = element_text(face = "bold", size = 14))
pE <- ec(E2, "E  TE of the other five libraries by D1 uIF group")
pF <- ec(E3, "F  TE of D1 by uIF group of the other five")
pG <- sc(C, "uIF_incl", "TE_paper", "G  inclusion-adjusted uIF", bquote(log[2]~"uIF / (%uORF+/100)"), bquote(log[2]~"mORF TE"))
pH <- ec(E4, "H  TE by inclusion-adjusted uIF group")
source("revision_modules/fig_compact.R")   # 6.5 x <= 9 in, Arial, >= 8 pt (2026-10-28)
fig <- ggpubr::annotate_figure(ggpubr::ggarrange(plotlist = lapply(list(pA, pB, pC, pD, pE, pF, pG, pH), compact, geom_scale = 0.5), ncol = 2, nrow = 4), top = fig_title("Not a manuscript figure (module A3). The uIF-TE relationship with independent libraries, its built-in component, and inclusion-adjusted uIF"))   # V40 (2026-11-11): this page had the title "Figure S7"; the manuscript figure of this module is A3_Figure_uIF_inclusion_adjusted.pdf (Supplemental Fig. S5)
fig_save(file.path(OUT_DIR, "A3_Figure_uIF_independence_and_inclusion.pdf"), fig, h = 9)
# V15 (2026-11-04, Larry): Supplemental Fig. S3 shows the inclusion-adjusted uIF only (the former panels G-H, relabelled A-B)
# V23 (2026-11-08): Supplemental Fig. S7 shows the adjustment against %uORF+ (its purpose: the inclusion share) and the TE pattern of the classes
scp <- function(d, ycol, title, yl) { d <- d %>% filter(is.finite(uIF_incl), is.finite(.data[[ycol]]), .data[[ycol]] > 0); v <- r2(d$pct_uORF_plus, log2(d[[ycol]]))
  ggplot(d, aes(x = pct_uORF_plus, y = log2(.data[[ycol]]))) + geom_point(colour = "#3498db", alpha = 0.4, size = 1.1) + geom_smooth(method = "lm", colour = "red", se = FALSE, linewidth = 1) +
    theme_classic() + labs(title = title, subtitle = sprintf("r = %.2f, rho = %.2f, n = %d", v["r"], v["rho"], v["n"]), x = "%uORF+ (CAGE-seq)", y = yl) + theme(plot.title = element_text(face = "bold", size = 14)) }
pS7A <- scp(C, "uIF_paper", "A  uIF versus %uORF+", bquote(log[2]~"uIF")); pS7B <- scp(C, "uIF_incl", "B  inclusion-adjusted uIF versus %uORF+", bquote(log[2]~"uIF / (%uORF+/100)"))
pS7C <- pH + labs(title = "C  TE by inclusion-adjusted uIF group")
# V30 (2026-11-09, Larry): titles wrapped for the three narrow panels (they ran off the panel); the legend of C at the right of its panel (no free room inside); aligned plotting areas (fig_grid)
cS7 <- list(compact(pS7A, geom_scale = 0.5, wrap_title = 22, wrap_sub = 32), compact(pS7B, geom_scale = 0.5, wrap_title = 22, wrap_sub = 32), legend_side(compact(pS7C, geom_scale = 0.5, wrap_title = 30), wrap_paren))
# V31 (2026-11-10, Larry): a second row with two tables. D = median uIF and median inclusion-adjusted uIF by %uORF+ group; E = the 5th and 95th percentile of genes
# (uIF and inclusion-adjusted uIF of all genes; uIF of the genes with at least 90% uORF+). Text 8 pt, titles 10 pt bold as in the panels above.
sig3 <- function(x) sub("\\.$", "", formatC(signif(x, 3), digits = 3, format = "fg", flag = "#"))                                  # three significant digits, so that the ratio can be recomputed from the table
tbl_panel <- function(df, title, widths) {
  th <- gridExtra::ttheme_minimal(base_size = FIG_PT_MIN, base_family = FIG_FONT, core = list(fg_params = list(hjust = 0, x = 0.06, lineheight = 0.9), padding = unit(c(2, 2.6), "mm")),
                                  colhead = list(fg_params = list(hjust = 0, x = 0.06, fontface = "bold", lineheight = 0.9), padding = unit(c(2, 2.2), "mm")))
  tg <- gridExtra::tableGrob(df, rows = NULL, theme = th); tg$widths <- unit(widths, "in")
  rule <- function() segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"), y0 = unit(0, "npc"), y1 = unit(0, "npc"), gp = gpar(lwd = 0.8))
  tg <- gtable::gtable_add_grob(tg, rule(), t = 1, b = 1, l = 1, r = ncol(tg)); tg <- gtable::gtable_add_grob(tg, rule(), t = nrow(tg), b = nrow(tg), l = 1, r = ncol(tg))   # rules under the header and under the last row
  tg$vp <- viewport(x = unit(8, "pt"), y = unit(1, "npc"), just = c("left", "top"), width = sum(tg$widths), height = sum(tg$heights))                                   # table at the top left, under its title
  ttl <- textGrob(title, x = unit(8, "pt"), y = unit(0.45, "npc"), just = c("left", "centre"), gp = gpar(fontsize = FIG_PT_MIN + 2, fontface = "bold", fontfamily = FIG_FONT, lineheight = 0.95))
  gridExtra::arrangeGrob(ttl, tg, ncol = 1, heights = unit.c(unit(0.2 + 0.15 * lengths(regmatches(title, gregexpr("\n", title))), "in"), unit(1, "null"))) }
dfD <- data.frame(a = sub("-", "–", as.character(tabD$pUP_bin), fixed = TRUE), b = trimws(format(tabD$n, big.mark = ",")), c = sprintf("%.2f", tabD$median_uIF), d = sprintf("%.2f", tabD$median_uIF_incl))
names(dfD) <- c("%uORF+", "Genes", "Median\nuIF", "Median\ninclusion-\nadjusted uIF")
dfE <- data.frame(a = c("uIF, all genes", "Inclusion-adjusted\nuIF, all genes", "uIF, genes with\n≥ 90% uORF+"), b = trimws(format(tabE$n, big.mark = ",")), c = sig3(tabE$p05), d = sig3(tabE$p95), e = sprintf("%.0f-fold", tabE$ratio_p95_p05))
names(dfE) <- c("Measure", "Genes", "5th\npercentile", "95th\npercentile", "95th / 5th\npercentile")
gD <- tbl_panel(dfD, "D medians by %uORF+ group", c(0.68, 0.48, 0.55, 0.88)); gE <- tbl_panel(dfE, "E differences among genes", c(1.2, 0.48, 0.64, 0.64, 0.68))
H1 <- 2.55; H2 <- 1.8                                                                                                              # row of the plots; row of the tables (in). V54 (2026-11-23, Larry: "plots A-C are too tall. I want them to be more square-shaped (e.g., Fig 2)"): 3.1 -> 2.55 in
figS3 <- gridExtra::arrangeGrob(fig_grid(cS7, ncol = 3, h = H1), gridExtra::arrangeGrob(gD, gE, ncol = 2, widths = c(2.72, 3.78)), ncol = 1, heights = c(H1, H2), top = fig_title("Figure S4. Inclusion-adjusted uIF"))
fig_save(file.path(OUT_DIR, "A3_Figure_uIF_inclusion_adjusted.pdf"), figS3, h = H1 + H2 + 0.2)
.add(""); .add(sprintf("A3 done in %.1f min. Outputs: A3_uIF_tests_summary.txt, A3_correlation_table.csv, A3_permutation_null.csv, A3_single_TuORF_cohort_uIF_variants.csv, A3_uIF_by_inclusion_group.csv, A3_uIF_range_among_genes.csv, A3_Figure_uIF_independence_and_inclusion.pdf, A3_Figure_uIF_inclusion_adjusted.pdf", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
writeLines(c(sprintf("A3 uIF tests — generated %s", format(Sys.time(), "%Y-%m-%d %H:%M")), .log), file.path(OUT_DIR, "A3_uIF_tests_summary.txt"))
