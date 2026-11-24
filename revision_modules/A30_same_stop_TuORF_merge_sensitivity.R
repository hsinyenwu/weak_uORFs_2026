################################################################################
# A30 (2026-11-08, V23; V22 full review) — TuORF calls that share a start or stop codon: how many distinct uORFs, and do the single-TuORF
# results change when such calls are merged into one uORF?
# Table S1A lists calls as separate TuORFs when they share gene, strand and stop codon but start at different AUGs (often one call by each
# caller, the CiPS start upstream), or share the start codon but end at different stop codons (isoform versions). Here the calls of a gene
# are joined into one uORF when they share a start or a stop codon (connected groups); the call with the most uORF P-sites represents the
# group (ties: the longer call). The single-TuORF cohort is rebuilt with the rules of the canonical script (exactly one uORF in the gene; no
# translated overlapping uORF; %uORF+ > 0; uIF > 0) and compared with the published cohort.
# Outputs: output/A30_Table_TuORF_groups.csv, output/A30_Table_cohort_comparison.csv, output/A30_summary.txt
# Run from Weak_uORF_submission/:  Rscript revision_modules/A30_same_stop_TuORF_merge_sensitivity.R
################################################################################
suppressPackageStartupMessages({ library(data.table); library(openxlsx) })
.log <- character(); .add <- function(...) { m <- sprintf(...); cat(m, "\n"); .log <<- c(.log, m) }
S1A <- as.data.table(read.xlsx("output/Table_S1_all_TuORFs.xlsx", sheet = "Table S1A (CAGE)"))
ou_genes <- unique(trimws(read.xlsx("Files/CiPS_ouORF_all_Final.xlsx")$gene_id))
S1A[, row := .I]
# connected groups of calls sharing a start or a stop codon within a gene (union-find over the two keys)
par <- seq_len(nrow(S1A)); fnd <- function(i) { while (par[i] != i) { par[i] <<- par[par[i]]; i <- par[i] }; i }
uni <- function(a, b) { ra <- fnd(a); rb <- fnd(b); if (ra != rb) par[rb] <<- ra }
for (k in list(S1A[, paste(gene_id, Strand, "start", uORF_start_genomic)], S1A[, paste(gene_id, Strand, "stop", uORF_stop_genomic)])) {
  for (idx in split(seq_len(nrow(S1A)), k)) if (length(idx) > 1) for (j in idx[-1]) uni(idx[1], j) }
S1A[, grp := vapply(row, fnd, integer(1))]; S1A[, n_in_group := .N, by = grp]
S1A[, shares_stop := .N > 1, by = .(gene_id, Strand, uORF_stop_genomic)]; S1A[, shares_start := .N > 1, by = .(gene_id, Strand, uORF_start_genomic)]
n_stop_groups <- S1A[shares_stop == TRUE, uniqueN(paste(gene_id, Strand, uORF_stop_genomic))]
.add("A30: %d TuORFs in %d genes; %d TuORFs share the stop codon with another call (%d stop codons), %d share the start codon; %d distinct uORFs after joining calls that share a start or stop codon (%d groups with more than one call)",
     nrow(S1A), uniqueN(S1A$gene_id), sum(S1A$shares_stop), n_stop_groups, sum(S1A$shares_start), uniqueN(S1A$grp), S1A[n_in_group > 1, uniqueN(grp)])
rt_only <- S1A[detection_method == "RiboTaper"]; ci_only <- S1A[detection_method == "CiPS"]
pairs <- merge(rt_only[, .(gene_id, Strand, uORF_stop_genomic, rt_start = uORF_start_genomic)], ci_only[, .(gene_id, Strand, uORF_stop_genomic, ci_start = uORF_start_genomic)], by = c("gene_id", "Strand", "uORF_stop_genomic"))
.add("  RiboTaper-only TuORFs that share the stop codon with a CiPS-only TuORF: %d of %d (the CiPS start codon lies upstream in %d of %d pairs)", uniqueN(pairs[, paste(gene_id, rt_start)]), nrow(rt_only),
     sum(ifelse(pairs$Strand == "+", pairs$ci_start < pairs$rt_start, pairs$ci_start > pairs$rt_start)), nrow(pairs))
fwrite(S1A[, .(gene_id, tx_id, uORF_genome_range, detection_method, uORF_length_aa, uORF_P_sites, pct_uORF_plus, group = paste0("G", grp), n_in_group, shares_start, shares_stop)], "output/A30_Table_TuORF_groups.csv")
# representative call per group, then the cohort rules
rep <- S1A[order(grp, -uORF_P_sites, -uORF_length_nt)][, .SD[1], by = grp]
build <- function(d) { cnt <- d[, .N, by = gene_id]; one <- cnt[N == 1, gene_id]
  c0 <- d[gene_id %in% one & !gene_id %in% ou_genes]; n_ou <- length(intersect(one, ou_genes))
  c1 <- c0[pct_uORF_plus > 0]; c1[, uIF_raw := ifelse(is.finite(mORF_TL) & mORF_TL > 0, uORF_TL / mORF_TL, NA_real_)]; c2 <- c1[is.finite(uIF_raw) & uIF_raw > 0]
  list(d = c2, steps = c(single = length(one), minus_ouORF = nrow(c0), minus_pUP0 = nrow(c1), final = nrow(c2))) }
pub <- build(S1A); mrg <- build(rep)
.add("  cohort rules on the published definition: %d genes with one TuORF -> %d without an overlapping uORF -> %d with %%uORF+ > 0 -> %d with uIF > 0 (published cohort: 2,996)", pub$steps[1], pub$steps[2], pub$steps[3], pub$steps[4])
.add("  merged definition: %d genes with one uORF -> %d -> %d -> %d (%d genes added, %d lost)", mrg$steps[1], mrg$steps[2], mrg$steps[3], mrg$steps[4], length(setdiff(mrg$d$gene_id, pub$d$gene_id)), length(setdiff(pub$d$gene_id, mrg$d$gene_id)))
summ <- function(d, nm) { g <- cut(d$pct_uORF_plus, c(-Inf, 30, 80, Inf), right = FALSE, labels = c("<30%", "30-80%", "80-100%")); te <- is.finite(d$mORF_TE) & d$mORF_TE > 0
  data.table(cohort = nm, genes = nrow(d), pct_below30 = 100 * mean(g == "<30%"), pct_30_80 = 100 * mean(g == "30-80%"), pct_80_100 = 100 * mean(g == "80-100%"), pct_below90 = 100 * mean(d$pct_uORF_plus < 90),
             r_pUP_log2uIF = cor(d$pct_uORF_plus, log2(d$uIF_raw)), genes_TE = sum(te), r_pUP_log2TE = cor(d$pct_uORF_plus[te], log2(d$mORF_TE[te])),
             medTE_below30 = median(d$mORF_TE[te & g == "<30%"]), medTE_30_80 = median(d$mORF_TE[te & g == "30-80%"]), medTE_80_100 = median(d$mORF_TE[te & g == "80-100%"])) }
cmp <- rbind(summ(pub$d, "published definition (calls listed separately)"), summ(mrg$d, "calls sharing a start or stop codon merged"))
for (i in seq_len(nrow(cmp))) with(cmp[i], .add("  %-48s n = %d: <30%% %.1f%%, 30-80%% %.1f%%, 80-100%% %.1f%%; below 90%% %.1f%%; r(%%uORF+, log2 uIF) = %.3f; r(%%uORF+, log2 TE) = %.3f (n = %d); median TE %.2f / %.2f / %.2f",
  cohort, genes, pct_below30, pct_30_80, pct_80_100, pct_below90, r_pUP_log2uIF, r_pUP_log2TE, genes_TE, medTE_below30, medTE_30_80, medTE_80_100))
fwrite(cmp, "output/A30_Table_cohort_comparison.csv")
writeLines(c(sprintf("A30 TuORF calls sharing a start or stop codon — generated %s", format(Sys.time(), "%Y-%m-%d %H:%M")), .log), "output/A30_summary.txt"); cat("A30 done\n")
