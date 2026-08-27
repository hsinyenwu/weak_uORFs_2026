################################################################################
# TSS heterogeneity and the uORF landscape in Arabidopsis
#
# Analysis code for:
#   "Transcription Start Site Heterogeneity Confounds the Landscape and
#    Functional Interpretation of uORFs" (Hsu lab)
#
# For each translated uORF (TuORF), the fraction of transcripts that include
# the uORF (%uORF+, called "pUP" in the code) is computed from TSS profiling
# data (CAGE-seq, or TSS-seq from Nielsen et al. 2019) and combined with
# Ribo-seq-based uORF identification (RiboTaper + CiPS), RNA-seq/Ribo-seq
# quantification, reporter assays, and public datasets to ask how transcription
# start site selection shapes uORF inclusion and function.
#
# See README.md for a section-by-section overview, the list of required input
# files, and instructions for running the script.
################################################################################
rm(list=ls())
gc()
options(stringsAsFactors = FALSE)

################################################################################
# 0) INPUT / OUTPUT LOCATIONS
################################################################################
## All input files are expected in DATA_DIR and all results are written to
## OUT_DIR_MAIN. Both are relative paths: run this script from the repository
## root (the folder containing this file), or setwd() to it first.
DATA_DIR     <- "Files"
OUT_DIR_MAIN <- "output"

## --- uORF identification and annotation -------------------------------------
PATH_RT_ORFS <- file.path(DATA_DIR, "ORFs_max_filt_expressed")       # RiboTaper ORF calls (expressed set)
PATH_CIPS    <- file.path(DATA_DIR, "CiPS_uORF_Final.xlsx")          # CiPS TuORF calls
PATH_GTF     <- file.path(DATA_DIR, "Araport11+CTRL_20181206.gtf")   # Araport11 annotation used throughout
PATH_GTF2    <- file.path(DATA_DIR, "Araport11_20181206_max2_expressed_isoform_TPM0.5_rm_diff_CDSstart.gtf")
                                                                     # expressed-isoform annotation (Figure 4 TxDb)
PATH_RIBO_P  <- file.path(DATA_DIR, "CTRL_expressed_P_sites_sort_count")  # Ribo-seq P-site positions/counts

## --- TSS data source used to compute %uORF+ (pUP) ---------------------------
##   "CAGE" = CAGE-seq  (basis of the manuscript figures)
##   "TSS"  = TSS-seq   (Nielsen et al. 2019; Figure S2 companion run)
TSS_SOURCE         <- "CAGE"
PATH_CAGE_CAGEseq  <- file.path(DATA_DIR, "wt_R123_genome_RiboPlotR.txt")
PATH_CAGE_TSSseq   <- file.path(DATA_DIR, "wt_TSS_Nielsen_2019.txt")
PATH_CAGE          <- if (TSS_SOURCE == "CAGE") PATH_CAGE_CAGEseq else PATH_CAGE_TSSseq

## --- expression, translation, and gene sets ---------------------------------
PATH_RSEM_RNA      <- file.path(DATA_DIR, "RSEM_Araport11_1st+2nd_Riboseq_RNA.genes.results_CDS_only")   # RNA-seq RSEM (gene level)
PATH_RSEM_RIBO     <- file.path(DATA_DIR, "RSEM_Araport11_1st+2nd_Riboseq_Ribo.genes.results_CDS_only")  # Ribo-seq RSEM (gene level)
PATH_5UTR_GENES    <- file.path(DATA_DIR, "genes_w_5UTR.tsv")             # genes with an annotated 5'UTR
PATH_CPUORF_MANUAL <- file.path(DATA_DIR, "CPuORF_collection_3_manully_fix.xlsx")  # curated CPuORF collection
PATH_FASTA_TAIR10  <- file.path(DATA_DIR, "TAIR10_chr_all_2.fas")         # TAIR10 genome sequence

## --- public datasets and reporter-assay measurements -------------------------
PATH_HALFLIFE_TSV  <- file.path(DATA_DIR, "RNA-half-lives-5EU.tsv")       # RNA half-lives (5-EU labeling dataset)
PATH_HALFLIFE_TXT  <- file.path(DATA_DIR, "RNA-half-lives-5EU.txt")       # same dataset, txt export
PATH_GOSLIM        <- file.path(DATA_DIR, "ATH_GO_GOSLIM.tsv")            # TAIR GO / GOslim annotations
PATH_RSEM_ISOFORMS <- file.path(DATA_DIR, "RSEM_isoforms_TPM_6samples.tsv")  # per-isoform RSEM TPM, 6 samples in one table
PATH_DLA_WEAK_UORF_UUORF  <- file.path(DATA_DIR, "DLA_weak_uORF_UuORFs.xlsx")               # dual-luciferase assays
PATH_QPCR_WEAK_UORF_UUORF <- file.path(DATA_DIR, "qPCR_weak_uORF_uuORF.xlsx")               # qPCR
PATH_QPCR_PERFECT_KOZAK   <- file.path(DATA_DIR, "qPCR_weak_uORF_uuORF_perfect_Kozak.xlsx") # qPCR (perfect-Kozak constructs)
################################################################################
# 0b) OUTPUT HELPERS
################################################################################
OUT_DIR <- OUT_DIR_MAIN
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

set_out_dir <- function(path) {
  OUT_DIR <<- path
  if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
  invisible(OUT_DIR)
}

save_pdf <- function(plot_obj, filename, w = 7, h = 5) {
  ggplot2::ggsave(
    filename = file.path(OUT_DIR, filename),
    plot = plot_obj,
    width = w, height = h
  )
}

################################################################################
# 0c) LIBRARIES (ONCE)
################################################################################
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
  
  library(ggplot2)
  library(ggpubr)
  library(rstatix)
  library(scales)
  
  library(openxlsx)
  
  library(GenomicFeatures)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
  library(S4Vectors)
  library(BiocGenerics)
  library(Biostrings)
  library(Rsamtools)
  library(ORFik)
  
  library(rcompanion)
  library(multcompView)
})

################################################################################
# 0d) SMALL FILE CHECKS
################################################################################
stopifnot(file.exists(PATH_RT_ORFS))
stopifnot(file.exists(PATH_CIPS))
stopifnot(file.exists(PATH_GTF))
stopifnot(file.exists(PATH_RIBO_P))
stopifnot(file.exists(PATH_CAGE))
stopifnot(file.exists(PATH_RSEM_RNA))
stopifnot(file.exists(PATH_RSEM_RIBO))
stopifnot(file.exists(PATH_5UTR_GENES))
stopifnot(file.exists(PATH_CPUORF_MANUAL))
stopifnot(file.exists(PATH_FASTA_TAIR10))
stopifnot(file.exists(PATH_HALFLIFE_TSV))
stopifnot(file.exists(PATH_HALFLIFE_TXT))
stopifnot(file.exists(PATH_GOSLIM))
stopifnot(file.exists(PATH_RSEM_ISOFORMS))
stopifnot(file.exists(PATH_DLA_WEAK_UORF_UUORF))
stopifnot(file.exists(PATH_QPCR_WEAK_UORF_UUORF))
stopifnot(file.exists(PATH_QPCR_PERFECT_KOZAK))

################################################################################
# 0e) HELPERS
################################################################################
get_gene_id_from_tx <- function(tx_id) substr(tx_id, 1, 9)

clean_peptide <- function(x) {
  x <- as.character(x)
  x <- toupper(x)
  x <- gsub("\\s+", "", x)
  x <- gsub("[\\*X]", "", x)
  x[nchar(x) == 0] <- NA_character_
  x
}

# ---------- KS letters for ECDF legends ----------
get_pairwise_ks_pvalues_clean <- function(data, group_col, val_col, id_map) {
  groups <- levels(data[[group_col]])
  n_g <- length(groups)
  clean_names <- id_map[groups]
  p_mat <- matrix(1, nrow = n_g, ncol = n_g, dimnames = list(clean_names, clean_names))
  if (n_g < 2) return(p_mat)
  
  combos <- utils::combn(seq_len(n_g), 2, simplify = FALSE)
  
  p_values <- vapply(combos, function(idx) {
    g1 <- groups[idx[1]]
    g2 <- groups[idx[2]]
    d1 <- as.numeric(stats::na.omit(data[[val_col]][data[[group_col]] == g1]))
    d2 <- as.numeric(stats::na.omit(data[[val_col]][data[[group_col]] == g2]))
    if (length(d1) > 0 && length(d2) > 0) stats::ks.test(d1, d2)$p.value else 1
  }, numeric(1))
  
  idx_mat <- do.call(rbind, combos)
  i <- idx_mat[, 1]; j <- idx_mat[, 2]
  p_mat[cbind(i, j)] <- p_values
  p_mat[cbind(j, i)] <- p_values
  p_mat
}

reorder_letters_by_appearance <- function(named_letters, ordered_groups) {
  current_sequence <- named_letters[ordered_groups]
  unique_codes <- unique(current_sequence)
  new_alphabet <- letters[1:length(unique_codes)]
  mapping <- setNames(new_alphabet, unique_codes)
  new_letters <- mapping[named_letters]
  names(new_letters) <- names(named_letters)
  new_letters
}

generate_legend_labels <- function(data_df, group_col, val_col) {
  group_levels <- levels(data_df[[group_col]])
  clean_ids <- paste0("G", seq_along(group_levels))
  names(clean_ids) <- group_levels
  
  p_mat <- get_pairwise_ks_pvalues_clean(data_df, group_col, val_col, clean_ids)
  
  p_vec <- p_mat[lower.tri(p_mat)]
  names_vec <- expand.grid(rownames(p_mat), colnames(p_mat))
  names_vec <- names_vec[lower.tri(p_mat), ]
  names(p_vec) <- paste(names_vec$Var1, names_vec$Var2, sep = "-")
  
  letters_list <- multcompView::multcompLetters(
    p_vec, compare = "<", threshold = 0.05, Letters = letters, reversed = FALSE
  )
  
  clean_letters <- letters_list$Letters
  clean_to_real <- setNames(names(clean_ids), clean_ids)
  initial_letters <- clean_letters
  names(initial_letters) <- clean_to_real[names(clean_letters)]
  
  final_letters <- reorder_letters_by_appearance(initial_letters, levels(data_df[[group_col]]))
  
  counts <- data_df %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")
  
  counts$letters <- final_letters[as.character(counts[[group_col]])]
  counts$Final_Label <- paste0(counts[[group_col]], " (n=", counts$n, ") ", counts$letters)
  counts
}

# ---------- Wilcoxon letters for boxplots ----------
get_wilcox_letters <- function(df, y_col, group_col, p_adjust = "BY") {
  df2 <- df %>%
    dplyr::filter(!is.na(.data[[group_col]])) %>%
    dplyr::filter(!is.na(.data[[y_col]]), is.finite(.data[[y_col]]))
  
  grps <- unique(as.character(df2[[group_col]]))
  if (length(grps) < 2) return(NULL)
  
  pw <- tryCatch(
    stats::pairwise.wilcox.test(df2[[y_col]], df2[[group_col]], p.adjust.method = p_adjust, exact = FALSE),
    error = function(e) NULL
  )
  if (is.null(pw) || is.null(pw$p.value)) return(NULL)
  
  pmat <- pw$p.value
  pmat[is.na(pmat)] <- 1
  full <- rcompanion::fullPTable(pmat)
  lets <- multcompView::multcompLetters(full, compare = "<", threshold = 0.05, Letters = letters, reversed = FALSE)
  
  data.frame(Group = names(lets$Letters), Letters = lets$Letters, stringsAsFactors = FALSE)
}

add_counts_labels <- function(df, group_col, group_order) {
  counts <- df %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(.ord = match(as.character(.data[[group_col]]), group_order)) %>%
    dplyr::arrange(.ord)
  
  counts$Group_Label <- paste0(counts[[group_col]], "\n(n=", counts$n, ")")
  list(counts = counts, label_map = setNames(counts$Group_Label, as.character(counts[[group_col]])))
}

# ---------- mapping helpers ----------
map_from_transcripts_to_grl <- function(gr_tx, exons_by_tx) {
  stopifnot(is(gr_tx, "GRanges"))
  stopifnot(is(exons_by_tx, "GRangesList"))
  
  mapped <- tryCatch(
    GenomicFeatures::mapFromTranscripts(gr_tx, exons_by_tx, ignore.strand = FALSE),
    error = function(e) GenomicFeatures::mapFromTranscripts(gr_tx, exons_by_tx)
  )
  
  if (is(mapped, "GRangesList")) {
    if (is.null(names(mapped))) names(mapped) <- names(gr_tx)[seq_along(mapped)]
    return(mapped)
  }
  
  if (is(mapped, "GRanges")) {
    mc <- S4Vectors::mcols(mapped)
    if (!("xHits" %in% colnames(mc))) stop("mapFromTranscripts returned GRanges without xHits.")
    sp <- split(mapped, mc$xHits)
    idx <- suppressWarnings(as.integer(names(sp)))
    if (any(is.na(idx))) stop("Could not interpret xHits indices.")
    names(sp) <- names(gr_tx)[idx]
    return(GenomicRanges::GRangesList(sp))
  }
  
  stop("Unsupported return type from mapFromTranscripts().")
}

span_df_from_grl <- function(grl) {
  if (length(grl) == 0) {
    return(data.frame(orf_id = character(0), chr = character(0), g_min = integer(0), g_max = integer(0)))
  }
  
  lens <- S4Vectors::elementNROWS(grl)
  orf_ids <- rep(names(grl), times = lens)
  flat <- unlist(grl, use.names = FALSE)
  
  if (length(flat) == 0) {
    return(data.frame(orf_id = character(0), chr = character(0), g_min = integer(0), g_max = integer(0)))
  }
  
  dt <- data.table::data.table(
    orf_id = orf_ids,
    chr    = as.character(GenomicRanges::seqnames(flat)),
    st     = as.integer(GenomicRanges::start(flat)),
    en     = as.integer(GenomicRanges::end(flat))
  )
  
  out <- dt[, .(
    chr   = chr[1],
    g_min = as.integer(min(st, na.rm = TRUE)),
    g_max = as.integer(max(en, na.rm = TRUE))
  ), by = orf_id]
  
  as.data.frame(out)
}

standardize_evidence_types <- function(df) {
  df %>%
    dplyr::mutate(
      gene_id       = as.character(gene_id),
      tx_id         = as.character(tx_id),
      Strand        = as.character(Strand),
      uORF_id        = as.character(uORF_id),
      uORF_id_short  = as.character(uORF_id_short),
      genome_range   = as.character(genome_range),
      uORF_Tx_range  = as.character(uORF_Tx_range),
      peptide_seq    = as.character(peptide_seq),
      method         = as.character(method),
      ORF_gstart     = as.integer(ORF_gstart),
      ORF_gstop      = as.integer(ORF_gstop),
      uORF_P_sites   = as.numeric(uORF_P_sites),
      uORF_length    = as.integer(uORF_length),
      uORF_length_nt = as.integer(uORF_length_nt),
      uORF_TL        = as.numeric(uORF_TL)
    )
}

normalize_chr_to_ref <- function(chr_vec, ref_seqlevels) {
  chr_vec <- as.character(chr_vec)
  if (length(ref_seqlevels) == 0) return(chr_vec)
  
  ref_has_prefix <- any(grepl("^(chr|Chr)", ref_seqlevels))
  in_has_prefix  <- any(grepl("^(chr|Chr)", chr_vec))
  
  if (ref_has_prefix && !in_has_prefix) {
    cand1 <- paste0("Chr", chr_vec)
    cand2 <- paste0("chr", chr_vec)
    return(ifelse(cand1 %in% ref_seqlevels, cand1,
                  ifelse(cand2 %in% ref_seqlevels, cand2, chr_vec)))
  }
  if (!ref_has_prefix && in_has_prefix) return(gsub("^(chr|Chr)", "", chr_vec))
  chr_vec
}

map_points_to_txcoord <- function(gr_points, tx_lookup, exonByTx) {
  stopifnot(is(gr_points, "GRanges"))
  stopifnot(length(gr_points) == length(tx_lookup))
  stopifnot(all(names(gr_points) == names(tx_lookup)))
  
  used_txs <- unique(unname(tx_lookup))
  used_txs <- used_txs[!is.na(used_txs)]
  used_txs <- used_txs[used_txs %in% names(exonByTx)]
  sub_exons <- exonByTx[used_txs]
  
  out <- setNames(rep(NA_integer_, length(gr_points)), names(gr_points))
  if (length(sub_exons) == 0 || length(gr_points) == 0) return(out)
  
  mapped <- GenomicFeatures::mapToTranscripts(gr_points, sub_exons, ignore.strand = FALSE)
  if (length(mapped) == 0) return(out)
  
  mc <- S4Vectors::mcols(mapped)
  if (!("xHits" %in% colnames(mc)) || !("transcriptsHits" %in% colnames(mc))) {
    mapped <- GenomicFeatures::mapToTranscripts(gr_points, sub_exons)
    if (length(mapped) == 0) return(out)
    mc <- S4Vectors::mcols(mapped)
    if (!("xHits" %in% colnames(mc)) || !("transcriptsHits" %in% colnames(mc))) {
      stop("mapToTranscripts() output missing xHits/transcriptsHits.")
    }
  }
  
  point_id <- names(gr_points)[mc$xHits]
  tx_mapped <- names(sub_exons)[mc$transcriptsHits]
  tx_expected <- tx_lookup[point_id]
  tx_pos <- as.integer(GenomicRanges::start(mapped))
  
  dt <- data.table::data.table(
    point_id = point_id,
    tx_mapped = tx_mapped,
    tx_expected = tx_expected,
    tx_pos = tx_pos
  )
  
  dt <- dt[tx_mapped == tx_expected]
  if (nrow(dt) == 0) return(out)
  
  dt <- dt[, .(tx_pos = as.integer(min(tx_pos))), by = point_id]
  out[dt$point_id] <- dt$tx_pos
  out
}

fill_uORF_Tx_range_by_endpoints <- function(df, exonByTx, ref_seqlevels) {
  if (!"uORF_Tx_range" %in% colnames(df)) df$uORF_Tx_range <- NA_character_
  df$uORF_Tx_range <- as.character(df$uORF_Tx_range)
  
  need <- which(is.na(df$uORF_Tx_range) | df$uORF_Tx_range == "" | df$uORF_Tx_range == "NA")
  if (length(need) == 0) return(df)
  
  ok <- need[
    !is.na(df$tx_id[need]) & df$tx_id[need] %in% names(exonByTx) &
      !is.na(df$chr[need]) &
      !is.na(df$ORF_gstart[need]) & !is.na(df$ORF_gstop[need]) &
      !is.na(df$Strand[need])
  ]
  if (length(ok) == 0) return(df)
  
  chr_use <- normalize_chr_to_ref(df$chr[ok], ref_seqlevels)
  sid <- paste0("S", ok); eid <- paste0("E", ok)
  
  grS <- GenomicRanges::GRanges(
    seqnames = chr_use,
    ranges   = IRanges::IRanges(start = as.integer(df$ORF_gstart[ok]), width = 1),
    strand   = as.character(df$Strand[ok])
  ); names(grS) <- sid
  
  grE <- GenomicRanges::GRanges(
    seqnames = chr_use,
    ranges   = IRanges::IRanges(start = as.integer(df$ORF_gstop[ok]), width = 1),
    strand   = as.character(df$Strand[ok])
  ); names(grE) <- eid
  
  tx_lookupS <- stats::setNames(as.character(df$tx_id[ok]), sid)
  tx_lookupE <- stats::setNames(as.character(df$tx_id[ok]), eid)
  
  txS <- map_points_to_txcoord(grS, tx_lookupS, exonByTx)[sid]
  txE <- map_points_to_txcoord(grE, tx_lookupE, exonByTx)[eid]
  
  lo <- pmin(txS, txE)
  hi <- pmax(txS, txE)
  
  new_rng <- ifelse(!is.na(lo) & !is.na(hi), paste0(lo, "-", hi), NA_character_)
  df$uORF_Tx_range[ok] <- ifelse(!is.na(new_rng), new_rng, df$uORF_Tx_range[ok])
  
  df
}

################################################################################
# 1) LOAD RT + CiPS + TxDb
################################################################################
set_out_dir(OUT_DIR_MAIN)

ORFs_max_filt <- read.delim(PATH_RT_ORFS, header = TRUE, sep = "\t", quote = "")
ORFs_max_filt_uORF <- ORFs_max_filt %>% dplyr::filter(category == "uORF")

CiPS <- openxlsx::read.xlsx(PATH_CIPS)

# normalize CiPS column names
if ("AA" %in% colnames(CiPS) && !"pept_length" %in% colnames(CiPS)) {
  CiPS <- CiPS %>% dplyr::rename(pept_length = AA)
}
if (!"uORF_Tx_range" %in% colnames(CiPS)) CiPS$uORF_Tx_range <- NA_character_

txdb <- GenomicFeatures::makeTxDbFromGFF(
  file = PATH_GTF, format = "gtf",
  dataSource = "Araport11",
  organism = "Arabidopsis"
)

# strand info
Strand_info <- rbind(
  ORFs_max_filt_uORF %>%
    dplyr::select(transcript_id, strand) %>%
    dplyr::rename(tx_id = transcript_id, Strand = strand),
  CiPS %>% dplyr::select(tx_id, Strand)
) %>% dplyr::distinct()

# detected gene set (for uuORF exclusion + Figure4 TE)
detected_genes <- unique(c(as.character(ORFs_max_filt_uORF$gene_id), as.character(CiPS$gene_id)))

################################################################################
# 2) mORF CDS metrics (P-sites / CDS length)
################################################################################
Ribo1 <- read.delim(PATH_RIBO_P, header = FALSE)
colnames(Ribo1) <- c("count", "chr", "start", "strand")
Ribo1$end <- Ribo1$start
RiboR <- GenomicRanges::makeGRangesFromDataFrame(Ribo1, keep.extra.columns = TRUE)

target_tx <- unique(c(ORFs_max_filt_uORF$transcript_id, as.character(CiPS$tx_id)))
cdsRanges <- GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)
cdsRanges <- cdsRanges[names(cdsRanges) %in% target_tx]

CDSLength <- vapply(cdsRanges, function(gr) sum(BiocGenerics::width(gr)), numeric(1))

cds_unlist <- BiocGenerics::unlist(cdsRanges, use.names = FALSE)
cds_unlist$tx_index <- rep(seq_along(cdsRanges), times = S4Vectors::elementNROWS(cdsRanges))

ol_cds <- GenomicRanges::findOverlaps(RiboR, cds_unlist)
agg_cds <- stats::aggregate(
  RiboR$count[S4Vectors::queryHits(ol_cds)] ~ cds_unlist$tx_index[S4Vectors::subjectHits(ol_cds)],
  FUN = sum
)

CDScountsInf <- numeric(length(cdsRanges))
if (nrow(agg_cds) > 0) CDScountsInf[agg_cds[, 1]] <- agg_cds[, 2]
names(CDScountsInf) <- names(cdsRanges)

mORFInfo <- data.frame(
  tx_id = names(CDScountsInf),
  mORF_Ribo_count = unname(CDScountsInf),
  mORF_length = unname(CDSLength),
  mORF_TL = unname(CDScountsInf / CDSLength)
)

################################################################################
# 2B) Merge RT + CiPS evidence (robust, consistent)
################################################################################
# format CiPS
CiPS <- CiPS %>%
  dplyr::mutate(
    gene_id    = as.character(gene_id),
    tx_id      = as.character(tx_id),
    Strand     = as.character(Strand),
    chr        = as.character(chr),
    ORF_gstart = as.integer(ORF_gstart),
    ORF_gstop  = as.integer(ORF_gstop),
    genome_range = paste0(chr, "_", pmin(ORF_gstart, ORF_gstop), "_", pmax(ORF_gstart, ORF_gstop)),
    uORF_id      = paste(tx_id, genome_range, sep = "_"),
    uORF_Tx_range = as.character(uORF_Tx_range),
    peptide_seq  = as.character(peptide_seq)
  )

tx_exons <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
ref_seqlevels <- tryCatch(GenomeInfoDb::seqlevels(BiocGenerics::unlist(tx_exons)),
                          error = function(e) character(0))

CiPS <- fill_uORF_Tx_range_by_endpoints(CiPS, tx_exons, ref_seqlevels)

# RT meta in tx coords
DF_RT_meta <- ORFs_max_filt_uORF %>%
  dplyr::transmute(
    gene_id = as.character(gene_id),
    tx_id   = as.character(transcript_id),
    Strand  = as.character(strand),
    
    uORF_id_short = as.character(ORF_id_tr),
    
    tx_start = as.integer(start_pos),
    tx_end_full = as.integer(stop_pos) + 2L,
    uORF_Tx_range = paste0(tx_start, "-", tx_end_full),
    
    pept_length = as.integer(nchar(clean_peptide(ORF_pept))),
    peptide_seq = paste0(clean_peptide(ORF_pept), "*"),
    
    uORF_P_sites = as.numeric(ORF_P_sites),
    uORF_length_nt = as.integer(ORF_length),
    method = "RT"
  ) %>%
  dplyr::filter(!is.na(tx_start), !is.na(tx_end_full)) %>%
  dplyr::filter(tx_id %in% names(tx_exons))

gr_tx <- GenomicRanges::GRanges(
  seqnames = DF_RT_meta$tx_id,
  ranges   = IRanges::IRanges(start = DF_RT_meta$tx_start, end = DF_RT_meta$tx_end_full),
  strand   = "*"
)
names(gr_tx) <- DF_RT_meta$uORF_id_short

RT_grl <- map_from_transcripts_to_grl(gr_tx, tx_exons)
RT_span_df <- span_df_from_grl(RT_grl) %>% dplyr::rename(uORF_id_short = orf_id)

DF_RT <- DF_RT_meta %>%
  dplyr::left_join(RT_span_df, by = "uORF_id_short") %>%
  dplyr::filter(!is.na(chr), !is.na(g_min), !is.na(g_max)) %>%
  dplyr::mutate(
    genome_range = paste0(chr, "_", g_min, "_", g_max),
    uORF_id = paste(tx_id, genome_range, sep = "_"),
    ORF_gstart = as.integer(ifelse(Strand == "+", g_min, g_max)),
    ORF_gstop  = as.integer(ifelse(Strand == "+", g_max, g_min)),
    uORF_length = as.integer(pept_length),
    uORF_TL = dplyr::if_else(uORF_length_nt > 0, uORF_P_sites / uORF_length_nt, NA_real_)
  ) %>%
  dplyr::select(
    gene_id, tx_id, Strand,
    uORF_id, uORF_id_short,
    genome_range,
    ORF_gstart, ORF_gstop,
    uORF_Tx_range,
    uORF_P_sites, uORF_length, uORF_length_nt,
    peptide_seq, uORF_TL,
    method
  )

DF_CiPS <- CiPS %>%
  dplyr::mutate(
    uORF_P_sites   = as.numeric(TotalC),
    uORF_length    = as.integer(pept_length),
    uORF_length_nt = as.integer(pept_length) * 3L,
    uORF_TL        = as.numeric(TotalC) / (as.integer(pept_length) * 3L),
    method         = "CiPS",
    uORF_id_short  = as.character(uORF_id),
    uORF_Tx_range  = as.character(uORF_Tx_range)
  ) %>%
  dplyr::select(
    gene_id, tx_id, Strand,
    uORF_id, uORF_id_short,
    genome_range,
    ORF_gstart, ORF_gstop,
    uORF_Tx_range,
    uORF_P_sites, uORF_length, uORF_length_nt,
    peptide_seq, uORF_TL,
    method
  )

DF_RT2   <- standardize_evidence_types(DF_RT)
DF_CiPS2 <- standardize_evidence_types(DF_CiPS)

# strong overlap key: gene + genome_range + strand
RT_key   <- paste(DF_RT2$gene_id, DF_RT2$genome_range, DF_RT2$Strand, sep = "|")
CiPS_key <- paste(DF_CiPS2$gene_id, DF_CiPS2$genome_range, DF_CiPS2$Strand, sep = "|")
overlap_keys <- intersect(unique(RT_key), unique(CiPS_key))

DF_RT3 <- DF_RT2 %>%
  dplyr::mutate(.key = paste(gene_id, genome_range, Strand, sep = "|"),
                method = dplyr::if_else(.key %in% overlap_keys, "Both", "RT")) %>%
  dplyr::select(-.key)

DF_CiPS_only <- DF_CiPS2 %>%
  dplyr::mutate(.key = paste(gene_id, genome_range, Strand, sep = "|")) %>%
  dplyr::filter(!(.key %in% overlap_keys)) %>%
  dplyr::select(-.key) %>%
  dplyr::mutate(method = "CiPS")

DF_All_Evidence <- dplyr::bind_rows(DF_RT3, DF_CiPS_only) %>%
  dplyr::mutate(
    chr  = vapply(strsplit(genome_range, "_"), `[`, 1, FUN.VALUE = character(1)),
    g_min = as.integer(vapply(strsplit(genome_range, "_"), `[`, 2, FUN.VALUE = character(1))),
    g_max = as.integer(vapply(strsplit(genome_range, "_"), `[`, 3, FUN.VALUE = character(1))),
    uORF_start_bio = ifelse(Strand == "+", g_min, g_max),
    uORF_stop_bio  = ifelse(Strand == "+", g_max, g_min)
  )

mORF_uORF_Info <- dplyr::inner_join(mORFInfo, DF_All_Evidence, by = "tx_id")

################################################################################
# 3) pUP (CAGE) for detected uORFs (TuORFs)
################################################################################
CAGE_raw <- read.delim(PATH_CAGE, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
colnames(CAGE_raw) <- c("count", "chr", "start", "strand")
CAGE_raw$end <- CAGE_raw$start

# normalize CAGE chr -> NCBI-like by removing chr prefix
CAGE_raw$chr <- gsub("^chr|^Chr", "", CAGE_raw$chr, ignore.case = TRUE)
CAGE_gr <- GenomicRanges::makeGRangesFromDataFrame(CAGE_raw, keep.extra.columns = TRUE)

all_tx <- GenomicFeatures::transcripts(txdb, use.names = TRUE)
tss_windows <- all_tx[names(all_tx) %in% mORF_uORF_Info$tx_id]

# strand-aware window: upstream 100 + downstream full transcript
tss_windows <- suppressWarnings(GenomicRanges::promoters(
  tss_windows, upstream = 100, downstream = BiocGenerics::width(tss_windows)
))
tss_windows <- GenomicRanges::trim(tss_windows)

# try to align seqlevels styles
try({
  GenomeInfoDb::seqlevelsStyle(tss_windows) <- "NCBI"
  GenomeInfoDb::seqlevelsStyle(CAGE_gr) <- GenomeInfoDb::seqlevelsStyle(tss_windows)
}, silent = TRUE)

common_lv <- intersect(GenomeInfoDb::seqlevels(CAGE_gr), GenomeInfoDb::seqlevels(tss_windows))
if (length(common_lv) > 0) {
  CAGE_gr <- GenomeInfoDb::keepSeqlevels(CAGE_gr, common_lv, pruning.mode = "coarse")
  tss_windows <- GenomeInfoDb::keepSeqlevels(tss_windows, common_lv, pruning.mode = "coarse")
}

overlaps <- GenomicRanges::findOverlaps(CAGE_gr, tss_windows)

dt_cage <- data.table::data.table(
  tx_id = names(tss_windows)[S4Vectors::subjectHits(overlaps)],
  cage_coord = BiocGenerics::start(CAGE_gr)[S4Vectors::queryHits(overlaps)],
  count = CAGE_gr$count[S4Vectors::queryHits(overlaps)]
)

dt_targets <- data.table::data.table(
  tx_id = mORF_uORF_Info$tx_id,
  uORF_start_bio = mORF_uORF_Info$uORF_start_bio,
  Strand = mORF_uORF_Info$Strand,
  unique_key = paste(mORF_uORF_Info$tx_id, mORF_uORF_Info$genome_range, sep = "|")
)

dt_join <- dt_cage[dt_targets, on = .(tx_id), allow.cartesian = TRUE]

OUT_pUP <- dt_join[, .(
  pUP = {
    total_c <- sum(count, na.rm = TRUE)
    if (total_c == 0 || is.na(total_c)) 0 else {
      Cb <- if (Strand[1] == "+") sum(count[cage_coord <= uORF_start_bio], na.rm = TRUE)
      else                        sum(count[cage_coord >= uORF_start_bio], na.rm = TRUE)
      (Cb / total_c) * 100
    }
  }
), by = .(unique_key)]

mORF_uORF_Info$unique_key <- paste(mORF_uORF_Info$tx_id, mORF_uORF_Info$genome_range, sep="|")
OUT_df5 <- dplyr::left_join(mORF_uORF_Info, OUT_pUP, by = "unique_key")
OUT_df5$pUP[is.na(OUT_df5$pUP)] <- 0
OUT_df5$pUP <- round(OUT_df5$pUP, 3)

# redundancy check: keep best per gene & exact coords
OUT_df5_clean <- OUT_df5 %>%
  dplyr::mutate(
    gene_id = substr(tx_id, 1, 9),
    method_rank = dplyr::case_when(
      method == "Both" ~ 1L,
      method == "RT"   ~ 2L,
      TRUE             ~ 3L
    )
  ) %>%
  dplyr::group_by(gene_id, chr, uORF_start_bio, uORF_stop_bio, Strand) %>%
  dplyr::arrange(method_rank, dplyr::desc(uORF_P_sites), dplyr::desc(uORF_TL), .by_group = TRUE) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::select(-method_rank)


################################################################################
# NEW OUTPUT TABLE: uORFs with pUP > 30%
# (uses the full cleaned TuORF set OUT_df5_clean)
################################################################################
uORFs_pUP_gt30 <- OUT_df5_clean %>%
  dplyr::filter(pUP > 30) %>%
  dplyr::mutate(
    uIF = round(uORF_TL / mORF_TL, 3),
    pUP_Group = dplyr::case_when(
      pUP >= 80 ~ "80~100",
      TRUE      ~ "30~80"
    )
  ) %>%
  dplyr::select(
    gene_id, tx_id, uORF_id, uORF_id_short,
    genome_range, Strand,
    pUP, pUP_Group,
    uORF_P_sites, uORF_length, uORF_length_nt,
    uORF_TL, mORF_TL, uIF,
    method, peptide_seq,
    ORF_gstart, ORF_gstop, uORF_Tx_range,
    chr, g_min = g_min, g_max = g_max   # already present in OUT_df5_clean
  ) %>%
  dplyr::arrange(dplyr::desc(pUP))

openxlsx::write.xlsx(
  list("uORFs_pUP_gt30" = uORFs_pUP_gt30),
  file = file.path(OUT_DIR, "Table_uORFs_pUP_gt30.xlsx"),
  overwrite = TRUE
)

message(sprintf("✓ Saved Table_uORFs_pUP_gt30.xlsx  (%d uORFs with pUP > 30%%)",
                nrow(uORFs_pUP_gt30)))

################################################################################
# 4) TE table (RSEM) + one-uORF genes + uIF
################################################################################
RNA <- read.delim(PATH_RSEM_RNA, header = TRUE, sep = "\t")
Ribo <- read.delim(PATH_RSEM_RIBO, header = TRUE, sep = "\t")
fiveUTRByGeneNames <- read.delim(PATH_5UTR_GENES, header = TRUE)$gene_id

TE_tbl <- data.frame(
  gene_id = as.character(RNA$gene_id),
  RNA = as.numeric(RNA$TPM),
  Ribo = as.numeric(Ribo$TPM),
  TE = as.numeric(Ribo$TPM / RNA$TPM),
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(gene_id %in% fiveUTRByGeneNames, RNA > 1, !is.na(TE) & is.finite(TE))

gene_counts <- table(OUT_df5_clean$gene_id)
single_uorf_genes <- names(gene_counts)[gene_counts == 1]

OUT_df6 <- dplyr::left_join(OUT_df5_clean, TE_tbl, by = "gene_id") %>%
  dplyr::filter(pUP > 0) %>%
  dplyr::filter(gene_id %in% single_uorf_genes) %>%
  dplyr::mutate(
    uIF_raw = uORF_TL / mORF_TL,
    uIF = round(uIF_raw, 3)
  ) %>%
  dplyr::filter(!is.na(uIF_raw) & is.finite(uIF_raw) & uIF_raw > 0)

################################################################################
# TABLE S2 EXPORT — genes containing a single TuORF
# The 3,213-gene analysis set used for Figures 1-2 ("We first focused on 3,213
# genes with a single TuORF (Table S2)"), with uIF groups as in Figure 1.
# Writes Table_S2_single_TuORF_genes.xlsx (data sheet + column legend) to OUT_DIR.
################################################################################
TableS2 <- OUT_df6 %>%
  dplyr::mutate(
    uIF_Group = dplyr::case_when(          # same cutoffs as Figures 1E/2C
      uIF > 1   ~ "Strong (uIF > 1)",
      uIF > 0.2 ~ "Medium (0.2 < uIF <= 1)",
      TRUE      ~ "Weak (uIF <= 0.2)"
    ),
    pUP_Group = dplyr::case_when(
      pUP >= 80 ~ "80~100",
      pUP >= 30 ~ "30~80",
      TRUE      ~ "0~30"
    ),
    method = dplyr::recode(method, RT = "RiboTaper", .default = method)
  ) %>%
  dplyr::transmute(
    gene_id, tx_id,
    chr, Strand,
    uORF_genome_range  = genome_range,
    uORF_start_genomic = uORF_start_bio,
    uORF_stop_genomic  = uORF_stop_bio,
    uORF_Tx_range,
    uORF_length_aa     = uORF_length,
    uORF_length_nt,
    peptide_seq,
    detection_method   = method,
    uORF_P_sites,
    uORF_TL,
    mORF_P_sites       = mORF_Ribo_count,
    mORF_length_nt     = mORF_length,
    mORF_TL,
    uIF, uIF_Group,
    pct_uORF_plus      = pUP,
    pUP_Group,
    RNA_TPM = RNA, Ribo_TPM = Ribo, mORF_TE = TE
  ) %>%
  dplyr::arrange(dplyr::desc(uIF))

stopifnot(nrow(TableS2) == dplyr::n_distinct(TableS2$gene_id))  # one row per gene

TableS2_legend <- data.frame(
  Column = c("gene_id", "tx_id", "chr", "Strand", "uORF_genome_range",
             "uORF_start_genomic", "uORF_stop_genomic", "uORF_Tx_range",
             "uORF_length_aa", "uORF_length_nt", "peptide_seq", "detection_method",
             "uORF_P_sites", "uORF_TL", "mORF_P_sites", "mORF_length_nt", "mORF_TL",
             "uIF", "uIF_Group", "pct_uORF_plus", "pUP_Group",
             "RNA_TPM", "Ribo_TPM", "mORF_TE"),
  Description = c(
    "Arabidopsis gene identifier (AGI)",
    "Transcript isoform carrying the TuORF",
    "Chromosome",
    "Genomic strand",
    "TuORF genomic range (chr_start_end)",
    "Genomic coordinate of the TuORF start codon (biological orientation)",
    "Genomic coordinate of the TuORF stop (biological orientation)",
    "TuORF range in transcript coordinates",
    "TuORF peptide length (amino acids)",
    "TuORF length (nucleotides)",
    "Predicted TuORF peptide sequence",
    "Method that identified the TuORF (RiboTaper, CiPS, or Both)",
    "Ribo-seq P-site count within the TuORF",
    "TuORF translation level: P-sites / TuORF length (nt)",
    "Ribo-seq P-site count within the mORF CDS",
    "mORF CDS length (nucleotides)",
    "mORF translation level: P-sites / CDS length (nt)",
    "uORF Impact Factor: uORF_TL / mORF_TL (rounded to 3 decimals)",
    "uIF group as in Figure 1: Strong (uIF > 1), Medium (0.2 < uIF <= 1), Weak (uIF <= 0.2)",
    paste0("%uORF+ (percentage of transcripts including the uORF) based on ",
           ifelse(TSS_SOURCE == "CAGE", "CAGE-seq", "TSS-seq (Nielsen et al. 2019)")),
    "%uORF+ category (0~30, 30~80, 80~100)",
    "mRNA abundance (RSEM TPM, RNA-seq)",
    "Ribosome footprint abundance (RSEM TPM, Ribo-seq)",
    "mORF translation efficiency: Ribo_TPM / RNA_TPM"
  ),
  stringsAsFactors = FALSE
)

wb_s2 <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb_s2, "Table S2")
openxlsx::addWorksheet(wb_s2, "Column legend")
.hdr_s2 <- openxlsx::createStyle(textDecoration = "bold")
openxlsx::writeData(wb_s2, "Table S2", TableS2, headerStyle = .hdr_s2)
openxlsx::writeData(wb_s2, "Column legend", TableS2_legend, headerStyle = .hdr_s2)
openxlsx::freezePane(wb_s2, "Table S2", firstRow = TRUE)
openxlsx::setColWidths(wb_s2, "Table S2", cols = seq_len(ncol(TableS2)), widths = "auto")
openxlsx::setColWidths(wb_s2, "Column legend", cols = 1:2, widths = c(24, 90))
openxlsx::saveWorkbook(wb_s2, file.path(OUT_DIR, "Table_S2_single_TuORF_genes.xlsx"), overwrite = TRUE)

message(sprintf("✓ Saved Table_S2_single_TuORF_genes.xlsx  (%d single-TuORF genes; %d Strong / %d Medium / %d Weak)",
                nrow(TableS2),
                sum(TableS2$uIF_Group == "Strong (uIF > 1)"),
                sum(TableS2$uIF_Group == "Medium (0.2 < uIF <= 1)"),
                sum(TableS2$uIF_Group == "Weak (uIF <= 0.2)")))

################################################################################
# TABLE S1 EXPORT — all TuORFs
# Table S1A: all unique TuORFs with %uORF+ from CAGE-seq   (Methods; Figure 2C)
# Table S1B: the same TuORFs with %uORF+ from TSS-seq      (Figure S2)
# Both %uORF+ columns are computed in this one run (the non-active TSS source
# is recomputed locally), so the two sheets always stay in sync.
# Writes Table_S1_all_TuORFs.xlsx (S1A + S1B + column legend) to OUT_DIR.
################################################################################

## pUP for an arbitrary TSS data source — mirrors section 3 exactly
.compute_pUP_from_source <- function(cage_path) {
  raw <- read.delim(cage_path, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  colnames(raw)[1:4] <- c("count", "chr", "start", "strand")
  raw$end <- raw$start
  raw$chr <- gsub("^chr|^Chr", "", raw$chr, ignore.case = TRUE)
  gr <- GenomicRanges::makeGRangesFromDataFrame(raw, keep.extra.columns = TRUE)
  try({ GenomeInfoDb::seqlevelsStyle(gr) <- GenomeInfoDb::seqlevelsStyle(tss_windows) }, silent = TRUE)
  lv <- intersect(GenomeInfoDb::seqlevels(gr), GenomeInfoDb::seqlevels(tss_windows))
  gr <- GenomeInfoDb::keepSeqlevels(gr, lv, pruning.mode = "coarse")
  tw <- GenomeInfoDb::keepSeqlevels(tss_windows, lv, pruning.mode = "coarse")
  ov <- GenomicRanges::findOverlaps(gr, tw)
  dtc <- data.table::data.table(
    tx_id      = names(tw)[S4Vectors::subjectHits(ov)],
    cage_coord = BiocGenerics::start(gr)[S4Vectors::queryHits(ov)],
    count      = gr$count[S4Vectors::queryHits(ov)]
  )
  dtt <- data.table::data.table(
    tx_id          = mORF_uORF_Info$tx_id,
    uORF_start_bio = mORF_uORF_Info$uORF_start_bio,
    Strand         = mORF_uORF_Info$Strand,
    unique_key     = paste(mORF_uORF_Info$tx_id, mORF_uORF_Info$genome_range, sep = "|")
  )
  dj <- dtc[dtt, on = .(tx_id), allow.cartesian = TRUE]
  dj[, .(
    pUP = {
      tc <- sum(count, na.rm = TRUE)
      if (tc == 0 || is.na(tc)) 0 else {
        Cb <- if (Strand[1] == "+") sum(count[cage_coord <= uORF_start_bio], na.rm = TRUE)
        else                        sum(count[cage_coord >= uORF_start_bio], na.rm = TRUE)
        (Cb / tc) * 100
      }
    }
  ), by = .(unique_key)]
}

## active source is already in OUT_pUP; compute only the other one fresh
pUP_CAGE_tbl <- if (TSS_SOURCE == "CAGE") OUT_pUP else .compute_pUP_from_source(PATH_CAGE_CAGEseq)
pUP_TSS_tbl  <- if (TSS_SOURCE == "TSS")  OUT_pUP else .compute_pUP_from_source(PATH_CAGE_TSSseq)

TableS1_base <- OUT_df5_clean %>%
  dplyr::left_join(pUP_CAGE_tbl %>% dplyr::rename(pUP_CAGE = pUP), by = "unique_key") %>%
  dplyr::left_join(pUP_TSS_tbl  %>% dplyr::rename(pUP_TSS  = pUP), by = "unique_key") %>%
  dplyr::mutate(
    pUP_CAGE = round(dplyr::coalesce(pUP_CAGE, 0), 3),
    pUP_TSS  = round(dplyr::coalesce(pUP_TSS,  0), 3),
    uIF      = round(uORF_TL / mORF_TL, 3),
    method   = dplyr::recode(method, RT = "RiboTaper", .default = method)
  ) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::mutate(n_TuORFs_in_gene = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(TE_tbl, by = "gene_id")

.make_S1_sheet <- function(df, pup_col) {
  df %>%
    dplyr::mutate(
      pct_uORF_plus = .data[[pup_col]],
      pUP_Group = dplyr::case_when(
        pct_uORF_plus >= 80 ~ "80~100",
        pct_uORF_plus >= 30 ~ "30~80",
        TRUE                ~ "0~30"
      )
    ) %>%
    dplyr::transmute(
      gene_id, tx_id,
      chr, Strand,
      uORF_genome_range  = genome_range,
      uORF_start_genomic = uORF_start_bio,
      uORF_stop_genomic  = uORF_stop_bio,
      uORF_Tx_range,
      uORF_length_aa     = uORF_length,
      uORF_length_nt,
      peptide_seq,
      detection_method   = method,
      n_TuORFs_in_gene,
      uORF_P_sites,
      uORF_TL,
      mORF_P_sites       = mORF_Ribo_count,
      mORF_length_nt     = mORF_length,
      mORF_TL,
      uIF,
      pct_uORF_plus, pUP_Group,
      RNA_TPM = RNA, Ribo_TPM = Ribo, mORF_TE = TE
    ) %>%
    dplyr::arrange(gene_id, uORF_start_genomic)
}

TableS1A <- .make_S1_sheet(TableS1_base, "pUP_CAGE")
TableS1B <- .make_S1_sheet(TableS1_base, "pUP_TSS")

stopifnot(nrow(TableS1A) == nrow(OUT_df5_clean), nrow(TableS1B) == nrow(TableS1A))

TableS1_legend <- data.frame(
  Column = c("gene_id", "tx_id", "chr", "Strand", "uORF_genome_range",
             "uORF_start_genomic", "uORF_stop_genomic", "uORF_Tx_range",
             "uORF_length_aa", "uORF_length_nt", "peptide_seq", "detection_method",
             "n_TuORFs_in_gene", "uORF_P_sites", "uORF_TL", "mORF_P_sites",
             "mORF_length_nt", "mORF_TL", "uIF", "pct_uORF_plus", "pUP_Group",
             "RNA_TPM", "Ribo_TPM", "mORF_TE"),
  Description = c(
    "Arabidopsis gene identifier (AGI)",
    "Transcript isoform carrying the TuORF",
    "Chromosome",
    "Genomic strand",
    "TuORF genomic range (chr_start_end)",
    "Genomic coordinate of the TuORF start codon (biological orientation)",
    "Genomic coordinate of the TuORF stop (biological orientation)",
    "TuORF range in transcript coordinates",
    "TuORF peptide length (amino acids)",
    "TuORF length (nucleotides)",
    "Predicted TuORF peptide sequence",
    "Method that identified the TuORF (RiboTaper, CiPS, or Both)",
    "Number of unique TuORFs in this gene (this table lists each on its own row)",
    "Ribo-seq P-site count within the TuORF",
    "TuORF translation level: P-sites / TuORF length (nt)",
    "Ribo-seq P-site count within the mORF CDS",
    "mORF CDS length (nucleotides)",
    "mORF translation level: P-sites / CDS length (nt)",
    "uORF Impact Factor: uORF_TL / mORF_TL (rounded to 3 decimals)",
    "%uORF+ (percentage of transcripts including the uORF); sheet S1A = CAGE-seq, sheet S1B = TSS-seq (Nielsen et al. 2019)",
    "%uORF+ category (0~30, 30~80, 80~100)",
    "mRNA abundance (RSEM TPM, RNA-seq); NA if the gene is not in the TE analysis set (5'UTR-containing, RNA TPM > 1)",
    "Ribosome footprint abundance (RSEM TPM, Ribo-seq); NA as above",
    "mORF translation efficiency: Ribo_TPM / RNA_TPM; NA as above"
  ),
  stringsAsFactors = FALSE
)

wb_s1 <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb_s1, "Table S1A (CAGE)")
openxlsx::addWorksheet(wb_s1, "Table S1B (TSS-seq)")
openxlsx::addWorksheet(wb_s1, "Column legend")
.hdr_s1 <- openxlsx::createStyle(textDecoration = "bold")
openxlsx::writeData(wb_s1, "Table S1A (CAGE)", TableS1A, headerStyle = .hdr_s1)
openxlsx::writeData(wb_s1, "Table S1B (TSS-seq)", TableS1B, headerStyle = .hdr_s1)
openxlsx::writeData(wb_s1, "Column legend", TableS1_legend, headerStyle = .hdr_s1)
openxlsx::freezePane(wb_s1, "Table S1A (CAGE)", firstRow = TRUE)
openxlsx::freezePane(wb_s1, "Table S1B (TSS-seq)", firstRow = TRUE)
openxlsx::setColWidths(wb_s1, "Table S1A (CAGE)", cols = seq_len(ncol(TableS1A)), widths = "auto")
openxlsx::setColWidths(wb_s1, "Table S1B (TSS-seq)", cols = seq_len(ncol(TableS1B)), widths = "auto")
openxlsx::setColWidths(wb_s1, "Column legend", cols = 1:2, widths = c(24, 100))
openxlsx::saveWorkbook(wb_s1, file.path(OUT_DIR, "Table_S1_all_TuORFs.xlsx"), overwrite = TRUE)

message(sprintf("✓ Saved Table_S1_all_TuORFs.xlsx  (%d TuORFs from %d genes; S1A = CAGE, S1B = TSS-seq)",
                nrow(TableS1A), dplyr::n_distinct(TableS1A$gene_id)))

################################################################################
# FIGURE 1 (C, D, E) — combined panel
################################################################################
median_uIF <- stats::median(OUT_df6$uIF, na.rm = TRUE)

pC <- ggplot2::ggplot(OUT_df6, ggplot2::aes(x = uIF)) +
  ggplot2::geom_histogram(binwidth = 0.1, fill = "#3498db", color = "white", size = 0.1) +
  ggplot2::geom_vline(xintercept = median_uIF, linetype = "dashed", color = "red", size = 1) +
  ggplot2::annotate("text", x = median_uIF + 0.2, y = Inf,
                    label = paste("Median =", round(median_uIF, 2)),
                    vjust = 2, hjust = 0, color = "red", fontface = "bold") +
  ggplot2::coord_cartesian(xlim = c(0, 8)) +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "C", x = "uIF", y = "Number of genes") +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16))

OUT_df6_log2_both <- OUT_df6 %>% dplyr::filter(is.finite(uIF_raw), uIF_raw > 0, is.finite(TE), TE > 0)

pD <- ggplot2::ggplot(OUT_df6_log2_both, ggplot2::aes(x = log2(uIF_raw), y = log2(TE))) +
  ggplot2::geom_point(color = "#3498db", alpha = 0.6, size = 1.5) +
  ggplot2::geom_smooth(method = "lm", color = "red", size = 1, se = FALSE) +
  ggpubr::stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 4) +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "D", x = bquote(log[2]~"(uIF)"), y = bquote(log[2]~"(mORF TE)")) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16))

TE_NoTuORF <- TE_tbl %>%
  dplyr::filter(!(gene_id %in% unique(OUT_df5_clean$gene_id))) %>%
  dplyr::transmute(gene_id, TE, Group_Label = "No TuORF")

PlotE_Data <- dplyr::bind_rows(
  OUT_df6 %>%
    dplyr::transmute(
      gene_id, TE,
      Group_Label = dplyr::case_when(
        uIF > 1   ~ "Strong uORF (>1)",
        uIF > 0.2 ~ "Mid uORF (>0.2 - 1)",
        TRUE      ~ "Weak uORF (0 - 0.2)"
      )
    ),
  TE_NoTuORF
) %>% dplyr::filter(is.finite(TE))

PlotE_Data$Group_Label <- factor(
  PlotE_Data$Group_Label,
  levels = c("Strong uORF (>1)", "Mid uORF (>0.2 - 1)", "Weak uORF (0 - 0.2)", "No TuORF")
)

statsE <- generate_legend_labels(PlotE_Data, "Group_Label", "TE")
PlotE_Data <- PlotE_Data %>%
  dplyr::left_join(statsE %>% dplyr::select(Group_Label, Final_Label), by = "Group_Label")

PlotE_Data$Final_Label <- factor(
  PlotE_Data$Final_Label,
  levels = statsE$Final_Label[order(match(statsE$Group_Label, levels(PlotE_Data$Group_Label)))]
)

pal3 <- scales::hue_pal()(3)
colsE <- c(pal3[1], pal3[2], pal3[3], "#999999")
names(colsE) <- levels(PlotE_Data$Final_Label)

pE <- ggplot2::ggplot(PlotE_Data %>% dplyr::filter(TE <= 3.5),
                      ggplot2::aes(x = TE, color = Final_Label)) +
  ggplot2::stat_ecdf(geom = "step", size = 1.2) +
  ggplot2::scale_color_manual(values = colsE) +
  ggplot2::coord_cartesian(xlim = c(0, 3.5)) +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "E", x = "mORF TE", y = "Fraction", color = "") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 16),
    legend.position = c(0.70, 0.4),
    legend.background = ggplot2::element_blank(),
    legend.key.size = grid::unit(0.5, "cm")
  )

Figure1_CDE <- ggpubr::ggarrange(pC, pD, pE, NULL, ncol = 2, nrow = 2)
print(Figure1_CDE)
save_pdf(Figure1_CDE, "Figure1_CDE.pdf", w = 10, h = 10)


################################################################################
# FIGURE 3 (B,C,D,E,F,G,H) — combined
################################################################################
df_main <- OUT_df6 %>%
  dplyr::mutate(
    uIF_Group = dplyr::case_when(
      uIF > 1   ~ "Strong uORF (>1)",
      uIF > 0.2 ~ "Mid uORF (>0.2 - 1)",
      TRUE      ~ "Weak uORF (0 - 0.2)"
    ),
    pUP_Group = dplyr::case_when(
      pUP >= 80 ~ "80~100",
      pUP >= 30 ~ "30~80",
      TRUE      ~ "0~30"
    )
  )
df_main$uIF_Group <- factor(df_main$uIF_Group,
                            levels = c("Strong uORF (>1)", "Mid uORF (>0.2 - 1)", "Weak uORF (0 - 0.2)"))
df_main$pUP_Group <- factor(df_main$pUP_Group, levels = c("80~100", "30~80", "0~30"))

# TE for "No TuORF"
df_other <- TE_NoTuORF %>% dplyr::select(gene_id, TE) %>% dplyr::mutate(pUP_Group = "No TuORF", pUP = NA_real_)
df_3F <- dplyr::bind_rows(
  df_main %>% dplyr::select(gene_id, TE, pUP_Group),
  df_other
)
df_3F$pUP_Group <- factor(df_3F$pUP_Group, levels = c("80~100", "30~80", "0~30", "No TuORF"))

# B: pUP histogram with 30/80 lines
total_n <- nrow(df_main)
stats_3B <- df_main %>% dplyr::count(pUP_Group, name="n") %>%
  dplyr::mutate(pct = round((n / total_n) * 100, 1))
label_low  <- paste0(stats_3B$pct[stats_3B$pUP_Group=="0~30"], "%\n(n=", stats_3B$n[stats_3B$pUP_Group=="0~30"], ")")
label_mid  <- paste0(stats_3B$pct[stats_3B$pUP_Group=="30~80"], "%\n(n=", stats_3B$n[stats_3B$pUP_Group=="30~80"], ")")
label_high <- paste0(stats_3B$pct[stats_3B$pUP_Group=="80~100"], "%\n(n=", stats_3B$n[stats_3B$pUP_Group=="80~100"], ")")

p3B <- ggplot2::ggplot(df_main, ggplot2::aes(x = pUP)) +
  ggplot2::geom_histogram(binwidth = 2, fill = "#3498db", color = "white", size = 0.1) +
  ggplot2::geom_vline(xintercept = c(30, 80), color = "black", size = 0.5) +
  ggplot2::annotate("text", x = 15, y = Inf, label = label_low,  vjust = 1.2, size = 3.5) +
  ggplot2::annotate("text", x = 55, y = Inf, label = label_mid,  vjust = 1.2, size = 3.5) +
  ggplot2::annotate("text", x = 90, y = Inf, label = label_high, vjust = 1.2, size = 3.5) +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "B", x = "%pUP (CAGE)", y = "Number of genes") +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16))

# C: ECDF pUP by uIF group (legend letters)
stats_3C <- generate_legend_labels(df_main, "uIF_Group", "pUP")
df_main_3C <- df_main %>% dplyr::left_join(stats_3C %>% dplyr::select(uIF_Group, Final_Label), by = "uIF_Group")
ordered_3C <- unique(stats_3C$Final_Label[order(match(stats_3C$uIF_Group, levels(df_main$uIF_Group)))])
df_main_3C$Final_Label <- factor(df_main_3C$Final_Label, levels = ordered_3C)
palC <- scales::hue_pal()(3)
cols_3C <- setNames(palC, ordered_3C)

p3C <- ggplot2::ggplot(df_main_3C, ggplot2::aes(x = pUP, color = Final_Label)) +
  ggplot2::stat_ecdf(geom = "step", size = 1.2) +
  ggplot2::scale_color_manual(values = cols_3C) +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "C", x = "%pUP (CAGE)", y = "Fraction", color = "") +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16),
                 legend.position = c(0.70, 0.4),
                 legend.background = ggplot2::element_blank())

# D: pUP vs log2(uIF_raw)
p3D <- ggplot2::ggplot(df_main, ggplot2::aes(x = pUP, y = log2(uIF_raw))) +
  ggplot2::geom_point(color = "#3498db", alpha = 0.5, size = 1) +
  ggplot2::geom_smooth(method = "lm", color = "red", size = 1, se = FALSE) +
  ggpubr::stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "D", x = "%pUP (CAGE)", y = bquote(log[2]~"(uIF)")) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16))

# E: pUP vs log2(TE)
p3E <- ggplot2::ggplot(df_main %>% dplyr::filter(is.finite(TE), TE > 0), ggplot2::aes(x = pUP, y = log2(TE))) +
  ggplot2::geom_point(color = "#3498db", alpha = 0.5, size = 1) +
  ggplot2::geom_smooth(method = "lm", color = "red", size = 1, se = FALSE) +
  ggpubr::stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "E", x = "%pUP (CAGE)", y = bquote(log[2]~"(mORF TE)")) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16))

# F: TE ECDF by pUP bins + No TuORF (legend letters)
stats_3F <- generate_legend_labels(df_3F, "pUP_Group", "TE")
df_3F_plot <- df_3F %>% dplyr::left_join(stats_3F %>% dplyr::select(pUP_Group, Final_Label), by = "pUP_Group")
ordered_3F <- unique(stats_3F$Final_Label[order(match(stats_3F$pUP_Group, levels(df_3F$pUP_Group)))])
df_3F_plot$Final_Label <- factor(df_3F_plot$Final_Label, levels = ordered_3F)
palF <- scales::hue_pal()(3)
cols_3F <- c(setNames(palF[1], ordered_3F[1]),
             setNames(palF[2], ordered_3F[2]),
             setNames(palF[3], ordered_3F[3]),
             setNames("#999999", ordered_3F[4]))

p3F <- ggplot2::ggplot(df_3F_plot %>% dplyr::filter(TE <= 3.5),
                       ggplot2::aes(x = TE, color = Final_Label)) +
  ggplot2::stat_ecdf(geom = "step", size = 1.2) +
  ggplot2::scale_color_manual(values = cols_3F) +
  ggplot2::coord_cartesian(xlim = c(0, 3.5)) +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "F", x = "mORF TE", y = "Fraction", color = "") +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16),
                 legend.position = c(0.70, 0.4),
                 legend.background = ggplot2::element_blank())

# G/H (mORF-start pUP)
cds_all_list <- GenomicFeatures::cdsBy(txdb, by="tx", use.names=TRUE)
fiveUTR_all_list <- GenomicFeatures::fiveUTRsByTranscript(txdb, use.names=TRUE)
valid_tx_ids <- intersect(names(cds_all_list), names(fiveUTR_all_list))

all_tx_map <- GenomicFeatures::transcripts(txdb, columns=c("TXNAME", "GENEID"))
df_tx_map <- data.frame(tx_id = all_tx_map$TXNAME, gene_id = unlist(all_tx_map$GENEID)) %>%
  dplyr::filter(tx_id %in% valid_tx_ids)

target_mORF_tx <- df_tx_map %>%
  dplyr::group_by(gene_id) %>%
  dplyr::arrange(tx_id) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

cds_target <- cds_all_list[names(cds_all_list) %in% target_mORF_tx$tx_id]
cds_flat <- BiocGenerics::unlist(cds_target, use.names = FALSE)
cds_flat$tx_id <- rep(names(cds_target), times = S4Vectors::elementNROWS(cds_target))

df_cds_coords <- as.data.frame(cds_flat) %>%
  dplyr::group_by(tx_id) %>%
  dplyr::summarise(
    strand = unique(strand),
    mORF_start_bio = ifelse(unique(strand) == "+", min(start), max(end)),
    .groups = "drop"
  )

all_tx_m <- GenomicFeatures::transcripts(txdb, use.names=TRUE)
all_tx_m <- all_tx_m[names(all_tx_m) %in% target_mORF_tx$tx_id]
tss_windows_m <- suppressWarnings(GenomicRanges::promoters(all_tx_m, upstream = 100, downstream = BiocGenerics::width(all_tx_m)))
tss_windows_m <- GenomicRanges::trim(tss_windows_m)

overlaps_m <- GenomicRanges::findOverlaps(CAGE_gr, tss_windows_m)
dt_cage_m <- data.table::data.table(
  tx_id = names(tss_windows_m)[S4Vectors::subjectHits(overlaps_m)],
  cage_coord = BiocGenerics::start(CAGE_gr)[S4Vectors::queryHits(overlaps_m)],
  count = CAGE_gr$count[S4Vectors::queryHits(overlaps_m)]
)

dt_targets_m <- data.table::data.table(
  tx_id = df_cds_coords$tx_id,
  mORF_start_bio = df_cds_coords$mORF_start_bio,
  Strand = df_cds_coords$strand
)

dt_join_m <- dt_cage_m[dt_targets_m, on = .(tx_id), nomatch=0]

OUT_pUP_mORF <- dt_join_m[, .(
  pUP = {
    total_c <- sum(count, na.rm = TRUE)
    if (total_c == 0) 0 else {
      Cb <- if (Strand[1] == "+") sum(count[cage_coord <= mORF_start_bio], na.rm = TRUE)
      else                        sum(count[cage_coord >= mORF_start_bio], na.rm = TRUE)
      (Cb / total_c) * 100
    }
  }
), by = .(tx_id)]

df_plot_G <- rbind(
  data.frame(pUP = df_main$pUP, Type = "uORF Start (1-uORF genes)"),
  data.frame(pUP = OUT_pUP_mORF$pUP, Type = "mORF Start (Protein-coding)")
)
df_plot_G$Type <- factor(df_plot_G$Type, levels = c("mORF Start (Protein-coding)", "uORF Start (1-uORF genes)"))

stats_3G <- generate_legend_labels(df_plot_G, "Type", "pUP")
df_plot_G2 <- df_plot_G %>% dplyr::left_join(stats_3G %>% dplyr::select(Type, Final_Label), by="Type")
df_plot_G2$Final_Label <- factor(df_plot_G2$Final_Label, levels = stats_3G$Final_Label)

# stable two colors
col_map_G <- setNames(c("#4D4D4D", "#E76F69"), levels(df_plot_G2$Final_Label))

p3G <- ggplot2::ggplot(df_plot_G2, ggplot2::aes(x = pUP, color = Final_Label)) +
  ggplot2::stat_ecdf(geom = "step", size = 1.2) +
  ggplot2::scale_color_manual(values = col_map_G) +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "G", x = "% TSS Upstream (CAGE)", y = "Fraction", color = "") +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16),
                 legend.position = c(0.6, 0.2),
                 legend.background = ggplot2::element_blank())

# H: mORF pUP histogram
mORF_hist_stats <- OUT_pUP_mORF %>%
  dplyr::mutate(
    pUP_Group = dplyr::case_when(
      pUP >= 80 ~ "80~100",
      pUP >= 30 ~ "30~80",
      TRUE ~ "0~30"
    )
  ) %>%
  dplyr::count(pUP_Group, name="n") %>%
  dplyr::mutate(pct = round((n / sum(n)) * 100, 1))

get_morf_label <- function(grp_name) {
  row <- mORF_hist_stats[mORF_hist_stats$pUP_Group == grp_name, ]
  if (nrow(row) == 0) return("0%\n(n=0)")
  paste0(row$pct, "%\n(n=", row$n, ")")
}

label_low_m  <- get_morf_label("0~30")
label_mid_m  <- get_morf_label("30~80")
label_high_m <- get_morf_label("80~100")

p3H <- ggplot2::ggplot(OUT_pUP_mORF, ggplot2::aes(x = pUP)) +
  ggplot2::geom_histogram(binwidth = 2, fill = "#3498db", color = "white", size = 0.1) +
  ggplot2::geom_vline(xintercept = c(30, 80), color = "black", size = 0.5) +
  ggplot2::annotate("text", x = 15, y = Inf, label = label_low_m,  vjust = 1.2, size = 3.5) +
  ggplot2::annotate("text", x = 55, y = Inf, label = label_mid_m,  vjust = 1.2, size = 3.5) +
  ggplot2::annotate("text", x = 90, y = Inf, label = label_high_m, vjust = 1.2, size = 3.5) +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "H", x = "% mORF pUP (CAGE)", y = "Number of genes") +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16))

Final_Fig3 <- ggpubr::ggarrange(
  p3B, p3H,
  p3C, p3G,
  p3D, p3E,
  p3F, NULL,
  ncol = 2, nrow = 4
)
print(Final_Fig3)
save_pdf(Final_Fig3, "Figure3_BCDEFGH.pdf", w = 10, h = 18)


################################################################################
# Figure 1 &3 ADD-ON: New analyses and combined figure (updated)
# (1) Correlations: log2(uIF) vs log2(half-life), log2(pUP) vs log2(half-life)
# (2) ECDF half-lives by uIF categories and by pUP categories (with 'Other genes')
# (3) New: Scatterplots for correlations by groups (removed old ECDF faceted)
# Combined in 4x6 inches figure
################################################################################

# Load half-lives data (path declared in section 0)
hl_path <- PATH_HALFLIFE_TSV
hl_df <- read.delim(hl_path, header = TRUE, sep = "\t")
colnames(hl_df)[1] <- "gene_id"  # Standardize Gene ID to gene_id
hl_df <- hl_df %>% dplyr::select(gene_id, Mean) %>% dplyr::rename(half_life = Mean)

# Join to OUT_df6 (single uORF genes)
df_hl_uorf <- dplyr::left_join(OUT_df6, hl_df, by = "gene_id") %>%
  dplyr::filter(is.finite(half_life), half_life > 0)

# Add categories (consistent with earlier script)
df_hl_uorf <- df_hl_uorf %>%
  dplyr::mutate(
    uIF_cat = dplyr::case_when(
      uIF > 1 ~ "Strong (>1)",
      uIF > 0.2 ~ "Mid (>0.2 - 1)",
      TRUE ~ "Weak (0 - 0.2)"
    ),
    pUP_cat = dplyr::case_when(
      pUP >= 80 ~ "80~100",
      pUP >= 30 ~ "30~80",
      TRUE ~ "0~30"
    )
  )
df_hl_uorf$uIF_cat <- factor(df_hl_uorf$uIF_cat, levels = c("Weak (0 - 0.2)", "Mid (>0.2 - 1)", "Strong (>1)"))
df_hl_uorf$pUP_cat <- factor(df_hl_uorf$pUP_cat, levels = c("0~30", "30~80", "80~100"))

# Other genes: half-lives for genes not in uORF genes
uorf_genes <- unique(df_hl_uorf$gene_id)
df_hl_other <- hl_df %>%
  dplyr::filter(!gene_id %in% uorf_genes) %>%
  dplyr::mutate(uIF_cat = "Other genes", pUP_cat = "Other genes") %>%
  dplyr::filter(is.finite(half_life), half_life > 0)

# (1) Correlation plots against half-lives (only on uORF genes)
df_cor <- df_hl_uorf %>%
  dplyr::filter(is.finite(uIF), uIF > 0, is.finite(pUP), pUP > 0, is.finite(half_life), half_life > 0)

p1a <- ggplot2::ggplot(df_cor, ggplot2::aes(x = log2(uIF), y = log2(half_life))) +
  ggplot2::geom_point(color = "#3498db", alpha = 0.6, size = 1) +
  ggplot2::geom_smooth(method = "lm", color = "red", size = 1, se = FALSE) +
  ggpubr::stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 4) +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "(1a) log2(uIF) vs log2(half-life)", x = bquote(log[2]~"(uIF)"), y = bquote(log[2]~"(half-life)"))

p1b <- ggplot2::ggplot(df_cor, ggplot2::aes(x = log2(pUP), y = log2(half_life))) +
  ggplot2::geom_point(color = "#3498db", alpha = 0.6, size = 1) +
  ggplot2::geom_smooth(method = "lm", color = "red", size = 1, se = FALSE) +
  ggpubr::stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 4) +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "(1b) log2(pUP) vs log2(half-life)", x = bquote(log[2]~"(pUP)"), y = bquote(log[2]~"(half-life)"))

# (2) ECDF plots for half-lives by categories (with KS letters, including Other genes)
# By uIF categories + Other
df_2a_all <- dplyr::bind_rows(
  df_hl_uorf %>% dplyr::select(gene_id, half_life, Group = uIF_cat),
  df_hl_other %>% dplyr::select(gene_id, half_life, Group = uIF_cat)
)
df_2a_all$Group <- factor(df_2a_all$Group, levels = c("Weak (0 - 0.2)", "Mid (>0.2 - 1)", "Strong (>1)", "Other genes"))

stats_2a <- generate_legend_labels(df_2a_all, "Group", "half_life")
df_2a <- df_2a_all %>% dplyr::left_join(stats_2a %>% dplyr::select(Group, Final_Label), by = "Group")
df_2a$Final_Label <- factor(df_2a$Final_Label, levels = stats_2a$Final_Label)
pal_2a <- c(scales::hue_pal()(3), "#999999")
cols_2a <- setNames(pal_2a, levels(df_2a$Final_Label))
p2a <- ggplot2::ggplot(df_2a, ggplot2::aes(x = half_life, color = Final_Label)) +
  ggplot2::stat_ecdf(geom = "step", size = 1.2) +
  ggplot2::scale_color_manual(values = cols_2a) +
  ggplot2::coord_cartesian(xlim = c(0, max(df_2a$half_life, na.rm = TRUE) * 1.1)) +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "(2a) ECDF half-lives by uIF cat", x = "mRNA half-life", y = "Fraction", color = "") +
  ggplot2::theme(legend.position = c(0.7, 0.3), legend.background = ggplot2::element_blank())

# By pUP categories + Other
df_2b_all <- dplyr::bind_rows(
  df_hl_uorf %>% dplyr::select(gene_id, half_life, Group = pUP_cat),
  df_hl_other %>% dplyr::select(gene_id, half_life, Group = pUP_cat)
)
df_2b_all$Group <- factor(df_2b_all$Group, levels = c("0~30", "30~80", "80~100", "Other genes"))

stats_2b <- generate_legend_labels(df_2b_all, "Group", "half_life")
df_2b <- df_2b_all %>% dplyr::left_join(stats_2b %>% dplyr::select(Group, Final_Label), by = "Group")
df_2b$Final_Label <- factor(df_2b$Final_Label, levels = stats_2b$Final_Label)
pal_2b <- c(scales::hue_pal()(3), "#999999")
cols_2b <- setNames(pal_2b, levels(df_2b$Final_Label))
p2b <- ggplot2::ggplot(df_2b, ggplot2::aes(x = half_life, color = Final_Label)) +
  ggplot2::stat_ecdf(geom = "step", size = 1.2) +
  ggplot2::scale_color_manual(values = cols_2b) +
  ggplot2::coord_cartesian(xlim = c(0, max(df_2b$half_life, na.rm = TRUE) * 1.1)) +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "(2b) ECDF half-lives by pUP cat", x = "mRNA half-life", y = "Fraction", color = "") +
  ggplot2::theme(legend.position = c(0.7, 0.3), legend.background = ggplot2::element_blank())

# (3) New scatterplots: correlations by groups (no half-lives, use OUT_df6 data)
df_scatter <- OUT_df6 %>%
  dplyr::filter(is.finite(uIF), uIF > 0, is.finite(pUP), pUP > 0) %>%
  dplyr::mutate(
    log_uIF = log2(uIF),
    log_pUP = log2(pUP),
    uIF_cat = dplyr::case_when(
      uIF > 1 ~ "Strong (>1)",
      uIF > 0.2 ~ "Mid (>0.2 - 1)",
      TRUE ~ "Weak (0 - 0.2)"
    ),
    pUP_cat = dplyr::case_when(
      pUP >= 80 ~ "80~100",
      pUP >= 30 ~ "30~80",
      TRUE ~ "0~30"
    )
  )
df_scatter$uIF_cat <- factor(df_scatter$uIF_cat, levels = c("Weak (0 - 0.2)", "Mid (>0.2 - 1)", "Strong (>1)"))
df_scatter$pUP_cat <- factor(df_scatter$pUP_cat, levels = c("0~30", "30~80", "80~100"))

# Compute correlations per group
cor_uIF_groups <- df_scatter %>%
  dplyr::group_by(uIF_cat) %>%
  dplyr::summarise(
    cor = cor(log_uIF, log_pUP, method = "pearson"),
    p = cor.test(log_uIF, log_pUP)$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(label = paste0("r = ", round(cor, 2), ", p = ", format.pval(p, digits = 2)))

cor_pUP_groups <- df_scatter %>%
  dplyr::group_by(pUP_cat) %>%
  dplyr::summarise(
    cor = cor(log_uIF, log_pUP, method = "pearson"),
    p = cor.test(log_uIF, log_pUP)$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(label = paste0("r = ", round(cor, 2), ", p = ", format.pval(p, digits = 2)))

pal_uIF <- scales::hue_pal()(3)
pal_pUP <- scales::hue_pal()(3)

# New 3A: log2(uIF) vs log2(pUP), colored by uIF groups, lines per group
p3a <- ggplot2::ggplot(df_scatter, ggplot2::aes(x = log_uIF, y = log_pUP, color = uIF_cat)) +
  ggplot2::geom_point(alpha = 0.6, size = 1) +
  ggplot2::geom_smooth(method = "lm", se = FALSE, size = 1) +
  ggplot2::scale_color_manual(values = pal_uIF) +
  ggplot2::geom_text(data = cor_uIF_groups, aes(x = -Inf, y = Inf, label = label), hjust = -0.1, vjust = 1.1, inherit.aes = FALSE, size = 3) +
  ggplot2::facet_wrap(~ uIF_cat, scales = "free") +  # Facet to separate annotations
  ggplot2::coord_cartesian(ylim = c(-8, 8)) +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "(3a) log2(uIF) vs log2(pUP) by uIF groups", x = bquote(log[2]~"(uIF)"), y = bquote(log[2]~"(pUP)"), color = "uIF cat")

# New 3B: log2(pUP) vs log2(uIF), colored by pUP groups, lines per group (swapped)
p3b <- ggplot2::ggplot(df_scatter, ggplot2::aes(x = log_pUP, y = log_uIF, color = pUP_cat)) +
  ggplot2::geom_point(alpha = 0.6, size = 1) +
  ggplot2::geom_smooth(method = "lm", se = FALSE, size = 1) +
  ggplot2::scale_color_manual(values = pal_pUP) +
  ggplot2::geom_text(data = cor_pUP_groups, aes(x = -Inf, y = Inf, label = label), hjust = -0.1, vjust = 1.1, inherit.aes = FALSE, size = 3) +
  ggplot2::facet_wrap(~ pUP_cat, scales = "free") +  # Facet to separate annotations
  ggplot2::coord_cartesian(ylim = c(-10, 10)) +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "(3b) log2(pUP) vs log2(uIF) by pUP groups", x = bquote(log[2]~"(pUP)"), y = bquote(log[2]~"(uIF)"), color = "pUP cat")

# Combine into 4x6 inches figure (portrait: 3 rows, 2 columns)
fig_add <- ggpubr::ggarrange(
  p1a, p1b,
  p2a, p2b,
  p3a, p3b,
  ncol = 2, nrow = 3,
  heights = c(1, 1, 1.5)  # Extra height for faceted plots
)
print(fig_add)
save_pdf(fig_add, "Figure_1FGH_combined_updated.pdf", w = 8, h = 12)


################################################################################
# FIGURE 4 — UuORF (uuORF) MODULE (FULL FIXED, NO SIDE-EFFECT NAMES)
# FIXES:
#  - Uses PATH_GTF2 explicitly (never reuses global txdb)
#  - Panel D % labels computed robustly (base table; not affected by dplyr grouping)
#  - Panel D is GENE-level (one uuORF representative per gene), y-axis = Number of genes
#  - Panel arrangement fixed: A B / C D
#  - Uses coord_cartesian() (no data dropping like ylim())
#  - All Figure4-created objects use *_4 suffix to avoid downstream name collisions
################################################################################

message("FIGURE 4: uuORF module starting (FULL FIXED; _4 names)...")

# ---- REQUIRED OBJECTS FROM UPSTREAM ----
stopifnot(exists("PATH_GTF2"), file.exists(PATH_GTF2))
stopifnot(exists("PATH_CAGE"), file.exists(PATH_CAGE))
stopifnot(exists("PATH_FASTA_TAIR10"), file.exists(PATH_FASTA_TAIR10))
stopifnot(exists("detected_genes"))
stopifnot(exists("fiveUTRByGeneNames"))
stopifnot(exists("RNA"), exists("Ribo"))
stopifnot(exists("save_pdf"))
stopifnot(exists("normalize_chr_to_ref"))
stopifnot(exists("generate_legend_labels"))
stopifnot(exists("add_counts_labels"))
stopifnot(exists("get_wilcox_letters"))

message("Figure 4 will use PATH_GTF2: ", PATH_GTF2)

# ------------------------------------------------------------------------------
# 1) Dedicated TxDb for Figure 4 from PATH_GTF2 (do NOT reuse txdb)
# ------------------------------------------------------------------------------
txdb_uu_4 <- GenomicFeatures::makeTxDbFromGFF(
  file = PATH_GTF2, format = "gtf",
  dataSource = "Araport11",
  organism = "Arabidopsis"
)

# ------------------------------------------------------------------------------
# 2) FASTA open (NO Rsamtools::open(); use open() generic)
# ------------------------------------------------------------------------------
tryCatch({ Rsamtools::indexFa(PATH_FASTA_TAIR10) }, error = function(e) NULL)
FA_UU_4 <- Rsamtools::FaFile(PATH_FASTA_TAIR10)
if (!Rsamtools::isOpen(FA_UU_4)) open(FA_UU_4)

fa_idx_4 <- tryCatch(Rsamtools::scanFaIndex(FA_UU_4), error = function(e) NULL)
if (is.null(fa_idx_4)) stop("Figure4: scanFaIndex() failed for FASTA.")
fa_levels_4 <- GenomeInfoDb::seqlevels(GenomeInfoDb::seqinfo(fa_idx_4))

getSeq_safe_4 <- function(fa, gr) Rsamtools::getSeq(fa, gr)

# ------------------------------------------------------------------------------
# 3) Get 5'UTRs + sequences; exclude detected genes
# ------------------------------------------------------------------------------
fiveUTRByTx_all_4 <- GenomicFeatures::fiveUTRsByTranscript(txdb_uu_4, use.names = TRUE)
fiveUTRByTx_all_4 <- fiveUTRByTx_all_4[S4Vectors::elementNROWS(fiveUTRByTx_all_4) > 0]

# harmonize 5'UTR seqlevels to FASTA style if needed
utr_lv_4 <- tryCatch(GenomeInfoDb::seqlevels(BiocGenerics::unlist(fiveUTRByTx_all_4)),
                     error = function(e) character(0))
if (length(utr_lv_4) > 0 && length(fa_levels_4) > 0) {
  new_lv_4 <- normalize_chr_to_ref(utr_lv_4, fa_levels_4)
  names(new_lv_4) <- utr_lv_4
  if (any(new_lv_4 != utr_lv_4)) {
    fiveUTRByTx_all_4 <- GenomeInfoDb::renameSeqlevels(fiveUTRByTx_all_4, new_lv_4)
  }
}

fiveUTR_seqs_all_4 <- GenomicFeatures::extractTranscriptSeqs(FA_UU_4, fiveUTRByTx_all_4)

utr_gene_ids_4 <- substr(names(fiveUTR_seqs_all_4), 1, 9)
keep_mask_4 <- !(utr_gene_ids_4 %in% detected_genes)

fiveUTRByTx_uu_4  <- fiveUTRByTx_all_4[keep_mask_4]
fiveUTR_seqs_uu_4 <- fiveUTR_seqs_all_4[keep_mask_4]

# ------------------------------------------------------------------------------
# 4) Find candidate uuORFs (AUG; min 8 aa; longest)
# ------------------------------------------------------------------------------
fiveUTR_ORFs_4 <- ORFik::findMapORFs(
  fiveUTRByTx_uu_4, fiveUTR_seqs_uu_4,
  startCodon = "ATG",
  longestORF = TRUE,
  groupByTx = FALSE,
  minimumLength = 8
)
if (length(fiveUTR_ORFs_4) == 0) stop("Figure4: No uuORFs found (findMapORFs returned 0).")

# keep tx with exactly one uuORF
orf_tx_ids_4 <- sub("_.+$", "", names(fiveUTR_ORFs_4))
tx_counts_4 <- table(orf_tx_ids_4)
single_uorf_tx_4 <- names(tx_counts_4)[tx_counts_4 == 1]

fiveUTR_ORFs_1_4 <- fiveUTR_ORFs_4[orf_tx_ids_4 %in% single_uorf_tx_4]
fiveUTRByTx_1_4  <- fiveUTRByTx_uu_4[names(fiveUTRByTx_uu_4) %in% single_uorf_tx_4]

# ------------------------------------------------------------------------------
# 5) Load CAGE and normalize chr to UTR/FASTA levels
# ------------------------------------------------------------------------------
CAGE_raw2_4 <- read.delim(PATH_CAGE, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
colnames(CAGE_raw2_4)[1:4] <- c("count", "chr", "start", "strand")
CAGE_raw2_4$end <- CAGE_raw2_4$start

chr0_4 <- gsub("^chr|^Chr", "", CAGE_raw2_4$chr, ignore.case = TRUE)
chr0_4[chr0_4 %in% c("M","m","mitochondria","MT","Mt")] <- "Mt"
chr0_4[chr0_4 %in% c("C","c","chloroplast","PT","Pt")]  <- "Pt"
CAGE_raw2_4$chr <- normalize_chr_to_ref(chr0_4, fa_levels_4)

CAGE_gr2_4 <- GenomicRanges::makeGRangesFromDataFrame(CAGE_raw2_4, keep.extra.columns = TRUE)

ref_levels_utr_4 <- tryCatch(GenomeInfoDb::seqlevels(BiocGenerics::unlist(fiveUTRByTx_all_4)),
                             error = function(e) character(0))
common_lv2_4 <- intersect(GenomeInfoDb::seqlevels(CAGE_gr2_4), ref_levels_utr_4)
common_lv2_4 <- intersect(common_lv2_4, fa_levels_4)
if (length(common_lv2_4) > 0) {
  CAGE_gr2_4 <- GenomeInfoDb::keepSeqlevels(CAGE_gr2_4, common_lv2_4, pruning.mode = "coarse")
}

# ------------------------------------------------------------------------------
# 6) Compute pUP per uuORF (NO subsetByOverlaps; use findOverlaps)
# ------------------------------------------------------------------------------
EXTEND_100_4 <- function(gr) {
  st <- as.character(GenomicRanges::strand(gr)[1])
  if (st == "+") {
    i <- which.min(GenomicRanges::start(gr))
    GenomicRanges::start(gr)[i] <- max(1L, GenomicRanges::start(gr)[i] - 100L)
  } else if (st == "-") {
    i <- which.max(GenomicRanges::end(gr))
    GenomicRanges::end(gr)[i] <- GenomicRanges::end(gr)[i] + 100L
  }
  gr
}
fiveUTRByTx_ext_4 <- GenomicRanges::GRangesList(lapply(fiveUTRByTx_1_4, EXTEND_100_4))

# ---- ensure names are preserved so [[tx_id]] works reliably ----
names(fiveUTRByTx_ext_4) <- names(fiveUTRByTx_1_4)

tx_ids_uu_4 <- names(fiveUTRByTx_1_4)
results_list_4 <- vector("list", length(tx_ids_uu_4))
pb <- txtProgressBar(min = 0, max = length(tx_ids_uu_4), style = 3)

for (i in seq_along(tx_ids_uu_4)) {
  setTxtProgressBar(pb, i)
  tx_id <- tx_ids_uu_4[i]
  
  gr_tx_ext <- fiveUTRByTx_ext_4[[tx_id]]
  tx_strand <- as.character(GenomicRanges::strand(gr_tx_ext)[1])
  
  # ---- avoid regex '.' in tx_id; match literal prefix "tx_id_" ----
  idx <- which(startsWith(names(fiveUTR_ORFs_1_4), paste0(tx_id, "_")))
  if (length(idx) == 0) next
  
  uu_name <- names(fiveUTR_ORFs_1_4)[idx[1]]
  uu_gr   <- fiveUTR_ORFs_1_4[[idx[1]]]
  
  uu_gr_ord <- if (tx_strand == "-") uu_gr[order(GenomicRanges::start(uu_gr), decreasing = TRUE)]
  else                               uu_gr[order(GenomicRanges::start(uu_gr), decreasing = FALSE)]
  
  uu_start_gen <- if (tx_strand == "+") min(GenomicRanges::start(uu_gr_ord))
  else                                  max(GenomicRanges::end(uu_gr_ord))
  
  hits <- GenomicRanges::findOverlaps(CAGE_gr2_4, gr_tx_ext, ignore.strand = FALSE)
  
  if (length(hits) == 0) {
    Ca <- 0; Cb <- 0
  } else {
    cage_hits <- CAGE_gr2_4[unique(S4Vectors::queryHits(hits))]
    if (tx_strand == "+") {
      Ca <- sum(cage_hits$count[GenomicRanges::start(cage_hits) > uu_start_gen], na.rm = TRUE)
      Cb <- sum(cage_hits$count[GenomicRanges::start(cage_hits) < uu_start_gen], na.rm = TRUE)
    } else {
      Ca <- sum(cage_hits$count[GenomicRanges::start(cage_hits) < uu_start_gen], na.rm = TRUE)
      Cb <- sum(cage_hits$count[GenomicRanges::start(cage_hits) > uu_start_gen], na.rm = TRUE)
    }
  }
  
  C_total <- as.numeric(Ca + Cb)
  pUP_val <- if (C_total > 0) (Cb / C_total * 100) else 0
  
  dna_set <- getSeq_safe_4(FA_UU_4, uu_gr_ord)
  dna_cat <- paste0(as.character(dna_set), collapse = "")
  aa_seq  <- as.character(Biostrings::translate(Biostrings::DNAString(dna_cat)))
  
  results_list_4[[i]] <- data.frame(
    gene_id         = substr(tx_id, 1, 9),
    tx_id           = tx_id,
    uuORF_id        = uu_name,
    pUP             = as.numeric(pUP_val),
    C_total         = as.numeric(C_total),
    aa_seq          = aa_seq,
    uuORF_start_gen = as.integer(uu_start_gen),
    Strand          = tx_strand,
    stringsAsFactors = FALSE
  )
}
close(pb)

UuORF_tx_df_4 <- dplyr::bind_rows(results_list_4) %>%
  dplyr::mutate(
    # ---- also catch NA (not only NaN/Inf) ----
    pUP     = dplyr::if_else(is.na(pUP) | is.nan(pUP) | !is.finite(pUP), 0, pUP),
    C_total = dplyr::if_else(is.na(C_total) | is.nan(C_total) | !is.finite(C_total), 0, C_total)
  )

# ------------------------------------------------------------------------------
# 7) GENE-level uuORF: pick 1 transcript per gene (prevents double counting)
#    rule: highest C_total, then highest pUP, then smallest tx_id
# ------------------------------------------------------------------------------
UuORF_gene_df_4 <- UuORF_tx_df_4 %>%
  dplyr::group_by(gene_id) %>%
  dplyr::arrange(dplyr::desc(C_total), dplyr::desc(pUP), tx_id, .by_group = TRUE) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    UU_Group = dplyr::case_when(
      pUP >= 80 ~ "80~100",
      pUP >= 30 ~ "30~80",
      TRUE      ~ "0~30"
    ),
    UU_Group = factor(UU_Group, levels = c("0~30", "30~80", "80~100"))
  )

# ------------------------------------------------------------------------------
# TABLE S3 EXPORT — UuORF genes
# Gene-level undetected-uORF (UuORF) table behind Figure S5D-E: one transcript
# per gene (highest CAGE count, then highest pUP), %uORF+ groups as plotted.
# Writes Table_S3_UuORF_genes.xlsx (data sheet + column legend) to OUT_DIR.
# ------------------------------------------------------------------------------
TableS3 <- UuORF_gene_df_4 %>%
  dplyr::left_join(TE_tbl, by = "gene_id") %>%
  dplyr::transmute(
    gene_id, tx_id,
    uuORF_id,
    Strand,
    uuORF_start_genomic = as.integer(uuORF_start_gen),
    uuORF_length_aa     = nchar(sub("\\*+$", "", aa_seq)),
    peptide_seq         = aa_seq,
    TSS_total_count     = C_total,
    pct_uORF_plus       = round(pUP, 3),
    pUP_Group           = as.character(UU_Group),
    RNA_TPM = RNA, Ribo_TPM = Ribo, mORF_TE = TE
  ) %>%
  dplyr::arrange(dplyr::desc(pct_uORF_plus))

stopifnot(nrow(TableS3) == dplyr::n_distinct(TableS3$gene_id))  # one row per gene

TableS3_legend <- data.frame(
  Column = c("gene_id", "tx_id", "uuORF_id", "Strand", "uuORF_start_genomic",
             "uuORF_length_aa", "peptide_seq", "TSS_total_count",
             "pct_uORF_plus", "pUP_Group", "RNA_TPM", "Ribo_TPM", "mORF_TE"),
  Description = c(
    "Arabidopsis gene identifier (AGI); genes with predicted but undetected uORFs (UuORFs) — not identified as TuORF genes by RiboTaper or CiPS",
    "Representative transcript (one per gene: highest total CAGE/TSS count, then highest %uORF+, then smallest tx_id)",
    "UuORF identifier (transcript_rank from ORF finding)",
    "Genomic strand",
    "Genomic coordinate of the UuORF start codon (biological orientation)",
    "UuORF peptide length in amino acids (stop codon excluded)",
    "Predicted UuORF peptide sequence (trailing * = stop codon)",
    paste0("Total ", ifelse(TSS_SOURCE == "CAGE", "CAGE-seq", "TSS-seq"),
           " read count within the annotated 5'UTR + 100 nt upstream (0 = no detectable signal; such genes are scored as 0% uORF+)"),
    paste0("%uORF+ (percentage of transcripts including the uORF) based on ",
           ifelse(TSS_SOURCE == "CAGE", "CAGE-seq", "TSS-seq (Nielsen et al. 2019)")),
    "%uORF+ category as in Figure S5D (0~30, 30~80, 80~100)",
    "mRNA abundance (RSEM TPM, RNA-seq); NA if the gene is not in the TE analysis set (5'UTR-containing, RNA TPM > 1)",
    "Ribosome footprint abundance (RSEM TPM, Ribo-seq); NA as above",
    "mORF translation efficiency: Ribo_TPM / RNA_TPM; NA as above"
  ),
  stringsAsFactors = FALSE
)

wb_s3 <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb_s3, "Table S3")
openxlsx::addWorksheet(wb_s3, "Column legend")
.hdr_s3 <- openxlsx::createStyle(textDecoration = "bold")
openxlsx::writeData(wb_s3, "Table S3", TableS3, headerStyle = .hdr_s3)
openxlsx::writeData(wb_s3, "Column legend", TableS3_legend, headerStyle = .hdr_s3)
openxlsx::freezePane(wb_s3, "Table S3", firstRow = TRUE)
openxlsx::setColWidths(wb_s3, "Table S3", cols = seq_len(ncol(TableS3)), widths = "auto")
openxlsx::setColWidths(wb_s3, "Column legend", cols = 1:2, widths = c(24, 110))
openxlsx::saveWorkbook(wb_s3, file.path(OUT_DIR, "Table_S3_UuORF_genes.xlsx"), overwrite = TRUE)

.s3_grp <- table(TableS3$pUP_Group)
message(sprintf("✓ Saved Table_S3_UuORF_genes.xlsx  (%d UuORF genes; 0~30: %d / 30~80: %d / 80~100: %d)",
                nrow(TableS3), .s3_grp[["0~30"]], .s3_grp[["30~80"]], .s3_grp[["80~100"]]))

# ------------------------------------------------------------------------------
# FIGURE S5D — %uORF+ distribution of UuORF genes, styled to match Figure 2B
# Two versions of the same data are written:
#   FigureS5D_UuORF_pctuORFplus_histogram_Fig2Bstyle.pdf  (2B-style histogram)
#   FigureS5D_UuORF_pctuORFplus_ECDF.pdf                  (cumulative version)
# Requires: UuORF_gene_df_4 (gene-level UuORF table), OUT_DIR, save_pdf()
# ------------------------------------------------------------------------------
df_S5D <- UuORF_gene_df_4 %>%
  dplyr::transmute(
    gene_id,
    pUP  = as.numeric(pUP),
    is100 = pUP >= 100          # genes whose uORF is in 100% of transcripts
  )

n_S5D    <- nrow(df_S5D)
n_low    <- sum(df_S5D$pUP <  30)
n_mid    <- sum(df_S5D$pUP >= 30 & df_S5D$pUP < 80)
n_high   <- sum(df_S5D$pUP >= 80)
n_100    <- sum(df_S5D$is100)

lab_low  <- paste0(round(100 * n_low  / n_S5D, 1), "%\n(n=", n_low,  ")")
lab_mid  <- paste0(round(100 * n_mid  / n_S5D, 1), "%\n(n=", n_mid,  ")")
lab_high <- paste0(round(100 * n_high / n_S5D, 1), "%\n(n=", n_high, ")")

## --- histogram (Figure 2B style: binwidth 2, 30/80 dividers, % + n labels) ----
BINW_S5D <- 2
.brk  <- seq(0, 100, by = BINW_S5D)
.cnt  <- graphics::hist(df_S5D$pUP, breaks = .brk, plot = FALSE)$counts
Y_CAP <- 250                     # the 0-2% bin (n = .cnt[1]) is clipped; labelled below
tallest <- max(.cnt)
bar_100 <- .cnt[length(.cnt)]    # height of the last bin (blue 98-<100 + red 100)

p_S5D_hist <- ggplot2::ggplot(df_S5D, ggplot2::aes(x = pUP, fill = is100)) +
  ggplot2::geom_histogram(binwidth = BINW_S5D, boundary = 0, color = "white", size = 0.1) +
  ggplot2::scale_fill_manual(values = c(`FALSE` = "#3498db", `TRUE` = "red"), guide = "none") +
  ggplot2::geom_vline(xintercept = c(30, 80), color = "black", size = 0.5) +
  ggplot2::annotate("text", x = 15, y = Inf, label = lab_low,  vjust = 1.2, size = 3.5) +
  ggplot2::annotate("text", x = 55, y = Inf, label = lab_mid,  vjust = 1.2, size = 3.5) +
  ggplot2::annotate("text", x = 90, y = Inf, label = lab_high, vjust = 1.2, size = 3.5) +
  ## clipped first bar: state its true height so the axis break is not misleading
  ggplot2::annotate("segment", x = 6.5, xend = 1.6, y = Y_CAP * 0.80, yend = Y_CAP * 0.985,
                    arrow = grid::arrow(length = grid::unit(0.14, "cm"), type = "closed"),
                    color = "#1F6FA9", linewidth = 0.35) +
  ggplot2::annotate("text", x = 7, y = Y_CAP * 0.80,
                    label = paste0("bar clipped\n(n=", tallest, ")"),
                    hjust = 0, vjust = 0.5, size = 2.8, fontface = "bold", color = "#1F6FA9") +
  ## black arrow marking the red 100% uORF+ bar (as in Figure 2B)
  ggplot2::annotate("segment", x = 100, xend = 100,
                    y = Y_CAP * 0.44, yend = bar_100 + Y_CAP * 0.03,
                    arrow = grid::arrow(length = grid::unit(0.18, "cm"), type = "closed"),
                    color = "black", linewidth = 0.4) +
  ggplot2::annotate("text", x = 98, y = Y_CAP * 0.47,
                    label = paste0("100% uORF+\n(n=", n_100, "; ",
                                   round(100 * n_100 / n_S5D, 1), "%)"),
                    hjust = 1, vjust = 0, size = 3, color = "red") +
  ggplot2::scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  ggplot2::coord_cartesian(xlim = c(-1.5, 101.5), ylim = c(0, Y_CAP)) +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "D", x = "%uORF+ (CAGE)", y = "Number of genes") +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16))

## --- ECDF of the same data (no binning, so nothing has to be clipped) --------
cum_30  <- n_low / n_S5D
cum_80  <- (n_low + n_mid) / n_S5D

p_S5D_ecdf <- ggplot2::ggplot(df_S5D, ggplot2::aes(x = pUP)) +
  ggplot2::stat_ecdf(geom = "step", linewidth = 1.1, color = "#3498db") +
  ggplot2::geom_vline(xintercept = c(30, 80), color = "black", size = 0.5) +
  ggplot2::annotate("segment", x = 0, xend = 30, y = cum_30, yend = cum_30,
                    linetype = "dashed", color = "grey40", linewidth = 0.4) +
  ggplot2::annotate("segment", x = 0, xend = 80, y = cum_80, yend = cum_80,
                    linetype = "dashed", color = "grey40", linewidth = 0.4) +
  ggplot2::annotate("point", x = c(30, 80), y = c(cum_30, cum_80), size = 1.8, color = "black") +
  ## cumulative values written against the y-axis (the top labels already give n)
  ggplot2::annotate("text", x = 1.5, y = cum_30 - 0.045, hjust = 0, size = 2.9, color = "grey25",
                    label = sprintf("%.3f", cum_30)) +
  ggplot2::annotate("text", x = 1.5, y = cum_80 + 0.035, hjust = 0, size = 2.9, color = "grey25",
                    label = sprintf("%.3f", cum_80)) +
  ggplot2::annotate("text", x = 15, y = Inf, label = lab_low,  vjust = 1.2, size = 3.5) +
  ggplot2::annotate("text", x = 55, y = Inf, label = lab_mid,  vjust = 1.2, size = 3.5) +
  ggplot2::annotate("text", x = 90, y = Inf, label = lab_high, vjust = 1.2, size = 3.5) +
  ## red marker for the 100% uORF+ genes (the final jump of the curve)
  ggplot2::annotate("segment", x = 100, xend = 100, y = 0.76, yend = 0.985,
                    arrow = grid::arrow(length = grid::unit(0.16, "cm"), type = "closed"),
                    color = "red", linewidth = 0.4) +
  ggplot2::annotate("text", x = 98, y = 0.73, hjust = 1, vjust = 1, size = 3, color = "red",
                    label = paste0("100% uORF+\n(n=", n_100, "; ",
                                   round(100 * n_100 / n_S5D, 1), "%)")) +
  ggplot2::scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  ggplot2::scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  ggplot2::coord_cartesian(xlim = c(0, 100), ylim = c(0, 1.10)) +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "D", x = "%uORF+ (CAGE)",
                y = "Cumulative fraction of UuORF genes") +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16))

save_pdf(p_S5D_hist, "FigureS5D_UuORF_pctuORFplus_histogram_Fig2Bstyle.pdf", w = 5, h = 4)
save_pdf(p_S5D_ecdf, "FigureS5D_UuORF_pctuORFplus_ECDF.pdf",                 w = 5, h = 4)

message(sprintf(
  "Figure S5D (Fig 2B style) saved: %d UuORF genes | 0~30 %d (%.1f%%) / 30~80 %d (%.1f%%) / 80~100 %d (%.1f%%) | 100%%: %d | tallest bin %d (clipped at %d)",
  n_S5D, n_low, 100*n_low/n_S5D, n_mid, 100*n_mid/n_S5D, n_high, 100*n_high/n_S5D,
  n_100, tallest, Y_CAP))

# ------------------------------------------------------------------------------
# 8) Panel D: Histogram + CORRECT % labels (ROBUST; NOT affected by grouping)
# ------------------------------------------------------------------------------
UU_levels_4 <- c("0~30", "30~80", "80~100")
uu_factor_4 <- factor(as.character(UuORF_gene_df_4$UU_Group), levels = UU_levels_4)
n_vec_4 <- as.integer(table(uu_factor_4))
total_n_4 <- sum(n_vec_4)
pct_vec_4 <- if (total_n_4 > 0) round(n_vec_4 / total_n_4 * 100, 1) else rep(0, length(n_vec_4))

hist_stats_4 <- data.frame(
  UU_Group = UU_levels_4,
  n = n_vec_4,
  pct = pct_vec_4,
  stringsAsFactors = FALSE
)

label_low_4  <- paste0(hist_stats_4$pct[1], "%\n(n=", hist_stats_4$n[1], ")")
label_mid_4  <- paste0(hist_stats_4$pct[2], "%\n(n=", hist_stats_4$n[2], ")")
label_high_4 <- paste0(hist_stats_4$pct[3], "%\n(n=", hist_stats_4$n[3], ")")

p4D_4 <- ggplot2::ggplot(UuORF_gene_df_4, ggplot2::aes(x = pUP)) +
  ggplot2::geom_histogram(binwidth = 1, fill = "#3498db", color = "white", size = 0.1) +
  ggplot2::geom_vline(xintercept = c(30, 80), color = "black", size = 0.6) +
  ggplot2::annotate("text", x = 15, y = Inf, label = label_low_4,  vjust = 1.2, size = 5) +
  ggplot2::annotate("text", x = 55, y = Inf, label = label_mid_4,  vjust = 1.2, size = 5) +
  ggplot2::annotate("text", x = 90, y = Inf, label = label_high_4, vjust = 1.2, size = 5) +
  ggplot2::scale_x_continuous(breaks = seq(0, 100, by = 20), limits = c(0, 100)) +
  ggplot2::coord_cartesian(ylim = c(0, 600)) +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "D", x = "%pUP (CAGE)", y = "Number of genes") +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16))

# ------------------------------------------------------------------------------
# 9) Panel A: ECDF pUP (GENE-level; KS letters)
# ------------------------------------------------------------------------------
statsA_4 <- generate_legend_labels(UuORF_gene_df_4, "UU_Group", "pUP")
UuORF_gene_df_A_4 <- UuORF_gene_df_4 %>%
  dplyr::left_join(statsA_4 %>% dplyr::select(UU_Group, Final_Label), by = "UU_Group")
UuORF_gene_df_A_4$Final_Label <- factor(UuORF_gene_df_A_4$Final_Label, levels = statsA_4$Final_Label)
colsA_4 <- setNames(scales::hue_pal()(length(statsA_4$Final_Label)), statsA_4$Final_Label)

p4A_4 <- ggplot2::ggplot(UuORF_gene_df_A_4, ggplot2::aes(x = pUP, color = Final_Label)) +
  ggplot2::stat_ecdf(geom = "step", linewidth = 1.1) +
  ggplot2::scale_color_manual(values = colsA_4) +
  ggplot2::coord_cartesian(xlim = c(0, 100)) +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "A", x = "pUP (%)", y = "ECDF", color = "") +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16),
                 legend.position = "bottom")

# ------------------------------------------------------------------------------
# 10) Panels B/C: TE & RNA boxplots (gene lists from GENE-level uuORF)
# ------------------------------------------------------------------------------
TE_tbl_0p1_4 <- data.frame(
  gene_id = as.character(RNA$gene_id),
  RNA     = as.numeric(RNA$TPM),
  Ribo    = as.numeric(Ribo$TPM),
  TE      = as.numeric(Ribo$TPM / RNA$TPM),
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(gene_id %in% fiveUTRByGeneNames) %>%
  dplyr::filter(is.finite(RNA), RNA > 0.1) %>%
  dplyr::filter(is.finite(TE))

uu_0_30_4   <- unique(UuORF_gene_df_4$gene_id[UuORF_gene_df_4$UU_Group == "0~30"])
uu_30_80_4  <- unique(UuORF_gene_df_4$gene_id[UuORF_gene_df_4$UU_Group == "30~80"])
uu_80_100_4 <- unique(UuORF_gene_df_4$gene_id[UuORF_gene_df_4$UU_Group == "80~100"])

TE_Detected_4 <- TE_tbl_0p1_4 %>% dplyr::filter(gene_id %in% detected_genes) %>% dplyr::mutate(Group = "Detected uORFs")
TE_0_30_4     <- TE_tbl_0p1_4 %>% dplyr::filter(gene_id %in% uu_0_30_4)       %>% dplyr::mutate(Group = "0~30")
TE_30_80_4    <- TE_tbl_0p1_4 %>% dplyr::filter(gene_id %in% uu_30_80_4)      %>% dplyr::mutate(Group = "30~80")
TE_80_100_4   <- TE_tbl_0p1_4 %>% dplyr::filter(gene_id %in% uu_80_100_4)     %>% dplyr::mutate(Group = "80~100")
TE_Other_4    <- TE_tbl_0p1_4 %>% dplyr::filter(!gene_id %in% c(UuORF_gene_df_4$gene_id, detected_genes)) %>% dplyr::mutate(Group = "Other")

group_order4_4 <- c("Detected uORFs", "0~30", "30~80", "80~100", "Other")
TE_Uu_Final_4 <- dplyr::bind_rows(TE_Detected_4, TE_0_30_4, TE_30_80_4, TE_80_100_4, TE_Other_4) %>%
  dplyr::mutate(Group = factor(Group, levels = group_order4_4))

lab_info_4 <- add_counts_labels(TE_Uu_Final_4, "Group", group_order4_4)
countsBC_4 <- lab_info_4$counts
label_map_4 <- lab_info_4$label_map

TE_Uu_Final_4 <- TE_Uu_Final_4 %>%
  dplyr::mutate(Group_Label = label_map_4[as.character(Group)]) %>%
  dplyr::mutate(Group_Label = factor(Group_Label, levels = countsBC_4$Group_Label))

letters_TE_4  <- get_wilcox_letters(TE_Uu_Final_4, "TE",  "Group", p_adjust = "BY")
letters_RNA_4 <- get_wilcox_letters(TE_Uu_Final_4, "RNA", "Group", p_adjust = "BY")

if (!is.null(letters_TE_4)) {
  letters_TE_4$Group_Label <- label_map_4[letters_TE_4$Group]
  letters_TE_4$Group_Label <- factor(letters_TE_4$Group_Label, levels = countsBC_4$Group_Label)
}
if (!is.null(letters_RNA_4)) {
  letters_RNA_4$Group_Label <- label_map_4[letters_RNA_4$Group]
  letters_RNA_4$Group_Label <- factor(letters_RNA_4$Group_Label, levels = countsBC_4$Group_Label)
}

p4B_4 <- ggpubr::ggboxplot(
  TE_Uu_Final_4, x = "Group_Label", y = "TE", color = "Group_Label",
  outlier.shape = NA, ylim = c(0, 3.2)
) +
  ggpubr::rotate_x_text(45) +
  ggplot2::theme_classic() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::labs(title = "B", x = NULL, y = "mORF TE") +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16))

if (!is.null(letters_TE_4)) {
  p4B_4 <- p4B_4 + ggplot2::geom_text(
    data = letters_TE_4,
    ggplot2::aes(x = Group_Label, y = 3.05, label = Letters),
    inherit.aes = FALSE, size = 5, fontface = "bold"
  )
}

p4C_4 <- ggpubr::ggboxplot(
  TE_Uu_Final_4, x = "Group_Label", y = "RNA", color = "Group_Label",
  outlier.shape = NA, ylim = c(0, 50)
) +
  ggpubr::rotate_x_text(45) +
  ggplot2::theme_classic() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::labs(title = "C", x = NULL, y = "RNA (TPM)") +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16))

if (!is.null(letters_RNA_4)) {
  p4C_4 <- p4C_4 + ggplot2::geom_text(
    data = letters_RNA_4,
    ggplot2::aes(x = Group_Label, y = 47, label = Letters),
    inherit.aes = FALSE, size = 5, fontface = "bold"
  )
}

# ------------------------------------------------------------------------------
# 11) Assemble Figure 4 (FIXED ORDER): A B / C D
# ------------------------------------------------------------------------------
Figure4_UuORF_4 <- ggpubr::ggarrange(
  p4A_4, p4B_4,
  p4C_4, p4D_4,
  ncol = 2, nrow = 2,
  align = "hv"
)

print(Figure4_UuORF_4)
save_pdf(Figure4_UuORF_4, "Figure4_UuORF.pdf", w = 12, h = 9)

openxlsx::write.xlsx(
  list(
    "UuORF_TranscriptLevel" = UuORF_tx_df_4,
    "UuORF_GeneLevel"       = UuORF_gene_df_4,
    "UuORF_pUP_HistStats"   = hist_stats_4,
    "UuORF_TE_Data"         = TE_Uu_Final_4,
    "Letters_TE"            = if (is.null(letters_TE_4)) data.frame() else letters_TE_4,
    "Letters_RNA"           = if (is.null(letters_RNA_4)) data.frame() else letters_RNA_4,
    "LegendLetters_A_ECDF"  = statsA_4
  ),
  file = file.path(OUT_DIR, "Figure4_UuORF_Data.xlsx"),
  overwrite = TRUE
)

try(close(FA_UU_4), silent = TRUE)

message("Figure 4 complete (_4 names): saved Figure4_UuORF.pdf and Figure4_UuORF_Data.xlsx in OUT_DIR = ", OUT_DIR)
################################################################################
################################################################################
# FIGURE 5 — CPuORF vs TuORF vs No TuORF (combined D/E/F)
# (Uses OUT_df5_clean evidence + manual CP list; robust peptide/length matching)
################################################################################
message("FIGURE 5: CPuORF module starting...")

CP_raw <- tryCatch(openxlsx::read.xlsx(PATH_CPUORF_MANUAL, sheet = 1),
                   error = function(e) openxlsx::read.xlsx(PATH_CPUORF_MANUAL))
CP_ref <- CP_raw %>%
  dplyr::mutate(
    gene_id   = trimws(as.character(gene_id)),
    pep_clean = clean_peptide(peptide),
    pep_len   = nchar(pep_clean)
  ) %>%
  dplyr::filter(!is.na(gene_id))

keys_strict <- paste0(CP_ref$gene_id, "_", CP_ref$pep_clean)
keys_length <- paste0(CP_ref$gene_id, "_", CP_ref$pep_len)

# use detected-uORF evidence with pUP already computed
df_uorf_all <- OUT_df5_clean %>%
  dplyr::mutate(
    pep_clean = clean_peptide(peptide_seq),
    pep_len   = nchar(pep_clean),
    key_strict = paste0(trimws(gene_id), "_", pep_clean),
    key_length = paste0(trimws(gene_id), "_", pep_len),
    ORF_Group  = ifelse(key_strict %in% keys_strict | key_length %in% keys_length, "CPuORF", "TuORF"),
    UI         = uORF_TL / mORF_TL
  )

# one representative per gene: prioritize CPuORF then higher TL
df_uorf_one <- df_uorf_all %>%
  dplyr::group_by(gene_id) %>%
  dplyr::arrange(dplyr::desc(ORF_Group == "CPuORF"), dplyr::desc(uORF_TL), dplyr::desc(uORF_P_sites), .by_group = TRUE) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

Plot_Data_uORF <- dplyr::inner_join(df_uorf_one %>% dplyr::select(gene_id, ORF_Group, pUP, UI),
                                    TE_tbl %>% dplyr::select(gene_id, TE),
                                    by = "gene_id")

Plot_Data_uORF$ORF_Group <- factor(Plot_Data_uORF$ORF_Group, levels = c("CPuORF", "TuORF"))

Plot_Data_NoTuORF <- TE_tbl %>%
  dplyr::filter(!(gene_id %in% unique(Plot_Data_uORF$gene_id))) %>%
  dplyr::transmute(gene_id, TE, ORF_Group = "No TuORF")

Plot_Data_All <- dplyr::bind_rows(
  Plot_Data_uORF %>% dplyr::select(gene_id, TE, ORF_Group),
  Plot_Data_NoTuORF
)
Plot_Data_All$ORF_Group <- factor(Plot_Data_All$ORF_Group, levels = c("CPuORF", "TuORF", "No TuORF"))

base_cols <- c("CPuORF"="#E76F69","TuORF"="#48B5C4","No TuORF"="#4D4D4D")
theme_clean <- ggplot2::theme_classic(base_size=14) +
  ggplot2::theme(plot.title = ggplot2::element_text(face="bold", size=16))

# D: ECDF pUP (CPuORF vs TuORF) with letters in legend
Plot_Data_uORF2 <- Plot_Data_uORF %>% dplyr::filter(is.finite(pUP))
Plot_Data_uORF2$ORF_Group <- factor(Plot_Data_uORF2$ORF_Group, levels=c("CPuORF","TuORF"))

stats_5D <- generate_legend_labels(Plot_Data_uORF2, "ORF_Group", "pUP")
Data_5D <- Plot_Data_uORF2 %>% dplyr::left_join(stats_5D %>% dplyr::select(ORF_Group, Final_Label), by="ORF_Group")
Data_5D$Final_Label <- factor(Data_5D$Final_Label, levels = stats_5D$Final_Label)

cols_5D <- setNames(base_cols[as.character(stats_5D$ORF_Group)], stats_5D$Final_Label)

p5D <- ggplot2::ggplot(Data_5D, ggplot2::aes(x=pUP, color=Final_Label)) +
  ggplot2::stat_ecdf(geom="step", size=1.2) +
  ggplot2::scale_color_manual(values=cols_5D) +
  ggplot2::coord_cartesian(xlim=c(0,100)) +
  theme_clean +
  ggplot2::labs(title="D", x="% pUP (CAGE)", y="Fraction", color="") +
  ggplot2::theme(legend.position=c(0.7,0.3), legend.background=ggplot2::element_blank())

# E: ECDF UI (CPuORF vs TuORF) with letters in legend
Plot_Data_uORF3 <- Plot_Data_uORF %>% dplyr::filter(is.finite(UI), UI <= 20)
Plot_Data_uORF3$ORF_Group <- factor(Plot_Data_uORF3$ORF_Group, levels=c("CPuORF","TuORF"))

stats_5E <- generate_legend_labels(Plot_Data_uORF3, "ORF_Group", "UI")
Data_5E <- Plot_Data_uORF3 %>% dplyr::left_join(stats_5E %>% dplyr::select(ORF_Group, Final_Label), by="ORF_Group")
Data_5E$Final_Label <- factor(Data_5E$Final_Label, levels = stats_5E$Final_Label)
cols_5E <- setNames(base_cols[as.character(stats_5E$ORF_Group)], stats_5E$Final_Label)

p5E <- ggplot2::ggplot(Data_5E, ggplot2::aes(x=UI, color=Final_Label)) +
  ggplot2::stat_ecdf(geom="step", size=1.2) +
  ggplot2::scale_color_manual(values=cols_5E) +
  ggplot2::coord_cartesian(xlim=c(0,20)) +
  theme_clean +
  ggplot2::labs(title="E", x="UI (uIF)", y="Fraction", color="") +
  ggplot2::theme(legend.position=c(0.7,0.3), legend.background=ggplot2::element_blank())

# F: ECDF TE (CPuORF/TuORF/NoTuORF) with letters in legend
Plot_Data_All2 <- Plot_Data_All %>% dplyr::filter(is.finite(TE), TE <= 4)
Plot_Data_All2$ORF_Group <- factor(Plot_Data_All2$ORF_Group, levels=c("CPuORF","TuORF","No TuORF"))

stats_5F <- generate_legend_labels(Plot_Data_All2, "ORF_Group", "TE")
Data_5F <- Plot_Data_All2 %>% dplyr::left_join(stats_5F %>% dplyr::select(ORF_Group, Final_Label), by="ORF_Group")
label_order_5F <- stats_5F$Final_Label[match(c("CPuORF","TuORF","No TuORF"), stats_5F$ORF_Group)]
Data_5F$Final_Label <- factor(Data_5F$Final_Label, levels=label_order_5F)

cols_5F <- setNames(base_cols[as.character(stats_5F$ORF_Group)], stats_5F$Final_Label)

p5F <- ggplot2::ggplot(Data_5F, ggplot2::aes(x=TE, color=Final_Label)) +
  ggplot2::stat_ecdf(geom="step", size=1.2) +
  ggplot2::scale_color_manual(values=cols_5F) +
  ggplot2::coord_cartesian(xlim=c(0,4)) +
  theme_clean +
  ggplot2::labs(title="F", x="mORF TE", y="Fraction", color="") +
  ggplot2::theme(legend.position=c(0.7,0.3), legend.background=ggplot2::element_blank())

Figure5 <- ggpubr::ggarrange(p5D, p5E, p5F, ncol=3, nrow=1)
print(Figure5)
save_pdf(Figure5, "Figure5_Robust_Matching_Stats.pdf", w=14, h=5)

# Close FASTA
try(close(FA_UU), silent = TRUE)

message("DONE. Outputs in: ", OUT_DIR)






################################################################################
# ======================== FULL FIXED CODE AFTER FIGURE 5 ======================
# Goals:
# - Consistent globals: PATH_*, OUT_DIR, save_pdf(), get_gene_id_from_tx()
# - No re-loading duplicate libraries (use pkg::fn for extra packages)
# - No re-defining helpers already defined above (unless missing)
# - FASTA open() generic (no Rsamtools::open())
################################################################################

# ------------------------------------------------------------------------------
# 0) Minimal sanity checks (objects from Figures 1–5 must already exist)
# ------------------------------------------------------------------------------
stopifnot(exists("OUT_DIR"))
stopifnot(exists("save_pdf"))
stopifnot(exists("get_gene_id_from_tx"))

stopifnot(exists("PATH_CAGE"), file.exists(PATH_CAGE))
stopifnot(exists("PATH_GTF"),  file.exists(PATH_GTF))
stopifnot(exists("PATH_FASTA_TAIR10"), file.exists(PATH_FASTA_TAIR10))

# txdb and OUT_df6 are created upstream (sections 1 and 4)
stopifnot(exists("txdb"))
stopifnot(exists("OUT_df6"))

# ------------------------------------------------------------------------------
# Helper: ensure a package exists (WITHOUT attaching)
# ------------------------------------------------------------------------------
.require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package not installed: ", pkg, "\nInstall it first (e.g., BiocManager::install or install.packages).")
  }
  invisible(TRUE)
}

################################################################################
# FIGURE 6. AUG frequency and plots
# (Chi-square + Cramer's V for all pairs) + FASTA-safe getSeq
################################################################################

message("FIGURE 6: starting...")

.require_pkg("ggseqlogo")
.require_pkg("forcats")

# ---- TuORF vs NoTuORF gene set (from OUT_df6) ----
uORF_gene_ids <- character(0)
if ("gene_id" %in% colnames(OUT_df6)) {
  uORF_gene_ids <- unique(toupper(as.character(OUT_df6$gene_id)))
} else if ("tx_id" %in% colnames(OUT_df6)) {
  uORF_gene_ids <- unique(toupper(substr(as.character(OUT_df6$tx_id), 1, 9)))
}

# ---- FASTA handle (NO Rsamtools::open(); use open() generic) ----
tryCatch({ Rsamtools::indexFa(PATH_FASTA_TAIR10) }, error = function(e) {
  message("NOTE: indexFa() failed (ignored): ", conditionMessage(e))
})

FA_6 <- Rsamtools::FaFile(PATH_FASTA_TAIR10)
if (!Rsamtools::isOpen(FA_6)) open(FA_6)

idx_6 <- tryCatch(Rsamtools::scanFaIndex(FA_6), error = function(e) NULL)
if (is.null(idx_6)) stop("scanFaIndex() failed for FASTA (cannot proceed).")

si_6 <- GenomeInfoDb::seqinfo(idx_6)
fasta_levels_6  <- GenomeInfoDb::seqlevels(si_6)
fasta_seqlens_6 <- GenomeInfoDb::seqlengths(si_6)

# FASTA-safe getSeq: drop transcripts that error (no crash)
safe_getSeq_drop_bad <- function(fa, gr) {
  seqs_try <- tryCatch(Rsamtools::getSeq(fa, gr), error = function(e) NULL)
  if (!is.null(seqs_try)) {
    names(seqs_try) <- names(gr)
    return(list(seqs = seqs_try, ok_names = names(gr)))
  }
  
  message("WARNING: batch getSeq() failed; using per-transcript getSeq() and dropping failures.")
  tx_names <- names(gr)
  
  char_seqs <- vapply(seq_along(gr), function(i) {
    tryCatch({ as.character(Rsamtools::getSeq(fa, gr[i]))[1] },
             error = function(e) NA_character_)
  }, character(1))
  
  ok <- !is.na(char_seqs) & nchar(char_seqs) > 0
  ok_names <- tx_names[ok]
  if (sum(!ok) > 0) message("Dropped ", sum(!ok), " transcripts due to FASTA retrieval errors (ignored).")
  
  seqs <- Biostrings::DNAStringSet(char_seqs[ok])
  names(seqs) <- ok_names
  list(seqs = seqs, ok_names = ok_names)
}

# Chi-square + Cramer's V helpers
cramers_v <- function(chi2, n, r, c) {
  denom <- n * min(r - 1, c - 1)
  if (denom <= 0) return(NA_real_)
  sqrt(as.numeric(chi2) / denom)
}

chisq_one <- function(tbl) {
  ct <- suppressWarnings(stats::chisq.test(tbl, correct = FALSE))
  n <- sum(tbl)
  v <- cramers_v(ct$statistic, n, nrow(tbl), ncol(tbl))
  list(chi2 = unname(ct$statistic), df = unname(ct$parameter), p = unname(ct$p.value), V = v, n = n)
}

# ------------------------------------------------------------------------------
# 6.0 LOAD CAGE + BUILD 1-EXON 5'UTR MODEL (reuses the txdb from section 1)
# ------------------------------------------------------------------------------
# CAGE as GRanges (consistent with earlier code: remove chr prefix)
CAGE_raw6 <- read.delim(PATH_CAGE, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
colnames(CAGE_raw6) <- c("count", "chr", "start", "strand")
CAGE_raw6$end <- CAGE_raw6$start

CAGE_raw6$chr <- gsub("^chr|^Chr", "", CAGE_raw6$chr, ignore.case = TRUE)
CAGE_raw6$chr[CAGE_raw6$chr %in% c("M","m","mitochondria","MT","Mt")] <- "Mt"
CAGE_raw6$chr[CAGE_raw6$chr %in% c("C","c","chloroplast","PT","Pt")]  <- "Pt"

CAGE_gr6 <- GenomicRanges::makeGRangesFromDataFrame(CAGE_raw6, keep.extra.columns = TRUE)

# 1-exon 5'UTRs (from txdb created in Fig1–5)
fiveUTR_6 <- GenomicFeatures::fiveUTRsByTranscript(txdb, use.names = TRUE)
fiveUTR_1exon <- fiveUTR_6[S4Vectors::elementNROWS(fiveUTR_6) == 1]

gr_flat <- BiocGenerics::unlist(fiveUTR_1exon)
expected_widths <- IRanges::width(gr_flat) + 100L
names(expected_widths) <- names(gr_flat)

# extend 100 nt upstream of TSS side (strand-aware)
BiocGenerics::start(gr_flat[BiocGenerics::strand(gr_flat) == "+"]) <- BiocGenerics::start(gr_flat[BiocGenerics::strand(gr_flat) == "+"]) - 100L
BiocGenerics::end(gr_flat[BiocGenerics::strand(gr_flat) == "-"])   <- BiocGenerics::end(gr_flat[BiocGenerics::strand(gr_flat) == "-"]) + 100L

# Harmonize styles + seqlevels intersection with FASTA
try({ GenomeInfoDb::seqlevelsStyle(gr_flat)  <- "NCBI" }, silent = TRUE)
try({ GenomeInfoDb::seqlevelsStyle(CAGE_gr6) <- "NCBI" }, silent = TRUE)

valid_levels <- intersect(GenomeInfoDb::seqlevels(gr_flat), GenomeInfoDb::seqlevels(CAGE_gr6))
valid_levels <- intersect(valid_levels, fasta_levels_6)

gr_flat  <- GenomeInfoDb::keepSeqlevels(gr_flat, valid_levels, pruning.mode = "coarse")
CAGE_gr6 <- GenomeInfoDb::keepSeqlevels(CAGE_gr6, valid_levels, pruning.mode = "coarse")

# Set seqlengths for safe trimming (from FASTA)
GenomeInfoDb::seqlengths(gr_flat) <- fasta_seqlens_6[GenomeInfoDb::seqlevels(gr_flat)]
gr_flat <- GenomicRanges::trim(gr_flat)

# Enforce expected widths post-trim
expected_widths2 <- expected_widths[names(gr_flat)]
gr_flat <- gr_flat[IRanges::width(gr_flat) == expected_widths2]

# FASTA-safe sequences + drop bad tx
seq_out <- safe_getSeq_drop_bad(FA_6, gr_flat)
ok_names <- seq_out$ok_names
if (length(ok_names) < 50) message("WARNING: very few transcripts left after FASTA-safe filtering: ", length(ok_names))

gr_flat      <- gr_flat[names(gr_flat) %in% ok_names]
fiveUTR_seqs <- seq_out$seqs[ok_names]
fiveUTR_grl  <- split(gr_flat, names(gr_flat))

message("Figure6: transcripts retained after FASTA-safe filtering: ", length(fiveUTR_grl))

# ------------------------------------------------------------------------------
# 6.0b MAP CAGE TO TRANSCRIPTS + BUILD dt_raw
# ------------------------------------------------------------------------------
mapped_gr_raw <- GenomicFeatures::mapToTranscripts(CAGE_gr6, fiveUTR_grl)

dt_raw <- data.table::as.data.table(BiocGenerics::as.data.frame(mapped_gr_raw))
dt_raw[, tx_name := as.character(seqnames)]
dt_raw[, count := S4Vectors::mcols(CAGE_gr6)$count[xHits]]
dt_raw <- dt_raw[!is.na(tx_name) & !is.na(start) & !is.na(count)]
dt_raw[, start := as.integer(start)]
dt_raw[, count := as.numeric(count)]
dt_raw <- dt_raw[, .(count = sum(count, na.rm = TRUE)), by = .(tx_name, start)]
data.table::setorder(dt_raw, tx_name, start)
dt_raw[, `:=`(cum_sum = cumsum(count), total_sum = sum(count, na.rm = TRUE)), by = tx_name]
dt_raw <- dt_raw[total_sum > 0.1]

# ------------------------------------------------------------------------------
# 6.1 THRESHOLD ANALYSIS (50% CAGE) + CHI-SQUARE/CRAMER'S V ALL PAIRS
# ------------------------------------------------------------------------------
run_threshold_analysis <- function(threshold_val, dt_source, seqs_source, uorf_gene_ids_vec, gr_lens, make_main_plots = FALSE) {
  
  # first transcript coord reaching threshold
  dt_point <- dt_source[cum_sum / total_sum >= threshold_val, .SD[1], by = tx_name]
  
  # tx lengths
  tx_len_vec <- vapply(gr_lens, function(x) sum(IRanges::width(x)), numeric(1))
  dt_point[, tx_len := tx_len_vec[tx_name]]
  dt_point <- dt_point[!is.na(tx_len)]
  
  # window fit: (pos-60 .. pos+62)
  dt_valid <- dt_point[start >= 61 & (tx_len - start) >= 62]
  if (nrow(dt_valid) < 10) return(NULL)
  
  current_seqs <- seqs_source[dt_valid$tx_name]
  win_start <- dt_valid$start - 60L
  win_end   <- dt_valid$start + 62L
  window_seqs <- Biostrings::subseq(current_seqs, start = win_start, end = win_end)
  
  # ATG positions
  atg_matches <- Biostrings::vmatchPattern("ATG", window_seqs)
  atg_starts  <- as.list(Biostrings::startIndex(atg_matches))
  atg_rel     <- lapply(atg_starts, function(x) as.integer(x) - 61L)  # -60..+60
  names(atg_rel) <- dt_valid$tx_name
  
  tx_ids   <- dt_valid$tx_name
  gene_ids <- toupper(get_gene_id_from_tx(tx_ids))
  is_tuorf_bool <- gene_ids %in% toupper(uorf_gene_ids_vec)
  
  n_t <- sum(is_tuorf_bool)
  n_n <- sum(!is_tuorf_bool)
  
  label_t <- paste0("TuORFs\n(n=", n_t, ")")
  label_n <- paste0("No TuORF\n(n=", n_n, ")")
  
  stats_cat_vec <- ifelse(is_tuorf_bool, "TuORF", "NoTuORF")
  plot_cat_vec  <- ifelse(is_tuorf_bool, label_t, label_n)
  
  # seqlogos (crop 21..71)
  s_t <- as.character(window_seqs[is_tuorf_bool])
  s_n <- as.character(window_seqs[!is_tuorf_bool])
  s_t_log <- substr(s_t, 21, 71)
  s_n_log <- substr(s_n, 21, 71)
  
  pC <- if (length(s_t) > 5) ggseqlogo::ggseqlogo(s_t_log) +
    ggplot2::ggtitle(paste0("TuORF (", threshold_val * 100, "%, n=", n_t, ")")) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
    ggplot2::ylim(0, 0.6) else ggplot2::ggplot() + ggplot2::theme_void()
  
  pD <- if (length(s_n) > 5) ggseqlogo::ggseqlogo(s_n_log) +
    ggplot2::ggtitle(paste0("No TuORF (", threshold_val * 100, "%, n=", n_n, ")")) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
    ggplot2::ylim(0, 0.6) else ggplot2::ggplot() + ggplot2::theme_void()
  
  if (!make_main_plots) return(list(pC = pC, pD = pD))
  
  # A: ATG frequency lines
  make_line_df <- function(hit_list, mask, name) {
    if (sum(mask) == 0) return(NULL)
    hits <- unlist(hit_list[mask], use.names = FALSE)
    tbl <- table(factor(hits, levels = -60:60))
    data.frame(X_Axis = -60:60, FreqPct = (as.integer(tbl) / sum(mask)) * 100, Category = name)
  }
  
  df_line <- dplyr::bind_rows(
    make_line_df(atg_rel,  is_tuorf_bool,  label_t),
    make_line_df(atg_rel, !is_tuorf_bool,  label_n)
  )
  df_line$Category <- factor(df_line$Category, levels = c(label_n, label_t))
  
  # Up/Down counts per transcript
  count_range <- function(v, lo, hi) if (length(v) == 0) 0L else sum(v >= lo & v <= hi)
  up_counts   <- vapply(atg_rel, count_range, integer(1), -60, -1)
  down_counts <- vapply(atg_rel, count_range, integer(1),   0, 60)
  
  df_counts <- data.frame(
    Stats_Cat  = stats_cat_vec,
    Plot_Cat   = plot_cat_vec,
    Upstream   = up_counts,
    Downstream = down_counts,
    stringsAsFactors = FALSE
  ) %>%
    tidyr::pivot_longer(cols = c(Upstream, Downstream), names_to = "Region", values_to = "Count") %>%
    dplyr::mutate(
      Region = factor(Region, levels = c("Upstream", "Downstream")),
      Count_Group = factor(dplyr::case_when(Count == 0 ~ "0", Count == 1 ~ "1", TRUE ~ "2+"),
                           levels = c("0", "1", "2+"), ordered = TRUE),
      Stats_Cat = factor(ifelse(Stats_Cat == "TuORF", "TuORF", "NoTuORF"),
                         levels = c("NoTuORF", "TuORF"))
    )
  
  # Bar plot summary
  df_bar_summary <- df_counts %>%
    dplyr::group_by(Region, Plot_Cat, Count_Group) %>%
    dplyr::summarise(N = dplyr::n(), .groups = "drop_last") %>%
    dplyr::mutate(Pct = N / sum(N) * 100) %>%
    dplyr::ungroup()
  
  # Pairwise chi-square + Cramer's V (all pairs among 4 groups)
  df_pairs <- df_counts %>%
    dplyr::mutate(Group4 = paste0(as.character(Region), "_", as.character(Stats_Cat)))
  
  group_levels <- c("Upstream_NoTuORF", "Upstream_TuORF", "Downstream_NoTuORF", "Downstream_TuORF")
  df_pairs$Group4 <- factor(df_pairs$Group4, levels = group_levels)
  
  pair_list <- utils::combn(group_levels, 2, simplify = FALSE)
  
  pair_res <- lapply(pair_list, function(gp) {
    dsub <- df_pairs %>% dplyr::filter(Group4 %in% gp)
    n_by <- table(dsub$Group4)
    tbl  <- table(factor(dsub$Count_Group, levels = c("0","1","2+")),
                  factor(dsub$Group4, levels = gp))
    
    out <- chisq_one(tbl)
    data.frame(
      group1 = gp[1], group2 = gp[2],
      n1 = as.integer(n_by[gp[1]]),
      n2 = as.integer(n_by[gp[2]]),
      chi2 = out$chi2, df = out$df, p = out$p, V = out$V,
      stringsAsFactors = FALSE
    )
  })
  
  pair_df <- dplyr::bind_rows(pair_res) %>%
    dplyr::mutate(p_adj_BH = stats::p.adjust(p, method = "BH"))
  
  utils::write.csv(
    pair_df,
    file.path(OUT_DIR, paste0("Figure6_pairwise_chisq_cramerv_", threshold_val*100, "pct.csv")),
    row.names = FALSE
  )
  
  cap_lines <- apply(pair_df, 1, function(r) {
    paste0(
      r[["group1"]], " vs ", r[["group2"]],
      " | χ²=", round(as.numeric(r[["chi2"]]), 2),
      " df=", as.integer(r[["df"]]),
      " p(BH)=", format.pval(as.numeric(r[["p_adj_BH"]]), digits = 2, eps = 1e-3),
      " V=", round(as.numeric(r[["V"]]), 3),
      " (n=", as.integer(r[["n1"]]), ",", as.integer(r[["n2"]]), ")"
    )
  })
  
  pairwise_caption <- paste0(
    "Pairwise χ² tests on Count_Group (0/1/2+) distributions + Cramer's V (effect size)\n",
    paste(cap_lines, collapse = "\n")
  )
  
  message("=== Figure6: Pairwise χ² + Cramer's V (BH-adjusted p) ===")
  print(pair_df)
  
  # Plot A
  pA <- ggplot2::ggplot(df_line, ggplot2::aes(x = X_Axis, y = FreqPct, color = Category)) +
    ggplot2::geom_line(size = 0.8) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
    ggplot2::scale_color_manual(values = c("#48B5C4", "#E76F69")) +
    ggplot2::theme_light() +
    ggplot2::labs(title = "A: ATG Frequency", x = "Position relative to threshold point", y = "% Transcripts") +
    ggplot2::theme(legend.position = c(0.85, 0.8), legend.background = ggplot2::element_blank())
  
  # Plot B
  pB <- ggplot2::ggplot(df_bar_summary, ggplot2::aes(x = Plot_Cat, y = Pct, fill = forcats::fct_rev(Count_Group))) +
    ggplot2::geom_bar(stat = "identity", position = "stack", width = 0.6) +
    ggplot2::facet_wrap(~Region) +
    ggplot2::scale_fill_manual(
      values = c("0" = "grey90", "1" = "#7FB3D5", "2+" = "#2471A3"),
      name = "# ATGs", breaks = c("2+", "1", "0")
    ) +
    ggplot2::theme_classic() +
    ggplot2::coord_cartesian(ylim = c(0, 115)) +
    ggplot2::labs(title = "B: ATG Count (Chi-square / Cramer's V)", y = "% Transcripts", x = "")
  
  list(pA = pA, pB = pB, pC = pC, pD = pD,
       pairwise_caption = pairwise_caption,
       pairwise_table = pair_df)
}

# ------------------------------------------------------------------------------
# 6.2 RUN 50% (AB + seqlogos) AND 30/80% (seqlogos only)
# ------------------------------------------------------------------------------
res_50 <- run_threshold_analysis(
  threshold_val = 0.50,
  dt_source = dt_raw,
  seqs_source = fiveUTR_seqs,
  uorf_gene_ids_vec = uORF_gene_ids,
  gr_lens = fiveUTR_grl,
  make_main_plots = TRUE
)

thresholds <- c(0.30, 0.80)
names(thresholds) <- paste0(thresholds * 100, "%")
extra_plots <- lapply(thresholds, function(x) {
  run_threshold_analysis(
    threshold_val = x,
    dt_source = dt_raw,
    seqs_source = fiveUTR_seqs,
    uorf_gene_ids_vec = uORF_gene_ids,
    gr_lens = fiveUTR_grl,
    make_main_plots = FALSE
  )
})

# ------------------------------------------------------------------------------
# 6.3 SAVE FIGURE 6 AB (50%) WITH PAIRWISE χ² + V CAPTION
# ------------------------------------------------------------------------------
if (!is.null(res_50)) {
  fig6_ab <- ggpubr::ggarrange(res_50$pA, res_50$pB, ncol = 2, widths = c(1.3, 1))
  fig6_ab_annotated <- ggpubr::annotate_figure(
    fig6_ab,
    top = ggpubr::text_grob("Figure 6: A and B (50% CAGE Threshold)", face = "bold", size = 16),
    bottom = ggpubr::text_grob(res_50$pairwise_caption, size = 8, hjust = 0)
  )
  print(fig6_ab_annotated)
  save_pdf(fig6_ab_annotated, "Figure6_AB_50pct_CAGE_threshold.pdf", w = 14, h = 7)
}

# ------------------------------------------------------------------------------
# 6.4 SAVE FIGURE 6 SEQLOGOS (30/50/80)
# ------------------------------------------------------------------------------
if (!is.null(res_50) && !is.null(extra_plots[["30%"]]) && !is.null(extra_plots[["80%"]])) {
  
  fig6_logos <- ggpubr::ggarrange(
    ggpubr::ggarrange(extra_plots[["30%"]]$pC, extra_plots[["30%"]]$pD, ncol = 2, labels = c("30% TuORF", "30% No TuORF")),
    ggpubr::ggarrange(res_50$pC, res_50$pD, ncol = 2, labels = c("50% TuORF", "50% No TuORF")),
    ggpubr::ggarrange(extra_plots[["80%"]]$pC, extra_plots[["80%"]]$pD, ncol = 2, labels = c("80% TuORF", "80% No TuORF")),
    nrow = 3
  )
  
  fig6_logos <- ggpubr::annotate_figure(
    fig6_logos,
    top = ggpubr::text_grob("Figure 6: All SeqLogos (y-axis 0 to 0.6)", face = "bold", size = 16)
  )
  
  print(fig6_logos)
  save_pdf(fig6_logos, "Figure6_All_SeqLogos.pdf", w = 12, h = 15)
}

# close FASTA handle for Figure 6
try(close(FA_6), silent = TRUE)

message("FIGURE 6 complete. Outputs in OUT_DIR = ", OUT_DIR)

################################################################################
# FIGURE 7. DLA and qPCR
################################################################################

message("FIGURE 7: starting...")

# --- PATHS (declared in section 0; re-checked here) ---
stopifnot(file.exists(PATH_DLA_WEAK_UORF_UUORF))
stopifnot(file.exists(PATH_QPCR_WEAK_UORF_UUORF))
stopifnot(file.exists(PATH_QPCR_PERFECT_KOZAK))

# ------------------------------------------------------------------------------
# 7.1 DLA for weak uORF and UuORFs (WT vs Mut only)
# ------------------------------------------------------------------------------
raw_dla_1 <- openxlsx::read.xlsx(PATH_DLA_WEAK_UORF_UUORF)

plot_data_7_1 <- raw_dla_1 %>%
  dplyr::rename(Lines = Plasmids) %>%
  dplyr::mutate(
    Repeat = as.factor(Repeat),
    Lines  = as.character(Lines)
  ) %>%
  dplyr::filter(Lines %in% c("WT", "Mut")) %>%
  dplyr::mutate(Lines = factor(Lines, levels = c("WT", "Mut"))) %>%
  tidyr::pivot_longer(
    cols = -c(Repeat, Lines),
    names_to = "Gene",
    values_to = "Expression"
  )

comparison_list_7 <- list(c("WT", "Mut"))

p7_1 <- ggplot2::ggplot(plot_data_7_1, ggplot2::aes(x = Lines, y = Expression, fill = Lines)) +
  ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  ggplot2::geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  ggplot2::facet_grid(Gene ~ Repeat, scales = "free_y") +
  ggpubr::stat_compare_means(
    comparisons = comparison_list_7,
    method = "t.test",
    label = "p.signif",
    hide.ns = FALSE,
    vjust = 0.5
  ) +
  ggplot2::scale_fill_manual(values = c("WT" = "#888888", "Mut" = "#E69F00")) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "DLA Assay Analysis (WT vs Mut)",
    subtitle = "Facet Rows: Genes | Facet Columns: Biological Repeats",
    x = "Plasmid / Line",
    y = "DLA Value"
  ) +
  ggplot2::theme(
    strip.background = ggplot2::element_rect(fill = "#f0f0f0"),
    strip.text = ggplot2::element_text(face = "bold", size = 10),
    legend.position = "none"
  )

print(p7_1)
save_pdf(p7_1, "DLA_Analysis_Plot_WuORFs_UuORFs.pdf", w = 10, h = 15)

# ------------------------------------------------------------------------------
# 7.2 qPCR for weak uORF and UuORFs (WT vs Mut only)
# ------------------------------------------------------------------------------
raw_qpcr_1 <- openxlsx::read.xlsx(PATH_QPCR_WEAK_UORF_UUORF)

plot_data_7_2 <- raw_qpcr_1 %>%
  dplyr::mutate(
    Repeat = as.factor(Repeat),
    Lines  = as.character(Lines)
  ) %>%
  dplyr::filter(Lines %in% c("WT", "Mut")) %>%
  dplyr::mutate(Lines = factor(Lines, levels = c("WT", "Mut"))) %>%
  tidyr::pivot_longer(
    cols = -c(Repeat, Lines, tech_repeat),
    names_to = "Gene",
    values_to = "Expression"
  )

p7_2 <- ggplot2::ggplot(plot_data_7_2, ggplot2::aes(x = Lines, y = Expression, fill = Lines)) +
  ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  ggplot2::geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  ggplot2::facet_grid(Gene ~ Repeat, scales = "free_y") +
  ggpubr::stat_compare_means(
    comparisons = comparison_list_7,
    method = "t.test",
    label = "p.signif",
    hide.ns = FALSE,
    vjust = 0.5
  ) +
  ggplot2::scale_fill_manual(values = c("WT" = "#888888", "Mut" = "#E69F00")) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "qPCR Expression Analysis per Gene (WT vs Mut)",
    subtitle = "Facet Rows: Genes | Facet Columns: Biological Repeats",
    x = "Line",
    y = "Relative Expression"
  ) +
  ggplot2::theme(
    strip.background = ggplot2::element_rect(fill = "#f0f0f0"),
    strip.text = ggplot2::element_text(face = "bold", size = 10),
    legend.position = "none"
  )

print(p7_2)
save_pdf(p7_2, "qPCR_Analysis_Plot_WuORFs_UuORFs.pdf", w = 10, h = 15)

# ------------------------------------------------------------------------------
# 7.3 DLA for perfect Kozak (Sheet 2)
# ------------------------------------------------------------------------------
raw_dla_2 <- openxlsx::read.xlsx(PATH_DLA_WEAK_UORF_UUORF, sheet = 2)

plot_data_7_3 <- raw_dla_2 %>%
  dplyr::rename(Lines = Plasmids) %>%
  dplyr::mutate(Repeat = as.factor(Repeat)) %>%
  dplyr::mutate(Lines = factor(Lines, levels = c("WT", "GK"))) %>%
  tidyr::pivot_longer(
    cols = -c(Repeat, Lines),
    names_to = "Gene",
    values_to = "Expression"
  )

comparison_list_7_3 <- list(c("WT", "GK"))

p7_3 <- ggplot2::ggplot(plot_data_7_3, ggplot2::aes(x = Lines, y = Expression, fill = Lines)) +
  ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  ggplot2::geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  ggplot2::facet_grid(Gene ~ Repeat, scales = "free_y") +
  ggpubr::stat_compare_means(
    comparisons = comparison_list_7_3,
    method = "t.test",
    label = "p.signif",
    hide.ns = FALSE,
    vjust = 0.5
  ) +
  ggplot2::scale_fill_manual(values = c("WT" = "#888888", "GK" = "#009E73")) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "DLA Assay Analysis (Perfect Kozak, Sheet 2)",
    subtitle = "Facet Rows: Genes | Facet Columns: Biological Repeats",
    x = "Plasmid / Line",
    y = "DLA Value"
  ) +
  ggplot2::theme(
    strip.background = ggplot2::element_rect(fill = "#f0f0f0"),
    strip.text = ggplot2::element_text(face = "bold", size = 10),
    legend.position = "none"
  )

print(p7_3)
save_pdf(p7_3, "DLA_Perfect_Kozak.pdf", w = 10, h = 15)

# ------------------------------------------------------------------------------
# 7.4 qPCR for perfect Kozak
# ------------------------------------------------------------------------------
raw_qpcr_2 <- openxlsx::read.xlsx(PATH_QPCR_PERFECT_KOZAK)

plot_data_7_4 <- raw_qpcr_2 %>%
  dplyr::mutate(Repeat = as.factor(Repeat)) %>%
  dplyr::mutate(Lines = factor(Lines, levels = c("WT", "KK"))) %>%
  tidyr::pivot_longer(
    cols = -c(Repeat, Lines, tech_repeat),
    names_to = "Gene",
    values_to = "Expression"
  )

comparison_list_7_4 <- list(c("WT", "KK"))

p7_4 <- ggplot2::ggplot(plot_data_7_4, ggplot2::aes(x = Lines, y = Expression, fill = Lines)) +
  ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  ggplot2::geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  ggplot2::facet_grid(Gene ~ Repeat, scales = "free_y") +
  ggpubr::stat_compare_means(
    comparisons = comparison_list_7_4,
    method = "t.test",
    label = "p.signif",
    hide.ns = FALSE,
    vjust = 0.5
  ) +
  ggplot2::scale_fill_manual(values = c("WT" = "#888888", "KK" = "#D55E00")) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "qPCR Expression Analysis: Perfect Kozak",
    subtitle = "Facet Rows: Genes | Facet Columns: Biological Repeats",
    x = "Line",
    y = "Relative Expression"
  ) +
  ggplot2::theme(
    strip.background = ggplot2::element_rect(fill = "#f0f0f0"),
    strip.text = ggplot2::element_text(face = "bold", size = 10),
    legend.position = "none"
  )

print(p7_4)
save_pdf(p7_4, "qPCR_perfect_Kozak.pdf", w = 10, h = 15)

message("FIGURE 7 complete.")

################################################################################
# FIGURE 8 alternative: uORF LENGTH ANALYSIS (VIOLIN + BOX OVERLAY)
# OUTPUT: Figure8_ABCDEFGH_Full_Violin.pdf
#         Figure8_uORF_length_before_after_TSSfilter.csv
# (panel H: uORF length before vs after TSS filtering,
#  in the previously empty slot next to panel E; output PDF renamed G -> GH)
################################################################################

message("FIGURE 8: starting...")

.require_pkg("RColorBrewer")

# Ensure no-uORF TE table exists (consistent with the earlier TE objects)
# - Preferred existing object: TE_NoTuORF from Fig1/3
# - If missing, rebuild from TE_tbl + OUT_df5_clean
if (!exists("TE_NoTuORF")) {
  if (exists("TE_tbl") && exists("OUT_df5_clean")) {
    TE_NoTuORF <- TE_tbl %>%
      dplyr::filter(!(gene_id %in% unique(OUT_df5_clean$gene_id))) %>%
      dplyr::transmute(gene_id, TE, Group_Label = "No TuORF")
  } else {
    stop("Figure8 needs TE_NoTuORF or (TE_tbl + OUT_df5_clean) to build it.")
  }
}

# Use TE_NoTuORF as the "No uORF" reference expected by the Figure 8 code
TE_CDS_Other <- TE_NoTuORF %>% dplyr::select(gene_id, TE)

# Sanity
stopifnot(exists("OUT_df6"))
req_uorf_cols <- c("gene_id", "TE", "uORF_length", "pUP", "uIF")
missing_uorf_cols <- setdiff(req_uorf_cols, colnames(OUT_df6))
if (length(missing_uorf_cols) > 0) {
  stop("OUT_df6 missing required columns for Figure 8: ", paste(missing_uorf_cols, collapse = ", "))
}

# Reference (No uORF)
df_no_uorf <- TE_CDS_Other %>%
  dplyr::select(gene_id, TE) %>%
  dplyr::mutate(
    uORF_length = NA_real_,
    Length_Bin  = "No uORF",
    pUP         = NA_real_,
    uIF         = NA_real_,
    log2_uIF    = NA_real_
  )

# uORF genes
df_uorf <- OUT_df6 %>%
  dplyr::filter(!is.na(uORF_length), uORF_length >= 1) %>%
  dplyr::mutate(
    Length_Bin = dplyr::case_when(
      uORF_length < 2   ~ "1",
      uORF_length <= 10 ~ "2-10",
      uORF_length <= 20 ~ "11-20",
      uORF_length <= 30 ~ "21-30",
      TRUE              ~ ">30"
    ),
    log2_uIF = dplyr::if_else(is.finite(uIF_raw) & uIF_raw > 0, log2(uIF_raw), NA_real_)
  ) %>%
  dplyr::select(gene_id, TE, uORF_length, Length_Bin, pUP, uIF, log2_uIF)

level_order_uorf <- c("1", "2-10", "11-20", "21-30", ">30")

create_fig8_violinbox <- function(data_df, x_col, y_col, plot_title, color_pal, y_limit = NULL) {
  
  data_df[[x_col]] <- as.character(data_df[[x_col]])
  has_no_uorf <- any(data_df[[x_col]] == "No uORF", na.rm = TRUE)
  
  level_order <- if (has_no_uorf) c("No uORF", level_order_uorf) else level_order_uorf
  data_df[[x_col]] <- factor(data_df[[x_col]], levels = level_order)
  
  data_df <- data_df %>%
    dplyr::filter(!is.na(.data[[y_col]])) %>%
    dplyr::filter(!is.na(.data[[x_col]]))
  
  n_groups <- length(unique(as.character(data_df[[x_col]])))
  do_stats <- (n_groups >= 2)
  
  counts_df <- data_df %>%
    dplyr::group_by(.data[[x_col]]) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")
  
  counts_df[[x_col]] <- as.character(counts_df[[x_col]])
  counts_df$Label_Text <- paste0(counts_df[[x_col]], "\n(n=", counts_df$n, ")")
  
  label_levels <- counts_df$Label_Text[match(level_order, counts_df[[x_col]])]
  label_levels <- label_levels[!is.na(label_levels)]
  
  if (do_stats) {
    stat_test <- stats::pairwise.wilcox.test(data_df[[y_col]], data_df[[x_col]], p.adjust.method = "BH", exact = FALSE)
    if (!is.null(stat_test$p.value)) {
      p_matrix <- stat_test$p.value
      p_matrix[is.na(p_matrix)] <- 1
      p_vec <- rcompanion::fullPTable(p_matrix)
      letters_list <- multcompView::multcompLetters(
        p_vec, compare = "<", threshold = 0.05, Letters = letters, reversed = FALSE
      )
      letter_df <- data.frame(
        bin = names(letters_list$Letters),
        Letters = letters_list$Letters,
        stringsAsFactors = FALSE
      )
      colnames(letter_df)[1] <- x_col
    } else {
      letter_df <- data.frame(bin = character(0), Letters = character(0), stringsAsFactors = FALSE)
      colnames(letter_df)[1] <- x_col
    }
  } else {
    letter_df <- data.frame(bin = character(0), Letters = character(0), stringsAsFactors = FALSE)
    colnames(letter_df)[1] <- x_col
  }
  
  df_final <- data_df %>% dplyr::left_join(counts_df, by = x_col)
  df_final$Label_Text <- factor(df_final$Label_Text, levels = label_levels)
  
  if (nrow(letter_df) > 0) {
    letter_df <- letter_df %>% dplyr::left_join(counts_df, by = x_col)
    letter_df$Label_Text <- factor(letter_df$Label_Text, levels = label_levels)
  }
  
  n_colors <- length(level_order)
  if (n_colors > 2) {
    cols <- if ("No uORF" %in% level_order) {
      c(
        "No uORF" = "#999999",
        setNames(
          colorRampPalette(RColorBrewer::brewer.pal(5, color_pal))(n_colors - 1),
          level_order[-1]
        )
      )
    } else {
      setNames(colorRampPalette(RColorBrewer::brewer.pal(5, color_pal))(n_colors), level_order)
    }
  } else {
    cols <- scales::hue_pal()(n_colors)
    names(cols) <- level_order
  }
  
  p <- ggplot2::ggplot(df_final, ggplot2::aes(x = Label_Text, y = .data[[y_col]], fill = .data[[x_col]])) +
    ggplot2::geom_violin(trim = TRUE, alpha = 0.7, scale = "width") +
    ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.9) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = plot_title, x = "uORF Length (AA)", y = y_col) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold", size = 10),
      axis.text.x = ggplot2::element_text(size = 9)
    )
  
  if (nrow(letter_df) > 0) {
    y_pos <- if (!is.null(y_limit)) y_limit[2] * 0.90 else max(df_final[[y_col]], na.rm = TRUE) * 1.05
    p <- p + ggplot2::geom_text(
      data = letter_df,
      ggplot2::aes(x = Label_Text, y = y_pos, label = Letters),
      inherit.aes = FALSE, size = 4, fontface = "bold"
    )
  }
  
  if (!is.null(y_limit)) p <- p + ggplot2::coord_cartesian(ylim = y_limit)
  p
}

p8A <- create_fig8_violinbox(
  data_df    = dplyr::bind_rows(df_uorf, df_no_uorf),
  x_col      = "Length_Bin",
  y_col      = "TE",
  plot_title = "A: All uORFs vs No uORF",
  color_pal  = "Blues",
  y_limit    = c(0, 3)
)

p8B <- create_fig8_violinbox(
  data_df    = dplyr::bind_rows(df_uorf %>% dplyr::filter(pUP > 30), df_no_uorf),
  x_col      = "Length_Bin",
  y_col      = "TE",
  plot_title = "B: Strong uORFs (pUP > 30) vs No uORF",
  color_pal  = "Greens",
  y_limit    = c(0, 3)
)

df_bar <- df_uorf %>%
  dplyr::mutate(Cutoff_Group = dplyr::if_else(pUP > 30, "Kept (>30% pUP)", "Removed (<=30% pUP)")) %>%
  dplyr::group_by(Length_Bin, Cutoff_Group) %>%
  dplyr::summarise(Count = dplyr::n(), .groups = "drop")

df_bar$Length_Bin <- factor(df_bar$Length_Bin, levels = level_order_uorf)
df_bar$Cutoff_Group <- factor(df_bar$Cutoff_Group, levels = c("Removed (<=30% pUP)", "Kept (>30% pUP)"))

p8C <- ggplot2::ggplot(df_bar, ggplot2::aes(x = Length_Bin, y = Count, fill = Cutoff_Group)) +
  ggplot2::geom_bar(stat = "identity", position = "stack", width = 0.7) +
  ggplot2::geom_text(
    ggplot2::aes(label = Count),
    position = ggplot2::position_stack(vjust = 0.5),
    color = "white", size = 3, fontface = "bold"
  ) +
  ggplot2::scale_fill_manual(values = c("Removed (<=30% pUP)" = "#E74C3C", "Kept (>30% pUP)" = "#2ECC71")) +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "C: Genes Removed vs Kept", x = "uORF Length (AA)", y = "Count", fill = "") +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 10), legend.position = "top")

df_perc <- df_uorf %>%
  dplyr::mutate(Is_Removed = dplyr::if_else(pUP <= 30, 1L, 0L)) %>%
  dplyr::group_by(Length_Bin) %>%
  dplyr::summarise(
    Total = dplyr::n(),
    Removed = sum(Is_Removed),
    Percent = round(Removed / Total * 100, 1),
    .groups = "drop"
  )

df_perc$Length_Bin <- factor(df_perc$Length_Bin, levels = level_order_uorf)

p8D <- ggplot2::ggplot(df_perc, ggplot2::aes(x = Length_Bin, y = Percent, fill = Length_Bin)) +
  ggplot2::geom_bar(stat = "identity", color = "black", alpha = 0.8) +
  ggplot2::geom_text(ggplot2::aes(label = paste0(Percent, "%")), vjust = -0.5, size = 3.5, fontface = "bold") +
  ggplot2::scale_fill_brewer(palette = "OrRd") +
  ggplot2::coord_cartesian(ylim = c(0, 100)) +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "D: % Removed by 30% pUP Cutoff", x = "uORF Length (AA)", y = "% Removed") +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 10), legend.position = "none")

df_uorf$Length_Bin <- factor(df_uorf$Length_Bin, levels = level_order_uorf)

p8E <- ggplot2::ggplot(df_uorf, ggplot2::aes(x = Length_Bin, y = pUP, fill = Length_Bin)) +
  ggplot2::geom_violin(trim = TRUE, alpha = 0.6, scale = "width") +
  ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  ggplot2::scale_fill_brewer(palette = "Purples") +
  ggplot2::coord_cartesian(ylim = c(0, 100)) +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "E: pUP Distribution", x = "uORF Length (AA)", y = "% pUP") +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 10), legend.position = "none")

p8F <- create_fig8_violinbox(
  data_df    = df_uorf,
  x_col      = "Length_Bin",
  y_col      = "log2_uIF",
  plot_title = "F: log2(uIF) vs Length (All uORFs)",
  color_pal  = "Oranges",
  y_limit    = c(-7, 8)
)

p8G <- create_fig8_violinbox(
  data_df    = df_uorf %>% dplyr::filter(pUP > 30),
  x_col      = "Length_Bin",
  y_col      = "log2_uIF",
  plot_title = "G: log2(uIF) vs Length (Strong uORFs >30% pUP)",
  color_pal  = "Oranges",
  y_limit    = c(-7, 8)
)

# ------------------------------------------------------------------------------
# Panel H: uORF length BEFORE vs AFTER TSS filtering
# "Before" = all single-TuORF uORFs (df_uorf); "After" = those retained at
# pUP > 30; "Removed" = pUP <= 30 shown for contrast. The Wilcoxon test compares
# Retained vs Removed only (disjoint sets; "All" overlaps both, so testing
# against it would be invalid).
# ------------------------------------------------------------------------------
df_len_H <- dplyr::bind_rows(
  df_uorf %>% dplyr::transmute(uORF_length, Filter_Group = "All TuORFs\n(before)"),
  df_uorf %>% dplyr::filter(pUP > 30)  %>% dplyr::transmute(uORF_length, Filter_Group = "Retained\n(>30% uORF+)"),
  df_uorf %>% dplyr::filter(pUP <= 30) %>% dplyr::transmute(uORF_length, Filter_Group = "Removed\n(<=30% uORF+)")
)
df_len_H$Filter_Group <- factor(df_len_H$Filter_Group,
  levels = c("All TuORFs\n(before)", "Retained\n(>30% uORF+)", "Removed\n(<=30% uORF+)"))

len_H_stats <- df_len_H %>%
  dplyr::group_by(Filter_Group) %>%
  dplyr::summarise(n = dplyr::n(),
                   median_len_aa = stats::median(uORF_length),
                   mean_len_aa   = round(mean(uORF_length), 1),
                   .groups = "drop")

wt_H <- stats::wilcox.test(df_uorf$uORF_length[df_uorf$pUP > 30],
                           df_uorf$uORF_length[df_uorf$pUP <= 30])

lab_H <- len_H_stats %>%
  dplyr::mutate(lab = paste0("n=", n, "\nmedian=", median_len_aa))

p8H <- ggplot2::ggplot(df_len_H, ggplot2::aes(x = Filter_Group, y = uORF_length, fill = Filter_Group)) +
  ggplot2::geom_violin(trim = TRUE, alpha = 0.7, scale = "width") +
  ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.9) +
  ggplot2::scale_y_log10(expand = ggplot2::expansion(mult = c(0.05, 0.22))) +
  ggplot2::scale_fill_manual(values = c("All TuORFs\n(before)" = "#999999",
                                        "Retained\n(>30% uORF+)" = "#2ECC71",
                                        "Removed\n(<=30% uORF+)" = "#E74C3C")) +
  ggplot2::geom_text(data = lab_H,
                     ggplot2::aes(x = Filter_Group, y = max(df_len_H$uORF_length) * 1.6, label = lab),
                     inherit.aes = FALSE, size = 2.8) +
  ggplot2::theme_classic() +
  ggplot2::labs(
    title = "H: uORF length before vs after TSS filtering",
    subtitle = sprintf("Retained vs Removed: Wilcoxon p = %.2e", wt_H$p.value),
    x = "", y = "uORF length (AA, log10 scale)") +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 10),
                 plot.subtitle = ggplot2::element_text(size = 8),
                 legend.position = "none",
                 axis.text.x = ggplot2::element_text(size = 8))

utils::write.csv(
  len_H_stats %>% dplyr::mutate(wilcoxon_retained_vs_removed_p = signif(wt_H$p.value, 3)),
  file.path(OUT_DIR, "Figure8_uORF_length_before_after_TSSfilter.csv"), row.names = FALSE)
message(sprintf(
  "Panel H: uORF length (AA) — all n=%d median=%g | retained (pUP>30) n=%d median=%g | removed n=%d median=%g | Wilcoxon retained vs removed p=%.2e",
  len_H_stats$n[1], len_H_stats$median_len_aa[1],
  len_H_stats$n[2], len_H_stats$median_len_aa[2],
  len_H_stats$n[3], len_H_stats$median_len_aa[3], wt_H$p.value))

Figure_8_Final <- ggpubr::ggarrange(
  ggpubr::ggarrange(p8A, p8B, ncol = 2),
  ggpubr::ggarrange(p8C, p8D, ncol = 2),
  ggpubr::ggarrange(p8E, p8H, ncol = 2),
  ggpubr::ggarrange(p8F, p8G, ncol = 2),
  nrow = 4
)

print(Figure_8_Final)
save_pdf(Figure_8_Final, "Figure8_ABCDEFGH_Full_Violin.pdf", w = 12, h = 20)

message("FIGURE 8 complete.")

################################################################################
# GO SLIM ENRICHMENT — BP/MF/CC bubble plots for TWO gene groupings
#   (A) by uIF     : 1-uORF Weak / Mid / Strong          + (>=1 strong uORF or >=2 uORFs)
#   (B) by %uORF+  : 1-uORF 0~30% / 30~80% / 80~100%     + (>=1 strong uORF or >=2 uORFs)
# Both groupings end with the SAME 4th column (a shared reference group), so the
# two figure sets are directly comparable; that column's numbers are identical.
#
# Outputs in OUT_DIR (filenames state the grouping and the number of groups):
#   GO_BP_by_uIF_4groups_Top15.pdf            GO_BP_by_pct_uORF_plus_4groups_Top15.pdf
#   GO_MF_by_uIF_4groups_Top15.pdf            GO_MF_by_pct_uORF_plus_4groups_Top15.pdf
#   GO_CC_by_uIF_4groups_Top15.pdf            GO_CC_by_pct_uORF_plus_4groups_Top15.pdf
#   GO_enrichment_results_by_uIF.csv          GO_enrichment_results_by_pct_uORF_plus.csv
# (These replace the older ambiguous names GO_BP/MF/CC_4groups_Top15.pdf.)
################################################################################

message("GO SLIM enrichment (uIF + %uORF+ groupings): starting...")

.require_pkg("GO.db")
.require_pkg("AnnotationDbi")
.require_pkg("stringr")

stopifnot(exists("mORF_uORF_Info"))
stopifnot(exists("OUT_df6"))

goslim_path <- path.expand(PATH_GOSLIM)  # declared in section 0
if (!file.exists(goslim_path)) stop("GO Slim file not found: ", goslim_path)

raw_go <- read.delim(goslim_path, header = FALSE, stringsAsFactors = FALSE, sep = "\t", quote = "")

gene_col_idx <- which.max(apply(raw_go[1:100, ], 2, function(x) mean(grepl("^AT[1-5MC]G", x, ignore.case = TRUE))))
go_col_idx   <- which.max(apply(raw_go[1:100, ], 2, function(x) mean(grepl("GO:", x))))

go_wide <- data.frame(
  gene_raw  = trimws(raw_go[, gene_col_idx]),
  go_string = trimws(raw_go[, go_col_idx]),
  stringsAsFactors = FALSE
)

go_long <- go_wide %>%
  tidyr::separate_rows(go_string, sep = ",") %>%
  dplyr::mutate(
    gene  = toupper(gsub("\\.\\d+$", "", gene_raw)),
    go_id = trimws(go_string)
  ) %>%
  dplyr::select(gene, go_id) %>%
  dplyr::filter(grepl("^GO:", go_id)) %>%
  dplyr::distinct()

go_keys <- unique(go_long$go_id)
go_info <- AnnotationDbi::select(GO.db::GO.db, keys = go_keys, columns = c("TERM", "ONTOLOGY"), keytype = "GOID")

go_full <- dplyr::left_join(go_long, go_info, by = c("go_id" = "GOID")) %>%
  dplyr::filter(!is.na(ONTOLOGY))

message("GO associations loaded: ", nrow(go_full))

df_all_uorfs_univ <- mORF_uORF_Info %>%
  dplyr::mutate(
    gene_id_u = toupper(substr(as.character(tx_id), 1, 9)),
    uIF = uORF_TL / mORF_TL
  ) %>%
  dplyr::filter(!is.na(uIF) & is.finite(uIF)) %>%
  dplyr::group_by(gene_id_u, uORF_start_bio, uORF_stop_bio) %>%
  dplyr::slice_max(order_by = uIF, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

all_uorf_genes <- unique(df_all_uorfs_univ$gene_id_u)
universe_genes <- intersect(all_uorf_genes, unique(go_full$gene))
universe_go <- go_full %>% dplyr::filter(gene %in% universe_genes)

message("Universe size (any-uORF genes with GO): ", length(universe_genes))
if (length(universe_genes) < 50) stop("Universe too small; check GO file / ID formats.")

OUT_df6_go <- OUT_df6 %>%
  dplyr::mutate(
    gene_id_go = toupper(if ("gene_id" %in% colnames(OUT_df6)) as.character(gene_id) else substr(as.character(tx_id), 1, 9)),
    uIF_num = as.numeric(uIF),
    pUP_num = as.numeric(pUP)      # %uORF+ (used for the second grouping)
  ) %>%
  dplyr::filter(!is.na(uIF_num) & is.finite(uIF_num)) %>%
  dplyr::select(gene_id_go, uIF_num, pUP_num) %>%
  dplyr::distinct()

## NOTE: boundary harmonized with Figures 1E/2C and Table S2
## (Weak: uIF <= 0.2; Mid: 0.2 < uIF <= 1). Previously uIF == 0.2 fell into Mid
## here (3 genes), making the GO group sizes 920/1075 instead of 923/1072.
genes_1weak   <- OUT_df6_go %>% dplyr::filter(uIF_num <= 0.2)                %>% dplyr::pull(gene_id_go) %>% unique()
genes_1mid    <- OUT_df6_go %>% dplyr::filter(uIF_num > 0.2 & uIF_num <= 1)  %>% dplyr::pull(gene_id_go) %>% unique()
genes_1strong <- OUT_df6_go %>% dplyr::filter(uIF_num > 1)                   %>% dplyr::pull(gene_id_go) %>% unique()

genes_1weak   <- intersect(genes_1weak, universe_genes)
genes_1mid    <- intersect(genes_1mid, universe_genes)
genes_1strong <- intersect(genes_1strong, universe_genes)

gene_summary <- df_all_uorfs_univ %>%
  dplyr::group_by(gene_id_u) %>%
  dplyr::summarise(
    n_uorfs = dplyr::n(),
    max_uIF = max(uIF, na.rm = TRUE),
    .groups = "drop"
  )

genes_group4 <- gene_summary %>%
  dplyr::filter(max_uIF > 1 | n_uorfs >= 2) %>%
  dplyr::pull(gene_id_u) %>%
  unique() %>%
  intersect(universe_genes)

GROUP4_NAME <- ">=1 strong uORF or >=2 uORFs"

## ---------------------------------------------------------------------------
## GROUPING (A): by uIF — TuORF translational strength   [original grouping]
## ---------------------------------------------------------------------------
group_lists_uIF <- list(
  "1-uORF Weak"   = genes_1weak,
  "1-uORF Mid"    = genes_1mid,
  "1-uORF Strong" = genes_1strong
)
group_lists_uIF[[GROUP4_NAME]] <- genes_group4
group_order_uIF <- c("1-uORF Weak", "1-uORF Mid", "1-uORF Strong", GROUP4_NAME)

## ---------------------------------------------------------------------------
## GROUPING (B): by %uORF+ (pUP) — TSS-determined uORF inclusion
## Same single-TuORF analysis set and the same three categories used in
## Figures 2B/2F (<30%, 30-80%, 80-100%), so group sizes match those panels
## (up to genes lacking any GO annotation, which the universe excludes).
## ---------------------------------------------------------------------------
genes_pup_low  <- OUT_df6_go %>% dplyr::filter(pUP_num <  30)                 %>% dplyr::pull(gene_id_go) %>% unique()
genes_pup_mid  <- OUT_df6_go %>% dplyr::filter(pUP_num >= 30 & pUP_num < 80)  %>% dplyr::pull(gene_id_go) %>% unique()
genes_pup_high <- OUT_df6_go %>% dplyr::filter(pUP_num >= 80)                 %>% dplyr::pull(gene_id_go) %>% unique()

genes_pup_low  <- intersect(genes_pup_low,  universe_genes)
genes_pup_mid  <- intersect(genes_pup_mid,  universe_genes)
genes_pup_high <- intersect(genes_pup_high, universe_genes)

group_lists_pUP <- list(
  "1-uORF 0~30% uORF+"   = genes_pup_low,
  "1-uORF 30~80% uORF+"  = genes_pup_mid,
  "1-uORF 80~100% uORF+" = genes_pup_high
)
## 4th group: the SAME gene set used as the 4th column of the uIF plots, so the
## two figure sets share an identical reference column and can be read side by
## side. Its enrichment results are therefore identical in both figures.
group_lists_pUP[[GROUP4_NAME]] <- genes_group4

group_order_pUP <- c("1-uORF 0~30% uORF+", "1-uORF 30~80% uORF+",
                     "1-uORF 80~100% uORF+", GROUP4_NAME)

message("Genes per group, grouping A (uIF), restricted to Universe:")
print(vapply(group_lists_uIF[group_order_uIF], length, integer(1)))
message("Genes per group, grouping B (%uORF+), restricted to Universe:")
print(vapply(group_lists_pUP[group_order_pUP], length, integer(1)))

run_fisher_enrichment <- function(group_genes, universe_go_df, all_bg_genes) {
  group_genes <- unique(group_genes)
  N <- length(all_bg_genes)
  n <- length(group_genes)
  if (n < 5) return(NULL)
  
  bg <- universe_go_df %>%
    dplyr::distinct(gene, go_id, TERM, ONTOLOGY) %>%
    dplyr::count(go_id, TERM, ONTOLOGY, name = "bg_count") %>%
    dplyr::filter(bg_count >= 3)
  
  if (nrow(bg) == 0) return(NULL)
  
  fg <- universe_go_df %>%
    dplyr::filter(gene %in% group_genes) %>%
    dplyr::distinct(gene, go_id) %>%
    dplyr::count(go_id, name = "fg_count")
  
  bg2 <- dplyr::left_join(bg, fg, by = "go_id") %>%
    dplyr::mutate(fg_count = dplyr::coalesce(fg_count, 0L)) %>%
    dplyr::filter(fg_count > 0)
  
  if (nrow(bg2) == 0) return(NULL)
  
  res <- lapply(seq_len(nrow(bg2)), function(i) {
    a <- as.integer(bg2$fg_count[i])
    bg_count <- as.integer(bg2$bg_count[i])
    
    b <- bg_count - a
    c <- n - a
    d <- (N - n) - b
    if (any(c(a, b, c, d) < 0)) return(NULL)
    
    ft <- stats::fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
    data.frame(
      GO_ID = bg2$go_id[i],
      Term = bg2$TERM[i],
      Ontology = bg2$ONTOLOGY[i],
      FG_Count = a,
      BG_Count = bg_count,
      p_value = ft$p.value,
      Odds_Ratio = unname(ft$estimate),
      stringsAsFactors = FALSE
    )
  })
  
  out <- dplyr::bind_rows(res)
  if (nrow(out) == 0) return(NULL)
  out
}

################################################################################
# GENERIC RUNNER: Fisher tests + BP/MF/CC bubble plots for
# ANY gene grouping. Called once per grouping (uIF, %uORF+), so both analyses
# use identical universe, statistics, term ordering and styling.
#   file_slug     -> becomes part of every output filename
#   grouping_desc -> printed in the plot subtitle so each PDF is self-identifying
################################################################################
run_go_for_grouping <- function(group_lists, group_order, file_slug, grouping_desc) {

  group_sizes <- vapply(group_lists[group_order], length, integer(1))
  message("--- GO enrichment, grouping: ", grouping_desc, " ---")

  all_results <- dplyr::bind_rows(lapply(group_order, function(grp) {
    res <- run_fisher_enrichment(group_lists[[grp]], universe_go, universe_genes)
    if (is.null(res) || nrow(res) == 0) return(NULL)
    res$Group <- grp
    res
  }))

  if (is.null(all_results) || nrow(all_results) == 0) {
    message("No enrichment results produced for grouping: ", grouping_desc)
    return(invisible(NULL))
  }

  ## Benjamini-Hochberg FDR within each group (reported in the CSV only; the
  ## plots keep ranking by raw p so they stay comparable to earlier versions)
  all_results <- all_results %>%
    dplyr::group_by(Group) %>%
    dplyr::mutate(p_adj_BH = stats::p.adjust(p_value, method = "BH")) %>%
    dplyr::ungroup()

  group_label_map <- setNames(
    paste0(group_order, "\n(n=", group_sizes[group_order], ")"),
    group_order
  )

  make_one_go_plot <- function(ont_code, ont_title) {

    plot_data <- all_results %>%
      dplyr::filter(Ontology == ont_code) %>%
      dplyr::group_by(Group) %>%
      dplyr::slice_min(order_by = p_value, n = 15, with_ties = FALSE) %>%
      dplyr::ungroup()

    if (nrow(plot_data) == 0) return(NULL)

    plot_data$Group <- factor(plot_data$Group, levels = group_order)
    plot_data$Group_Label <- factor(
      as.character(plot_data$Group),
      levels = group_order,
      labels = group_label_map[group_order]
    )

    plot_data$Term_Label <- paste0(stringr::str_trunc(plot_data$Term, 55), " (", plot_data$GO_ID, ")")

    ## --- TERM ORDERING -------------------------------------------------------
    ## Sort bubbles into blocks that follow the x-axis groups left to right:
    ##   top block    = terms plotted ONLY in the leftmost group, best p first
    ##   next blocks  = terms plotted only in group 2, then only 3, then only 4
    ##   bottom block = terms plotted in >1 group ("overlapping"); 2-group overlaps
    ##                  first, the most widely shared terms at the very bottom
    ## (previous behavior: single ordering by best p across all groups)
    term_membership <- plot_data %>%
      dplyr::group_by(Term_Label) %>%
      dplyr::summarise(
        n_groups    = dplyr::n_distinct(Group),
        first_group = min(as.integer(Group)),
        best_p      = min(p_value),
        .groups     = "drop"
      ) %>%
      dplyr::mutate(
        block = dplyr::if_else(n_groups == 1L, first_group, length(group_order) + 1L)
      ) %>%
      dplyr::arrange(block, n_groups, first_group, best_p)   # top-to-bottom order

    ## ggplot draws the FIRST factor level at the BOTTOM of a discrete y-axis,
    ## so reverse the top-to-bottom order when setting levels:
    plot_data$Term_Label <- factor(plot_data$Term_Label, levels = rev(term_membership$Term_Label))

    ggplot2::ggplot(plot_data, ggplot2::aes(x = Group_Label, y = Term_Label)) +
      ggplot2::geom_point(ggplot2::aes(size = FG_Count, color = -log2(p_value))) +
      ggplot2::scale_color_gradient(low = "#56B4E9", high = "#D55E00", name = "-log2(p)") +
      ggplot2::scale_x_discrete(drop = FALSE) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = ont_title,
        subtitle = paste0("Genes grouped by ", grouping_desc,
                          "\nTop 15 GO Slim terms per group (Fisher's exact test)"),
        x = "", y = "", size = "Gene count"
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = ggplot2::element_text(size = 8),
        plot.title = ggplot2::element_text(face = "bold")
      )
  }

  n_grp <- length(group_order)
  for (ont in c("BP", "MF", "CC")) {
    ont_title <- switch(ont,
                        BP = "GO Enrichment: Biological Process (BP)",
                        MF = "GO Enrichment: Molecular Function (MF)",
                        CC = "GO Enrichment: Cellular Component (CC)")
    p <- make_one_go_plot(ont, ont_title)
    fn <- sprintf("GO_%s_by_%s_%dgroups_Top15.pdf", ont, file_slug, n_grp)
    if (!is.null(p)) {
      print(p)
      save_pdf(p, fn, w = 7, h = 8)
      message("  saved ", fn)
    } else {
      message("  ", ont, ": no terms.")
    }
  }

  ## full enrichment table for this grouping (all tested terms, not just top 15)
  csv_fn <- sprintf("GO_enrichment_results_by_%s.csv", file_slug)
  utils::write.csv(
    all_results %>%
      dplyr::mutate(Group = factor(Group, levels = group_order),
                    Group_size = group_sizes[as.character(Group)]) %>%
      dplyr::arrange(Group, Ontology, p_value) %>%
      dplyr::select(Group, Group_size, Ontology, GO_ID, Term,
                    FG_Count, BG_Count, Odds_Ratio, p_value, p_adj_BH),
    file.path(OUT_DIR, csv_fn), row.names = FALSE
  )
  message("  saved ", csv_fn, "  (", nrow(all_results), " group x term tests)")

  invisible(all_results)
}

## ---- run both groupings -----------------------------------------------------
go_res_uIF <- run_go_for_grouping(
  group_lists_uIF, group_order_uIF,
  file_slug     = "uIF",
  grouping_desc = "uORF impact factor (uIF): Weak <= 0.2 < Mid <= 1 < Strong"
)

go_res_pUP <- run_go_for_grouping(
  group_lists_pUP, group_order_pUP,
  file_slug     = "pct_uORF_plus",
  grouping_desc = paste0("%uORF+ (uORF inclusion, ",
                         ifelse(TSS_SOURCE == "CAGE", "CAGE-seq", "TSS-seq"), ")")
)

message("GO SLIM enrichment complete (both groupings).")







################################################################################
# FIGURE 9 — mRNA half-life by uIF groups (9A/B) and pUP groups (9C/D)
# UPDATED: 9B now uses log2(uIF) on x-axis (exactly like Figure 1D / 3D)
# ECDF groups remain on linear uIF (Strong >1, Mid >0.2) for biological consistency
################################################################################
################################################################################
# FIGURE 9 — mRNA half-life by uIF groups (9A/B/C) and pUP groups (9D/E/F)
# FINAL ADJUSTMENTS:
#   - Legends in A and D moved to RIGHT side (c(0.78, 0.25))
#   - Legend text reduced (size = 9) + smaller key size
#   - Figure saved at 7.5 × 5 inches (compact, publication-ready)
################################################################################
message("FIGURE 9: half-life analysis starting...")

# ------------------------------------------------------------------------------
# 1) Load & clean half-lives
# ------------------------------------------------------------------------------
if (!exists("half_lives")) {
  half_lives <- read.delim(
    file = PATH_HALFLIFE_TXT,   # declared in section 0
    header = TRUE, stringsAsFactors = FALSE, sep = "\t"
  )
}
half_lives <- half_lives %>%
  dplyr::mutate(
    gene_id = as.character(Gene.ID),
    half_life_mean = as.numeric(Mean)
  ) %>%
  dplyr::filter(!is.na(half_life_mean), half_life_mean <= 15, is.finite(half_life_mean))

# ------------------------------------------------------------------------------
# 2) uIF groups + half-life (for 9A/B/C)
# ------------------------------------------------------------------------------
df_uif_hl <- OUT_df6 %>%
  dplyr::inner_join(
    half_lives %>% dplyr::select(gene_id, half_life_mean),
    by = "gene_id"
  ) %>%
  dplyr::filter(uIF < 15) %>%
  dplyr::mutate(
    uIF_Group = dplyr::case_when(
      uIF > 1   ~ "Strong uORF (>1)",
      uIF > 0.2 ~ "Mid uORF (>0.2 - 1)",
      TRUE      ~ "Weak uORF (0 - 0.2)"
    )
  )

df_no_hl <- half_lives %>%
  dplyr::filter(!gene_id %in% unique(OUT_df6$gene_id)) %>%
  dplyr::transmute(gene_id, half_life_mean, Group_Label = "No TuORF")

# 9A: ECDF — legend adjusted (right side + smaller)
plot9A_data <- dplyr::bind_rows(
  df_uif_hl %>% dplyr::transmute(half_life_mean, Group_Label = uIF_Group),
  df_no_hl
) %>% dplyr::filter(is.finite(half_life_mean))

plot9A_data$Group_Label <- factor(
  plot9A_data$Group_Label,
  levels = c("Strong uORF (>1)", "Mid uORF (>0.2 - 1)", "Weak uORF (0 - 0.2)", "No TuORF")
)

stats9A <- generate_legend_labels(plot9A_data, "Group_Label", "half_life_mean")
plot9A_data <- plot9A_data %>%
  dplyr::left_join(stats9A %>% dplyr::select(Group_Label, Final_Label), by = "Group_Label")

plot9A_data$Final_Label <- factor(
  plot9A_data$Final_Label,
  levels = stats9A$Final_Label[order(match(stats9A$Group_Label, levels(plot9A_data$Group_Label)))]
)

pal_uif <- scales::hue_pal()(3)
cols9A <- c(pal_uif[1], pal_uif[2], pal_uif[3], "#999999")
names(cols9A) <- levels(plot9A_data$Final_Label)

p9A <- ggplot2::ggplot(plot9A_data, ggplot2::aes(x = half_life_mean, color = Final_Label)) +
  ggplot2::stat_ecdf(geom = "step", size = 1.2) +
  ggplot2::scale_color_manual(values = cols9A) +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "A", x = "mRNA half-life (h)", y = "Fraction", color = "") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 16),
    legend.position = c(0.78, 0.25),           # RIGHT side
    legend.background = ggplot2::element_blank(),
    legend.key.size = grid::unit(0.35, "cm"),  # smaller
    legend.text = ggplot2::element_text(size = 9)  # reduced font
  )

# 9B: log2(uIF) scatter — ALL genes
p9B <- ggplot2::ggplot(
  df_uif_hl %>% dplyr::filter(is.finite(uIF), uIF > 0),
  ggplot2::aes(x = log2(uIF), y = half_life_mean)
) +
  ggplot2::geom_point(color = "#3498db", alpha = 0.6, size = 1.5) +
  ggplot2::geom_smooth(method = "lm", color = "red", size = 1, se = FALSE) +
  ggpubr::stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 4) +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "B", x = bquote(log[2]~"(uIF)"), y = "mRNA half-life (h)") +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16))

# 9C: log2(uIF) scatter — MID + STRONG only
p9C <- ggplot2::ggplot(
  df_uif_hl %>% 
    dplyr::filter(uIF_Group %in% c("Mid uORF (>0.2 - 1)", "Strong uORF (>1)"),
                  is.finite(uIF), uIF > 0),
  ggplot2::aes(x = log2(uIF), y = half_life_mean)
) +
  ggplot2::geom_point(color = "#3498db", alpha = 0.6, size = 1.5) +
  ggplot2::geom_smooth(method = "lm", color = "red", size = 1, se = FALSE) +
  ggpubr::stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 4) +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "C: log2(uIF) vs half-life (mid + strong uORF genes only)",
                x = bquote(log[2]~"(uIF)"), y = "mRNA half-life (h)") +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16))

# ------------------------------------------------------------------------------
# 4) pUP groups + half-life (for 9D/E/F)
# ------------------------------------------------------------------------------
df_pup_hl <- OUT_df6 %>%
  dplyr::inner_join(
    half_lives %>% dplyr::select(gene_id, half_life_mean),
    by = "gene_id"
  ) %>%
  dplyr::mutate(
    pUP_Group = dplyr::case_when(
      pUP >= 80 ~ "80~100",
      pUP >= 30 ~ "30~80",
      TRUE      ~ "0~30"
    )
  )

# 9D: ECDF — legend adjusted (right side + smaller)
plot9D_data <- dplyr::bind_rows(
  df_pup_hl %>% dplyr::transmute(half_life_mean, Group_Label = pUP_Group),
  df_no_hl %>% dplyr::select(half_life_mean, Group_Label)
) %>% dplyr::filter(is.finite(half_life_mean))

plot9D_data$Group_Label <- factor(
  plot9D_data$Group_Label,
  levels = c("80~100", "30~80", "0~30", "No TuORF")
)

stats9D <- generate_legend_labels(plot9D_data, "Group_Label", "half_life_mean")
plot9D_data <- plot9D_data %>%
  dplyr::left_join(stats9D %>% dplyr::select(Group_Label, Final_Label), by = "Group_Label")

plot9D_data$Final_Label <- factor(
  plot9D_data$Final_Label,
  levels = stats9D$Final_Label[order(match(stats9D$Group_Label, levels(plot9D_data$Group_Label)))]
)

pal_pup <- scales::hue_pal()(3)
cols9D <- c(pal_pup[1], pal_pup[2], pal_pup[3], "#999999")
names(cols9D) <- levels(plot9D_data$Final_Label)

p9D <- ggplot2::ggplot(plot9D_data, ggplot2::aes(x = half_life_mean, color = Final_Label)) +
  ggplot2::stat_ecdf(geom = "step", size = 1.2) +
  ggplot2::scale_color_manual(values = cols9D) +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "D", x = "mRNA half-life (h)", y = "Fraction", color = "") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 16),
    legend.position = c(0.78, 0.25),           # RIGHT side
    legend.background = ggplot2::element_blank(),
    legend.key.size = grid::unit(0.35, "cm"),  # smaller
    legend.text = ggplot2::element_text(size = 9)  # reduced font
  )

# 9E: pUP scatter — ALL genes
p9E <- ggplot2::ggplot(
  df_pup_hl %>% dplyr::filter(is.finite(pUP)),
  ggplot2::aes(x = pUP, y = half_life_mean)
) +
  ggplot2::geom_point(color = "#3498db", alpha = 0.5, size = 1) +
  ggplot2::geom_smooth(method = "lm", color = "red", size = 1, se = FALSE) +
  ggpubr::stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "E", x = "%pUP (CAGE)", y = "mRNA half-life (h)") +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16))

# 9F: pUP scatter — MID + STRONG only
p9F <- ggplot2::ggplot(
  df_pup_hl %>% dplyr::filter(pUP_Group %in% c("30~80", "80~100")),
  ggplot2::aes(x = pUP, y = half_life_mean)
) +
  ggplot2::geom_point(color = "#3498db", alpha = 0.5, size = 1) +
  ggplot2::geom_smooth(method = "lm", color = "red", size = 1, se = FALSE) +
  ggpubr::stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "F: %pUP vs half-life (mid + strong uORF genes only)",
                x = "%pUP (CAGE)", y = "mRNA half-life (h)") +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16))

# ------------------------------------------------------------------------------
# 5) Combine + SAVE (compact 7.5 × 5 inches)
# ------------------------------------------------------------------------------
Figure9_ABCDEF <- ggpubr::ggarrange(
  p9A, p9B, p9C,
  p9D, p9E, p9F,
  ncol = 3, nrow = 2
)
print(Figure9_ABCDEF)
save_pdf(Figure9_ABCDEF, "Figure9_ABCDEF.pdf", w = 7.5, h = 5)

# === ANALYSIS: half-life means per group ===
cat("\n=== Half-life mean analysis (3 uIF groups + No TuORF) ===\n")
print(
  dplyr::bind_rows(
    df_uif_hl %>%
      dplyr::group_by(uIF_Group) %>%
      dplyr::summarise(n = dplyr::n(), mean_hl = mean(half_life_mean), median_hl = median(half_life_mean), .groups = "drop") %>%
      dplyr::rename(Group = uIF_Group),
    df_no_hl %>%
      dplyr::summarise(n = dplyr::n(), mean_hl = mean(half_life_mean), median_hl = median(half_life_mean)) %>%
      dplyr::mutate(Group = "No TuORF")
  )
)

cat("\n=== Half-life mean analysis (3 pUP groups + No TuORF) ===\n")
print(
  dplyr::bind_rows(
    df_pup_hl %>%
      dplyr::group_by(pUP_Group) %>%
      dplyr::summarise(n = dplyr::n(), mean_hl = mean(half_life_mean), median_hl = median(half_life_mean), .groups = "drop") %>%
      dplyr::rename(Group = pUP_Group),
    df_no_hl %>%
      dplyr::summarise(n = dplyr::n(), mean_hl = mean(half_life_mean), median_hl = median(half_life_mean)) %>%
      dplyr::mutate(Group = "No TuORF")
  )
)

message("FIGURE 9 complete. Saved: Figure9_ABCDEF.pdf (compact 7.5×5, legends right + smaller)")






################################################################################
# FIGURE 11: uORF ORDER ANALYSIS (log2 uIF) - UPDATED RANGES
################################################################################

message("FIGURE 11: starting...")

stopifnot(exists("mORF_uORF_Info"))
stopifnot(exists("OUT_pUP"))

df_all_uorfs_clean <- mORF_uORF_Info %>%
  dplyr::mutate(
    gene_id = substr(tx_id, 1, 9),
    uIF = uORF_TL / mORF_TL,
    uORF_start_bio = as.numeric(uORF_start_bio),
    uORF_stop_bio  = as.numeric(uORF_stop_bio)
  ) %>%
  dplyr::filter(!is.na(uIF) & is.finite(uIF) & uIF > 0) %>%
  dplyr::mutate(log2_uIF = log2(uIF)) %>%
  dplyr::group_by(gene_id, uORF_start_bio, uORF_stop_bio) %>%
  dplyr::slice_max(order_by = uIF, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

df_ranked_all <- df_all_uorfs_clean %>%
  dplyr::group_by(gene_id) %>%
  dplyr::arrange(
    ifelse(Strand == "+", uORF_start_bio, -uORF_start_bio),
    .by_group = TRUE
  ) %>%
  dplyr::mutate(
    rank_val = dplyr::row_number(),
    total_uorfs = dplyr::n()
  ) %>%
  dplyr::ungroup()

df_ranked_all$unique_key <- paste(df_ranked_all$tx_id, df_ranked_all$genome_range, sep="|")

if ("unique_key" %in% colnames(OUT_pUP)) {
  df_ranked_final <- dplyr::left_join(df_ranked_all, OUT_pUP, by="unique_key")
} else {
  message("Joining pUP by tx_id (fallback)...")
  OUT_pUP_temp <- OUT_pUP
  if (colnames(OUT_pUP_temp)[1] != "tx_id") colnames(OUT_pUP_temp)[1] <- "tx_id"
  df_ranked_final <- dplyr::left_join(df_ranked_all, OUT_pUP_temp, by="tx_id")
}

df_ranked_final$pUP[is.na(df_ranked_final$pUP)] <- 0

generate_uorf_order_fig <- function(data, n_uorfs, file_suffix) {
  
  df_subset <- data %>%
    dplyr::filter(total_uorfs == n_uorfs) %>%
    dplyr::mutate(
      Rank_Label = dplyr::case_when(
        rank_val == 1 ~ "1st",
        rank_val == 2 ~ "2nd",
        rank_val == 3 ~ "3rd",
        rank_val == 4 ~ "4th",
        TRUE ~ as.character(rank_val)
      )
    )
  
  possible_ranks <- c("1st", "2nd", "3rd", "4th", "5th")
  rank_levels <- possible_ranks[1:n_uorfs]
  df_subset$Rank_Label <- factor(df_subset$Rank_Label, levels = rank_levels)
  
  gene_count <- length(unique(df_subset$gene_id))
  plot_title_base <- paste0("(Genes w/ ", n_uorfs, " uORFs)")
  message("Generating Figure 11 for ", plot_title_base, " | n=", gene_count)
  
  if (gene_count < 5) { message("  - Too few genes. Skipping."); return(invisible(NULL)) }
  
  get_stat_letters <- function(d, val_col) {
    stat.test <- stats::pairwise.wilcox.test(d[[val_col]], d$Rank_Label, p.adjust.method = "BH", exact = FALSE)
    if (all(is.nan(stat.test$p.value)) || all(is.na(stat.test$p.value))) return(NULL)
    p_mat <- rcompanion::fullPTable(stat.test$p.value)
    lets <- multcompView::multcompLetters(p_mat, compare="<", threshold=0.05, Letters=letters, reversed=FALSE)
    data.frame(Rank_Label=names(lets$Letters), Letters=lets$Letters)
  }
  
  pA <- ggplot2::ggplot(df_subset, ggplot2::aes(x = log2_uIF, color = Rank_Label)) +
    ggplot2::stat_ecdf(geom = "step", size = 1.2) +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::coord_cartesian(xlim = c(-7, 8)) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste0("A: log2(uIF) ECDF ", plot_title_base, "\n(n = ", gene_count, " genes)"),
      x="log2(uIF)", y="Fraction", color="Order"
    ) +
    ggplot2::theme(plot.title = ggplot2::element_text(face="bold", size=10))
  
  stats_uif <- get_stat_letters(df_subset, "log2_uIF")
  pB <- ggplot2::ggplot(df_subset, ggplot2::aes(x = Rank_Label, y = log2_uIF, fill = Rank_Label)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    ggplot2::coord_cartesian(ylim = c(-7, 8)) +
    ggplot2::scale_fill_brewer(palette = "Set1") +
    ggplot2::theme_classic() +
    ggplot2::labs(title = "B: log2(uIF) Boxplot", x="", y="log2(uIF)") +
    ggplot2::theme(plot.title = ggplot2::element_text(face="bold", size=10), legend.position="none")
  
  if (!is.null(stats_uif)) {
    pB <- pB + ggplot2::geom_text(data=stats_uif, ggplot2::aes(x=Rank_Label, y=7, label=Letters),
                                  inherit.aes=FALSE, size=5)
  }
  
  pC <- ggplot2::ggplot(df_subset, ggplot2::aes(x = Rank_Label, y = log2_uIF, fill = Rank_Label)) +
    ggplot2::geom_violin(trim = TRUE, alpha = 0.7, scale = "width") +
    ggplot2::geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    ggplot2::coord_cartesian(ylim = c(-7, 8)) +
    ggplot2::scale_fill_brewer(palette = "Set1") +
    ggplot2::theme_classic() +
    ggplot2::labs(title = "C: log2(uIF) Distribution", x="", y="log2(uIF)") +
    ggplot2::theme(plot.title = ggplot2::element_text(face="bold", size=10), legend.position="none")
  
  pD <- ggplot2::ggplot(df_subset, ggplot2::aes(x = pUP, color = Rank_Label)) +
    ggplot2::stat_ecdf(geom = "step", size = 1.2) +
    ggplot2::scale_color_brewer(palette = "Dark2") +
    ggplot2::theme_classic() +
    ggplot2::labs(title = paste("D: pUP ECDF", plot_title_base), x="% pUP", y="Fraction", color="Order") +
    ggplot2::theme(plot.title = ggplot2::element_text(face="bold", size=10))
  
  stats_pup <- get_stat_letters(df_subset, "pUP")
  pE <- ggplot2::ggplot(df_subset, ggplot2::aes(x = Rank_Label, y = pUP, fill = Rank_Label)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    ggplot2::coord_cartesian(ylim = c(0, 110)) +
    ggplot2::scale_fill_brewer(palette = "Dark2") +
    ggplot2::theme_classic() +
    ggplot2::labs(title = "E: pUP Boxplot", x="", y="% pUP") +
    ggplot2::theme(plot.title = ggplot2::element_text(face="bold", size=10), legend.position="none")
  
  if (!is.null(stats_pup)) {
    pE <- pE + ggplot2::geom_text(data=stats_pup, ggplot2::aes(x=Rank_Label, y=105, label=Letters),
                                  inherit.aes=FALSE, size=5)
  }
  
  pF <- ggplot2::ggplot(df_subset, ggplot2::aes(x = Rank_Label, y = pUP, fill = Rank_Label)) +
    ggplot2::geom_violin(trim = TRUE, alpha = 0.7, scale = "width") +
    ggplot2::geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    ggplot2::scale_fill_brewer(palette = "Dark2") +
    ggplot2::theme_classic() +
    ggplot2::labs(title = "F: pUP Distribution", x="", y="% pUP") +
    ggplot2::theme(plot.title = ggplot2::element_text(face="bold", size=10), legend.position="none")
  
  Fig_combined <- ggpubr::ggarrange(pA, pB, pC, pD, pE, pF, ncol=3, nrow=2, labels="AUTO")
  print(Fig_combined)
  
  final_filename <- paste0("Figure11_log2uIF_", file_suffix, ".pdf")
  save_pdf(Fig_combined, final_filename, w = 15, h = 10)
  message("Saved: ", final_filename)
}

generate_uorf_order_fig(df_ranked_final, 2, "2_uORFs")
generate_uorf_order_fig(df_ranked_final, 3, "3_uORFs")
generate_uorf_order_fig(df_ranked_final, 4, "4_uORFs")

# ------------------------------------------------------------------------------
# %uORF+ BY uORF ORDER — summary numbers behind Figure S3
# For genes with exactly N TuORFs (N = 2, 3, 4; the same uIF>0-filtered set the
# Figure 11 / Figure S3 panels use), summarize %uORF+ (pUP) at each positional
# rank (1st = 5'-most / mORF-distal ... Nth = mORF-proximal), and run a PAIRED
# Wilcoxon test of 1st vs last rank within the same genes.
# Writes Figure11_pUP_by_uORF_rank_summary.csv to OUT_DIR.
# ------------------------------------------------------------------------------
rank_pup_summary <- dplyr::bind_rows(lapply(2:4, function(N) {
  df_ranked_final %>%
    dplyr::filter(total_uorfs == N) %>%
    dplyr::group_by(rank_val) %>%
    dplyr::summarise(
      n_genes         = dplyr::n(),
      median_pUP      = round(stats::median(pUP, na.rm = TRUE), 1),
      mean_pUP        = round(mean(pUP, na.rm = TRUE), 1),
      pct_low_0_30    = round(100 * mean(pUP <  30, na.rm = TRUE), 1),
      pct_mid_30_80   = round(100 * mean(pUP >= 30 & pUP < 80, na.rm = TRUE), 1),
      pct_high_80_100 = round(100 * mean(pUP >= 80, na.rm = TRUE), 1),
      .groups = "drop"
    ) %>%
    dplyr::mutate(n_uORFs_in_gene = N, .before = 1)
}))

## paired comparison within genes: 1st (5'-most) vs last (mORF-proximal) rank
rank_pup_pairs <- dplyr::bind_rows(lapply(2:4, function(N) {
  w <- df_ranked_final %>%
    dplyr::filter(total_uorfs == N, rank_val %in% c(1L, N)) %>%
    dplyr::mutate(which = dplyr::if_else(rank_val == 1L, "first", "last")) %>%
    dplyr::select(gene_id, which, pUP) %>%
    tidyr::pivot_wider(names_from = which, values_from = pUP)
  wt <- stats::wilcox.test(w$first, w$last, paired = TRUE)
  data.frame(
    n_uORFs_in_gene   = N,
    n_genes           = nrow(w),
    median_first_pUP  = round(stats::median(w$first), 1),
    median_last_pUP   = round(stats::median(w$last), 1),
    pct_first_lt_last = round(100 * mean(w$first < w$last), 1),
    paired_wilcoxon_p = signif(wt$p.value, 3)
  )
}))

for (i in seq_len(nrow(rank_pup_pairs))) {
  message(sprintf(
    "%d-uORF genes (n=%d): median %%uORF+ 1st (5'-most) = %.1f vs last (mORF-proximal) = %.1f; 1st < last in %.1f%% of genes; paired Wilcoxon p = %.2e",
    rank_pup_pairs$n_uORFs_in_gene[i], rank_pup_pairs$n_genes[i],
    rank_pup_pairs$median_first_pUP[i], rank_pup_pairs$median_last_pUP[i],
    rank_pup_pairs$pct_first_lt_last[i], rank_pup_pairs$paired_wilcoxon_p[i]))
}

utils::write.csv(rank_pup_summary,
                 file.path(OUT_DIR, "Figure11_pUP_by_uORF_rank_summary.csv"), row.names = FALSE)
utils::write.csv(rank_pup_pairs,
                 file.path(OUT_DIR, "Figure11_pUP_first_vs_last_paired.csv"), row.names = FALSE)
message("Saved Figure11_pUP_by_uORF_rank_summary.csv and Figure11_pUP_first_vs_last_paired.csv")

message("FIGURE 11 complete.")

################################################################################
# NOTE:
# Remaining modules: transcriptome impact, length bins, Kozak, correlation
# were already largely consistent with Figures 1–5, but had repeated library blocks
# and occasional redundant helper definitions. If you want, paste the remainder you
# actually still use (after Figure 11) and I’ll collapse it the same way (no dup libs,
# strict PATH_* vars, save_pdf everywhere) in one clean block.
################################################################################









# Table 1. %transcripts affected by uORFs
# TRANSCRIPTOME-WIDE uORF IMPACT (MULTI-ISOFORM AWARE)
# - Expression axis uses log2
################################################################################
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(knitr)
  library(scales)
  library(readr)
  library(GenomicFeatures)
  library(GenomicRanges)
})

message(">>> Starting Transcript-Level Impact Analysis (Multi-Isoform Aware) <<<")

# ==============================================================================
# 1. LOAD & PROCESS RSEM ISOFORM DATA
# ==============================================================================
message("Reading per-isoform RSEM TPM table (6 samples in one file)...")

# One table with columns: transcript_id, gene_id, TPM_D1 ... TPM_D7
# (merged from the six RSEM *.isoforms.results files; TPM columns only)
df_rsem_iso <- readr::read_tsv(PATH_RSEM_ISOFORMS, show_col_types = FALSE)

# Average TPM across the 6 samples per transcript
df_transcripts <- df_rsem_iso %>%
  dplyr::mutate(TPM = rowMeans(dplyr::across(dplyr::starts_with("TPM_")), na.rm = TRUE)) %>%
  dplyr::select(transcript_id, gene_id, TPM)

message(paste("Total Unique Transcripts Loaded:", nrow(df_transcripts)))

# ==============================================================================
# 2. LOAD ANNOTATION (Required for Coordinate Matching)
# ==============================================================================
message("Loading GTF to map transcripts to genome coordinates...")
if (!exists("txdb")) {
  txdb <- GenomicFeatures::makeTxDbFromGFF(
    file = PATH_GTF, format = "gtf",
    dataSource = "Araport11", organism = "Arabidopsis"
  )
}

# Get Exons grouped by Transcript
exons_by_tx <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)

# Filter for transcripts actually present in RSEM data
common_tx <- intersect(names(exons_by_tx), df_transcripts$transcript_id)
exons_target <- exons_by_tx[common_tx]

# ==============================================================================
# 3. PREPARE uORF DATA & EXPAND TO ISOFORMS
# ==============================================================================
if (!exists("mORF_uORF_Info")) stop("mORF_uORF_Info missing. Run Step 2 first.")
if (!exists("OUT_pUP")) stop("OUT_pUP missing. Run Step 3 first.")

# ==============================================================================
# Multi-isoform mapping note: uORFs are assigned to an isoform by requiring
# BOTH uORF endpoints to fall within that isoform's exons (spliced-safe);
# a simple bounding-box "within" test would under-map spliced uORFs.
#
# ASSUMES ALREADY EXISTS IN THIS MODULE:
#   - df_transcripts  (avg TPM per transcript)
#   - exons_target    (exonsBy(txdb) filtered to transcripts present in df_transcripts)
#   - mORF_uORF_Info
#   - OUT_pUP
#   - fiveUTRByGeneNames (or PATH_5UTR_GENES)
#   - OUT_DIR + save_pdf()
# ==============================================================================

# ---- The impact analysis below runs for BOTH TSS sources ----
# The whole computation below is now wrapped in run_transcriptome_impact(), so
# the same table/plots are produced once with CAGE-seq pUP and once with
# TSS-seq pUP (Nielsen et al. 2019), saved separately:
#   Table_Transcriptome_Impact_MultiIsoform_CAGE.csv    / _TSSseq.csv
#   Figure_Transcriptome_Impact_MultiIsoform_CAGE.pdf   / _TSSseq.pdf
# The table now also reports BOTH denominator conventions:
#   (a) nuclear transcripts with annotated 5'UTR  -> manuscript Figure 3C (21.6% / 9.6%)
#   (b) ALL nucleus-encoded transcripts (with + without annotated 5'UTR)
# plus expressed-gene counts backing the "~20% of expressed genes" sentence.

if (!exists(".compute_pUP_from_source"))
  stop(".compute_pUP_from_source missing. Run the Table S1 module helper first.")
if (!exists("fiveUTRByGeneNames")) {
  if (file.exists(PATH_5UTR_GENES)) {
    fiveUTRByGeneNames <- read.delim(PATH_5UTR_GENES, header = TRUE)$gene_id
  } else stop("fiveUTRByGeneNames missing and PATH_5UTR_GENES not found.")
}
## pUP tables for both sources (reuse the Table S1 module's if already computed)
if (!exists("pUP_CAGE_tbl"))
  pUP_CAGE_tbl <- if (TSS_SOURCE == "CAGE") OUT_pUP else .compute_pUP_from_source(PATH_CAGE_CAGEseq)
if (!exists("pUP_TSS_tbl"))
  pUP_TSS_tbl  <- if (TSS_SOURCE == "TSS")  OUT_pUP else .compute_pUP_from_source(PATH_CAGE_TSSseq)

run_transcriptome_impact <- function(pUP_tbl, src_tag) {

  message(">>> Transcriptome-wide TuORF impact [", src_tag, "] <<<")

  # ---- unique uORFs + pUP per genomic location ----
  df_uorf_unique <- mORF_uORF_Info %>%
    dplyr::select(uORF_start_bio, uORF_stop_bio, chr, Strand, tx_id, genome_range) %>%
    dplyr::mutate(unique_key = paste(tx_id, genome_range, sep = "|")) %>%
    dplyr::left_join(
      if ("unique_key" %in% colnames(pUP_tbl)) pUP_tbl else dplyr::rename(pUP_tbl, unique_key = tx_id),
      by = "unique_key"
    ) %>%
    dplyr::mutate(
      g_start = pmin(uORF_start_bio, uORF_stop_bio),
      g_end   = pmax(uORF_start_bio, uORF_stop_bio)
    ) %>%
    dplyr::filter(!is.na(g_start) & !is.na(g_end)) %>%
    dplyr::group_by(chr, g_start, g_end, Strand) %>%
    dplyr::summarise(pUP = max(pUP, na.rm = TRUE), .groups = "drop")

  # ---- endpoint-in-exon mapping (spliced-safe) ----
  gr_uorf_start <- GenomicRanges::GRanges(
    seqnames = df_uorf_unique$chr,
    ranges   = IRanges::IRanges(start = df_uorf_unique$g_start, width = 1),
    strand   = df_uorf_unique$Strand
  )
  gr_uorf_end <- GenomicRanges::GRanges(
    seqnames = df_uorf_unique$chr,
    ranges   = IRanges::IRanges(start = df_uorf_unique$g_end, width = 1),
    strand   = df_uorf_unique$Strand
  )
  hits_s <- GenomicRanges::findOverlaps(gr_uorf_start, exons_target, type = "within")
  hits_e <- GenomicRanges::findOverlaps(gr_uorf_end,   exons_target, type = "within")
  df_map_s <- data.frame(uorf_idx = S4Vectors::queryHits(hits_s),
                         transcript_id = names(exons_target)[S4Vectors::subjectHits(hits_s)],
                         stringsAsFactors = FALSE)
  df_map_e <- data.frame(uorf_idx = S4Vectors::queryHits(hits_e),
                         transcript_id = names(exons_target)[S4Vectors::subjectHits(hits_e)],
                         stringsAsFactors = FALSE)
  df_mapping <- unique(merge(df_map_s, df_map_e, by = c("uorf_idx", "transcript_id")))
  df_mapping$pUP <- df_uorf_unique$pUP[df_mapping$uorf_idx]

  # best uORF per transcript (max pUP = fraction of transcripts carrying >=1 uORF)
  df_best_uorf_expanded <- df_mapping %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::slice_max(order_by = pUP, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::rename(Max_pUP = pUP)

  message("[", src_tag, "] transcripts mapped to uORFs (endpoint-safe): ",
          nrow(df_best_uorf_expanded))

  # ---- per-transcript TPM, Max_pUP, category ----
  df_calc_tx <- df_transcripts %>%
    dplyr::left_join(df_best_uorf_expanded, by = "transcript_id") %>%
    dplyr::mutate(
      Has_uORF      = !is.na(Max_pUP),   # transcript carries a mapped TuORF (annotation view)
      Max_pUP       = ifelse(is.na(Max_pUP), 0, Max_pUP),
      Has_uORF_pup0 = Max_pUP > 0,       # original definition behind the ms 21.6% (requires pUP > 0)
      Affected_TPM  = TPM * (Max_pUP / 100),
      Category = dplyr::case_when(
        grepl("^ATCG", gene_id, ignore.case = TRUE) ~ "Chloroplast",
        grepl("^ATMG", gene_id, ignore.case = TRUE) ~ "Mitochondrial",
        gene_id %in% fiveUTRByGeneNames ~ "Nuclear_with_5UTR",
        TRUE ~ "Nuclear_no_5UTR"
      )
    )

  tpm_sums <- df_calc_tx %>%
    dplyr::group_by(Category) %>%
    dplyr::summarise(Total_TPM = sum(TPM, na.rm = TRUE)) %>%
    dplyr::ungroup()
  .val <- function(cat) { v <- tpm_sums$Total_TPM[tpm_sums$Category == cat]
                          if (length(v) == 0) 0 else v }
  val_chl     <- .val("Chloroplast")
  val_mito    <- .val("Mitochondrial")
  val_nuc_utr <- .val("Nuclear_with_5UTR")
  val_nuc_no  <- .val("Nuclear_no_5UTR")
  val_nuc_all <- val_nuc_utr + val_nuc_no

  df_nuc_tx <- df_calc_tx %>%
    dplyr::filter(Category %in% c("Nuclear_with_5UTR", "Nuclear_no_5UTR"))
  total_affected_tpm        <- sum(df_nuc_tx$Affected_TPM, na.rm = TRUE)
  total_uorf_containing_tpm <- sum(df_nuc_tx$TPM[df_nuc_tx$Has_uORF], na.rm = TRUE)
  total_uorf_pup0_tpm       <- sum(df_nuc_tx$TPM[df_nuc_tx$Has_uORF_pup0], na.rm = TRUE)

  pct1 <- function(num, den, digits = 2)
    ifelse(den > 0, paste0(round(num / den * 100, digits), "%"), "NA")
  fmt_n_pct <- function(n, denom, digits = 2) {
    if (is.na(denom) || denom == 0) return("NA")
    paste0(round(n / denom * 100, digits), "% (n=", n, "/", denom, ")")
  }

  # ---- transcript / gene counts (5'UTR-containing nuclear set) ----
  df_tx_5utr <- df_calc_tx %>% dplyr::filter(Category == "Nuclear_with_5UTR")
  n_tx_5utr_total    <- nrow(df_tx_5utr)
  n_tx_5utr_uorf     <- sum(df_tx_5utr$Max_pUP > 0,  na.rm = TRUE)
  n_tx_5utr_uorf_p30 <- sum(df_tx_5utr$Max_pUP > 30, na.rm = TRUE)
  df_gene_5utr <- df_tx_5utr %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(Has_uORF     = any(Max_pUP > 0,  na.rm = TRUE),
                     Has_uORF_p30 = any(Max_pUP > 30, na.rm = TRUE), .groups = "drop")
  n_gene_5utr_total    <- nrow(df_gene_5utr)
  n_gene_5utr_uorf     <- sum(df_gene_5utr$Has_uORF)
  n_gene_5utr_uorf_p30 <- sum(df_gene_5utr$Has_uORF_p30)

  # ---- expressed-gene context (manuscript: "~20% of expressed genes") ----
  df_gene_all <- df_calc_tx %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(TPM = sum(TPM, na.rm = TRUE),
                     Has_TuORF_tx = any(Has_uORF), .groups = "drop")
  n_genes_expressed  <- sum(df_gene_all$TPM > 0, na.rm = TRUE)
  n_genes_TuORF_expr <- sum(df_gene_all$Has_TuORF_tx & df_gene_all$TPM > 0, na.rm = TRUE)

  # ---- expressed-gene count as in the manuscript (gene-level RSEM, ref 24) ----
  n_expr_ms <- NA_integer_; n_tu_expr_ms <- NA_integer_; n_tu_universe <- NA_integer_
  try({
    rna_g <- read.delim(path.expand(PATH_RSEM_RNA), header = TRUE, stringsAsFactors = FALSE)
    gid_expr <- unique(rna_g$gene_id[!is.na(rna_g$TPM) & rna_g$TPM > 0])
    n_expr_ms <- length(gid_expr)
    tu_genes <- unique(substr(as.character(OUT_df5_clean$gene_id), 1, 9))
    n_tu_universe <- length(tu_genes)
    n_tu_expr_ms <- length(intersect(tu_genes, gid_expr))
  }, silent = TRUE)

  results_table <- data.frame(
    Metric = c(
      "0. TSS source for %uORF+ (pUP)",
      "1. Chloroplast Transcripts (TPM)",
      "2. Mitochondrial Transcripts (TPM)",
      "3. Nuclear Transcripts with 5' UTR (TPM)",
      "4. Nuclear Transcripts without 5' UTR (TPM)",
      "5. All Nuclear Transcripts (TPM; row 3 + row 4)",
      "----------------------------------",
      "6. Nuclear transcripts with 5' UTR (N)",
      "   Nuclear 5'UTR transcripts with uORFs (pUP>0)",
      "   Nuclear 5'UTR transcripts with uORFs (pUP>30%)",
      "7. 5'UTR-containing nuclear genes (N)",
      "   5'UTR genes regulated by uORFs (pUP>0)",
      "   5'UTR genes with uORFs (pUP>30%)",
      "8. Genes with RNA-seq TPM > 0 (2-isoform RSEM set)",
      "   TuORF-containing genes among them (mapped TuORF)",
      "   % TuORF genes of expressed genes (2-isoform set)",
      "9. Genes with RNA-seq TPM > 0 (gene-level RSEM, as in ms)  [ms: 23,862]",
      "   TuORF gene universe (all TuORF genes)  [ms: 4,775]",
      "   TuORF genes among expressed  [ms: ~20% of expressed genes]",
      "----------------------------------",
      "A. Total TPM of Transcripts with mapped TuORFs (Unscaled, any pUP)",
      "   % Mass Containing TuORFs (of 5' UTR Transcripts)",
      "   % Mass Containing TuORFs (of ALL Nuclear Transcripts)",
      "A2. Total TPM of Transcripts with TuORFs and pUP > 0 (original definition)",
      "   % Mass (of 5' UTR Transcripts)  [ms Figure 3C: 21.6%]",
      "   % Mass (of ALL Nuclear Transcripts)",
      "B. Total TPM Affected by TuORFs (Scaled by Max pUP)",
      "   % Mass with TuORF(s) after TSS scaling (of 5' UTR Transcripts)  [ms Figure 3C: 9.6%]",
      "   % Mass with TuORF(s) after TSS scaling (of ALL Nuclear Transcripts)"
    ),
    Value = c(
      ifelse(src_tag == "CAGE", "CAGE-seq", "TSS-seq (Nielsen et al. 2019)"),
      round(val_chl, 2),
      round(val_mito, 2),
      round(val_nuc_utr, 2),
      round(val_nuc_no, 2),
      round(val_nuc_all, 2),
      "----------------",
      n_tx_5utr_total,
      fmt_n_pct(n_tx_5utr_uorf,     n_tx_5utr_total),
      fmt_n_pct(n_tx_5utr_uorf_p30, n_tx_5utr_total),
      n_gene_5utr_total,
      fmt_n_pct(n_gene_5utr_uorf,     n_gene_5utr_total),
      fmt_n_pct(n_gene_5utr_uorf_p30, n_gene_5utr_total),
      n_genes_expressed,
      n_genes_TuORF_expr,
      pct1(n_genes_TuORF_expr, n_genes_expressed),
      ifelse(is.na(n_expr_ms), "NA (RSEM gene file not read)", n_expr_ms),
      ifelse(is.na(n_tu_universe), "NA", n_tu_universe),
      ifelse(is.na(n_tu_expr_ms) | is.na(n_expr_ms) | n_expr_ms == 0, "NA",
             paste0(n_tu_expr_ms, " (", round(n_tu_expr_ms / n_expr_ms * 100, 2), "%)")),
      "----------------",
      round(total_uorf_containing_tpm, 2),
      pct1(total_uorf_containing_tpm, val_nuc_utr),
      pct1(total_uorf_containing_tpm, val_nuc_all),
      round(total_uorf_pup0_tpm, 2),
      pct1(total_uorf_pup0_tpm, val_nuc_utr),
      pct1(total_uorf_pup0_tpm, val_nuc_all),
      round(total_affected_tpm, 2),
      pct1(total_affected_tpm, val_nuc_utr),
      pct1(total_affected_tpm, val_nuc_all)
    ),
    stringsAsFactors = FALSE
  )

  print(knitr::kable(results_table,
        caption = paste0("Transcriptome-wide TuORF impact (multi-isoform corrected; TSS source = ",
                         src_tag, ")")))
  out_csv <- file.path(OUT_DIR, paste0("Table_Transcriptome_Impact_MultiIsoform_", src_tag, ".csv"))
  write.csv(results_table, out_csv, row.names = FALSE)
  message("[", src_tag, "] wrote ", out_csv)

  # ---- plots ----
  df_plot_data <- df_nuc_tx %>%
    dplyr::filter(Max_pUP > 0) %>%
    dplyr::mutate(pUP_Bin = cut(Max_pUP, breaks = seq(0, 100, by = 5),
                                include.lowest = TRUE, labels = seq(5, 100, 5)))

  pA <- ggplot2::ggplot(df_plot_data, ggplot2::aes(x = Max_pUP)) +
    ggplot2::geom_histogram(binwidth = 5, fill = "#2C3E50", color = "white", alpha = 0.9) +
    ggplot2::coord_cartesian(xlim = c(0, 100)) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = paste0("A: Transcripts Affected (", src_tag, ")"),
                  x = "% pUP (max per transcript)", y = "Count") +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  pB <- ggplot2::ggplot(df_plot_data, ggplot2::aes(x = pUP_Bin, y = TPM)) +
    ggplot2::geom_boxplot(fill = "#E67E22", outlier.size = 0.5, outlier.alpha = 0.3) +
    ggplot2::scale_y_continuous(trans = "log2", labels = scales::comma) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = paste0("B: Expression vs pUP (", src_tag, ")"),
                  x = "% pUP Bin", y = "TPM (log2)") +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))

  Fig_Impact_Multi <- ggpubr::ggarrange(pA, pB, ncol = 2, nrow = 1)
  print(Fig_Impact_Multi)
  save_pdf(Fig_Impact_Multi,
           paste0("Figure_Transcriptome_Impact_MultiIsoform_", src_tag, ".pdf"),
           w = 12, h = 6)

  invisible(list(table = results_table, per_transcript = df_calc_tx))
}

Impact_CAGE   <- run_transcriptome_impact(pUP_CAGE_tbl, "CAGE")
Impact_TSSseq <- run_transcriptome_impact(pUP_TSS_tbl,  "TSSseq")

message("Success! Multi-Isoform aware analysis complete for CAGE-seq AND TSS-seq.")





################################################################################
# uORF LENGTH BIN PLOTS (RT vs CiPS; "Both" included in BOTH groups)
# - Makes ONE bar plot:
#   (#) # of uORFs in each length bin (within each detection set)
#
# REQUIREMENTS:
#   - DF_All_Evidence exists (from the RT+CiPS merge section), with columns:
#       uORF_id, method, uORF_length (AA)
#   - OUT_DIR and save_pdf() exist; if not, this script defines safe defaults.
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

# ----------------------------
# 0) Output helpers (safe defaults)
# ----------------------------
if (!exists("OUT_DIR")) {
  OUT_DIR <- getwd()
  message("OUT_DIR not found; using working directory: ", OUT_DIR)
}
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

if (!exists("save_pdf")) {
  save_pdf <- function(plot_obj, filename, w = 8, h = 5) {
    ggplot2::ggsave(filename = file.path(OUT_DIR, filename), plot = plot_obj, width = w, height = h)
  }
}

# ----------------------------
# 1) Sanity checks
# ----------------------------
if (!exists("DF_All_Evidence")) stop("DF_All_Evidence missing. Run the RT+CiPS merge section first.")
need_cols <- c("uORF_id", "method", "uORF_length")
missing_cols <- setdiff(need_cols, colnames(DF_All_Evidence))
if (length(missing_cols) > 0) stop("DF_All_Evidence missing: ", paste(missing_cols, collapse = ", "))

# ----------------------------
# 2) Build length base + STANDARDIZE method labels
# ----------------------------
df_len_base <- DF_All_Evidence %>%
  dplyr::select(uORF_id, method, uORF_length) %>%
  dplyr::mutate(
    uORF_id = as.character(uORF_id),
    method_raw = as.character(method),
    method_raw = trimws(method_raw),
    method_std = dplyr::case_when(
      tolower(method_raw) %in% c("rt", "ribotaper", "ribo taper") ~ "RT",
      tolower(method_raw) %in% c("cips", "ci-ps")                ~ "CiPS",
      tolower(method_raw) %in% c("both")                         ~ "Both",
      TRUE                                                       ~ method_raw
    ),
    uORF_length = suppressWarnings(as.numeric(uORF_length))
  ) %>%
  dplyr::filter(!is.na(uORF_length), is.finite(uORF_length), uORF_length > 0)

# (Optional diagnostics)
message("=== DIAGNOSTIC: method_std counts (rows) ===")
print(table(df_len_base$method_std, useNA = "ifany"))

# ----------------------------
# 3) TWO inclusive detection sets (Both duplicated into each)
# ----------------------------
# Force method_std to a plain character vector
stopifnot("method_std" %in% names(df_len_base))

if (is.list(df_len_base$method_std)) {
  df_len_base$method_std <- vapply(df_len_base$method_std, function(z) as.character(z)[1], character(1))
} else {
  df_len_base$method_std <- as.character(df_len_base$method_std)
}
df_len_base$method_std <- trimws(df_len_base$method_std)

# Now filter will return TRUE/FALSE cleanly
df_len_two <- dplyr::bind_rows(
  df_len_base %>%
    dplyr::filter(.data$method_std == "RT" | .data$method_std == "Both") %>%
    dplyr::mutate(DetectedBy = "RiboTaper"),
  df_len_base %>%
    dplyr::filter(.data$method_std == "CiPS" | .data$method_std == "Both") %>%
    dplyr::mutate(DetectedBy = "CiPS")
) %>%
  dplyr::mutate(DetectedBy = factor(.data$DetectedBy, levels = c("RiboTaper", "CiPS")))
# Legend labels with unique-uORF counts
n_two <- df_len_two %>%
  dplyr::distinct(DetectedBy, uORF_id) %>%
  dplyr::count(DetectedBy, name = "n_uORFs")

label_map <- setNames(
  paste0(as.character(n_two$DetectedBy), "\n(n=", n_two$n_uORFs, ")"),
  as.character(n_two$DetectedBy)
)

# ----------------------------
# 4) Length bins + counts
# ----------------------------
bin_levels <- c("1", "2-10", "11-20", "21-30", ">30")

df_bins <- df_len_two %>%
  dplyr::mutate(
    uORF_length_int = as.integer(uORF_length),
    LenBin = dplyr::case_when(
      uORF_length_int == 1 ~ "1",
      uORF_length_int >= 2  & uORF_length_int <= 10 ~ "2-10",
      uORF_length_int >= 11 & uORF_length_int <= 20 ~ "11-20",
      uORF_length_int >= 21 & uORF_length_int <= 30 ~ "21-30",
      uORF_length_int > 30 ~ ">30",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(LenBin)) %>%
  dplyr::distinct(DetectedBy, uORF_id, LenBin) %>%     # count each uORF once per set per bin
  dplyr::count(DetectedBy, LenBin, name = "n_uORFs") %>%
  tidyr::complete(
    DetectedBy,
    LenBin = bin_levels,
    fill = list(n_uORFs = 0)
  ) %>%
  dplyr::mutate(
    DetectedBy_Label = factor(
      label_map[as.character(DetectedBy)],
      levels = unname(label_map[c("RiboTaper", "CiPS")])
    )
  )

df_bins$LenBin <- factor(df_bins$LenBin, levels = bin_levels)

# ----------------------------
# 5) Plot: counts (# uORFs)
# ----------------------------
p_bins_n <- ggplot2::ggplot(df_bins, ggplot2::aes(x = LenBin, y = n_uORFs, fill = DetectedBy_Label)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.75) +
  ggplot2::geom_text(
    ggplot2::aes(label = n_uORFs),
    position = ggplot2::position_dodge(width = 0.8),
    vjust = -0.25, size = 3
  ) +
  ggplot2::theme_classic() +
  ggplot2::labs(
    title = "uORF length bins: uORF counts in each bin\n(Both included in BOTH RT and CiPS)",
    x = "uORF length bin (AA)",
    y = "# uORFs (unique uORF_id)",
    fill = ""
  ) +
  ggplot2::theme(legend.position = "top") +
  ggplot2::scale_y_continuous(expand = expansion(mult = c(0, 0.10)))

print(p_bins_n)
save_pdf(p_bins_n, "LengthDist_LengthBins_Counts_withNinLegend.pdf", w = 8.5, h = 5)

message("Done: saved bin count plot to OUT_DIR = ", OUT_DIR)

# =============================================================================
# ADDITIONAL PLOT: 3 exclusive groups (CiPS only / RT only / Both)
#   - uses df_len_base (uORF_id, method_std, uORF_length)
#   - counts unique uORF_id once per bin within each group
#   - saves: LengthDist_LengthBins_Counts_3groups_Exclusive.pdf
# =============================================================================

stopifnot(exists("df_len_base"))
need_cols3 <- c("uORF_id", "method_std", "uORF_length")
miss3 <- setdiff(need_cols3, colnames(df_len_base))
if (length(miss3) > 0) stop("df_len_base missing: ", paste(miss3, collapse=", "))

# Ensure method_std is a plain character vector (handles list-columns)
if (is.list(df_len_base$method_std)) {
  df_len_base$method_std <- vapply(df_len_base$method_std, function(z) as.character(z)[1], character(1))
} else {
  df_len_base$method_std <- as.character(df_len_base$method_std)
}
df_len_base$method_std <- trimws(df_len_base$method_std)

# Keep only the 3 expected classes (RT / CiPS / Both)
df_len_3 <- df_len_base %>%
  dplyr::mutate(
    Group3 = dplyr::case_when(
      method_std == "CiPS" ~ "CiPS only",
      method_std == "RT"   ~ "RT only",
      method_std == "Both" ~ "Both",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(Group3)) %>%
  dplyr::mutate(
    uORF_length = suppressWarnings(as.numeric(uORF_length))
  ) %>%
  dplyr::filter(!is.na(uORF_length), is.finite(uORF_length), uORF_length > 0)

group_order3 <- c("CiPS only", "RT only", "Both")
df_len_3$Group3 <- factor(df_len_3$Group3, levels = group_order3)

# Legend labels with unique-uORF counts per group
n3 <- df_len_3 %>%
  dplyr::distinct(Group3, uORF_id) %>%
  dplyr::count(Group3, name = "n_uORFs")

label_map3 <- setNames(
  paste0(group_order3, "\n(n=", n3$n_uORFs[match(group_order3, as.character(n3$Group3))], ")"),
  group_order3
)

# Length bins (same binning as the length-bin module)
bin_levels <- c("1", "2-10", "11-20", "21-30", ">30")

df_bins_3 <- df_len_3 %>%
  dplyr::mutate(
    uORF_length_int = as.integer(uORF_length),
    LenBin = dplyr::case_when(
      uORF_length_int == 1 ~ "1",
      uORF_length_int >= 2  & uORF_length_int <= 10 ~ "2-10",
      uORF_length_int >= 11 & uORF_length_int <= 20 ~ "11-20",
      uORF_length_int >= 21 & uORF_length_int <= 30 ~ "21-30",
      uORF_length_int > 30 ~ ">30",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(LenBin)) %>%
  dplyr::distinct(Group3, uORF_id, LenBin) %>%   # each uORF counted once per bin per group
  dplyr::count(Group3, LenBin, name = "n_uORFs") %>%
  tidyr::complete(
    Group3 = factor(group_order3, levels = group_order3),
    LenBin = bin_levels,
    fill = list(n_uORFs = 0)
  ) %>%
  dplyr::mutate(
    LenBin = factor(LenBin, levels = bin_levels),
    Group3_Label = factor(label_map3[as.character(Group3)], levels = unname(label_map3[group_order3]))
  )

# Plot (dodged bars)
p_bins_n3 <- ggplot2::ggplot(df_bins_3, ggplot2::aes(x = LenBin, y = n_uORFs, fill = Group3_Label)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.75) +
  ggplot2::geom_text(
    ggplot2::aes(label = n_uORFs),
    position = ggplot2::position_dodge(width = 0.8),
    vjust = -0.25, size = 3
  ) +
  ggplot2::theme_classic() +
  ggplot2::labs(
    title = "uORF length bins: exclusive detection classes",
    subtitle = "CiPS only vs RT only vs Both (no duplication of 'Both')",
    x = "uORF length bin (AA)",
    y = "# uORFs (unique uORF_id)",
    fill = ""
  ) +
  ggplot2::theme(legend.position = "top") +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.10)))

print(p_bins_n3)
save_pdf(p_bins_n3, "LengthDist_LengthBins_Counts_3groups_Exclusive.pdf", w = 8.5, h = 5)






















################################################################################
# MODULE: KOZAK SCORE + ECDF + SEQLOGO  (FULL, FIXED)
# - uORFs: 3 uIF groups (Weak/Mid/Strong)  [OPTIONAL filters: pUP>=MIN_PUP, RNA>=MIN_RNA_TPM]
# - mORFs: 5 groups:
#     (1) mORF of Weak-uIF genes (single-uORF genes; pUP>=MIN_PUP)
#     (2) mORF of Mid-uIF genes  (single-uORF genes; pUP>=MIN_PUP)
#     (3) mORF of Strong-uIF genes (single-uORF genes; pUP>=MIN_PUP)
#     (4) mORF of multi-uORF genes
#     (5) mORF of no-uORF genes
#
# OUTPUTS (in OUT_DIR):
#   Kozak_Score_uORF_and_mORF_Violin.pdf
#   Kozak_ECDF_uORF_uIFgroups.pdf
#   Kozak_ECDF_mORF_groups.pdf
#   Kozak_SeqLogo_uORF_and_mORF.pdf
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(ggpubr)
  library(ggseqlogo)
  library(Biostrings)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
  library(Rsamtools)
  library(S4Vectors)
  library(rcompanion)
  library(multcompView)
  library(scales)
})

# ----------------------------
# 0) Checks / globals
# ----------------------------
if (!exists("OUT_DIR")) OUT_DIR <- getwd()
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
if (!exists("save_pdf")) {
  save_pdf <- function(plot_obj, filename, w = 8, h = 6) {
    ggplot2::ggsave(filename = file.path(OUT_DIR, filename), plot = plot_obj, width = w, height = h)
  }
}

if (!exists("PATH_GTF")) stop("PATH_GTF missing.")
if (!exists("PATH_FASTA_TAIR10")) stop("PATH_FASTA_TAIR10 missing.")
if (!exists("OUT_df6")) stop("OUT_df6 missing.")
if (!exists("mORF_uORF_Info")) stop("mORF_uORF_Info missing.")

# Optional but recommended
if (!exists("TE_CDS")) {
  message("NOTE: TE_CDS not found. 'no-uORF genes' will be defined from TxDb 5'UTR transcripts minus uORF genes.")
}

UIF_WEAK_MAX   <- 0.5
UIF_STRONG_MIN <- 1.5

MIN_PUP <- 30          # pUP>=30 filter for uORF uIF groups
MIN_RNA_TPM <- 0.1     # OPTIONAL RNA filter (only used if RNA info exists)
USE_RNA_FILTER <- TRUE # set FALSE to disable

SEQLOGO_YMAX <- 0.5    # requested
MAX_LOGO_SEQS <- 5000  # cap to keep runtime/memory stable

# ----------------------------
# 1) Build genome resources + transcript sequences
# ----------------------------
message("Kozak module: loading TxDb + FASTA...")
txdb_k <- GenomicFeatures::makeTxDbFromGFF(PATH_GTF, format = "gtf")

tryCatch({ Rsamtools::indexFa(PATH_FASTA_TAIR10) }, error = function(e) NULL)
fa_k <- Rsamtools::FaFile(PATH_FASTA_TAIR10)

tx_seqs_k     <- GenomicFeatures::extractTranscriptSeqs(fa_k, txdb_k, use.names = TRUE)
exons_by_tx_k <- GenomicFeatures::exonsBy(txdb_k, by = "tx", use.names = TRUE)
fiveUTRs_k    <- GenomicFeatures::fiveUTRsByTranscript(txdb_k, use.names = TRUE)

# ----------------------------
# 2) Reference Kozak matrix
#    - uses CDS start context from transcripts (5'UTR length defines CDS start)
# ----------------------------
message("Kozak module: building reference Kozak matrix...")
common_tx <- intersect(names(tx_seqs_k), names(fiveUTRs_k))

utr5_len_vec <- vapply(fiveUTRs_k[common_tx], function(gr) sum(IRanges::width(gr)), integer(1))
nt_seqs_ref  <- as.character(tx_seqs_k[common_tx])

ref_df <- data.frame(
  tx_id = common_tx,
  fUTR_length = as.integer(utr5_len_vec),
  nt_seq = nt_seqs_ref,
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    nt_len = nchar(nt_seq),
    ok = !is.na(fUTR_length) & fUTR_length >= 9 & (fUTR_length + 5) <= nt_len,
    Kozak = ifelse(ok, substr(nt_seq, fUTR_length - 8, fUTR_length + 5), NA_character_)
  ) %>%
  dplyr::filter(!is.na(Kozak), nchar(Kozak) == 14, substr(Kozak, 10, 12) == "ATG")

if (nrow(ref_df) < 200) {
  message("WARNING: reference Kozak set is small (n=", nrow(ref_df), "). Check TxDb/FASTA consistency.")
}

bases <- c("A","C","G","T")
full_counts <- Biostrings::consensusMatrix(Biostrings::DNAStringSet(ref_df$Kozak), as.prob = FALSE)[bases, , drop = FALSE] + 0.01
full_consensus <- prop.table(full_counts, 2)
letter_freq <- full_consensus[bases, , drop = FALSE]
aF <- Biostrings::alphabetFrequency(Biostrings::DNAStringSet(ref_df$Kozak), as.prob = TRUE, collapse = TRUE)[bases]
natlogratio <- log(letter_freq / aF)

# zero the ATG columns (10:12)
if (ncol(natlogratio) == 14) natlogratio[, 10:12] <- 0

# ----------------------------
# 3) Helpers
# ----------------------------
extract_last_two_ints <- function(x) {
  if (is.na(x)) return(c(NA_integer_, NA_integer_))
  nums <- stringr::str_extract_all(as.character(x), "\\d+")[[1]]
  if (length(nums) < 2) return(c(NA_integer_, NA_integer_))
  as.integer(tail(nums, 2))
}

order_exons_tx <- function(ex_gr) {
  if (length(ex_gr) == 0) return(ex_gr)
  st <- as.character(GenomicRanges::strand(ex_gr))[1]
  if (st == "-") ex_gr[order(GenomicRanges::start(ex_gr), decreasing = TRUE)]
  else           ex_gr[order(GenomicRanges::start(ex_gr), decreasing = FALSE)]
}

map_genomic_pos_to_tx <- function(ex_gr_ord, pos) {
  if (is.na(pos)) return(NA_integer_)
  s <- GenomicRanges::start(ex_gr_ord); e <- GenomicRanges::end(ex_gr_ord)
  idx <- which(pos >= s & pos <= e)
  if (length(idx) == 0) return(NA_integer_)
  idx <- idx[1]
  lens <- IRanges::width(ex_gr_ord)
  cum_before <- c(0L, cumsum(lens)[-length(lens)])
  st <- as.character(GenomicRanges::strand(ex_gr_ord))[1]
  if (st == "+") as.integer(cum_before[idx] + (pos - s[idx] + 1L))
  else           as.integer(cum_before[idx] + (e[idx] - pos + 1L))
}

compute_orf_tx_start <- function(tid, tx_range, gstart, gstop, strand_char, tx_seqs, exons_by_tx) {
  # 1) from Tx_range if possible
  ab <- extract_last_two_ints(tx_range)
  a <- ab[1]; b <- ab[2]
  if (!is.na(tid) && tid %in% names(tx_seqs)) {
    tx_len <- length(tx_seqs[[tid]])
    if (!is.na(a) && !is.na(b) && a <= tx_len && b <= tx_len) {
      return(min(a, b))
    }
  }
  
  # 2) map from genome coordinate
  if (is.na(tid) || !tid %in% names(exons_by_tx) || is.na(gstart) || is.na(gstop)) return(NA_integer_)
  ex_gr <- order_exons_tx(exons_by_tx[[tid]])
  st_g <- if (!is.na(strand_char) && nzchar(strand_char)) strand_char else as.character(GenomicRanges::strand(ex_gr))[1]
  g1 <- as.integer(gstart); g2 <- as.integer(gstop)
  start_g <- if (st_g == "-") max(g1, g2) else min(g1, g2)
  map_genomic_pos_to_tx(ex_gr, start_g)
}

calc_kozak <- function(seqs, tid, tx_start, matrix14) {
  if (is.na(tid) || !tid %in% names(seqs) || is.na(tx_start)) return(c(NA_character_, NA_real_))
  w_start <- tx_start - 9L
  w_end   <- tx_start + 4L
  if (w_start < 1L || w_end > length(seqs[[tid]])) return(c(NA_character_, NA_real_))
  k_seq <- as.character(Biostrings::subseq(seqs[[tid]], start = w_start, end = w_end))
  ch <- strsplit(k_seq, "")[[1]]
  score <- NA_real_
  if (length(ch) == 14 && all(ch %in% rownames(matrix14))) {
    score <- sum(diag(matrix14[ch, 1:14, drop = FALSE]))
  }
  c(k_seq, score)
}

# ----------------------------
# 4) KS-letter helper
# ----------------------------
get_pairwise_ks_pvalues_clean <- function(data, group_col, val_col, id_map) {
  groups <- levels(data[[group_col]])
  n_g <- length(groups)
  clean_names <- id_map[groups]
  p_mat <- matrix(1, nrow = n_g, ncol = n_g, dimnames = list(clean_names, clean_names))
  
  combos <- utils::combn(seq_len(n_g), 2, simplify = FALSE)
  
  p_values <- vapply(combos, function(idx) {
    g1 <- groups[idx[1]]
    g2 <- groups[idx[2]]
    d1 <- as.numeric(stats::na.omit(data[[val_col]][data[[group_col]] == g1]))
    d2 <- as.numeric(stats::na.omit(data[[val_col]][data[[group_col]] == g2]))
    if (length(d1) > 0 && length(d2) > 0) stats::ks.test(d1, d2)$p.value else 1
  }, numeric(1))
  
  idx_mat <- do.call(rbind, combos)
  i <- idx_mat[, 1]
  j <- idx_mat[, 2]
  
  p_mat[cbind(i, j)] <- p_values
  p_mat[cbind(j, i)] <- p_values
  p_mat
}

reorder_letters_by_appearance <- function(named_letters, ordered_groups) {
  current_sequence <- named_letters[ordered_groups]
  unique_codes <- unique(current_sequence)
  new_alphabet <- letters[1:length(unique_codes)]
  mapping <- setNames(new_alphabet, unique_codes)
  new_letters <- mapping[named_letters]
  names(new_letters) <- names(named_letters)
  new_letters
}

generate_legend_labels <- function(data_df, group_col, val_col) {
  group_levels <- levels(data_df[[group_col]])
  clean_ids <- paste0("G", seq_along(group_levels))
  names(clean_ids) <- group_levels
  
  p_mat <- get_pairwise_ks_pvalues_clean(data_df, group_col, val_col, clean_ids)
  
  p_vec <- p_mat[lower.tri(p_mat)]
  names_vec <- expand.grid(rownames(p_mat), colnames(p_mat))
  names_vec <- names_vec[lower.tri(p_mat), ]
  names(p_vec) <- paste(names_vec$Var1, names_vec$Var2, sep = "-")
  
  letters_list <- multcompView::multcompLetters(p_vec, compare="<", threshold=0.05, Letters=letters, reversed = FALSE)
  
  clean_letters <- letters_list$Letters
  clean_to_real <- setNames(names(clean_ids), clean_ids)
  initial_letters <- clean_letters
  names(initial_letters) <- clean_to_real[names(clean_letters)]
  
  final_letters <- reorder_letters_by_appearance(initial_letters, levels(data_df[[group_col]]))
  
  counts <- data_df %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")
  
  counts$letters <- final_letters[as.character(counts[[group_col]])]
  counts$Final_Label <- paste0(counts[[group_col]], " (n=", counts$n, ") ", counts$letters)
  counts
}

# ----------------------------
# 5) Build uORF set (3 uIF groups) + compute uORF Kozak
# ----------------------------
message("Kozak module: computing uORF Kozak for 3 uIF groups...")

df_uorf_base <- OUT_df6 %>%
  dplyr::mutate(pUP = suppressWarnings(as.numeric(pUP))) %>%
  dplyr::filter(!is.na(pUP), is.finite(pUP), pUP >= MIN_PUP) %>%
  dplyr::mutate(
    gene_id = as.character(gene_id),
    tx_id   = as.character(tx_id),
    uIF_val = dplyr::coalesce(
      suppressWarnings(as.numeric(.data[["uIF_raw"]])),
      suppressWarnings(as.numeric(.data[["uIF"]])),
      NA_real_
    )
  ) %>%
  dplyr::filter(!is.na(tx_id), !is.na(gene_id), !is.na(uIF_val), is.finite(uIF_val)) %>%
  dplyr::mutate(
    uIF_Group3 = dplyr::case_when(
      uIF_val > UIF_STRONG_MIN ~ "Strong",
      uIF_val >= UIF_WEAK_MAX  ~ "Mid",
      TRUE ~ "Weak"
    ),
    uIF_Group3 = factor(uIF_Group3, levels = c("Weak","Mid","Strong")),
    GroupLabel = factor(
      as.character(uIF_Group3),
      levels = c("Weak","Mid","Strong"),
      labels = c("uORF Weak (uIF<0.5)", "uORF Mid (0.5–1.5)", "uORF Strong (>1.5)")
    )
  )

# OPTIONAL RNA filter (only if RNA exists)
if (USE_RNA_FILTER) {
  if ("RNA" %in% names(df_uorf_base)) {
    df_uorf_base <- df_uorf_base %>% dplyr::filter(!is.na(RNA), is.finite(RNA), RNA >= MIN_RNA_TPM)
  } else if (exists("TE_CDS") && all(c("gene_id","RNA") %in% names(TE_CDS))) {
    te_rna <- TE_CDS %>% dplyr::select(gene_id, RNA) %>% dplyr::distinct()
    df_uorf_base <- df_uorf_base %>%
      dplyr::left_join(te_rna, by = "gene_id") %>%
      dplyr::filter(!is.na(RNA), is.finite(RNA), RNA >= MIN_RNA_TPM)
  } else {
    message("NOTE: RNA TPM filter requested but RNA not found in OUT_df6 or TE_CDS; skipping RNA filter.")
  }
}

# Ensure required fields exist (else rely on Tx_range mapping)
if (!"uORF_Tx_range" %in% names(df_uorf_base)) df_uorf_base$uORF_Tx_range <- NA_character_
if (!"ORF_gstart" %in% names(df_uorf_base)) df_uorf_base$ORF_gstart <- NA_integer_
if (!"ORF_gstop"  %in% names(df_uorf_base)) df_uorf_base$ORF_gstop  <- NA_integer_
if (!"Strand"     %in% names(df_uorf_base)) df_uorf_base$Strand     <- NA_character_

df_uorf_base <- df_uorf_base %>%
  dplyr::mutate(
    tx_start_calc = mapply(
      compute_orf_tx_start,
      tid        = tx_id,
      tx_range   = uORF_Tx_range,
      gstart     = ORF_gstart,
      gstop      = ORF_gstop,
      strand_char = Strand,
      MoreArgs   = list(tx_seqs = tx_seqs_k, exons_by_tx = exons_by_tx_k)
    )
  )

k_res_u <- mapply(
  FUN = function(tid, txs) calc_kozak(tx_seqs_k, tid, txs, natlogratio),
  df_uorf_base$tx_id,
  df_uorf_base$tx_start_calc
)

df_uorf_kozak <- df_uorf_base %>%
  dplyr::mutate(
    Kozak_Seq   = k_res_u[1, ],
    Kozak_Score = as.numeric(k_res_u[2, ])
  ) %>%
  dplyr::filter(!is.na(Kozak_Score), is.finite(Kozak_Score), !is.na(Kozak_Seq), nchar(Kozak_Seq) == 14)

# Reference line for uORF plots: Strong uORF median
ref_uorf_strong <- suppressWarnings(stats::median(df_uorf_kozak$Kozak_Score[df_uorf_kozak$uIF_Group3 == "Strong"], na.rm = TRUE))
if (!is.finite(ref_uorf_strong)) ref_uorf_strong <- NA_real_

# ----------------------------
# 6) Build mORF groups + compute mORF Kozak
# ----------------------------
message("Kozak module: computing mORF Kozak for 5 groups...")

all_uorf_genes <- unique(toupper(substr(as.character(mORF_uORF_Info$tx_id), 1, 9)))

gene_uorf_n <- mORF_uORF_Info %>%
  dplyr::mutate(gene_id = toupper(substr(as.character(tx_id), 1, 9))) %>%
  dplyr::count(gene_id, name = "n_uorfs")
multi_uorf_genes <- gene_uorf_n$gene_id[gene_uorf_n$n_uorfs >= 2]

weak_genes   <- unique(toupper(df_uorf_base$gene_id[df_uorf_base$uIF_Group3 == "Weak"]))
mid_genes    <- unique(toupper(df_uorf_base$gene_id[df_uorf_base$uIF_Group3 == "Mid"]))
strong_genes <- unique(toupper(df_uorf_base$gene_id[df_uorf_base$uIF_Group3 == "Strong"]))

# no-uORF gene universe (prefer TE_CDS; optional RNA filter applied here too)
if (exists("TE_CDS") && "gene_id" %in% names(TE_CDS)) {
  if (USE_RNA_FILTER && "RNA" %in% names(TE_CDS)) {
    gene_univ <- unique(toupper(as.character(TE_CDS$gene_id[is.finite(TE_CDS$RNA) & TE_CDS$RNA >= MIN_RNA_TPM])))
  } else {
    gene_univ <- unique(toupper(as.character(TE_CDS$gene_id)))
  }
} else {
  # fallback: all genes implied by TxDb 5'UTR transcript names
  gene_univ <- unique(toupper(gsub("\\.\\d+$", "", names(fiveUTRs_k))))
}
no_uorf_genes <- setdiff(gene_univ, all_uorf_genes)

# Representative tx: gene+".1"
make_tx1 <- function(g) paste0(g, ".1")

morf_groups <- list(
  "mORF (Weak-uIF genes)"    = make_tx1(weak_genes),
  "mORF (Mid-uIF genes)"     = make_tx1(mid_genes),
  "mORF (Strong-uIF genes)"  = make_tx1(strong_genes),
  "mORF (multi-uORF genes)"  = make_tx1(multi_uorf_genes),
  "mORF (no-uORF genes)"     = make_tx1(no_uorf_genes)
)
morf_order <- names(morf_groups)

valid_tx_morf <- function(tx_vec) {
  tx_vec <- unique(tx_vec)
  intersect(intersect(tx_vec, names(fiveUTRs_k)), names(tx_seqs_k))
}

morf_df <- dplyr::bind_rows(lapply(morf_order, function(lbl) {
  txs <- valid_tx_morf(morf_groups[[lbl]])
  if (length(txs) == 0) return(NULL)
  data.frame(tx_id = txs, MorfGroup = lbl, stringsAsFactors = FALSE)
}))
morf_df$MorfGroup <- factor(morf_df$MorfGroup, levels = morf_order)

# CDS start in transcript coords = spliced 5'UTR length + 1
morf_df <- morf_df %>%
  dplyr::mutate(
    tx_start_calc = vapply(tx_id, function(tid) {
      if (!tid %in% names(fiveUTRs_k)) return(NA_integer_)
      as.integer(sum(IRanges::width(fiveUTRs_k[[tid]]))) + 1L
    }, integer(1))
  )

k_res_m <- mapply(
  FUN = function(tid, txs) calc_kozak(tx_seqs_k, tid, txs, natlogratio),
  morf_df$tx_id,
  morf_df$tx_start_calc
)

df_morf_kozak <- morf_df %>%
  dplyr::mutate(
    Kozak_Seq   = k_res_m[1, ],
    Kozak_Score = as.numeric(k_res_m[2, ])
  ) %>%
  dplyr::filter(!is.na(Kozak_Score), is.finite(Kozak_Score), !is.na(Kozak_Seq), nchar(Kozak_Seq) == 14)

# Reference line for mORF plots: no-uORF median
ref_morf_nouorf <- suppressWarnings(stats::median(df_morf_kozak$Kozak_Score[df_morf_kozak$MorfGroup == "mORF (no-uORF genes)"], na.rm = TRUE))
if (!is.finite(ref_morf_nouorf)) ref_morf_nouorf <- NA_real_

# ----------------------------
# 7) Violin/box + significance letters (+ horizontal ref line)
# ----------------------------
message("Kozak module: plotting violin/box with significance letters...")

plot_violin_with_letters <- function(df, group_col, score_col, title, ylab, hline_val = NA_real_, hline_lab = NULL) {
  df <- df %>% dplyr::filter(!is.na(.data[[score_col]]), is.finite(.data[[score_col]]))
  
  # KS letters
  stat_df <- generate_legend_labels(df, group_col, score_col)
  label_map <- setNames(stat_df$Final_Label, as.character(stat_df[[group_col]]))
  df$Final_Label <- factor(label_map[as.character(df[[group_col]])], levels = stat_df$Final_Label)
  
  y_top <- max(df[[score_col]], na.rm = TRUE)
  y_bot <- min(df[[score_col]], na.rm = TRUE)
  y_pos <- y_top + 0.06 * (y_top - y_bot + 1e-9)
  
  letter_df <- stat_df %>%
    dplyr::mutate(y = y_pos) %>%
    dplyr::mutate(Final_Label = factor(Final_Label, levels = stat_df$Final_Label))
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Final_Label, y = .data[[score_col]])) +
    ggplot2::geom_violin(aes(fill = Final_Label), trim = FALSE, alpha = 0.7) +
    ggplot2::geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA) +
    ggplot2::geom_text(
      data = letter_df,
      aes(x = Final_Label, y = y, label = letters),
      inherit.aes = FALSE,
      size = 5,
      fontface = "bold"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::labs(title = title, x = "", y = ylab, fill = "") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold")
    )
  
  if (is.finite(hline_val)) {
    p <- p +
      ggplot2::geom_hline(yintercept = hline_val, linetype = "dashed", linewidth = 0.8) +
      ggplot2::annotate("text", x = Inf, y = hline_val, label = hline_lab %||% sprintf("ref=%.2f", hline_val),
                        hjust = 1.05, vjust = -0.4, size = 3.5)
  }
  p
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

p_vio_uorf <- plot_violin_with_letters(
  df = df_uorf_kozak %>% dplyr::mutate(Group = GroupLabel),
  group_col = "Group",
  score_col = "Kozak_Score",
  title = sprintf("uORF Kozak scores by uIF group (pUP>=%.0f%s)", MIN_PUP,
                  if (USE_RNA_FILTER) sprintf(", RNA>=%.2f", MIN_RNA_TPM) else ""),
  ylab = "Kozak score",
  hline_val = ref_uorf_strong,
  hline_lab = if (is.finite(ref_uorf_strong)) sprintf("Strong median=%.2f", ref_uorf_strong) else NULL
)

p_vio_morf <- plot_violin_with_letters(
  df = df_morf_kozak %>% dplyr::mutate(Group = MorfGroup),
  group_col = "Group",
  score_col = "Kozak_Score",
  title = "mORF Kozak scores by gene group",
  ylab = "Kozak score",
  hline_val = ref_morf_nouorf,
  hline_lab = if (is.finite(ref_morf_nouorf)) sprintf("No-uORF median=%.2f", ref_morf_nouorf) else NULL
)

p_violin_both <- ggpubr::ggarrange(p_vio_uorf, p_vio_morf, ncol = 1, heights = c(1, 1))
print(p_violin_both)
save_pdf(p_violin_both, "Kozak_Score_uORF_and_mORF_Violin.pdf", w = 11, h = 10)

# ----------------------------
# 8) ECDF plots + significance letters in legend (+ vertical ref line)
# ----------------------------
message("Kozak module: plotting ECDF with significance letters...")

plot_ecdf_with_letters <- function(df, group_col, score_col, title, xlab, vline_val = NA_real_, vline_lab = NULL) {
  df <- df %>% dplyr::filter(!is.na(.data[[score_col]]), is.finite(.data[[score_col]]))
  
  stat_df <- generate_legend_labels(df, group_col, score_col)
  stat_df_join <- stat_df %>% dplyr::select(dplyr::all_of(group_col), Final_Label)
  df <- dplyr::left_join(df, stat_df_join, by = group_col)
  df$Final_Label <- factor(df$Final_Label, levels = stat_df$Final_Label)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[score_col]], color = Final_Label)) +
    ggplot2::stat_ecdf(geom = "step", linewidth = 1.1) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::labs(title = title, x = xlab, y = "Fraction", color = "") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )
  
  if (is.finite(vline_val)) {
    p <- p +
      ggplot2::geom_vline(xintercept = vline_val, linetype = "dashed", linewidth = 0.8) +
      ggplot2::annotate("text", x = vline_val, y = 0.05, label = vline_lab %||% sprintf("ref=%.2f", vline_val),
                        angle = 90, vjust = -0.4, hjust = 0, size = 3.2)
  }
  p
}

p_ecdf_uorf <- plot_ecdf_with_letters(
  df = df_uorf_kozak %>% dplyr::mutate(Group = GroupLabel),
  group_col = "Group",
  score_col = "Kozak_Score",
  title = "ECDF: uORF Kozak scores (uIF groups)",
  xlab = "Kozak score",
  vline_val = ref_uorf_strong,
  vline_lab = if (is.finite(ref_uorf_strong)) sprintf("Strong median=%.2f", ref_uorf_strong) else NULL
)
print(p_ecdf_uorf)
save_pdf(p_ecdf_uorf, "Kozak_ECDF_uORF_uIFgroups.pdf", w = 10, h = 6)

p_ecdf_morf <- plot_ecdf_with_letters(
  df = df_morf_kozak %>% dplyr::mutate(Group = MorfGroup),
  group_col = "Group",
  score_col = "Kozak_Score",
  title = "ECDF: mORF Kozak scores (5 groups)",
  xlab = "Kozak score",
  vline_val = ref_morf_nouorf,
  vline_lab = if (is.finite(ref_morf_nouorf)) sprintf("No-uORF median=%.2f", ref_morf_nouorf) else NULL
)
print(p_ecdf_morf)
save_pdf(p_ecdf_morf, "Kozak_ECDF_mORF_groups.pdf", w = 10, h = 6)

# ----------------------------
# 9) Seqlogos (probability), y 0–0.5, SAME HEIGHT per panel
#    Layout: 5 rows x 2 cols
#      Left col: uORF Weak/Mid/Strong then 2 blanks
#      Right col: 5 mORF groups
# ----------------------------
message("Kozak module: plotting seqlogos (bits; y=0..0.5; equal panel heights)...")

# downsample helper
cap_seqs <- function(x, max_n = MAX_LOGO_SEQS) {
  x <- x[!is.na(x) & nchar(x) == 14]
  if (length(x) > max_n) x <- sample(x, max_n)
  x
}

# Build per-group sequence vectors
uorf_levels <- levels(df_uorf_base$GroupLabel)
df_uorf_logo <- df_uorf_kozak %>%
  dplyr::mutate(Group = GroupLabel) %>%
  dplyr::filter(!is.na(Kozak_Seq), nchar(Kozak_Seq) == 14)
df_uorf_logo$Group <- factor(df_uorf_logo$Group, levels = uorf_levels)

morf_levels <- morf_order
df_morf_logo <- df_morf_kozak %>%
  dplyr::mutate(Group = MorfGroup) %>%
  dplyr::filter(!is.na(Kozak_Seq), nchar(Kozak_Seq) == 14)
df_morf_logo$Group <- factor(df_morf_logo$Group, levels = morf_levels)

get_group_seqs <- function(df, grp, seq_col = "Kozak_Seq") {
  cap_seqs(df %>% dplyr::filter(.data[["Group"]] == grp) %>% dplyr::pull(.data[[seq_col]]))
}

make_logo_plot <- function(seqs, title_txt) {
  if (length(seqs) < 20) {
    return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::ggtitle(paste0(title_txt, "\n(n<20)")))
  }
  ggseqlogo::ggseqlogo(seqs, method = "bits") +
    ggplot2::coord_cartesian(ylim = c(0, SEQLOGO_YMAX)) +
    ggplot2::labs(title = title_txt) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 10),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    )
}

# uORF panels (3 + 2 blanks)
u_titles <- as.character(uorf_levels)
u_plots <- lapply(u_titles, function(tt) make_logo_plot(get_group_seqs(df_uorf_logo, tt), tt))
while (length(u_plots) < 5) u_plots[[length(u_plots) + 1]] <- ggplot2::ggplot() + ggplot2::theme_void()

# mORF panels (5)
m_titles <- as.character(morf_levels)
m_plots <- lapply(m_titles, function(tt) make_logo_plot(get_group_seqs(df_morf_logo, tt), tt))

# arrange row-wise: (u1,m1,u2,m2,...)
plotlist <- list(
  u_plots[[1]], m_plots[[1]],
  u_plots[[2]], m_plots[[2]],
  u_plots[[3]], m_plots[[3]],
  u_plots[[4]], m_plots[[4]],
  u_plots[[5]], m_plots[[5]]
)

p_logo_grid <- ggpubr::ggarrange(
  plotlist = plotlist,
  ncol = 2, nrow = 5,
  labels = NULL,
  align = "v"
)

p_logo_grid <- ggpubr::annotate_figure(
  p_logo_grid,
  top = ggpubr::text_grob("Kozak seqlogos (probability; y=0–0.5)\nLeft: uORFs (uIF groups) | Right: mORFs (5 groups)",
                          face = "bold", size = 12)
)

print(p_logo_grid)
save_pdf(p_logo_grid, "Kozak_SeqLogo_uORF_and_mORF.pdf", w = 12, h = 12)

message("Kozak module complete. Outputs written to OUT_DIR = ", OUT_DIR)






################################################################################
# MODULE: CORRELATION PLOTS (Before vs After pUP>=30 & RNA>=0.1)
#  (1) dist(uORF_STOP -> CDS_start) vs TE
#  (2) dist(50% CAGE -> CDS_start) vs TE           [tx coords; no introns]
#  (3) uORF Kozak score vs TE
#  (4) dist(50% CAGE -> uORF_START) vs TE          [NEW]
#  (5) NO-uORF genes: dist(50% CAGE -> CDS_start) vs TE
#
# Pearson + Spearman + p-values via cor.test() + annotate()
#      (avoid ggpubr::stat_cor parse/grob errors)
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(ggpubr)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
  library(Rsamtools)
  library(Biostrings)
  library(S4Vectors)
})

# ----------------------------
# 0) Parameters
# ----------------------------
MIN_PUP_COR <- 30
MIN_RNA_COR <- 0.1

# ----------------------------
# 1) Sanity checks
# ----------------------------
if (!exists("OUT_df5_clean")) stop("OUT_df5_clean missing (run the pipeline up to Step 4).")
if (!exists("mORF_uORF_Info")) stop("mORF_uORF_Info missing (run the RT+CiPS merge step).")
if (!exists("RNA") || !exists("Ribo")) stop("RNA/Ribo missing (run Step 5).")
if (!exists("fiveUTRByGeneNames")) {
  if (file.exists(PATH_5UTR_GENES)) {
    fiveUTRByGeneNames <- read.delim(PATH_5UTR_GENES, header = TRUE)$gene_id
  } else {
    stop("fiveUTRByGeneNames missing and PATH_5UTR_GENES not found.")
  }
}
if (!exists("save_pdf")) stop("save_pdf() missing.")
if (!exists("OUT_DIR")) stop("OUT_DIR missing.")
if (!"uORF_Tx_range" %in% colnames(OUT_df5_clean)) stop("OUT_df5_clean must contain uORF_Tx_range.")
if (!exists("PATH_GTF")) stop("PATH_GTF missing.")
if (!exists("PATH_CAGE")) stop("PATH_CAGE missing.")

# ----------------------------
# 2) Build TE table WITHOUT RNA>1 restriction (so RNA>=0.1 filter matters)
# ----------------------------
TE_CDS_cor <- data.frame(
  gene_id = as.character(RNA$gene_id),
  RNA = as.numeric(RNA$TPM),
  Ribo = as.numeric(Ribo$TPM),
  TE = as.numeric(Ribo$TPM / RNA$TPM),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(gene_id = as.character(gene_id)) %>%
  dplyr::filter(gene_id %in% fiveUTRByGeneNames) %>%
  dplyr::filter(!is.na(RNA) & is.finite(RNA) & RNA > 0) %>%
  dplyr::filter(!is.na(TE) & is.finite(TE))

# ----------------------------
# 3) Define single-uORF genes (match the OUT_df6 logic)
# ----------------------------
gene_counts_cor <- table(as.character(OUT_df5_clean$gene_id))
single_uorf_genes_cor <- names(gene_counts_cor)[gene_counts_cor == 1]

# ----------------------------
# 4) Helper: parse BOTH endpoints from uORF_Tx_range like "123-456"
# ----------------------------
parse_tx_start_stop_from_range <- function(x) {
  if (is.na(x) || !nzchar(x)) return(c(NA_integer_, NA_integer_))
  nums <- regmatches(as.character(x), gregexpr("\\d+", as.character(x)))[[1]]
  if (length(nums) < 2) return(c(NA_integer_, NA_integer_))
  a <- as.integer(nums[1]); b <- as.integer(nums[2])
  c(min(a, b), max(a, b))
}

safe_range <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) == 0) return(NULL)
  range(v)
}

# ----------------------------
# 5) Build TxDb once for CDS start in transcript coords
# ----------------------------
txdb_cor <- if (exists("txdb_k")) txdb_k else if (exists("txdb")) txdb else {
  GenomicFeatures::makeTxDbFromGFF(
    file = PATH_GTF, format = "gtf",
    dataSource = "Araport11",
    organism = "Arabidopsis"
  )
}

fiveUTR_by_tx_cor <- GenomicFeatures::fiveUTRsByTranscript(txdb_cor, use.names = TRUE)
utr_len_cor <- vapply(fiveUTR_by_tx_cor, function(gr) sum(IRanges::width(gr)), integer(1))
cds_start_tx_cor <- utr_len_cor + 1L  # CDS start in transcript coords

# ----------------------------
# 6) Build uORF correlation base dataset (SINGLE-uORF genes, with TE/RNA/pUP)
#     (1) uses uORF STOP -> CDS start
# ----------------------------
tmp_ss <- t(vapply(as.character(OUT_df5_clean$uORF_Tx_range),
                   parse_tx_start_stop_from_range,
                   integer(2)))
colnames(tmp_ss) <- c("uORF_start_tx", "uORF_stop_tx")

df_uorf_cor_base <- OUT_df5_clean %>%
  dplyr::mutate(
    gene_id = as.character(gene_id),
    tx_id   = as.character(tx_id),
    pUP     = as.numeric(pUP),
    uORF_Tx_range = as.character(uORF_Tx_range),
    uORF_start_tx = as.integer(tmp_ss[, "uORF_start_tx"]),
    uORF_stop_tx  = as.integer(tmp_ss[, "uORF_stop_tx"]),
    CDS_start_tx  = as.integer(cds_start_tx_cor[tx_id]),
    dist_uORFstop_to_CDS_nt = as.numeric(CDS_start_tx - uORF_stop_tx)
  ) %>%
  dplyr::filter(gene_id %in% single_uorf_genes_cor) %>%
  dplyr::left_join(TE_CDS_cor, by = "gene_id") %>%
  dplyr::filter(!is.na(TE) & is.finite(TE)) %>%
  dplyr::filter(!is.na(dist_uORFstop_to_CDS_nt) & is.finite(dist_uORFstop_to_CDS_nt)) %>%
  # keep only uORFs upstream of CDS (avoid negative distances)
  dplyr::filter(!is.na(uORF_stop_tx) & !is.na(CDS_start_tx) & uORF_stop_tx <= CDS_start_tx)

df_uorf_cor_before <- df_uorf_cor_base
df_uorf_cor_after  <- df_uorf_cor_base %>%
  dplyr::filter(!is.na(pUP) & is.finite(pUP) & pUP >= MIN_PUP_COR) %>%
  dplyr::filter(!is.na(RNA) & is.finite(RNA) & RNA >= MIN_RNA_COR)

# ----------------------------
# 7) Compute CAGE 50% position in TRANSCRIPT COORDS (no introns)
#     returns per-tx: cage50_pos_tx and dist to CDS start
# ----------------------------
compute_cage50_to_cds_txcoords <- function(tx_ids, txdb_obj, cds_start_tx_vec, path_cage) {
  tx_ids <- unique(as.character(tx_ids))
  tx_ids <- tx_ids[!is.na(tx_ids)]
  tx_ids <- tx_ids[tx_ids %in% names(cds_start_tx_vec)]
  if (length(tx_ids) == 0) {
    return(data.frame(tx_id = character(0), cage50_pos_tx = integer(0), dist_cage50_to_CDS_nt = numeric(0)))
  }
  
  CAGE_rawX <- read.delim(path_cage, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  colnames(CAGE_rawX) <- c("count", "chr", "start", "strand")
  CAGE_rawX$end <- CAGE_rawX$start
  CAGE_rawX$chr <- gsub("^chr|^Chr", "", CAGE_rawX$chr, ignore.case = TRUE)
  CAGE_grX <- GenomicRanges::makeGRangesFromDataFrame(CAGE_rawX, keep.extra.columns = TRUE)
  
  ex_by_tx <- GenomicFeatures::exonsBy(txdb_obj, by = "tx", use.names = TRUE)
  ex_by_tx <- ex_by_tx[tx_ids]
  
  GenomeInfoDb::seqlevelsStyle(CAGE_grX) <- "NCBI"
  tryCatch({ GenomeInfoDb::seqlevelsStyle(ex_by_tx) <- "NCBI" }, error = function(e) {})
  
  ex_flat <- BiocGenerics::unlist(ex_by_tx, use.names = FALSE)
  common_levels <- intersect(GenomeInfoDb::seqlevels(CAGE_grX), GenomeInfoDb::seqlevels(ex_flat))
  if (length(common_levels) > 0) {
    CAGE_grX <- GenomeInfoDb::keepSeqlevels(CAGE_grX, common_levels, pruning.mode = "coarse")
    ex_by_tx <- GenomeInfoDb::keepSeqlevels(ex_by_tx, common_levels, pruning.mode = "coarse")
  }
  
  mapped <- tryCatch(
    GenomicFeatures::mapToTranscripts(CAGE_grX, ex_by_tx, ignore.strand = FALSE),
    error = function(e) GenomicFeatures::mapToTranscripts(CAGE_grX, ex_by_tx)
  )
  if (length(mapped) == 0) {
    return(data.frame(tx_id = character(0), cage50_pos_tx = integer(0), dist_cage50_to_CDS_nt = numeric(0)))
  }
  
  mc <- S4Vectors::mcols(mapped)
  if (!("xHits" %in% colnames(mc))) stop("mapToTranscripts() output missing xHits; cannot recover CAGE counts.")
  
  dt <- data.table::data.table(
    tx_id = as.character(GenomicRanges::seqnames(mapped)),
    pos_tx = as.integer(GenomicRanges::start(mapped)),
    count = as.numeric(S4Vectors::mcols(CAGE_grX)$count[mc$xHits])
  )
  dt <- dt[!is.na(tx_id) & !is.na(pos_tx) & !is.na(count)]
  dt[, cds_start_tx := as.integer(cds_start_tx_vec[tx_id])]
  dt <- dt[!is.na(cds_start_tx)]
  dt <- dt[pos_tx >= 1 & pos_tx <= cds_start_tx]
  
  if (nrow(dt) == 0) {
    return(data.frame(tx_id = character(0), cage50_pos_tx = integer(0), dist_cage50_to_CDS_nt = numeric(0)))
  }
  
  dt <- dt[, .(count = sum(count, na.rm = TRUE)), by = .(tx_id, pos_tx, cds_start_tx)]
  data.table::setorder(dt, tx_id, pos_tx)
  dt[, cum_sum := cumsum(count), by = tx_id]
  dt[, total_sum := sum(count, na.rm = TRUE), by = tx_id]
  dt <- dt[total_sum > 0]
  
  dt50 <- dt[cum_sum / total_sum >= 0.5, .SD[1], by = tx_id]
  dt50[, dist_cage50_to_CDS_nt := as.numeric(cds_start_tx - pos_tx)]
  
  as.data.frame(dt50[, .(tx_id, cage50_pos_tx = pos_tx, dist_cage50_to_CDS_nt)])
}

# uORF genes: CAGE50 info
df_cage50_uorf <- compute_cage50_to_cds_txcoords(
  tx_ids = df_uorf_cor_before$tx_id,
  txdb_obj = txdb_cor,
  cds_start_tx_vec = cds_start_tx_cor,
  path_cage = PATH_CAGE
)

# (2) 50% CAGE -> CDS start
df_uorf_cage_before <- df_uorf_cor_before %>%
  dplyr::left_join(df_cage50_uorf, by = "tx_id") %>%
  dplyr::filter(!is.na(dist_cage50_to_CDS_nt) & is.finite(dist_cage50_to_CDS_nt))

df_uorf_cage_after <- df_uorf_cor_after %>%
  dplyr::left_join(df_cage50_uorf, by = "tx_id") %>%
  dplyr::filter(!is.na(dist_cage50_to_CDS_nt) & is.finite(dist_cage50_to_CDS_nt))

# (4) NEW: 50% CAGE -> uORF start
df_uorf_cage_to_uorfstart_before <- df_uorf_cor_before %>%
  dplyr::left_join(df_cage50_uorf, by = "tx_id") %>%
  dplyr::mutate(dist_cage50_to_uORFstart_nt = as.numeric(uORF_start_tx - cage50_pos_tx)) %>%
  dplyr::filter(!is.na(dist_cage50_to_uORFstart_nt) & is.finite(dist_cage50_to_uORFstart_nt)) %>%
  dplyr::filter(dist_cage50_to_uORFstart_nt >= 0)

df_uorf_cage_to_uorfstart_after <- df_uorf_cor_after %>%
  dplyr::left_join(df_cage50_uorf, by = "tx_id") %>%
  dplyr::mutate(dist_cage50_to_uORFstart_nt = as.numeric(uORF_start_tx - cage50_pos_tx)) %>%
  dplyr::filter(!is.na(dist_cage50_to_uORFstart_nt) & is.finite(dist_cage50_to_uORFstart_nt)) %>%
  dplyr::filter(dist_cage50_to_uORFstart_nt >= 0)

# ----------------------------
# 8) Kozak score for uORF start vs TE (build natlogratio if missing)
# ----------------------------
build_natlogratio_if_missing <- function() {
  if (exists("natlogratio") && is.matrix(natlogratio) && ncol(natlogratio) == 14) return(invisible(TRUE))
  
  if (!exists("PATH_FASTA_TAIR10")) stop("PATH_FASTA_TAIR10 missing.")
  if (!file.exists(PATH_FASTA_TAIR10)) stop("FASTA not found: ", PATH_FASTA_TAIR10)
  
  tryCatch({ Rsamtools::indexFa(PATH_FASTA_TAIR10) }, error = function(e) NULL)
  fa <- Rsamtools::FaFile(PATH_FASTA_TAIR10)
  
  tx_seqs_tmp <- GenomicFeatures::extractTranscriptSeqs(fa, txdb_cor, use.names = TRUE)
  fiveUTR_tmp <- GenomicFeatures::fiveUTRsByTranscript(txdb_cor, use.names = TRUE)
  common_tx <- intersect(names(tx_seqs_tmp), names(fiveUTR_tmp))
  
  utr_len <- vapply(fiveUTR_tmp[common_tx], function(gr) sum(IRanges::width(gr)), integer(1))
  nt_seqs <- as.character(tx_seqs_tmp[common_tx])
  
  ref_df <- data.frame(
    tx_id = common_tx,
    fUTR_length = as.integer(utr_len),
    nt_seq = nt_seqs,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      nt_len = nchar(nt_seq),
      ok = !is.na(fUTR_length) & fUTR_length >= 9 & (fUTR_length + 5) <= nt_len,
      Kozak = ifelse(ok, substr(nt_seq, fUTR_length - 8, fUTR_length + 5), NA_character_)
    ) %>%
    dplyr::filter(!is.na(Kozak), nchar(Kozak) == 14, substr(Kozak, 10, 12) == "ATG")
  
  bases <- c("A","C","G","T")
  full_counts <- Biostrings::consensusMatrix(Biostrings::DNAStringSet(ref_df$Kozak), as.prob = FALSE)[bases, , drop = FALSE] + 0.01
  full_consensus <- prop.table(full_counts, 2)
  letter_freq <- full_consensus[bases, , drop = FALSE]
  aF <- Biostrings::alphabetFrequency(Biostrings::DNAStringSet(ref_df$Kozak), as.prob = TRUE, collapse = TRUE)[bases]
  nat <- log(letter_freq / aF)
  if (ncol(nat) == 14) nat[, 10:12] <- 0
  
  natlogratio <<- nat
  tx_seqs_k <<- tx_seqs_tmp
  invisible(TRUE)
}

build_natlogratio_if_missing()

if (!exists("tx_seqs_k")) {
  tryCatch({ Rsamtools::indexFa(PATH_FASTA_TAIR10) }, error = function(e) NULL)
  fa <- Rsamtools::FaFile(PATH_FASTA_TAIR10)
  tx_seqs_k <- GenomicFeatures::extractTranscriptSeqs(fa, txdb_cor, use.names = TRUE)
}

kozak_score_one <- function(tid, start_pos, seqs, mat14) {
  if (is.na(tid) || is.na(start_pos)) return(NA_real_)
  if (!tid %in% names(seqs)) return(NA_real_)
  w_start <- as.integer(start_pos) - 9L
  w_end   <- as.integer(start_pos) + 4L
  if (w_start < 1L) return(NA_real_)
  if (w_end > length(seqs[[tid]])) return(NA_real_)
  kseq <- as.character(Biostrings::subseq(seqs[[tid]], start = w_start, end = w_end))
  if (nchar(kseq) != 14) return(NA_real_)
  ch <- strsplit(kseq, "")[[1]]
  if (length(ch) != 14) return(NA_real_)
  if (!all(ch %in% rownames(mat14))) return(NA_real_)
  sum(diag(mat14[ch, 1:14, drop = FALSE]))
}

df_uorf_kozak_before <- df_uorf_cor_before %>%
  dplyr::mutate(
    Kozak_Score_uORF = mapply(kozak_score_one, tx_id, uORF_start_tx,
                              MoreArgs = list(seqs = tx_seqs_k, mat14 = natlogratio))
  ) %>%
  dplyr::filter(!is.na(Kozak_Score_uORF) & is.finite(Kozak_Score_uORF))

df_uorf_kozak_after <- df_uorf_kozak_before %>%
  dplyr::filter(!is.na(pUP) & is.finite(pUP) & pUP >= MIN_PUP_COR) %>%
  dplyr::filter(!is.na(RNA) & is.finite(RNA) & RNA >= MIN_RNA_COR)

# ----------------------------
# 9) NO-uORF genes set + CAGE50 distance  (plot 5)
# ----------------------------
genes_with_uorf_any <- unique(toupper(substr(as.character(mORF_uORF_Info$tx_id), 1, 9)))
genes_univ_te <- unique(toupper(as.character(TE_CDS_cor$gene_id)))
genes_no_uorf <- setdiff(genes_univ_te, genes_with_uorf_any)

tx_map <- GenomicFeatures::transcripts(txdb_cor, columns = c("TXNAME","GENEID"))
df_tx_map <- data.frame(
  tx_id = as.character(tx_map$TXNAME),
  gene_id = as.character(unlist(tx_map$GENEID)),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    gene_id_u = toupper(substr(gene_id, 1, 9)),
    tx_id = as.character(tx_id)
  ) %>%
  dplyr::filter(gene_id_u %in% genes_no_uorf) %>%
  dplyr::filter(tx_id %in% names(cds_start_tx_cor)) %>%
  dplyr::arrange(gene_id_u, tx_id) %>%
  dplyr::group_by(gene_id_u) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

df_cage50_nouorf <- compute_cage50_to_cds_txcoords(
  tx_ids = df_tx_map$tx_id,
  txdb_obj = txdb_cor,
  cds_start_tx_vec = cds_start_tx_cor,
  path_cage = PATH_CAGE
)

df_no_uorf_plot <- df_tx_map %>%
  dplyr::left_join(df_cage50_nouorf, by = "tx_id") %>%
  dplyr::left_join(TE_CDS_cor %>% dplyr::mutate(gene_id_u = toupper(gene_id)),
                   by = "gene_id_u") %>%
  dplyr::filter(!is.na(dist_cage50_to_CDS_nt) & is.finite(dist_cage50_to_CDS_nt)) %>%
  dplyr::filter(!is.na(TE) & is.finite(TE))

# ----------------------------
# 10) Plot helper (Pearson + Spearman + p-values)
# ----------------------------
fmt_p <- function(p) {
  if (!is.finite(p)) return("NA")
  format.pval(p, digits = 2, eps = 1e-3)
}
fmt_num <- function(x) {
  if (!is.finite(x)) return("NA")
  sprintf("%.2f", x)
}

make_cor_plot <- function(df, xcol, ycol, title_txt, xlab_txt, ylab_txt, xlim_use = NULL, ylim_use = NULL) {
  df2 <- df %>% dplyr::filter(is.finite(.data[[xcol]]), is.finite(.data[[ycol]]))
  n <- nrow(df2)
  
  # compute correlations (robust)
  lab <- paste0("n=", n)
  if (n >= 3) {
    x <- df2[[xcol]]
    y <- df2[[ycol]]
    
    pear <- tryCatch(stats::cor.test(x, y, method = "pearson"), error = function(e) NULL)
    spea <- tryCatch(stats::cor.test(x, y, method = "spearman", exact = FALSE), error = function(e) NULL)
    
    r_val   <- if (!is.null(pear) && is.finite(unname(pear$estimate))) unname(pear$estimate) else NA_real_
    p_val   <- if (!is.null(pear)) pear$p.value else NA_real_
    rho_val <- if (!is.null(spea) && is.finite(unname(spea$estimate))) unname(spea$estimate) else NA_real_
    p_rho   <- if (!is.null(spea)) spea$p.value else NA_real_
    
    lab <- paste0(
      "Pearson r = ", fmt_num(r_val), " (p=", fmt_p(p_val), ")\n",
      "Spearman ρ = ", fmt_num(rho_val), " (p=", fmt_p(p_rho), ")\n",
      "n=", n
    )
  }
  
  # text position inside plot window
  x_min <- if (!is.null(xlim_use)) xlim_use[1] else min(df2[[xcol]], na.rm = TRUE)
  x_max <- if (!is.null(xlim_use)) xlim_use[2] else max(df2[[xcol]], na.rm = TRUE)
  y_min <- if (!is.null(ylim_use)) ylim_use[1] else min(df2[[ycol]], na.rm = TRUE)
  y_max <- if (!is.null(ylim_use)) ylim_use[2] else max(df2[[ycol]], na.rm = TRUE)
  if (!is.finite(x_min) || !is.finite(x_max)) { x_min <- 0; x_max <- 1 }
  if (!is.finite(y_min) || !is.finite(y_max)) { y_min <- 0; y_max <- 1 }
  if (x_max == x_min) x_max <- x_min + 1
  if (y_max == y_min) y_max <- y_min + 1
  
  x_text <- x_min + 0.02 * (x_max - x_min)
  y_text <- y_max - 0.02 * (y_max - y_min)
  
  p <- ggplot2::ggplot(df2, ggplot2::aes(x = .data[[xcol]], y = .data[[ycol]])) +
    ggplot2::geom_point(alpha = 0.35, size = 1) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.8) +
    ggplot2::annotate("text", x = x_text, y = y_text, label = lab, hjust = 0, vjust = 1, size = 3.2) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::labs(title = title_txt, x = xlab_txt, y = ylab_txt)
  
  if (!is.null(xlim_use) || !is.null(ylim_use)) {
    p <- p + ggplot2::coord_cartesian(xlim = xlim_use, ylim = ylim_use)
  } else {
    p <- p + ggplot2::coord_cartesian()
  }
  
  p
}

TE_XLIM <- c(0, 4)

yl1 <- safe_range(df_uorf_cor_before$dist_uORFstop_to_CDS_nt)
yl2 <- safe_range(df_uorf_cage_before$dist_cage50_to_CDS_nt)
yl3 <- safe_range(df_uorf_kozak_before$Kozak_Score_uORF)
yl4 <- safe_range(df_uorf_cage_to_uorfstart_before$dist_cage50_to_uORFstart_nt)

# ----------------------------
# 11) Build plots (Before vs After)
# ----------------------------

# (1) uORF stop -> CDS start
p1_before <- make_cor_plot(
  df_uorf_cor_before, "TE", "dist_uORFstop_to_CDS_nt",
  "1) uORF stop → mORF CDS start (Before)",
  "mORF TE", "Distance (nt)", xlim_use = TE_XLIM, ylim_use = yl1
)
p1_after <- make_cor_plot(
  df_uorf_cor_after, "TE", "dist_uORFstop_to_CDS_nt",
  sprintf("1) uORF stop → mORF CDS start (After: pUP≥%d & RNA≥%.1f)", MIN_PUP_COR, MIN_RNA_COR),
  "mORF TE", "Distance (nt)", xlim_use = TE_XLIM, ylim_use = yl1
)

# (2) 50% CAGE -> CDS start
p2_before <- make_cor_plot(
  df_uorf_cage_before, "TE", "dist_cage50_to_CDS_nt",
  "2) 50% CAGE → mORF CDS start (Before)",
  "mORF TE", "Distance (nt, transcript coords)", xlim_use = TE_XLIM, ylim_use = yl2
)
p2_after <- make_cor_plot(
  df_uorf_cage_after, "TE", "dist_cage50_to_CDS_nt",
  sprintf("2) 50%% CAGE → mORF CDS start (After: pUP≥%d & RNA≥%.1f)", MIN_PUP_COR, MIN_RNA_COR),
  "mORF TE", "Distance (nt, transcript coords)", xlim_use = TE_XLIM, ylim_use = yl2
)

# (3) uORF Kozak score
p3_before <- make_cor_plot(
  df_uorf_kozak_before, "TE", "Kozak_Score_uORF",
  "3) uORF Kozak score (Before)",
  "mORF TE", "uORF Kozak score", xlim_use = TE_XLIM, ylim_use = yl3
)
p3_after <- make_cor_plot(
  df_uorf_kozak_after, "TE", "Kozak_Score_uORF",
  sprintf("3) uORF Kozak score (After: pUP≥%d & RNA≥%.1f)", MIN_PUP_COR, MIN_RNA_COR),
  "mORF TE", "uORF Kozak score", xlim_use = TE_XLIM, ylim_use = yl3
)

# (4) NEW: 50% CAGE -> uORF start
p4_before <- make_cor_plot(
  df_uorf_cage_to_uorfstart_before, "TE", "dist_cage50_to_uORFstart_nt",
  "4) 50% CAGE → uORF start (Before)",
  "mORF TE", "Distance (nt, transcript coords)", xlim_use = TE_XLIM, ylim_use = yl4
)
p4_after <- make_cor_plot(
  df_uorf_cage_to_uorfstart_after, "TE", "dist_cage50_to_uORFstart_nt",
  sprintf("4) 50%% CAGE → uORF start (After: pUP≥%d & RNA≥%.1f)", MIN_PUP_COR, MIN_RNA_COR),
  "mORF TE", "Distance (nt, transcript coords)", xlim_use = TE_XLIM, ylim_use = yl4
)

# (5) NO uORF genes: 50% CAGE -> CDS start
p5_nouorf <- make_cor_plot(
  df_no_uorf_plot %>% dplyr::filter(RNA >= MIN_RNA_COR),
  "TE", "dist_cage50_to_CDS_nt",
  sprintf("5) NO uORF genes: 50%% CAGE → mORF CDS start (RNA≥%.1f)", MIN_RNA_COR),
  "mORF TE", "Distance (nt, transcript coords)", xlim_use = TE_XLIM, ylim_use = NULL
)

# ----------------------------
# 12) Arrange + Save
# ----------------------------
blank <- ggplot2::ggplot() + ggplot2::theme_void()

Fig_Cor_All <- ggpubr::ggarrange(
  p1_before, p1_after,
  p2_before, p2_after,
  p3_before, p3_after,
  p4_before, p4_after,
  p5_nouorf, blank,
  ncol = 2, nrow = 5,
  align = "hv"
)

Fig_Cor_All <- ggpubr::annotate_figure(
  Fig_Cor_All,
  top = ggpubr::text_grob(
    sprintf("Correlation panels (Before vs After pUP≥%d & RNA≥%.1f)\nPearson r and Spearman ρ with p-values are shown on each panel. TE xlim=0–4.",
            MIN_PUP_COR, MIN_RNA_COR),
    face = "bold", size = 12
  )
)

print(Fig_Cor_All)
save_pdf(
  Fig_Cor_All,
  "Correlation_uORF_BeforeAfter_pUP30_RNA0.1_PLUS_CAGEtoUORFstart_PLUS_NoTuORF_WITH_Pearson_Spearman.pdf",
  w = 14, h = 22
)

message("Saved PDF -> ", OUT_DIR)
################################################################################









################################################################################
# RECONCILIATION MODULE
# Internal consistency checks: recomputes every headline count reported in the
# manuscript (TuORF totals, UuORF groups, uORFs per gene, length bins) so the
# published numbers can be traced to the data in one place.
# Writes:  Reconciliation_A_numbers.txt (human-readable summary)
#          Reconcile_A1_RT_CiPS_union.csv
#          Reconcile_A2_UuORF_groups.csv
#          Reconcile_A3_uORFs_per_gene.csv
#          Reconcile_A4_length_bins.csv
# All into OUT_DIR. Uses the TSS_SOURCE selected in section 0 (default CAGE,
# the data source behind the disputed manuscript numbers).
################################################################################
message("RECONCILIATION MODULE: starting...")

.recon_lines <- c(
  "================================================================================",
  paste0("RECONCILIATION NUMBERS  (TSS_SOURCE = ", TSS_SOURCE,
         "; generated ", format(Sys.time(), "%Y-%m-%d %H:%M"), ")"),
  "Resolves Section A of Weak_uORF_issue_report_2026_08_12.docx",
  "================================================================================")
.add <- function(...) .recon_lines <<- c(.recon_lines, paste0(...))

## ---------------------------------------------------------------------------
## A1) RiboTaper vs CiPS: counts, overlap, union
## ---------------------------------------------------------------------------
tryCatch({
  stopifnot(exists("OUT_df5_clean"), exists("DF_RT2"), exists("DF_CiPS2"),
            exists("overlap_keys"), exists("DF_All_Evidence"))

  n_rt_raw    <- nrow(DF_RT2)
  n_cips_raw  <- nrow(DF_CiPS2)
  n_keys_both <- length(overlap_keys)
  n_union_raw <- nrow(DF_All_Evidence)

  meth <- table(factor(OUT_df5_clean$method, levels = c("RT", "Both", "CiPS")))
  n_rt_only   <- as.integer(meth["RT"])
  n_both      <- as.integer(meth["Both"])
  n_cips_only <- as.integer(meth["CiPS"])
  n_union     <- nrow(OUT_df5_clean)
  n_rt_tot    <- n_rt_only + n_both
  n_cips_tot  <- n_cips_only + n_both
  n_genes_tot <- dplyr::n_distinct(OUT_df5_clean$gene_id)

  .add("")
  .add("--- A1. RiboTaper + CiPS TuORF accounting -----------------------------")
  .add("Raw evidence rows            : RT = ", n_rt_raw, " | CiPS = ", n_cips_raw,
       " | overlapping keys = ", n_keys_both, " | union rows = ", n_union_raw)
  .add("FINAL deduplicated TuORF set (OUT_df5_clean, basis of Table S1A):")
  .add("  TuORFs identified by RiboTaper (incl. shared)  : ", n_rt_tot)
  .add("  TuORFs identified by CiPS      (incl. shared)  : ", n_cips_tot)
  .add("  Identified by BOTH (overlap)                   : ", n_both)
  .add("  RiboTaper-only                                 : ", n_rt_only)
  .add("  CiPS-only                                      : ", n_cips_only)
  .add("  TOTAL unique TuORFs (union)                    : ", n_union)
  .add("  Check: RT_total + CiPS_total - overlap = ",
       n_rt_tot + n_cips_tot - n_both, "  (must equal union)")
  .add("  TOTAL TuORF genes                              : ", n_genes_tot)
  .add("  Manuscript claim: 2113 (RT) + 6633 (CiPS) = total 7051, implying overlap 1695.")

  utils::write.csv(
    data.frame(metric = c("RT_total_incl_shared", "CiPS_total_incl_shared",
                          "Overlap_Both", "RT_only", "CiPS_only",
                          "Union_total_TuORFs", "Total_TuORF_genes",
                          "raw_RT_rows", "raw_CiPS_rows", "raw_overlap_keys", "raw_union_rows"),
               value  = c(n_rt_tot, n_cips_tot, n_both, n_rt_only, n_cips_only,
                          n_union, n_genes_tot,
                          n_rt_raw, n_cips_raw, n_keys_both, n_union_raw)),
    file.path(OUT_DIR, "Reconcile_A1_RT_CiPS_union.csv"), row.names = FALSE)
}, error = function(e) .add("A1 FAILED: ", conditionMessage(e)))

## ---------------------------------------------------------------------------
## A2) UuORF gene categories (0~30 / 30~80 / 80~100), incl. no-CAGE-signal genes
## ---------------------------------------------------------------------------
tryCatch({
  stopifnot(exists("UuORF_gene_df_4"), exists("hist_stats_4"))

  n_total_uu <- nrow(UuORF_gene_df_4)
  n_zero     <- sum(UuORF_gene_df_4$C_total == 0, na.rm = TRUE)
  zero_grp   <- table(UuORF_gene_df_4$UU_Group[UuORF_gene_df_4$C_total == 0])

  .add("")
  .add("--- A2. UuORF gene categories (Figure S5D) ----------------------------")
  for (i in seq_len(nrow(hist_stats_4))) {
    .add(sprintf("  %-7s : n = %5d  (%.1f%%)", hist_stats_4$UU_Group[i],
                 hist_stats_4$n[i], hist_stats_4$pct[i]))
  }
  .add("  TOTAL UuORF genes                        : ", n_total_uu)
  .add("  Genes with NO CAGE/TSS signal (C_total=0): ", n_zero,
       "  (", sprintf("%.1f", 100 * n_zero / n_total_uu), "% of total)")
  .add("  -> zero-signal genes fall in group(s): ",
       paste(names(zero_grp), zero_grp, sep = "=", collapse = ", "),
       "  (pUP is set to 0 when C_total = 0, so they are counted in 0~30)")
  .add("  Percent excluding zero-signal genes from 0~30: ",
       sprintf("%.1f%%", 100 * (hist_stats_4$n[1] - n_zero) / n_total_uu))
  .add("  Manuscript claim: 2,400 + 126 + 147 = 2,673 of 2,713 (40 genes unaccounted; 2,400/2,713 = 88.5%, not 90%).")

  a2 <- hist_stats_4
  a2$total <- n_total_uu
  a2$n_zero_CAGE_signal_total <- n_zero
  utils::write.csv(a2, file.path(OUT_DIR, "Reconcile_A2_UuORF_groups.csv"), row.names = FALSE)
}, error = function(e) .add("A2 FAILED: ", conditionMessage(e)))

## ---------------------------------------------------------------------------
## A3) Genes by number of TuORFs (Figure S3 / Table S1A accounting)
## ---------------------------------------------------------------------------
tryCatch({
  stopifnot(exists("OUT_df5_clean"))

  cnt_exh <- table(table(OUT_df5_clean$gene_id))          # exhaustive, deduped set
  ks <- as.integer(names(cnt_exh)); ns <- as.integer(cnt_exh)
  exh <- data.frame(n_uORFs_per_gene = ks, n_genes = ns,
                    n_TuORFs = ks * ns, set = "exhaustive_OUT_df5_clean")

  .add("")
  .add("--- A3. Genes by number of TuORFs -------------------------------------")
  .add("EXHAUSTIVE (deduplicated TuORF set; sums to the Table S1A totals):")
  for (i in seq_along(ks))
    .add(sprintf("  genes with %d TuORF(s): %5d   (=%5d TuORFs)", ks[i], ns[i], ks[i]*ns[i]))
  .add("  TOTAL genes  : ", sum(ns))
  .add("  TOTAL TuORFs : ", sum(ks * ns))

  filt <- NULL
  if (exists("df_ranked_final")) {
    per_gene <- df_ranked_final[!duplicated(df_ranked_final$gene_id),
                                c("gene_id", "total_uorfs")]
    cnt_fig <- table(per_gene$total_uorfs)
    kf <- as.integer(names(cnt_fig)); nf <- as.integer(cnt_fig)
    filt <- data.frame(n_uORFs_per_gene = kf, n_genes = nf,
                       n_TuORFs = kf * nf, set = "Figure_S3_filtered_uIF_gt0")
    .add("AS PLOTTED IN FIGURE S3 (after uIF>0/finite filter + per-coord dedup;")
    .add("this is the filtered subset the panels actually show):")
    for (i in seq_along(kf))
      .add(sprintf("  genes with %d uORF(s): %5d", kf[i], nf[i]))
    .add("  Manuscript S3 legend: n = 1020 (2 uORFs), 320 (3), 109 (4).")
  }

  utils::write.csv(rbind(exh, filt),
                   file.path(OUT_DIR, "Reconcile_A3_uORFs_per_gene.csv"), row.names = FALSE)
}, error = function(e) .add("A3 FAILED: ", conditionMessage(e)))

## ---------------------------------------------------------------------------
## A4) uORF length bins (manuscript Figure 5): n, %removed, medians, 1-aa group
## ---------------------------------------------------------------------------
tryCatch({
  stopifnot(exists("df_uorf"), exists("df_perc"))

  med_all <- df_uorf %>%
    dplyr::group_by(Length_Bin) %>%
    dplyr::summarise(n = dplyr::n(), median_TE = stats::median(TE, na.rm = TRUE),
                     median_log2_uIF = stats::median(log2_uIF, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(panel = "A_all_uORFs")
  med_kept <- df_uorf %>%
    dplyr::filter(pUP > 30) %>%
    dplyr::group_by(Length_Bin) %>%
    dplyr::summarise(n = dplyr::n(), median_TE = stats::median(TE, na.rm = TRUE),
                     median_log2_uIF = stats::median(log2_uIF, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(panel = "B_pUP_gt30")

  .add("")
  .add("--- A4. uORF length bins (manuscript Figure 5) ------------------------")
  .add("Bins used by the script: 1, 2-10, 11-20, 21-30, >30 AA  (a 1-aa bin EXISTS).")
  .add("% removed by the <=30% uORF+ cutoff, per bin:")
  for (i in seq_len(nrow(df_perc)))
    .add(sprintf("  %-6s : total = %5d | removed = %5d | %% removed = %.1f%%",
                 as.character(df_perc$Length_Bin[i]), df_perc$Total[i],
                 df_perc$Removed[i], df_perc$Percent[i]))
  .add("  Manuscript claim: 17/32/48/54%% excluded for 2-10/11-20/21-30/>30 AA;")
  .add("  the Results text does not mention the 1-aa bin, but the Figure 5A-B legend does.")
  .add("Median TE per bin (panel A, all uORFs):")
  for (i in seq_len(nrow(med_all)))
    .add(sprintf("  %-6s : n = %5d | median TE = %.3f | median log2(uIF) = %.3f",
                 as.character(med_all$Length_Bin[i]), med_all$n[i],
                 med_all$median_TE[i], med_all$median_log2_uIF[i]))
  .add("Median TE per bin (panel B, pUP > 30 only):")
  for (i in seq_len(nrow(med_kept)))
    .add(sprintf("  %-6s : n = %5d | median TE = %.3f | median log2(uIF) = %.3f",
                 as.character(med_kept$Length_Bin[i]), med_kept$n[i],
                 med_kept$median_TE[i], med_kept$median_log2_uIF[i]))

  utils::write.csv(dplyr::bind_rows(
      med_all, med_kept,
      df_perc %>% dplyr::transmute(Length_Bin, n = Total, Removed, Percent,
                                   panel = "CD_removed_by_30pct_cutoff")),
    file.path(OUT_DIR, "Reconcile_A4_length_bins.csv"), row.names = FALSE)
}, error = function(e) .add("A4 FAILED: ", conditionMessage(e)))

## ---------------------------------------------------------------------------
## Supporting numbers for the surrounding Results text (items 2, 6 context)
## ---------------------------------------------------------------------------
tryCatch({
  stopifnot(exists("OUT_df6"), exists("OUT_df5_clean"))
  gene_counts_all  <- table(OUT_df5_clean$gene_id)
  n_single_all     <- sum(gene_counts_all == 1)
  n6      <- nrow(OUT_df6)
  n_100   <- sum(OUT_df6$pUP == 100, na.rm = TRUE)
  n_lt30  <- sum(OUT_df6$pUP < 30, na.rm = TRUE)
  n_30_80 <- sum(OUT_df6$pUP >= 30 & OUT_df6$pUP < 80, na.rm = TRUE)
  n_80_100<- sum(OUT_df6$pUP >= 80, na.rm = TRUE)
  n_uIFlt1<- sum(OUT_df6$uIF_raw < 1, na.rm = TRUE)

  .add("")
  .add("--- Supporting numbers (single-TuORF analysis set, OUT_df6) -----------")
  .add("  Genes with exactly 1 TuORF (deduped set, no TE filter): ", n_single_all)
  .add("  Single-TuORF genes analyzed (OUT_df6; TE data + pUP>0 + uIF>0): ", n6)
  .add("  (Manuscript uses 3,213 for this analysis set)")
  .add(sprintf("  Genes with 100%% uORF inclusion : %d / %d  (%.2f%%)  [manuscript: 31/3,213 <1%%]",
               n_100, n6, 100 * n_100 / n6))
  .add(sprintf("  pUP <30 / 30-80 / 80-100        : %d (%.0f%%) / %d (%.0f%%) / %d (%.0f%%)  [manuscript: 30/26/44%%]",
               n_lt30, 100*n_lt30/n6, n_30_80, 100*n_30_80/n6, n_80_100, 100*n_80_100/n6))
  .add(sprintf("  uIF < 1                          : %d / %d  (%.0f%%; median uIF = %.2f)  [manuscript: 62%%, 1,995/3,213; median 0.61]",
               n_uIFlt1, n6, 100 * n_uIFlt1 / n6, stats::median(OUT_df6$uIF_raw, na.rm = TRUE)))
}, error = function(e) .add("Supporting numbers FAILED: ", conditionMessage(e)))

.recon_lines <- c(.recon_lines, "", "================================================================================")
writeLines(.recon_lines, file.path(OUT_DIR, "Reconciliation_A_numbers.txt"))
cat(paste(.recon_lines, collapse = "\n"), "\n")
message("RECONCILIATION MODULE complete -> ", file.path(OUT_DIR, "Reconciliation_A_numbers.txt"))
