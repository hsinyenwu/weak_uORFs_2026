#!/usr/bin/env Rscript
################################################################################
# psite_extract.R — P-site counts per genomic position from a Ribo-seq BAM       #
# (R port of psite_extract.py, 2026-10-28; same rules, same output)             #
#                                                                               #
# usage: Rscript psite_extract.R BAM RULE OUT.tsv[.gz] [--roi roi.bed]          #
#                                 [--mapq Q] [--yield N]                        #
#   RULE  A = 24-28 nt, offsets 8-12 (RiboTaper settings used in this study)    #
#         B = Ribo-seQC cutoffs for the merged libraries (23-30 nt)              #
#         C = length - 16 for 21-36 nt (3'-anchored)      D = fixed 12 (20-36)   #
#   --roi  keep alignments overlapping the BED regions (like samtools view -L)  #
#   --mapq keep alignments with MAPQ >= Q (e.g. 255 = STAR unique mappers)      #
#   default: all primary alignments (secondary and supplementary removed)       #
# The offset is walked along the aligned bases from the read's 5' end; spliced   #
# gaps (N) and deletions (D) are skipped, insertions and clips do not count.    #
# Output columns (no header): chr, pos (1-based P-site), strand, count.         #
################################################################################
suppressPackageStartupMessages({ library(GenomicAlignments); library(Rsamtools); library(GenomicRanges); library(data.table) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: Rscript psite_extract.R BAM RULE OUT.tsv[.gz] [--roi roi.bed] [--mapq Q] [--yield N]")
bam <- args[1]; rule <- args[2]; out <- args[3]
opt <- function(flag, default) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1] }
roi_path <- opt("--roi", NA); mapq_min <- as.numeric(opt("--mapq", NA)); yield <- as.integer(opt("--yield", "5000000"))
RULES <- list(A = c(`24` = 8, `25` = 9, `26` = 10, `27` = 11, `28` = 12),
              B = c(`23` = 7, `24` = 8, `25` = 9, `26` = 10, `27` = 11, `28` = 12, `29` = 12, `30` = 13),
              C = setNames(21:36 - 16, 21:36), D = setNames(rep(12, 17), 20:36))
off <- RULES[[rule]]; if (is.null(off)) stop("unknown rule: ", rule)
roi <- NULL
if (!is.na(roi_path)) { b <- fread(roi_path, header = FALSE); roi <- reduce(GRanges(as.character(b$V1), IRanges(b$V2 + 1L, b$V3)), ignore.strand = TRUE) }

psites_chunk <- function(ga, off) {
  offv <- off[as.character(qwidth(ga))]; keep <- !is.na(offv); if (!any(keep)) return(NULL)
  ga <- ga[keep]; offv <- as.integer(offv[keep])
  bl <- grglist(ga, drop.D.ranges = TRUE)                       # aligned blocks (M/=/X); N and D split blocks
  ul <- unlist(bl, use.names = FALSE)
  d <- data.table(rid = rep(seq_along(bl), elementNROWS(bl)), chr = as.character(seqnames(ul)), s = start(ul), e = end(ul), strand = as.character(strand(ul)))
  d[, w := e - s + 1L]; d[, ord := ifelse(strand == "+", s, -e)]; setorder(d, rid, ord)   # blocks in 5'->3' order of the read
  d[, cb := cumsum(w) - w, by = rid]; d[, o := offv[rid]]
  hit <- d[o >= cb & o < cb + w]
  hit[, pos := ifelse(strand == "+", s + (o - cb), e - (o - cb))]
  hit[, .(N = .N), by = .(chr, pos, strand)]
}

bf <- BamFile(bam, yieldSize = yield); open(bf)
param <- ScanBamParam(what = c("flag", "mapq"))
acc <- list(); n_in <- 0; n_kept <- 0; n_used <- 0; t0 <- Sys.time()
repeat {
  ga <- readGAlignments(bf, param = param); if (length(ga) == 0) break
  n_in <- n_in + length(ga)
  fl <- mcols(ga)$flag; keep <- bitwAnd(fl, 0x904L) == 0L                       # mapped, primary, not supplementary
  if (!is.na(mapq_min)) keep <- keep & mcols(ga)$mapq >= mapq_min
  ga <- ga[keep]
  if (!is.null(roi)) ga <- ga[overlapsAny(granges(ga), roi, ignore.strand = TRUE)]
  n_kept <- n_kept + length(ga)
  r <- psites_chunk(ga, off); if (!is.null(r)) { acc[[length(acc) + 1L]] <- r; n_used <- n_used + sum(r$N) }
  message(sprintf("  %s records read, %s kept, %s P-sites (%.1f min)", format(n_in, big.mark = ","), format(n_kept, big.mark = ","), format(n_used, big.mark = ","), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
close(bf)
res <- rbindlist(acc)[, .(count = sum(N)), by = .(chr, pos, strand)]
setorder(res, chr, strand, pos)
fwrite(res[, .(chr, pos, strand, count)], out, sep = "\t", col.names = FALSE)
message(sprintf("reads in: %s  kept: %s  used: %s  positions: %s  -> %s", format(n_in, big.mark = ","), format(n_kept, big.mark = ","), format(n_used, big.mark = ","), format(nrow(res), big.mark = ","), out))
