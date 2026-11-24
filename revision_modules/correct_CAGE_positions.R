################################################################################
# correct_CAGE_positions.R — position register of the CAGE-seq input file (2026-11-22)                                      #
#                                                                                                                            #
# The CAGE-seq signal file used until V52 (Files/wt_R123_genome_RiboPlotR.txt, made from the pooled track of the three      #
# wild-type replicates of Thieffry et al. 2020 by taking the end of each interval) places every position 1 nt DOWNSTREAM    #
# of the first transcribed nucleotide, in the direction of transcription, on both strands. The TSS-seq file                 #
# (Files/wt_TSS_Nielsen_2019.txt) is in the correct register. Evidence (printed again below):                               #
#   (1) a transcription start site has a pyrimidine at -1 and a purine at +1; in the TSS-seq file this pair is at the       #
#       reported position, in the CAGE-seq file it is one nucleotide upstream of the reported position;                     #
#   (2) the two files share most signal when the CAGE-seq positions are moved 1 nt toward the 5' end.                       #
# With the uncorrected file a transcript that starts exactly at the A of a uORF AUG was counted as lacking the uORF with    #
# CAGE-seq and as containing it with TSS-seq.                                                                               #
#                                                                                                                            #
# This script keeps the file as received under the name Files/wt_R123_genome_RiboPlotR_before_position_correction.txt and   #
# writes the corrected file under the name that all scripts read, Files/wt_R123_genome_RiboPlotR.txt:                       #
#   plus strand: position - 1; minus strand: position + 1. Signal values and the order of the lines are unchanged.          #
# It is safe to run again: the corrected file is always made from the file as received.                                     #
# Output: Files/wt_R123_genome_RiboPlotR.txt (corrected), output/CAGE_position_correction_summary.txt                       #
# Run from Weak_uORF_submission/, before the main script:  Rscript revision_modules/correct_CAGE_positions.R                #
################################################################################
suppressPackageStartupMessages({ library(data.table); library(Biostrings); library(tools) })
F_USED <- "Files/wt_R123_genome_RiboPlotR.txt"; F_RAW <- "Files/wt_R123_genome_RiboPlotR_before_position_correction.txt"
F_TSS <- "Files/wt_TSS_Nielsen_2019.txt"; F_FA <- "Files/TAIR10_chr_all_2.fas"
MD5_RAW <- "d2fc87faddba5adb6f407d45c5e01aa2"   # the file as received (used until V52)
.log <- character(0); .add <- function(...) { s <- sprintf(...); message(s); .log <<- c(.log, s) }
.add("CAGE position correction — %s", format(Sys.time(), "%Y-%m-%d %H:%M"))
if (!file.exists(F_RAW)) { stopifnot(md5sum(F_USED) == MD5_RAW); stopifnot(file.copy(F_USED, F_RAW)); .add("the file as received was copied to %s", F_RAW) }
stopifnot(md5sum(F_RAW) == MD5_RAW)
raw <- fread(F_RAW, header = FALSE, sep = "\t", colClasses = c("character", "character", "integer", "character"), col.names = c("count", "chr", "pos", "strand"))
stopifnot(all(raw$strand %in% c("+", "-")), !anyNA(raw$pos))
new <- copy(raw); new[, pos := ifelse(strand == "+", pos - 1L, pos + 1L)]
fwrite(new, F_USED, sep = "\t", col.names = FALSE, quote = FALSE)
.add("%s positions; plus strand %s (position - 1), minus strand %s (position + 1); signal values and line order unchanged", format(nrow(new), big.mark = ","), format(sum(new$strand == "+"), big.mark = ","), format(sum(new$strand == "-"), big.mark = ","))
.add("written: %s (md5 %s); file as received: %s (md5 %s)", F_USED, md5sum(F_USED), F_RAW, md5sum(F_RAW))

## evidence: base composition around the reported positions, and agreement with TSS-seq, before and after
g <- readDNAStringSet(F_FA); names(g) <- sub("^[Cc]hr", "", sub(" .*", "", names(g)))
tss <- fread(F_TSS, header = FALSE, sep = "\t", col.names = c("count", "chr", "pos", "strand"))
prep <- function(d) { d <- copy(d); d[, count := as.numeric(count)]; d[, chr := sub("^[Cc]hr", "", as.character(chr))]; d[chr %in% as.character(1:5)] }
base_at <- function(d, off) { out <- character(nrow(d))
  for (ch in unique(d$chr)) { s <- g[[ch]]; i <- which(d$chr == ch); p <- ifelse(d$strand[i] == "+", d$pos[i] + off, d$pos[i] - off); ok <- p >= 1 & p <= length(s)
    b <- rep(NA_character_, length(i)); b[ok] <- strsplit(as.character(s[p[ok]]), "")[[1]]; m <- d$strand[i] == "-"; b[m] <- chartr("ACGT", "TGCA", b[m]); out[i] <- b }
  out }
share <- function(d, off, set) { b <- base_at(d, off); sum(d$count[b %in% set], na.rm = TRUE) / sum(d$count[!is.na(b)]) }
shared <- function(a, b) { a <- copy(a); b <- copy(b); a[, f := count / sum(count)]; b[, f := count / sum(count)]; m <- merge(a[, .(chr, pos, strand, fa = f)], b[, .(chr, pos, strand, fb = f)], by = c("chr", "pos", "strand")); sum(pmin(m$fa, m$fb)) }
R <- prep(raw); N <- prep(new); T <- prep(tss)
.add("")
.add("Share of the signal with a pyrimidine 1 nt upstream of the reported position and a purine at the reported position (transcript orientation):")
for (x in list(list("TSS-seq", T), list("CAGE-seq as received", R), list("CAGE-seq corrected", N)))
  .add("  %-22s pyrimidine at -1: %.1f%%   purine at the position: %.1f%%   (purine at -1: %.1f%%)", x[[1]], 100 * share(x[[2]], -1, c("C", "T")), 100 * share(x[[2]], 0, c("A", "G")), 100 * share(x[[2]], -1, c("A", "G")))
.add("Signal shared with the TSS-seq file (sum over positions of the smaller of the two normalised signals): as received %.3f, corrected %.3f", shared(R, T), shared(N, T))
writeLines(.log, "output/CAGE_position_correction_summary.txt")
