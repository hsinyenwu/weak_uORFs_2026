################################################################################
# compare_output_runs.R — compare two output folders of the canonical script    #
# (e.g. the run before and after the class renaming). Reports, per file:         #
#   xlsx: sheets, dimensions, columns whose values differ (numeric: max |diff|;   #
#         character: number of differing cells and an example pair)             #
#   csv/txt: identical or not (and the first differing lines)                   #
#   pdf: pdftotext word sets that differ (labels only, expected after renaming) #
# usage: Rscript revision_modules/compare_output_runs.R OLD_DIR NEW_DIR [pattern]#
################################################################################
suppressPackageStartupMessages({ library(openxlsx) })
args <- commandArgs(trailingOnly = TRUE); OLD <- args[1]; NEW <- args[2]; PAT <- if (length(args) >= 3) args[3] else NULL
stopifnot(dir.exists(OLD), dir.exists(NEW))
files <- intersect(list.files(OLD), list.files(NEW)); if (!is.null(PAT)) files <- files[grepl(PAT, files)]
files <- files[!grepl("^A[0-9]+_|^FigS2_|^FigureS2_|^Fig3C|^GO_FDR|^A9_|_run\\.log$", files)]     # revision-module outputs are not produced by the canonical script
only_old <- setdiff(list.files(OLD), list.files(NEW)); only_new <- setdiff(list.files(NEW), list.files(OLD))
only_old <- only_old[!grepl("^A[0-9]+_|^FigS2_|^FigureS2_|^Fig3C|^GO_FDR|^A9_", only_old)]
cat(sprintf("files compared: %d | only in OLD: %s | only in NEW: %s\n\n", length(files), paste(only_old, collapse = ", "), paste(only_new, collapse = ", ")))
cmp_xlsx <- function(a, b) {
  sa <- getSheetNames(a); sb <- getSheetNames(b); out <- character(0)
  if (!identical(sa, sb)) out <- c(out, sprintf("  sheets differ: %s vs %s", paste(sa, collapse = "/"), paste(sb, collapse = "/")))
  for (sh in intersect(sa, sb)) {
    x <- read.xlsx(a, sheet = sh); y <- read.xlsx(b, sheet = sh)
    if (!identical(dim(x), dim(y))) { out <- c(out, sprintf("  [%s] dimensions %s vs %s", sh, paste(dim(x), collapse = "x"), paste(dim(y), collapse = "x"))); next }
    if (!identical(names(x), names(y))) out <- c(out, sprintf("  [%s] column names differ: %s", sh, paste(setdiff(union(names(x), names(y)), intersect(names(x), names(y))), collapse = ", ")))
    for (cn in intersect(names(x), names(y))) {
      u <- x[[cn]]; v <- y[[cn]]
      if (is.numeric(u) && is.numeric(v)) { d <- suppressWarnings(max(abs(u - v), na.rm = TRUE)); if (is.finite(d) && d > 1e-9) out <- c(out, sprintf("  [%s] %s: numeric max |diff| = %.3g", sh, cn, d)); if (!identical(is.na(u), is.na(v))) out <- c(out, sprintf("  [%s] %s: NA pattern differs", sh, cn)) }
      else { same <- (is.na(u) & is.na(v)) | (!is.na(u) & !is.na(v) & as.character(u) == as.character(v)); if (!all(same)) { i <- which(!same)[1]; out <- c(out, sprintf("  [%s] %s: %d cells differ, e.g. '%s' -> '%s'", sh, cn, sum(!same), substr(as.character(u[i]), 1, 60), substr(as.character(v[i]), 1, 60))) } }
    }
  }
  if (length(out) == 0) "  identical content" else out
}
cmp_text <- function(a, b) { x <- readLines(a, warn = FALSE); y <- readLines(b, warn = FALSE)
  if (identical(x, y)) return("  identical"); i <- which(x != y | is.na(x) | is.na(y))[1]
  c(sprintf("  %d vs %d lines; first difference at line %s:", length(x), length(y), i), sprintf("    OLD: %s", substr(x[i], 1, 150)), sprintf("    NEW: %s", substr(y[i], 1, 150))) }
cmp_pdf <- function(a, b) {
  w <- function(f) { t <- tryCatch(system2("pdftotext", c("-layout", shQuote(f), "-"), stdout = TRUE, stderr = FALSE), error = function(e) character(0)); table(unlist(strsplit(gsub("\\s+", " ", paste(t, collapse = " ")), " "))) }
  ta <- w(a); tb <- w(b); ka <- names(ta); kb <- names(tb)
  d_old <- setdiff(ka, kb); d_new <- setdiff(kb, ka); same <- intersect(ka, kb); cnt <- same[ta[same] != tb[same]]
  if (!length(d_old) && !length(d_new) && !length(cnt)) return("  text identical")
  c(sprintf("  words only in OLD: %s", paste(head(d_old, 12), collapse = " | ")), sprintf("  words only in NEW: %s", paste(head(d_new, 12), collapse = " | ")), if (length(cnt)) sprintf("  words with changed counts: %s", paste(head(cnt, 10), collapse = " | ")) else NULL) }
for (f in files) {
  a <- file.path(OLD, f); b <- file.path(NEW, f); ext <- tolower(tools::file_ext(f))
  res <- if (ext == "xlsx") cmp_xlsx(a, b) else if (ext %in% c("csv", "txt", "tsv")) cmp_text(a, b) else if (ext == "pdf") cmp_pdf(a, b) else sprintf("  (%s not compared; sizes %d vs %d)", ext, file.size(a), file.size(b))
  cat(f, "\n", paste(res, collapse = "\n"), "\n", sep = "")
}
