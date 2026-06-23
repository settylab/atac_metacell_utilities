args <- commandArgs(trailingOnly = TRUE)
peaks.file <- args[1]
outfile    <- args[2]
span       <- as.numeric(args[3])
org        <- args[4]

suppressPackageStartupMessages({
  library(Biostrings)
  library(GenomicRanges)
})

load.bsgenome <- function(genome) {
  if (genome == "hg38") {
    suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
    return(Hsapiens)
  } else if (genome == "hg19") {
    suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
    return(Hsapiens)
  } else if (genome == "mm10") {
    suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm10))
    return(Mmusculus)
  } else {
    stop(paste("Unsupported genome:", genome))
  }
}

regions <- read.table(peaks.file, stringsAsFactors = FALSE, header = FALSE)
colnames(regions) <- c("chrom", "chromStart", "chromEnd", "summit", "score", "name")
all.regions <- GRanges(
  regions[, "chrom"],
  IRanges(regions[, "chromStart"], regions[, "chromEnd"]),
  score  = regions[, "score"],
  summit = regions[, "summit"],
  name   = regions[, "name"]
)

# Identical coord transformation to seq_gl.R lines 57-59:
start(all.regions) <- end(all.regions) <- start(all.regions) + all.regions$summit - 1
all.regions <- resize(all.regions, fix = "center", width = span)

genome <- load.bsgenome(org)
seqs <- getSeq(genome, all.regions)
names(seqs) <- all.regions$name
writeXStringSet(seqs, outfile)
cat(sprintf("wrote %d sequences to %s\n", length(seqs), outfile))
