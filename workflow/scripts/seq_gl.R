args <-commandArgs(trailingOnly=TRUE)
library(SeqGL)

span <- as.numeric(args[3])
org <- args[4]

# Peaks directory
peaks.file <- args[1]


# Load genome function 
load.bsgenome <- function (genome) {

  if (!genome %in% c('hg38', 'hg19', 'hg18', 'mm9', 'mm10', 'pFal')) {
    stop ('The specified organism is not supported. Install the necessary Bioconductor package and update the script seq_gl.R')
  }

  if (genome == 'hg38') {
    library (BSgenome.Hsapiens.UCSC.hg38)
    org <- Hsapiens
  }

  if (genome == 'hg19') {
    library (BSgenome.Hsapiens.UCSC.hg19)
    org <- Hsapiens
  }

  if (genome == 'hg18') {
    library (BSgenome.Hsapiens.UCSC.hg18)
    org <- Hsapiens
  }

  if (genome == 'mm9') {
    library (BSgenome.Mmusculus.UCSC.mm9)
    org <- Mmusculus
  }

  if (genome == 'mm10') {
    library (BSgenome.Mmusculus.UCSC.mm10)
    org <- Mmusculus
  }

  if (genome == 'pFal'){
    library (BSgenome.Pfalciparum.PlasmoDB.plasFalc3D76)
    org <- Pfalciparum
  }

  return (org)
}

# Load peaks 
regions <- read.table(peaks.file, stringsAsFactors=FALSE, header=FALSE)
colnames(regions)<-c('chrom', 'chromStart', 'chromEnd', 'summit', 'score', 'name')
all.regions <- GRanges(regions[,'chrom'], IRanges (regions[,'chromStart'], regions[,'chromEnd']),
			    score=regions[,'score'], summit=regions[,'summit'], name=regions[,'name'])

start(all.regions) <- end(all.regions) <- start(all.regions) + all.regions$summit - 1

all.regions <- resize(all.regions, fix='center', span)

# Identify sequences and flag any sequences with N
seqs <- SeqGL:::get.seqs(load.bsgenome(org), all.regions)
names(seqs) <- all.regions$name

# Save 
writeXStringSet(seqs, args[2])
