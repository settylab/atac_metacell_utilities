args <-commandArgs(trailingOnly=TRUE)
library(SeqGL)

span <- args[3]
org <- args[4]

# Peaks directory
peaks.file <- args[1]

# Load peaks 
regions <- read.table(peaks.file, stringsAsFactors=FALSE, header=TRUE)
all.regions <- GRanges(regions[,'chrom'], IRanges (regions[,'chromStart'], regions[,'chromEnd']),
			    score=regions[,'score'], summit=regions[,'summit'], name=regions[,'name'])

start(all.regions) <- end(all.regions) <- start(all.regions) + all.regions$summit - 1

all.regions <- resize(all.regions, fix='center', span)

# Identify sequences and flag any sequences with N
seqs <- SeqGL:::get.seqs(SeqGL:::load.bsgenome(org), all.regions)
names(seqs) <- all.regions$name

# Save 
writeXStringSet(seqs, args[2])
