args <-commandArgs(trailingOnly=TRUE)

library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)

peak_file <- args[1]
chromvar_indir <- args[2]
chromvar_outdir <- args[3]

peaks <- getPeaks(peak_file, sort_peaks=TRUE, extra_cols=c('summit'=4, 'score'=5,'name'=6))
count_matrix <- Matrix::readMM(sprintf("%s/sc_atac_counts.mtx", chromvar_indir))
count_matrix <- t(count_matrix)

nfrags <- readr::read_csv(sprintf("%s/sc_nfrags.csv", chromvar_indir))
col_data <- DataFrame(nfrags$nFrags, row.names=nfrags$x)
names(col_data) <- 'depth'

fragment_counts <- SummarizedExperiment(assays=list(counts=count_matrix), rowRanges=peaks, colData=col_data)
fragment_counts <- addGCBias(fragment_counts, genome=BSgenome.Hsapiens.UCSC.hg38)
counts_filtered <- filterPeaks(fragment_counts, non_overlapping=TRUE)

tf_names <- readr::read_csv(sprintf("%s/chromvar_tf_names.csv", chromvar_indir))
insilico_chip_mat <- Matrix::readMM(sprintf("%s/bin_ins_chip.mtx", chromvar_indir))

anno <- getAnnotations(insilico_chip_mat, rowRanges=rowRanges(counts_filtered))
colnames(anno) <-tf_names$tf_name

bg <- getBackgroundPeaks(object=counts_filtered)
dev <- computeDeviations(object=counts_filtered, annotations=anno, background_peaks=bg)

# export output
dir.create(sprintf("%s/export", chromvar_dir))
deviations <- assay(dev, "deviations")
write.csv(deviations, sprintf("%s/deviations.csv", chromvar_outdir))
zs_deviations <- assay(dev, "z")
write.csv(zs_deviations, sprintf("%s/deviations_zs.csv", chromvar_outdir))