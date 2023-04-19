args <-commandArgs(trailingOnly=TRUE)
library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)

peak_file <- args[1]
sc_counts <- args[2]
n_frag_file <- args[3]
bin_ins_chip_file <- args[4]
tf_name_file <- args[5]
outdir <- args[6]

peaks <- getPeaks(peak_file, sort_peaks=TRUE, extra_cols=c('summit'=4, 'score'=5,'name'=6))
count_matrix <- Matrix::readMM(sc_counts)
count_matrix <- t(count_matrix)

nfrags <- readr::read_csv(n_frag_file)
col_data <- DataFrame(nfrags$nFrags, row.names=nfrags$x)
names(col_data) <- 'depth'

fragment_counts <- SummarizedExperiment(assays= list(counts=count_matrix), rowRanges=peaks, colData = col_data)

fragment_counts <- addGCBias(fragment_counts, genome= BSgenome.Hsapiens.UCSC.hg38)
#counts_filtered <- filterSamples(fragment_counts, shiny=FALSE)
counts_filtered <- filterPeaks(counts_filtered, non_overlapping = TRUE)

# motifs <- getJasparMotifs()
# motif_ix <- matchMotifs(motifs, counts_filtered, genome=BSgenome.Hsapiens.UCSC.hg38)

tf_names <- readr::read_csv(tf_name_file)
insilico_chip_mat <- Matrix::readMM(bin_ins_chip_file)
anno <- getAnnotations(insilico_chip_mat, rowRanges = rowRanges(counts_filtered))
colnames(anno) <-tf_names$tf_name

bg <- getBackgroundPeaks(object = counts_filtered)
dev <- computeDeviations(object = counts_filtered, annotations = anno, background_peaks=bg)

# export output
deviations <- assay(dev, "deviations")
write.csv(deviations,sprintf("%s/cv_dev.csv", outdir))
zs_deviations <- assay(dev, "z")
write.csv(zs_deviations,sprintf("%s/cv_dev_zs.csv", outdir))

# Run vanilla version
dev_vanilla <- computeDeviations(object = counts_filtered, annotations = motif_ix, background_peaks=bg)
van_deviations <- assay(dev_vanilla, "deviations")
write.csv(van_deviations, "/fh/fast/setty_m/user/cdien/insilico_chip/output/van_cv_dev.csv")

van_zs_deviations <- assay(dev_vanilla, "z")
write.csv(van_zs_deviations, "/fh/fast/setty_m/user/cdien/insilico_chip/output/van_cv_dev_zs.csv")
