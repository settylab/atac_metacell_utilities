library(tidyverse)
library(data.table)
library(Matrix)
library(SummarizedExperiment)
library(edgeR)
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3){
    stop("Please supply the directory containing the required data files, cell type variable name, and a vector of tuples of cell types to compare")
} else {
    data_dir <- args[1] #'celltype_combined'
    to_compare <- args[2]  #c("Mono",  'HSC')
    group_variable <- args[3]
}
to_compare = strsplit(to_compare, as.character('/'))[[1]]
to_compare = strsplit(to_compare, as.character(','), fixed = TRUE)

print(to_compare)
fname <- paste(data_dir, '/meta_atac_metadata.csv',  sep = "")
print(fname)
atac_metadata <- fread(paste(data_dir, '/meta_atac_metadata.csv' , sep = ""))

peak_names <- fread(paste(data_dir, '/meta_atac_peaks.csv',  sep = ""))
metacell_mtx <- readMM(paste(data_dir, '/meta_atac_X.mtx', sep = ""))
# No user input below this line
#####################################
peak_names <- unlist(peak_names$V2)
peak_names <- peak_names[! peak_names %in% c('0')]


metacell_mtx <- as(metacell_mtx, "dgCMatrix")
metacell_mtx <- t(metacell_mtx)


metacell_se <- SummarizedExperiment(assays=metacell_mtx)
colnames(metacell_se) <- atac_metadata$V1
diff_acc_edgeR <- function(metacell_mtx, atac_metadata, group_variable, to_compare, out_dir){
    sub_atac_md <- atac_metadata %>% 
      .[,group:=eval(as.name(group_variable))] %>% 
      .[group%in%to_compare] %>%
      .[,group:=factor(group, levels=to_compare)] %>% 
      setorder(group)
    sub_se <- metacell_se[,sub_atac_md$V1]
    colData(sub_se) <- sub_atac_md %>% tibble::column_to_rownames("V1") %>% DataFrame
    rownames(sub_se) <- peak_names
    assayNames(sub_se)[1] <- "counts"
    assay(sub_se,"logcounts") <- log(1e6*(sweep(assay(sub_se,"counts"),2,colSums(assay(sub_se,"counts")),"/"))+1)


    acc.dt <- data.table(feature = rownames(sub_se),
                         detection_rate_groupA = rowMeans(assay(sub_se,"counts")[,sub_se$celltype_combined==to_compare[1]]>0) %>% round(2),
                         detection_rate_groupB = rowMeans(assay(sub_se,"counts")[,sub_se$celltype_combined==to_compare[2]]>0) %>% round(2),
                         mean_groupA = rowMeans(assay(sub_se[,sub_se$celltype_combined==to_compare[1]],"logcounts")) %>% round(2),
                         mean_groupB = rowMeans(assay(sub_se[,sub_se$celltype_combined==to_compare[2]],"logcounts")) %>% round(2)
    )

    features.to.use <- acc.dt[detection_rate_groupA>=0.3 | detection_rate_groupB>=0.3,feature]

    atac_metacells.dge <- DGEList(assay(sub_se[features.to.use,],"counts"))


    cdr <- colMeans(assay(sub_se,"counts")>0)
    design <- model.matrix(~cdr+sub_se$group)

    # Estimate dispersions
    atac_metacells.dge  <- estimateDisp(atac_metacells.dge,design)

    # Fit GLM
    fit <- glmQLFit(atac_metacells.dge,design)

    # Likelihood ratio test
    lrt <- glmQLFTest(fit)

    out <- topTags(lrt, n=nrow(lrt))$table %>% as.data.table(keep.rownames=T) %>%
      setnames(c("feature","logFC","logCPM","LR","p.value","padj_fdr")) %>%
      .[,c("logCPM","LR","p.value"):=NULL] %>%
      .[,c("groupA_N","groupB_N"):=list(table(sub_atac_md$group)[1],table(sub_atac_md$group)[2])]%>% 
      .[,c("padj_fdr","logFC"):=list(signif(padj_fdr,digits=3),round(logFC,3))] %>%
      merge(acc.dt[,c("feature","mean_groupA","mean_groupB")], by="feature", all.y=TRUE) %>%
      .[is.na(logFC),c("logFC","padj_fdr"):=list(0,1)] %>%
      setorder(padj_fdr, na.last=T)

    fwrite(out, sprintf("%s/%s_%s_diff_acc.tsv", data_dir, to_compare[1], to_compare[2]), sep="\t", na="NA", quote=F)


    sigA_peaks <- out[out$logFC <= -1.25 & out$padj_fdr <=0.01,"feature"] %>% unlist()
    sigB_peaks <- out[out$logFC >= 1.25 & out$padj_fdr <=0.01,"feature"]  %>% unlist()


    pdf(file=sprintf("%s/%s_%s_MA.pdf", data_dir, to_compare[1], to_compare[2]))

    with(lrt$table, plot(logCPM,logFC,pch=16,cex=0.2, col="gray"))

    with(lrt$table[sigA_peaks,], points(logCPM,logFC,pch=16, cex=0.5,col="red"))
    with(lrt$table[sigB_peaks,], points(logCPM,logFC,pch=16, cex=0.5,col="dodgerblue"))
    legend("bottomleft",legend=c(sprintf("%s peaks", to_compare[1]),sprintf("%s peaks", to_compare[2])), pch=16,col=c("red","dodgerblue"))

    dev.off()
    }

for (comparison in to_compare){
    diff_acc_edgeR(metacell_mtx, atac_metadata, group_variable, as.vector(comparison), data_dir)
}

