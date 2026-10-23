#args <-commandArgs(trailingOnly=TRUE)
libloc <- snakemake@input[["libloc"]]
install.packages("renv", lib=libloc, repos='http://cran.us.r-project.org')
