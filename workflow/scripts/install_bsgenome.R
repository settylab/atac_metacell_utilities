args <-commandArgs(trailingOnly=TRUE)
package <- paste0("bioc::BSgenome.Hsapiens.UCSC.", args[1])
message(package)
renv::install(package)
