#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

# Install missing Bioconductor packages
if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager", repos="https://cran.rstudio.com")
}

# List of Bioconductor packages needed
biocPkgs <- c("DiffBind", "ggplot2", "rtracklayer")

# Install any that are missing
for (pkg in biocPkgs) {
  if (!requireNamespace(pkg, quietly=TRUE)) {
    BiocManager::install(pkg, ask=FALSE, update=FALSE)
  }
}
suppressPackageStartupMessages({
  library(DiffBind)
  library(ggplot2)
  library(rtracklayer)
})

#==============================================================================#
# Input validation
if(length(args) < 2) {
  stop("Usage: Rscript diffbind_3.R <out_csv> <out_dir> [fdr_threshold]")
}
metadata_path <- args[1]
output_dir <- args[2]
fdr_thresh <- ifelse(length(args) >= 3, as.numeric(args[3]), 0.05)
threads <- ifelse(length(args) >= 4, as.integer(args[4]), 8)
if (is.na(threads) || threads <= 0) {
  threads <- 8
}
options(mc.cores=threads)
message("Using ", threads, " threads for parallel processing.")


if(!file.exists(metadata_path)) {
  stop(paste("Metadata file not found:", metadata_path))
}

metadata <- read.csv(metadata_path)


# Count how many replicates per condition
rep_counts <- table(metadata$Condition)
if (length(rep_counts) < 2) {
  stop("Need at least two conditions in metadata")
}
min_reps <- min(rep_counts)
message("Detected replicate counts per condition:\n",
        paste(names(rep_counts), rep_counts, sep="=", collapse=", "),
        "\nUsing minMembers = ", min_reps)


# DiffBind analysis
tryCatch({
  dbObj <- dba(sampleSheet = metadata)
  
  # Skip greylist and blacklist analysis (blacklist removal preferably done during peak calling)
  dbObj <- dba.blacklist(dbObj,
                         blacklist = FALSE,
                         greylist  = FALSE)
  
  dbObj <- dba.count(dbObj)
  
  # do we add a normalization step here? (DiffBind performs its default normalization (TMM))
  
  # Auto-create exactly one contrast, allowing groups down to min_reps
  dbObj <- dba.contrast(dbObj,
                        categories = DBA_CONDITION,
                        minMembers = min_reps)
  
  dbObj <- dba.analyze(dbObj, bBlacklist = FALSE, bGreylist = FALSE)
  
  # Save results
  results <- dba.report(dbObj, th = fdr_thresh)
  write.csv(as.data.frame(results), 
            file.path(output_dir, "diffbind_results.csv"))
  
  
# ============================================================================
  # Generate plots
  pdf(file.path(output_dir, "diffbind_enhanced_plots.pdf"))
  
  # 1. Correlation Heatmap
  try(dba.plotHeatmap(dbObj, correlations=TRUE))
  
  # 2. PCA Plot
  try(dba.plotPCA(dbObj, attributes=DBA_CONDITION))
  
  # 3. MA Plot
  try(dba.plotMA(dbObj))
  
  # 4. Volcano Plot
  try(dba.plotVolcano(dbObj))
  
  # 5. Boxplots
  try(dba.plotBox(dbObj))
  
  # 6. Binding Affinity Heatmap
  try(dba.plotHeatmap(dbObj, contrast=1, correlations=FALSE, 
                      scale="row", colScheme = colorRampPalette(c("blue", "white", "red"))(256)))
  
  dev.off()
  
  # Additional Outputs --------------------------------------------------------
  # Save normalized read counts
  norm_counts <- dba.peakset(dbObj, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)
  write.csv(norm_counts, file.path(output_dir, "normalized_counts.csv"))
  
  # Save DB sites as BED files
  gr <- dba.report(dbObj, th=fdr_thresh, DataType=DBA_DATA_GRANGES)
  rtracklayer::export(gr[gr$Fold > 0], file.path(output_dir, "gain_sites.bed"))
  rtracklayer::export(gr[gr$Fold < 0], file.path(output_dir, "loss_sites.bed"))
  
  # Save session info
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
  
}, error = function(e) {
  stop(paste("DiffBind analysis failed:", e$message))
})
