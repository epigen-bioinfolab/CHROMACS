#!/usr/bin/env Rscript
# ===============================
# R Script: diffbind_3.R
# ===============================

args <- commandArgs(trailingOnly=TRUE)

# =============================================================================#
if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager", repos="https://cran.rstudio.com")
}

biocPkgs <- c("DiffBind", "ggplot2", "rtracklayer")

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

if(length(args) < 2) {
  stop("Usage: Rscript diffbind_3.R <out_csv> <out_dir> [fdr_threshold]")
}
metadata_path <- args[1]
output_dir <- args[2]
fdr_thresh <- ifelse(length(args) >= 3, as.numeric(args[3]), 0.05)
threads <- ifelse(length(args) >= 4, as.integer(args[4]), 4)
if (is.na(threads) || threads <= 0) {
  threads <- 4
}
options(mc.cores=threads)
message("Using ", threads, " threads for parallel processing.")


if(!file.exists(metadata_path)) {
  stop(paste("Metadata file not found:", metadata_path))
}

metadata <- read.csv(metadata_path)


#count replicates per condn
rep_counts <- table(metadata$Condition)
if (length(rep_counts) < 2) {
  stop("Need at least two conditions in metadata")
}
min_reps <- min(rep_counts)
message("Detected replicate counts per condition:\n",
        paste(names(rep_counts), rep_counts, sep="=", collapse=", "),
        "\nUsing minMembers = ", min_reps)


#diffbind analysis
tryCatch({
  dbObj <- dba(sampleSheet = metadata)
  
  # Skip greylist and blacklist analysis (blacklist removal preferably done during peak calling)
  dbObj <- dba.blacklist(dbObj,
                         blacklist = FALSE,
                         greylist  = FALSE)
  
  dbObj <- dba.count(dbObj)
  
  # consensus peaks
  consensus_gr <- dba.peakset(dbObj, bRetrieve=TRUE, DataType=DBA_DATA_GRANGES)
  names(consensus_gr) <- as.character(seq_along(consensus_gr))
  rtracklayer::export(consensus_gr, file.path(output_dir, "consensus_peaks_for_fimo.bed"))
  
  
  dbObj <- dba.contrast(dbObj,
                        categories = DBA_CONDITION,
                        minMembers = min_reps,
                        reorderMeta = list(Condition="baseline")) # ensures the baseline
  
  dbObj <- dba.analyze(dbObj, bBlacklist = FALSE, bGreylist = FALSE)
  
  results_gr <- dba.report(dbObj, th = fdr_thresh)
  
  #peakID
  consensus_gr <- dba.peakset(dbObj, bRetrieve=TRUE, DataType=DBA_DATA_GRANGES)
  consensus_df <- as.data.frame(consensus_gr)
  results_df <- as.data.frame(results_gr)
  coord_key <- function(gr) paste0(seqnames(gr), ":", start(gr), "-", end(gr))
  consensus_coords <- coord_key(consensus_gr)
  results_coords <- coord_key(results_gr)
  peak_id_map <- match(results_coords, consensus_coords)
  results_df <- cbind(PeakID = as.character(peak_id_map), results_df)
  
  # Save
  write.csv(results_df, file.path(output_dir, "diffbind_results.csv"), row.names = FALSE)
  
  
# ============================================================================
  #plots
  pdf(file.path(output_dir, "diffbind_enhanced_plots.pdf"))
  
  #correlation heatmap
  try(dba.plotHeatmap(dbObj, correlations=TRUE))
  
  #pca
  try(dba.plotPCA(dbObj, attributes=DBA_CONDITION))
  
  #ma
  try(dba.plotMA(dbObj))
  
  #volcano
  try(dba.plotVolcano(dbObj))
  
  #box
  try(dba.plotBox(dbObj))
  
  #affinity hm
  try(dba.plotHeatmap(dbObj, contrast=1, correlations=FALSE, 
                      scale="row", colScheme = colorRampPalette(c("blue", "white", "red"))(256)))
  
  dev.off()
  
  # extra outputs
  norm_counts <- dba.peakset(dbObj, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)
  write.csv(norm_counts, file.path(output_dir, "normalized_counts.csv"))
  
  gr <- dba.report(dbObj, th=fdr_thresh, DataType=DBA_DATA_GRANGES)
  rtracklayer::export(gr[gr$Fold > 0], file.path(output_dir, "gain_sites.bed"))
  rtracklayer::export(gr[gr$Fold < 0], file.path(output_dir, "loss_sites.bed"))
  
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
  
}, error = function(e) {
  stop(paste("DiffBind analysis failed:", e$message))
})
