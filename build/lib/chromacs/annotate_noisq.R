#!/usr/bin/env Rscript
# ===============================
# R Script: annotate_noisq.R
# ===============================

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(openxlsx)
  library(tools)
  library(data.table)
  library(rtracklayer)
  library(clusterProfiler)
  library(grid)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript annotate_noisq.R <noisq_results.xlsx> <peak_counts.txt> <assembly> <ref_dir>")
}

noisq_xlsx <- args[1]
peak_counts_file <- args[2]
assembly <- args[3]
ref_dir <- args[4]

# ===============================
# Genome mapping and OrgDb logic
# ===============================

genome_mapping <- list(
  GRCh38 = list(sp="homo_sapiens", cap="Homo_sapiens", org_db="org.Hs.eg.db", tax_id=9606),
  GRCm39 = list(sp="mus_musculus", cap="Mus_musculus", org_db="org.Mm.eg.db", tax_id=10090),
  mRatBN7.2 = list(sp="rattus_norvegicus", cap="Rattus_norvegicus", org_db="org.Rn.eg.db", tax_id=10116),
  `ARS-UCD1.3` = list(sp="bos_taurus", cap="Bos_taurus", org_db="org.Bt.eg.db", tax_id=9913),
  Sscrofa11.1 = list(sp="sus_scrofa", cap="Sus_scrofa", org_db="org.Ss.eg.db", tax_id=9825),
  GRCg7b = list(sp="gallus_gallus", cap="Gallus_gallus", org_db="org.Gg.eg.db", tax_id=208526),
  Pan_tro_3.0 = list(sp="pan_troglodytes", cap="Pan_troglodytes", org_db="org.Pt.eg.db", tax_id=9598),
  ROS_Cfam_1.0 = list(sp="canis_lupus_familiaris", cap="Canis_lupus_familiaris", org_db="org.Cf.eg.db", tax_id=9615),
  ARS1 = list(sp="capra_hircus", cap="Capra_hircus", org_db="org.Che.eg.db", tax_id=9925),
  CVASU_BBG_1.0 = list(sp="capra_hircus", cap="Capra_hircus", org_db="org.Che.eg.db", tax_id=9925),
  OryCun2.0 = list(sp="oryctolagus_cuniculus", cap="Oryctolagus_cuniculus", org_db="org.Ocu.eg.db", tax_id=9986),
  gorGor4 = list(sp="gorilla_gorilla", cap="Gorilla_gorilla", org_db="org.Gor.eg.db", tax_id=9593),
  Mmul_10 = list(sp="macaca_mulatta", cap="Macaca_mulatta", org_db="org.Mmu.eg.db", tax_id=9544),
  GRCz11 = list(sp="danio_rerio", cap="Danio_rerio", org_db="org.Dr.eg.db", tax_id=7955),
  UCB_Xtro_10.0 = list(sp="xenopus_tropicalis", cap="Xenopus_tropicalis", org_db="org.Xt.eg.db", tax_id=8364),
  Ssal_v3.1 = list(sp="salmo_salar", cap="Salmo_salar", org_db="org.Ssa.eg.db", tax_id=8030),
  BDGP6 = list(sp="drosophila_melanogaster", cap="Drosophila_melanogaster", org_db="org.Dm.eg.db", tax_id=7227),
  WBcel235 = list(sp="caenorhabditis_elegans", cap="Caenorhabditis_elegans", org_db="org.Ce.eg.db", tax_id=6239)
)

info <- genome_mapping[[assembly]]
if (is.null(info)) stop("No mapping for assembly: ", assembly)

txdb_file <- file.path(ref_dir, paste0("TxDb_", assembly, ".sqlite"))
if (!file.exists(txdb_file)) stop("TxDb not found. Run peak annotation first.")
txdb <- loadDb(txdb_file)

org_pkg <- info$org_db
org_version <- "1.0"
org_tar <- file.path(ref_dir, paste0(org_pkg, "_", org_version, ".tar.gz"))

if (requireNamespace(org_pkg, quietly=TRUE)) {
  suppressPackageStartupMessages(library(org_pkg, character.only=TRUE))
} else {
  if (requireNamespace("BiocManager", quietly=TRUE)) {
    tryCatch({
      BiocManager::install(org_pkg, ask=FALSE, update=FALSE)
      suppressPackageStartupMessages(library(org_pkg, character.only=TRUE))
    }, error = function(e) {
      message("Bioconductor installation failed for OrgDb: ", org_pkg)
    })
  }
  if (!requireNamespace(org_pkg, quietly=TRUE) && file.exists(org_tar)) {
    install.packages(org_tar, repos=NULL, type="source")
    suppressPackageStartupMessages(library(org_pkg, character.only=TRUE))
  }
  if (!requireNamespace(org_pkg, quietly=TRUE)) {
    message("Proceeding without OrgDb package.")
  }
}

# ===============================
# Setup
# ===============================

annot_dir <- file.path(dirname(noisq_xlsx), "Annotated_NOISeq")
dir.create(annot_dir, showWarnings=FALSE, recursive=TRUE)

# Load peak counts
lines <- readLines(peak_counts_file)
lines <- lines[!grepl("^#", lines)]
peak_df <- fread(paste(lines, collapse="\n"), sep="\t", header=TRUE)
required_cols <- c("Geneid", "Chr", "Start", "End")
if (!all(required_cols %in% colnames(peak_df))) {
  stop("Missing required columns in peak_counts.txt. Needed: Geneid, Chr, Start, End")
}
setnames(peak_df, "Geneid", "peak_id")

# ===============================
# Function to annotate and save
# ===============================

annotate_from_table <- function(df, tag) {
  merged <- merge(peak_df, df, by = "peak_id")
  if (nrow(merged) == 0) {
    message("No overlapping peaks for ", tag)
    return()
  }
  gr <- GRanges(seqnames = merged$Chr,
                ranges   = IRanges(start = merged$Start, end = merged$End))
  ann <- annotatePeak(gr, TxDb=txdb, annoDb=org_pkg, tssRegion=c(-3000,3000))
  
  df_ann <- as.data.frame(ann)
  df_ann[] <- lapply(df_ann, function(x) {
    if (is.character(x)) {
      x <- gsub("[\r\n\t]+", " ", x)
      x <- gsub(" +", " ", x)
      trimws(x)
    } else x
  })
  
  out_prefix <- file.path(annot_dir, tag)
  write.table(df_ann, file=paste0(out_prefix, "_annotated.tsv"),
              sep="\t", quote=TRUE, row.names=FALSE)
  
  pdf(paste0(out_prefix, "_summary.pdf"), width=8, height=10)
  vennpie(ann)
  plotAnnoPie(ann)
  
  grid.newpage()
  grid.text("Annotation Summary", x=0.5, y=0.95,
            gp=gpar(fontsize=14, fontface="bold"))
  summary_text <- capture.output(print(ann))
  for (i in seq_along(summary_text)) {
    grid.text(summary_text[i], x=0.02, y=0.9 - (i-1)*0.025,
              just="left", gp=gpar(fontsize=9, family="mono"))
  }
  dev.off()
  message("[✓] Annotated: ", tag)
}

# ===============================
# Main Annotation Tasks
# ===============================

# 1. Annotate main noisq_results.xlsx
main_df <- read.xlsx(noisq_xlsx)
if (!"peak_id" %in% colnames(main_df)) stop("Missing 'peak_id' column in main NOISeq results.")
annotate_from_table(main_df, "noisq")

# 2. Annotate gain/loss files (if present)
gain_file <- file.path(dirname(noisq_xlsx), "NOISeq_gain_sites.tsv")
loss_file <- file.path(dirname(noisq_xlsx), "NOISeq_loss_sites.tsv")

if (file.exists(gain_file)) {
  gain_df <- fread(gain_file)
  if ("peak_id" %in% colnames(gain_df)) {
    annotate_from_table(gain_df, "gain_sites")
  } else {
    message("gain_sites file exists but lacks 'peak_id' column.")
  }
} else {
  message("NOISeq_gain_sites.tsv not found.")
}

if (file.exists(loss_file)) {
  loss_df <- fread(loss_file)
  if ("peak_id" %in% colnames(loss_df)) {
    annotate_from_table(loss_df, "loss_sites")
  } else {
    message("loss_sites file exists but lacks 'peak_id' column.")
  }
} else {
  message("NOISeq_loss_sites.tsv not found.")
}

message("\n✅ All annotation steps complete. Results in: ", annot_dir)
