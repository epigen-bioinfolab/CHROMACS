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

#genome_mapping
genome_mapping <- list(
  GRCh38 = list(sp="homo_sapiens", cap="Homo_sapiens", org_db="org.Hs.eg.db", tax_id=9606),
  GRCm39 = list(sp="mus_musculus", cap="Mus_musculus", org_db="org.Mm.eg.db", tax_id=10090),
  mRatBN7.2 = list(sp="rattus_norvegicus", cap="Rattus_norvegicus", org_db="org.Rn.eg.db", tax_id=10116),
  `ARS-UCD1.3` = list(sp="bos_taurus", cap="Bos_taurus", org_db="org.Bt.eg.db", tax_id=9913),
  Sscrofa11.1 = list(sp="sus_scrofa", cap="Sus_scrofa", org_db="org.Ss.eg.db", tax_id=9825),
  GRCg7b = list(sp="gallus_gallus", cap="Gallus_gallus", org_db="org.Gg.eg.db", tax_id=208526),
  Pan_tro_3.0 = list(sp="pan_troglodytes", cap="Pan_troglodytes", org_db="org.Pt.eg.db", tax_id=9598),
  ROS_Cfam_1.0 = list(sp="canis_lupus_familiaris", cap="Canis_lupus_familiaris", org_db="org.Cf.eg.db", tax_id=9615),
  ARS1 = list(sp="capra_hircus", cap="Capra_hircus", org_db="org.Che.eg.db", tax_id=9925), # no such org_db found
  CVASU_BBG_1.0 = list(sp="capra_hircus", cap="Capra_hircus", org_db="org.Che.eg.db", tax_id=9925), # no such org_db found
  OryCun2.0 = list(sp="oryctolagus_cuniculus", cap="Oryctolagus_cuniculus", org_db="org.Ocu.eg.db", tax_id=9986), # no such org_db found
  gorGor4 = list(sp="gorilla_gorilla", cap="Gorilla_gorilla", org_db="org.Gor.eg.db", tax_id=9593), # no such org_db found
  Mmul_10 = list(sp="macaca_mulatta", cap="Macaca_mulatta", org_db="org.Mmu.eg.db", tax_id=9544),
  
  # Fish & Amphibians
  GRCz11 = list(sp="danio_rerio", cap="Danio_rerio", org_db="org.Dr.eg.db", tax_id=7955),
  UCB_Xtro_10.0 = list(sp="xenopus_tropicalis", cap="Xenopus_tropicalis", org_db="org.Xt.eg.db", tax_id=8364), # no such org_db found
  Ssal_v3.1 = list(sp="salmo_salar", cap="Salmo_salar", org_db="org.Ssa.eg.db", tax_id=8030), # no such org_db found
  
  # Invertebrates
  BDGP6 = list(sp="drosophila_melanogaster", cap="Drosophila_melanogaster", org_db="org.Dm.eg.db", tax_id=7227),
  WBcel235 = list(sp="caenorhabditis_elegans", cap="Caenorhabditis_elegans", org_db="org.Ce.eg.db", tax_id=6239)
)

info <- genome_mapping[[assembly]]
if (is.null(info)) stop("No mapping for assembly: ", assembly)

# Load TxDb
txdb_file <- file.path(ref_dir, paste0("TxDb_", assembly, ".sqlite"))
if (!file.exists(txdb_file)) stop("TxDb not found. Run peak annotation first.")
txdb <- loadDb(txdb_file)

# OrgDb handling logic
org_pkg <- info$org_db
org_version <- "1.0"
org_tar <- file.path(ref_dir, paste0(org_pkg, "_", org_version, ".tar.gz"))

# Check if OrgDb is installed
if (requireNamespace(org_pkg, quietly=TRUE)) {
  message("Loading installed OrgDb: ", org_pkg)
  suppressPackageStartupMessages(library(org_pkg, character.only=TRUE))
} else {
  # Try installing OrgDb from Bioconductor
  if (requireNamespace("BiocManager", quietly=TRUE)) {
    tryCatch({
      message("Attempting to install OrgDb from Bioconductor: ", org_pkg)
      BiocManager::install(org_pkg, ask=FALSE, update=FALSE)
      suppressPackageStartupMessages(library(org_pkg, character.only=TRUE))
    }, error = function(e) {
      message("Bioconductor installation failed for OrgDb: ", org_pkg)
    })
  }
  
  # If installation from Bioconductor fails, attempt to install from tarball if available
  if (!requireNamespace(org_pkg, quietly=TRUE) && file.exists(org_tar)) {
    message("Installing cached OrgDb from: ", org_tar)
    install.packages(org_tar, repos=NULL, type="source")
    suppressPackageStartupMessages(library(org_pkg, character.only=TRUE))
  }
  
  # If still not available, try to build the OrgDb from GTF using AnnotationForge
  if (!requireNamespace(org_pkg, quietly=TRUE)) {
    message("Building OrgDb from GTF using AnnotationForge…")
    if (!requireNamespace("AnnotationForge", quietly=TRUE)) {
      BiocManager::install("AnnotationForge")
    }
    library(AnnotationForge)
    
    # Import GTF and extract gene entries
    gtf_file <- file.path(ref_dir, paste0(assembly, ".gtf"))
    if (!file.exists(gtf_file)) stop("GTF file not found at: ", gtf_file)
    
    gtf <- rtracklayer::import(gtf_file)
    genes <- gtf[gtf$type == "gene"]
    
    # Create a gene_info data.frame
    gene_info <- data.frame(
      GID = genes$gene_id,
      SYMBOL = genes$gene_name,
      BIOTYPE = genes$gene_biotype,
      stringsAsFactors = FALSE
    )
    
    # Remove duplicate GIDs and genes with missing SYMBOLs
    gene_info <- gene_info[!duplicated(gene_info$GID), ]
    gene_info <- gene_info[!is.na(gene_info$SYMBOL), ]
    
    # Build OrgDb package using AnnotationForge
    gs <- strsplit(info$sp, "_")[[1]]
    genus <- tools::toTitleCase((gs[1]))
    species <- gs[2]
    
    AnnotationForge::makeOrgPackage(
      gene_info = gene_info,          # Data.frame with GID as first column
      version = org_version,
      maintainer = "KI <KI@epigenbioinfo.org>",
      author = "KI",
      outputDir = ref_dir,
      tax_id = info$tax_id,
      genus = genus,
      species = species,
      goTable = NULL,               
      verbose = TRUE
    )
    
    # Try installing the newly built OrgDb from the tarball
    install.packages(org_tar, repos=NULL, type="source")
    suppressPackageStartupMessages(library(org_pkg, character.only=TRUE))
  }
  
  # If still no OrgDb available, proceed without it
  if (!requireNamespace(org_pkg, quietly=TRUE)) {
    message("Proceeding without OrgDb package.")
  }
}


# ... [previous code remains unchanged] ...

# Output directory
annot_dir <- file.path(dirname(noisq_xlsx), "Annotated_NOISeq")
dir.create(annot_dir, showWarnings=FALSE, recursive=TRUE)

# Load data with proper handling
de_df <- read.xlsx(noisq_xlsx)

# Read peak counts file as tab-delimited, skipping comment lines
lines <- readLines(peak_counts_file)
lines <- lines[!grepl("^#", lines)]  # remove comment lines
peak_df <- fread(paste(lines, collapse = "\n"), sep = "\t", header = TRUE)

# Check required columns
required_cols <- c("Geneid", "Chr", "Start", "End")
if (!all(required_cols %in% colnames(peak_df))) {
  stop("Missing required columns in peak_counts.txt. Needed: Geneid, Chr, Start, End")
}

# Rename Geneid to peak_id for proper merging
setnames(peak_df, "Geneid", "peak_id")

# Merge by peak_id (inner join)
merged <- merge(peak_df, de_df, by = "peak_id")

# Build GRanges from merged data
gr <- GRanges(
  seqnames = merged$Chr,
  ranges   = IRanges(start = merged$Start, end = merged$End)
)

# Annotate peaks
ann <- annotatePeak(gr, TxDb=txdb, annoDb=org_pkg, tssRegion=c(-3000, 3000))

# Clean annotation table
df_ann <- as.data.frame(ann)
df_ann[] <- lapply(df_ann, function(x) {
  if (is.character(x)) {
    x <- gsub("[\r\n\t]+", " ", x)   # Remove tabs/newlines
    x <- gsub(" +", " ", x)          # Collapse extra spaces
    trimws(x)                        # Trim leading/trailing
  } else x
})

# Save as TSV
write.table(df_ann,
            file = file.path(annot_dir, "noisq_annotated.tsv"),
            sep = "\t", quote = TRUE, row.names = FALSE)


# Generate annotation summary plots - CORRECTED SECTION
pdf(file.path(annot_dir, "NOISeq_annotation_summary.pdf"), width=8, height=10)
vennpie(ann)
plotAnnoPie(ann)

# New page for summary text
grid.newpage()
grid.text("Annotation Summary", x = 0.5, y = 0.95,
          gp = gpar(fontsize = 14, fontface = "bold"))

# Capture and print summary to PDF
summary_text <- capture.output(print(ann))
if (length(summary_text) > 0) {
  for (i in seq_along(summary_text)) {
    grid.text(summary_text[i],
              x    = 0.02,
              y    = 0.9 - (i - 1) * 0.025,
              just = "left",
              gp   = gpar(fontsize = 9, family = "mono"))
  }
}

# Also save summary as a text file
writeLines(summary_text, file.path(annot_dir, "noisq_annotation_summary.txt"))

dev.off()

message("[✓] Annotation complete: ", annot_dir)
