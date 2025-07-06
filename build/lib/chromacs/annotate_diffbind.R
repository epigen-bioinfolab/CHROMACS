#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ChIPseeker)
  library(clusterProfiler)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(grid)
  library(gridExtra)
  library(tools)
  library(R.utils)
  library(rtracklayer)
  if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript annotate_diffbind.R <diffbind_dir> <assembly> <ref_dir>")
}
diffbind_dir <- args[1]
assembly <- args[2]
ref_dir  <- args[3]

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

# Load existing TxDb
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
# output directory
annotated_dir <- file.path(diffbind_dir, "Annotated_DiffBind")
dir.create(annotated_dir, recursive=TRUE, showWarnings=FALSE)

# Helper to annotate and plot a GRanges object
annotate_gr <- function(gr, name) {
  # 1) for the annotation
  ann <- annotatePeak(
    gr,
    TxDb      = txdb,
    annoDb    = org_pkg,
    tssRegion = c(-3000, 3000)
  )
  
  # 2) clean and export annotation table
  df_ann <- as.data.frame(ann)
  df_ann[] <- lapply(df_ann, function(x) {
    if (is.character(x)) {
      x <- gsub("[\r\n\t]+", " ", x)   # Remove all newline/tab characters
      x <- gsub(" +", " ", x)          # Collapse multiple spaces
      trimws(x)                        # Trim leading/trailing spaces
    } else {
      x
    }
  })

  out_txt <- file.path(annotated_dir, paste0(name, "_annotated.tsv"))
  write.table(df_ann,
              file      = out_txt,
              sep       = "\t",
              quote     = TRUE,
              row.names = FALSE,
              col.names = TRUE)
  
  # 3) opens a combined PDF
  out_pdf <- file.path(annotated_dir,
                       paste0("AnnotVis_", name, "_combined.pdf"))
  pdf(out_pdf, width=8, height=10)
  
  # 3a) venn + pie
  vennpie(ann)
  
  # 3b) new page for summary text
  grid.newpage()
  grid.text("Annotation Summary",
            x = 0.5, y = 0.95,
            gp = gpar(fontsize = 14, fontface = "bold"))
  
  # capture and draw the summary of `ann`  
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
  
  # 3c) final annotation pie
  plotAnnoPie(ann)
  
  dev.off()
  
  message("[✓] Annotated & plotted: ", name)
}

# 1) Annotate diffbind_results.csv (handles both "Chr"/"Start"/"End" and "seqnames"/"start"/"end") # seqnames/start/end is our preferred one
csv_file <- file.path(diffbind_dir, "diffbind_results.csv")
if (file.exists(csv_file)) {
  df <- read.csv(csv_file, stringsAsFactors=FALSE)
  cols_lower <- tolower(colnames(df))

  # detect which naming scheme is present
  if (all(c("chr","start","end") %in% cols_lower) ||
      all(c("seqnames","start","end") %in% cols_lower)) {

    # pick the right column names
    if ("seqnames" %in% cols_lower) {
      chr_col   <- colnames(df)[ which(cols_lower == "seqnames")[1] ]
    } else {
      chr_col   <- colnames(df)[ which(cols_lower == "chr")[1] ]
    }
    start_col <- colnames(df)[ which(cols_lower == "start")[1] ]
    end_col   <- colnames(df)[ which(cols_lower == "end")[1] ]

    # build the GRanges
    gr_csv <- GRanges(
      seqnames = df[[chr_col]],
      ranges   = IRanges(
        start = df[[start_col]],
        end   = df[[end_col]]
      )
    )
    annotate_gr(gr_csv, "diffbind_results")

  } else {
    message("Skipping diffbind_results.csv: neither Chr/Start/End nor seqnames/start/end found")
  }

} else {
  message("diffbind_results.csv not found in ", diffbind_dir)
}

# 2) Annotate each BED file in the directory
bed_files <- list.files(diffbind_dir, pattern="\\.bed$", full.names=TRUE)
if (length(bed_files) == 0) {
  message("No .bed files found in ", diffbind_dir)
} else {
  for (bf in bed_files) {
    name <- file_path_sans_ext(basename(bf))
    gr   <- import(bf, format="BED")
    annotate_gr(gr, name)
  }
}

message("All requested annotations complete. Results are in:\n  ", annotated_dir)


