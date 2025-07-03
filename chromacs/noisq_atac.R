#!/usr/bin/env Rscript

library(NOISeq)
library(openxlsx)
library(data.table)
library(ggplot2)
library(ggrepel)

args <- commandArgs(TRUE)

if (length(args) != 4) {
  stop("Usage: Rscript noisq_atac.R <peak_counts.csv> <metadata.csv> <output.xlsx> <qval>")
}

counts_file  <- args[1]
metadata_file <- args[2]
output_file   <- args[3]
q_threshold   <- as.numeric(args[4])

if (is.na(q_threshold) || q_threshold < 0 || q_threshold > 1) {
  stop("❌ Invalid q-value threshold. Must be between 0 and 1.")
}

# ===================== Load Data =====================
counts_dt <- fread(counts_file)
metadata <- fread(metadata_file)
pid <- counts_dt$peak_id
sample_ids <- metadata$SampleID
conditions <- metadata$Condition
if (!all(sample_ids %in% colnames(counts_dt))) {
  stop("❌ SampleIDs in metadata do not match count matrix.")
}
counts_dt <- counts_dt[, sample_ids, with = FALSE]

# Subset and convert
counts <- as.matrix(counts_dt)
rownames(counts) <- pid

# ===================== Factors =====================
factors <- data.frame(condition = conditions)
rownames(factors) <- sample_ids

# ===================== NOISeq Object & Normalize =====================
mydata <- readData(data = counts, factors = factors)
myTMM <- tmm(exprs(mydata), long = 1000, lc = 0)
myLog2 <- log2(myTMM + 1)

# Filter out rows with all 0s
keep <- rowSums(counts) > 0
counts <- counts[keep, ]
myTMM <- myTMM[keep, ]
myLog2 <- myLog2[keep, ]

# Updated NOISeq object
mydata_filtered <- readData(data = counts, factors = factors)

# ===================== Run NOISeq =====================
myresults <- noiseq(
  mydata_filtered,
  factor = "condition",
  k = NULL,
  norm = "n",
  pnr = 0.2,
  nss = 5,
  v = 0.02,
  lc = 0,
  replicates = "no"
)

# Try to get DE genes with q >= 0.9
res_all <- myresults@results[[1]]
de_res <- tryCatch({
  degenes(myresults, q = 0.9, M = NULL)
}, error = function(e) {
  warning("⚠️ degenes() failed; using all features instead.")
  res_all
})

if (nrow(de_res) == 0) {
  warning("⚠️ No DE features detected (q >= 0.9); using all results.")
  de_res <- res_all
}

# ============ Final Output Assembly ================
colnames(myTMM) <- paste0("TMM:", colnames(myTMM))
colnames(myLog2) <- paste0("log2TMM:", colnames(myLog2))

filtered_rows <- rownames(de_res)
filtered_rows <- intersect(filtered_rows, rownames(counts))  # Prevent out-of-bounds error

out_df <- data.frame(
  peak_id = filtered_rows,
  counts[filtered_rows, , drop = FALSE],
  myTMM[filtered_rows, , drop = FALSE],
  myLog2[filtered_rows, , drop = FALSE],
  res_all[filtered_rows, , drop = FALSE]
)

# ============  Generate Plots ==============
message("Generating NOISeq diagnostic + volcano + PCA plots...")
plot_file <- file.path(dirname(output_file), "NOISeq_plots.pdf")

pdf(plot_file, width = 10, height = 8)

# --- Built-in NOISeq plots ---
DE.plot(myresults, q = q_threshold, graphic = "expr", log.scale = TRUE)
DE.plot(myresults, q = q_threshold, graphic = "MD")

# --- Volcano plot (safeguarded) ---
if (!all(c("M", "probability") %in% colnames(res_all))) {
  warning("⚠️ Missing 'M' or 'probability' in NOISeq results. Skipping volcano plot.")
} else {
  volcano_df <- data.frame(
    M = res_all$M,
    prob = res_all$probability,
    peak_id = rownames(res_all)
  )
  volcano_df$pval <- 1 - volcano_df$prob
  volcano_df$signif <- with(volcano_df, ifelse(pval < 0.05 & abs(M) > 1, "Significant", "NS"))

  print(
    ggplot(volcano_df, aes(x = M, y = -log10(pval), color = signif)) +
      geom_point(alpha = 0.7) +
      scale_color_manual(values = c("NS" = "grey60", "Significant" = "blue")) +
      labs(title = "NOISeq Volcano Plot", x = "log2 Fold Change (M)", y = "-log10(p-value)") +
      theme_minimal()
  )
}

# --- PCA plot ---
pca <- prcomp(t(myLog2), scale. = TRUE)
pca_df <- data.frame(pca$x, condition = factors$condition)
ggplot(pca_df, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  labs(title = "PCA of TMM-normalized log2 counts", x = "PC1", y = "PC2") +
  theme_minimal()

dev.off()
message(" Plots saved to: ", plot_file)

# ============ Final Output ============
write.xlsx(out_df, output_file, rowNames = FALSE)
cat(paste(" NOISeq results written to:", output_file, "\n"))


