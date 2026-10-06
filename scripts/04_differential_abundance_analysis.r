# =========================================================
# Project      : 16S rRNA Sequencing - Microbiome Analysis Pipeline
# Dataset      : PRJNA759642
# Samples      : Human gut microbiome, 16S rRNA V4 amplicon sequencing
# Script       : Differential Abundance Analysis
# Description  : Load the genus-level phyloseq object and sample metadata,
# filter low-abundance and low-prevalence families,
# normalize sequencing counts,
# identify differentially abundant families between clinical groups,
# and evaluate disease-associated taxonomic shifts.
# =========================================================


# =========================================================
# Libraries & Setup
# =========================================================

source("Setup_Environment/00_paths.R")
source("Setup_Environment/01_environment.R")
source("Setup_Environment/02_io_helpers.R")
source("Setup_Environment/03_seed.R")


# =========================================================
# Load and Prepare Data
# =========================================================

message("Load genus-level phyloseq object...")

ps_obj_unfiltered <- readRDS(file.path(paths$obj, "phyloseq_object_unfiltered.rds"))

ps_family <- tax_glom(
  ps_obj_unfiltered,
  taxrank = "Family",
  NArm = FALSE
)

message("[QC] Number of samples: ", nsamples(ps_family))
message("[QC] Number of families: ", (ntaxa(ps_family)))

meta_da <- data.frame(sample_data(ps_family))

message("[QC] Sample metadata loaded for ", nrow(meta_da), " samples.")

# ---------------------------------------------------------
# Validate clinical metadata
# ---------------------------------------------------------

message("[CHECK] Validating clinical group information...")

stopifnot("SampleType" %in% colnames(meta_da))

meta_da$SampleType <- factor(
  meta_da$SampleType,
  levels = c("Healthy", "UC")
)

stopifnot(!anyNA(meta_da$SampleType))

message("[QC] Sample groups:")
print(table(meta_da$SampleType))


# =========================================================
# Prepare Family-Level Count Matrix
# =========================================================

message("[ANALYSIS] Preparing family-level count matrix...")

counts_family <- as(
  otu_table(ps_family),
  "matrix"
)

if (taxa_are_rows(ps_family)) {
  counts_family <- counts_family
} else {
  counts_family <- t(counts_family)
}

stopifnot(identical(colnames(counts_family), sample_names(ps_family)))

message(
  "[QC] Count matrix dimensions: ",
  nrow(counts_family),
  " families x ",
  ncol(counts_family),
  " samples."
)


# =========================================================
# Filter Low-Information Families
# =========================================================

message("[ANALYSIS] Filtering low-information families...")

# ---------------------------------------------------------
# Remove genera with >= 90% zero counts
# ---------------------------------------------------------

zero_fraction <- rowMeans(counts_family == 0)

keep_zero <- zero_fraction < 0.90

message(
  "[QC] Families retained after zero-count filtering: ",
  sum(keep_zero),
  " / ",
  length(keep_zero)
)

counts_filtered <- counts_family[keep_zero, , drop = FALSE]


# ---------------------------------------------------------
# Remove families with very low variance
# ---------------------------------------------------------

feature_variance <- apply(
  counts_family,
  1,
  var
)

variance_threshold <- 0.5 * median(feature_variance)

keep_variance <- feature_variance >= variance_threshold

keep_combined <- keep_zero & keep_variance

message(
  "[QC] Families retained after variance filtering: ",
  sum(keep_combined),
  " / ",
  length(keep_combined),
  " (variance threshold = ",
  signif(variance_threshold, 4),
  ")"
)

counts_filtered <- counts_family[keep_combined, , drop = FALSE]

message(
  "[QC] Final count matrix: ",
  nrow(counts_filtered),
  " families x ",
  ncol(counts_filtered),
  " samples."
)


# =========================================================
# Create edgeR DGEList
# =========================================================

message("[ANALYSIS] Creating edgeR DGEList object...")

dge <- DGEList(
  counts = counts_filtered,
  group = meta_da$SampleType
)

group <- factor(
  as.character(sample_data(ps_family)$SampleType[
    match(colnames(dge), rownames(sample_data(ps_family)))
  ]),
  levels = c("Healthy", "UC")
)

stopifnot(identical(colnames(dge$counts), rownames(meta_da)))
stopifnot(!anyNA(group))
stopifnot(length(group) == ncol(dge))

message(
  "[QC] edgeR library sizes - min: ",
  min(dge$samples$lib.size),
  ", median: ",
  median(dge$samples$lib.size),
  ", max: ",
  max(dge$samples$lib.size)
)

# ------------------------------------------------
#         Quality Control (QC)
# ------------------------------------------------ 

# Plot library sizes and samples before normalization 

message("[QC] Generating library size plot before TMM normalization...")

open_png(filename = "library_size_before_TMM.png", dir = paths$plot_diff_abu,
         width = 1800, height = 1000)

barplot(
  dge$samples$lib.size * 1e-6,
  names.arg = colnames(dge),
  las = 2,
  cex.names = 0.7,
  ylab = "Library size (millions)",
  main = "Library size before TMM normalization"
)

close_png()


# =========================================================
# TMM Normalization
# =========================================================

# ---------------------------------------------------------
# Boxplot before normalization
# ---------------------------------------------------------

open_png(filename = "log2_CPM_before_TMM.png", dir = paths$plot_diff_abu,
         width = 1800, height = 1000)

boxplot(
  cpm(dge, log = TRUE, prior.count = 1), 
  las = 2,
  xaxt = "n",
  main = "Log2 CPM before TMM normalization",
  ylab = "log2(CPM)"
)

axis(
  1,
  at = seq_along(colnames(dge)),
  labels = colnames(dge),
  las = 2,
  cex.axis = 0.7
)

close_png()

message("[ANALYSIS] Calculating TMM normalization factors...")

dge <- calcNormFactors(
  dge,
  method = "TMM"
)

stopifnot(all(is.finite(dge$samples$norm.factors)))

message("[QC] TMM normalization factors:")
print(dge$samples[, c("lib.size", "norm.factors")])

# ---------------------------------------------------------
# Boxplot afetr normalization
# ---------------------------------------------------------

open_png(filename = "log2_CPM_after_TMM.png", dir = paths$plot_diff_abu,
         width = 1800, height = 1000)

boxplot(
  cpm(dge, log = TRUE, prior.count = 1), 
  las = 2,
  xaxt = "n",
  main = "Log2 CPM after TMM normalization",
  ylab = "log2(CPM)"
)

axis(
  1,
  at = seq_along(colnames(dge)),
  labels = colnames(dge),
  las = 2,
  cex.axis = 0.7
)

close_png()

# ------------------------------------------------
#         Exploratory Data Analysis (EDA)
# ------------------------------------------------ 

# ---------------------------------------------------------
# MDS plot
# ---------------------------------------------------------

message("[ANALYSIS] Generating MDS plot after TMM normalization...")

group_colors <- c(
  "Healthy" = "#0072B2",
  "UC" = "#D55E00"
)

open_png(filename = "MDS_plot_TMM.png", dir = paths$plot_diff_abu,
         width = 1800, height = 1000)

plotMDS(
  dge,
  col = group_colors[as.character(group)],
  pch = 16,
  main = "MDS plot after TMM normalization"
)

legend(
  "topright",
  legend = names(group_colors),
  col = group_colors,
  pch = 16,
  bty = "n"
)

close_png()

# ---------------------------------------------------------
# PCA
# ---------------------------------------------------------

logCPM <- cpm(dge, log = TRUE, prior.count = 1)
pca <- prcomp(t(logCPM), scale. = TRUE)
pca_var <- pca$sdev^2
pca_var_percent <- round(100 * pca_var / sum(pca_var), 1)

metadata <- data.frame(
  Sample = colnames(logCPM),
  Group = group
)

pca_df <- data.frame(pca$x[, 1:2], metadata)
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.8, size = 3) +
  xlab(paste0("PC1 (", pca_var_percent[1], "%)")) +
  ylab(paste0("PC2 (", pca_var_percent[2], "%)")) +
  theme_bw() +
  scale_color_manual(values = group_colors) +
  ggtitle("PCA of normalized log2 CPM values")

save_plot(pca_plot, filename = "PCA_normalized_log2CPM.png", dir = paths$plot_diff_abu)

# ---------------------------------------------------------
# Sample distance heatmap
# ---------------------------------------------------------

open_png(filename = "sample_distance_heatmap.png", dir = paths$plot_diff_abu,
         width = 1800, height = 1000)

distance <- dist(t(logCPM))
pheatmap(as.matrix(distance),
         main = "Heatmap of sample-to-sample distances",
         fontsize_row = 10,
         fontsize_col = 10)

close_png()


# =========================================================
# Estimate Negative Binomial Dispersion
# =========================================================

message("[ANALYSIS] Estimating negative binomial dispersion...")

design <- model.matrix(
  ~ group,
  data = dge$samples
)

dge <- estimateDisp(dge, design)

message(
  "[QC] Common dispersion estimate: ",
  signif(dge$common.dispersion, 4)
)

open_png(filename = "BCV_and_meanVariance_plots.png", dir = paths$plot_diff_abu,
  width = 1200, height = 600)

par(
  mfrow = c(1, 2),
  mar = c(5, 4, 4, 2)
)

plotBCV(dge)
plotMeanVar(dge)

par(mfrow = c(1, 1))

close_png()


# =========================================================
# Differential Abundance Testing
# =========================================================

message("[ANALYSIS] Performing edgeR exact test: UC vs Healthy...")

exact_test <- exactTest(
  dge,
  pair = c("Healthy", "UC")
)

message("[QC] Differential abundance testing completed.")


# =========================================================
# Extract Differential Abundance Results
# =========================================================

message("[ANALYSIS] Extracting differential abundance results...")

da_results <- topTags(
  exact_test,
  n = Inf,
  sort.by = "PValue"
)$table

da_results$FeatureID <- rownames(da_results)

# ---------------------------------------------------------
# Add taxonomic annotations to differential abundance results
# ---------------------------------------------------------

message("[ANALYSIS] Adding taxonomic annotations to differential abundance results...")

tax_family <- as.data.frame(tax_table(ps_family)[, "Family"])
tax_family$FeatureID <- rownames(tax_family)

da_results <- merge(
  da_results,
  tax_family,
  by = "FeatureID",
  all.x = TRUE,
  sort = FALSE
)

message(
  "[CHECK] Taxonomic annotations added for ",
  sum(!is.na(da_results$Family)),
  " of ",
  nrow(da_results),
  " differential abundance features."
)

da_results$Significant <- da_results$FDR < 0.25

# ---------------------------------------------------------
# Restore statistical ranking
# ---------------------------------------------------------

da_results <- da_results[
  order(da_results$FDR),
  ,
  drop = FALSE
]

message(
  "[STAT] Families with BH-FDR < 0.25: ",
  sum(da_results$Significant)
)


# =========================================================
# Differential abudance
# =========================================================

# ---------------------------------------------------------
# MD plot
# ---------------------------------------------------------

open_png(filename  = "MD_plot_exactTest.png", dir = paths$plot_diff_abu,
         width = 1000, height = 800)

plotMD(exact_test, main = "MD plot of exact-test Data")
abline(h = c(-1, 1), col = "blue", lty = 2)

close_png()

# ---------------------------------------------------------
# Volcano plot
# ---------------------------------------------------------

message("[ANALYSIS] Generating differential abundance volcano plot...")

da_results$Category <- "Not significant"

da_results$Category[
  da_results$FDR < 0.25 & da_results$logFC > 0
] <- "Enriched in UC"

da_results$Category[
  da_results$FDR < 0.25 & da_results$logFC < 0
] <- "Depleted in UC"

volcano_plot <- ggplot(
  da_results,
  aes(
    x = logFC,
    y = -log10(FDR),
    color = Category
  )
) +
  geom_point(
    alpha = 0.7,
    size = 2
  ) +
  geom_hline(
    yintercept = -log10(0.25),
    linetype = "dashed"
  ) +
  scale_color_manual(
    values = c(
      "Enriched in UC" = "#D55E00",
      "Depleted in UC" = "#0072B2",
      "Not significant" = "grey70"
    )
  ) +
  labs(
    title = "Differential abundance of families",
    x = "Log2 fold change",
    y = "-Log10(BH-FDR)",
    color = "Association"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(face = "bold")
  )

save_plot(volcano_plot, filename = "differential_abundance_volcano.png", dir = paths$plot_diff_abu)

# ---------------------------------------------------------
# Heatmap of differentially abundant families
# ---------------------------------------------------------

message("[ANALYSIS] Generating heatmap of differentially abundant families (ComplexHeatmap)...")

# ---------------------------------------------------------
# Select and rank significantly differentially abundant families
# ---------------------------------------------------------

da_significant <- da_results[da_results$Significant, , drop = FALSE]

# Order families by increasing BH-FDR,

da_significant <- da_significant[order(da_significant$FDR), , drop = FALSE]

da_features <- da_significant$FeatureID

# ---------------------------------------------------------
# Extract normalized abundance values for significant families
# ---------------------------------------------------------

# Extract log2-CPM values for the differentially abundant families
logCPM_da <- logCPM[da_features, , drop = FALSE]

# Standardize abundance values within each families to obtain row-wise z-scores

logCPM_da_scaled <- t(scale(t(logCPM_da)))

# ---------------------------------------------------------
# Assign readable family labels
# ---------------------------------------------------------

tax_sig <- as.data.frame(
  tax_table(ps_family)[da_features, , drop = FALSE],
  stringsAsFactors = FALSE
)

rank_cols <- rank_names(ps_family)
rank_cols <- rank_cols[seq_len(match("Family", rank_cols))]

label_feature <- function(row) {
  fam <- row[["Family"]]
  if (!is.na(fam) && nzchar(fam)) return(fam)
  for (r in rev(setdiff(rank_cols, "Family"))) {
    val <- row[[r]]
    if (!is.na(val) && nzchar(val)) {
      return(paste0("Unclassified (", tolower(r), " ", val, ")"))
    }
  }
  "Unclassified (no assignment)"
}

family_labels <- vapply(
  seq_len(nrow(tax_sig)),
  function(i) label_feature(tax_sig[i, , drop = FALSE]),
  character(1)
)

message("[CHECK] Unclassified features among significant families: ",
        sum(grepl("^Unclassified", family_labels)), " of ", length(family_labels))

# Direction of change (used to split rows and colour the FDR bars)
direction <- factor(
  ifelse(da_significant$logFC > 0, "Enriched in UC", "Depleted in UC"),
  levels = c("Enriched in UC", "Depleted in UC")
)

# Ensure unique row names when multiple features belong to the same family
rownames(logCPM_da_scaled) <- make.unique(family_labels)

if (nrow(logCPM_da_scaled) > 0) {
  
  # -------------------------------------------------------
  # Define heatmap color scale
  # -------------------------------------------------------
  
  col_fun <- colorRamp2(
    c(min(logCPM_da_scaled), 0, max(logCPM_da_scaled)),
    c("#F4C2C2", "white", "#5B2C6F")
  )
  
  # -------------------------------------------------------
  # Add clinical group annotation
  # -------------------------------------------------------
  
  col_annotation <- HeatmapAnnotation(
    SampleType = group,
    col = list(
      SampleType = c(
        "Healthy" = "#0072B2",
        "UC" = "#D55E00"
      )
    )
  )
  
  # -------------------------------------------------------
  # Add FDR annotation
  # -------------------------------------------------------
  
  # Display -log10(BH-FDR) for each family as a bar, coloured by direction.
  
  fdr_annotation <- rowAnnotation(
    "-log10(FDR)" = anno_barplot(
      -log10(da_significant$FDR),
      gp = gpar(
        fill = ifelse(direction == "Enriched in UC", "#D55E00", "#0072B2"),
        col = NA
      ),
      bar_width = 0.8,
      border = FALSE,
      width = unit(3, "cm"),
      axis_param = list(side = "bottom", gp = gpar(fontsize = 9))
    ),
    annotation_name_gp = gpar(fontsize = 10)
  )
  
  # -------------------------------------------------------
  # Construct heatmap
  # -------------------------------------------------------
  
  ht <- Heatmap(
    logCPM_da_scaled,
    name = "z-score\n(log2 CPM)",
    col = col_fun,
    top_annotation = col_annotation,
    left_annotation = fdr_annotation,
    
    # Split columns according to clinical group
    
    column_split = group,
    cluster_column_slices = FALSE,
    
    # Split rows by direction of change; keep the FDR ordering within each block
    
    row_split = direction,
    cluster_row_slices = FALSE,
    row_title = "%s",
    row_title_rot = 0,
    row_title_gp = gpar(fontsize = 12, fontface = "bold"),
    
    # Display family names and hide individual sample names
    
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 11),
    row_names_max_width = unit(9, "cm"),
    show_column_names = FALSE,
    
    # Preserve the statistical ordering of families
    
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    
    column_title = "Differentially abundant families (BH-FDR < 0.25)",
    column_title_gp = gpar(fontsize = 13, fontface = "bold"),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 10),
      labels_gp = gpar(fontsize = 9)
    )
  )
  
  # -------------------------------------------------------
  # Save heatmap
  # -------------------------------------------------------
  
  open_png(filename = "differentially_abundant_families_heatmap.png", dir = paths$plot_diff_abu,
           width = 2600, height = 800 + 80 * nrow(logCPM_da_scaled))
  
  # Combine the FDR dot plot and the abundance heatmap
  
  draw(ht)
  
  close_png()
  
} else {
  
  # Report when no family passes the predefined BH-FDR threshold
  message(
    "[QC] No families passed the BH-FDR < 0.25 threshold; ",
    "heatmap was not generated."
  )
}


# =========================================================
# Save Analysis Objects
# =========================================================

message("[OUTPUT] Saving Differential Abudance analysis objects...")

save_rds(dge, filename = "edgeR_dgelist_family.rds", dir = paths$obj)

save_rds(exact_test, filename = "edgeR_exact_test_family.rds", dir = paths$obj)

save_rds(da_results, filename = "differential_abundance_results_family.rds", dir = paths$obj)

save_csv(da_results, filename = "differential_abundance_results_family.csv", dir = paths$obj)

message("[OUTPUT] Differential abundance analysis objects saved successfully.")


# =========================================================
# Save session info
# =========================================================

message("[OUTPUT] Saving session information...")

save_session_info(
  filename = "sessionInfo_16S_PRJNA759642_differential_abundance.txt",
  dir = paths$logs,
  label = "16S rRNA sequencing - Differential abundance analysis (PRJNA759642)"
)

message("[OUTPUT] Session information saved to: ", paths$logs)

# =========================================================
# Final pipeline message
# =========================================================

message("=================================================")
message("[PIPELINE] Differential abundance analysis completed successfully for PRJNA759642.")
message("=================================================")


