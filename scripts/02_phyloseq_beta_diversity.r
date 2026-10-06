# =========================================================
# Project      : 16S rRNA Sequencing - Microbiome Analysis Pipeline
# Dataset      : PRJNA759642
# Samples      : Human gut microbiome, 16S rRNA V4 amplicon sequencing
# Script       : Phyloseq Construction, Taxonomic Composition and Beta-Diversity Analysis
# Description  : Load processed DADA2 outputs and sample metadata,
#                validate sample and ASV identifiers,
#                construct the phyloseq object,
#               assess sequencing depth and phylum-level composition,
#               perform Bray-Curtis beta-diversity analysis,
#               visualize samples using PCoA,
#               and test group differences using PERMANOVA and beta-dispersion.
# =========================================================


# =========================================================
#                  Libraries & Setup
# =========================================================

source("Setup_Environment/00_paths.R")
source("Setup_Environment/01_environment.R")
source("Setup_Environment/02_io_helpers.R")
source("Setup_Environment/03_seed.R")


# =========================================================
# Load processed DADA2 outputs and sample metadata
# =========================================================

message("[LOAD] Loading DADA2 outputs and sample metadata...")

taxa <- readRDS(file.path(paths$obj, "taxa.rds"))

seqtab.nochim <- readRDS(file.path(paths$obj, "seqtab_nochim.rds"))

meta_raw <- read_excel(file.path(paths$metadata, "meta_clean.xlsx"), sheet = 1)  

meta <- as.data.frame(meta_raw)

message(
  "[LOAD] Loaded ",
  nrow(seqtab.nochim),
  " samples and ",
  ncol(seqtab.nochim),
  " non-chimeric ASVs."
)

message(
  "[LOAD] Metadata contains ",
  nrow(meta),
  " samples and ",
  ncol(meta),
  " variables."
)


# =========================================================
# Validate sample and ASV identifiers
# =========================================================

message("[CHECK] Validating sample identifiers...")

# Check required sample identifier
stopifnot(all(rownames(seqtab.nochim) %in% meta[, "Run"]))

# Check uniqueness
stopifnot(!anyDuplicated(rownames(seqtab.nochim)))
stopifnot(!anyDuplicated(meta$Run))

# Check completeness of sample matching
stopifnot(all(rownames(seqtab.nochim) %in% meta$Run))
stopifnot(all(meta$Run %in% rownames(seqtab.nochim)))

# Reorder metadata to exactly match the ASV table

meta <- meta[match(rownames(seqtab.nochim), meta[, "Run"]), ]

rownames(meta) <- meta[, "Run"] 

# Convert to data.frame explicitly
meta <- as.data.frame(meta)

# Final exact matching check
stopifnot(identical(rownames(seqtab.nochim), rownames(meta)))

message(
  "[CHECK] Sample identifiers validated: ",
  nrow(meta),
  " samples matched exactly."
)

# Validate ASV identifiers

message("[CHECK] Validating ASV identifiers...")

stopifnot(identical(
  colnames(seqtab.nochim),
  rownames(taxa)
))

message(
  "[CHECK] ASV identifiers validated: ",
  ncol(seqtab.nochim),
  " ASVs matched between sequence table and taxonomy table."
)

# =========================================================
# Construct phyloseq object
# =========================================================

message("[ANALYSIS] Constructing phyloseq object...")

otu <- otu_table(seqtab.nochim, taxa_are_rows = FALSE)

tax <- tax_table(taxa)

samples <- sample_data(meta)

ps_obj <- phyloseq::phyloseq(otu, tax, samples)

message("[CHECK] Validating phyloseq object...")

stopifnot(nsamples(ps_obj) == nrow(seqtab.nochim))

stopifnot(ntaxa(ps_obj) == ncol(seqtab.nochim))

stopifnot(identical(sample_names(ps_obj), rownames(meta)))

message(
  "[ANALYSIS] Phyloseq object constructed successfully: ",
  nsamples(ps_obj),
  " samples × ",
  ntaxa(ps_obj),
  " ASVs."
)

ps_obj_unfiltered <- ps_obj

# =========================================================
# Prevalence exploration 
# =========================================================

message("[ANALYSIS] Computing ASV prevalence across samples...")

prevdf <- apply(otu_table(ps_obj), 
                MARGIN = ifelse(taxa_are_rows(ps_obj), 1, 2),
                FUN = function(x) sum(x > 0))

prevdf <- data.frame(
  Prevalence = prevdf,
  TotalAbundance = taxa_sums(ps_obj),
  tax_table(ps_obj)
)

message(
  "[QC] Median ASV prevalence: ",
  median(prevdf$Prevalence),
  " samples out of ",
  nsamples(ps_obj),
  "."
)

# ---------------------------------------------------------
# Prevalence vs abundance plot
# ---------------------------------------------------------

message("[ANALYSIS] Generating ASV prevalence-abundance plot...")

p_prevalence <- ggplot(
  prevdf,
  aes(x = TotalAbundance, y = Prevalence / nsamples(ps_obj), color = Phylum)
) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey50") +
  geom_point(size = 2, alpha = 0.6) +
  scale_x_log10() +
  labs(
    title = "ASV prevalence vs total abundance",
    x = "Total abundance (log10)",
    y = "Prevalence (fraction of samples)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  ) +
  facet_wrap(~ Phylum, scales = "free_y")

save_plot(p_prevalence, filename = "asv_prevalence_vs_abundance.png", 
          dir = paths$plot_phy, width = 12, height = 8)



# =========================================================
# Sequencing depth assessment
# =========================================================

message("[ANALYSIS] Assessing sequencing depth across samples...")

sample_sum_df <- data.frame(
  sample = names(sample_sums(ps_obj)),
  sum = as.numeric(sample_sums(ps_obj))
)


smin <- min(sample_sums(ps_obj))
smean <- mean(sample_sums(ps_obj))
smax <- max(sample_sums(ps_obj))

message(
  "[QC] Sequencing depth - minimum: ",
  format(smin, big.mark = ","),
  " reads; mean: ",
  format(round(smean), big.mark = ","),
  " reads; maximum: ",
  format(smax, big.mark = ","),
  " reads."
)

# ---------------------------------------------------------
# Sequencing depth distribution plot
# ---------------------------------------------------------

sum_reads <- ggplot(sample_sum_df, aes(x = sum)) +
  geom_histogram(
    bins = 15,
    color = "black",
    fill = "grey70"
  ) +
  labs(
    title = "Distribution of sequencing depth",
    x = "Sequencing depth (reads)",
    y = "Number of samples"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

save_plot(sum_reads, filename = "sequencing_depth_distribution.png", dir = paths$plot_phy)


# =========================================================
# Phylum-level taxonomic composition
# =========================================================

message("[ANALYSIS] Preparing Phylum-level taxonomic composition...")

# ---------------------------------------------------------
# Aggregate ASVs at Phylum level
# ---------------------------------------------------------

message("[ANALYSIS] Aggregating ASVs at Phylum level...")

ps_phylum <- tax_glom(
  ps_obj,
  taxrank = "Phylum",
  NArm = FALSE
)

# ---------------------------------------------------------
# Convert counts to relative abundance
# ---------------------------------------------------------

message("[ANALYSIS] Converting phylum counts to relative abundance...")

ps_phylum_rel <- transform_sample_counts(
  ps_phylum,
  function(x) x / sum(x)
)

# ---------------------------------------------------------
# Convert phyloseq object to plotting data frame
# ---------------------------------------------------------

phylum_df <- psmelt(ps_phylum_rel)

# ---------------------------------------------------------
# Handle missing taxonomy
# ---------------------------------------------------------

phylum_df$Phylum <- as.character(phylum_df$Phylum)

phylum_df$Phylum[
  is.na(phylum_df$Phylum) | phylum_df$Phylum == ""
] <- "Unclassified"

# ---------------------------------------------------------
# Identify dominant Phyla
# ---------------------------------------------------------

phylum_abundance <- phylum_df |>
  group_by(Phylum) |>
  summarise(
    mean_abundance = mean(Abundance),
    .groups = "drop"
  )

dominant_phyla <- phylum_abundance |>
  filter(mean_abundance >= 0.02) |>
  pull(Phylum)

message(
  "[QC] Phyla with mean relative abundance >= 2%: ",
  length(dominant_phyla)
)

# ---------------------------------------------------------
# Group low-abundance phyla as "Other"
# ---------------------------------------------------------

phylum_df <- phylum_df |>
  mutate(
    Phylum_plot = if_else(
      Phylum %in% dominant_phyla,
      Phylum,
      "Other"
    )
  )

# ---------------------------------------------------------
# Order samples by clinical group
# ---------------------------------------------------------

phylum_df$SampleType <- factor(
  phylum_df$SampleType,
  levels = c("Healthy", "UC")
)

message("[QC] Sample distribution by clinical group:")

print(table(phylum_df$SampleType))

sample_order <- phylum_df |>
  distinct(Sample, SampleType) |>
  arrange(SampleType, Sample) |>
  pull(Sample)

phylum_df$Sample <- factor(
  phylum_df$Sample,
  levels = sample_order
)

# ---------------------------------------------------------
# Order phyla by mean abundance
# ---------------------------------------------------------

phylum_order <- phylum_df |>
  group_by(Phylum_plot) |>
  summarise(
    mean_abundance = mean(Abundance),
    .groups = "drop"
  ) |>
  arrange(mean_abundance) |>
  pull(Phylum_plot)

phylum_df$Phylum_plot <- factor(
  phylum_df$Phylum_plot,
  levels = phylum_order
)

# ---------------------------------------------------------
# Generate phylum-level composition plot
# ---------------------------------------------------------

p_phylum <- ggplot(
  phylum_df,
  aes(
    x = Sample,
    y = Abundance,
    fill = Phylum_plot
  )
) +
  geom_col(width = 0.9) +
  facet_grid(
    ~ SampleType,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  labs(
    title = "Phylum-level microbial composition",
    x = NULL,
    y = "Relative abundance",
    fill = "Phylum"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 7
    ),
    axis.title.y = element_text(
      face = "bold"
    ),
    legend.title = element_text(
      face = "bold"
    ),
    legend.text = element_text(
      size = 9
    )
  )

save_plot(p_phylum, filename = "phylum_composition.png", dir = paths$plot_phy,
          width = 12, height = 7)

message("[OUTPUT] Phylum-level composition plot saved.")


# =========================================================
# Beta diversity - Genus-level Bray-Curtis and PCoA
# =========================================================

message("[ANALYSIS] Performing genus-level beta-diversity analysis...")

# ---------------------------------------------------------
# Aggregate ASVs at genus level before relative-abundance transformation
# ---------------------------------------------------------

message("[ANALYSIS] Aggregating ASVs at genus level...")

ps_genus <- tax_glom(
  ps_obj_unfiltered,
  taxrank = "Genus",
  NArm = FALSE
)

message(
  "[QC] Number of genera after taxonomic aggregation: ",
  ntaxa(ps_genus)
)

# ---------------------------------------------------------
# Convert genus-level counts to relative abundance
# ---------------------------------------------------------

message("[ANALYSIS] Converting genus-level counts to relative abundance...")

ps_genus_rel <- transform_sample_counts(
  ps_genus,
  function(x) x / sum(x)
)

# ---------------------------------------------------------
# Calculate Bray-Curtis distance
# ---------------------------------------------------------

message("[ANALYSIS] Calculating genus-level Bray-Curtis distance...")

bray_dist_genus <- phyloseq::distance(
  ps_genus_rel,
  method = "bray"
)

# ---------------------------------------------------------
# Principal Coordinates Analysis
# ---------------------------------------------------------

message("[ANALYSIS] Performing genus-level Bray-Curtis PCoA...")

pcoa_bray_genus <- ordinate(
  ps_genus_rel,
  method = "PCoA",
  distance = bray_dist_genus
)

pcoa_variance_genus <- pcoa_bray_genus$values$Relative_eig[1:2]

message(
  "[QC] PCoA variance explained - Axis 1: ",
  round(pcoa_variance_genus[1] * 100, 2),
  "%; Axis 2: ",
  round(pcoa_variance_genus[2] * 100, 2),
  "%."
)

# ---------------------------------------------------------
# Prepare sample metadata for PERMANOVA
# ---------------------------------------------------------

message("[ANALYSIS] Preparing metadata for genus-level PERMANOVA...")

sample_df_genus <- data.frame(
  sample_data(ps_genus_rel)
)

sample_df_genus$SampleType <- factor(
  sample_df_genus$SampleType,
  levels = c("Healthy", "UC")
)

sample_df_genus$Gender <- factor(
  sample_df_genus$Gender
)

sample_df_genus$Age <- as.numeric(
  sample_df_genus$Age
)

# ---------------------------------------------------------
# PERMANOVA
# ---------------------------------------------------------

message("[ANALYSIS] Running genus-level PERMANOVA with 999 permutations...")

set.seed(1234)

permanova_genus <- vegan::adonis2(
  bray_dist_genus ~ SampleType + Age + Gender,
  data = sample_df_genus,
  permutations = 999,
  by = "margin"
)

print(permanova_genus)

# ---------------------------------------------------------
# PCoA plot
# ---------------------------------------------------------

message("[ANALYSIS] Generating genus-level Bray-Curtis PCoA plot...")

p_pcoa_genus <- plot_ordination(
  ps_genus_rel,
  pcoa_bray_genus,
  color = "SampleType"
) +
  geom_point(
    size = 4,
    alpha = 0.8
  ) +
  labs(
    title = "PCoA of Bray-Curtis distances (genus level)",
    color = "Sample type"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )

save_plot(p_pcoa_genus, filename = "pcoa_bray_curtis_genus.png", dir = paths$plot_beta)

message("[OUTPUT] Genus-level Bray-Curtis PCoA plot saved.")

# ---------------------------------------------------------
# Homogeneity of multivariate dispersions
# ---------------------------------------------------------

message("[ANALYSIS] Testing homogeneity of multivariate dispersions...")

dispersion_genus <- vegan::betadisper(
  bray_dist_genus,
  sample_df_genus$SampleType
)

set.seed(1234)

dispersion_test_genus <- vegan::permutest(
  dispersion_genus,
  permutations = 999
)

print(dispersion_test_genus)

png(file.path(paths$plot_beta, "betadisper_boxplot.png"), width = 1600, height = 1200, res = 200)
boxplot(dispersion_genus,
        main = "Distance to group centroid",
        xlab = "Sample type",
        ylab = "Distance to centroid")
dev.off()


# =========================================================
# Multiple Regression on distance Matrices (MRM)
# =========================================================

message("[ANALYSIS] Preparing pairwise distance matrices for MRM...")

# ---------------------------------------------------------
# Extract sample labels and align metadata to Bray-Curtis distance matrix
# ---------------------------------------------------------

labels <- attr(bray_dist_genus, "Labels")
n <- length(labels)

status <- setNames(as.character(sample_df_genus$SampleType), rownames(sample_df_genus))[labels]
gender <- setNames(as.character(sample_df_genus$Gender), rownames(sample_df_genus))[labels]
age <- setNames(as.numeric(sample_df_genus$Age), rownames(sample_df_genus))[labels]

# ---------------------------------------------------------
# Construct pairwise dummy variables following the model specification
# ---------------------------------------------------------

pair_idx <- combn(n, 2)

d1 <- d2 <- g1 <- g2 <- age_diff <- numeric(ncol(pair_idx))

for (k in seq_len(ncol(pair_idx))) {
  i <- pair_idx[1, k]
  j <- pair_idx[2, k]
  
  si <- status[i]
  sj <- status[j]
  
  d1[k] <- as.numeric(
    (si == "Healthy" && sj == "UC") ||
      (si == "UC" && sj == "Healthy")
  )
  
  d2[k] <- as.numeric(
    si == "UC" && sj == "UC"
  )
  
  # Healthy-Healthy is the reference category: d1 = 0, d2 = 0
  
  gi <- gender[i]
  gj <- gender[j]
  
  g1[k] <- as.numeric(
    (gi == "m" && gj == "f") ||
      (gi == "f" && gj == "m")
  )
  
  g2[k] <- as.numeric(
    gi == "f" && gj == "f"
  )
  
  # Male-Male is the reference category: g1 = 0, g2 = 0
  
  age_diff[k] <- abs(age[i] - age[j])
}

# ---------------------------------------------------------
# Validate pairwise matrix alignment
# ---------------------------------------------------------

message("[CHECK] Validating pairwise distance matrix alignment...")

stopifnot(all(labels == rownames(sample_df_genus)[match(labels, rownames(sample_df_genus))]))
stopifnot(!anyNA(status), !anyNA(gender), !anyNA(age))
stopifnot(length(d1) == length(bray_dist_genus))

# ---------------------------------------------------------
# Convert predictor vectors to distance objects
# ---------------------------------------------------------

make_dist <- function(x, template) {
  stopifnot(length(x) == length(template))
  d <- x
  attributes(d) <- attributes(template)
  d
}

d1_dist <- make_dist(d1, bray_dist_genus)
d2_dist <- make_dist(d2, bray_dist_genus)
g1_dist <- make_dist(g1, bray_dist_genus)
g2_dist <- make_dist(g2, bray_dist_genus)
age_dist <- make_dist(age_diff, bray_dist_genus)

message("[CHECK] Pairwise predictor distance matrices constructed successfully.")

# ---------------------------------------------------------
# Multiple Regression on distance Matrices
# ---------------------------------------------------------

message("[ANALYSIS] Running MRM with 999 permutations...")

set.seed(1234)

mrm_genus <- ecodist::MRM(
  bray_dist_genus ~ d1_dist + d2_dist + age_dist + g1_dist + g2_dist,
  nperm = 999
)

print(mrm_genus)

message("[OUTPUT] Genus-level MRM analysis completed successfully.")


# =========================================================
# Save analysis objects
# =========================================================

message("Saving  results")

save_rds(ps_obj_unfiltered, filename = "phyloseq_object_unfiltered.rds", dir = paths$obj)
save_rds(bray_dist_genus, filename = "bray_curtis_distance_genus.rds", dir = paths$obj)
save_rds(permanova_genus, filename = "permanova_bray_curtis_genus.rds", dir = paths$obj)
save_rds(dispersion_test_genus, filename = "betadisper_test_bray_curtis_genus.rds", dir = paths$obj)
save_rds(pcoa_bray_genus, filename = "pcoa_bray_curtis_genus.rds", dir = paths$obj)
save_rds(prevdf, filename = "prevalence_table.rds", dir = paths$obj)
save_rds(ps_genus_rel, filename = "phyloseq_genus_relative.rds", dir = paths$obj)
save_rds(mrm_genus, filename = "mrm_genus.rds", dir = paths$obj)

message("Objects saved successfully.")


# =========================================================
# Phyloseq and beta-diversity analysis completed
# =========================================================

message("==============================================")
message("Phyloseq and beta-diversity analysis completed successfully.")
message("Samples analyzed: ", nsamples(ps_obj))
message("ASVs analyzed: ", ntaxa(ps_obj))
message("==============================================")


# =========================================================
# Save session info
# =========================================================

message("[OUTPUT] Saving session information...")

save_session_info(
  filename = "sessionInfo_16S_PRJNA759642_phyloseq.txt",
  dir = paths$logs,
  label = "16S rRNA sequencing - Phyloseq and beta-diversity analysis (PRJNA759642)"
)

message("[OUTPUT] Session information saved to: ", paths$logs)

# =========================================================
# Final pipeline message
# =========================================================

message("=================================================")
message("[PIPELINE] Phyloseq construction and beta-diversity analysis completed successfully for PRJNA759642.")
message("=================================================")

