# =========================================================
# Project      : 16S rRNA Sequencing - Microbiome Analysis Pipeline
# Dataset      : PRJNA759642
# Samples      : Human gut microbiome, 16S rRNA V4 amplicon sequencing
# Script       : Alpha-Diversity Analysis
# Description  : Load the phyloseq object and sample metadata,
#                aggregate ASV counts at genus level,
#                calculate Shannon alpha diversity,
#                assess the relationship between sequencing depth
#                and alpha diversity,
#                evaluate the joint contributions of disease status,
#                age, and gender using linear regression,
#                and perform rarefaction-based sensitivity analysis
#                following the published study.
# =========================================================


# =========================================================
#                  Libraries & Setup
# =========================================================

source("Setup_Environment/00_paths.R")
source("Setup_Environment/01_environment.R")
source("Setup_Environment/02_io_helpers.R")
source("Setup_Environment/03_seed.R")


# =========================================================
# Load and Prepare Data
# =========================================================

message("Load phyloseq object...")

ps_obj_unfiltered <- readRDS(file.path(paths$obj, "phyloseq_object_unfiltered.rds"))

# ---------------------------------------------------------
# Aggregate ASVs at genus level
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
# Calculate library size per sample
# ---------------------------------------------------------

counts_genus <- as(otu_table(ps_genus), "matrix")
if (taxa_are_rows(ps_genus)) counts_genus <- t(counts_genus)
lib_sizes <- rowSums(counts_genus)
stopifnot(setequal(names(lib_sizes), sample_names(ps_genus)))

meta_alpha <- data.frame(sample_data(ps_genus))

message("[QC] Genus-table library size - min: ", min(lib_sizes),
        ", median: ", median(lib_sizes), ", max: ", max(lib_sizes))

# =========================================================
# Shannon Alpha Diversity
# =========================================================

# ---------------------------------------------------------
# Calculate Shannon alpha diversity
# ---------------------------------------------------------

message("[ANALYSIS] Calculating Shannon alpha diversity at genus level...")

alpha_diversity <- phyloseq::estimate_richness(
  ps_genus,
  measures = "Shannon"
)

message("[QC] Shannon alpha diversity calculated for ",
        nrow(alpha_diversity),
        " samples.")

alpha_diversity$sample <- rownames(alpha_diversity)

# ---------------------------------------------------------
# Validate and align sample identifiers
# ---------------------------------------------------------

message("[CHECK] Validating and aligning sample identifiers...")

stopifnot(setequal(alpha_diversity$sample, sample_names(ps_genus)))

stopifnot(!anyDuplicated(alpha_diversity$sample))

alpha_diversity <- alpha_diversity[
  match(sample_names(ps_genus), alpha_diversity$sample),
]

stopifnot(identical(sample_names(ps_genus), alpha_diversity$sample))

alpha_diversity$SampleType <- sample_data(ps_genus)$SampleType

stopifnot(!anyNA(alpha_diversity$SampleType))

message(
  "[CHECK] Sample identifiers and metadata aligned successfully for ",
  nrow(alpha_diversity),
  " samples."
)

# ---------------------------------------------------------
# Add Shannon diversity to phyloseq sample metadata
# ---------------------------------------------------------

sample_data(ps_genus)$Shannon_diversity <- alpha_diversity$Shannon

message("[OUTPUT] Shannon diversity added to phyloseq sample metadata.")


# =========================================================
# Shannon Diversity Distribution and Group Comparison
# =========================================================

# ---------------------------------------------------------
# Shannon diversity distribution
# ---------------------------------------------------------

shannon_hist <- ggplot(
  as.data.frame(sample_data(ps_genus)),
  aes(x = Shannon_diversity)
) +
  geom_histogram(
    bins = 20,
    fill = "#7A8FA6",
    color = "black"
  ) +
  labs(
    title = "Distribution of Shannon diversity",
    x = "Shannon index",
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

save_plot(shannon_hist, filename = "shannon_diversity_distribution.png", dir = paths$plot_alpha)

# ---------------------------------------------------------
# Wilcoxon test between clinical groups
# ---------------------------------------------------------

message("[ANALYSIS] Comparing Shannon diversity between Healthy and UC samples...")

wilcoxon_shannon <- wilcox.test(
  Shannon ~ SampleType,
  data = alpha_diversity,
  correct = FALSE
)

message(
  "[STAT] Wilcoxon rank-sum test p-value: ",
  signif(wilcoxon_shannon$p.value, 4)
)

# ---------------------------------------------------------
# Shannon diversity by clinical group
# ---------------------------------------------------------

message("[ANALYSIS] Generating Shannon diversity boxplot...")

shannon_box <- ggplot(
  as.data.frame(sample_data(ps_genus)),
  aes(
    x = SampleType,
    y = Shannon_diversity,
    fill = SampleType
  )
) +
  geom_boxplot(
    color = "black",
    width = 0.6
  ) +
  geom_jitter(
    width = 0.1,
    size = 2,
    alpha = 0.7
  ) +
  scale_fill_manual(
    values = c(
      "Healthy" = "#0072B2",
      "UC" = "#D55E00"
    )
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(sample_data(ps_genus)$Shannon_diversity) * 1.05,
    label = paste0("Wilcoxon p = ", format.pval(wilcoxon_shannon$p.value, digits = 3)),
    size = 4
  ) +
  labs(
    title = "Shannon diversity by clinical group",
    x = "Sample type",
    y = "Shannon index"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "none"
  )

save_plot(shannon_box, filename = "shannon_diversity_by_sample_type.png", dir = paths$plot_alpha)


# =========================================================
# Linear Regression Analysis
# =========================================================

# ---------------------------------------------------------
# Prepare regression dataset
# ---------------------------------------------------------

message("[ANALYSIS] Preparing data for Shannon diversity linear regression...")

metadata_regression <- data.frame(sample_data(ps_genus))

metadata_regression$sample <- rownames(metadata_regression)

alpha_regression <- merge(
  alpha_diversity[, c("sample", "Shannon")],
  metadata_regression[, c("sample", "SampleType", "Age", "Gender")],
  by = "sample",
  sort = FALSE
)

alpha_regression$SampleType <- factor(
  alpha_regression$SampleType,
  levels = c("Healthy", "UC")
)

alpha_regression$Gender <- factor(alpha_regression$Gender)

alpha_regression$Age <- as.numeric(alpha_regression$Age)

alpha_regression$log_Shannon <- log(alpha_regression$Shannon)

stopifnot(!anyNA(alpha_regression$log_Shannon))
stopifnot(!anyNA(alpha_regression$SampleType))
stopifnot(!anyNA(alpha_regression$Age))
stopifnot(!anyNA(alpha_regression$Gender))

message(
  "[CHECK] Regression dataset prepared: ",
  nrow(alpha_regression),
  " samples."
)

# ---------------------------------------------------------
# Fit linear regression model
# ---------------------------------------------------------

message("[ANALYSIS] Fitting linear regression model...")

shannon_lm <- lm(
  log_Shannon ~ SampleType + Age + Gender,
  data = alpha_regression
)

print(summary(shannon_lm))
print(round(confint(shannon_lm), 3))

message("[STAT] UC vs Healthy, Shannon fold-change: ",
        round(exp(coef(shannon_lm)["SampleTypeUC"]), 3))

# ---------------------------------------------------------
# Assess sequencing depth effects on alpha diversity
# ---------------------------------------------------------

message("[ANALYSIS] Assessing the association between sequencing depth and alpha diversity...")

alpha_regression$LibrarySize     <- unname(lib_sizes[alpha_regression$sample])
alpha_regression$Observed_genera <- unname(rowSums(counts_genus > 0)[alpha_regression$sample])

stopifnot(!anyNA(alpha_regression$LibrarySize))
stopifnot(!anyNA(alpha_regression$Observed_genera))

message("[STAT] Testing correlation between library size and Shannon diversity...")

print(cor.test(alpha_regression$LibrarySize, alpha_regression$Shannon,
               method = "spearman", exact = FALSE))

message("[STAT] Testing correlation between library size and observed genera...")

print(cor.test(alpha_regression$LibrarySize, alpha_regression$Observed_genera,
               method = "spearman", exact = FALSE))

# ---------------------------------------------------------
# Sensitivity analysis: sequencing depth as a covariate
# ---------------------------------------------------------

message("[SENSITIVITY] Refitting Shannon model with sequencing depth as a covariate...")

shannon_lm_depth <- lm(log_Shannon ~ SampleType + Age + Gender + log10(LibrarySize),
                       data = alpha_regression)

print(summary(shannon_lm_depth))

# ---------------------------------------------------------
# Generate genus-level rarefaction curves
# ---------------------------------------------------------

message("[QC] Generating genus-level rarefaction curves...")

png(file.path(paths$plot_phy, "rarefaction_curves_genus.png"),
    width = 1800, height = 1200, res = 200)
vegan::rarecurve(counts_genus, step = 500, xlim = c(0, 2e5), label = FALSE,
                 xlab = "Number of sequences", ylab = "Number of genera",
                 col = ifelse(meta_alpha[rownames(counts_genus), "SampleType"] == "UC",
                              "#D55E00", "#0072B2"))
dev.off()


# =========================================================
# Rarefaction Sensitivity Analysis
# =========================================================

rarefaction_levels <- c(1000, 5000, 10000, 30000, 50000, 100000)
n_rarefactions <- 100
drop_low_depth <- FALSE

message("[ANALYSIS] Assessing alpha diversity across sequencing depths...")

for (depth in rarefaction_levels) {
  message(
    "[QC] Depth ", depth, ": ",
    sum(lib_sizes < depth), " of ", length(lib_sizes),
    " samples below the level (retained with all reads)"
  )
}

# ---------------------------------------------------------
# Rarefy samples at a given sequencing depth
# ---------------------------------------------------------

rarefy_keep_low_depth <- function(counts, depth, seed) {
  set.seed(seed)
  
  out <- counts
  enough <- rowSums(counts) >= depth
  
  if (any(enough)) {
    out[enough, ] <- vegan::rrarefy(
      counts[enough, , drop = FALSE],
      sample = depth
    )
  }
  
  out
}

# ---------------------------------------------------------
# Rarefy repeatedly and average the Shannon index
# ---------------------------------------------------------

message(
  "[ANALYSIS] Calculating Shannon diversity on ",
  n_rarefactions,
  " rarefactions per sequencing depth..."
)

rarefaction_shannon <- do.call(
  rbind,
  lapply(rarefaction_levels, function(depth) {
    
    # samples x repetitions matrix of Shannon values
    shannon_reps <- vapply(
      seq_len(n_rarefactions),
      function(i) {
        vegan::diversity(
          rarefy_keep_low_depth(counts_genus, depth, seed = 1234 + i),
          index = "shannon"
        )
      },
      numeric(nrow(counts_genus))
    )
    
    data.frame(
      sample = rownames(counts_genus),
      Shannon = rowMeans(shannon_reps),
      Shannon_SD = apply(shannon_reps, 1, sd),
      Depth = depth,
      stringsAsFactors = FALSE
    )
  })
)

if (drop_low_depth) {
  keep <- unname(lib_sizes[rarefaction_shannon$sample]) >= rarefaction_shannon$Depth
  rarefaction_shannon <- rarefaction_shannon[keep, , drop = FALSE]
}

rarefaction_shannon$SampleType <- factor(
  meta_alpha[rarefaction_shannon$sample, "SampleType"],
  levels = c("Healthy", "UC")
)

stopifnot(!anyNA(rarefaction_shannon$SampleType))

# ---------------------------------------------------------
# Compare Shannon diversity between groups
# ---------------------------------------------------------

message("[STAT] Comparing Shannon diversity between Healthy and UC samples...")

wilcoxon_rarefaction <- do.call(
  rbind,
  lapply(
    split(rarefaction_shannon, rarefaction_shannon$Depth),
    function(df) {
      
      depth <- unique(df$Depth)
      
      data.frame(
        Depth = depth,
        P_value = wilcox.test(
          Shannon ~ SampleType,
          data = df,
          correct = FALSE
        )$p.value,
        N_samples = nrow(df),
        N_below_depth = sum(lib_sizes < depth)
      )
    }
  )
)

print(wilcoxon_rarefaction)

wilcoxon_rarefaction$P_label <- paste0(
  "Wilcoxon p = ",
  formatC(
    wilcoxon_rarefaction$P_value,
    digits = 3,
    format = "g"
  )
)

# ---------------------------------------------------------
# Generate rarefaction sensitivity plot
# ---------------------------------------------------------

message("[ANALYSIS] Generating rarefaction sensitivity plot...")

p_alpha_rarefaction <- ggplot(
  rarefaction_shannon,
  aes(x = SampleType, y = Shannon, fill = SampleType)
) +
  geom_boxplot(color = "black", width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.7) +
  geom_text(
    data = wilcoxon_rarefaction,
    aes(x = 1.5, y = Inf, label = P_label),
    inherit.aes = FALSE, vjust = 1.5, size = 3.5
  ) +
  facet_wrap(
    ~ Depth, nrow = 2, scales = "fixed",
    labeller = labeller(Depth = c(
      "1000" = "1K reads", "5000" = "5K reads", "10000" = "10K reads",
      "30000" = "30K reads", "50000" = "50K reads", "1e+05" = "100K reads"
    ))
  ) +
  scale_fill_manual(values = c("Healthy" = "#0072B2", "UC" = "#D55E00")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(title = "Shannon diversity across sequencing depths",
       subtitle = paste0("Mean Shannon index over ", n_rarefactions, " rarefactions per depth"),
       x = "Sample type", y = "Shannon index") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

save_plot(p_alpha_rarefaction, filename = "shannon_rarefaction_sensitivity.png", dir = paths$plot_alpha)


# =========================================================
# Save Analysis Objects
# =========================================================

message("[OUTPUT] Saving alpha-diversity analysis objects...")

save_rds(ps_genus, filename = "phyloseq_object_genus.rds", dir = paths$obj)

save_rds(alpha_diversity, filename = "alpha_diversity_shannon.rds", dir = paths$obj)

save_rds(wilcoxon_shannon, filename = "wilcoxon_shannon.rds", dir = paths$obj)

save_rds(shannon_lm, filename = "shannon_linear_regression.rds", dir = paths$obj)

save_rds(shannon_lm_depth, filename = "shannon_linear_regression_depth.rds", dir = paths$obj)

save_rds(rarefaction_shannon, filename = "shannon_rarefaction_results.rds", dir = paths$obj)

save_rds(wilcoxon_rarefaction, filename = "wilcoxon_shannon_rarefaction.rds", dir = paths$obj)

message("[OUTPUT] Alpha-diversity analysis objects saved successfully.")


# =========================================================
# Save session info
# =========================================================

message("[OUTPUT] Saving session information...")

save_session_info(
  filename = "sessionInfo_16S_PRJNA759642_alpha_diversity.txt",
  dir = paths$logs,
  label = "16S rRNA sequencing - Alpha-diversity analysis (PRJNA759642)"
)

message("[OUTPUT] Session information saved to: ", paths$logs)


# =========================================================
# Final pipeline message
# =========================================================

message("=================================================")
message("[PIPELINE] Alpha-diversity analysis completed successfully for PRJNA759642.")
message("=================================================")


