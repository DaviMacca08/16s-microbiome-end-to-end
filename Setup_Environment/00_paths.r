# =========================================================
# Project paths (single source of truth)
# =========================================================

base_dir <- "/Users/davidemaccarrone/Desktop/Bioinformatics/MyProjects🤞🏻/Metagenomics/16srRNA/v1"

paths <- list(
  
  base = base_dir,
  
  # data 
  data_raw_16sfastq = file.path(base_dir, "data/raw_fastq_16S_main_cohort"),
  data_16sfiltered = file.path(base_dir, "data/raw_fastq_16s_main_cohort_filtered"),
  metadata = file.path(base_dir, "data/metadata"),
  
  # plot
  plot_qc = file.path(base_dir, "results/plot/qc/"),
  plot_qc_fw = file.path(base_dir, "results/plot/qc/forward"),
  plot_qc_rv = file.path(base_dir, "results/plot/qc/reverse"),
  plot_phy = file.path(base_dir, "results/plot/phyloseq_beta_alpha"),
  plot_beta = file.path(base_dir, "results/plot/phyloseq_beta_alpha/beta"),
  plot_alpha = file.path(base_dir, "results/plot/phyloseq_beta_alpha/alpha"),
  plot_diff_abu = file.path(base_dir, "results/plot/differential_abudance"),
  
  # objects
  obj = file.path(base_dir, "results/objects"),
  
  # log
  logs = file.path(base_dir, "results/log"),
  
  # references
  ref = file.path(base_dir, "references")
)

# Ensure directories exist
invisible(lapply(paths, function(x) {
  if (!dir.exists(x)) dir.create(x, recursive = TRUE, showWarnings = FALSE)
}))