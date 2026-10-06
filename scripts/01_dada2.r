# =========================================================
# Project      : 16S rRNA Sequencing - Microbiome Analysis Pipeline
# Dataset      : PRJNA759642
# Samples      : Human gut microbiome, 16S rRNA V4 amplicon sequencing
# Script       : DADA2 Processing and Taxonomic Assignment
# Description  : Load raw paired-end FASTQ files, inspect read quality,
#                filter and trim reads, learn error rates, infer ASVs,
#                merge paired-end reads, remove chimeric sequences,
#                and assign taxonomy using the SILVA 138.2 reference database.
#                and assign genus- and species-level taxonomy using
#                the SILVA 138.2 reference databases.
# =========================================================


# =========================================================
#                  Libraries & Setup
# =========================================================

source("Setup_Environment/00_paths.R")
source("Setup_Environment/01_environment.R")
source("Setup_Environment/02_io_helpers.R")
source("Setup_Environment/03_seed.R")


# =========================================================
# Load raw sequencing data
# =========================================================

message("Loading raw sequencing data")

raw_fastq <- paths$data_raw_16sfastq

fns <- sort(list.files(raw_fastq, full.names = TRUE))

forward <- fns[grepl("_1\\.fastq\\.gz$", fns)]
reverse <- fns[grepl("_2\\.fastq\\.gz$", fns)]

message("Forward FASTQ files: ", length(forward))
message("Reverse FASTQ files: ", length(reverse))

if (length(forward) != length(reverse)) {
  message("ERROR: Forward and reverse FASTQ counts do not match.")
  stop()
}

message("Paired-end samples identified: ", length(forward))


# =========================================================
# Inspect read quality profiles
# =========================================================

message("Inspecting read quality profiles")

list_qc_forward <- list()
list_qc_reverse <- list()

# -------------------------
# Forward reads
# -------------------------

message("Generating forward-read quality profiles")

for (i in seq_along(forward)) {
  
  message("QC forward: ", basename(forward[i]))
  
  p <- plotQualityProfile(forward[i])
  
  list_qc_forward[[i]] <- p
  
  run_name <- sub("_1\\.fastq\\.gz$", "", basename(forward[i]))
  
  save_plot(p, filename = paste0(run_name, "_forward_QC.png"), dir = paths$plot_qc_fw,
         width = 8, height = 6)
  
  rm(p)
}


# -------------------------
# Reverse reads
# -------------------------

message("Generating reverse-read quality profiles")

for (i in seq_along(reverse)) {
  
  message("QC reverse: ", basename(reverse[i]))
  
  p <- plotQualityProfile(reverse[i])
  
  list_qc_reverse[[i]] <- p
  
  run_name <- sub("_2\\.fastq\\.gz$", "", basename(reverse[i]))
  
  save_plot(p, filename = paste0(run_name, "_reverse_QC.png"), dir = paths$plot_qc_rv,
            width = 8, height = 6)
  
  rm(p)
}

# -------------------------
# Aggregate quality profiles
# -------------------------

message("Generating aggregate quality profiles")

qc_fw_aggregate <- plotQualityProfile(forward, aggregate = TRUE)
qc_rv_aggregate <- plotQualityProfile(reverse, aggregate = TRUE)

save_plot(plot = qc_fw_aggregate, filename = "QC_forward_aggregate.png", dir = paths$plot_qc)
save_plot(plot = qc_rv_aggregate, filename = "QC_reverse_aggregate.png", dir = paths$plot_qc)

message("Read quality inspection completed.")


# =========================================================
# Filter and trim reads
# ========================================================= 

message("Filtering and trimming reads")

filtFs <- file.path(paths$data_16sfiltered, paste0(sub("_1\\.fastq\\.gz$", "", basename(forward)), "_F_filt.fastq.gz"))
filtRs <- file.path(paths$data_16sfiltered, paste0(sub("_2\\.fastq\\.gz$", "", basename(reverse)), "_R_filt.fastq.gz"))

names(filtFs) <- sub("_1\\.fastq\\.gz$", "", basename(forward))
names(filtRs) <- sub("_2\\.fastq\\.gz$", "", basename(reverse))

sample.names <- names(filtFs)

message("Samples to be filtered: ", length(sample.names))

out <- filterAndTrim(
  fwd = forward, filt = filtFs,
  rev = reverse, filt.rev = filtRs,
  truncLen = c(140, 130),
  maxEE = c(2, 2),
  truncQ = 2,
  maxN = 0,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE,  
  verbose = TRUE
)

message("Filtering and trimming completed.")


# =========================================================
# Learn DADA2 error rates
# =========================================================

message("Learning DADA2 error rates")

errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

message("Error models successfully learned.")

# -------------------------
# Plot error models
# -------------------------

message("Generating DADA2 error model plots")

errF_plot <- plotErrors(errF, nominalQ = TRUE) +
  ggtitle("DADA2 Error Model — Forward Reads")

errR_plot <- plotErrors(errR, nominalQ = TRUE) +
  ggtitle("DADA2 Error Model — Reverse Reads")


save_plot(errF_plot, filename = "dada2_error_model_forward.png", dir = paths$plot_qc)
save_plot(errR_plot, filename = "dada2_error_model_reverse.png", dir = paths$plot_qc)


# =========================================================
# Infer amplicon sequence variants
# =========================================================

message("Inferring ASVs with DADA2")

# Standard independent sample inference
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

message("ASV inference completed.")


# =========================================================
# Merge paired-end reads
# =========================================================

message("Merging paired-end reads")

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

message("Paired-end read merging completed.")


# =========================================================
# Construct ASV sequence table
# =========================================================

message("Constructing ASV sequence table")

seqtab <- makeSequenceTable(mergers)

message(
  "Sequence table dimensions: ",
  nrow(seqtab),
  " samples × ",
  ncol(seqtab),
  " ASVs"
)

# -------------------------
# Inspect ASV length distribution
# -------------------------

seq_lengths <- table(
  nchar(getSequences(seqtab))
)

print(seq_lengths)


# =========================================================
# Remove chimeric sequences
# =========================================================

message("Removing chimeric sequences")

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread=TRUE, verbose=TRUE)

message(
  "Non-chimeric sequence table dimensions: ",
  nrow(seqtab.nochim),
  " samples × ",
  ncol(seqtab.nochim),
  " ASVs"
)

chimera_retention <- sum(seqtab.nochim) / sum(seqtab)

message(
  "Proportion of reads retained after chimera removal: ",
  round(chimera_retention, 4)
)


# =========================================================
# Track reads through the DADA2 pipeline
# =========================================================

message("Tracking reads through the DADA2 pipeline")

getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

message("Read tracking table generated.")


# =========================================================
# Assign taxonomy
# =========================================================

message("Assigning taxonomy")

taxa <- assignTaxonomy(seqtab.nochim, refFasta = file.path(paths$ref, "silva_nr99_v138.2_toGenus_trainset.fa.gz"), multithread = TRUE)

message("Genus-level taxonomy assignment completed.")

# -------------------------
# Species-level assignment
# -------------------------

message("Assigning species by exact sequence matching")

taxa <- addSpecies(taxa, refFasta = file.path(paths$ref, "silva_v138.2_assignSpecies.fa.gz"))

taxa.print <- taxa 

rownames(taxa.print) <- NULL

message("Species-level assignment completed.")


# =========================================================
# Save DADA2 results
# =========================================================

message("Saving DADA2 results")

save_rds(taxa, filename = "taxa.rds", dir = paths$obj)
save_rds(seqtab.nochim, filename = "seqtab_nochim.rds", dir = paths$obj)
save_rds(track, filename = "track_table.rds", dir = paths$obj)

message("DADA2 objects saved successfully.")


# =========================================================
# DADA2 processing completed
# =========================================================

message("==============================================")
message("DADA2 processing completed successfully.")
message("Samples processed: ", nrow(seqtab.nochim))
message("Non-chimeric ASVs: ", ncol(seqtab.nochim))
message("==============================================")


# =========================================================
# Save session info
# =========================================================

message("[OUTPUT] Saving session information...")

save_session_info(
  filename = "sessionInfo_16S_PRJNA759642_DADA2.txt",
  dir = paths$logs,
  label = "16S rRNA sequencing - DADA2 processing and taxonomic assignment (PRJNA759642)"
)

message("[OUTPUT] Session information saved to: ", paths$logs)


# =========================================================
# Final pipeline message
# =========================================================

message("=================================================")
message("[PIPELINE] DADA2 processing and taxonomic assignment completed successfully for PRJNA759642.")
message("=================================================")

