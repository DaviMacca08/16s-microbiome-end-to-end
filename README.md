End-to-End 16S rRNA Microbiome Analysis Framework: From Raw Sequencing
Data to Differential Abundance
================

![R](https://img.shields.io/badge/R-4.6.0-blue)
![16S-rRNA](https://img.shields.io/badge/16S--rRNA-DADA2-purple)
![Microbiome](https://img.shields.io/badge/Microbiome-phyloseq-green)
![Status](https://img.shields.io/badge/status-completed-success)

------------------------------------------------------------------------

# 💼 Bioinformatics Service Demonstration

This repository presents a modular **end-to-end 16S rRNA microbiome
analysis pipeline** designed for the systematic characterization of
microbial communities from amplicon sequencing data.

The workflow covers the main analytical stages required for a
reproducible 16S rRNA sequencing project, from raw sequencing quality
assessment and denoising to taxonomic profiling, microbial diversity
analysis and differential abundance testing.

The project is structured as a practical bioinformatics workflow
suitable for:

- exploratory analysis of 16S rRNA amplicon sequencing datasets
- characterization of microbial community composition
- comparison of microbiome profiles across biological conditions
- assessment of within-sample and between-sample microbial diversity
- identification of taxa associated with experimental or clinical groups
- reproducible analysis of microbiome sequencing data
- development of client-oriented microbiome bioinformatics workflows

**Core analytical functionalities implemented in the pipeline include:**

- raw FASTQ quality assessment using FastQC
- aggregated quality reporting using MultiQC
- paired-end read quality assessment, including verification of
  primer/adapter absence
- amplicon sequence variant (ASV) inference using DADA2
- sequencing-error modeling and denoising
- paired-end read merging
- chimera removal
- construction of an ASV feature table
- taxonomic classification against SILVA 138.2 (genus level, plus
  species-level exact matching)
- integration of sequencing data and sample metadata using phyloseq
- phylum-level taxonomic composition profiling
- genus-level alpha diversity analysis (Shannon index), including
  covariate-adjusted linear regression
- genus-level beta diversity analysis (Bray-Curtis), PCoA ordination,
  PERMANOVA, dispersion testing and Multiple Regression on distance
  Matrices (MRM)
- sensitivity analysis of alpha diversity to sequencing depth
  (library-size covariate and repeated rarefaction)
- family-level differential abundance analysis using a negative-binomial
  framework
- TMM normalization using edgeR
- multiple-testing correction using Benjamini–Hochberg FDR

**📤 Deliverables**

- Quality-control summaries
- Processed ASV feature tables
- Taxonomic assignments
- Microbiome composition profiles
- Alpha and beta diversity results
- Differential abundance results
- Reproducible analytical scripts
- Structured visualization and result tables
- Analytical report documenting the complete workflow

------------------------------------------------------------------------

# 📊 Case Study

This repository uses a publicly available human gut microbiome dataset
as a demonstration case for the implementation of a reproducible 16S
rRNA microbiome analysis workflow.

**Organism:** Homo sapiens

**Biological system:** Human gut microbiome (fecal samples)

**Cohort:** Pediatric cohort (ages 7–21 years), 42 samples: 19
ulcerative colitis (UC) and 23 healthy controls

**Experimental groups:** Healthy and ulcerative colitis (UC)

**Sequencing approach:** 16S rRNA gene amplicon sequencing

**Target region:** V4 region (515F/806R primers)

**Sequencing:** Paired-end Illumina MiSeq, 2 × 150 bp

**Public dataset:** NCBI Sequence Read Archive (SRA), PRJNA759642

**Source study:** Zuo et al. (2022), *Scientific Reports*, DOI:
10.1038/s41598-022-07995-7

**Primary feature type:** Amplicon Sequence Variants (ASVs)

The dataset contains human gut microbiome samples generated using 16S
rRNA V4 amplicon sequencing and provides a suitable case study for
demonstrating an end-to-end microbial community analysis workflow. The
original study also generated shotgun metagenomic data; only the 16S
rRNA data are analyzed here.

------------------------------------------------------------------------

# 🧬 Biological Context

16S rRNA gene amplicon sequencing provides a targeted approach for
characterizing bacterial community composition by sequencing a conserved
phylogenetic marker containing variable regions that enable taxonomic
discrimination.

The analysis of microbiome sequencing data requires consideration of
several characteristics that distinguish these datasets from
conventional bulk transcriptomic data.

In particular:

- sequencing depth can vary substantially between samples
- microbial abundance tables are sparse
- microbial sequencing data are compositional
- taxonomic features can differ substantially in prevalence and
  abundance
- technical variation can influence observed community structure
- statistical inference must account for the count-based nature of
  sequencing data

The pipeline therefore separates **data processing, ecological
characterization and statistical testing** into distinct analytical
stages.

The workflow begins with raw sequencing quality assessment and proceeds
through denoising and ASV inference, taxonomic classification and
construction of an integrated microbiome data object.

The resulting feature table is subsequently used to characterize:

- microbial taxonomic composition
- within-sample diversity
- between-sample community structure
- condition-associated differences in microbial abundance

The use of ASVs rather than fixed similarity-based OTUs follows modern
amplicon-sequencing workflows in which sequencing errors are modeled
explicitly and exact biological sequence variants can be resolved. DADA2
was specifically developed to model Illumina amplicon sequencing errors
and infer high-resolution sequence variants.

------------------------------------------------------------------------

# 📌 Analytical Framework

The complete workflow can be summarized as:

**raw FASTQ files → quality control → DADA2 filtering and denoising →
ASV inference → paired-end merging → chimera removal → taxonomy
assignment → phyloseq integration → community composition → alpha
diversity → beta diversity → differential abundance**

Each stage is designed to generate structured outputs that can be used
by subsequent analytical modules.

## Sequencing Quality Control

Initial sequencing quality is assessed on the 42 main-cohort runs using
**FastQC**, followed by aggregation of sample-level reports with
**MultiQC** (`qc.sh`).

The QC stage is used to evaluate:

- per-base sequence quality
- sequence length distributions
- adapter/primer-related sequence content
- duplication patterns
- general sequencing quality
- consistency between sequencing samples

Quality assessment is performed before downstream denoising to identify
potential technical problems and to guide filtering decisions.

Some FastQC modules (per-base sequence content, sequence duplication,
overrepresented sequences, GC content) are flagged as failed or warning
for this dataset. These flags are expected for amplicon data, where all
reads originate from the same short, low-diversity locus, and do not
indicate technical problems.

## Amplicon Processing and ASV Inference

The preprocessing stage operates on paired-end FASTQ data.

Primers had already been removed by the data submitters before SRA
deposition, and the absence of adapter/primer contamination was
confirmed at the QC stage; no additional primer-trimming step was
therefore required.

DADA2 is used for:

- quality filtering and trimming (reads truncated at 140 bp forward /
  130 bp reverse, based on the per-sample and aggregate quality
  profiles; maximum 2 expected errors per read; no ambiguous bases; PhiX
  removal)
- error-rate learning, performed separately for forward and reverse
  reads
- independent, per-sample denoising (no sample pooling)
- inference of exact sequence variants
- paired-end read merging
- chimera detection and removal (consensus method)

The resulting ASV table provides the primary microbial feature matrix
for downstream analysis.

DADA2 was selected because it explicitly models sequencing errors and
enables inference of exact amplicon sequence variants rather than
relying exclusively on arbitrary similarity thresholds.

## Taxonomic Classification

The inferred ASVs are assigned taxonomic identities with the DADA2 naive
Bayesian classifier (`assignTaxonomy`) trained on the SILVA 138.2
non-redundant (99%) reference database, followed by species-level
assignment by exact sequence matching (`addSpecies`).

Taxonomic annotations are subsequently integrated with the ASV abundance
table to generate a structured microbiome object suitable for downstream
ecological and statistical analyses.

## Microbiome Data Integration

Processed ASV counts, taxonomy and sample metadata are integrated using
**phyloseq**.

The phyloseq framework provides a structured representation of
microbiome count data, taxonomy and sample information and supports
filtering, agglomeration, diversity analysis, ordination and
visualization within a reproducible R workflow.

Sample and ASV identifiers are validated between the sequence table,
taxonomy and metadata before the phyloseq object is built.

Taxonomic aggregation is performed at the level appropriate to each
analysis: phylum level for composition profiling, genus level for alpha
and beta diversity, and family level for differential abundance.
Aggregation is applied to the unfiltered ASV table, so no prevalence
filter is applied before the diversity analyses.

Quality checks at this stage include the sequencing-depth distribution
and an ASV prevalence-versus-abundance overview. Phylum-level
composition is shown as relative abundance per sample, with phyla below
2% mean relative abundance grouped as “Other”.

------------------------------------------------------------------------

# 📈 Diversity Analysis

## Alpha Diversity

Alpha diversity is used to characterize microbial diversity **within
individual samples**.

The analysis uses the **Shannon diversity index**, which incorporates
both feature richness and relative abundance, computed on
genus-aggregated counts.

The analysis is designed to assess whether the observed within-sample
microbial diversity differs between Healthy and UC samples, and
includes:

- a group comparison using the Wilcoxon rank-sum test
- a linear regression of log-transformed Shannon diversity on disease
  status, age and gender, reporting the UC-versus-Healthy fold-change
  and 95% confidence intervals
- an assessment of the relationship between sequencing depth and
  diversity (see *Sequencing-Depth Sensitivity Analysis* below)

## Beta Diversity

Beta diversity is used to characterize differences in microbial
community composition **between samples**.

Beta diversity is computed on genus-level relative abundances using
**Bray-Curtis** dissimilarity. The analysis includes:

- Principal Coordinates Analysis (PCoA) to visualize relationships among
  samples
- **PERMANOVA** (`vegan::adonis2`, 999 permutations, marginal effects)
  testing disease status while adjusting for age and gender
- a test of homogeneity of multivariate dispersion (`betadisper` and
  `permutest`), to verify that PERMANOVA results are not driven by
  differences in within-group dispersion
- **Multiple Regression on distance Matrices (MRM)** (`ecodist::MRM`,
  999 permutations), with pairwise predictors for disease-status
  combination, gender combination and age difference

Beta-diversity analyses are interpreted at the community level rather
than as direct measurements of individual taxon abundance.

## Sequencing-Depth Sensitivity Analysis

Because sequencing depth can influence observed diversity estimates, the
Shannon analysis is complemented by the following checks:

- Spearman correlation between library size and both Shannon diversity
  and the number of observed genera
- refitting of the Shannon regression model with log10 library size as
  an additional covariate
- genus-level rarefaction curves
- repeated rarefaction at 1,000, 5,000, 10,000, 30,000, 50,000 and
  100,000 reads (100 rarefactions per depth, averaged into a mean count
  table), with a Wilcoxon comparison between Healthy and UC at each
  depth; samples with fewer reads than the target depth are retained
  with all their reads

This sensitivity analysis follows the approach of the original
publication.

Importantly, rarefaction is used here as a **sensitivity analysis for
sequencing-depth effects**, rather than as the general normalization
strategy for all downstream statistical analyses. This distinction is
important because rarefaction discards sequencing observations and has
been criticized as an inefficient normalization strategy for
quantitative microbiome analyses.

------------------------------------------------------------------------

# 📊 Differential Abundance Analysis

Differential abundance analysis is performed to identify microbial taxa
whose observed abundance differs between the biological conditions (UC
versus Healthy). Testing is performed at the **family level**.

The implemented statistical framework is based on **edgeR** and
includes:

- aggregation of ASV counts at family level
- filtering of low-information families (zero counts in at least 90% of
  samples, or count variance below half of the median variance across
  families)
- construction of count-based feature matrices and an edgeR `DGEList`
- TMM normalization using `calcNormFactors()`
- exploratory analysis of the normalized data (log2-CPM distributions,
  MDS, PCA, sample-distance heatmap)
- estimation of dispersion parameters using `estimateDisp()`
- negative-binomial modeling and two-group testing using `exactTest()`
- Benjamini–Hochberg false-discovery-rate correction, with families at
  FDR \< 0.25 reported as significant
- volcano plot and heatmap of the significant families

The analysis is performed on count data rather than relying exclusively
on relative-abundance comparisons.

TMM normalization is a robust scaling approach designed to account for
differences in library composition and sequencing depth and has been
incorporated into several microbiome differential-abundance workflows.

The compositional nature of microbiome sequencing data is explicitly
considered when interpreting differential-abundance results.
Sequencing-based microbiome datasets represent relative observations
constrained by the total sequencing output, meaning that changes in one
taxon can influence the observed relative abundances of other taxa.

Therefore, differential-abundance results are interpreted as
**condition-associated differences in observed microbial abundance**,
rather than direct measurements of absolute microbial biomass.

Multiple-testing correction is performed using the **Benjamini–Hochberg
procedure** to control the false discovery rate across tested microbial
features. Families are considered significant at BH-FDR \< 0.25; the
rationale for this threshold is discussed in the report.

The exact test compares two groups and does not adjust for covariates:
age and gender are accounted for in the alpha- and beta-diversity
models, but not in the differential abundance test.

------------------------------------------------------------------------

# 🗂 Repository Contents

The repository is organized into six independent modules, executed in
order: two Bash scripts for data acquisition and raw-read quality
control, followed by four R scripts for the analysis. Each module reads
the files or objects produced by the previous ones and writes its own
outputs; the R scripts also save the session information of their run.

| Script | Analysis | Main outputs |
|:---|:---|:---|
| `download_16S_fastq.sh` | Download of the paired-end FASTQ files of the 42 runs of the main cohort from the ENA Portal API (parallel downloads, interrupted downloads are resumed) | `*_1.fastq.gz` / `*_2.fastq.gz` files (not versioned) and the filtered ENA run report |
| `qc.sh` | FastQC on every run and aggregation with MultiQC | FastQC reports and MultiQC report |
| `01_dada2.r` | Read-quality inspection, filtering and trimming, error learning, ASV inference, merging, chimera removal, SILVA taxonomic assignment | Quality profiles (per sample and aggregate), error-model plots, ASV table (`seqtab_nochim.rds`), taxonomy table (`taxa.rds`) |
| `02_phyloseq_beta_diversity.r` | phyloseq object construction and validation, ASV prevalence, sequencing depth, phylum composition, genus-level Bray-Curtis beta diversity | Prevalence, depth and composition plots, PCoA, dispersion boxplot; phyloseq objects, distance matrix, PERMANOVA, dispersion and MRM results |
| `03_phyloseq_alpha_diversity.r` | Genus-level Shannon diversity, group comparison, covariate-adjusted regression, sequencing-depth sensitivity analysis | Shannon distribution and group plots, rarefaction curves, rarefaction sensitivity plot; alpha-diversity table, Wilcoxon, regression and rarefaction results |
| `04_differential_abundance_analysis.r` | Family-level differential abundance with edgeR | Normalization, MDS, PCA, dispersion and sample-distance plots, volcano plot, heatmap; edgeR objects and result tables (`.rds` and `.csv`) |

The overall flow is:

**raw sequencing data → quality control → ASV inference → taxonomy →
microbiome profiling → diversity analysis → differential abundance**

and, within each module:

**analysis code → processed objects → result tables → visualization
outputs → report**

## Outputs

- **Figures** are saved in `results/plot/`, grouped by analytical
  module.
- **Processed objects** (ASV table, taxonomy, phyloseq objects, distance
  matrices, statistical results and edgeR objects) are saved as RDS
  files in `results/object/`; the differential abundance table is also
  exported as CSV.
- **Logs** with the session information of each script are saved in
  `results/log/`.

## Report

The results, figures and their biological interpretation are presented
in the analytical report in the `report/` directory, together with the
methodological decisions (filtering thresholds, significance criteria,
sensitivity analyses). This README documents what the project contains
and how to run it; it does not report results.

## Reproducibility

The workflow is designed around reproducible execution and transparent
methodological tracking. The project environment is managed using **R
and renv** (exact package versions are recorded in `renv.lock`), while
shared setup and utility functions are centralized within the
`Setup_Environment/` directory. Random seeds are fixed for permutation
tests and rarefactions.

------------------------------------------------------------------------

# ▶️ Reproducing the Analysis

1.  Clone the repository.
2.  Restore the R environment with `renv::restore()`.
3.  Download the raw FASTQ files of the 42 main-cohort runs of
    PRJNA759642 with `download_16S_fastq.sh` (requires Bash, `curl`,
    `awk` and `xargs`). Files are expected as `*_1.fastq.gz` /
    `*_2.fastq.gz`.
4.  Place the SILVA 138.2 files
    (`silva_nr99_v138.2_toGenus_trainset.fa.gz` and
    `silva_v138.2_assignSpecies.fa.gz`) in the reference directory and
    the cleaned metadata (`meta_clean.xlsx`) in `data/metadata/`. All
    locations are defined in `Setup_Environment/00_paths.R`.
5.  Run the raw-read quality control with `qc.sh`, from an environment
    providing FastQC and MultiQC (e.g. a conda environment).
6.  Run the R scripts in `scripts/` in numerical order (see *Repository
    Contents*).

> ⚠️ **Check the paths before running.** The Bash scripts use relative
> paths and directory names defined at the top of each file
> (e.g. `OUTDIR`, `RAW_DIR`, `FASTQC_DIR`, `MULTIQC_DIR`), and should be
> launched from the `scripts/` directory. These must be adapted to your
> local setup and must be consistent with each other: the directory
> where the FASTQ files are downloaded must be the one read by `qc.sh`
> and by the R scripts (defined in `Setup_Environment/00_paths.R`). The
> `conda activate` line in `qc.sh` must also be adapted to the name of
> your own environment.

Scripts 03 and 04 load the phyloseq object saved by script 02, which
must therefore be run first.

------------------------------------------------------------------------

# 📁 Project Structure

``` text
Metagenomics_16S_project/
├── data/
│   ├── metadata/
│   └── raw FASTQ files (not versioned)
│
├── results/
│   ├── fastqc/
│   ├── multiqc/
│   ├── plot/
│   │   └── Analysis-specific visualization outputs
│   │
│   ├── object/
│   │   └── Generated analytical summaries and RDS objects
│   │
│   └── log/
│       └── Execution logs and reproducibility information
│
├── report/
│   └── Analysis report
│
├── scripts/
│   ├── download_16S_fastq.sh
│   ├── qc.sh
│   ├── 01_dada2.r
│   ├── 02_phyloseq_beta_diversity.r
│   ├── 03_phyloseq_alpha_diversity.r
│   └── 04_differential_abundance_analysis.r
│
├── Setup_Environment/
│   ├── 00_paths.R
│   ├── 01_environment.R
│   ├── 02_io_helpers.R
│   └── 03_seed.R
│
├── README.md
└── LICENSE.txt
```

------------------------------------------------------------------------

# 🔬 Methodological Principles

The workflow follows several principles relevant to reproducible
microbiome analysis.

### ASV-based inference

The pipeline uses DADA2-based ASV inference rather than conventional
fixed-threshold OTU clustering.

ASVs provide high-resolution sequence features and allow sequencing
errors to be modeled during inference.

### Count-based statistical analysis

Where appropriate, the workflow retains count information for
statistical modeling rather than converting the complete feature table
to relative abundance before differential testing.

### Sequencing-depth assessment

Sequencing depth is explicitly evaluated during quality control, and its
influence on alpha diversity is tested directly through a library-size
covariate and repeated rarefaction.

### Compositionality-aware interpretation

Microbiome sequencing data are compositional by construction.
Consequently, observed relative abundances should not automatically be
interpreted as absolute microbial abundance.

### Multiple-testing control

Differential-abundance testing involves simultaneous statistical testing
across multiple microbial features. Benjamini–Hochberg FDR correction is
therefore applied to control the expected proportion of false
discoveries among significant results.

### Reproducibility

The project separates:

- environment configuration
- path management
- input/output utilities
- analytical scripts
- intermediate objects
- result tables
- visualization outputs
- analytical reports

This structure facilitates reproducibility, auditing and adaptation of
the workflow to new datasets.

------------------------------------------------------------------------

# 📚 References

- Zuo W, Wang B, Bai X, Luan Y, Fan Y, Michail S, Sun F. (2022). **16S
  rRNA and metagenomic shotgun sequencing data revealed consistent
  patterns of gut microbiome signature in pediatric ulcerative
  colitis.** *Scientific Reports*. DOI: 10.1038/s41598-022-07995-7.

- Callahan BJ, McMurdie PJ, Rosen MJ, et al. (2016). **DADA2:
  High-resolution sample inference from Illumina amplicon data.**
  *Nature Methods*, 13, 581–583. DOI: 10.1038/nmeth.3869.

- Nearing JT, Douglas GM, Comeau AM, Langille MGI. (2018). **Denoising
  the Denoisers: an independent evaluation of microbiome sequence
  error-correction approaches.** *PeerJ*, 6:e5364. DOI:
  10.7717/peerj.5364.

- Quast C, Pruesse E, Yilmaz P, et al. (2013). **The SILVA ribosomal RNA
  gene database project: improved data processing and web-based tools.**
  *Nucleic Acids Research*, 41(D1):D590–D596.

- McMurdie PJ, Holmes S. (2013). **phyloseq: An R Package for
  Reproducible Interactive Analysis and Graphics of Microbiome Census
  Data.** *PLoS ONE*, 8(4): e61217. DOI: 10.1371/journal.pone.0061217.

- McMurdie PJ, Holmes S. (2014). **Waste Not, Want Not: Why Rarefying
  Microbiome Data Is Inadmissible.** *PLoS Computational Biology*,
  10(4): e1003531. DOI: 10.1371/journal.pcbi.1003531.

- Robinson MD, McCarthy DJ, Smyth GK. (2010). **edgeR: a Bioconductor
  package for differential expression analysis of digital gene
  expression data.** *Bioinformatics*, 26(1):139–140.

- Robinson MD, Oshlack A. (2010). **A scaling normalization method for
  differential expression analysis of RNA-seq data.** *Genome Biology*,
  11:R25.

- Anderson MJ. (2001). **A new method for non-parametric multivariate
  analysis of variance.** *Austral Ecology*, 26:32–46.

- Lichstein JW. (2007). **Multiple regression on distance matrices: a
  multivariate spatial analysis tool.** *Plant Ecology*, 188:117–131.

- Gloor GB, Macklaim JM, Pawlowsky-Glahn V, Egozcue JJ. (2017).
  **Microbiome Datasets Are Compositional: And This Is Not Optional.**
  *Frontiers in Microbiology*, 8:2224.

- Nearing JT, Douglas GM, Hayes MG, et al. (2022). **Microbiome
  differential abundance methods produce different results across 38
  datasets.** *Nature Communications*, 13:342.

------------------------------------------------------------------------

# 🛠 Tools and Resources

**Programming environment**

- R v4.6.0
- renv v1.2.4
- Bash (data download and raw-read QC)

**Quality control**

- FastQC v0.12.1
- MultiQC v1.33

**Amplicon processing**

- DADA2 v1.40.0

**Taxonomic reference**

- SILVA v138.2 (non-redundant 99% training set and species assignment
  file)

**Microbiome analysis**

- phyloseq v1.56.0
- vegan (PERMANOVA, dispersion, rarefaction), ecodist (MRM)

**Statistical analysis**

- edgeR v4.10.5

**Data resources**

- NCBI Sequence Read Archive (SRA)
- European Nucleotide Archive (ENA) Portal API, used to retrieve the
  FASTQ files
- BioProject: PRJNA759642

**Visualization**

- ggplot2, ComplexHeatmap, pheatmap

------------------------------------------------------------------------

# 📌 Note

This repository is structured as a **reproducible end-to-end microbiome
bioinformatics workflow** rather than a single analysis script.

The pipeline emphasizes:

- modularity
- reproducibility
- ASV-based microbial profiling
- explicit sequencing-quality assessment
- structured microbiome data integration
- diversity analysis
- sequencing-depth sensitivity assessment
- count-based differential abundance analysis
- transparent statistical methodology
- biological interpretability
- suitability for real-world 16S rRNA microbiome projects

The workflow is designed as a **service-oriented demonstration
project**, illustrating how a complete 16S rRNA analysis can be
structured from raw sequencing data to interpretable biological results.

------------------------------------------------------------------------

# 📬 Contact

For questions, collaborations, or bioinformatics consulting inquiries:

**Davide Maccarrone**

- GitHub: <https://github.com/DaviMacca08>
- LinkedIn: <https://www.linkedin.com/in/davidemaccarrone>
- Email: <davide_maccarrone@icloud.com>
