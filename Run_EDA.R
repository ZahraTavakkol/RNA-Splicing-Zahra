source("/home/zahrath/projects/def-hsn/zahrath/code/splicing_eda_functions.R")

# Load required libraries
library(Matrix)
library(data.table)
library(ggplot2)
library(Seurat)
library(Biostrings)
library(GenomicRanges)
library(Rsamtools)
library(stringr)
library(patchwork)

Sample_Name <- c(
  "SRR26455873", "SRR26455874", "SRR26455875",
  "SRR26455967", "SRR26455733", "SRR26456012",
  "SRR26528474", "SRR26528508",
  "SRR26528542", "SRR26528531",
  "SRR26528349", "SRR26528347",
  "SRR26455910", "SRR26455911",
  "SRR26528385", "SRR26528396", "SRR26528343"
)

# Set paths
Sample_Name <- 'SRR26455875'
INPUT_DIR <- paste0("/home/zahrath/scratch/Splicing-Pipeline/Pipeline_Output/", Sample_Name)
GENOME_FA <- "/home/zahrath/scratch/GRCm39.primary_assembly.genome.fa"
OUT_DIR <- paste0("/home/zahrath/scratch/Splicing-Pipeline/EDA_Results/", Sample_Name)
print(paste0("/home/zahrath/scratch/data/AON/", Sample_Name, "/Solo.out/Gene/filtered/"))
DATA_PATH <- paste0("/home/zahrath/scratch/data/AON/", Sample_Name, "/Solo.out/Gene/filtered/")

FLANK <- 50
KMER_SIZE <- 9

# Load data
expr <- readRDS(file.path(INPUT_DIR, "gene_expression.rds"))
m1 <- readRDS(file.path(INPUT_DIR, "m1_matrix.rds"))
event_data <- readRDS(file.path(INPUT_DIR, "m1_event_data.rds"))

# Ensure column/row names exist
colnames(expr) <- colnames(m1) <- colnames(expr)
rownames(expr) <- rownames(expr) %||% paste0("G", seq_len(nrow(expr)))
rownames(m1) <- rownames(m1) %||% event_data$row_names_mtx

# generate_summary_csv(expr, m1, event_data, OUT_DIR)

# Run EDA analyses
# expression_eda(expr, file.path(OUT_DIR, "expression_EDA"))
# expression_eda_from_path(DATA_PATH, file.path(OUT_DIR, "expression_EDA_seurat"))

# m1_eda(m1, file.path(OUT_DIR, "M1_EDA"))
# plot_junctions_per_cell_by_chr(m1, event_data, file.path(OUT_DIR, "M1_EDA"))
# plot_detected_junctions_vs_chr_length(m1, event_data, OUT_DIR)
# plot_junctions_vs_chr_length(event_data, OUT_DIR)
# rank_junctions(m1, event_data, OUT_DIR)

# Sequence motif analysis
seqs <- extract_sequences(event_data, GENOME_FA, KMER_SIZE)
# motif_summary <- analyze_motifs(seqs, KMER_SIZE)
# dir.create(file.path(OUT_DIR, "sequence_EDA"), recursive = TRUE, showWarnings = FALSE)
# write.csv(motif_summary, file.path(OUT_DIR, "sequence_EDA", "motif_summary.csv"), row.names = FALSE)
# gc_content_profile(seqs, FLANK)
# base_freq_plot(seqs, FLANK)

# Save top kmers
dir.create(file.path(OUT_DIR, "sequence_EDA"), recursive = TRUE, showWarnings = FALSE)
# top_kmers <- head(kmer_freqs(seqs, k = KMER_SIZE), 50)
# write.csv(as.data.frame(top_kmers), file.path(OUT_DIR, "sequence_EDA", "top_kmers.csv"))

top_kmers <- head(kmer_freqs_fast(seqs, k = KMER_SIZE), 50)
top_kmers_df <- data.frame(
  kmer = names(top_kmers),
  count = as.numeric(top_kmers),
  frequency = as.numeric(top_kmers) / sum(top_kmers)
)
write.csv(top_kmers_df, file.path(OUT_DIR, "sequence_EDA", "top_kmers_9.csv"), row.names = FALSE)
