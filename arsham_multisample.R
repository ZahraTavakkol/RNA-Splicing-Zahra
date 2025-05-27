# Zahra's Splicing EDA Pipeline in R (Multi-sample version)

# -------------------- Libraries --------------------
library(Matrix)
library(data.table)
library(ggplot2)
library(Seurat)
library(Biostrings)
library(GenomicRanges)
library(Rsamtools)
library(stringr)
library(patchwork)

# -------------------- Settings --------------------
DATA_BASE <- "/scratch/zahrath/data/"
RESULT_DIR <- "/home/zahrath/scratch/Splicing-Pipeline/Pipeline_Output/MutiSample/"

# -------------------- Sample Information --------------------
Sample_Name <- c(
  "SRR26455873", "SRR26455874", "SRR26455875",
  "SRR26455967", "SRR26455733", "SRR26456012",
  "SRR26528474", "SRR26528508",
  "SRR26528542", "SRR26528531",
  "SRR26528349", "SRR26528347",
  "SRR26455910", "SRR26455911",
  "SRR26528385", "SRR26528396", "SRR26528343"
)

Regions <- c(
  "AON", "AON", "AON",
  "OLFB", "OLFB", "OLFB",
  "PAL", "PAL",
  "MOp", "MOp",
  "STR-STRv", "STR-STRv",
  "ORB", "ORB",
  "striatum", "striatum", "striatum"
)

SJ_dirs <- file.path(DATA_BASE, Regions, Sample_Name, "Solo.out/SJ")
Gene_dirs <- file.path(DATA_BASE, Regions, Sample_Name, "Solo.out/Gene")

# -------------------- Load pipeline function --------------------
source("/home/zahrath/projects/def-hsn/zahrath/Arsham-s-Splicing-Pipeline-main/R/gedi_input_matrices.R")

# -------------------- Run GEDI functions on all samples --------------------
SJ_object <- multigedi_make_junction_ab(
  STARsolo_SJ_dirs = SJ_dirs,
  sample_ids = Sample_Name,
  use_internal_whitelist = TRUE
)

m1_obj <- multigedi_make_m1(junction_ab_object = SJ_object)
m1 <- m1_obj$m1_inclusion_matrix

saveRDS(m1, file.path(RESULT_DIR, "m1_matrix.rds"))
saveRDS(m1_obj$event_data, file.path(RESULT_DIR, "m1_event_data.rds"))

multigedi_countsplit(
  m1_inclusion_matrix = m1,
  folds = 2,
  epsilon = c(0.5, 0.5),
  object_names = "m1"
)

m2_test <- multigedi_make_m2(m1_inclusion_matrix = m1_test, eventdata = m1_obj$event_data)
m2_train <- multigedi_make_m2(m1_inclusion_matrix = m1_train, eventdata = m1_obj$event_data)

saveRDS(m2_test, file.path(RESULT_DIR, "m2_test.rds"))
saveRDS(m2_train, file.path(RESULT_DIR, "m2_train.rds"))

gene_expression <- multigedi_make_gene(
  expression_dirs = Gene_dirs,
  sample_ids = Sample_Name,
  use_internal_whitelist = TRUE
)

saveRDS(gene_expression, file.path(RESULT_DIR, "gene_expression.rds"))
