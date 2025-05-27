# Load required function
source("/home/zahrath/projects/def-hsn/zahrath/Arsham-s-Splicing-Pipeline-main/R/gedi_input_matrices.R")

# Define full list of sample IDs manually (from output of find)
Sample_Name <- c(
  "SRR26455873", "SRR26455874", "SRR26455875",
  "SRR26455967", "SRR26455733", "SRR26456012",
  "SRR26528474", "SRR26528508",
  "SRR26528542", "SRR26528531",
  "SRR26528349", "SRR26528347",
  "SRR26455910", "SRR26455911",
  "SRR26528385", "SRR26528396", "SRR26528343"
)

# Define base input/output paths
Data_dir <- "/scratch/zahrath/data/"
Result_dir <- "/home/zahrath/scratch/Splicing-Pipeline/Pipeline_Output/"

# Loop over samples
for (sample in Sample_Name) {
  message("Processing sample: ", sample)

  # Detect region (e.g., AON, MOp, etc.)
  region_path <- system(paste0("find ", Data_dir, " -type d -name ", sample), intern = TRUE)
  region <- basename(dirname(region_path))

  # Build input paths
  sj_path <- file.path(Data_dir, region, sample, "Solo.out", "SJ")
  gene_path <- file.path(Data_dir, region, sample, "Solo.out", "Gene")
  out_path <- file.path(Result_dir, sample)
  dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

  # Run pipeline
  SJ_object <- multigedi_make_junction_ab(
    STARsolo_SJ_dirs = sj_path,
    sample_ids = sample,
    use_internal_whitelist = TRUE
  )

  m1_obj <- multigedi_make_m1(junction_ab_object = SJ_object)
  m1 <- m1_obj$m1_inclusion_matrix

  saveRDS(m1, file.path(out_path, "m1_matrix.rds"))
  saveRDS(m1_obj$event_data, file.path(out_path, "m1_event_data.rds"))

  multigedi_countsplit(
    m1_inclusion_matrix = m1,
    folds = 2,
    epsilon = c(0.5, 0.5),
    object_names = "m1"
  )

  m2_test <- multigedi_make_m2(m1_inclusion_matrix = m1_test, eventdata = m1_obj$event_data)
  m2_train <- multigedi_make_m2(m1_inclusion_matrix = m1_train, eventdata = m1_obj$event_data)

  saveRDS(m2_test, file.path(out_path, "m2_test.rds"))
  saveRDS(m2_train, file.path(out_path, "m2_train.rds"))

  gene_expression <- multigedi_make_gene(
    expression_dirs = gene_path,
    sample_ids = sample,
    use_internal_whitelist = TRUE
  )

  saveRDS(gene_expression, file.path(out_path, "gene_expression.rds"))
}
