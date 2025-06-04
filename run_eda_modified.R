source("/home/zahrath/projects/def-hsn/zahrath/functions_modified.R")

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

DATA_BASE <- "/scratch/zahrath/data/"
RESULT_DIR <- "/home/zahrath/scratch/Splicing-Pipeline/Pipeline_Output"
genome_fa <- "/home/zahrath/scratch/GRCm39.primary_assembly.genome.fa"

# Create output directory
GLOBAL_OUT_DIR <- "/home/zahrath/scratch/Splicing-Pipeline/EDA_Results_Mod"
dir.create(GLOBAL_OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Define parameters
FLANK <- 50
KMER_SIZE <- 6

for (i in seq_along(Sample_Name)) {
    print(Sample_Name[i])
    sample <- Sample_Name[i]
    region <- Regions[i]

    INPUT_DIR <- file.path(RESULT_DIR, sample)
    OUT_DIR <- file.path(GLOBAL_OUT_DIR, sample)
    DATA_PATH <- file.path(DATA_BASE, region, sample, "Solo.out/Gene/filtered/")
    dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

    # Load data
    expr <- readRDS(file.path(INPUT_DIR, "gene_expression.rds"))
    m1 <- readRDS(file.path(INPUT_DIR, "m1_matrix.rds"))
    event_data <- readRDS(file.path(INPUT_DIR, "m1_event_data.rds"))

    colnames(expr) <- colnames(m1) <- colnames(expr)
    rownames(expr) <- rownames(expr) %||% paste0("G", seq_len(nrow(expr)))
    rownames(m1) <- rownames(m1) %||% event_data$row_names_mtx

  # ---- Summary CSV ----
    generate_summary_csv(expr, m1, event_data, OUT_DIR)

  # ---- Expression EDA ----
    expression_eda(expr, OUT_DIR)

  # ---- M1 EDA ----
    m1_eda(m1, OUT_DIR)

  # ---- Motif & GC ----

  seqs_list <- extract_donor_acceptor_sequences(event_data, genome_fa, flank = FLANK)

  kmer_freqs(seqs = seqs_list$donor, out_path = OUT_DIR, k = KMER_SIZE, label = "donor")
  kmer_freqs(seqs = seqs_list$acceptor, out_path = OUT_DIR, k = KMER_SIZE, label = "acceptor")

  # ---- Junction Graph ----
  make_junction_graph(m1, OUT_DIR)

  # ---- Junction Stats ----
  rank_junctions(m1, event_data, OUT_DIR, top_n = 20)
  plot_junction_length_per_chr(event_data, OUT_DIR)
  plot_junction_count_per_chr(m1, OUT_DIR)


  gtf_file <- "/home/zahrath/scratch/Mus_musculus.GRCm39.110.gtf"
  txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
  gene_gr <- genes(txdb)
  gene_chr_map <- data.frame(
    gene = names(gene_gr),  # Ensembl IDs
    chr = as.character(seqnames(gene_gr)),
    stringsAsFactors = FALSE
  )

  expr_mat <- ReadMtx(
    mtx = file.path(DATA_PATH, "matrix.mtx"),
    features = file.path(DATA_PATH, "features.tsv"),
    cells = file.path(DATA_PATH, "barcodes.tsv"),
    feature.column = 1
  )

  rownames(expr_mat) <- sub("\\..*", "", rownames(expr_mat))

  print(row.names(expr_mat)[1:10])

  gene_chr_vec <- gene_chr_map$chr[match(rownames(expr_mat), gene_chr_map$gene)]

  plot_std_boxplot_per_chr(expr_mat, gene_chr_vec, OUT_DIR, label = "Gene")
  plot_std_boxplot_per_chr(m1, event_data$chr, OUT_DIR, label = "Junction")

}
