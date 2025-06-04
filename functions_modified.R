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
library(igraph)
library(reshape2)
library(GenomicFeatures)

# -------------------- Overal Summary --------------------

generate_summary_csv <- function(expr, m1, event_data, out_path) {
  dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
  
  cell_count <- ncol(expr)
  gene_count <- nrow(expr)
  junction_count <- nrow(m1)
  junctions_per_chr <- table(event_data$chr)

  summary_df <- data.frame(
    metric = c("cells", "genes", "junctions", paste0("junctions_chr_", names(junctions_per_chr))),
    value = c(cell_count, gene_count, junction_count, as.numeric(junctions_per_chr))
  )
  write.csv(summary_df, file.path(out_path, "summary.csv"), row.names = FALSE)
}

make_junction_graph <- function(m1, out_path, min_weight = 1) {
  dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

  # Parse junction name cleanly without suffixes
  parse_coords <- function(jn) {
    jn_clean <- sub("_[SE]$", "", jn)  # Remove _S or _E
    parts <- unlist(strsplit(jn_clean, "[:-]"))
    if (length(parts) < 3) return(NULL)
    strand_info <- ifelse(grepl("_S$", jn), "S", ifelse(grepl("_E$", jn), "E", NA))
    if (is.na(strand_info)) return(NULL)
    list(chr = parts[1], start = as.integer(parts[2]), end = as.integer(parts[3]), dir = strand_info)
  }

  edge_list <- lapply(rownames(m1), function(jn) {
    parsed <- parse_coords(jn)
    if (is.null(parsed)) {
      message("Skipping invalid junction format: ", jn)
      return(NULL)
    }
    from <- paste0(parsed$chr, ":", ifelse(parsed$dir == "S", parsed$start, parsed$end))
    to   <- paste0(parsed$chr, ":", ifelse(parsed$dir == "S", parsed$end, parsed$start))
    weight <- sum(m1[jn, ])
    if (weight < min_weight) return(NULL)
    data.frame(from = from, to = to, weight = weight, stringsAsFactors = FALSE)
  })

  edge_df <- do.call(rbind, edge_list)
  if (nrow(edge_df) == 0) stop("No valid junctions found for graph")

  library(igraph)
  g <- graph_from_data_frame(edge_df, directed = TRUE)

  V(g)$degree <- degree(g)
  E(g)$weight <- edge_df$weight

  edge_table <- data.frame(
    start_node = edge_df$from,
    end_node = edge_df$to,
    weight = edge_df$weight,
    start_node_degree = V(g)[edge_df$from]$degree,
    end_node_degree = V(g)[edge_df$to]$degree
  )
  write.csv(edge_table, file.path(out_path, "junction_graph_edges.csv"), row.names = FALSE)

  png(file.path(out_path, "junction_graph_sparse.png"), width = 1600, height = 1200)
  plot(g,
       layout = layout_with_fr(g),
       vertex.size = V(g)$degree * 1.5,
       edge.width = E(g)$weight / max(E(g)$weight) * 5,
       vertex.label = NA,
       edge.arrow.size = 0.3,
       vertex.color = "orange",
       main = "Junction Graph (Sparse Layout)")
  dev.off()
}


# -------------------- Updated extract_sequences --------------------
extract_donor_acceptor_sequences <- function(event_df, genome_fa, flank = 50, out_path = "sequence_EDA") {
  dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

  message("Opening genome FASTA...")
  fa <- FaFile(genome_fa)
  open(fa)

  message("Preparing GRanges for donor and acceptor...")
  gr_both <- GRanges(
    seqnames = c(event_df$chr, event_df$chr),
    ranges = IRanges(
      start = c(pmax(event_df$start - flank, 1), pmax(event_df$end - flank, 1)),
      end =   c(event_df$start + flank,          event_df$end + flank)
    ),
    strand = rep(ifelse(event_df$strand == 1, "+", "-"), 2)
  )
  names(gr_both) <- c(
    paste0(event_df$chr, ":", event_df$start, "_donor"),
    paste0(event_df$chr, ":", event_df$end, "_acceptor")
  )

  message("Fetching sequences from genome...")
  all_seqs <- getSeq(fa, gr_both)
  message("Done fetching sequences.")
  close(fa)
  message("Genome FASTA closed.")

  # Split into donor and acceptor
  donor_seqs <- all_seqs[1:nrow(event_df)]
  acceptor_seqs <- all_seqs[(nrow(event_df) + 1):(2 * nrow(event_df))]

  # # Save to CSV
  # write.csv(
  #   data.frame(name = names(donor_seqs), sequence = as.character(donor_seqs)),
  #   file = file.path(out_path, "donor_sequences.csv"),
  #   row.names = FALSE
  # )
  # write.csv(
  #   data.frame(name = names(acceptor_seqs), sequence = as.character(acceptor_seqs)),
  #   file = file.path(out_path, "acceptor_sequences.csv"),
  #   row.names = FALSE
  # )

  return(list(donor = donor_seqs, acceptor = acceptor_seqs))
}

# -------------------- Expression EDA --------------------
expression_eda <- function(expr_mat, outdir) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  gene_counts <- Matrix::rowSums(expr_mat > 0)
  cell_counts <- Matrix::colSums(expr_mat > 0)

  png(file.path(outdir, "gene_cell_hist.png"))
  hist(gene_counts, breaks = 50, main = "Cells per Gene", xlab = "# Cells")
  dev.off()

  png(file.path(outdir, "cell_gene_hist.png"))
  hist(cell_counts, breaks = 50, main = "Genes per Cell", xlab = "# Genes")
  dev.off()

  expr_filtered <- expr_mat[gene_counts > 0, cell_counts > 0]
  gene_counts <- Matrix::rowSums(expr_filtered > 0)
  cell_counts <- Matrix::colSums(expr_filtered > 0)

  png(file.path(outdir, "gene_cell_hist_filtered.png"))
  hist(gene_counts, breaks = 50, main = "Cells per Gene", xlab = "# Cells")
  dev.off()

  png(file.path(outdir, "cell_gene_hist_filtered.png"))
  hist(cell_counts, breaks = 50, main = "Genes per Cell", xlab = "# Genes")
  dev.off()

  seurat_obj <- CreateSeuratObject(counts = expr_mat)
  
  mito_genes <- grep("^MT-", rownames(seurat_obj), value = TRUE)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, features = mito_genes)

  seurat_obj$depthLog10 <- log10(seurat_obj$nCount_RNA + 1)

  png(file.path(outdir, "box_UMI_counts.png"))
  boxplot(seurat_obj$nCount_RNA, main = "UMI Counts per Cell", ylab = "UMIs")
  dev.off()

  png(file.path(outdir, "box_gene_counts.png"))
  boxplot(seurat_obj$nFeature_RNA, main = "Genes Detected per Cell", ylab = "# Genes")
  dev.off()

  png(file.path(outdir, "box_mito_percent.png"))
  boxplot(seurat_obj$percent.mt, main = "Mitochondrial Content", ylab = "% Mito")
  dev.off()

  seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.2, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:20, verbose = FALSE)

  png(file.path(outdir, "umap_clusters.png"))
  print(DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters") + ggtitle("UMAP - Clusters"))
  dev.off()

  png(file.path(outdir, "umap_percent_mt.png"))
  print(FeaturePlot(seurat_obj, features = "percent.mt") + ggtitle("UMAP - % Mito"))
  dev.off()

  png(file.path(outdir, "umap_depthLog10.png"))
  print(FeaturePlot(seurat_obj, features = "depthLog10") + ggtitle("UMAP - Log10 Depth"))
  dev.off()
}

# -------------------- M1 EDA --------------------
m1_eda <- function(m1_mat, outdir) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  junction_counts <- Matrix::rowSums(m1_mat > 0)
  cell_counts <- Matrix::colSums(m1_mat > 0)

  png(file.path(outdir, "junction_cell_hist.png"))
  hist(junction_counts, breaks = 50, main = "Cells per Junction", xlab = "# Cells")
  dev.off()

  png(file.path(outdir, "cell_junction_hist.png"))
  hist(cell_counts, breaks = 50, main = "Junctions per Cell", xlab = "# Junctions")
  dev.off()
}

# -------------------- rank_junctions with plot --------------------
rank_junctions <- function(m1, event_data, out_path, top_n = 20) {
  dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
  presence_count <- Matrix::rowSums(m1 > 0)
  total_values <- Matrix::rowSums(m1)

  rank_df <- data.frame(
    junction = event_data$row_names_mtx,
    presence_in_cells = presence_count,
    total_UMI = total_values
  )
  rank_df <- rank_df[order(-rank_df$presence, -rank_df$total), ]
  write.csv(rank_df, file.path(out_path, "ranked_junctions.csv"), row.names = FALSE)
}

# -------------------- GC content and base freq plots --------------------

gc_content_profile <- function(seqs, out_path, label = "", flank = 50, window = 10) {
  dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

  seq_mat <- sapply(as.character(seqs), function(s) strsplit(s, "")[[1]])
  seq_mat <- do.call(rbind, seq_mat)

  gc_profile <- sapply(1:(ncol(seq_mat) - window + 1), function(i) {
    window_seq <- seq_mat[, i:(i + window - 1)]
    rowMeans(window_seq %in% c("G", "C"))
  })
  gc_avg <- colMeans(gc_profile)

  df <- data.frame(
    position = seq(-flank, flank - window + 1),
    gc = gc_avg
  )

  png(file.path(out_path, paste0("gc_content_plot_", label, ".png")), width = 1000, height = 600)
  print(
    ggplot(df, aes(x = position, y = gc)) +
        geom_line(color = "darkgreen") +
        labs(title = paste("GC Content Meta-profile", label), x = "Position around Junction", y = "GC%")
  )
  dev.off()
}

base_freq_plot <- function(seqs, out_path, label = "", flank = 50) {
  dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

  seq_mat <- sapply(as.character(seqs), function(s) strsplit(s, "")[[1]])
  seq_mat <- do.call(rbind, seq_mat)

  base_freq <- sapply(1:ncol(seq_mat), function(i) {
    table(factor(seq_mat[, i], levels = c("A", "C", "G", "T"))) / nrow(seq_mat)
  })
  base_df <- as.data.frame(t(base_freq))
  base_df$position <- seq(-flank, flank)

  base_df_long <- reshape2::melt(base_df, id.vars = "position", variable.name = "base", value.name = "frequency")

  png(file.path(out_path, paste0("base_freq_plot_", label, ".png")), width = 1000, height = 600)
  print(
    ggplot(base_df_long, aes(x = position, y = frequency, color = base)) +
        geom_line() +
        labs(title = paste("Base Frequency Around Splice Sites", label), y = "Frequency")
  )
  dev.off()
}


# -------------------- Std deviation per chromosome (Log-scaled) --------------------
plot_std_boxplot_per_chr <- function(mat, chr_vec, out_path, label = "") {
  stopifnot(length(chr_vec) == nrow(mat))

  # Define canonical chromosomes
  canonical_chr <- paste0("chr", c(1:19, "X", "Y"))

  # Compute full SD
  std_all <- apply(mat, 1, sd)

  # Compute filtered SD: only non-zero entries per row
  std_nonzero <- apply(mat, 1, function(x) if (sum(x != 0) > 1) sd(x[x != 0]) else NA)

  df <- data.frame(
    chr = chr_vec,
    std = std_all,
    std_filtered = std_nonzero
  )

  # Ensure chr is a factor and sort it
  df$chr <- factor(df$chr, levels = sort(unique(df$chr)))

  # Plot helpers
  save_plot <- function(data, col, name_suffix, title_suffix) {
    ggsave(
      file.path(out_path, paste0("std_boxplot_", label, "_", name_suffix, ".png")),
      plot = ggplot(data, aes(x = chr, y = .data[[col]])) +
        geom_boxplot(outlier.size = 0.5, fill = "steelblue") +
        scale_y_log10() +
        labs(
          title = paste("Standard Deviation per Chromosome -", label, title_suffix),
          x = "Chromosome", y = "Standard Deviation (log10)"
        ) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
      width = 10, height = 6
    )
  }

  # 1. All chromosomes, full data
  save_plot(df, "std", "all_chr", "(all chromosomes)")

  # 2. Canonical only, full data
  save_plot(subset(df, chr %in% canonical_chr), "std", "canonical_chr", "(canonical only)")

  # 3. All chromosomes, filtered non-zero SD
  save_plot(subset(df, !is.na(std_filtered)), "std_filtered", "nonzero_all_chr", "(non-zero only)")

  # 4. Canonical only, filtered non-zero SD
  save_plot(subset(df, chr %in% canonical_chr & !is.na(std_filtered)), "std_filtered", "nonzero_canonical_chr", "(non-zero + canonical)")
}


# -------------------- Junction count per chromosome --------------------
plot_junction_count_per_chr <- function(m1, out_path) {
  chr <- sub(":.*", "", rownames(m1))
  count_df <- as.data.frame(table(chr))
  png(file.path(out_path, "junctions_per_chr_scatter.png"), width = 800, height = 500)
  print(
    ggplot(count_df, aes(x = chr, y = Freq)) +
      geom_point(size = 3, color = "darkred") +
      labs(title = "Junctions per Chromosome", x = "Chromosome", y = "# Junctions") +
      theme_minimal()
  )
  dev.off()
}

# -------------------- Junction length per chromosome --------------------
plot_junction_length_per_chr <- function(event_data, out_path) {
  dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
  event_data$length <- event_data$end - event_data$start + 1

  library(ggplot2)
  library(gtools)

  make_plot <- function(data, fname, title_suffix) {
    data$chr <- factor(data$chr, levels = mixedsort(unique(data$chr)))
    png(file.path(out_path, fname), width = 1200, height = 700)
    print(
      ggplot(data, aes(x = chr, y = length)) +
        geom_boxplot(outlier.size = 0.4, fill = "darkgreen") +
        scale_y_log10() +
        labs(title = paste0("Junction Length Distribution ", title_suffix),
             x = "Chromosome", y = "Junction Length (log10 bp)") +
        theme_minimal(base_size = 12) +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
              plot.margin = margin(10, 10, 20, 10))
    )
    dev.off()
  }

  # Plot with all chromosomes
  make_plot(event_data, "junction_length_boxplot_all_chromosomes.png", "(All Chromosomes)")

  # Plot with canonical chromosomes only
  canonical_data <- subset(event_data, grepl("^chr([0-9]{1,2}|X|Y)$", chr))
  make_plot(canonical_data, "junction_length_boxplot_canonical_chromosomes.png", "(Canonical Chromosomes)")
}


# -------------------- Kmer Frequency --------------------
kmer_freqs <- function(seqs, out_path, k = 6, top_n = 50, label = "") {
  dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

  message("Calculating k-mer frequencies...")
  freqs_mat <- oligonucleotideFrequency(seqs, width = k, step = 1)
  freqs_total <- colSums(freqs_mat)
  freqs_sorted <- sort(freqs_total, decreasing = TRUE)

  df <- data.frame(
    kmer = names(freqs_sorted),
    frequency = as.numeric(freqs_sorted)
  )

  top_df <- head(df, top_n)

  # Save CSV
  csv_file <- file.path(out_path, paste0("top_kmers_k", k, if (label != "") paste0("_", label), ".csv"))
  write.csv(top_df, csv_file, row.names = FALSE)
  message("Saved top k-mers to: ", csv_file)

  # Optional plot
  png(file.path(out_path, paste0("top_kmers_k", k, if (label != "") paste0("_", label), ".png")), width = 1000, height = 600)
  print(
    ggplot(top_df, aes(x = reorder(kmer, -frequency), y = frequency)) +
      geom_bar(stat = "identity", fill = "darkblue") +
      labs(title = paste("Top", top_n, paste0(k, "-mers"), if (label != "") paste0("(", label, ")")),
           x = "k-mer", y = "Frequency") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  )
  dev.off()

  return(df)
}