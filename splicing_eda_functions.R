
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

# -------------------- Expression EDA --------------------
expression_eda <- function(expr_mat, outdir) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  gene_counts <- Matrix::rowSums(expr_mat > 0)
  cell_counts <- Matrix::colSums(expr_mat > 0)

  # Plot histograms BEFORE filtering
  png(file.path(outdir, "gene_cell_hist.png"))
  hist(gene_counts, breaks = 50, main = "Cells per Gene", xlab = "# Cells")
  dev.off()

  png(file.path(outdir, "cell_gene_hist.png"))
  hist(cell_counts, breaks = 50, main = "Genes per Cell", xlab = "# Genes")
  dev.off()

  # Filter genes and cells
  expr_filtered <- expr_mat[gene_counts > 0, cell_counts > 0]

  if (nrow(expr_filtered) >= 2 && ncol(expr_filtered) >= 2) {
    expr_subset <- expr_filtered[, 1:min(100, ncol(expr_filtered))]
    mat <- t(as.matrix(expr_subset))

    vars <- apply(mat, 2, var)
    mat <- mat[, vars > 0, drop = FALSE]

    if (ncol(mat) >= 2) {
      expr_pca <- prcomp(mat, scale. = TRUE)
      png(file.path(outdir, "pca_expression.png"))
      plot(expr_pca$x[, 1:2], pch = 20, col = "blue", main = "PCA - Expression")
      dev.off()
    } else {
      message("Not enough non-constant features for PCA.")
    }
  } else {
    message("Not enough non-zero genes/cells for PCA.")
  }

  gene_counts <- Matrix::rowSums(expr_filtered > 0)
  cell_counts <- Matrix::colSums(expr_filtered > 0)

  # Plot histograms BEFORE filtering
  png(file.path(outdir, "gene_cell_hist_filtered.png"))
  hist(gene_counts, breaks = 50, main = "Cells per Gene", xlab = "# Cells")
  dev.off()

  png(file.path(outdir, "cell_gene_hist_filtered.png"))
  hist(cell_counts, breaks = 50, main = "Genes per Cell", xlab = "# Genes")
  dev.off()
}

expression_eda_from_path <- function(data_path, outdir, umi_min = 500) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  expr_mat <- ReadMtx(
    mtx = file.path(data_path, "matrix.mtx"),
    features = file.path(data_path, "features.tsv"),
    cells = file.path(data_path, "barcodes.tsv")
  )
  seurat_obj <- CreateSeuratObject(counts = expr_mat)
  
  mito_genes <- grep("^MT-", rownames(seurat_obj), value = TRUE)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, features = mito_genes)

  seurat_obj$depthLog10 <- log10(seurat_obj$nCount_RNA + 1)

  # seurat_obj <- subset(seurat_obj, subset = percent.mt < 12 & nCount_RNA > umi_min & nCount_RNA < 40000)

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
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE)
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

  return(seurat_obj)
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

  # PCA on cells
  if (ncol(m1_mat) >= 2) {
    m1_subset <- m1_mat[, 1:min(100, ncol(m1_mat))]
    mat <- t(as.matrix(m1_subset))
    vars <- apply(mat, 2, var)
    mat <- mat[, vars > 0, drop = FALSE]

    if (ncol(mat) >= 2) {
      m1_pca <- prcomp(mat, scale. = TRUE)
      png(file.path(outdir, "pca_splicing.png"))
      plot(m1_pca$x[, 1:2], pch = 20, col = "darkgreen", main = "PCA - Splicing Junctions")
      dev.off()
    }
  }
}

plot_junctions_per_cell_by_chr <- function(m1_mat, event_data, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  stopifnot(nrow(m1_mat) == nrow(event_data))
  chr_vec <- event_data$chr

  chr_levels <- unique(chr_vec)
  chr_levels <- chr_levels[order(chr_levels)]

  chr_counts_list <- lapply(chr_levels, function(chr) {
    idx <- which(chr_vec == chr)
    counts <- Matrix::colSums(m1_mat[idx, , drop = FALSE])
    data.frame(
      cell = names(counts),
      chromosome = chr,
      count = as.numeric(counts)
    )
  })

  junction_counts_df <- do.call(rbind, chr_counts_list)

  png(file.path(outdir, "junctions_per_cell_by_chr.png"), width = 1200, height = 600)
  print(
    ggplot(junction_counts_df, aes(x = chromosome, y = count)) +
      geom_boxplot(outlier.size = 0.5) +
      theme_bw(base_size = 12) +
      labs(title = "Splice Junctions per Cell by Chromosome", x = "Chromosome", y = "# Junctions") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
  dev.off()
}

plot_junctions_vs_chr_length <- function(event_data, out_path) {
  chr_lengths <- tapply(event_data$end - event_data$start + 1, event_data$chr, sum)
  junction_counts <- table(event_data$chr)
  junction_lengths <- tapply(event_data$end - event_data$start + 1, event_data$chr, function(x) x)

  df1 <- data.frame(chr = names(chr_lengths), length = as.numeric(chr_lengths), count = as.numeric(junction_counts))

  png(file.path(out_path, "box_junctions_vs_chr_length.png"))
  print(
    ggplot(df1, aes(x = length, y = count)) +
      geom_boxplot() +
      labs(title = "Junction Count vs. Chromosome Length", x = "Chromosome Length", y = "# Junctions")
  )
  dev.off()

  df2 <- stack(junction_lengths)
  colnames(df2) <- c("length", "chr")
  df2$chr_len <- chr_lengths[df2$chr]

  png(file.path(out_path, "box_junction_length_vs_chr_length.png"))
  print(
    ggplot(df2, aes(x = chr_len, y = length)) +
      geom_boxplot() +
      labs(title = "Junction Length vs. Chromosome Length", x = "Chromosome Length", y = "Junction Length")
  )
  dev.off()
}

plot_detected_junctions_vs_chr_length <- function(m1, event_data, out_path) {
  dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

  junction_detected <- Matrix::rowSums(m1 > 0) > 0
  detected_events <- event_data[junction_detected, ]

  junction_length <- detected_events$end - detected_events$start + 1
  chr <- detected_events$chr

  chr_junction_count <- as.data.frame(table(chr))

  chr_junction_length_list <- split(junction_length, chr)
  chr_junction_length_df <- stack(chr_junction_length_list)
  colnames(chr_junction_length_df) <- c("length", "chr")

  chr_detected_range <- tapply(detected_events$end, chr, max) - tapply(detected_events$start, chr, min)
  chr_length_df <- data.frame(chr = names(chr_detected_range), chr_length = as.numeric(chr_detected_range))

  count_df <- merge(chr_junction_count, chr_length_df, by = "chr")
  colnames(count_df) <- c("chr", "junction_count", "chr_length")

  # ---- Plot 1: Boxplot of # junctions per chromosome vs. chr length ----
  png(file.path(out_path, "box_detected_junctions_vs_chr_length.png"), width = 1000, height = 600)
  print(
    ggplot(count_df, aes(x = chr_length, y = junction_count)) +
      geom_boxplot(aes(group = cut_width(chr_length, 1e7)), outlier.size = 0.8) +
      labs(title = "Detected Junctions vs. Chromosome Length",
           x = "Chromosome Length (bp)",
           y = "# Detected Junctions") +
      theme_bw()
  )
  dev.off()

  # ---- Plot 2: Boxplot of junction length vs. chr length ----
  chr_junction_length_df$chr_length <- chr_detected_range[chr_junction_length_df$chr]

  png(file.path(out_path, "box_junction_length_vs_chr_length_detected.png"), width = 1000, height = 600)
  print(
    ggplot(chr_junction_length_df, aes(x = chr_length, y = length)) +
      geom_boxplot(aes(group = cut_width(chr_length, 1e7)), outlier.size = 0.8) +
      labs(title = "Junction Length vs. Chromosome Length (Detected)",
           x = "Chromosome Length (bp)",
           y = "Junction Length") +
      theme_bw()
  )
  dev.off()
}

rank_junctions <- function(m1, event_data, out_path) {
  dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
  presence_count <- Matrix::rowSums(m1 > 0)
  total_values <- Matrix::rowSums(m1)

  rank_df <- data.frame(
    junction = event_data$row_names_mtx,
    presence = presence_count,
    total = total_values
  )
  rank_df <- rank_df[order(-rank_df$presence, -rank_df$total), ]
  write.csv(rank_df, file.path(out_path, "ranked_junctions.csv"), row.names = FALSE)
}

# -------------------- Sequence + Motif --------------------
extract_sequences <- function(event_df, genome_fa, flank = 50) {
  fa <- FaFile(genome_fa)
  open(fa)
  gr <- GRanges(
    seqnames = event_df$chr,
    ranges = IRanges(start = pmax(event_df$start - flank, 1), end = event_df$end + flank),
    strand = ifelse(event_df$strand == 1, "+", "-")
  )
  names(gr) <- paste0(event_df$chr, ":", event_df$start, "-", event_df$end)
  seqs <- getSeq(fa, gr)
  close(fa)
  return(seqs)
}

gc_content_profile <- function(seqs, flank = 50, window = 10) {
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

  ggplot(df, aes(x = position, y = gc)) +
    geom_line(color = "darkgreen") +
    labs(title = "GC Content Meta-profile", x = "Position around Junction", y = "GC%")
}

base_freq_plot <- function(seqs, flank = 50) {
  seq_mat <- sapply(as.character(seqs), function(s) strsplit(s, "")[[1]])
  seq_mat <- do.call(rbind, seq_mat)

  base_freq <- sapply(1:ncol(seq_mat), function(i) {
    table(factor(seq_mat[, i], levels = c("A", "C", "G", "T"))) / nrow(seq_mat)
  })
  base_df <- as.data.frame(t(base_freq))
  base_df$position <- seq(-flank, flank)

  base_df_long <- reshape2::melt(base_df, id.vars = "position", variable.name = "base", value.name = "frequency")

  ggplot(base_df_long, aes(x = position, y = frequency, color = base)) +
    geom_line() +
    labs(title = "Base Frequency Around Splice Sites", y = "Frequency")
}

analyze_motifs <- function(seqs, flank = 50) {
  donor <- substr(seqs, flank + 1, flank + 2)
  acceptor <- substr(seqs, width(seqs) - flank - 1, width(seqs) - flank)
  branch_region <- subseq(seqs, start = width(seqs) - flank*2, end = width(seqs) - flank)
  ppt_region <- subseq(seqs, start = width(seqs) - flank - 10, end = width(seqs) - flank)

  data.frame(
    donor_GT_GC = donor %in% c("GT", "GC"),
    acceptor_AG = acceptor == "AG",
    has_branch_point = vcountPattern("A", branch_region) > 0,
    has_ppt = sapply(as.character(ppt_region), function(s) {
      sum(strsplit(s, "")[[1]] %in% c("C", "T")) >= 8
    })
  )
}

kmer_freqs <- function(seqs, k = 6) {
  kmers <- unlist(lapply(as.character(seqs), function(s) {
    sapply(1:(nchar(s) - k + 1), function(i) substr(s, i, i + k - 1))
  }))
  sort(table(kmers), decreasing = TRUE)
}

kmer_freqs_fast <- function(seqs, k = 6) {
  freqs_mat <- oligonucleotideFrequency(seqs, width = k, step = 1)
  freqs_total <- colSums(freqs_mat)
  freqs_sorted <- sort(freqs_total, decreasing = TRUE)
  return(freqs_sorted)
}
