setwd('/home/saberi/projects/csg/splicing-group/RNA-Splicing-Zahra')

library(data.table)
library(stringr)
library(magrittr)
library(purrr)
library(ggplot2)

graph_csv <- 'SRR26455873/junction_graph_edges.csv'
graph <- fread(graph_csv)
graph <- graph[start_node %like% 'chr']
graph[, chr := str_extract(start_node, 'chr.*[:]') %>% str_remove('[:]')]
graph[, chr := factor(chr, levels = str_c('chr', c(1:22, 'X', 'Y', 'M')), ordered = TRUE)]
setorder(graph, chr)

thresh <- graph[, quantile(weight, 0.95)]
graph[weight > thresh, .(total_reads = sum(weight)), chr][order(-total_reads)]
graph <- graph[weight < thresh]

stats <- graph[, .(num_junctions = .N, num_reads = sum(weight)), chr]
stats[, chr := factor(chr, levels = str_c('chr', c(1:22, 'X', 'Y', 'M')), ordered = TRUE)]
setorder(stats, chr)

stats[, frac_reads := (num_reads / sum(num_reads) * 100) %>% round(2)]
stats[, reads_per_junction := num_reads / num_junctions]
stats[, frac_reads_per_junction := (reads_per_junction / sum(reads_per_junction) * 100) %>% round(2)]

stats

ggplot(stats, aes(x = chr, y = num_junctions)) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = mean(stats$num_junctions), color = 'red', linetype = 'dashed', linewidth = 1) +
  theme_minimal() +
  labs(x = 'Chromosome', y = 'Number of Junctions')

ggsave('plots/chr_num_junctions.png', width = 10, height = 5)

ggplot(stats, aes(x = chr, y = frac_reads_per_junction)) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = mean(stats$frac_reads_per_junction), color = 'red', linetype = 'dashed', linewidth = 1) +
  theme_minimal() +
  labs(x = 'Chromosome', y = 'Fraction of Reads per Junction')

ggsave('plots/chr_frac_reads_per_junction.png', width = 10, height = 5)

pdf('plots/weight_hist.pdf')
graph[, hist(weight)]
dev.off()

ggplot(graph) +
    aes(x = weight) +
    geom_density() +
    theme_minimal() +
    labs(x = 'Weight', y = 'Density')

ggsave('plots/chr_weight.png', width = 10, height = 5)