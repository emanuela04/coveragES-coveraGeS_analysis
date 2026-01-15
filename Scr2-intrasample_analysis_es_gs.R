suppressPackageStartupMessages({
  library(optparse)
  library(GenomicRanges)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(car)
  library(DescTools)
  library(pheatmap)
})

#Rscript analysis_LCR.R \
#--es /path/pairwise-ES.txt \
#--gs /path/pairwise-GS.txt \
#--outdir /path/output \
#--prefix LCR \
#--dpi 600

option_list <- list(
  make_option(c("--es"), type="character", help="Path to pairwise-ES.txt (tab-separated)", metavar="FILE"),
  make_option(c("--gs"), type="character", help="Path to pairwise-GS.txt (tab-separated)", metavar="FILE"),
  make_option(c("--omim"), type="character", default=NULL, help="Path to OMIM file (tsv)", metavar="FILE"),
  make_option(c("--sysndd"), type="character", default=NULL, help="Path to SysNDD file (tsv)", metavar="FILE"),
  make_option(c("--outdir"), type="character", default=".", help="Output directory [default: %default]", metavar="DIR"),
  make_option(c("--prefix"), type="character", default="results", help="Output prefix [default: %default]", metavar="STR"),
  make_option(c("--dpi"), type="integer", default=600, help="DPI for saved figures [default: %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Controlli minimi
if (is.null(opt$es) || is.null(opt$gs)) {
  stop("Devi fornire almeno --es e --gs. Esempio: Rscript analysis_LCR.R --es pairwise-ES.txt --gs pairwise-GS.txt")
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

### load linput bed file selected from mosdepth
es <- read.table(opt$es, header=TRUE, sep="\t", stringsAsFactors=FALSE)
gs <- read.table(opt$gs, header=TRUE, sep="\t", stringsAsFactors=FALSE)


# Create GRanges for overlap

GS_gr <- makeGRangesFromDataFrame(gs,
                                  seqnames.field = "chr",
                                  start.field    = "start",
                                  end.field      = "end",
                                  keep.extra.columns = TRUE)

ES_gr <- makeGRangesFromDataFrame(es,
                                  seqnames.field = "chr",
                                  start.field    = "start",
                                  end.field      = "end",
                                  keep.extra.columns = TRUE)
#Compute overlap per sample


all_samples <- intersect(unique(gs$sample), unique(es$sample))

results <- lapply(all_samples, function(s) {
  gs_s <- GS_gr[GS_gr$sample == s]
  es_s <- ES_gr[ES_gr$sample == s]
  
  ov   <- findOverlaps(gs_s, es_s)
  
  data.frame(
    sample = s,
    GS_low = length(gs_s),
    ES_low = length(es_s),
    overlap = length(ov),
    only_GS = length(setdiff(gs_s, es_s)),
    only_ES = length(setdiff(es_s, gs_s))
  )
})

results <- bind_rows(results)
results$sample <- as.character(results$sample)


df_melt <- results %>%
  pivot_longer(cols = c(overlap, only_GS, only_ES),
               names_to = "category",
               values_to = "count")


#  Plot regions per sample


ggplot(df_melt, aes(sample, count, fill = category)) +
  geom_col(position = "dodge") +
  theme_bw() 

#  Paired test ES vs GS per sample
results$diff_ES_GS <- results$only_ES - results$only_GS

wilcox_paired <- wilcox.test(
  results$only_ES,
  results$only_GS,
  paired = TRUE,
  alternative = "greater"
)

print(wilcox_paired)
#V = 86, p-value = 0.01764 --> ES has more LCRs compared to GS


#  Fraction test: ES-specific > 50% 
results$ES_fraction_specific <- results$only_ES / (results$only_ES + results$only_GS)

wilcox_frac <- wilcox.test(
  results$ES_fraction_specific,
  mu = 0.5,
  alternative = "greater"
)

print(wilcox_frac)

# Levene test 
# Group = platform
platform <- factor(rep(c("ES","GS"), each=nrow(results)))
values   <- c(results$only_ES, results$only_GS)

leveneTest(values ~ platform)


#  Correlation ES-only vs GS-only
# the number of ES-only LCRs per individual was not correlated with the number of 
#GS-only LCRs (Spearman ρ = 0.08; Figure S3), indicating that samples with 
#high counts of ES-specific failures did not necessarily exhibit correspondingly 
#high counts of GS-specific failures.

cor.test(results$only_ES, results$only_GS, method = "spearman")


#jaccard per sample + concordance test: To what extent do ES and GS fail in the same regions?
#Low-coverage regions do not completely overlap, so part of the LCRs is approach-specific rather than intrinsic to the genome.

jaccard <- sapply(all_samples, function(s) {
  gs_s <- GS_gr[GS_gr$sample == s]
  es_s <- ES_gr[ES_gr$sample == s]
  #numero di intervalli LCR condivisi tra ES e GS nello stesso individuo
  inter <- length(intersect(gs_s, es_s))
  #numero totale di intervalli LCR osservati in almeno uno dei due approcci
  uni   <- length(union(gs_s, es_s))
  
  if (uni == 0) return(NA)
  inter / uni
})

results$jaccard <- jaccard

# Test: jaccard << 1
wilcox.test(results$jaccard, mu = 1, alternative = "less")


ggplot(results, aes(sample, jaccard)) +
  geom_col(fill="steelblue") +
  theme_bw() +
  labs(title="Jaccard similarity ES–GS (per sample)")
df_melt$sample_id <- factor(df_melt$sample,
                            labels = paste0("Sample ", seq_along(unique(df_melt$sample))))

df_melt$sample_index <- factor(as.numeric(factor(df_melt$sample)))

mm<- ggplot(df_melt, aes(sample_index, count, fill = category)) +
  geom_col(position = "dodge") +
  scale_fill_manual(
    values = c("only_ES" = "#1f78b4",
               "only_GS" = "#33a02c",
               "overlap" = "gold"),
    labels = c("ES-only", "GS-only", "Shared")
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank()
  ) +
  labs(
    x = "Samples (n=14)",
    y = "Number of LCRs",
    title = "ES vs GS LCRs per sample"
  )


ggsave(
  filename = "Figure_S1.png",
  plot = mm,
  width = 10,
  height = 6,
  units = "in",
  dpi = 600
)


# PCA + HEATMAP OF LOW-COVERAGE REGIONS (ES vs GS)


library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)

ES_gr$sample <- as.character(ES_gr$sample)
GS_gr$sample <- as.character(GS_gr$sample)

samples <- sort(intersect(unique(ES_gr$sample), unique(GS_gr$sample)))

all_intervals <- reduce(c(ES_gr, GS_gr))

#  SAMPLE × INTERVAL MATRIX

mat <- matrix(0,
              nrow = length(samples),
              ncol = length(all_intervals),
              dimnames = list(samples, paste0("interval_", seq_len(length(all_intervals)))))

for (s in samples) {
  message("sample: ", s)
  
  es_s <- ES_gr[ES_gr$sample == s]
  gs_s <- GS_gr[GS_gr$sample == s]
  present <- reduce(c(es_s, gs_s))
  
  # Find overlaps between all_intervals and sample intervals
  hits <- findOverlaps(all_intervals, present)
  mat[s, unique(queryHits(hits))] <- 1
}





# PCA con prcomp che crea la decomposizione ai valori singolari per fare la pca
pca <- prcomp(mat, center = TRUE, scale. = FALSE)

pca_df <- data.frame(
  sample = rownames(mat),
  PC1 = pca$x[,1],
  PC2 = pca$x[,2]
)


library(pheatmap)


ES_gr$sample <- as.character(ES_gr$sample)
GS_gr$sample <- as.character(GS_gr$sample)

samples <- sort(intersect(unique(ES_gr$sample), unique(GS_gr$sample)))

# Universe of all genomic intervals seen in ES or GS
all_intervals <- reduce(c(ES_gr, GS_gr))

#I profili di low coverage sono simili o diversi in modo sistematico?similarità complessiva dei profili binari (vettori di migliaia di intervalli)
#ES e GS producono profili di low coverage distinti o simili, considerando l’intero genoma target?

profiles <- list()
rownames_list <- c()

for (s in samples) {
  message("Processing sample ", s)
  
  # ES profile
  es_s <- ES_gr[ES_gr$sample == s]
  hits_es <- findOverlaps(all_intervals, es_s)
  vec_es <- rep(0, length(all_intervals))
  vec_es[unique(queryHits(hits_es))] <- 1
  
  # GS profile
  gs_s <- GS_gr[GS_gr$sample == s]
  hits_gs <- findOverlaps(all_intervals, gs_s)
  vec_gs <- rep(0, length(all_intervals))
  vec_gs[unique(queryHits(hits_gs))] <- 1
  
  profiles[[paste0(s, "_ES")]] <- vec_es
  profiles[[paste0(s, "_GS")]] <- vec_gs
}

mat2 <- do.call(rbind, profiles)
colnames(mat2) <- paste0("interval_", seq_len(ncol(mat2)))

# Metadata for PCA plotting
meta <- data.frame(
  sample = gsub("_(ES|GS)", "", rownames(mat2)),
  platform = ifelse(grepl("ES$", rownames(mat2)), "ES", "GS")
)


pca2 <- prcomp(mat2, center = TRUE, scale. = FALSE)

pca_df <- cbind(
  rownames = rownames(mat2),
  meta,
  PC1 = pca2$x[,1],
  PC2 = pca2$x[,2]
)


#  PCA PLOT 

ggplot(pca_df, aes(PC1, PC2, color = platform, label = sample)) +
  geom_point(size = 3) +
  #geom_text(vjust = -0.7, size = 3) +
  scale_color_manual(values = c("ES" = "#1f78b4", "GS" = "#e31a1c")) +
  theme_bw(base_size = 14) +
  labs(
    title = "PCA of low-coverage profiles (ES vs GS separate)",
    x = paste0("PC1 (", round(summary(pca2)$importance[2,1]*100,1), "%)"),
    y = paste0("PC2 (", round(summary(pca2)$importance[2,2]*100,1), "%)")
  )


#  HEATMAP ES vs GS 
ann <- data.frame(platform = meta$platform)
rownames(ann) <- rownames(mat2)

ann_colors <- list(
  platform = c(
    ES = "#1F77B4",   
    GS = "#D62728"    
  )
)

pheatmap(mat2,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_row = ann,
         annotation_colors = ann_colors,   
         color = c("white", "#1f78b4"),   
         main = "Heatmap of low-coverage regions (ES vs GS separated)")



# Prepare a data frame that pairs ES and GS points
pca_df$pair_id <- pca_df$sample

# Create segment data for connecting ES–GS pairs
segments_df <- pca_df %>%
  tidyr::pivot_wider(
    id_cols = pair_id,
    names_from = platform,
    values_from = c(PC1, PC2)
  ) %>%
  na.omit()

library(ggplot2)

ggplot() +
  # Lines connecting ES–GS pairs
  geom_segment(data = segments_df,
               aes(x = PC1_ES, y = PC2_ES,
                   xend = PC1_GS, yend = PC2_GS),
               color = "grey5", linewidth = 0.7, alpha = 0.8) +
  
  # Points (ES / GS)
  geom_point(data = pca_df,
             aes(PC1, PC2, color = platform),
             size = 3) +
  
  scale_color_manual(values = c("ES" = "#1f78b4", "GS" = "#e31a1c")) +
  theme_bw(base_size = 14) +
  labs(
    title = "PCA of LCRs profiles (paired ES–GS)",
    x = paste0("PC1 (", round(summary(pca2)$importance[2,1]*100,1), "%)"),
    y = paste0("PC2 (", round(summary(pca2)$importance[2,2]*100,1), "%)"),
    color = "Sequencing methods"
  )


p <- ggplot() +
  geom_segment(
    data = segments_df,
    aes(x = PC1_ES, y = PC2_ES,
        xend = PC1_GS, yend = PC2_GS),
    color = "grey20", linewidth = 0.5, alpha = 0.8
  ) +
  geom_point(
    data = pca_df,
    aes(PC1, PC2, color = platform),
    size = 3.5
  ) +
  scale_color_manual(values = c(
    "ES" = "#1f78b4",
    "GS" = "#e31a1c"
  )) +
  theme_bw(base_size = 14) +
  labs(
    title = "PCA of LCR profiles (paired ES–GS)",
    x = paste0("PC1 (", round(summary(pca2)$importance[2,1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca2)$importance[2,2] * 100, 1), "%)"),
    color = "Sequencing method"
  )


#################

omim <- read.csv(opt$omim, sep="\t", stringsAsFactors=FALSE, quote="", fill=TRUE)
SYS_NDD <- read.csv(opt$sysndd, sep="\t", stringsAsFactors=FALSE, quote="", fill=TRUE)

omim_ndd <- omim %>% full_join(SYS_NDD, by = "SYMBOL")

db <- omim_ndd %>%
  dplyr::select(SYMBOL, Phenotypes) %>%
  left_join(
    SYS_NDD %>% dplyr::select(SYMBOL, disease_ontology_name, disease_ontology_id_version, ndd_phenotype, category),
    by = "SYMBOL"
  )

db_flags <- db %>%
  transmute(
    gene = SYMBOL,
    is_OMIM = !is.na(Phenotypes) & Phenotypes != "",
    is_NDD  = (ndd_phenotype == 1) & category %in% c("Definitive", "Moderate")
  ) %>%
  distinct(gene, .keep_all = TRUE)

TOT_NDD <- sum(db_flags$is_NDD, na.rm = TRUE)









