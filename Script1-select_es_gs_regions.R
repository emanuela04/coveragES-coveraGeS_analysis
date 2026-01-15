
#_______LCR ANALYSIS PIPELINE — TWIST ES VS GS
#________Coverage metrics, LCR definition, ES–GS overlap

library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(ggplot2)

# 1. LOAD MOSDEPTH REGIONS FOR ES & GS
#
path_gs <- "./coverage/wgs_mosdepth/"  
path_es <-'./coverage/twist/results/'

read_mosdepth <- function(file) {
  df <- read_tsv(file, col_names = FALSE)
  colnames(df) <- c("chr", "start", "end", "name", "depth")
  df$length <- df$end - df$start
  df$sample <- tools::file_path_sans_ext(basename(file))
  return(df)
}

read_mosdepth <- function(file) {
  df <- read_tsv(file, col_names = FALSE)
  colnames(df) <- c("chr", "start", "end", "name", "depth")
  df$length <- df$end - df$start
  df$sample <- tools::file_path_sans_ext(basename(file))
  return(df)
}


files_es <- list.files(path_es, pattern = "*regions.bed$", recursive = TRUE, full.names = TRUE)
files_gs <- list.files(path_gs, pattern="regions.bed", full.names = TRUE)

ES <- map_df(files_es, read_mosdepth)
GS <- map_df(files_gs, read_mosdepth)
ES <- ES[ES$chr != "chrX",]
GS <- GS%>% filter(!chr %in% c("chrX", "chrY","chrM"))




paired <- read.table("samples_es_gs.txt", header=F, col.names = "sample", stringsAsFactors = F)
ES$sample <- gsub(".regions","",ES$sample)


# 2. CALCULATE MEAN COVERAGE PER SAMPLE AND SELECTE THE PAIRED ES SAMPLE WITH GS 


ES_paired <- ES %>% filter(sample %in% paired$sample)



coverage_stats_regions <- function(df) {
  df %>%
    group_by(sample) %>%
    summarise(
      n_regions = n(),
      n_regions_low20 = sum(depth <= 20),
      pct_regions_low20 = 100 * n_regions_low20 / n_regions,
      pct_regions_gt20  = 100 * (n_regions - n_regions_low20) / n_regions,
      weighted_cov = sum(depth * length) / sum(length),
      .groups = "drop"
    )
}
cov_GS=coverage_stats_regions(GS)
cov_ES=coverage_stats_regions(ES)
write.table(cov_GS, file ="stats_cov_GS", row.names=F, col.names=T, quote=F, sep="\t")
write.table(cov_ES, file ="stats_cov_ES", row.names=F, col.names=T, quote=F, sep="\t")




# 3. COMPUTE Z-SCORES FOR EACH INTERVAL


batch_info <- read.table("twist/results/batch.txt", header=F, col.names = c("sample","batch"), sep ="\t", stringsAsFactors = F) #LOAD BATCH INFO SOTRED IN SEPARED FILE
batch_info$sample<- as.character(batch_info$sample)
ES_paired<- ES_paired %>%
  left_join(batch_info, by="sample")
GS$batch <- "GS"

compute_z <- function(df) {
  df %>%
    group_by(batch, chr, start, end) %>%
    mutate(
      mean_depth = mean(depth),
      sd_depth = sd(depth),
      z = (depth - mean_depth) / sd_depth
    ) %>%
    ungroup()
}

ES_z <- compute_z(ES_paired)
GS_z <- compute_z(GS)


# 4. DEFINE LCR (depth <20)


ES_LCR <- ES_z %>% filter(depth < 20)
GS_LCR <- GS_z %>% filter(depth < 20)

write.table(ES_LCR, "ES_paired.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(GS_LCR, "GS_paired.txt", quote=FALSE, row.names=FALSE, sep="\t")

