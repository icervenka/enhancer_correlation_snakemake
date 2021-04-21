suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))

source('snakemake/scripts/common.R')

# path_base = "~/server/programming/enhancer_pipeline/output/"
# 
# raw=read.table(paste0(path_base, "enhancer/07_annotated/enhancers.filtered.bed"), stringsAsFactors = F)
# annotation=read.table(paste0(path_base, "enhancer/07_annotated/enhancers.annotated.bed"), stringsAsFactors = F)
# directionality=read.table(paste0(path_base, "counts/strand_enhancer_counts.txt"), stringsAsFactors = F, header = T)
# gene_counts=read.table(paste0(path_base, "counts/gene_counts.txt"), stringsAsFactors = F, header = T)
# enhancer_counts=read.table(paste0(path_base, "counts/enhancer_counts.txt"), stringsAsFactors = F, header = T)
# alignment_counts=read.table(paste0(path_base, "counts/alignment_counts.txt"), stringsAsFactors = F, header = T)
# rnaseq_metadata=read.table("~/server/programming/enhancer_pipeline/rnaseq_metadata.tsv", stringsAsFactors = F, header = T)

### load inputs and parameters from snakemake ###
annotation = read.table(snakemake@input[['enhancer_annotation']], stringsAsFactors = F)
directionality = read.table(snakemake@input[['enhancer_directionality']], stringsAsFactors = F, header = T)
gene_counts = read.table(snakemake@input[['gene_counts']], stringsAsFactors = F, header = T)
enhancer_counts = read.table(snakemake@input[['enhancer_counts']], stringsAsFactors = F, header = T)
alignment_counts = read.table(snakemake@input[['alignment_counts']], stringsAsFactors = F, header = T)

rnaseq_metadata = read.table(snakemake@params[['rnaseq_metadata']], stringsAsFactors = F, header = T)
normalization = snakemake@params[['normalize_counts']]
correlation_cutoff = as.numeric(snakemake@params[['correlation_cutoff']])

### create a mapping df between featurecounts columns and sample IDs ###
fc_sample_map = featurecount_sample_name_map(names(gene_counts) %>% unique, rnaseq_metadata) %>%
  left_join(alignment_counts, by = "sample")

### Add new columns to enhancer-gene annotation table (enhancer center, TSS, distance etc.) ###
annotation = annotation %>%
  setNames(expand.grid(bed_header, c("enh", "gene")) %>% mutate(h = paste(Var1, Var2, sep = "_")) %>% pull(h)) %>%
  mutate(center_enh = floor((chromEnd_enh - chromStart_enh)/2 + chromStart_enh)) %>%
  mutate(TSS = ifelse(strand_gene == "+", chromStart_gene, chromEnd_gene)) %>% # check if strandedness is assigned correctly
  mutate(position_enh = case_when(center_enh - TSS < 0 & strand_gene == "+" ~ "upstream",
                                  center_enh - TSS < 0 & strand_gene == "-" ~ "downstream",
                                  center_enh - TSS >= 0 & strand_gene == "+" ~ "downstream",
                                  center_enh - TSS >= 0 & strand_gene == "-" ~ "upstream")) %>%
  left_join(directionality %>% select(id, direction_enh = direction),
            by = c("name_enh" = "id")) %>%
  mutate(distance = abs(TSS - center_enh)) %>%  
  group_by(name_gene) %>% 
  arrange(distance) %>% 
  mutate(distance_rank = row_number()) %>%
  select(matches("_enh"), matches("_gene"), everything()) # reorder columns

### Combine gene-enhancer combinations with enhancer and gene expression data  ###
expression_df = annotation %>%
  select(name_enh, name_gene) %>%
  left_join(enhancer_counts %>%
              normalize_counts(fc_sample_map, normalization = normalization) %>%
              rownames_to_column(var = "name_enh") %>%
              pivot_longer(-name_enh, names_to = "sample", values_to = "value") %>%
              select(name_enh, sample, value_enh = value), by = "name_enh") %>%
  left_join(gene_counts %>%
              normalize_counts(fc_sample_map, normalization = normalization) %>%
              rownames_to_column(var = "name_gene") %>%
              pivot_longer(-name_gene, names_to = "sample", values_to = "value") %>%
              select(name_gene, sample, value_gene = value), by = c("name_gene", "sample")) %>% 
  left_join(rnaseq_metadata %>% select(sample, group), by = "sample") %>%
  select(name_enh, value_enh, name_gene, value_gene, sample, group) # reorder columns

### Create annotation file  ###
write.table(annotation,
            snakemake@output[[1]],
            quote = F,
            col.names = T,
            sep = "\t",
            row.names = F)

### Export enhancer-gene expression output file  ###
write.table(expression_df,
            snakemake@output[[2]],
            quote = F,
            col.names = T,
            sep = "\t",
            row.names = F)
