suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))

source('snakemake/scripts/common.R')

### load inputs and parameters from snakemake ###
gene_counts = read.table(snakemake@input[['gene_counts']], stringsAsFactors = F, header = T)
enhancer_counts = read.table(snakemake@input[['enhancer_counts']], stringsAsFactors = F, header = T)
annotation_df = read.table(snakemake@input[['annotation']], stringsAsFactors = F, header = T)
expression_df = read.table(snakemake@input[['expression']], stringsAsFactors = F, header = T)

rnaseq_metadata = read.table(snakemake@params[['rnaseq_metadata']], stringsAsFactors = F, header = T)
outdir = snakemake@params[['outdir']]
contrast = snakemake@params[['contrast']]
reference_levels = snakemake@params[['reference_levels']]
lfc_shrink = snakemake@params[['lfc_shrink']]
fdr = snakemake@params[['fdr']]
correlation_cutoff = snakemake@params[['enhancer_gene_correlation_cutoff']]

### calculate differential expression for genes and enhancers ###
gene_diffexp = deseq_diffexp(gene_counts, 
                             rnaseq_metadata, 
                             ref_levels = reference_levels, 
                             fdr = fdr)

enhancer_diffexp = deseq_diffexp(enhancer_counts, 
                                 rnaseq_metadata, 
                                 ref_levels = reference_levels, 
                                 fdr = fdr)

# add enhancer diffexp to annotation
ec = map(enhancer_diffexp, function(x, annotation) {
  annotation %>%
    left_join(x %>%
                as.data.frame() %>%
                rownames_to_column(var = "name_enh"),
              by = "name_enh")
}, annotation = annotation_df)

# add gene diffexp to annotation
gec = map2(ec, gene_diffexp, function(x, y) {
  y = y %>%
    as.data.frame() %>%
    rownames_to_column(var = "name_gene")
  x %>%
    left_join(y, by = "name_gene", suffix = c("_enh", "_gene")) %>% 
    drop_na()
})

# combine into large data frame for easier filtering
# calculate contrast-wise correlation
contrast_correlation = map_dfr(names(gec), function(x, df, expression_df, reference_levels) {
  use_groups = create_ref_level_df(reference_levels) %>% 
    filter(name == x) %>% 
    select(V1, V2) %>% unlist() %>% 
    unname
  
  contrast_corr_df = expression_df %>% 
    filter(group %in% use_groups) %>%
    group_by(name_enh, name_gene) %>%
    summarise(contrast_correlation = cor(value_enh, value_gene),
           contrast_mean_expression_enh = mean(value_enh),
           contrast_mean_expression_gene = mean(value_gene),
           .groups = "drop") %>%
    drop_na()
  
  df[[x]] %>%
    mutate(contrast = x) %>%
    left_join(contrast_corr_df, by = c("name_enh", "name_gene"))
}, df = gec, expression_df = expression_df, reference_levels = reference_levels)

# calculate global correlation
global_correlation = expression_df %>% group_by(name_enh, name_gene) %>%
  dplyr::summarise(global_correlation = cor(value_enh, value_gene),
                   global_mean_expression_enh = mean(value_enh),
                   global_mean_expression_gene = mean(value_gene),
                   .groups = "drop") %>%
  drop_na()

# combine all together
all_res = contrast_correlation %>%
  left_join(global_correlation, by = c("name_enh", "name_gene"))

# organize column order
all_res = all_res %>%
  select(contrast, name_enh, name_gene, TSS, center_enh, position_enh, direction_enh, distance, distance_rank,
         contrast_correlation, contrast_mean_expression_enh, contrast_mean_expression_gene,
         global_correlation, global_mean_expression_enh, global_mean_expression_gene,
         matches("_enh"), matches("_gene"))

# filtering
all_res_filtered = all_res %>%
  filter(padj_enh <= fdr) %>%
  filter(padj_gene <= fdr) %>%
  filter(abs(global_correlation) >= correlation_cutoff | 
           abs(contrast_correlation) >= correlation_cutoff)

# export all results in a single file
write.table(all_res,
            snakemake@output[[1]],
            sep = "\t", 
            quote = F, 
            row.names = F)
            
# export filtered results separated into contrasts
all_res_filtered %>%
  group_by(contrast) %>%
  group_walk(~ write.table(.x, 
                           paste0(outdir, "corr_filtered_", .y$contrast, ".txt"), 
                           sep = "\t", 
                           quote = F, 
                           row.names = F))
