suppressPackageStartupMessages(library(dplyr))
source('snakemake/scripts/common.R')

walk(c(1,2), function(x) {
  print(snakemake@input[[x]])
  f = read.table(snakemake@input[[x]], 
             sep = "\t", 
             header = T,
             stringsAsFactors = F) %>%
    mutate(row_sums = (.) %>% select(-one_of(featurecounts_params)) %>% rowSums()) %>%
    filter(row_sums >= snakemake@params[[x]]) %>%
    featurecounts2bed() %>%
    arrange(chrom, chromStart)
    
  write.table.bed(f, snakemake@output[[x]])
})