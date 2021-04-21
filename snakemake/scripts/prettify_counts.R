suppressPackageStartupMessages(library(tidyverse))

counts = read.table(snakemake@input[[1]], 
                    sep = "\t", 
                    header = T, 
                    stringsAsFactors = F)

counts = counts %>%
  dplyr::select(-c(Chr, Start, End, Strand, Length))

header_samples = data.frame(colname = names(counts)[-1], stringsAsFactors = F) %>%
  extract(colname, 
          into = c("sample", "strand"), 
          regex = ".*\\.([[:alnum:]]+)\\.accepted_hits_([[:alnum:]]{1}).bam",
          remove = F) %>%
  dplyr::mutate(sample_strand = paste0(sample, "_", strand))

names(counts) = c("id", header_samples$sample_strand)

counts_strand = data.frame(
  id = counts$id,
  counts_f = counts %>%
    dplyr::select(matches("_f")) %>%
    rowSums,
  counts_r = counts %>%
    dplyr::select(matches("_r")) %>%
    rowSums,
  stringsAsFactors = F
)

counts_strand = counts_strand %>%
  dplyr::mutate(total = counts_f + counts_r,
                difference = abs(counts_f - counts_r),
                perc = (difference/total) * 100, 
                direction = ifelse(perc >= snakemake@params[['cutoff']], "uni", "bi"))

write.table(counts_strand,
            snakemake@output[[1]], 
            quote = F,
            col.names = T,
            sep = "\t", 
            row.names = F)