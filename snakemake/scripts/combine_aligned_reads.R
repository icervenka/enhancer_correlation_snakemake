suppressPackageStartupMessages(library(dplyr))

read_count_files = snakemake@input
read_counts = purrr::map_dfr(read_count_files, function(x) {
  data.frame(file = x, alignment_count = readr::read_lines(x, n_max=1))
}) %>%
  tidyr::extract(file, into = "sample", regex = "/([[:alnum:]]+)_aligned_reads") %>%
  select(sample, alignment_count)

write.table(read_counts,
            snakemake@output[[1]], 
            quote = F,
            col.names = T,
            sep = "\t", 
            row.names = F)