suppressPackageStartupMessages(library(dplyr))

counts = read.table(snakemake@input[[1]], 
                    sep = "\t", 
                    header = T, 
                    stringsAsFactors = F)

counts = counts %>%
  dplyr::select(-c(Chr, Start, End, Strand, Length))

counts_strand = data.frame(
  id = counts$Geneid,
  counts_fwd = counts %>%
    dplyr::select(matches(".fwd.")) %>%
    rowSums,
  counts_rev = counts %>%
    dplyr::select(matches(".rev.")) %>%
    rowSums,
  stringsAsFactors = F
)

counts_strand = counts_strand %>%
  dplyr::mutate(total = counts_fwd + counts_rev,
                difference = abs(counts_fwd - counts_rev),
                perc = (difference/total) * 100, 
                direction = ifelse(perc >= snakemake@params[['cutoff']], "uni", "bi"))

write.table(counts_strand,
            snakemake@output[[1]], 
            quote = F,
            col.names = T,
            sep = "\t", 
            row.names = F)