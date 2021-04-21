suppressPackageStartupMessages(library(dplyr))

enh = read.table(snakemake@input[[1]],
                 sep = "\t",
                 header = F,
                 stringsAsFactors = F) %>%
  setNames(c("seqname", "start", "end"))

enh = enh %>%
  dplyr::mutate(start = start + 1) %>%
  dplyr::mutate(source = "ENSEMBL") %>%
  dplyr::mutate(feature = "gene") %>%
  dplyr::mutate(score = ".") %>%
  dplyr::mutate(strand = ".") %>%
  dplyr::mutate(frame = ".") %>%
  dplyr::mutate(attribute = paste0('gene_id "', seqname, ":", floor((end - start)/2 + start),'";')) %>%
  dplyr::select(seqname, source, feature, start, end, score, strand, frame, attribute)

gtf_header = paste("##description: putative enhancers identified by combination of bed files", sep="\n")
gtf_header = paste(gtf_header, "##provider: Igor Cervenka", sep="\n")
gtf_header = paste(gtf_header, "##format: gtf", sep="\n")

write.table(gtf_header,
            snakemake@output[[1]], 
            quote = F,
            col.names = F,
            sep = "", 
            row.names = F)

write.table(enh, 
            snakemake@output[[1]], 
            quote = F,
            col.names = F,
            sep = "\t", 
            row.names = F,
            append = T)
