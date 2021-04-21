suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(DESeq2))

### functions ###
# chrom_to_UCSC = function(arr) {
#   map_chr(arr, function(x) {
#     if(!grepl("^chr", x)) {
#       paste0("chr", x)
#     } else {
#       x
#     }
#   })
# }
# 
# chrom_to_ENSEMBL = function(arr){
#   map_chr(arr, function(x) {
#     if(grepl("^chr", x)) {
#       stringr::str_replace(x, "chr", "")
#     } else {
#       x
#     }
#   })
# }

gtf2bed = function(input, extract_attribute = NULL, unify_seqid = F, FUN = NULL) {
  gtf_header = c("seqid", "source", "feature", "start",  "end", "score", "strand", "frame", "attribute")
  names(input) = gtf_header
  
  bedfile = input %>%
    mutate(start = start + 1)
  
  if(unify_seqid) {
    bedfile = bedfile %>%
      mutate(seqid = FUN(seqid))
  } 
  
  if(!is.null(extract_attribute)) {
    bedfile = bedfile %>%
      tidyr::extract(attribute,
                     into = c("name"), 
                     regex = paste0(extract_attribute, " ([[:graph:]]+);"),
                     remove = T)
  } else {
    bedfile = bedfile %>%
      mutate(name = "none")
  }
  
  bedfile = bedfile %>%
    select(chrom = seqid, 
                  chromStart = start, 
                  chromEnd = end, 
                  name, 
                  score, 
                  strand)
  return(bedfile)
}

featurecounts2bed = function(input) {
  input = input %>%
    mutate(chrom = Chr) %>%
    mutate(score = ".") %>%
    mutate(Start = as.numeric(Start) - 1,
           End = as.numeric(End)) %>%
    select(chrom,
           chromStart = Start,
           chromEnd = End,
           name = Geneid,
           score,
           strand = Strand)
  return(input)
}

# this is only a workaround, the sample mapping should be created somewhere in a rule to be consistent
featurecount_sample_name_map = function(fc_samples, metadata, sample_colname = "sample") {
  sample_ind = map_int(fc_samples, function(x) {
    # the dots are created when converted from directory slashes, 
    # need to be present as pattern and string delimiters
    ind = which(str_detect(paste0(".", x, "."), 
                           paste0("\\.", rnaseq_metadata[["sample"]], "\\.")))
    
    if(!purrr::is_empty(ind)) {
      return(ind)
    } else {
      return(NA)
    }
  })
  
  data.frame(column = seq(length(sample_ind)),
             sample_n = sample_ind) %>%
    drop_na() %>%
    left_join(metadata %>% 
                rownames_to_column(var = "sample_n") %>%
                mutate(sample_n = as.numeric(sample_n)), 
              by = "sample_n")
}

fc_counts_to_matrix = function(count_file, fc_sample_map) {
  m = count_file %>% 
    select(fc_sample_map$column) %>% 
    `rownames<-`(count_file$Geneid) %>%
    `names<-`(fc_sample_map$sample)
}

normalize_counts = function(count_file, fc_sample_map, normalization = "none") {
  df = fc_counts_to_matrix(count_file, fc_sample_map)
  
  if(normalization == "deseq2") {
    dds = DESeqDataSetFromMatrix(countData = df,
                                 colData = fc_sample_map,
                                 design = ~ group)
    dds = DESeq(dds)
    count_norm = counts(dds, normalized = TRUE)
  } else if(normalization == "rpm") {
    count_norm = df / fc_sample_map$alignment_count
  } else {
    count_norm = df
  }
  return(count_norm)
}

name_contrast = function(c1, c2, sep = "_vs_") {
  return(paste0(c1, sep, c2))
}


create_ref_level_df = function(x) {
  combn(x, 2) %>% 
    t %>% 
    as.data.frame(stringsAsFactors = F) %>%
    dplyr::mutate(cond = contrast) %>%
    dplyr::select(cond, V2, V1) %>%
    dplyr::mutate(name = name_contrast(V2, V1))
}

deseq_diffexp = function(count_file, rnaseq_metadata, ref_levels, fdr, contrast = "group") {
  # create a mapping df between featurecounts columns and sample IDs
  fc_sample_map = featurecount_sample_name_map(names(count_file) %>% unique, rnaseq_metadata)
  
  m = fc_counts_to_matrix(count_file, fc_sample_map)
  dds = DESeqDataSetFromMatrix(countData = m,
                               colData = fc_sample_map,
                               design = as.formula(paste0("~ ", contrast)))

  if(length(ref_levels) > 1) {
    dds[[contrast]] = factor(dds[[contrast]], levels = ref_levels)
  } else {
    dds[[contrast]] = relevel(dds[[contrast]], ref_levels)
  }
  
  dds = DESeq(dds)
  
  selected_contrasts = create_ref_level_df(levels(dds[[contrast]]))
  
  map_contrasts = selected_contrasts %>% 
    dplyr::select(-name) %>% 
    t %>% 
    as.data.frame(stringsAsFactors = F)
  
  # adapted from rnaseq_pipeline
  result_array = map(map_contrasts, function(x) { 
    sc = x 
    res = results(dds, contrast=sc, alpha = fdr)
    if(lfc_shrink == "apeglm") {
      suppressMessages(library(apeglm))
      #TODO
      lfc_res = lfcShrink(dds, contrast=sc, res=res, type = "ashr")
    } else if(lfc_shrink == "ashr") {
      suppressMessages(library(ashr))
      lfc_res = lfcShrink(dds, contrast=sc, res=res, type = "ashr")
    } else if(lfc_shrink == "normal" & contrast_type == "C"){
      message("When selecting numeric contrast and normal lfc shrinkage need following number of contrasts specified:")
      print(resultsNames(dds))
      lfc_res = res
    } else {
      lfc_res = lfcShrink(dds, contrast=sc, res=res, type = "normal")
    }
    lfc_res[order(lfc_res$padj),]
  })
  
  names(result_array) = selected_contrasts$name
  return(result_array)
}

write.table.bed = function(x, file) {
  write.table(x,
              file, 
              quote = F,
              col.names = F,
              sep = "\t", 
              row.names = F)
}

write.table.gtf = function(x, file) {
  write.table(x,
              file, 
              quote = F,
              col.names = T,
              sep = "\t", 
              row.names = F)
}

### default gtf file header ###
gtf_header = c("seqid", "source", "feature", "start",  "end", "score", "strand", "frame", "attribute")

### default bed file header ###
bed_header = c("chrom", "chromStart", "chromEnd", "name",  "score", "strand")

### default featurecounts file columns other than sample names ###
featurecounts_params = c("Geneid", "Chr", "Start", "End", "Strand", "Length")
