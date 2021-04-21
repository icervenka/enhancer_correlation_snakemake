suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
source('snakemake/scripts/common.R')

bedfile = read.table(snakemake@input[[1]], stringsAsFactors = F, header = F)

species = snakemake@params[['species']]
bam_chrom_style = snakemake@params[['bam_chrom_style']]
chrom_conversion_table = snakemake@params[['chrom_conversion_table']]

if(length(chrom_conversion_table) > 1) {
  conversion_table = read.table(chrom_conversion_table,
                                stringsAsFactors = F,
                                header = T)
} else {
  conversion_table = read.table(paste0("snakemake/data/",
                                       species,
                                       ".chrom.names"),
                                stringsAsFactors = F,
                                header = T)
}

if(length(conversion_table) != 2) {
  stop("Conversion table has too many columns, only two are allowed.")
}
if(!(bam_chrom_style %in% names(conversion_table))) {
  stop("Chromosome style column in not present in conversion table.")
}

match_perecentage = map_dbl(1:2, function(x, chr_vec) {
  chr_vec %in% conversion_table[[x]] %>% sum() / length(chr_vec)
}, chr_vec = bedfile[[1]][sample(1:dim(bedfile)[1], 1000)])

if(max(match_perecentage) < 0.8) {
  stop("Chromosome style in input file match frequency is insufficient. Are the chromosome names consistent?")
}

source_col_ind = which.max(match_perecentage)
dest_col_ind = which(names(conversion_table) == bam_chrom_style)

if(source_col_ind != dest_col_ind) {
  new_chrom = merge(bedfile,conversion_table, 
                    by.x=names(bedfile)[1], # chrom column
                    by.y=names(conversion_table)[source_col_ind]) %>%
    mutate(V1 = !!as.symbol(names(conversion_table)[dest_col_ind])) %>%
    select(-!!as.symbol(names(conversion_table)[dest_col_ind]))
} else {
  new_chrom = bedfile
}

write.table.bed(new_chrom, snakemake@output[[1]])
