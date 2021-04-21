import sys

### metadata selector functions ###
def match_row(metadata, column, equals_to):
    v = metadata.loc[metadata[column] == equals_to, ].dropna()
    return v

def get_merge_distance(wildcards):
    v = match_row(enhancer_metadata, "file", wildcards.file)
    v = v['merge_distance'].iloc[0]
    return v

def get_min_frag_size(wildcards):
    v = match_row(enhancer_metadata, "file", wildcards.file)
    v = v['min_frag_size'].iloc[0]
    return v

def get_max_frag_size(wildcards):
    v = match_row(enhancer_metadata, "file", wildcards.file)
    v = v['max_frag_size'].iloc[0]
    return v

def get_extend_left(wildcards):
    v = match_row(enhancer_metadata, "file", wildcards.file)
    return v["extend_left"].iloc[0]

def get_extend_right(wildcards):
    v = match_row(enhancer_metadata, "file", wildcards.file)
    return v["extend_right"].iloc[0]

def get_extended_beds(wildcards):
    bed = match_row(enhancer_metadata, "group", wildcards.group)
    bed = bed["file"].to_list()
    bed = [OUTPUT_EXTENDED + x for x  in bed]
    return bed

def get_exclude(wildcards):
    v = match_row(enhancer_metadata, "operation", "exclude")
    v = v['group'].to_list()
    v = [OUTPUT_GROUPED + x + ".bed" for x in v ]
    return v

def get_include(wildcards):
    v = match_row(enhancer_metadata, "operation", "include")
    v = v['group'].to_list()
    v = [OUTPUT_GROUPED + x + ".bed" for x  in v]
    return v

def get_bam_samples(wildcards):
    v = match_row(rnaseq_metadata, "sample", wildcards.sample)
    v = v['file'].to_list()
    v = [INPUT_RNASEQ_DIR + wildcards.sample + "/" + x  for x  in v]
    return v

### other ###
def bedtools_preserve_strand(s):
    if s == "yes":
        return "-s"
    else:
        return ""

def featurecounts_stranded(s):
    def stranded_switch(x):
        select = {
            "no": "-s 0 ",
            "yes": "-s 1 ",
            "reverse": "-s 2 "
        }
        return(select.get(x, "-s 0"))

    # v = match_row(rnaseq_metadata, "sample", wildcards.sample)
    # v = v['stranded'].to_list()[0]
    return stranded_switch(s)

def bam_paired(wildcards):
    v = match_row(rnaseq_metadata, "sample", wildcards.sample)
    v = v['file'].to_list()[0]
    if(v == "yes"):
        return 1
    else:
        return 0


### functions generating output files for top leve snakemake rule ###

def log_files():
    files = []
    files.append(LOG_DIR + "gene_counts.txt.summary")
    files.append(LOG_DIR + "enhancer_counts.txt.summary")
    files.append(LOG_DIR + "split_enhancer_counts.txt.summary")
    return(files)
