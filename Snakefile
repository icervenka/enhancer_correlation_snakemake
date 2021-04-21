import pandas as pd
from snakemake.shell import shell
from snakemake.utils import validate, min_version
from snakemake.io import load_configfile

##### set minimum snakemake version #####
min_version("5.7.0")

##### load common constants and functions #####
include: "snakemake/rules/common.smk"
include: "snakemake/rules/functions.smk"

##### load config file #####
configfile: "config.yaml"
#validate(config, schema="snakemake/schema/config.schema.yaml")

##### sample/file metadata information #####
enhancer_metadata = pd.read_table(config["enhancer_metadata"], dtype=str)
#validate(Metadata, schema="snakemake/schema/metadata.schema.yaml")

rnaseq_metadata = pd.read_table(config["rnaseq_metadata"], dtype=str)
#validate(Metadata, schema="snakemake/schema/metadata.schema.yaml")
Samples = list(rnaseq_metadata["sample"])
Files = list(rnaseq_metadata["file"])

##### top level snakemake rule #####
rule all:
    input:
        expand(OUTPUT_BAMDIR + "{sample}/{strand}.bam", sample=Samples, strand = ["fwd", "rev"]),
        OUTPUT_COMBINED + "include.bed",
        OUTPUT_COMBINED + "enhancers.bed",
        OUTPUT_COMBINED + "enhancers.gtf",
        OUTPUT_ANNOTATED + "enhancers.filtered.bed",
        OUTPUT_ANNOTATED + "enhancers.filtered.extended.bed",
        OUTPUT_ANNOTATED + "genes.filtered.bed",
        OUTPUT_ANNOTATED + "enhancers.genes.bed",
        OUTPUT_CORRELATED + "enhancer_gene_annotation.txt",
        OUTPUT_CORRELATED + "enhancer_gene_expression.txt",
        OUTPUT_CORRELATED + "enhancer_gene_correlation.txt",
        OUTPUT_COUNTDIR + "gene_counts.txt",
        OUTPUT_COUNTDIR + "enhancer_counts.txt",
        OUTPUT_COUNTDIR + "enhancer_counts_strand.txt",
        OUTPUT_COUNTDIR + "alignment_counts.txt",
        log_files()

##### load remaining pipleline rules #####
include: "snakemake/rules/process_bed.smk"
include: "snakemake/rules/process_rnaseq.smk"
include: "snakemake/rules/identify_enhancers.smk"
include: "snakemake/rules/diffexp_correlation.smk"
