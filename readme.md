# Snakemake pipeline for enhancer identification and quantification

This pipeline will overlay enhancer information from genomic regions supplied as
BED files with RNASeq data to identify correlations between gene and enhancer
expression.

---

## Input folders

### enhancer

BED files with genomic regions to include as potential enhancers. Pipeline will
combine data from several types of experiments such as ATACSeq, GROSeq, CHIPSeq
of various histone marks as well as removing ones that fall within coding regions.

Each of the supplied BED files needs to have an entry in the **enhancer_metadata**
file with following information specified.
- **merge_distance** - in the first step, all the regions within this distance are merged
- **min_frag_size, max_frag_size** - minimal and maximal length of the region in bp.
Regions below and above the threshold will be discarded
- **extend_left, extend_right** - Number of bp the genomic region should be extended
to the left and right, eg. broaden the ATACSeq peaks or to account for polyA readthrough.
If strand_specific is set to 'yes', extend_left will extend 5' and extend_right will
extend 3' direction of the region.
- **strand_specific** - If the BED files contains strand information, you can extend your
region in a strand specific manner
- **group** - grouping of the bed files either based on method or other variable. Currently
has no effect on the outcome, implemented for future extensions.
- **operation**
  - include - BED files which contain regions of putative enhancers
  - exclude - BED files which contain regions to be excluded, such as gene coding
  sequences

### rnaseq

BAM files containing both gene and enhancer expression data in format \<sample\>/\<file\>
that are reflected in the **rnaseq_metadata** file. This structure can be obtained
by using our snakemake bulk rnaseq pipeline in the 'align only' mode.

Following information needs to be provided:
- **sample** - name of the sample
- **file** - name of the BAM file
- **paired** - Were the the corresponding fastq files paired? To quantify directionality of enhancer transcription reads are split into forward and reverse. Current implementation sets this correctly for Illumina sequencing.
- **group** - group the sample belongs to. For comparing differential expression of genes and enhancers.

---
## Output folders

### enhancers
Contains all intermediate processing steps of supplied bed files as well as
identified and annotated enhancers. The folder are arranged in order of processing
with following structure:

- **01_converted** - BED files are converted to chromosome names corresponding
to BAM files
- **02_merged** - Regions of BED files within merge distance are merged together
- **03_filtered** - Regions in BED files are filtered based on min and max allowed
fragment size
- **04_extended** - BED regions are extended to the left and right based on the
metadata file
- **05_grouped** - BED files are grouped according to metadata file and overlapping
regions are merged
- **06_combined** - Grouped BED files to be included are combined into one and overlapping
regions are merged (creates include.bed). Regions in files that are marked as 'exclude'
are overlapped and removed (creates enhancer.bed and enhancer.gtf)
- **07_annotated** - Enhancers from previous step are extended based on **max_annotation_distance** parameter and combined with genes from supplied **gtf** file, selecting the ones that overlap
with particular enhancer
- **08_correlated** - Differentially expressed genes and enhancers are calculated by DESeq2 package based on provided contrasts and correlation is calculated both globally and on contrast level for for each gene-enhancer pair from annotation

### counts
Contains counts for genes, enhancers and total aligned reads.

### strand_split_bam
Contains strand specific BAM files to count the directionality of the enhancer transcription.
To save space, only the reads that overlap with regions from gtf file containing
putative enhancers are stored.

---
## Config

Following parameters need to be specified in the config.yaml file:

### _General_

#### experiment_name
Experiment identifier

#### threads
Maximum number of threads to use for calculations

#### species
Species that the samples come from

#### enhancer_metadata
Location of metadata file containing information about genomic regions of putative
enhancers

#### rnaseq_metadata
Location of metadata file containing information about BAM files for RNASeq

### _RNASeq related_

#### bam_chrom_style
Chromosome nomenclature from genome assembly. See **chrom_conversion_table**

#### gtf
Location of gtf file for counting of alignments

#### is_rnaseq_stranded
Strandedness of supplied BAM files, one of: 'yes', 'no', 'reverse'. Current
implementation only supports same strandedness for all files.

### _BED file related_

#### chrom_sizes
Location of BED file with chromosome sizes used for extending features
can be created from genome fasta file by:
```bash
samtools faidx genome.fa
awk '{print $1"\t"$2}' genome.fa.fai > chrom_sizes.genome
```
#### chrom_conversion_table
BED files and RNASeq files might be created using assemblies with different
chromosome nomenclatures. In order to convert between chromosome names either specify species and choromosome name style of BAM files or a location of conversion table currently supported conversions for assemblies:
- **ucsc** - begins with 'chr' string
- **ensembl** - only number wihtout 'chr' string

for species:

- **mouse**
- **human**

for other combinations you need to generate your own conversion table. Please
refer to 'snakemake/data/*.chrom.names' files for expected format.

#### preserve_strand
Whether the strandedness of the features in BED files for enhancer identification
should be preserved. Requires that the strand information in present in all
supplied BED files

### _Enhancer identification related_

#### combine
How to combine enhancer fragments from individual BED files that are to be included. One of:
- **merge** -  uses bedtools merge
- **intersect** - uses bedtools intersect (not yet implemented)
- **no** - all fragments are preserved

#### directionality_cutoff
Percentage threshold of reads that align to one strand that will determine transcription from enhancer is considered uni- or bidirectional.

#### max_annotation_distance
Max distance of gene from enhancer center to be included for annotation for
particular enhancer

#### gene_read_cutoff
Minimum total raw reads per gene

#### enhancer_read_cutoff
Minimum total raw reads per enhancer

#### normalize_counts
Normalize enhancer and gene counts per million reads. One of:
- **none** - do not normalize read numbers
- **rpm** - reads per million (not yet implemented)
- **deseq2** - DESeq2 normalization (not yet implemented)

### _Differential expression and correlation related_

#### contrast
Column in the **rnaseq_metadata**, which will be the basis of the group comparison.

#### reference_levels
List of item in contrast column in the order you want the comparisons constructed,
eg. ["ctrl","ko"]. The baseline condition is always to the left.

#### lfc_shrink
Log2Fold change shrinkage type applied by DESeq2. Can be either 'normal' for legacy
shrinking or 'ashr'. See [DESeq2 vignette](https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#alternative-shrinkage-estimators) for more information.

#### fdr
False discovery rate cutoff to determine differentially expressed genes and enhancers.

#### enhancer_gene_correlation_cutoff
Minimum Pearson correlation coefficient of gene-enhancer expression. Included are the
genes that pass the threshold either in global or in contrast specific correlation.
