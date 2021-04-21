rule filter_expression:
    input:
        rules.enhancer_counts.output.counts,
        rules.rnaseq_counts.output.counts
    output:
        OUTPUT_ANNOTATED + "enhancers.filtered.bed",
        OUTPUT_ANNOTATED + "genes.filtered.bed"
    params:
        enhancer_read_cutoff=config['enhancer_read_cutoff'],
        gene_read_cutoff=config['gene_read_cutoff']
    script:
        "../scripts/filter_expression.R"

rule extend_enhancers:
    input:
        rules.filter_expression.output[0]
    output:
        OUTPUT_ANNOTATED + "enhancers.filtered.extended.bed"
    params:
        chrom_sizes=config['chrom_sizes'],
        extend_by=config['max_annotation_distance']
    run:
        shell(
            "bedtools slop "
            "-g {params.chrom_sizes} "
            "-b {params.extend_by} "
            "-i {input} > {output} "
        )

rule annotate_enhancers:
    input:
        enhancers=rules.extend_enhancers.output,
        genes=rules.filter_expression.output[1]
    output:
        OUTPUT_ANNOTATED + "enhancers.genes.bed"
    run:
        shell(
            "bedtools intersect "
            "-sorted "
            "-wa -wb "
            "-a {input.enhancers} "
            "-b {input.genes} > {output} "
        )
