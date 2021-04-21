rule combine_expression:
    input:
        enhancer_annotation=rules.annotate_enhancers.output,
        enhancer_directionality=rules.prettify_strand_count_file.output,
        enhancer_counts=rules.enhancer_counts.output.counts,
        gene_counts=rules.rnaseq_counts.output.counts,
        alignment_counts=rules.combine_aligned_reads.output
    output:
        annotation=OUTPUT_CORRELATED + "enhancer_gene_annotation.txt",
        expression=OUTPUT_CORRELATED + "enhancer_gene_expression.txt"
    params:
        rnaseq_metadata=config["rnaseq_metadata"],
        normalize_counts=config['normalize_counts']
    script:
        "../scripts/combine_expression.R"

rule deseq_diffexp_correlation:
    input:
        enhancer_counts=rules.enhancer_counts.output.counts,
        gene_counts=rules.rnaseq_counts.output.counts,
        annotation=rules.combine_expression.output.annotation,
        expression=rules.combine_expression.output.expression
    output:
        OUTPUT_CORRELATED + "enhancer_gene_correlation.txt"
    params:
        rnaseq_metadata=config["rnaseq_metadata"],
        outdir=OUTPUT_CORRELATED,
        contrast=config["contrast"],
        reference_levels=config["reference_levels"],
        lfc_shrink=config["lfc_shrink"],
        fdr=config["fdr"],
        enhancer_gene_correlation_cutoff=config["enhancer_gene_correlation_cutoff"]
    script:
        "../scripts/deseq_diffexp_correlation.R"
