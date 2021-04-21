rule count_aligned_reads:
    input:
        get_bam_samples
    output:
        temp(OUTPUT_COUNTDIR + "{sample}_aligned_reads.txt")
    threads:
        config['threads']
    shell:
        "samtools view -@ {threads} -c -F 260 {input} >> {output}"

rule combine_aligned_reads:
    input:
        expand(rules.count_aligned_reads.output, sample=Samples)
    output:
        OUTPUT_COUNTDIR + "alignment_counts.txt"
    script:
        "../scripts/combine_aligned_reads.R"

rule separate_read_strands:
    input:
        bam=get_bam_samples,
        enh_overlap=rules.combine_all.output
    output:
        # TODO maybe change to temp files
        f=OUTPUT_BAMDIR + "{sample}/fwd.bam",
        r=OUTPUT_BAMDIR + "{sample}/rev.bam"
    params:
        outdir=OUTPUT_BAMDIR + "{sample}/",
        is_paired=bam_paired,
    threads:
        config['threads']
    run:
        if(params.is_paired):
            shell(
                "snakemake/scripts/split_strands_paired.sh "
                "{input.bam} "
                "{params.outdir} "
                "{output.f} "
                "{output.r} "
                "{input.enh_overlap} "
                "{threads}"
            )
        else:
            shell(
                "samtools view "
                "-bh "
                "-F 20 "
                "-L {input.enh_overlap} "
                "-@ {threads} "
                "-o {output.f} {input.bam};"
                "samtools view "
                "-bh "
                "-f 16 "
                "-L {input.enh_overlap} "
                "-@ {threads} "
                "-o {output.r} {input.bam};"
                "samtools index -@ {threads} {output.f};"
                "samtools index -@ {threads} {output.r};"
            )

rule rnaseq_counts:
    input:
        bam=expand(INPUT_RNASEQ_DIR + "{sample}/{file}", zip, sample=Samples, file=Files),
        gtf=config['gtf']
    output:
        counts=OUTPUT_COUNTDIR + "gene_counts.txt",
        summary=OUTPUT_COUNTDIR + "gene_counts.txt.summary"
    log:
        LOG_DIR + "gene_counts.log"
    params:
    # TODO it is possible to implement a chromosome lookup table via -A parameter
    # TODO it is possible to implement a comma separated list of values for individual files
    # TODO overlapping regions mess up parsing of featurecounts file to bed
        stranded=featurecounts_stranded(config['is_rnaseq_stranded'])
    threads:
        config['threads']
    run:
        shell(
            "featureCounts "
            "--largestOverlap "
            "{params.stranded} "
            "-T {threads} "
            "-t gene "
            "-a {input.gtf} "
            "-o {output.counts} {input.bam} "
            "2> {log} "
        )

rule enhancer_counts:
    input:
        bam=expand(INPUT_RNASEQ_DIR + "{sample}/{file}", zip, sample=Samples, file=Files),
        gtf=rules.enhancers2gtf.output
    output:
        counts=OUTPUT_COUNTDIR + "enhancer_counts.txt",
        summary=temp(OUTPUT_COUNTDIR + "enhancer_counts.txt.summary")
    log:
        LOG_DIR + "enhancer_counts.log"
    params:
        stranded=featurecounts_stranded("no")
    threads:
        config['threads']
    run:
        # TODO if the preserve_strand is not specified the counting only works in unstraded manner
        shell(
            "featureCounts "
            "--largestOverlap "
            "{params.stranded} "
            "-T {threads} "
            "-t gene "
            "-a {input.gtf} "
            "-o {output.counts} {input.bam} "
            "2> {log} "
        )

rule strand_enhancer_counts:
    input:
        bam=expand(rules.separate_read_strands.output, sample=Samples),
        gtf=rules.enhancers2gtf.output
    output:
        counts=OUTPUT_COUNTDIR + "split_enhancer_counts.txt",
        summary=OUTPUT_COUNTDIR + "split_enhancer_counts.txt.summary"
    log:
        LOG_DIR + "split_enhancer_counts.log"
    params:
        stranded=featurecounts_stranded("no")
    threads:
        config['threads']
    run:
        shell(
            "featureCounts "
            "--largestOverlap "
            "{params.stranded} "
            "-T {threads} "
            "-t gene "
            "-a {input.gtf} "
            "-o {output.counts} {input.bam} "
            "2> {log} "
        )

rule prettify_strand_count_file:
    input:
        rules.strand_enhancer_counts.output.counts
    output:
        OUTPUT_COUNTDIR + "enhancer_counts_strand.txt"
    params:
        cutoff=config['directionality_cutoff']
    script:
        "../scripts/prettify_strand_enhancer_counts.R"

rule move_count_logs:
    input:
        rules.rnaseq_counts.output.summary,
        rules.enhancer_counts.output.summary,
        rules.strand_enhancer_counts.output.summary,
    output:
        LOG_DIR + "gene_counts.txt.summary",
        LOG_DIR + "enhancer_counts.txt.summary",
        LOG_DIR + "split_enhancer_counts.txt.summary"
    params:
        outdir=directory(LOG_DIR)
    run:
        input = " ".join(input)
        shell(
            "cp {input} {params.outdir}"
        )
