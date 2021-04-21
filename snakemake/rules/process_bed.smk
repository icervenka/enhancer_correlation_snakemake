ps = bedtools_preserve_strand(config['preserve_strand'])

rule covert_chrom_names:
    input:
        INPUT_ENHANCER_DIR + "{file}"
    output:
        OUTPUT_CONVERTED + "{file}"
    params:
        species=config['species'],
        bam_chrom_style=config['bam_chrom_style'],
        chrom_conversion_table=config['chrom_conversion_table']
    threads:
        1
    script:
        "../scripts/convert_chrom_names.R"

rule merge_individual_bed:
    input:
        rules.covert_chrom_names.output
    output:
        tmp=temp(OUTPUT_MERGED + "{file}.tmp"),
        bed=OUTPUT_MERGED + "{file}"
    params:
        dist=get_merge_distance,
        preserve_strand=ps
    threads:
        1
    shell:
        """
        sort -k 1,1 -k 2,2n {input} > {output.tmp};
        bedtools merge {params.preserve_strand} -d {params.dist} -i {output.tmp} > {output.bed}
        """

rule filter_bed:
    input:
        rules.merge_individual_bed.output.bed
    output:
        OUTPUT_FILTERED + "{file}"
    params:
        min=get_min_frag_size,
        max=get_max_frag_size
    threads:
        1
    shell:
        """
        awk -F"\t" "BEGIN {{ OFS = FS }} {{ if((\$3-\$2 > {params.min}) && (\$3-\$2 < {params.max})) print }}" {input} > {output};
        """

rule extend_bed:
    input:
        rules.filter_bed.output
    output:
        OUTPUT_EXTENDED + "{file}"
    params:
        chrom_sizes=config['chrom_sizes'],
        extend_left=get_extend_left,
        extend_right=get_extend_right,
        preserve_strand=ps
    threads:
        1
    run:
        shell(
            "bedtools slop "
            "-g {params.chrom_sizes} "
            "{params.preserve_strand} "
            "-l {params.extend_left} "
            "-r {params.extend_right} "
            "-i {input} > {output} "
        )

rule group_bed:
    input:
        get_extended_beds
    output:
        tmp=temp(OUTPUT_GROUPED + "{group}.tmp"),
        bed=OUTPUT_GROUPED + "{group}.bed"
    params:
        preserve_strand=ps
    threads:
        1
    run:
        input = " ".join(input)
        shell(
            "cat {input} | sort -k 1,1 -k 2,2n > {output.tmp}; "
            "bedtools merge {params.preserve_strand} -i {output.tmp} > {output.bed}"
        )

# TODO finish implementing the different combine
rule combine_include:
    input:
        get_include
    output:
        tmp=temp(OUTPUT_COMBINED + "include.tmp"),
        bed=OUTPUT_COMBINED + "include.bed"
    params:
        combine=config['combine']
    threads:
        1
    run:
        if params.combine == "intersect":
            shell(
                "cat {input} | sort -k 1,1 -k 2,2n > {output.tmp}; "
                "bedtools merge -i {output.tmp} > {output.bed}"
            )
        elif params.combine == "merge":
            shell(
                "cat {input} | sort -k 1,1 -k 2,2n > {output.tmp}; "
                "bedtools merge -i {output.tmp} > {output.bed}"
            )
        elif params.combine == "no":
            shell(
                "cat {input} | sort -k 1,1 -k 2,2n > {output.tmp}"
            )
        else:
            import sys
            sys.exit("Incorrect combine value for processing enhancer bed file was selected.")

rule combine_all:
    input:
        include=rules.combine_include.output.bed,
        exclude=get_exclude
    output:
        OUTPUT_COMBINED + "enhancers.bed"
    run:
        exclude = " ".join(input.exclude)
        shell(
            "bedtools intersect -wa -v -a {input.include} -b {exclude} > {output}"

        )

rule enhancers2gtf:
    input:
        rules.combine_all.output
    output:
        OUTPUT_COMBINED + "enhancers.gtf"
    script:
        "../scripts/enhancers2gtf.R"
