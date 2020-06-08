rule Metilene_function:
    input:
        MethylMatrix_chr = 'references/reference_methylomes/MethylMatrix/{chr}/MethylMatrix_{chr}'
    output:
        compared_groups = 'references/reference_methylomes/MethylMatrix/{chr}/tissue_comparisons/{comp_group}_{chr}'
    log: 'logs/Metilene/{comp_group}_{chr}'
    params:
        groupA= lambda wildcards: string_split(wildcards.comp_group, '_',1)[0],
        groupB= lambda wildcards: string_split(wildcards.comp_group, '_',1)[1],
        M = '100', #originally 100
        m = '10', #originally 10
        d = '0.2' #originally 0.2
    shell:
        """
        ({METILENE} -M {params.M} -m {params.m} -d {params.d} -a {params.groupA}_ -b {params.groupB}_ {input.MethylMatrix_chr} |
            sort -V -k1,1 -k2,2n > {output.compared_groups}) 2> {log}
        """

################################################################################
# GOLDEN MARKERS METHYLMATRIX
################################################################################
rule golden_markers:
    input:
        expand('references/reference_methylomes/MethylMatrix/{chr}/tissue_comparisons/{comp_group}_{chr}', comp_group=comparing_groups, chr=AUTOSOMALCHROMO.split(' '))
    output:
        tissue_chr_DMR = 'references/reference_methylomes/{tissue}/{tissue}_{chr}_DMR.bed'
    threads: 1
    params:
        min_groups='15'
    shell:
        """
        files_to_compare=$(find references/reference_methylomes/MethylMatrix/{wildcards.chr}/tissue_comparisons/ -name "*{wildcards.tissue}*")
        files_to_compare=$(echo $files_to_compare)
        bedtools multiinter -i $files_to_compare | awk '{{if ($4>={params.min_groups}) print $0}}' > {output.tissue_chr_DMR}
        """

rule aggregate_regions:
    input:
        expand('references/reference_methylomes/{tissue}/{tissue}_{chr}_DMR.bed', tissue = tissues, chr = AUTOSOMALCHROMO.split(' '))
    output:
        'references/reference_methylomes/processed_refs/golden_markers.tmp.bed'
    shell:
        """
        cat {input} | bedtools sort -i - | bedtools merge -i - > {output[0]}
        """

rule binning_MM_on_golden_markers:
    input:
        MethylMatrix_body = 'references/reference_methylomes/MethylMatrix/MethylMatrix.body',
        bedfile= 'references/reference_methylomes/processed_refs/golden_markers.tmp.bed'
    output:
        MethylMatrix_sbp = 'references/reference_methylomes/processed_refs/MethylMatrix.sbp.bed',
        MethylMatrix_tmp = 'references/reference_methylomes/processed_refs/MM.tmp',
        MethylMatrix_binned = 'references/reference_methylomes/processed_refs/MethylMatrix_binned',
        gm = 'references/reference_methylomes/processed_refs/golden_markers.bed'
    shell:
        """
        grep -v "-" {input.MethylMatrix_body} > {output.MethylMatrix_sbp}

        bedtools intersect -a {input.bedfile} -b {input.MethylMatrix_body} -sorted -wa -wb > {output.MethylMatrix_tmp}

        Rscript scripts/reference_methylomes/aggregate_methylmatrix.R {output.MethylMatrix_tmp} {output.MethylMatrix_binned}
        cut -f1,2,3 {output.MethylMatrix_binned} > {output.gm}
        """
