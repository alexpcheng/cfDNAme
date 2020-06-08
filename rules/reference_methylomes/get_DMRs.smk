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
##### IMPROVING HERE DOWN
################################################################################
# GOLDEN MARKERS METHYLMATRIX
################################################################################
rule golden_markers:
    input:
        expand('references/reference_methylomes/MethylMatrix/{chr}/tissue_comparisons/{comp_group}_{chr}', comp_group=comparing_groups, chr=AUTOSOMALCHROMO.split(' '))
    output:
        golden_markers = 'references/reference_methylomes/MethylMatrix/golden_markers/golden_markers.bed.tmp'
    threads: 30
    shell:
        """
        python scripts/reference_methylomes/golden_markers.py --out {output.golden_markers} --threads {threads}
        """

rule binning_MM_on_golden_markers:
    input:
        MethylMatrix_body = 'references/reference_methylomes/MethylMatrix/MethylMatrix.body',
        bedfile= 'references/reference_methylomes/MethylMatrix/golden_markers/golden_markers.bed.tmp'
    output:
        MethylMatrix_tmp = 'references/reference_methylomes/MethylMatrix/golden_markers/MM.tmp',
        MethylMatrix_binned = 'references/reference_methylomes/MethylMatrix/golden_markers/MethylMatrix_binned',
        gm = 'references/reference_methylomes/MethylMatrix/golden_markers/golden_markers.bed'
    shell:
        """
        bedtools intersect -wo -a {input.bedfile} -b {input.MethylMatrix_body} -sorted |
            awk '$6-$5==1 {{print $0}}' | awk 'NF{{NF-=1}};1' > {output.MethylMatrix_tmp}
        Rscript Bin/aggregate_methylmatrix.R {output.MethylMatrix_tmp} {output.MethylMatrix_binned}
        cut -f1,2,3 {output.MethylMatrix_binned} > {output.gm}
        """
