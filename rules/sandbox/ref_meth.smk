rule bin_on_regions:
    input:
        expand('reference_methylomes/MethylMatrix/{chr}/tissue_comparisons/{comp_group}_{chr}', comp_group=comparing_groups, chr=AUTOSOMALCHROMO.split(' ')),
        MethylMatrix_body='reference_methylomes/MethylMatrix/MethylMatrix.body'
    output:
        all_DMR_tmp = 'sandbox/reference_methylomes/all_DMRs.tmp',
        all_DMR = 'sandbox/reference_methylomes/all_DMRs.txt',
        binned='sandbox/reference_methylomes/binned.bed'
    shell:
        """
        find ./reference_methylomes/MethylMatrix/ -name "*_*_chr*" | while read F; do cat ${{F}} >> {output.all_DMR_tmp}; done
        bedtools sort -i {output.all_DMR_tmp} | uniq > {output.all_DMR}
        bedtools intersect -wa -wb -a {output.all_DMR} -b {input.MethylMatrix_body} -sorted | cut -f1,2,3,14- > {output.binned}
        awk '{{print $0 >> "binned."$1".bed"}}' {output.binned}
        Rscript aggregate.R binned.chr1.bed ../../reference_methylomes.tsv
        """

rule merge_chr:
    input:
        expand('reference_methylomes/MethylMatrix/{chr}/tissue_comparisons/{comp_group}_{chr}', comp_group=comparing_groups, chr=AUTOSOMALCHROMO.split(' '))
    output:
        all_chr_tmp = temp('sandbox/reference_methylomes/MethylMatrix/all_chr/tissue_comparisons/{comp_group}_all_chromosomes.tmp'),
        all_chr = 'sandbox/reference_methylomes/MethylMatrix/all_chr/tissue_comparisons/{comp_group}_all_chromosomes'
    shell:
        """
        find ./reference_methylomes/MethylMatrix/ -name "*{wildcards.comp_group}*" | while read F; do cat ${{F}} >> {output.all_chr_tmp}; done
        bedtools sort -i {output.all_chr_tmp} | uniq > {output.all_chr}
        """

rule tissue_markers:
    input:
        expand('sandbox/reference_methylomes/MethylMatrix/all_chr/tissue_comparisons/{comp_group}_all_chromosomes', comp_group=comparing_groups),
        MM_body='reference_methylomes/MethylMatrix/MethylMatrix.body'
    output:
        hyper='sandbox/reference_methylomes/hyper/{tissue}',
        hypo='sandbox/reference_methylomes/hypo/{tissue}',
        DMR='sandbox/reference_methylomes/MethylMatrix/{tissue}/{tissue}_DMR.bed',
        MM_on_tissue='sandbox/reference_methylomes/MethylMatrix/{tissue}/MM_on_{tissue}_regions.txt'
    shell:
        """
        ls sandbox/reference_methylomes/MethylMatrix/all_chr/tissue_comparisons/*{wildcards.tissue}* | while read file
        do
            group1=$(echo ${{file##*/}} | cut -d'_' -f1)
            group2=$(echo ${{file##*/}} | cut -d'_' -f2)
            comp=$group1"_"$group2
            if [[ {wildcards.tissue} == $group1 ]]
            then
                cat $file | awk '$5>=0' >> {output.hyper}.$comp
                cat $file | awk '$5<=0' >> {output.hypo}.$comp
            fi
            if [[ {wildcards.tissue} == $group2 ]]
            then
                cat $file | awk '$5<=0' >> {output.hyper}.$comp
                cat $file | awk '$5>=0' >> {output.hypo}.$comp
            fi
        done
        touch {output.hyper} {output.hypo}
        bedtools multiinter -i sandbox/reference_methylomes/hy*/*{wildcards.tissue}* | awk '$4>=12' | cut -f1,2,3 > {output.DMR}
        bedtools intersect -wa -wb -a {output.DMR} -b {input.MM_body} -sorted | cut -f1,2,3,7- > {output.MM_on_tissue}
        """
rule agg:
    input:
        expand('sandbox/reference_methylomes/MethylMatrix/{tissue}/{tissue}_DMR.bed', tissue = tissues)
    output:
        'all.txt'
    shell:
        """
        touch {output}
        """
#awk '{print $0 >> $1".bed"}' example.bed
