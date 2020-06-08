rule mapping_stats:
    input:
        original_r1 = CLUSTER+config['DATA']+'samples/{sample}_R1.fastq.gz',
        trimmed_r1 = 'sample_output/trim/{sample}_R1_trim.fastq',
        raw_bam = 'sample_output/aligned/raw_aligned/{sample}.bam',
        deduped_bam = 'sample_output/aligned/all_chr/{sample}_mapped_all_chr.bam',
        autosomal_bam = 'sample_output/aligned/autosomal/{sample}_mapped_autosomal.bam',
        unmapped_r1 = 'sample_output/aligned/unmapped/{sample}_pe_unmapped_R1.fastq.gz',
        decon_r1 = 'sample_output/nonhuman_fasta/R1/{sample}_R1.fa'
    output:
        stats = 'sample_output/stats/{sample}_mapping_stats.txt'
    params:
        mappable_hg19='2684573069'
    shell:
        """
        original_reads=$(pigz -dc {input.original_r1} | wc -l)
        original_reads=$((original_reads / 4))
        trimmed_reads=$(wc -l {input.trimmed_r1} | cut -d' ' -f1)
        trimmed_reads=$((trimmed_reads / 4))
        mapped_reads=$(samtools view -c {input.raw_bam})
        mapped_reads=$((mapped_reads / 2))
        deduped_reads=$(samtools view -c {input.deduped_bam})
        deduped_reads=$((deduped_reads / 2))
        mapping_eff=$(echo "scale=2;$mapped_reads/$trimmed_reads" | bc)
        deduped_frac=$(echo "scale=2;$deduped_reads/$mapped_reads" | bc)
        depth_var=$(samtools depth {input.autosomal_bam} | awk '{{sum+=$3}} END {{print sum/{params.mappable_hg19}" "NR/{params.mappable_hg19}}}')
        depth=$(echo $depth_var | cut -f1 -d' ')
        bp_frac=$(echo $depth_var | cut -f2 -d' ')

        unmapped=$(pigz -dc {input.unmapped_r1} | wc -l)
        unmapped=$((unmapped / 4))
        decon=$(wc -l {input.decon_r1} | cut -d' ' -f1)
        decon=$((decon / 4))

        echo -e "NUM_READS\tREADS_AFTER_TRIM\tALIGNED\tMAPPING_EFFICIENCY\tDEDUPED_READS\tDEDUP_FRAC\tDEPTH\tFRACTION_PB_BP\tUNMAPPED_READS\tDECONTAMINATED_READS" > {output.stats}
        echo -e "$original_reads\t$trimmed_reads\t$mapped_reads\t$mapping_eff\t$deduped_reads\t$deduped_frac\t$depth\t$bp_frac\t$unmapped\t$decon" >> {output.stats}
        """
