rule total_abundance: #allow for qubit scaling too
    input:
        mapping_stats = ancient('sample_output/stats/{sample}_mapping_stats.txt'),
        unmapped_r1 = ancient('sample_output/aligned/unmapped/{sample}_pe_unmapped_R1.fastq.gz'),
        unmapped_r2 = ancient('sample_output/aligned/unmapped/{sample}_pe_unmapped_R2.fastq.gz')
    output:
        bam = 'sample_output/total_abundance/bam/{sample}_spikeins.bam',
        total_abundance = 'sample_output/total_abundance/{sample}.txt'
    threads: 2
    params:
        spikeinref = 'references/total_abundance_controls/ucsf_spikeins/',
        outdir = 'sample_output/total_abundance/bam/',
        prep_and_seq_type = get_seq_type #first is prep, second is seq
    log: 'logs/total_abundance/{sample}.log'
    shell:
        """
        seq_type={params.prep_and_seq_type[1]}

        if [[ $seq_type == 2x* ]]; then
            {BISMARK} --local \
                    --quiet \
                    --parallel {threads} \
                    -o {params.outdir} \
                    --genome {params.spikeinref} \
                    -1 {input.unmapped_r1} -2 {input.unmapped_r2}

            cp {params.outdir}{wildcards.sample}_pe_unmapped_R1_bismark_bt2_pe.bam {output.bam}
		    mv {params.outdir}{wildcards.sample}_pe_unmapped_R1_bismark_bt2_PE_report.txt {log}

        elif [[ $seq_type == 1x* ]]; then
            {BISMARK} --local \
                    --quiet \
                    --parallel {threads} \
                    -o {params.outdir} \
                    --genome {params.spikeinref} \
                    {input.unmapped_r1}

            cp {params.outdir}{wildcards.sample}_pe_unmapped_R1_bismark_bt2.bam {output.bam}
		    mv {params.outdir}{wildcards.sample}_pe_unmapped_R1_bismark_bt2_SE_report.txt {log}
        fi

        Rscript scripts/total_abundance/spikein.R {output.bam} {input.mapping_stats} {output.total_abundance} {wildcards.sample}
        """
