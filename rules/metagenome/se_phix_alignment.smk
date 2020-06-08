rule phix_alignment:
    input:
        R1 = 'sample_output/aligned/unmapped/{sample}_pe_unmapped_R1.fastq.gz'
    output:
        bam = temp('sample_output/phix/{sample}.bam'),
        unmapped_R1 = temp('sample_output/phix/unmapped/{sample}_unmapped_R1.fastq')
    params:
        genome_prefix = 'references/phix/phix',
        out_files = 'sample_output/phix/unmapped/{sample}_unmapped_R%.fastq'
    shell:
        """
        bowtie2 -x {params.genome_prefix} \
                -U {input.R1} \
                --local --very-sensitive-local \
                --un {output.unmapped_R1} | samtools view -bh - > {output.bam}
        """
