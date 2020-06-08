rule references_for_decontamination: #rule is necessary because bbduk only looks at either Forward strand or F and R (not just Reverse)
    input:
        hg19_CT = 'references/hg19_Bismark/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa',
        hg19_GA = 'references/hg19_Bismark/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa',
        mic_ctl_CT = 'references/controls_Bismark/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa',
        mic_ctl_GA = 'references/controls_Bismark/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa',
        phix = 'references/phix/phix.fa'
    output:
        CT = 'references/for_decontamination/CT.fa',
        GA = 'references/for_decontamination/GA.fa',
        both = 'references/for_decontamination/both.fa'
    shell:
        """
        cp {input.phix} {output.both}
        cat {input.hg19_CT} {input.mic_ctl_CT} > {output.CT}
        cat {input.hg19_GA} {input.mic_ctl_GA} | {SEQTK} seq -r - > {output.GA}
        """
rule decontaminate:
    input:
        unmapped_R1 = 'sample_output/phix/unmapped/{sample}_unmapped_R1.fastq',
        CT = 'references/for_decontamination/CT.fa',
        GA = 'references/for_decontamination/GA.fa',
        both = 'references/for_decontamination/both.fa',
        mic_ctl_CT = 'references/controls_Bismark/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa',
        mic_ctl_GA = 'references/controls_Bismark/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa'
    output:
        r1_CT = temp('sample_output/decontaminate/{sample}.CT_R1.fastq'),
        r1_pass1 = temp('sample_output/decontaminate/{sample}.pass1_R1.fastq'),
        r1_pass2 = temp('sample_output/decontaminate/{sample}.pass2_R1.fastq'),
        decon_r1 = temp('sample_output/decontaminate/{sample}_decontaminated_R1.fastq'),
        nonhumanfa_r1 = temp('sample_output/nonhuman_fasta/R1/{sample}_R1.fa')
    resources:
        mem_mb=80000
    shell:
        """
        python scripts/metagenome/insilico_conversion.py {input.unmapped_R1} {output.r1_CT} R1

        {BBDUK} in={output.r1_CT} \
                out={output.r1_pass1} \
                -Xmx{resources.mem_mb}m \
                prealloc=t rcomp=f \
                ref={input.GA} k=50
        {BBDUK} in={output.r1_pass1} \
                out={output.r1_pass2} \
                -Xmx{resources.mem_mb}m \
                prealloc=t rcomp=f \
                ref={input.CT} k=50
        {BBDUK} in={output.r1_pass2} \
                out={output.decon_r1} \
                -Xmx{resources.mem_mb}m \
                prealloc=t \
                ref={input.both} k=20
        fastq_to_fasta -Q33 -i {output.decon_r1} -o {output.nonhumanfa_r1}
        """
