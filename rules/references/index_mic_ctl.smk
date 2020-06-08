rule get_and_index_microbial_ctl:
    input:
        controls = 'references/controls_Bismark/mic.fa'
    output:
        expand('references/controls_Bismark/Bisulfite_Genome/{NUC}_conversion/BS_{NUC}.{ext}.bt2', NUC=['CT', 'GA'], ext=['1', '2', '3', '4', 'rev.1', 'rev.2']),
        expand('references/controls_Bismark/Bisulfite_Genome/{NUC}_conversion/genome_mfa.{NUC}_conversion.fa', NUC=['CT', 'GA']),
    threads: indexing_threads
    params:
        genome_path = 'references/controls_Bismark'
    shell:
        """
        {BISMARKINDEX} --bowtie2 --parallel {threads} {params.genome_path}
        """
