rule get_and_index_phix:
    output:
        expand('references/phix/phix.{ext}.bt2', ext=['1', '2', '3', '4', 'rev.1', 'rev.2']),
        phix = 'references/phix/phix.fa'
    threads: indexing_threads
    params:
        genome_path = 'references/phix/phix'
    shell:
        """
        wget -O {output.phix} https://www.ncbi.nlm.nih.gov/search/api/sequence/NC_001422.1/?report=fasta
        bowtie2-build {output.phix} {params.genome_path}
        """
