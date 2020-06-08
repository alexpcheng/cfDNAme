rule get_and_index_hg19:
    output:
        expand('references/hg19_Bismark/Bisulfite_Genome/{NUC}_conversion/BS_{NUC}.{ext}.bt2', NUC=['CT', 'GA'], ext=['1', '2', '3', '4', 'rev.1', 'rev.2']),
        expand('references/hg19_Bismark/Bisulfite_Genome/{NUC}_conversion/genome_mfa.{NUC}_conversion.fa', NUC=['CT', 'GA']),
        hg19 = 'references/hg19_Bismark/hg19.fa',
        hg19_index = 'references/hg19_Bismark/hg19.fa.fai',
        hg19_chr_lengths = 'references/hg19.chr.lengths'
    threads: indexing_threads
    params:
        genome_path = 'references/hg19_Bismark'
    shell:
        """
        wget -O {output.hg19_chr_lengths} http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
        wget -O - http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz | gunzip > {output.hg19}
        {BISMARKINDEX} --bowtie2 --parallel {threads} {params.genome_path}
        samtools faidx {output.hg19}
        """
