rule estimate_BSconversion:
	input:
		mapped_autosomal = 'sample_output/aligned/autosomal/{sample}_mapped_autosomal.bam',
		hg19 = 'references/hg19_Bismark/hg19.fa',
		mic_ctl = 'references/controls_Bismark/mic.fa'
	output:
		mr = temp('sample_output/{sample}.mr'),
        bsrate = 'sample_output/conversion_rates/{sample}.bsrate.txt',
		sample_rate = 'sample_output/conversion_rates/{sample}.bsconversion'
	shell:
		"""
        {METHPIPETOMR} -m general -o {output.mr} -L 500 {input.mapped_autosomal}
		if [[ {wildcards.sample} == *MCB* ]]
		then
			{METHPIPEBSRATE} -c {input.mic_ctl} -o {output.bsrate} {output.mr}
		else
			{METHPIPEBSRATE} -c {input.hg19} -o {output.bsrate} {output.mr}
		fi
        X=$(head -1 {output.bsrate})
        echo -e "{wildcards.sample}\t$X" > {output.sample_rate}
        """
