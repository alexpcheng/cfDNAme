rule grammy:
	input:
		nonhumanfa_r1 = 'sample_output/nonhuman_fasta/R1/{sample}_R1.fa',
		nonhumanfa_r2 = 'sample_output/nonhuman_fasta/R2/{sample}_R2.fa',
		tblat1 = 'sample_output/grammy/{sample}/{sample}.tblat.1'
	output:
		nonhumanfa_gz = temp('sample_output/grammy/{sample}/{sample}.fa.gz'),
		nonhumanfasta_gz = temp('sample_output/grammy/{sample}/{sample}.fasta.gz'),
		rdt = 'sample_output/grammy/{sample}/{sample}.rdt',
		mtx = 'sample_output/grammy/{sample}/{sample}.mtx',
		lld = 'sample_output/grammy/{sample}/{sample}.lld',
		btp = 'sample_output/grammy/{sample}/{sample}.btp',
		est = 'sample_output/grammy/{sample}/{sample}.est',
		gra = 'sample_output/grammy/{sample}/{sample}.gra',
		avl = 'sample_output/grammy/{sample}/{sample}.avl'
	resources: mem_mb=1
	shell:
		"""
		prep_type=$(grep {wildcards.sample} {SEQUENCING_PREP_TABLE} | cut -f2)
		if [[ $prep_type == "SRSLY_SSLP" ]]
		then
			sed '/>/s/$/-1/' {input.nonhumanfa_r1} | gzip -1 > {output.nonhumanfa_gz}
			sed '/>/s/$/-2/' {input.nonhumanfa_r2} | gzip -1 >> {output.nonhumanfa_gz}
		else
			cat {input.nonhumanfa_r1} {input.nonhumanfa_r2} | sed 's/_1.*/-1/g' |  sed 's/_2.*/-2/g' | gzip -1 > {output.nonhumanfa_gz}
		fi
		cd sample_output/grammy/{wildcards.sample}
		python2.7 {GRAMMY_RDT} -t illumina . .
		python2.7 {GRAMMY_PRE} -q "40,75,-5" {wildcards.sample} {GRAMMY_REF_FASTA}
		python2.7 {GRAMMY_EM} -c L -b 5 -t .00001 -n 100 {wildcards.sample}.mtx
		python2.7 {GRAMMY_POST} {wildcards.sample}.est {GRAMMY_REF_FASTA} {wildcards.sample}.btp
		cd ../../../
		"""

rule annotate_grammy:
	input:
		rdt = 'sample_output/grammy/{sample}/{sample}.rdt',
		mtx = 'sample_output/grammy/{sample}/{sample}.mtx',
		lld = 'sample_output/grammy/{sample}/{sample}.lld',
		btp = 'sample_output/grammy/{sample}/{sample}.btp',
		est = 'sample_output/grammy/{sample}/{sample}.est',
		gra = 'sample_output/grammy/{sample}/{sample}.gra',
		avl = 'sample_output/grammy/{sample}/{sample}.avl',
		stats = 'sample_output/stats/{sample}_mapping_stats.txt',
		tblat1 = 'sample_output/grammy/{sample}/{sample}.tblat.1',
		LUT = 'LUTGrammy/taxids_names_lengths_tax.tab'
	output:
		tab='sample_output/grammy/{sample}/{sample}.tab',
		anno = 'sample_output/grammy/{sample}/{sample}.grammy.tab'
	params:
		DIR='sample_output/',
		DB='grammy/{sample}/',
		LUT='LUTGrammy/taxids_names_lengths_tax.tab'
	shell:
		"""
		Rscript scripts/metagenome/filter_gra_file.R {input.gra} {output.tab} {wildcards.sample}
		Rscript scripts/metagenome/annotate_grammy_apc.R {output.tab} {input.tblat1} {input.stats} {input.LUT} {output.anno}
		"""

rule cleanup:
	input:
		anno = 'sample_output/grammy/{sample}/{sample}.grammy.tab',
		decon_r1 = 'sample_output/decontaminate/{sample}_decontaminated_R1.fastq',
		decon_r2 = 'sample_output/decontaminate/{sample}_decontaminated_R2.fastq',
		nonhumanfa_r1 = 'sample_output/nonhuman_fasta/R1/{sample}_R1.fa',
		nonhumanfa_r2 = 'sample_output/nonhuman_fasta/R2/{sample}_R2.fa'
	output:
		decon_r1_gz = 'sample_output/decontaminate/{sample}_decontaminated_R1.fastq.gz',
		decon_r2_gz = 'sample_output/decontaminate/{sample}_decontaminated_R2.fastq.gz',
		nonhumanfa_r1_gz = 'sample_output/nonhuman_fasta/R1/{sample}_R1.fa.gz',
		nonhumanfa_r2_gz = 'sample_output/nonhuman_fasta/R2/{sample}_R2.fa.gz'
	shell:
		"""
		gzip -c {input.decon_r1} > {output.decon_r1_gz}
		gzip -c {input.decon_r2} > {output.decon_r2_gz}
		gzip -c {input.nonhumanfa_r1} > {output.nonhumanfa_r1_gz}
		gzip -c {input.nonhumanfa_r2} > {output.nonhumanfa_r2_gz}
		"""
