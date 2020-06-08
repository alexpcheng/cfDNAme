rule alignment:
	input:
		r1p = 'sample_output/trim/{sample}_R1_trim.fastq',
		r2p = 'sample_output/trim/{sample}_R2_trim.fastq',
		hg19 = 'references/hg19_Bismark/hg19.fa',
		mic = 'references/controls_Bismark/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa'
	output:
		bam = 'sample_output/aligned/raw_aligned/{sample}.bam',
		unmapped_R1 = 'sample_output/aligned/unmapped/{sample}_pe_unmapped_R1.fastq.gz',
		unmapped_R2 = 'sample_output/aligned/unmapped/{sample}_pe_unmapped_R2.fastq.gz'
	log: 'logs/alignment/{sample}.alignment.log'
	threads: alignment_threads
	params:
		outdir = 'sample_output/aligned/raw_aligned/',
		prep_and_seq_type = get_seq_type #first is prep, second is seq
	shell:
		"""
		seq_type={params.prep_and_seq_type[1]}
		threads=$(( {threads} / 2 ))


		if [[ {wildcards.sample} == *MCB* ]]
		then
			METHREF={CTLMETH}
		else
			METHREF={HG19METH}
		fi


		if [[ $seq_type == 2x* ]]
		then
			{BISMARK} --genome $METHREF \
					--parallel $threads \
  					--quiet \
                    --unmapped \
  					-o {params.outdir} \
  					-1 {input.r1p} \
  					-2 {input.r2p}
  			mv {params.outdir}{wildcards.sample}_R1_trim_bismark_bt2_pe.bam {output.bam}
  			mv {params.outdir}{wildcards.sample}_R1_trim_bismark_bt2_PE_report.txt {log}

			mv {params.outdir}{wildcards.sample}_R1_trim.fastq_unmapped_reads_1.fq.gz {output.unmapped_R1}
  			mv {params.outdir}{wildcards.sample}_R2_trim.fastq_unmapped_reads_2.fq.gz {output.unmapped_R2}
		elif [[ $seq_type == 1x* ]]
		then
			{BISMARK} --genome $METHREF \
  					--parallel $threads \
  					--quiet \
                    --unmapped \
  					-o {params.outdir} \
  					{input.r1p}

  			mv {params.outdir}{wildcards.sample}_R1_trim_bismark_bt2.bam {output.bam}
			mv {params.outdir}{wildcards.sample}_R1_trim_bismark_bt2_SE_report.txt {log}

			mv {params.outdir}{wildcards.sample}_R1_trim.fastq_unmapped_reads.fq.gz {output.unmapped_R1}
			echo "This was processed as a single-end dataset so this is just a fun placeholder file" > {output.unmapped_R2}
		fi
		"""


rule filter_bam:
	input:
		bam ='sample_output/aligned/raw_aligned/{sample}.bam'
	output:
		bismark_dup = temp('sample_output/aligned/raw_aligned/{sample}.deduplicated.bam'),
		mapped_all_chr=temp('sample_output/aligned/all_chr/{sample}_mapped_all_chr.bam'),
		mapped_all_chr_bai = temp('sample_output/aligned/all_chr/{sample}_mapped_all_chr.bam.bai'),
		mapped_autosomal=temp('sample_output/aligned/autosomal/{sample}_mapped_autosomal.bam'),
		mapped_autosomal_bai = temp('sample_output/aligned/autosomal/{sample}_mapped_autosomal.bam.bai'),
		name_sorted = temp('sample_output/aligned/autosomal/{sample}_mapped_autosomal_namesorted.bam')
	threads: filter_bam_threads
	log: 'logs/deduplication/{sample}.log'
	params:
		mapQ='10',
		outdir = 'sample_output/aligned/raw_aligned/'
	shell:
		"""
		if [[ {wildcards.sample} == *MCB* ]]
		then
			cp {input.bam} {output.bismark_dup}
			cp {output.bismark_dup} {output.mapped_all_chr}
			samtools index {output.mapped_all_chr}
			cp {output.mapped_all_chr} {output.mapped_autosomal}
			samtools index {output.mapped_autosomal}
			samtools sort -@ {threads} -n {output.mapped_autosomal} -o {output.name_sorted}
		else
			if [[ {wildcards.sample} == *L004* ]] || [[ {wildcards.sample} == *L003* ]]
			then
				extra_params=""
			else
				extra_params="--barcode"
			fi
			{RMDUPS} -o {wildcards.sample} --output_dir {params.outdir} --bam {input.bam} $extra_params &> {log}
			samtools sort -@ {threads} {output.bismark_dup} -o {output.mapped_all_chr}
			samtools index {output.mapped_all_chr}

			samtools view -@ {threads} -b -q {params.mapQ} {output.mapped_all_chr} {AUTOSOMALCHROMO} -o - | samtools sort -@ {threads} - -o {output.mapped_autosomal}
			samtools index {output.mapped_autosomal}

			samtools sort -@ {threads} -n {output.mapped_autosomal} -o {output.name_sorted}
		fi
		"""
