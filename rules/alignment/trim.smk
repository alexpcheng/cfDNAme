#rule decompress:
#	input:
#		data = CLUSTER + config['DATA']+'samples/{sample}.spring'
#	output:
#		r1 = CLUSTER+config['DATA']+'samples/{sample}_R1.fastq.gz',
#		r2 = CLUSTER+config['DATA']+'samples/{sample}_R2.fastq.gz'
#	shell:
#		"""
#		{SPRING} -d -i {input.data} -o {output.r1} {output.r2} -g
#		"""

rule trim:
	input:
		adapter_file = get_adapter_file,
		reads = get_fastq_reads
	output:
		r1p = temp('sample_output/trim/{sample}_R1_trim.fastq'),
		r2p = temp('sample_output/trim/{sample}_R2_trim.fastq')
	threads: trim_threads
	log: 'logs/trim/{sample}.trim.log'
	params:
		min_avg_phred='10',
		min_entropy='0.25',
		prep_and_seq_type = get_seq_type, #first is prep, second is seq
		mem_mb=1000
	shell:
		"""
		prep_type={params.prep_and_seq_type[0]}
		seq_type={params.prep_and_seq_type[1]}

		if [[ $seq_type == 1x* ]]
		then
			bash scripts/trim/BBDUK_WRAP.sh --bbduk_path={BBDUK} --seq_type=$seq_type --prep_type=$prep_type \
											--adapter_file={input.adapter_file} --mem_mb={params.mem_mb} --threads={threads} \
											--maq={params.min_avg_phred} --entropy={params.min_entropy} --log_file={log} \
											--R1={input.reads[0]} --R1_trim={output.r1p}
		elif [[ $seq_type == 2x* ]]
		then
			bash scripts/trim/BBDUK_WRAP.sh --bbduk_path={BBDUK} --seq_type=$seq_type --prep_type=$prep_type \
											--adapter_file={input.adapter_file} --mem_mb={params.mem_mb} --threads={threads} \
											--maq={params.min_avg_phred} --entropy={params.min_entropy} --log_file={log} \
											--R1={input.reads[0]} --R1_trim={output.r1p} \
											--R2={input.reads[1]} --R2_trim={output.r2p}
		fi
		"""
