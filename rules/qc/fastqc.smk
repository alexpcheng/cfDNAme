rule fqc:
	input:
		get_fastq_reads
	output:
		'sample_output/fastqc/{sample}_R1_fastqc.html',
		'sample_output/fastqc/{sample}_R2_fastqc.html'
	threads: fastqc_threads
	params:
		outdir = 'sample_output/fastqc/',
		prep_and_seq_type = get_seq_type, #first is prep, second is seq
	log: 'logs/fastqc/{sample}.fqc.log'
	shell:
		"""
		seq_type={params.prep_and_seq_type[1]}
		if [[ $seq_type == 1x* ]]
		then
			fastqc {input[0]} -t {threads} --outdir {params.outdir} &>{log}
			touch {output[1]}
		elif [[ $seq_type == 2x* ]]
		then
			fastqc {input[0]} {input[1]} -t {threads} --outdir {params.outdir} &>{log}
		fi
		"""
