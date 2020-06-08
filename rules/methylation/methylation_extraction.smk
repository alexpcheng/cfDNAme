rule methylation_extraction:
	input:
		mapped_autosomal= 'sample_output/aligned/autosomal/{sample}_mapped_autosomal.bam',
		name_sorted = 'sample_output/aligned/autosomal/{sample}_mapped_autosomal_namesorted.bam'
	output:
		CpG_bg='sample_output/methylation_extraction{dir}/bedGraph/{sample}.bedGraph.gz',
		CpG_bg_gz = temp('sample_output/methylation_extraction{dir}/{sample}_mapped_autosomal_namesorted.bedGraph.gz'),
		CpG_bismark = 'sample_output/methylation_extraction{dir}/bismark/{sample}.bismark.cov.gz',
		mbias = 'sample_output/methylation_extraction{dir}/mbias/{sample}.M-bias.txt',
		CHGOB = 'sample_output/methylation_extraction{dir}/CHG/{sample}_CHG_OB.txt.gz',
		CHGOT = 'sample_output/methylation_extraction{dir}/CHG/{sample}_CHG_OT.txt.gz',
		CHHOB = 'sample_output/methylation_extraction{dir}/CHH/{sample}_CHH_OB.txt.gz',
		CHHOT = 'sample_output/methylation_extraction{dir}/CHH/{sample}_CHH_OT.txt.gz',
		CpGOB = 'sample_output/methylation_extraction{dir}/CpG/{sample}_CpG_OB.txt.gz',
		CpGOT = 'sample_output/methylation_extraction{dir}/CpG/{sample}_CpG_OT.txt.gz',
		log = 'logs/methylation_extraction{dir}/{sample}.methylation_extraction.log'
	threads: methylation_extraction_threads
	params:
		outdir = 'sample_output/methylation_extraction{dir}/',
		prep_and_seq_type = get_seq_type
	shell:
		"""
		seq_type={params.prep_and_seq_type[1]}
		extra_params="--gzip"
		if [[ {wildcards.dir} == "_meth_trim" ]]; then

			extra_params=$extra_params" --ignore 10"

			if [[ $seq_type == 2x* ]]; then
				extra_params=$extra_params" --ignore_r2 5"
			fi

		fi
		echo $extra_params
		{METHEXT} --parallel {threads} \
					--bedGraph \
					-o {params.outdir} \
					$extra_params \
					{input.name_sorted}
		gunzip -c {params.outdir}{wildcards.sample}_mapped_autosomal_namesorted.bedGraph.gz | sort-bed - | gzip > {output.CpG_bg}
		mv {params.outdir}{wildcards.sample}_mapped_autosomal_namesorted.bismark.cov.gz {output.CpG_bismark}
		mv {params.outdir}{wildcards.sample}_mapped_autosomal_namesorted_splitting_report.txt {output.log}

		mv {params.outdir}{wildcards.sample}_mapped_autosomal_namesorted.M-bias.txt {output.mbias}
		mv {params.outdir}CHG_OT_{wildcards.sample}_mapped_autosomal_namesorted.txt.gz {output.CHGOT}
		mv {params.outdir}CHG_OB_{wildcards.sample}_mapped_autosomal_namesorted.txt.gz {output.CHGOB}
		mv {params.outdir}CHH_OT_{wildcards.sample}_mapped_autosomal_namesorted.txt.gz {output.CHHOT}
		mv {params.outdir}CHH_OB_{wildcards.sample}_mapped_autosomal_namesorted.txt.gz {output.CHHOB}
		mv {params.outdir}CpG_OT_{wildcards.sample}_mapped_autosomal_namesorted.txt.gz {output.CpGOT}
		mv {params.outdir}CpG_OB_{wildcards.sample}_mapped_autosomal_namesorted.txt.gz {output.CpGOB}
		"""

rule mbias:
	input:
		mbias_untrimmed = 'sample_output/methylation_extraction_meth_untrim/mbias/{sample}.M-bias.txt',
		mbias_trimmed = 'sample_output/methylation_extraction_meth_trim/mbias/{sample}.M-bias.txt'
	output:
		pdf= 'sample_output/mbias/{sample}.pdf'
	params:
		prep_and_seq_type = get_seq_type
	shell:
		"""
		seq_type={params.prep_and_seq_type[0]}
		cat {input.mbias_trimmed} {input.mbias_untrimmed} | \
			grep -P '\t' | grep -v position | \
			Rscript scripts/methylation_extraction/plot_mbias.R $seq_type {output.pdf}
		"""
