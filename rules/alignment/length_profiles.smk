rule length_profiles:
	input:
		mapped_autosomal= 'sample_output/aligned/autosomal/{sample}_mapped_autosomal.bam'
	output:
		lengths='sample_output/Lengths/{sample}_aligned.lengths.gz',
		length_counts = 'sample_output/Lengths/{sample}_aligned.lengths.counts',
		pdf='sample_output/Lengths/{sample}_FragsHistogram.pdf'
	params:
		prep_and_seq_type = get_seq_type
	shell:
		"""
		seq_type={params.prep_and_seq_type[1]}

		if [[ $seq_type == 2x* ]]; then
			samtools view {input.mapped_autosomal} | awk -v OFS='\t' '{{if ($9>0 && $9<1000) print $1,$9}}' | gzip -9 > {output.lengths}
        	zcat {output.lengths} | cut -f2 | sort | uniq -c > {output.length_counts}
			Rscript scripts/length_profiles/HistogramFragmentLengths.R {output.length_counts} {output.pdf}
		elif [[ $seq_type == 1x* ]]; then
			touch {output.lengths} {output.length_counts} {output.pdf}
		fi
		"""
