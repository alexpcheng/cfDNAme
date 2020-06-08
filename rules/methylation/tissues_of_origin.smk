rule binned_methylation:
	input:
		CpG_bg='sample_output/methylation_extraction_meth{dir}/bedGraph/{sample}.bedGraph.gz',
		#golden_bed = 'references/reference_methylomes/golden_markers.bed' USED TO BE
		golden_bed = ancient('references/reference_methylomes/processed_refs/golden_markers.bed')
	output:
		CpG_bg_tmp=temp('sample_output/aligned/autosomal/{sample}_mapped_autosomal_CpG.bedGraph{dir}.tmp'),
		intersected_golden_sample_tmp = temp('sample_output/{sample}.intersect.golden{dir}.tmp'),
		binned_CpG_DMR_golden = 'sample_output/binned_samples{dir}/golden_markers/{sample}'
	shell:
		"""
		zcat {input.CpG_bg} | tail -n +2 | bedtools sort -i - > {output.CpG_bg_tmp}

		bedtools intersect -wo -a {input.golden_bed} -b {output.CpG_bg_tmp} -sorted |
			awk '$6-$5==1 {{print $0}}' | awk 'NF{{NF-=1}};1' > {output.intersected_golden_sample_tmp}
		Rscript scripts/tissues_of_origin/aggregate_over_regions.R {output.intersected_golden_sample_tmp} {output.binned_CpG_DMR_golden}
		"""

rule tissue_of_origin:
	input:
		binned_CpG_DMR_golden = 'sample_output/binned_samples{dir}/golden_markers/{sample}',
		#reference_methylomes =  'references/reference_methylomes/MethylMatrix/MethylMatrix_binned_tissues.bed',
		reference_methylomes = ancient('references/reference_methylomes/processed_refs/MethylMatrix_binned')
	output:
		mp = 'sample_output/classic_tissues_of_origin{dir}/{sample}.tsv'
	params:
		sum_to_one='FALSE',
		other='TRUE',
		group_by_celltype='FALSE',
		removals='colon_5/hsc_5/hsc_2'
	shell:
		"""
		Rscript scripts/tissues_of_origin/tissues_of_origin_v2.R \
			{input.binned_CpG_DMR_golden} \
			{input.reference_methylomes} \
			{output.mp} \
			{METHYLOME_TABLE} \
			{wildcards.sample} \
			{params.sum_to_one} \
			{params.other} \
			{params.group_by_celltype} \
			{params.removals}
		"""

rule RODEM_methylmatrix_prep:
	input:
		MM = 'references/reference_methylomes/MethylMatrix/MethylMatrix.txt'
	output:
		MM_curated = 'references/reference_methylomes/MethylMatrix/MethylMatrix_nomissing_or_outliers.txt'
	run:
		removals = ['bcell_3', 'bcell_7', 'colon_15', 'colon_17', 'bladder_1']
		with open(input.MM) as f:
			header = f.readline()
			a = [(x,y) for x,y in enumerate(header.strip().split('\t'), start=1) if y not in removals]
			new_header = [x[1] for x in a]
			indexes = [str(x[0]) for x in a]
			index_string = ','.join(indexes)
		command = f'cut -f{index_string} {input.MM} | grep -v "-" > {output.MM_curated}'
		os.system(command)

rule RODEM:
	input:
		rodem_MM = 'references/reference_methylomes/MethylMatrix/MethylMatrix_nomissing_or_outliers.txt',
		sample = 'sample_output/methylation_extraction_meth{dir}/bismark/{sample}.bismark.cov.gz'
	output:
		too = 'sample_output/updated_tissues_of_origin{dir}/{sample}.tsv'
	shell:
		"""
		python scripts/tissues_of_origin/RODEM.py --MethylMatrix {input.rodem_MM} --sample {input.sample} --output {output.too} --prefix {wildcards.sample}_{wildcards.dir}
		"""
