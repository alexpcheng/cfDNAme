rule XY_fraction:
	input:
		mapped_all_chr='sample_output/aligned/all_chr/{sample}_mapped_all_chr.bam',
		mapped_all_chr_bai = 'sample_output/aligned/all_chr/{sample}_mapped_all_chr.bam.bai',
		bins='references/donor_fraction/bins_to_remove.txt',
		mapfile ='references/donor_fraction/500.hg19.subset.map.wig',
		gcfile ='references/donor_fraction/500.hg19.subset.gc.wig',
		removals = 'references/donor_fraction/bins_to_remove.txt'
	output:
		readcounts = temp('sample_output/donor_fraction/{sample}.auto.readcounts.wig'),
		XY = 'sample_output/donor_fraction/{sample}.XY'
	params:
		mapQ = '10',
	shell:
		"""
		samtools index {input.mapped_all_chr}
		{HMMCOPY_READCOUNTER} -w 500 -q {params.mapQ} \
			-c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
			{input.mapped_all_chr} > {output.readcounts}
		Rscript scripts/donor_fraction/XY_fraction.R {wildcards.sample} {output.readcounts} {input.mapfile} {input.gcfile} {output.XY} {input.removals}
		"""
