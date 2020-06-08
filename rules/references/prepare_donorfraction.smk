rule synthetic_reads_per_chr:
	input:
		hg19CT = 'references/hg19_Bismark/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa',
		hg19GA = 'references/hg19_Bismark/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa'
	output:
		region_fileCT = temp('references/donor_fraction/{chr}/{chr}CT.txt'),
		region_fileGA = temp('references/donor_fraction/{chr}/{chr}GA.txt'),
		chr_fasta = temp('references/donor_fraction/{chr}/{chr}.fa'),
		fq1 = temp('references/donor_fraction/{chr}/synthetic_reads_{chr}_R1.fq'),
		fq2 = temp('references/donor_fraction/{chr}/synthetic_reads_{chr}_R2.fq')
	shell:
		"""
		echo {wildcards.chr}_CT_converted > {output.region_fileCT}
		echo {wildcards.chr}_GA_converted > {output.region_fileGA}
		{SEQTK} subseq {input.hg19CT} {output.region_fileCT} > {output.chr_fasta}
		{SEQTK} subseq {input.hg19GA} {output.region_fileGA} >> {output.chr_fasta}
		wgsim -d 167 -s 15 -N 10000000 -1 75 -2 75 -S 1 {output.chr_fasta} {output.fq1} {output.fq2}
		"""

rule map_synthetic_reads:
	input:
		fq1 = 'references/donor_fraction/{chr}/synthetic_reads_{chr}_R1.fq',
		fq2 = 'references/donor_fraction/{chr}/synthetic_reads_{chr}_R2.fq'
	output:
		bam = temp('references/donor_fraction/{chr}/synthetic_reads_{chr}.bam'),
		sorted_bam = 'references/donor_fraction/{chr}/synthetic_reads_{chr}_sorted.bam',
		log = 'logs/donor_fraction/{chr}.df.alignment.log'
	params:
		outdir = 'references/donor_fraction/{chr}/'
	shell:
		"""
		{BISMARK} --genome {HG19METH} \
					--parallel {threads} \
					--quiet \
					--local \
					-o {params.outdir} \
					-1 {input.fq1} \
					-2 {input.fq2}
		mv {params.outdir}synthetic_reads_{wildcards.chr}_R1_bismark_bt2_pe.bam {output.bam}
		mv {params.outdir}synthetic_reads_{wildcards.chr}_R1_bismark_bt2_PE_report.txt {output.log}
		samtools sort -@ {threads} {output.bam} > {output.sorted_bam}
		samtools index {output.sorted_bam}
		"""
rule bin_chromosomes:
	input:
		chr_lengths = 'references/hg19.chr.lengths'
	output:
		chr_length = temp('references/donor_fraction/{chr}/{chr}.length'),
		chr_binned = 'references/donor_fraction/{chr}/{chr}.binned'
	params:
		bin='500'
	shell:
		"""
		grep -w '{wildcards.chr}' {input.chr_lengths} | awk -F'\t' '{{print $1,"1",$2}}' > {output.chr_length}
		bedops --chop 499 --stagger {params.bin} -x {output.chr_length} > {output.chr_binned}
		"""
rule autosomal_binned:
	input:
		expand('references/donor_fraction/{autosomal}/{autosomal}.binned', autosomal = config['AUTOSOMALCHROMO'].split(' '))
	output:
		'references/donor_fraction/autosomal.binned'
	shell:
		"""
		cat {input} | sort-bed - > {output[0]}
		"""

rule chrY_mapping:
	input:
		sorted_bam = 'references/donor_fraction/{chr}/synthetic_reads_{chr}_sorted.bam',
		chr_binned = 'references/donor_fraction/{chr}/{chr}.binned',
		autosomal_binned = 'references/donor_fraction/autosomal.binned'
	output:
		depth_file = 'references/donor_fraction/{chr}/{chr}_depth.txt',
		depth_binned = 'references/donor_fraction/{chr}/{chr}_depth_binned.txt',
		depth_binned_cov = 'references/donor_fraction/{chr}/{chr}_depth_binned_cov.txt'
	params:
		chrY_binned = 'references/donor_fraction/chrY/chrY.binned'
	shell:
		"""
		if [ {wildcards.chr} != "chrY" ]; then
			samtools view -h {input.sorted_bam} 'chrY' | samtools depth - | awk -F'\t' '{{print $1,$2,$2+1,NR,$3}}' | sort-bed - > {output.depth_file}
			bedmap --echo --fraction-map 1 --sum --delim '\t' {params.chrY_binned} {output.depth_file} > {output.depth_binned}
			grep -v 'NAN' {output.depth_binned} > {output.depth_binned_cov} || true
		fi

		if [ {wildcards.chr} == "chrY" ]; then
			samtools view -h {input.sorted_bam} {AUTOSOMALCHROMO} | samtools sort - | samtools depth - | awk -F'\t' '{{print $1,$2,$2+1,NR,$3}}' | sort-bed - > {output.depth_file}
			bedmap --echo --fraction-map 1 --sum --delim '\t' {input.autosomal_binned} {output.depth_file} > {output.depth_binned}
			grep -v 'NAN' {output.depth_binned} > {output.depth_binned_cov} || true
		fi
		"""

rule get_removal_bins:
	input:
		expand('references/donor_fraction/{chr}/{chr}_depth_binned_cov.txt', chr = KNOWNCHROMO)
	output:
		'references/donor_fraction/bins_to_remove.txt'
	params:
		file_path = 'references/donor_fraction/',
		max_cov = '10'
	shell:
		"""
		Rscript scripts/donor_fraction/donor_fraction_bins_to_remove.R {params.file_path} {output[0]} {params.max_cov}
		"""

rule HMMcopy_reference_files:
	input:
		hg19 = 'references/hg19_Bismark/hg19.fa'
	output:
		hg19 = 'references/donor_fraction/hg19.fa',
		bw='references/donor_fraction/hg19.bw',
		mapfile ='references/donor_fraction/500.hg19.full.map.wig',
		gcfile ='references/donor_fraction/500.hg19.full.gc.wig',
		subsetmap = 'references/donor_fraction/500.hg19.subset.map.wig',
		subsetgc = 'references/donor_fraction/500.hg19.subset.gc.wig'
	shell:
		"""
		cp {input.hg19} {output.hg19}
		bowtie-build {output.hg19} {output.hg19}
		{HMMCOPY_FA_BW} -o {output.bw} {output.hg19}
		{HMMCOPY_MAP} -w 500 {output.bw} > {output.mapfile}
		{HMMCOPY_GC} -w 500 {input.hg19} > {output.gcfile}
		python scripts/donor_fraction/filter_wigs.py {output.mapfile} {output.subsetmap} {KNOWNCHROMO}
		python scripts/donor_fraction/filter_wigs.py {output.gcfile} {output.subsetgc} {KNOWNCHROMO}
		"""
		#cp {input.hg19} {output.hg19}
		#bowtie-build {output.hg19} {output.hg19}
		#{HMMCOPY_FA_BW} -o {output.bw} {output.hg19}
		#{HMMCOPY_MAP} -w 500 {output.bw} > {output.mapfile}
		#{HMMCOPY_GC} -w 500 {input.hg19} > {output.gcfile}
		#python scripts/donor_fraction/filter_wigs.py {output.mapfile} {output.subsetmap} {KNOWNCHROMO}
		#python scripts/donor_fraction/filter_wigs.py {output.gcfile} {output.subsetgc} {KNOWNCHROMO}
