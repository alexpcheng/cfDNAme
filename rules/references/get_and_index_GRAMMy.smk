rule convert_reference:
	input:
		reference_fasta = 'databases/GenomeDB/NCBIGenomes06.fna'
	output:
		CT_ref_fasta = 'databases/blast/CT_conversion/NCBIGenomes06_CT.fna',
		GA_ref_fasta = 'databases/blast/GA_conversion/NCBIGenomes06_GA.fna',
		read_ids = 'databases/GenomeDB/NCBIGenomes06_readsIDs.txt'
	shell:
		"""
		python scripts/grammy_reference/reference_conversion_wgbs.py {input.reference_fasta} {output.CT_ref_fasta} {output.GA_ref_fasta} {output.read_ids}
		"""

rule gi_to_taxid:
	input:
		database_read_ids='databases/GenomeDB/NCBIGenomes06_readsIDs.txt',
		dmp='databases/GenomeDB/gi_taxid_nucl.dmp'
	output:
		taxids='databases/blast/NCBIGenomes06.gis.taxids',
		just_gi = 'databases/blast/NCBIGenomes06.gis'
	shell:
		"""
		cat {input.database_read_ids} | while read gi
		do
			taxid=$({SGREP} -n $gi {input.dmp} | cut -f 2,3 || true)
			if [ ! -z "$taxid" ]
			then
				echo -e "$gi\t$taxid" >> {output.taxids}
			fi
		done
		cut -f1 {output.taxids} > {output.just_gi}
		"""

rule index_reference:
	input:
		gi_taxid = 'databases/blast/NCBIGenomes06.gis.taxids',
		reference_fasta = 'databases/GenomeDB/NCBIGenomes06.fna'
	output:
		fai = 'databases/GenomeDB/NCBIGenomes06.fna.fai',
		lengths = 'databases/GenomeDB/taxids_lengths.txt'
	shell:
		"""
		samtools faidx {input.reference_fasta}
		Rscript scripts/grammy_reference/taxid_genome_sizes.R {output.fai} {input.gi_taxid} {output.lengths}
		"""

rule create_blast_database:
	input:
		taxids='databases/blast/NCBIGenomes06.gis.taxids',
		genome='databases/blast/{conversion}_conversion/NCBIGenomes06_{conversion}.fna',
		just_gi = 'databases/blast/NCBIGenomes06.gis'
	output:
		counts = 'databases/blast/{conversion}_conversion/NCBIGenomes06_{conversion}.fna.counts',
		obinary = 'databases/blast/{conversion}_conversion/NCBIGenomes06_{conversion}.fna.counts.obinary',
		hs_blast = 'databases/blast/{conversion}_conversion/NCBIGenomes06_{conversion}.fna.header'
	#log:
	#	'logs/databases/{conversion}_database.log'
	params:
		out_path = 'databases/blast/{conversion}_conversion/',
		genome = 'NCBIGenomes06_{conversion}.fna',
		out_prefix='NCBIGenomes06_{conversion}',
		windowmasker_f1 = 'NCBIGenomes06_{conversion}.fna.counts',
		windowmasker_f2 = 'NCBIGenomes06_{conversion}.fna.counts.obinary'
	shell:
		"""
		makeblastdb -in {input.genome} -out {params.out_path}{params.out_prefix}.fna -dbtype nucl -parse_seqids -taxid_map {input.taxids}
		cd databases/blast/{wildcards.conversion}_conversion/
		echo "Blast aliastoo"
		blastdb_aliastool -db {params.out_prefix}.fna -gilist ../../../{input.just_gi} -dbtype nucl -out {params.out_prefix}.curated
		echo "windowmasker"
		windowmasker -in {params.genome} -infmt blastdb -mk_counts -out {params.windowmasker_f1}
		windowmasker -in {params.windowmasker_f1} -sformat obinary -out {params.windowmasker_f2} -convert
		cd ../../../
		{HSBLASTN} index {input.genome}
		"""

#		"""
#		makeblastdb -in {input.genome} -out {params.out_path}{params.out_prefix} -dbtype nucl -parse_seqids -taxid_map {input.taxids} &>{log}
#		cd databases/blast/{wildcards.conversion}_conversion/
#		blastdb_aliastool -db {params.out_prefix} -gilist ../../../{input.just_gi} -dbtype nucl -out {params.out_prefix}.curated
#		cd ../../../
#		{HSBLAST} index {input.genome} &>>{log}
#		"""

rule create_grammy_database:
	input:
		genome='databases/blast/{conversion}_conversion/NCBIGenomes06_{conversion}.fna',
		taxids='databases/blast/NCBIGenomes06.gis.taxids',
		hs_blast = 'databases/blast/{conversion}_conversion/NCBIGenomes06_{conversion}.fna.header'
	output:
		prompt='logs/grammy/grammy_prompts_{conversion}',
		taxids='grefs/{conversion}/gid_tid.dmp'
	threads: indexing_threads
	params:
		genome_prefix='databases/blast/{conversion}_conversion/NCBIGenomes06_{conversion}',
		genome='NCBIGenomes06_{conversion}',
		curated='databases/blast/{conversion}_conversion/NCBIGenomes06_{conversion}.curated'
	shell:
		"""
		rm -f {output.prompt}
		mkdir -p grefs
		mkdir -p grefs/{wildcards.conversion}
		blastdbcmd -db {params.genome_prefix}.curated -entry all -outfmt "%T" | sort | uniq | while read taxid
		do
			echo "./scripts/grammy_reference/build.grefs.sh $taxid {input.taxids} {wildcards.conversion} {params.curated}" >> {output.prompt}
		done
		perl_fork_univ.pl {output.prompt} {threads}

		cp {input.taxids} {output.taxids}
		TIDS=$(cat {input.taxids} | cut -f2 | sort -u | tr '\n' ',')
		TAXIDS=${{TIDS::-1}}

		python2.7 {GRAMMY_GDT} -p 200000 -d {output.taxids} -r grefs/{wildcards.conversion} {params.genome} $TAXIDS
		mkdir -p grammy
		mv {params.genome}.gdt grammy
		"""
