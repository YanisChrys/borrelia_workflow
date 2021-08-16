#5 Extract each sample's SNPS from the vcf
rule extract_filtered_sample_snps:
	input:
		ref = REF_GENOM,
		my_sample="data/output/called/{sample}.g.vcf.gz",   #used for the sample name wildcard
		allvcf="data/output/called/filtered/all.filtered.snps.vcf"
	output:
		"data/output/called/extracted_samples/{sample}.snp.vcf"
	threads:
		config["n_cores"]
	log:
		"logs/called/filtered/extract_snps_{sample}.log"
	shell: """
		gatk SelectVariants \
		-R {input.ref} \
		-V {input.allvcf} \
		-O {output} \
		--sample-name {wildcards.sample}
	"""



# 6) create fasta consensus genome for each sample (for each plasmid and the chromosome)
# remove ">1" from chromosome name - which is added by gatk
# and replace ":" with space so the chromosome name is recognised
rule create_snp_fasta:
	input:
		ref = REF_GENOM,
		my_sample="data/output/called/extracted_samples/{sample}.snp.vcf"
	output:
		my_fasta=temp("data/output/called/fasta/{sample}.snp.ref.fasta")
	threads:
		config["n_cores"]
	shell: """
		wait $(gatk FastaAlternateReferenceMaker \
		-R {input.ref} \
		-V {input.my_sample} \
		-O {output.my_fasta}) &&		
		printf '%s\n' '%s/>[0-9]./>/g' 'x' | ex {output.my_fasta}  &&
		printf '%s\n' '%s/:/ /g' 'x' | ex {output.my_fasta}
	"""



#6.5
# picard doesn't allow to overwrite dict files and FastaAlternateReferenceMaker
# automatically generates a wrong one
# so we move the fasta to another folder without taking the dictionaries with us and index in the new folder
rule move_fasta:
	input:
		"data/output/called/fasta/{sample}.snp.ref.fasta"
	output:
		protected("data/output/called/newref/{sample}.snp.ref.fasta")
	run:
		os.renames(input[0], output[0])


# 7) index fasta
rule index_fasta:
	input:
		my_fasta="data/output/called/newref/{sample}.snp.ref.fasta"
	output:
		"data/output/called/newref/{sample}.snp.ref.fasta.fai"
	shell: """
	samtools faidx {input.my_fasta} -o {output}
	"""

# 8) create dictionary for fasta
rule create_dict_for_fasta:
	input:
		my_fasta="data/output/called/newref/{sample}.snp.ref.fasta"
	output:
		"data/output/called/newref/{sample}.snp.ref.dict"
	shell: """
		picard CreateSequenceDictionary -R {input.my_fasta} -O {output}
	"""

#9) use new reference with indel file to extract indels using the new reference:
rule extract_samples_indels:
	input:
		ref = "data/output/called/newref/{sample}.snp.ref.fasta",
		dict="data/output/called/newref/{sample}.snp.ref.dict",
		fai="data/output/called/newref/{sample}.snp.ref.fasta.fai",
		allvcf="data/output/called/filtered/all.filtered.indels.vcf"
	output:
		"data/output/called/extracted_samples/{sample}.snp.indel.vcf"
	threads:
		config["n_cores"]
	log:
		"logs/called/filtered/extract_indels_{sample}.log"
	shell: """
		gatk SelectVariants \
		-R {input.ref} \
		-V {input.allvcf} \
		-O {output} \
		--sample-name {wildcards.sample}
	"""
	
	
# 10)
rule create_indel_fasta:
	input:
		ref = "data/output/called/newref/{sample}.snp.ref.fasta",
		my_sample="data/output/called/extracted_samples/{sample}.snp.indel.vcf"
	output:
		"data/output/called/newref/{sample}.snp.indel.ref.fasta"
	threads:
		config["n_cores"]
	shell: """
		gatk FastaAlternateReferenceMaker \
		-R {input.ref} \
		-V {input.my_sample} \
		-O {output}
	"""
	
	

	