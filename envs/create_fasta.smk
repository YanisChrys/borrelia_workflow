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
rule create_snp_fasta:
	input:
		ref = REF_GENOM,
		my_sample="data/output/called/extracted_samples/{sample}.snp.vcf"
	output:
		"data/output/called/fasta/{sample}.snp.ref.fasta"
	threads:
		config["n_cores"]
	log:
		"logs/called/filtered/create_snp_{sample}.fasta.log"
	shell: """
		gatk FastaAlternateReferenceMaker \
		-R {input.ref} \
		-V {input.my_sample} \
		-O {output}
	"""
	
# 7) index fasta
rule index_fasta:
	input:
		"data/output/called/fasta/{sample}.snp.ref.fasta"
	output:
		"data/output/called/fasta/{sample}.snp.ref.dict"
	shell:
		"picard CreateSequenceDictionary R={input} O={output}"



#8) use new reference with indel file:
rule extract_samples_indels:
	input:
		ref = "data/output/called/fasta/{sample}.snp.ref.fasta",
		dict="data/output/called/fasta/{sample}.snp.ref.dict",
		allvcf="data/output/called/filtered/all.filtered.indels.vcf"
	output:
		"data/output/called/extracted_samples/{sample}.vcf"
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
	
	
	
rule create_indel_fasta:
	input:
		ref = "data/output/called/fasta/{sample}.snp.ref.fasta",
		my_sample="data/output/called/extracted_samples/{sample}.vcf"
	output:
		"data/output/called/fasta/{sample}.snp.indel.ref.fasta"
	threads:
		config["n_cores"]
	shell: """
		gatk FastaAlternateReferenceMaker \
		-R {input.ref} \
		-V {input.my_sample} \
		-O {output}
	"""
	
	

	