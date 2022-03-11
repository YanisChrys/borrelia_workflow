# Lenient filter consensus genomes where:
# no filter based on DP (more invariants called)


# 1) new mask:
rule create_lenient_mask_files:
	input:
		ref=REF_GENOM,
		my_sample="data/output/called/{sample}.g.vcf.gz",   #used for the sample name wildcard
		allvcf="data/output/called/allsites.vcf.gz"
	output:
		"data/output/called/lenientmasks/{sample}.lenientmask.vcf"
	threads:
		n_cores
	shell: """
		gatk SelectVariants \
		-R {input.ref} \
		-V {input.allvcf} \
		-O {output} \
		--remove-unused-alternates \
		--sample-name "{wildcards.sample}" \
		-select 'vc.getGenotype("{wildcards.sample}").isNoCall() || ( QUAL < 30.0 && QUAL > -1.0 )  || ( QD < 2.0 && QD > -1.0 ) || vc.isMixed()  || (vc.isSNP() &&  FS > 60.0 ) || (vc.isSNP() && MQ < 40.0 ) || (vc.isSNP() && SOR > 3.0 ) || (vc.isSNP() && MQRankSum < -12.5 ) || (vc.isSNP() && ReadPosRankSum < -8.0 ) || (vc.isIndel() && FS > 200.0 ) || (vc.isIndel() && ReadPosRankSum < -20.0 )'
	"""

# 2)
rule create_lenientmask:
	input:
		"data/output/called/lenientmasks/{sample}.lenientmask.vcf"
	output:
		bed_mask="data/output/called/lenientmasks/{sample}.editedmask.bed"
	shell: """
		vcf2bed < {input} > {output} 
	"""

# 3) create and edit the consensus sequence
# it is likely the file will have deletions marked with - spanning deletions with *
# for downstream analyses any non-bases could be a problem so we will edit those out 
rule create_lenientfasta:
	input:
		ref=REF_GENOM,
		mask="data/output/called/lenientmasks/{sample}.editedmask.bed",
		vcf="data/output/called/vcf_extracted_SNP_INDELS/{sample}.variant.vcf.gz"
	output:
		initial_fasta="data/output/called/lenientfasta/initial/{sample}.lenient.fasta",
		edited_fasta="data/output/called/lenientfasta/edited/{sample}.lenient.fasta"
	shell: """
	bcftools consensus -m {input.mask} -s {wildcards.sample} -o {output.initial_fasta} -f {input.ref} {input.vcf}
	samtools faidx {output.initial_fasta} Pbr > {output.edited_fasta}
	printf '%s\n' '%s/\*/-/g' 'x' | ex {output.edited_fasta}
	printf '%s\n' '%s/-//g' 'x' | ex {output.edited_fasta}
	sed -i "1s/.*/>Pbr_{wildcards.sample}/"  {output.edited_fasta}
	"""


# 4) create dictionary for fasta
rule create_dict_for_lenientfasta:
	input:
		"data/output/called/lenientfasta/edited/{sample}.lenient.fasta"
	output:
		"data/output/called/lenientfasta/edited/{sample}.lenient.dict"
	shell: """
		picard CreateSequenceDictionary -R {input} -O {output}
	"""


#### EXTRA #####
# CREATE FASTA WHERE ONLY CHANGES ARE
# REPLACED ON THE REFERENCE

# 5) create and edit the consensus sequence
# it is likely the file will have deletions marked with - spanning deletions with *
# for downstream analyses any non-bases could be a problem so we will edit those out 
rule create_simplefasta:
	input:
		ref=REF_GENOM,
		vcf="data/output/called/vcf_extracted_SNP_INDELS/{sample}.variant.vcf.gz"
	output:
		initial_fasta="data/output/called/simple_fasta/initial/{sample}.simple.fasta",
		edited_fasta="data/output/called/simple_fasta/edited/{sample}.simple.fasta"
	shell: """
	bcftools consensus -s {wildcards.sample} -o {output.initial_fasta} -f {input.ref} {input.vcf}
	samtools faidx {output.initial_fasta} Pbr > {output.edited_fasta}
	printf '%s\n' '%s/\*/-/g' 'x' | ex {output.edited_fasta}
	printf '%s\n' '%s/-//g' 'x' | ex {output.edited_fasta}
	sed -i "1s/.*/>Pbr_{wildcards.sample}/"  {output.edited_fasta}
	"""

