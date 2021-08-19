#5 Extract each sample's unique SNPS (exclude ref calls) from the vcf
#used for creating fasta with mask
rule extract_filtered_sample_snps:
	input:
		ref = REF_GENOM,
		my_sample="data/output/called/{sample}.g.vcf.gz",   #used for the sample name wildcard
		allvcf="data/output/called/filtered/allsites.filtered.snps.vcf"
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
		- select 'vc.getGenotype({wildcards.sample}).isHomVar() && vc.getGenotype({wildcards.sample}).isCalled()' \
		--use-iupac-sample \
		--sample-name {wildcards.sample}
	"""



#create mask with uncalled sites, all filters and low depth sites
#anything else will be a good quality variant
# 	-select 'DP < 2 && vc.getGenotype("RB-NU5H-MAY19-167").isHomVar()'
rule create_mask_files:
	input:
		ref = REF_GENOM,
		my_sample="data/output/called/{sample}.g.vcf.gz",   #used for the sample name wildcard
		allvcf="data/output/called/allsites.vcf"
	output:
		"data/output/called/snpmasks/{sample}.mask.vcf"
	threads:
		config["n_cores"]
	log:
		"logs/called/snpmasks/extract_snps_{sample}.log"
	shell: """
		gatk SelectVariants \ 
		-R {input.ref} \
		-V {input.allvcf} \
		-O {output} \
		--sample-name {wildcards.sample} \
		-select 'vc.getGenotype({wildcards.sample}).isNoCall() || DP < 2.0 || ((QD < 2.0 && vc.isSNP()) || (FS > 60.0 && vc.isSNP()) || (MQ < 40.0 && vc.isSNP()) || (SOR > 3.0 && vc.isSNP()) || (MQRankSum < -12.5 && vc.isSNP()) || (ReadPosRankSum < -8.0 && vc.isSNP()) || (QUAL < 30.0 && vc.isSNP()) ) || 
		((QD < 2.0  && vc.isIndel() ) || (FS > 200.0  && vc.isIndel() ) || (ReadPosRankSum < -20.0  && vc.isIndel()) || (QUAL < 30.0  && vc.isIndel()) )'
	"""

# 6) create fasta consensus genome for each sample (for each plasmid and the chromosome)
#use mask to turn filtered sites to N and give priority to the mask instead of the file
# remove ">1" from chromosome name - which is added by gatk
# and replace ":" with space so the chromosome name is recognised
rule create_snp_fasta:
	input:
		ref = REF_GENOM,
		my_sample="data/output/called/extracted_samples/{sample}.snp.vcf",
		my_mask="data/output/called/snpmasks/{sample}.mask.vcf"
	output:
		my_fasta="data/output/called/fasta/{sample}.snp.ref.fasta"
	threads:
		config["n_cores"]
	shell: """
		wait $(gatk FastaAlternateReferenceMaker \
		-R {input.ref} \
		-V {input.my_sample} \
		--line-width 80 \
		--snp-mask {input.my_mask} \
		--snp-mask-priority \
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

#9) use new reference with indel file to extract indels that are called as alts:
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
		- select 'vc.getGenotype({wildcards.sample}).isHomVar() && vc.getGenotype({wildcards.sample}).isCalled()' \
		--use-iupac-sample \
		--sample-name {wildcards.sample}
	"""


#create mask with uncalled sites, all filters and low depth sites
#anything else will be a good quality variant
# 	-select 'DP < 2 && vc.getGenotype("RB-NU5H-MAY19-167").isHomVar()'
rule create_mask_files:
	input:
		ref = "data/output/called/newref/{sample}.snp.ref.fasta",
		my_sample="data/output/called/{sample}.g.vcf.gz",   #used for the sample name wildcard
		allvcf="data/output/called/allsites.vcf"
	output:
		"data/output/called/indelmasks/{sample}.mask.vcf"
	threads:
		config["n_cores"]
	log:
		"logs/called/masks/extract_snps_{sample}.log"
	shell: """
		gatk SelectVariants \ 
		-R {input.ref} \
		-V {input.allvcf} \
		-O {output} \
		--sample-name {wildcards.sample} \
		-select 'vc.getGenotype({wildcards.sample}).isNoCall() || DP < 2.0 || 
		((QD < 2.0 && vc.isSNP()) || (FS > 60.0 && vc.isSNP()) || (MQ < 40.0 && vc.isSNP()) || (SOR > 3.0 && vc.isSNP()) || (MQRankSum < -12.5 && vc.isSNP()) || (ReadPosRankSum < -8.0 && vc.isSNP()) || (QUAL < 30.0 && vc.isSNP()) ) || ((QD < 2.0  && vc.isIndel() ) || (FS > 200.0  && vc.isIndel() ) || (ReadPosRankSum < -20.0  && vc.isIndel()) || (QUAL < 30.0  && vc.isIndel()) )'
	"""


# 10)
rule create_indel_fasta:
	input:
		ref = "data/output/called/newref/{sample}.snp.ref.fasta",
		my_sample="data/output/called/extracted_samples/{sample}.snp.indel.vcf",
		my_mask="data/output/called/indelmasks/{sample}.mask.vcf"
	output:
		"data/output/called/newref/{sample}.snp.indel.ref.fasta"
	threads:
		config["n_cores"]
	shell: """
		gatk FastaAlternateReferenceMaker \
		-R {input.ref} \
		-V {input.my_sample} \
		--line-width 80 \
		--snp-mask {input.my_mask} \
		--snp-mask-priority \
		-O {output}  &&		
		printf '%s\n' '%s/>[0-9]./>/g' 'x' | ex {output}  &&
		printf '%s\n' '%s/:/ /g' 'x' | ex {output}
	"""

# 11) move second reference
rule move_snpindel_ref:
	input:
		"data/output/called/newref/{sample}.snp.indel.ref.fasta"
	output:
		protected("data/output/called/snp_indel_ref/{sample}.snp.indel.ref.fasta")
	run:
		os.renames(input[0], output[0])
	
# 7) index 2nd fasta
rule index_2nd_fasta:
	input:
		my_fasta="data/output/called/snp_indel_ref/{sample}.snp.indel.ref.fasta"
	output:
		"data/output/called/snp_indel_ref/{sample}.snp.indel.ref.fasta.fai"
	shell: """
		samtools faidx {input.my_fasta} -o {output}
	"""

# 8) create dictionary for 2nd fasta
rule create_dict_for_2nd_fasta:
	input:
		my_fasta="data/output/called/snp_indel_ref/{sample}.snp.indel.ref.fasta"
	output:
		"data/output/called/snp_indel_ref/{sample}.snp.indel.ref.dict"
	shell: """
		picard CreateSequenceDictionary -R {input.my_fasta} -O {output}
	"""
	