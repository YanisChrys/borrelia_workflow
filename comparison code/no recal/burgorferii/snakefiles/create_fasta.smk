# This subprocess filters SNPs and INDELS separately, merges the files and uses a "mask" file to mark "bad" sites as N
# then creates the consensus genome fasta.
# Filtering is strict, same as the stringent_bed subprocess

# 1) Extract each samples called SNPs and Indels
rule extract_sample_snp_indels:
	input:
		ref=REF_GENOM,
		my_sample="data/output/called/{sample}.g.vcf.gz",   #used for the sample name wildcard
		allvcf="data/output/called/filtered/allsites.filtered.vcf"
	output:
		"data/output/called/vcf_extracted_SNP_INDELS/{sample}.variant.vcf.gz"
	threads:
		n_cores
	log:
		"logs/called/vcf_extracted_SNP_INDELS/extract_snps_{sample}.log"
	shell: """
		gatk SelectVariants \
		-R {input.ref} \
		-V {input.allvcf} \
		-O {output} \
		--exclude-filtered \
		--remove-unused-alternates \
		--sample-name "{wildcards.sample}" \
		-select 'vc.isSNP() || vc.isIndel()'
	"""



# 2) create stringent version mask with: 
# uncalled sites, 
# mixed sites (not handled well - or at all which means they become false positive REF calls)
# all sites(SNP or Indel) not passing any one of our filters and 
# low depth sites (AD of ALT-called site is less than 2 or DP of REF-called site is less than 2)
# anything else will be a good quality, called site
rule create_mask_files:
	input:
		ref=REF_GENOM,
		my_sample="data/output/called/{sample}.g.vcf.gz",   #used for the sample name wildcard
		allvcf="data/output/called/allsites.vcf.gz"
	output:
		"data/output/called/stringent_masks/{sample}.mask.vcf"
	threads:
		n_cores
	log:
		"logs/called/snpmasks/extract_snps_{sample}.log"
	shell: """
		gatk SelectVariants \
		-R {input.ref} \
		-V {input.allvcf} \
		-O {output} \
		--remove-unused-alternates \
		--set-filtered-gt-to-nocall \
		--sample-name "{wildcards.sample}" \
		-select 'vc.getGenotype("{wildcards.sample}").isNoCall() || DP < 2.0  || ( QUAL < 30.0 && QUAL > -1.0 )  || ( QD < 2.0 && QD > -1.0 ) || vc.isMixed() || (vc.getGenotype("{wildcards.sample}").isHomRef() && vc.getGenotype("{wildcards.sample}").hasAD() && vc.getGenotype("{wildcards.sample}").getAD().0 < 2.0 )  || (vc.getGenotype("{wildcards.sample}").isHomVar() && vc.getGenotype("{wildcards.sample}").hasAD() && vc.getGenotype("{wildcards.sample}").getAD().1 < 2.0 )  || (vc.isSNP() &&  FS > 60.0 ) || (vc.isSNP() && MQ < 40.0 ) || (vc.isSNP() && SOR > 3.0 ) || (vc.isSNP() && MQRankSum < -12.5 ) || (vc.isSNP() && ReadPosRankSum < -8.0 ) || (vc.isIndel() && FS > 200.0 ) || (vc.isIndel() && ReadPosRankSum < -20.0 )'
	"""
	
# 3) Edit mask files
# add n to mask for bcftools and create bed file used by bcftools
# NOTE: awk is for snakemake, not directly portable to terminal

rule edit_mask_files:
	input:
		"data/output/called/stringent_masks/{sample}.mask.vcf"
	output:
		vcf_mask="data/output/called/stringent_masks/editedmasks/{sample}.editedmask.vcf",
		bed_mask="data/output/called/stringent_masks/editedmasks/{sample}.editedmask.bed"
	threads:
		n_cores
	shell: """
		awk 'BEGIN{{OFS="\t"}} $1 ~ /^#/ {{print $0;next}} {{ $5 = "N";print }}' {input} > {output.vcf_mask} &&
		gatk IndexFeatureFile -I {output.vcf_mask}
		vcf2bed < {output.vcf_mask} > {output.bed_mask} 
	"""

# 4)
rule create_fasta:
	input:
		ref=REF_GENOM,
		mask="data/output/called/stringent_masks/editedmasks/{sample}.editedmask.bed",
		vcf="data/output/called/vcf_extracted_SNP_INDELS/{sample}.variant.vcf.gz"
	output:
		"data/output/called/stringent_fasta/initial/{sample}.fasta"
	shell: """
		bcftools consensus -m {input.mask} -o {output} -f {input.ref} {input.vcf}
	"""

# 5) exclude all plasmids and write sample name to the chromosome name
#also replace asterisks with dashes
rule edit_fasta:
	input:
		"data/output/called/stringent_fasta/initial/{sample}.fasta"
	output:
		"data/output/called/stringent_fasta/edited/{sample}.fasta"
	shell: """
		samtools faidx {input} NC_011728.1 > {output}
		printf '%s\n' '%s/\*/-/g' 'x' | ex {output}
		printf '%s\n' '%s/-//g' 'x' | ex {output}
		sed -i "1s/.*/>NC_011728.1_{wildcards.sample}/"  {output}
	"""

# 6) create dictionary for fasta
rule create_dict_for_fasta:
	input:
		"data/output/called/stringent_fasta/edited/{sample}.fasta"
	output:
		"data/output/called/stringent_fasta/edited/{sample}.dict"
	shell: """
		picard CreateSequenceDictionary -R {input} -O {output}
	"""

# 7) save all sites of a sample to check the fasta later
rule vcf_to_check_fasta:
	input:
		ref=REF_GENOM,
		my_sample="data/output/called/{sample}.g.vcf.gz",   #used for the sample name wildcard
		allvcf="data/output/called/allsites.vcf.gz"
	output:
		"data/output/called/vcf_allsites/{sample}.allsites.vcf"
	threads:
		n_cores
	shell: """
		gatk SelectVariants \
		-R {input.ref} \
		-V {input.allvcf} \
		-O {output} \
		--remove-unused-alternates \
		--sample-name "{wildcards.sample}" \
	"""


#### EXTRA #####
# CREATE FASTA WHERE ONLY CHANGES ARE
# REPLACED ON THE REFERENCE

# only keep "good" variants to make the fasta file
rule extract_stringent_snp_indels:
	input:
		ref=REF_GENOM,
		my_sample="data/output/called/{sample}.g.vcf.gz",   #used for the sample name wildcard
		vcf="data/output/called/vcf_extracted_SNP_INDELS/{sample}.variant.vcf.gz"
	output:
		"data/output/called/stringent_snpsindels/{sample}.str_snps.vcf.gz"
	threads:
		n_cores
	log:
		"logs/called/stringent_snpsindels/str_snps_{sample}.log"
	shell: """
		gatk SelectVariants \
		-R {input.ref} \
		-V {input.vcf} \
		-O {output} \
		--remove-unused-alternates \
		--sample-name "{wildcards.sample}" \
		-select '!(vc.getGenotype("{wildcards.sample}").isNoCall() || DP < 2.0  || ( QUAL < 30.0 && QUAL > -1.0 )  || ( QD < 2.0 && QD > -1.0 ) || vc.isMixed() || (vc.getGenotype("{wildcards.sample}").isHomRef() && vc.getGenotype("{wildcards.sample}").hasAD() && vc.getGenotype("{wildcards.sample}").getAD().0 < 2.0 )  || (vc.getGenotype("{wildcards.sample}").isHomVar() && vc.getGenotype("{wildcards.sample}").hasAD() && vc.getGenotype("{wildcards.sample}").getAD().1 < 2.0 )  || (vc.isSNP() &&  FS > 60.0 ) || (vc.isSNP() && MQ < 40.0 ) || (vc.isSNP() && SOR > 3.0 ) || (vc.isSNP() && MQRankSum < -12.5 ) || (vc.isSNP() && ReadPosRankSum < -8.0 ) || (vc.isIndel() && FS > 200.0 ) || (vc.isIndel() && ReadPosRankSum < -20.0 ) )'
	"""

rule create_simple_string_fasta:
	input:
		ref=REF_GENOM,
		vcf="data/output/called/stringent_snpsindels/{sample}.str_snps.vcf.gz"
	output:
		initial_fasta="data/output/called/simple_string_fasta/initial/{sample}.simple.fasta",
		edited_fasta="data/output/called/simple_string_fasta/edited/{sample}.simple.fasta"
	shell: """
	bcftools consensus -s {wildcards.sample} -o {output.initial_fasta} -f {input.ref} {input.vcf}
	samtools faidx {output.initial_fasta} NC_011728.1 > {output.edited_fasta}
	printf '%s\n' '%s/\*/-/g' 'x' | ex {output.edited_fasta}
	printf '%s\n' '%s/-//g' 'x' | ex {output.edited_fasta}
	sed -i "1s/.*/>NC_011728.1_{wildcards.sample}/"  {output.edited_fasta}
	"""
