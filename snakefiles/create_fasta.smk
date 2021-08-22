#1) Extract each samples called SNPs and Indels to use with FastaAlternateReferenceMaker
# the good quality ones will be replaced in the reference
# and the bad quality ones will be masked by the mask file containing all 
# unwanted sites
# the snp and indel based filtering has been done, now remove sites with AD<2
rule extract_sample_snp_indels:
	input:
		ref=REF_GENOM,
		my_sample="data/output/called/{sample}.g.vcf.gz",   #used for the sample name wildcard
		allvcf="data/output/called/filtered/allsites.filtered.vcf"
	output:
		"data/output/called/vcf_extracted_SNP_INDELS/{sample}.variant.vcf"
	threads:
		config["n_cores"]
	log:
		"logs/called/vcf_extracted_SNP_INDELS/extract_snps_{sample}.log"
	shell: """
		gatk SelectVariants \
		-R {input.ref} \
		-V {input.allvcf} \
		-O {output} \
		--remove-unused-alternates \
		--exclude-filtered \
		--sample-name "{wildcards.sample}" \
		-select 'vc.getGenotype("{wildcards.sample}").hasAD() && vc.isVariant() && (vc.getGenotype("{wildcards.sample}").getAD().1 > 2.0)'
	"""



# 2) create mask with: 
# uncalled sites, 
# all sites(SNP or Indel) not passing any one of our filters and 
# low depth sites (AD of ALT-called site is less than 2 or DP of REF-called site is less than 2)
# anything else will be a good quality, called site
# also, get anything that isn't a SNP or an indel because they will be ignored by the tool and be entered as a reference or worse
rule create_mask_files:
	input:
		ref=REF_GENOM,
		my_sample="data/output/called/{sample}.g.vcf.gz",   #used for the sample name wildcard
		allvcf="data/output/called/allsites.vcf"
	output:
		"data/output/called/masks/{sample}.mask.vcf"
	threads:
		config["n_cores"]
	log:
		"logs/called/snpmasks/extract_snps_{sample}.log"
	shell: """
		gatk SelectVariants \
		-R {input.ref} \
		-V {input.allvcf} \
		-O {output} \
		--remove-unused-alternates \
		--sample-name "{wildcards.sample}" \
		-select 'vc.getGenotype("{wildcards.sample}").isNoCall() || DP < 2.0 || vc.isMixed() || (vc.getGenotype("{wildcards.sample}").isCalled() && !vc.isVariant() && vc.getGenotype("{wildcards.sample}").hasAD() && vc.getGenotype("{wildcards.sample}").getAD().0 < 2.0 ) || (vc.getGenotype("{wildcards.sample}").isHomVar() && vc.getGenotype("{wildcards.sample}").getAD().1 < 2.0 ) || ( vc.isSNP() && ( (QD < 2.0 ) || (FS > 60.0 ) || (MQ < 40.0 ) || (SOR > 3.0 ) || (MQRankSum < -12.5 ) || (ReadPosRankSum < -8.0 ) || (QUAL < 30.0 ) )) || ( vc.isIndel() && ( (QD < 2.0 ) || (FS > 200.0 ) || (ReadPosRankSum < -20.0 ) || (QUAL < 30.0 ) ))'
	"""

# 3) Edit mask files
# add "N" to all ALT positions in mask file to trick FastaAlternateReferenceMaker
# into thinking they are variants and actually mask non-variant called or uncalled sites 
# otherwise it will only mask unwanted variant sites
# NOTE: awk is for snakemake, not directly portable to terminal
rule edit_mask_files:
	input:
		"data/output/called/masks/{sample}.mask.vcf"
	output:
		"data/output/called/editedmasks/{sample}.editedmask.vcf"
	threads:
		config["n_cores"]
	shell: """
	awk 'BEGIN{{OFS="\t"}} $1 ~ /^#/ {{print $0;next}} {{ $5 = "N";print }}' {input} > {output} &&
	gatk IndexFeatureFile -I {output}
	"""



# 4) Create fasta consensus genome for each sample (for each plasmid and the chromosome)
# for SNPs and simple INDELS simulataneously:
# use mask to turn unwanted sites to N and give priority to the mask instead of the file
# so if there is a variant in the vcf file that is unwanted (very likely, vcf file is not filtered)
# that site will also be replaced with "N"
# then:
# printf: remove ">1" from chromosome name - which is added by gatk
# and replace ":" with space so the chromosome name is recognised
rule create_fasta:
	input:
		ref=REF_GENOM,
		my_sample="data/output/called/vcf_extracted_SNP_INDELS/{sample}.variant.vcf",
		my_mask="data/output/called/editedmasks/{sample}.editedmask.vcf"
	output:
		"data/output/called/fasta/{sample}.newref.fasta"
	threads:
		config["n_cores"]
	shell: """
		gatk FastaAlternateReferenceMaker \
		-R {input.ref} \
		-V {input.my_sample} \
		-O {output} \
		--line-width 80 \
		--snp-mask {input.my_mask} \
		--snp-mask-priority &&
		printf '%s\n' '%s/>[0-9]./>/g' 'x' | ex {output}  &&
		printf '%s\n' '%s/:/ /g' 'x' | ex {output}
	"""



#5) Move fasta
# picard doesn't allow to overwrite dict files and FastaAlternateReferenceMaker
# automatically generates a "wrong" one
# so we move the fasta to another folder without taking the dictionaries with us and index in the new folder
rule move_fasta:
	input:
		"data/output/called/fasta/{sample}.newref.fasta"
	output:
		protected("data/output/called/newref/{sample}.newref.fasta")
	run:
		os.renames(input[0], output[0])


# 6) Index fasta
rule index_fasta:
	input:
		"data/output/called/newref/{sample}.newref.fasta"
	output:
		"data/output/called/newref/{sample}.newref.fasta.fai"
	shell: """
		samtools faidx {input} -o {output}
	"""

# 7) create dictionary for fasta
rule create_dict_for_fasta:
	input:
		"data/output/called/newref/{sample}.newref.fasta"
	output:
		"data/output/called/newref/{sample}.newref.dict"
	shell: """
		picard CreateSequenceDictionary -R {input} -O {output}
	"""


	