
# new referefences where:
# no filter based on DP (more invariants called)
# and note deletions and insertions

# new mask:
rule create_lenient_mask_files:
	input:
		ref=REF_GENOM,
		my_sample="data/output/called/{sample}.g.vcf.gz",   #used for the sample name wildcard
		allvcf="data/output/called/allsites.vcf.gz"
	output:
		"data/output/called/leneintmasks/{sample}.lenientmask.vcf"
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
		-select 'vc.getGenotype("{wildcards.sample}").isNoCall() || ( QUAL < 30.0 && QUAL > -1.0 )  || ( QD < 2.0 && QD > -1.0 ) || vc.isMixed()  || (vc.isSNP() &&  FS > 60.0 ) || (vc.isSNP() && MQ < 40.0 ) || (vc.isSNP() && SOR > 3.0 ) || (vc.isSNP() && MQRankSum < -12.5 ) || (vc.isSNP() && ReadPosRankSum < -8.0 ) || (vc.isIndel() && FS > 200.0 ) || (vc.isIndel() && ReadPosRankSum < -20.0 )'
	"""


rule bed_lenientmask:
	input:
		"data/output/called/leneintmasks/{sample}.lenientmask.vcf"
	output:
		bed_mask="data/output/called/leneintmasks/{sample}.editedmask.bed"
	threads:
		config["n_cores"]
	shell: """
		vcf2bed < {input} > {output} 
	"""

rule create_lenientfasta:
	input:
		ref=REF_GENOM,
		mask="data/output/called/leneintmasks/{sample}.editedmask.bed",
		vcf="data/output/called/vcf_extracted_SNP_INDELS/{sample}.variant.vcf.gz"
	output:
		"data/output/called/leneintfasta/{sample}.lenient.fasta"
	shell: """
	bcftools consensus --mark-del "-" --mark-ins lc -s {wildcards.sample} -m {input.mask} -o {output} -f {input.ref} {input.vcf}
	"""


# 7) create dictionary for fasta
rule create_dict_for_lenientfasta:
	input:
		"data/output/called/leneintfasta/{sample}.lenient.fasta"
	output:
		"data/output/called/leneintfasta/{sample}.lenient.dict"
	shell: """
		picard CreateSequenceDictionary -R {input} -O {output}
	"""
