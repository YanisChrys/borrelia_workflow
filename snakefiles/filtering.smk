
rule extract_SNPS:
	input:
		allvcf="data/output/called/allsites.vcf"
	output:
		"data/output/called/allsites.snp.vcf"
	threads:
		config["n_cores"]
	shell: """
		gatk SelectVariants \
		-V {input.allvcf} \
		-O {output} \
		-select-type SNP
	"""
	
	
rule extract_INDELS:
	input:
		allvcf="data/output/called/allsites.vcf"
	output:
		"data/output/called/allsites.indel.vcf"
	threads:
		config["n_cores"]
	shell: """
		gatk SelectVariants \
		-V {input.allvcf} \
		-O {output} \
		-select-type INDEL
	"""

#4.a) filter SNPs

rule SNP_filtration:
	input:
		ref = REF_GENOM,
		allvcf = "data/output/called/allsites.snp.vcf"
	output:
		protected("data/output/called/filtered/allsites.filtered.snps.vcf")
	threads:
		config["n_cores"]
	shell: """
		gatk VariantFiltration \
		-R {input.ref} \
		-V {input.allvcf} \
		-O {output} \
		--filter-expression "QD < 2.0" \
		--filter-name "QD2" \
		--filter-expression "FS > 60.0" \
		--filter-name "FS60" \
		--filter-expression "MQ < 40.0" \
		--filter-name "MQ40" \
		--filter-expression "SOR > 3.0" \
		--filter-name "SOR3" \
		--filter-expression "MQRankSum < -12.5" \
		--filter-name "MQRS12.5" \
		--filter-expression "ReadPosRankSum < -8.0" \
		--filter-name "RPRS8" \
		--filter-expression "QUAL < 30.0" \
		--filter-name "QUAL30"
		"""


#indels

rule indel_filtration:
	input:
		ref = REF_GENOM,
		allvcf = "data/output/called/allsites.indel.vcf"
	output:
		protected("data/output/called/filtered/allsites.filtered.indels.vcf")
	threads:
		config["n_cores"]
	shell: """
		gatk VariantFiltration \
		-R {input.ref} \
		-V {input.allvcf} \
		-O {output} \
		--filter-expression "QD < 2.0" \
		--filter-name "QD2" \
		--filter-expression "FS > 200.0" \
		--filter-name "FS200" \
		--filter-expression "ReadPosRankSum < -20.0" \
		--filter-name "RPRS20" \
		--filter-expression "QUAL < 30.0" \
		--filter-name "QUAL30"
		"""



