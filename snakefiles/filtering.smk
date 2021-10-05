# This subprocess filters SNPs and INDELs separately and remerges them into a single file which is used for the consensus sequences

# 1)
rule extract_SNPS:
	input:
		allvcf="data/output/called/allsites.vcf.gz"
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

# 2)
rule extract_INDELS:
	input:
		allvcf="data/output/called/allsites.vcf.gz"
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

# 3.a) filter SNPs
# based on the GATK best practises recomendations for SNP
rule SNP_filtration:
	input:
		ref = REF_GENOM,
		allvcf = "data/output/called/allsites.snp.vcf"
	output:
		"data/output/called/filtered/allsites.filtered.snps.vcf"
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


# 3.b) filter INDELS
# based on the GATK best practises recomendations for INDELS
rule INDEL_filtration:
	input:
		ref = REF_GENOM,
		allvcf = "data/output/called/allsites.indel.vcf"
	output:
		"data/output/called/filtered/allsites.filtered.indels.vcf"
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

# 4) merge vcfs with same number of samples but different sites
rule merge_filtered_vcfs:
    input:
        filt_snp="data/output/called/filtered/allsites.filtered.snps.vcf",
        filt_indel="data/output/called/filtered/allsites.filtered.indels.vcf"
    output:
        "data/output/called/filtered/allsites.filtered.vcf"
    threads:
        config["n_cores"]
    shell: """
    picard MergeVcfs -I {input.filt_snp} -I {input.filt_indel} -O {output}
    """
	
