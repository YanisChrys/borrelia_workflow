# This subprocess creates a plink bedfile from a vcf file with a more strict filter on called variants:
# we don't want to trust calls made with less than 2 reads

# 1) extract the chromosomal sequence for each sample in a separate file
# verbosity; only print to screen if there is an error
rule extract_sample_chrom:
	input:
		ref=REF_GENOM,
		my_sample="data/output/called/{sample}.g.vcf.gz",   #used for the sample name wildcard
		allvcf="data/output/called/allsites.vcf.gz"
	output:
		"data/output/called/stringent_vcf/{sample}.chrom.vcf"
	shell: """
		gatk SelectVariants \
		-R {input.ref} \
		-V {input.allvcf} \
		-O {output} \
		--remove-unused-alternates \
		-L NC_017238.1 \
		-verbosity ERROR \
		--sample-name "{wildcards.sample}"
	"""

# 2) apply GATK Best Practises filters and depth filters and set filtered sites to uncalled
rule stringent_filtration:
	input:
		ref = REF_GENOM,
		vcf = "data/output/called/stringent_vcf/{sample}.chrom.vcf"
	output:
		"data/output/called/stringent_vcf/{sample}.chrom_stringent.vcf.gz"
	shell: """
		gatk VariantFiltration \
		-R {input.ref} \
		-V {input.vcf} \
		-O {output} \
		--set-filtered-genotype-to-no-call \
		-verbosity ERROR \
		--filter-expression "QD < 2.0 && QD > -1.0 " \
		--filter-name "QD2" \
		--filter-expression "vc.isSNP() &&  FS > 60.0 " \
		--filter-name "FS60snp" \
		--filter-expression "vc.isIndel() && FS > 200.0 " \
		--filter-name "FS200indel" \
		--filter-expression "vc.isSNP() && MQ < 40.0" \
		--filter-name "MQ40snp" \
		--filter-expression "vc.isSNP() && SOR > 3.0" \
		--filter-name "SOR3snp" \
		--filter-expression "vc.isSNP() && MQRankSum < -12.5 " \
		--filter-name "MQRS12.5snp" \
		--filter-expression "vc.isIndel() && ReadPosRankSum < -20.0 " \
		--filter-name "RPRS20indel" \
		--filter-expression "vc.isSNP() && ReadPosRankSum < -8.0 " \
		--filter-name "RPRS8" \
		--filter-expression "QUAL < 30.0 && QUAL > -1.0" \
		--filter-name "QUAL30" \
		--filter-expression "DP < 2.0 " \
		--filter-name "DP2" \
		--filter-expression "vc.isMixed()" \
		--filter-name "mixed" \
		--filter-expression 'vc.getGenotype("{wildcards.sample}").isHomRef() && vc.getGenotype("{wildcards.sample}").hasAD() && vc.getGenotype("{wildcards.sample}").getAD().0 < 2.0' \
		--filter-name "AD2rf" \
		--filter-expression 'vc.getGenotype("{wildcards.sample}").isHomVar() && vc.getGenotype("{wildcards.sample}").hasAD() && vc.getGenotype("{wildcards.sample}").getAD().1 < 2.0' \
		--filter-name "AD2snp" 
	"""

# 3) merge all samples into a multi-sample vcf file and index it
rule merge_stringent vcf:
	input:
		expand("data/output/called/stringent_vcf/{sample}.chrom_stringent.vcf.gz", sample=FINAL_SAMPLES)
	output:
		"data/output/called/stringent_vcf/allchrom_stringent.vcf.gz"
	shell: """
		bcftools merge {input} -o {output}
		tabix -p vcf {output}
	"""


# 4) extract only SNPs from file
# remove the three non-afzelii samples that pass genotype filters
rule extract_stringent_SNPS:
	input:
		"data/output/called/stringent_vcf/allchrom_stringent.vcf.gz"
	output:
		"data/output/called/stringent_vcf/all_chrom_SNP_stringent.vcf.gz"
	shell: """
		gatk SelectVariants \
		-V {input} \
		-O {output} \
		-select-type SNP \
		-xl-sn RB-K41 \
		-xl-sn RB-K42 \
		-xl-sn RB-N2
	"""

# 5) create bedile:
# tell it we only have 1 chromosome, we don't have x or y chromosomes and to allow for non-standard chromosomes
# or it will not recognise the name of the chromosome
# filter for minimum allele frequency to exclude vey rare variants
# keep only sites with 100% coverage
# and samples with less than 50% missing data
rule stringent_bedfile:
	input:
		"data/output/called/stringent_vcf/all_chrom_SNP_stringent.vcf.gz"
	output:
		"data/output/called/stringent_bed/all_chrom_SNP_stringent"
	shell: """
		plink --vcf {input} --make-bed --chr-set 1 no-xy --allow-extra-chr --geno 0 --mind 0.5 --maf 0.05 --out {output}
	"""
	