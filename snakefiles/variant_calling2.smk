# call variants and filter

#workflow:
#haplotypecaller gvcf
#combinegvcfs
#genotypegvcfs
#VariantFiltration

# 1) call variants on each sample
rule call_variants:
	input:
		sample = "data/output/final_alignment/{run_id}/{sample}.bam",
		ref = REF_GENOM
	output:
		protected("data/output/called/{run_id}/{sample}.g.vcf.gz")
	shell:
		"gatk HaplotypeCaller --sample-ploidy 1 -ERC GVCF -R {input.ref} -I {input.sample} -O {output}"

# 2) combine all vcf files into one
# <needs all samples as input>
rule combine_gvcfs:
	input:
		gvcfs = "data/output/called/{run_id}/{sample}.g.vcf.gz", 
		ref = REF_GENOM
	output:
		"data/output/called/{run_id}/all.g.vcf"
	shell:
		"gatk CombineGVCFs -R {input.ref} --variant {input.gvcfs} -O {output}"

# 3) call joint variants on all samples	
rule joint_variant_calling:
	input:
		"data/output/called/{run_id}/all.g.vcf"
	output:
		"data/output/called/{run_id}/all.vcf"
	shell:
		"gatk GenotypeGVCFs -V {input} -O {output}"

#4) filter variants

rule joint_variant_filtration:
	input:
		"data/output/called/{run_id}/all.vcf"
	output:
		"data/output/called/{run_id}/all.filtered.vcf"
	shell: """ 
		gatk VariantFiltration 
		--filter expression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)""QD < 2.0 || DP < 40 || FS > 60.0 || MQ < 40.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0\"
		-R genome.fa -V all.vcf -O all.filtered.vcf
	"""







