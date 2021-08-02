# call variants and filter

#workflow:
#haplotypecaller gvcf
#combinegvcfs
#genotypegvcfs
#VariantFiltration

# 1) call variants on each sample
rule call_variants:
    input:
        sample="data/output/final_alignment/{run_id}/{sample}.bam",
		ref_genome=""
    output:
        protected("data/output/called/{sample}.g.vcf.gz")
    shell:
        "gatk HaplotypeCaller --sample-ploidy 1 -ERC GVCF -R {input.ref_genome} -I {input.sample} -O {output} "

# 2) combine all vcf files into one
# <needs all samples as input>
rule combine_gvcfs:
	input:
		ref_genome="data/ref_genome/{ref_genom}",
        gvcfs="data/output/called/{sample}.g.vcf.gz",
	output:
		"data/output/called/all.g.vcf"
	shell:
		"gatk CombineGVCFs -R {input.ref_genom} --variant {input.gvcfs} -O {output}"

# 3) call joint variants on all samples	
rule joint_variant_calling:
	input:
		"data/output/called/all.g.vcf"
	output:
		"data/output/called/all.vcf"
	shell:
		"gatk GenotypeGVCFs -V {input} -O {output}"

#4) filter variants

rule joint_variant_calling:
	input:
		"data/output/called/all.vcf"
	output:
		"data/output/called/all.filtered.vcf"
	shell: """ 
		gatk VariantFiltration 
		--filter expression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)""QD < 2.0 || DP < 40 || FS > 60.0 || MQ < 40.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0\"
		-R genome.fa -V all.vcf -O all.filtered.vcf
	"""








