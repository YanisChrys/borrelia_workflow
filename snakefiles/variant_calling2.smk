# call variants and filter

#workflow:
#haplotypecaller gvcf
#combinegvcfs
#genotypegvcfs
#VariantFiltration

# 0) merge alignments of the same sample across multiple runs
# find all run_ids of each sample
# - from here onwards, run_id gets dropped from names (but should still be available in readgroup tags)
def find_sample_alignments(sample, runid_lookup = RUN_ID_DICT):
    runids = runid_lookup[sample]
    filenames = []

    for rid in runids:
        filenames.append("data/output/final_alignment/{0}/{1}.bam".format(rid, sample))

    return filenames

# if a sample has more than 1 runs then merge them, otherwise it stays the same
rule merge_resequenced:
	input: 
		lambda wildcards: find_sample_alignments(wildcards.sample)
	output:
		"data/output/run_ids_merged/{sample}.bam"
	shell:
		"samtools cat {input} > {output}"

# - Index merged bam file
# re-index the files
rule samtools_index_merged:
	input:
		"data/output/run_ids_merged/{sample}.bam"
	output:
		"data/output/run_ids_merged/{sample}.bam.bai"
	shell:
		"samtools index {input}"


# 1) call variants on each sample
rule call_variants:
	input:
		sample = "data/output/run_ids_merged/{sample}.bam",
		ref = REF_GENOM
	output:
		protected("data/output/called/{sample}.g.vcf.gz")
	shell:
		"gatk HaplotypeCaller --sample-ploidy 1 -ERC GVCF -R {input.ref} -I {input.sample} -O {output}"

# 2) combine all vcf files into one
# needs all samples as input
rule combine_gvcfs:
	input:
		gvcfs = expand("data/output/called/{sample}.g.vcf.gz", sample=FINAL_SAMPLES),
		ref = REF_GENOM
	output:
		"data/output/called/all.g.vcf"
	shell:
		"gatk CombineGVCFs -R {input.ref} --variant {input.gvcfs} -O {output}"

# 3) call joint variants on all samples	
rule joint_variant_calling:
	input:
		"data/output/called/all.g.vcf"
	output:
		"data/output/called/all.vcf"
	shell:
		"gatk GenotypeGVCFs -V {input} -O {output}"

#4) filter variants

rule joint_variant_filtration:
	input:
		"data/output/called/all.vcf"
	output:
		"data/output/called/all.filtered.vcf"
	shell: """ 
		gatk VariantFiltration 
		--filter expression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)""QD < 2.0 || DP < 40 || FS > 60.0 || MQ < 40.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0\"
		-R genome.fa -V all.vcf -O all.filtered.vcf
	"""







