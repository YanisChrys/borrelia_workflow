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
		protected("data/output/run_ids_merged/{sample}.bam")
	shell:
		"samtools cat {input} > {output}"

# - Re-index merged bam files
rule index_merged:
	input:
		"data/output/run_ids_merged/{sample}.bam"
	output:
		protected("data/output/run_ids_merged/{sample}.bam.bai")
	shell:
		"picard BuildBamIndex -I {input} -O {output}"


# 1) call variants on each sample
rule call_variants:
	input:
		sample = "data/output/run_ids_merged/{sample}.bam",
		index = "data/output/run_ids_merged/{sample}.bam.bai",
		ref = REF_GENOM
	output:
		"data/output/called/{sample}.g.vcf.gz"
	benchmark:
		"benchmarks/{sample}.HapCaller.benchmark.txt"
	threads:
		config["n_cores"]
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
	benchmark:
		"benchmarks/combinegvcfs.benchmark.txt"
	threads:
		config["n_cores"]
	shell:
		"gatk CombineGVCFs -R {input.ref} $(echo {input.gvcfs} | sed 's/data/ -V data/g') -O {output}"


#"gatk GenomicsDBImport  -R {input.ref} $(echo {input.gvcfs} | sed 's/data/ -V data/g') -O {output} --genomicsdb-workspace-path data/output/my_dbi_database"


# 3) call joint variants on all samples	
rule joint_variant_calling:
	input:
		ref = REF_GENOM,
		gvcfs = "data/output/called/all.g.vcf"
	output:
		"data/output/called/allsites.vcf"
	benchmark:
		"benchmarks/genotypegvcfs.benchmark.txt"
	threads:
		config["n_cores"]
	shell:
		"gatk GenotypeGVCFs -R {input.ref} -V {input.gvcfs} -O {output} -all-sites"


