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
		protected("data/output/called/{sample}.g.vcf.gz")
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

# 3) call joint variants on all samples	
rule joint_variant_calling:
	input:
		ref = REF_GENOM,
		gvcfs = "data/output/called/all.g.vcf"
	output:
		"data/output/called/all.vcf"
	benchmark:
		"benchmarks/genotypegvcfs.benchmark.txt"
	threads:
		config["n_cores"]
	shell:
		"gatk GenotypeGVCFs -R {input.ref} -V {input.gvcfs} -O {output}"

#4) filter variants

rule joint_variant_filtration:
	input:
		ref = REF_GENOM,
		allvcf = "data/output/called/all.vcf"
	output:
		protected("data/output/called/all.filtered.vcf")
	threads:
		config["n_cores"]
	shell: """
		gatk VariantFiltration \
		-R {input.ref} \
		-V {input.allvcf} \
		-O {output} \
		--filter-name "my_filters" \
		--filter-expression "QD < 2.0 || DP < 40 || FS > 60.0 || MQ < 40.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0\"
		"""


#5 Extract each sample from the vcf
rule extract_samples:
	input:
		ref = REF_GENOM,
		my_sample="data/output/called/{sample}.g.vcf.gz",
		allvcf="data/output/called/all.filtered.vcf"
	output:
		"data/output/consensus/{sample}.vcf"
	threads:
		config["n_cores"]
	shell: """
		gatk SelectVariants \
		-R {input.ref} \
		-V {input.allvcf} \
		-O {output} \
		--sample-name {wildcards.sample}
	"""
		
# 6) create fasta consensus genome for each sample for each plasmid and the chromosome
rule create_fasta:
	input:
		ref = REF_GENOM,
		my_sample="data/output/consensus/{sample}.vcf"
	output:
		"data/output/consensus/fasta/{sample}.fasta"
	threads:
		config["n_cores"]
	shell: """
		gatk FastaAlternateReferenceMaker \
		-R {input.ref} \
		-V {input.my_sample} \
		-O {output}
	"""
	
# 7) index fasta
rule index_fasta:
	input:
		"data/output/consensus/fasta/{sample}.fasta"
	output:
		"data/output/consensus/fasta/{sample}.fasta.fai"
	shell:
		"samtools faidx {input} -o {output}"

	









