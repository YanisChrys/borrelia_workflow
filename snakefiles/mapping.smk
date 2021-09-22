# This subprocess extracts paired and singleton fastq reads from bamfiles, maps them to a reference, processes the output and merges samples from differen runs (re-requencing)

# 1.a) extract paired reads:
#paired (0x1) and exclude:
#read unmapped (0x4)
#not primary alignment (0x100)
#read fails platform/vendor quality checks (0x200)
#read is PCR or optical duplicate (0x400)
#supplementary alignment (0x800)
rule extract_paired_reads_from_bam:
	input: 
		"data/input/multi_refs/{run_id}/{sample}.bam"
	output:
		"data/input/reads/{run_id}/{sample}_paired.fastq"
	log:
		"logs/{run_id}/extract_reads_from_bam_{sample}_paired.log"
	shell:
		"samtools fastq -n -f 1 -F 3844 {input} > {output}" 


# 1.b ) extract singletons by excluding:
#read paired (0x1)
#read unmapped (0x4)
#not primary alignment (0x100)
#read fails platform/vendor quality checks (0x200)
#read is PCR or optical duplicate (0x400)
#supplementary alignment (0x800)
rule extract_singleton_reads_from_bam:
	input: 
		"data/input/multi_refs/{run_id}/{sample}.bam"
	output:
		"data/input/reads/{run_id}/{sample}_singleton.fastq"
	log:
		"logs/{run_id}/extract_reads_from_bam_{sample}_singleton.log"
	shell:
		"samtools fastq -n -F 3845 {input} > {output}"



# 2) map paired and unpaired reads separately
rule bwa_map_paired_reads:
	input:
		ref = REF_GENOM,
		reads = "data/input/reads/{run_id}/{sample}_paired.fastq",
		rg_source = "data/input/multi_refs/{run_id}/{sample}.bam"
	output:
		"data/output/mapped_reads/{run_id}/{sample}_paired.bam"
	threads:
		config["n_cores"]
	params:
		rg_group="'@RG\\tSM:{sample}\\tID:{run_id}_{sample}\\tPL:Illumina'"
	log:
		"logs/bwa_map_paired/{run_id}/{sample}_paired.log"
	shell:
		"""
		(bwa mem -M -R {params.rg_group} -t {threads} {input.ref} {input.reads} | samtools view -Sb - > {output})
		&> {log}
		"""

# 3)
rule bwa_map_singleton_reads:
	input:
		ref = REF_GENOM,
		reads = "data/input/reads/{run_id}/{sample}_singleton.fastq",
		rg_source = "data/input/multi_refs/{run_id}/{sample}.bam"
	output:
		"data/output/mapped_reads/{run_id}/{sample}_singleton.bam"
	threads:
		config["n_cores"]
	log:
		"logs/bwa_map_paired/{run_id}/{sample}_singleton.log"
	params:
		rg_group="'@RG\\tSM:{sample}\\tID:{run_id}_{sample}\\tPL:Illumina'"
	shell:
		"""
		(bwa mem -M -R {params.rg_group} -t {threads} {input.ref} {input.reads} | samtools view -Sb - > {output})
		&> {log}
		"""

# 4) merge paired and unpaired mapped reads
rule merge_alignments:
	input:
		paired_files = "data/output/mapped_reads/{run_id}/{sample}_paired.bam",
		singleton_files = "data/output/mapped_reads/{run_id}/{sample}_singleton.bam"
	output:
		temp("data/output/merged_alignments/{run_id}/{sample}.bam")
	shell:
		"samtools cat {input.paired_files} {input.singleton_files} > {output}"


# 5) extract good mapped reads
#exclude:
#read unmapped (0x4)
#not primary alignment (0x100)
#read fails platform/vendor quality checks (0x200)
#read is PCR or optical duplicate (0x400)
#supplementary alignment (0x800)
rule extract_only_mapped_reads:
	input:
		"data/output/merged_alignments/{run_id}/{sample}.bam"
	output:
		temp("data/output/extracted_merged_alignments/{run_id}/{sample}.bam")
	shell:
		"samtools view -b -F 3844 {input} > {output}"

		
# Alignment post-processing steps

# 6) Sort alignments
rule sort:
	input:
		"data/output/extracted_merged_alignments/{run_id}/{sample}.bam"
	output:
		"data/output/sorted_alignment/{run_id}/{sample}.bam"
	shell:
		"picard SortSam -I {input} -O {output} -SORT_ORDER coordinate "

# 7) Index alignments
rule index:
	input:
		"data/output/sorted_alignment/{run_id}/{sample}.bam"
	output:
		"data/output/sorted_alignment/{run_id}/{sample}.bam.bai"
	shell:
		"picard BuildBamIndex -I {input} -O {output}"


# 8) Mark remaining duplicates
# this step ensures the GATK variant caller doesn't use PCR duplicates
rule mark_duplicates:
	input:
		bam="data/output/sorted_alignment/{run_id}/{sample}.bam",
		bai="data/output/sorted_alignment/{run_id}/{sample}.bam.bai"
	output:
		bam = "data/output/duplicates_marked/{run_id}/{sample}.bam",
		metrics = protected("data/output/duplicates_marked/{run_id}/{sample}.dup_metrics.txt")
	log:
		"logs/mark_duplicates/duplicates_{run_id}_{sample}.log"
	threads:
		config["n_cores"]
	shell:
		"picard MarkDuplicates -I {input} -O {output.bam} -M {output.metrics} --REMOVE_DUPLICATES T"
		"&> {log}"


# 9) Re-index final alignment and move to final output location:
rule finalise_alignment:
	input:
		"data/output/duplicates_marked/{run_id}/{sample}.bam",
	output:
		alignment = protected("data/output/final_alignment/{run_id}/{sample}.bam"),
		index = protected("data/output/final_alignment/{run_id}/{sample}.bam.bai")
	run:
		os.renames(input[0], output.alignment)
		shell("samtools index {output.alignment}")  # Re-index


# 10) merge alignments of the same sample across multiple runs
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
		"picard MergeSamFiles $(echo {input} | sed 's/data/ -I data/g') -O {output}"

# 11) Re-index merged bam files
rule index_merged:
	input:
		"data/output/run_ids_merged/{sample}.bam"
	output:
		protected("data/output/run_ids_merged/{sample}.bam.bai")
	shell:
		"picard BuildBamIndex -I {input} -O {output}"