
# 1) extract paired and singleton reads from bam files.
# Each read is uniquely mapped to a reference, so we extract them all before remapping them.
# -t: copy RG, BC and QT tags to the FASTQ header line
rule extract_paired_reads_from_bam:
	input: 
		"data/input/multi_refs_positive/{run_id}/{sample}.bam"
	output:
		"data/input/reads/{run_id}/{sample}_paired.fastq"
	log:
		"logs/extract_reads_from_bam_{sample}_paired.log"
	shell:
		"samtools fastq -ntf 13 {input} > {output}"


rule extract_singleton_reads_from_bam:
	input: 
		"data/input/multi_refs_positive/{run_id}/{sample}.bam"
	output:
		"data/input/reads/{run_id}/{sample}_singleton.fastq"
	log:
		"logs/extract_reads_from_bam_{sample}_singleton.log"
	shell:
		"samtools fastq -ntf 4 {input} > {output}"	


# 2) map paired and unpaired reads separately
rule bwa_map_paired_reads:
	input:
		ref = REF_GENOM,
		reads = "data/input/reads/{run_id}/{sample}_paired.fastq",
	output:
		"data/output/mapped_reads/{run_id}/{sample}_paired.bam"
	threads:
		config["n_cores"]
	log:
		"logs/bwa_map_paired/{sample}_paired.log"
	shell:
		"""
		RG_GROUP="$(samtools view -H {input.reads} | grep '^@RG' | sed 's/\t/\\t/g')"
		(bwa mem -M -R $RG_GROUP -t {threads} {input.ref} {input.reads} | samtools view -Sb - > {output})
		&> {log}
		"""


rule bwa_map_singleton_reads:
	input:
		ref = REF_GENOM,
		reads = "data/input/reads/{run_id}/{sample}_singleton.fastq",
	output:
		"data/output/mapped_reads/{run_id}/{sample}_singleton.bam"
	threads:
		config["n_cores"]
	log:
		"logs/bwa_map_paired/{sample}_singleton.log"
	shell:
		"""
		RG_GROUP="$(samtools view -H {input.reads} | grep '^@RG' | sed 's/\t/\\t/g')"
		(bwa mem -M -R $RG_GROUP -t {threads} {input.ref} {input.reads} | samtools view -Sb - > {output})
		&> {log}
		"""

# 3) merge paired and unpaired mapped reads

rule merge_alignments:
	input:
		paired_files = lambda wildcards: ["data/output/mapped_reads/{run_id}/{0}_paired.bam".format(stem) for stem in SAMPLE_FILES[wildcards.run_id][wildcards.sample]],
		singleton_files = lambda wildcards: ["data/output/mapped_reads/{run_id}/{0}_singletons.bam".format(stem) for stem in SAMPLE_FILES[wildcards.run_id][wildcards.sample]]
	output:
		temp("data/output/merged_alignments/{run_id}/{sample}.bam")
	shell:
		"samtools cat {input.paired_files} {input.singleton_files} > {output}"
		
		
# Alignment post processing steps
# 4) Sort alignments
rule samtools_sort:
	input:
		"data/output/merged_alignments/{run_id}/{sample}.bam"
	output:
		temp("data/output/sorted_alignment/{run_id}/{sample}.bam")
	shell:
		"picard SortSam INPUT={input} OUTPUT={output} SORT_ORDER=coordinate "

# Index alignments
rule samtools_index:
	input:
		"data/output/sorted_alignment/{run_id}/{sample}.bam"
	output:
		temp("data/output/sorted_alignment/{run_id}/{sample}.bam.bai")
	shell:
		"samtools index {input}"


# Mark remaining duplicates
rule mark_duplicates:
	input:
		"data/output/sorted_alignment/{run_id}/{sample}.bam"
	output:
		bam = temp("data/output/duplicates_marked/{run_id}/{sample}.bam"),
		metrics = protected("data/output/duplicates_marked/{run_id}/{sample}.dup_metrics.txt")
	log:
		"logs/mark_duplicates/duplicates_{run_id}_{sample}.log"
	threads:
		config["n_cores"]
	shell:
		"picard -XX:ParallelGCThreads={threads} -Xmx5g MarkDuplicates I={input} O={output.bam} M={output.metrics} "
		"&> {log}"


# Re-index final alignment and move to final output location:
#  - Also see definition of FINAL_ALIGNMENTS constant in setup.smk, which makes this a terminal node
rule finalise_aligment:
	input:
		"data/output/duplicates_marked/{run_id}/{sample}.bam",
	output:
		alignment = protected("data/output/final_alignment/{run_id}/{sample}.bam"),
		index = protected("data/output/final_alignment/{run_id}/{sample}.bam.bai")
	run:
		os.renames(input[0], output.alignment)
		shell("samtools index {output.alignment}")  # Re-index
