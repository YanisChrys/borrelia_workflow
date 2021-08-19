# 1.a) extract paired (0x1) mapped (0x2) reads
rule extract_paired_reads_from_bam:
	input: 
		"data/input/multi_refs/{run_id}/{sample}.bam"
	output:
		"data/input/reads/{run_id}/{sample}_paired.fastq"
	log:
		"logs/{run_id}/extract_reads_from_bam_{sample}_paired.log"
	shell:
		"samtools fastq -n -f 3 {input} > {output}"


# 1.b ) extract singletons by excluding paired read (0x1) and unmapped reads (0x4)
rule extract_singleton_reads_from_bam:
	input: 
		"data/input/multi_refs/{run_id}/{sample}.bam"
	output:
		"data/input/reads/{run_id}/{sample}_singleton.fastq"
	log:
		"logs/{run_id}/extract_reads_from_bam_{sample}_singleton.log"
	shell:
		"samtools fastq -n -G 2053 {input} > {output}"



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

# 3) merge paired and unpaired mapped reads

rule merge_alignments:
	input:
		paired_files = "data/output/mapped_reads/{run_id}/{sample}_paired.bam",
		singleton_files = "data/output/mapped_reads/{run_id}/{sample}_singleton.bam"
	output:
		temp("data/output/merged_alignments/{run_id}/{sample}.bam")
	shell:
		"samtools cat {input.paired_files} {input.singleton_files} > {output}"


# 4) extract mapped reads
rule extract_mapped_reads:
	input: 
		"data/output/merged_alignments/{run_id}/{sample}.bam"
	output:
		temp("data/output/extracted_merged_alignments/{run_id}/{sample}.bam")
	log:
		"logs/{run_id}/extract_mapped_reads_from_bam_{sample}.log"
	shell:
		"samtools view -b -F 2052 {input} > {output}"
		
# Alignment post processing steps
# 5) Sort alignments
rule samtools_sort:
	input:
		"data/output/extracted_merged_alignments/{run_id}/{sample}.bam"
	output:
		"data/output/sorted_alignment/{run_id}/{sample}.bam"
	shell:
		"picard SortSam -I {input} -O {output} -SORT_ORDER coordinate "

# 6) Index alignments
rule samtools_index:
	input:
		"data/output/sorted_alignment/{run_id}/{sample}.bam"
	output:
		"data/output/sorted_alignment/{run_id}/{sample}.bam.bai"
	shell:
		"picard BuildBamIndex -I {input} -O {output}"


# 7) Mark remaining duplicates
rule mark_duplicates:
	input:
		"data/output/sorted_alignment/{run_id}/{sample}.bam"
	output:
		bam = "data/output/duplicates_marked/{run_id}/{sample}.bam",
		metrics = protected("data/output/duplicates_marked/{run_id}/{sample}.dup_metrics.txt")
	log:
		"logs/mark_duplicates/duplicates_{run_id}_{sample}.log"
	threads:
		config["n_cores"]
	shell:
		"picard MarkDuplicates -I {input} -O {output.bam} -M {output.metrics} "
		"&> {log}"


# 9) Re-index final alignment and move to final output location:
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
