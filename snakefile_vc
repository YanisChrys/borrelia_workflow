# Pipeline for mapping reads, calling variants and creating consensus genomes in fasta format.
# Also produce plink bed files of variants

# This script creates a variable with the samples to be processed, creates the "all" rule 
# and calls the sub-processes of the pipeline

import csv
import warnings

# Random seed:
random_seed: 23102020

# save reference genome as a variable
REF_GENOM = "data/ref_genom/GCF_000222835.1_ASM22283v1_genomic.fasta"

## Samples to process
# Find all sample files in data directory of the form "data/input/multi_refs/run_id/sample.bam"
# and save names as wildcards
AVAILABLE_FILES = glob_wildcards("data/input/multi_refs/{run_id}/{sample}.bam")


# Check qPCR results so as to only use borrelia positive samples
known_positives = []

# open metadate file with qPCR information and if the samples is positive take a note of it
with open("data/metadata/borrelia_metadata2.csv", mode="r", encoding="utf-8-sig") as handle:
	sample_file = csv.DictReader(handle)
	
	for line in sample_file:
		if line["pcr result"].lower() == "pos":
			known_positives.append(line["sample"])
			

# Check that all positives have data files available:
# Ignoring missing files with a warning
for sample in known_positives:
    if sample not in AVAILABLE_FILES.sample:
        warnings.warn("Input file not found for sample {}. Sample skipped.".format(sample))

# Build a list of samples to process and where to find them
# A single sample might be associated with multiple run_ids that can be later merged
FINAL_SAMPLES = []
RUN_ID_DICT = {} 

for i in range(len(AVAILABLE_FILES.run_id)):
	sample = AVAILABLE_FILES.sample[i]
	run_id = AVAILABLE_FILES.run_id[i]

	if sample in known_positives:
		if sample not in FINAL_SAMPLES:
			FINAL_SAMPLES.append(sample)
			RUN_ID_DICT[sample] = [run_id]
		else:
			RUN_ID_DICT[sample].append(run_id)



# Sub-modules, each focusing on a specific task:
include: "snakefiles/mapping.smk"
include: "snakefiles/recalibration.smk"
include: "snakefiles/variant_calling.smk"
include: "snakefiles/filtering.smk"
include: "snakefiles/create_fasta.smk"
include: "snakefiles/newfasta.smk"
include: "snakefiles/create_stringent_bed.smk"



# Terminal node:
rule all:
	input:
		expand("data/output/called/vcf_allsites/{sample}.allsites.vcf", sample=FINAL_SAMPLES),
		expand("data/output/called/stringent_fasta/edited/{sample}.fasta", sample=FINAL_SAMPLES),
		expand("data/output/called/lenientfasta/edited/{sample}.lenient.dict", sample=FINAL_SAMPLES),
		expand("triphasic_covariate_recalibration_diagnostics/{sample}.pdf", sample=FINAL_SAMPLES),
		"data/output/called/stringent_bed/all_chrom_SNP_stringent.bed"
		
	
	
