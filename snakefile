# Tick diet sequence assemblies: Overall workflow
import csv
import warnings
#import pysam

## SETUP:
configfile: "config.yaml"

REF_GENOM = "data/ref_genom/GCF_000222835.1_ASM22283v1_genomic.fasta"

## Samples to process
# Find all sample files in data directory:
AVAILABLE_FILES = glob_wildcards("data/input/multi_refs/{run_id}/{sample}.bam")


# Check qPCR results
known_positives = []

with open("data/metadata/borrelia_metadata2.csv", mode="r", encoding="utf-8-sig") as handle:
	sample_file = csv.DictReader(handle)
	
	for line in sample_file:
		if line["pcr result"].lower() == "pos":
			known_positives.append(line["sample"])
			
#			# fetch run_id using the sample if it exists in the filesystem:
#			if line["sample"] in AVAILABLE_FILES.sample:
#				run_id_index=AVAILABLE_FILES.sample.index(line["sample"])
#				my_run_id=AVAILABLE_FILES.run_id[run_id_index]
#				
#				my_sample="data/input/multi_refs/" + my_run_id + "/" + line["sample"] + ".bam"
#				
#				#get number of mapped reads, remove newline and turn to integer
#				num_mapped=int(pysam.view("-c", "-F", "4", my_sample).strip())
#				# if there are too few mapped reads write the sample name to a file to mark it as a potential false positive
#				if num_mapped<=500:
#					with open('false_positive_borrelia.txt', 'a+') as my_false_positives:
#						my_false_positives.write(line["sample"])
			
			



# Check that all positives have data files available:
# - Ignoring missing files with a warning
for sample in known_positives:
    if sample not in AVAILABLE_FILES.sample:
        warnings.warn("Input file not found for sample {}. Sample skipped.".format(sample))

# Build a list of samples to process and where to find them
# - A single sample might be associated with multiple run_ids
FINAL_SAMPLES = []
RUN_ID_DICT = {}  # Format is {sample: [run_id, run_id]}

for i in range(len(AVAILABLE_FILES.run_id)):
	sample = AVAILABLE_FILES.sample[i]
	run_id = AVAILABLE_FILES.run_id[i]

	if sample in known_positives:
		if sample not in FINAL_SAMPLES:
			FINAL_SAMPLES.append(sample)
			RUN_ID_DICT[sample] = [run_id]
		else:
			RUN_ID_DICT[sample].append(run_id)



## Sub-workflows (partially independent branches):
include: "snakefiles/mapping2.smk"
include: "snakefiles/recalibration.smk"
include: "snakefiles/variant_calling2.smk"
include: "snakefiles/filtering.smk"
include: "snakefiles/create_fasta.smk"
include: "snakefiles/newfasta.smk"



# Terminal node:
rule all:
	input:
#		"data/output/called/allsites.vcf.gz"
#		"data/output/my_dbi_database"
#		"data/output/called/allsites.vcf.gz"
		expand("data/output/called/vcf_allsites/{sample}.allsites.vcf", sample=FINAL_SAMPLES),
		expand("data/output/called/fasta/{sample}.newref.dict", sample=FINAL_SAMPLES),
		expand("data/output/called/leneintfasta/{sample}.lenient.dict", sample=FINAL_SAMPLES),
		expand("data/output/called/leneintfasta/{sample}.lenient.fasta", sample=FINAL_SAMPLES)
		
	
	
