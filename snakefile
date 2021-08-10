# Tick diet sequence assemblies: Overall workflow
import csv
import warnings

## SETUP:
#configfile: "config.yml"

REF_GENOM = "data/ref_genome/ref.fa"

## Samples to process
# Find all sample files in data directory:
AVAILABLE_FILES = glob_wildcards("data/input/multi_refs/{run_id}/{sample}.bam")

# Check qPCR results
known_positives = []

with open("data/metadata/SampleMetadata.csv", mode="r", encoding="utf-8-sig") as handle:
    sample_file = csv.DictReader(handle)
    
    for line in sample_file:
        if line["pcr result"].lower() == "pos":
            known_positives.append(line["sample"])

# <add code to check false positives here>

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
include: "snakefiles/variant_calling2.smk"

# Terminal node:
rule all: 
    input: "data/output/called/all.filtered.vcf"
