# Tick diet sequence assemblies: Overall workflow
import os
from re import search

## SETUP:
configfile: "config.yml"

# Find all samples to process:
RAW_READS = glob_wildcards("data/reads/{run_id}/{filestem}_{samplenum}_{lane}_{read}_001.fastq.gz")

## Parse names:
include: "makefiles/setup.smk"

## Sub-workflows (partially independent branches):
include: "snakefiles/mapping2.smk"
include: "snakefiles/variant_calling2.smk"


