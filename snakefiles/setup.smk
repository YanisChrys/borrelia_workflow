# Setup pipeline constants by parsing filenames of raw read data
#   Requires: a glob object named RAW_READS, containing fields
#               run_id, filestem, samplenum, lane, read
#   This is part of a larger pipeline, more fully described in the main makefile (../Snakemake)


# Parse file names and get unique sample IDs:
sampleRegex = config["sample_name_regex"]
SAMPLE_FILES = {}
RUN_ID_DICT = {}
SAMPLE_NAMES = []  # Linear versions of all sampleID / runID combo's (suitable for zip())
RUN_IDS = []

for i in range(len(RAW_READS.run_id)):
    runID = RAW_READS.run_id[i]
    fileStem = RAW_READS.filestem[i]
    sampleID = search(sampleRegex, fileStem).group(1)
    samplenum = RAW_READS.samplenum[i]
    lane = RAW_READS.lane[i]
    read = RAW_READS.read[i]

    # Ensure all sample names use the same delimiter (so resequences can be detected):
    sampleID = sampleID.replace('_', '-')

    # Initialise sample file name dictionary (if needed):
    if runID not in SAMPLE_FILES.keys():
        SAMPLE_FILES[runID] = {}

    if sampleID not in SAMPLE_FILES[runID].keys():
        SAMPLE_FILES[runID][sampleID] = []

    # Store all file name for this sample / runID combo
    fullstem = "{runid}/{stem}_{num}_{lane}".format(runid = runID,
                                                    stem = fileStem,
                                                    num = samplenum,
                                                    lane = lane)

    if fullstem not in SAMPLE_FILES[runID][sampleID]:  # Each stem could occur twice (for paired end reads)
        SAMPLE_FILES[runID][sampleID].append(fullstem)

    # Store all unique sample / runID combo's
    if sampleID not in RUN_ID_DICT.keys():
        RUN_ID_DICT[sampleID] = [runID]
        SAMPLE_NAMES.append(sampleID)
        RUN_IDS.append(runID)
    elif runID not in RUN_ID_DICT[sampleID]:
        RUN_ID_DICT[sampleID].append(runID)
        SAMPLE_NAMES.append(sampleID)
        RUN_IDS.append(runID)


# Expected output given these input files:
FINAL_ALIGNMENTS = expand("output/final_alignment/{run_id}/{sample}.bam", zip, run_id = RUN_IDS, sample = SAMPLE_NAMES)
