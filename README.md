# borrelia_workflow

Pipeline for mapping and calling Borrelia samples.

The find_borrelia_positive_samples.sh script takes a text file(false_positive_borrelia.txt) with all the qPCR positive files and looks for them in a folder system of the type: data/input/multi_refs/{run_id}/
(none of those terminal directories can be empty or it will cause errors) and moves them to the folder data/input/multi_refs_positive/{run_id}/
which will be the input of the analysis. The script also checks for false positive samples by finding how many mapped reads they have. 
If it's positive and it has more than 500 mapped reads it's marked as a false positive inside the text file <xxx>.

The reference genome needs to be unzipped, indexed and placed inside "data/ref_genome/"

The mapping uses bwa mem and the variant calling uses HaplotypeCaller in gVCF mode.
