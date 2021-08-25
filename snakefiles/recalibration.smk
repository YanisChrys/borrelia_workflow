#0) merge alignments of the same sample across multiple runs
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

# - Re-index merged bam files
rule index_merged:
	input:
		"data/output/run_ids_merged/{sample}.bam"
	output:
		protected("data/output/run_ids_merged/{sample}.bam.bai")
	shell:
		"picard BuildBamIndex -I {input} -O {output}"

# create list of intervals to use for the DB: all "chromosomes":
rule get_chromosome_names:
	input:
		"data/ref_genom/GCF_000222835.1_ASM22283v1_genomic.fasta.fai"
	output:
		"data/ref_genom/interval.list"
	shell:
		"awk ' {{print $1}}' {input} > {output}"



#Base recalibration:
# use our own data to perform a sort of bootstrap recalibration
# 3 times and correct for the different batch and sequencing technilogy biases in the quality/probability scores

# call variants
rule call_variants1:
	input:
		sample = "data/output/run_ids_merged/{sample}.bam",
		index = "data/output/run_ids_merged/{sample}.bam.bai",
		ref = REF_GENOM
	output:
		temp("data/output/recalibration1/myvcf/{sample}.g.vcf.gz")
	threads:
		config["n_cores"]
	shell:
		"gatk HaplotypeCaller --sample-ploidy 1 -ERC GVCF -R {input.ref} -I {input.sample} -O {output}"
		
# create db
rule genomicsdbimport1:
	input:
		ref=REF_GENOM,
		mylist="data/ref_genom/interval.list",
		gvcfs = expand("data/output/recalibration1/myvcf/{sample}.g.vcf.gz", sample=FINAL_SAMPLES)
	output:
		directory("data/output/recalibration1/my_dbi_database")
	shell: """
		mkdir -p data/output/recalibration1/dbimporttempdir/ && chmod a+rw data/output/recalibration1/dbimporttempdir &&
		gatk --java-options "-Xmx1g -Xms1g" GenomicsDBImport \
		$(echo {input.gvcfs} | sed 's/data/ -V data/g') \
		-R {input.ref} \
		--genomicsdb-shared-posixfs-optimizations \
		--genomicsdb-workspace-path {output} \
		--tmp-dir data/output/recalibration1/dbimporttempdir  \
		--batch-size 3 \
		-L {input.mylist}
	"""

# joint call
# call joint variants on all samples	
rule joint_variant_calling1:
	input:
		ref = REF_GENOM,
		my_db="data/output/recalibration1/my_dbi_database"
	output:
		temp("data/output/recalibration1/allsites.vcf.gz")
	threads:
		config["n_cores"]
	shell: """
		mkdir -p data/output/recalibration1/genotypegvcfstempdir && chmod a+rw data/output/recalibration1/genotypegvcfstempdir &&
		gatk GenotypeGVCFs \
		-R {input.ref} \
		-V gendb://{input.my_db}\
		-O {output} \
		--tmp-dir data/output/recalibration1/genotypegvcfstempdir/ \
		-stand-call-conf 45
	"""

# recalibrate bases of bam files
rule recalibrate1:
	input:
		ref=REF_GENOM,
		reads="data/output/run_ids_merged/{sample}.bam",
		dbsnp="data/output/recalibration1/allsites.vcf.gz"
	output:
		temp("data/output/recalibration1/{sample}.recalc.table")
	threads:
		config["n_cores"]
	shell: """
		gatk BaseRecalibrator \
		-I {input.reads} \
		-R {input.ref} \
		--known-sites {input.dbsnp} \
		-max-cycle 700 \
		-O {output}
	"""

# apply recalibration

rule apply_bqsr1:
    input:
        ref=REF_GENOM,
        reads="data/output/run_ids_merged/{sample}.bam",
        recalc_tb="data/output/recalibration1/{sample}.recalc.table"
    output:
        bam="data/output/recalibration1/mybam/{sample}.bam"
    threads:
        config["n_cores"]
    shell: """
        gatk ApplyBQSR \
        -R {input.ref}  \
        -I {input.reads}  \
		--bqsr-recal-file {input.recalc_tb} \
		--create-output-bam-index \
		-O {output}
	"""


	
##### 2


# call variants 2
rule call_variants2:
	input:
		sample = "data/output/recalibration1/mybam/{sample}.bam",
		ref = REF_GENOM
	output:
		temp("data/output/recalibration2/myvcf/{sample}.g.vcf.gz")
	threads:
		config["n_cores"]
	shell:
		"gatk HaplotypeCaller --sample-ploidy 1 -ERC GVCF -R {input.ref} -I {input.sample} -O {output}"	
        
# create db
rule genomicsdbimport2:
	input:
		ref=REF_GENOM,
		mylist="data/ref_genom/interval.list",
		gvcfs = expand("data/output/recalibration2/myvcf/{sample}.g.vcf.gz", sample=FINAL_SAMPLES)
	output:
		directory("data/output/recalibration2/my_dbi_database")
	shell: """
		mkdir -p data/output/recalibration2/dbimporttempdir/ && chmod a+rw data/output/recalibration2/dbimporttempdir &&
		
		gatk --java-options "-Xmx1g -Xms1g" GenomicsDBImport \
		$(echo {input.gvcfs} | sed 's/data/ -V data/g') \
		-R {input.ref} \
		--genomicsdb-shared-posixfs-optimizations \
		--genomicsdb-workspace-path {output} \
		--tmp-dir data/output/recalibration2/dbimporttempdir  \
		--batch-size 3 \
		-L {input.mylist}
	"""

# joint call
rule joint_variant_calling2:
	input:
		ref = REF_GENOM,
		my_db="data/output/recalibration2/my_dbi_database"
	output:
		temp("data/output/recalibration2/allsites.vcf.gz")
	threads:
		config["n_cores"]
	shell: """
		mkdir -p data/output/recalibration2/genotypegvcfstempdir && chmod a+rw data/output/recalibration2/genotypegvcfstempdir &&
		gatk GenotypeGVCFs \
		-R {input.ref} \
		-V gendb://{input.my_db} \
		-O {output} \
		--tmp-dir data/output/recalibration2/genotypegvcfstempdir/ \
		-stand-call-conf 45
	"""
    
# recalibrate bases of bam files
rule recalibrate2:
	input:
		ref=REF_GENOM,
		reads="data/output/recalibration1/mybam/{sample}.bam",
		dbsnp="data/output/recalibration2/allsites.vcf.gz"
	output:
		temp("data/output/recalibration2/{sample}.recalc.table")
	threads:
		config["n_cores"]
	shell: """
		gatk BaseRecalibrator \
		-I {input.reads} \
		-R {input.ref} \
		--known-sites {input.dbsnp} \
		-max-cycle 700 \
		-O {output}
	"""
    
# apply recalibration
rule apply_bqsr2:
	input:
		ref=REF_GENOM,
		reads="data/output/recalibration1/mybam/{sample}.bam",
		recalc_tb="data/output/recalibration2/{sample}.recalc.table"
	output:
		"data/output/recalibration2/mybam/{sample}.bam"
	threads:
		config["n_cores"]
	shell: """
		gatk ApplyBQSR \
		-R {input.ref}  \
		-I {input.reads}  \
		--bqsr-recal-file {input.recalc_tb} \
		--create-output-bam-index \
		-O {output}
	"""


##### 3


# call variants
rule call_variants3:
	input:
		sample = "data/output/recalibration2/mybam/{sample}.bam",
		ref = REF_GENOM
	output:
		temp("data/output/recalibration3/myvcf/{sample}.g.vcf.gz")
	threads:
		config["n_cores"]
	shell:
		"gatk HaplotypeCaller --sample-ploidy 1 -ERC GVCF -R {input.ref} -I {input.sample} -O {output}"	
        
# create db
rule genomicsdbimport3:
	input:
		ref=REF_GENOM,
		mylist="data/ref_genom/interval.list",
		gvcfs = expand("data/output/recalibration3/myvcf/{sample}.g.vcf.gz", sample=FINAL_SAMPLES)
	output:
		directory("data/output/recalibration3/my_dbi_database")
	shell: """
		mkdir -p data/output/recalibration3/dbimporttempdir && chmod a+rw data/output/recalibration3/dbimporttempdir &&
		gatk --java-options "-Xmx1g -Xms1g" GenomicsDBImport \
		$(echo {input.gvcfs} | sed 's/data/ -V data/g') \
		-R {input.ref} \
		--genomicsdb-shared-posixfs-optimizations \
		--genomicsdb-workspace-path {output} \
		--tmp-dir data/output/recalibration3/dbimporttempdir  \
		--batch-size 3 \
		-L {input.mylist}
    """
        
# joint call
rule joint_variant_calling3:
	input:
		ref = REF_GENOM,
		my_db="data/output/recalibration3/my_dbi_database"
	output:
		temp("data/output/recalibration3/allsites.vcf.gz")
	threads:
		config["n_cores"]
	shell: """
		mkdir -p data/output/recalibration3/genotypegvcfstempdir && chmod a+rw data/output/recalibration3/genotypegvcfstempdir &&
		gatk GenotypeGVCFs \
		-R {input.ref} \
		-V gendb://{input.my_db} \
		-O {output} \
		--tmp-dir data/output/recalibration3/genotypegvcfstempdir/ \
		-stand-call-conf 45
	"""
    
# recalibrate bases of bam files
rule recalibrate3:
	input:
		ref=REF_GENOM,
		reads="data/output/recalibration2/mybam/{sample}.bam",
		dbsnp="data/output/recalibration3/allsites.vcf.gz"
	output:
		temp("data/output/recalibration3/{sample}.recalc.table")
	threads:
		config["n_cores"]
	shell: """
		gatk BaseRecalibrator \
		-I {input.reads} \
		-R {input.ref} \
		--known-sites {input.dbsnp} \
		-max-cycle 700 \
		-O {output}
	"""

# apply final recalibration
rule apply_bqsr3:
	input:
		ref=REF_GENOM,
		reads="data/output/recalibration2/mybam/{sample}.bam",
		recalc_tb="data/output/recalibration3/{sample}.recalc.table"
	output:
		protected("data/output/RECAL_BAM/{sample}.bam")
	threads:
		config["n_cores"]
	shell: """
		gatk ApplyBQSR \
		-R {input.ref}  \
		-I {input.reads}  \
		--bqsr-recal-file {input.recalc_tb} \
		--create-output-bam-index \
		-O {output}
	"""