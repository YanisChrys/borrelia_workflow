# Use our own samples as bootstraps to recalibrate the quality scores (for model organisms ready-made dbSNP files are used instead).
# Sequencers are usually overconfident about their calls, or different sequencers have systematic biases
# in the way they assign Q scores. With samples from different runs, this step is important to remove biases and to 
# create more reliable quality scores that are consistent throughout all the samples.

# use reference genome index file to create a list of intervals to use for GATK DBImport: all "chromosomes":
rule get_chromosome_names:
	input:
		"data/ref_genom/GCF_000172275.2_ASM17227v2_genomic.fasta.fai"
	output:
		"data/ref_genom/interval.Bb.list"
	shell:
		"awk ' {{print $1}}' {input} > {output}"


######################### Begin 3 rounds of Base recalibration:
## the recalibration is complete if the quality scores converge.
## we check for Q score convergence with the rule "level_1_recalibration_diagnostic"
## which assembles all 3 rounds and graphically summarises the results

# call variants using HaplotypeCaller in GVCF mode:
# use a ploidy of 1 and don't allow the caller to note spanning deletions, limiting the number of "*" bases in the output
rule call_variants1:
	input:
		sample = "data/output/final_alignment/{sample}.bam",
		index = "data/output/final_alignment/{sample}.bam.bai",
		ref = REF_GENOM
	output:
		temp("data/output/recalibration1/myvcf/{sample}.g.vcf.gz")
	threads:
		n_cores
	shell:
		"gatk HaplotypeCaller --sample-ploidy 1 --disable-spanning-event-genotyping -ERC GVCF -R {input.ref} -I {input.sample} -O {output}"
		
# create database used for combining the GVCFS
# genomicsdb-shared-posixfs-optimizations and --batch-size 2 are used to increase the speed of the operation
rule genomicsdbimport1:
	input:
		ref=REF_GENOM,
		mylist="data/ref_genom/interval.Bb.list",
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
		--batch-size 2 \
		-L {input.mylist}
	"""

# joint call
# call joint variants on all samples
# use a temporary directory to reduce RAM use
# set minimum quality score for calling a variant to 45 
rule joint_variant_calling1:
	input:
		ref = REF_GENOM,
		my_db="data/output/recalibration1/my_dbi_database"
	output:
		temp("data/output/recalibration1/allsites.vcf.gz")
	threads:
		n_cores
	shell: """
		mkdir -p data/output/recalibration1/genotypegvcfstempdir && chmod a+rw data/output/recalibration1/genotypegvcfstempdir &&
		gatk GenotypeGVCFs \
		-R {input.ref} \
		-V gendb://{input.my_db}\
		-O {output} \
		--tmp-dir data/output/recalibration1/genotypegvcfstempdir/ \
		-stand-call-conf 45
	"""

# recalibrate bases of bam files by using the called variants above as a teaching sample
rule recalibrate1:
	input:
		ref=REF_GENOM,
		reads="data/output/final_alignment/{sample}.bam",
		dbsnp="data/output/recalibration1/allsites.vcf.gz"
	output:
		temp("data/output/recalibration1/{sample}.recalc.table")
	threads:
		n_cores
	shell: """
		gatk BaseRecalibrator \
		-I {input.reads} \
		-R {input.ref} \
		--known-sites {input.dbsnp} \
		-max-cycle 700 \
		-O {output}
	"""

# apply recalibration to the original bam files
rule apply_bqsr1:
    input:
        ref=REF_GENOM,
        reads="data/output/final_alignment/{sample}.bam",
        recalc_tb="data/output/recalibration1/{sample}.recalc.table"
    output:
        bam="data/output/recalibration1/mybam/{sample}.bam"
    threads:
        n_cores
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
		n_cores
	shell:
		"gatk HaplotypeCaller --sample-ploidy 1 --disable-spanning-event-genotyping -ERC GVCF -R {input.ref} -I {input.sample} -O {output}"	
        
# create db
rule genomicsdbimport2:
	input:
		ref=REF_GENOM,
		mylist="data/ref_genom/interval.Bb.list",
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
		--batch-size 2 \
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
		n_cores
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
		n_cores
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
		n_cores
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
		n_cores
	shell:
		"gatk HaplotypeCaller --sample-ploidy 1 --disable-spanning-event-genotyping -ERC GVCF -R {input.ref} -I {input.sample} -O {output}"	
        
# create db
rule genomicsdbimport3:
	input:
		ref=REF_GENOM,
		mylist="data/ref_genom/interval.Bb.list",
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
		--batch-size 2 \
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
		n_cores
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
		n_cores
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
		n_cores
	shell: """
		gatk ApplyBQSR \
		-R {input.ref}  \
		-I {input.reads}  \
		--bqsr-recal-file {input.recalc_tb} \
		--create-output-bam-index \
		-O {output}
	"""
	
	
# compare recalibration tables:
rule level_1_recalibration_diagnostic:
	input:
		table1="data/output/recalibration1/{sample}.recalc.table",
		table2="data/output/recalibration2/{sample}.recalc.table",
		table3="data/output/recalibration3/{sample}.recalc.table"
	output:
		"triphasic_covariate_recalibration_diagnostics/{sample}.pdf"
	shell: """
		gatk AnalyzeCovariates \
		-bqsr {input.table1} \
		-before {input.table2}  \
		-after {input.table3}  \
		-plots {output}
	"""

