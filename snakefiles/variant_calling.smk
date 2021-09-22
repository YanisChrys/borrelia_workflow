# call variants 

# 1) call variants on each sample
rule call_variants:
	input:
		sample="data/output/RECAL_BAM/{sample}.bam",
		ref = REF_GENOM
	output:
		"data/output/called/{sample}.g.vcf.gz"
	threads:
		config["n_cores"]
	shell:
		"gatk HaplotypeCaller --sample-ploidy 1 --disable-spanning-event-genotyping -ERC GVCF -R {input.ref} -I {input.sample} -O {output}"


# 2) create db of samples
#additional options resposible for performance
rule genomicsdbimport:
	input:
		ref=REF_GENOM,
		mylist="data/ref_genom/interval.list",
		gvcfs = expand("data/output/called/{sample}.g.vcf.gz", sample=FINAL_SAMPLES)
	output:
		directory("data/output/my_dbi_database")
	shell: """
		mkdir -p data/output/dbimporttempdir/ && chmod a+rw data/output/dbimporttempdir
			
		gatk --java-options "-Xmx1g -Xms1g" GenomicsDBImport \
		$(echo {input.gvcfs} | sed 's/data/ -V data/g') \
		-R {input.ref} \
		--genomicsdb-shared-posixfs-optimizations \
		--genomicsdb-workspace-path {output} \
		--tmp-dir data/output/dbimporttempdir  \
		--batch-size 3 \
		-L {input.mylist}
	"""


# 3) call joint variants on all samples	
rule joint_variant_calling:
	input:
		ref = REF_GENOM,
		dir="data/output/my_dbi_database"
	output:
		"data/output/called/allsites.vcf.gz"
	threads:
		config["n_cores"]
	shell: """
		mkdir -p data/output/genotypegvcfstempdir && chmod a+rw data/output/genotypegvcfstempdir 
		
		gatk GenotypeGVCFs \
		-R {input.ref} \
		-V gendb://{input.dir} \
		-O {output} \
		--tmp-dir data/output/genotypegvcfstempdir/ \
		-all-sites
	"""
	