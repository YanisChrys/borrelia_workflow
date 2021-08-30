#9) use new reference with indel file to extract indels that are called as alts:
rule extract_samples_indels:
	input:
		ref = "data/output/called/newref/{sample}.snp.ref.fasta",
		dict="data/output/called/newref/{sample}.snp.ref.dict",
		fai="data/output/called/newref/{sample}.snp.ref.fasta.fai",
		allvcf="data/output/called/filtered/allsites.filtered.indels.vcf"
	output:
		"data/output/called/vcf_extracted_SNP_INDELS/{sample}.snp.indel.vcf"
	threads:
		config["n_cores"]
	log:
		"logs/called/filtered/extract_indels_{sample}.log"
	shell: """
		gatk SelectVariants \
		-R {input.ref} \
		-V {input.allvcf} \
		-O {output} \
		--sample-name "{wildcards.sample}" \
		--remove-unused-alternates \
		-select 'vc.getGenotype("{wildcards.sample}").isHomVar() && vc.getGenotype("{wildcards.sample}").isCalled()'
	"""


#create mask with uncalled sites, all filters and low depth sites
#anything else will be a good quality variant
# 	-select 'DP < 2 && vc.getGenotype("RB-NU5H-MAY19-167").isHomVar()'
rule create_snp_indel_mask_files:
	input:
		ref = "data/output/called/newref/{sample}.snp.ref.fasta",
		my_sample="data/output/called/{sample}.g.vcf.gz",   #used for the sample name wildcard
		allvcf="data/output/called/allsites.vcf"
	output:
		"data/output/called/indelmasks/{sample}.mask.vcf"
	threads:
		config["n_cores"]
	log:
		"logs/called/masks/extract_snps_{sample}.log"
	shell: """
		gatk SelectVariants \ 
		-R {input.ref} \
		-V {input.allvcf} \
		-O {output} \
		--remove-unused-alternates \
		--sample-name "{wildcards.sample}" \
		-select 'vc.getGenotype("{wildcards.sample}").isNoCall() || DP < 2.0 || (vc.getGenotype("{wildcards.sample}").getAD.0 < 2.0 && vc.getGenotype("{wildcards.sample}").isHomRef() ) || (vc.getGenotype("{wildcards.sample}").getAD.1 < 2.0 && vc.getGenotype("{wildcards.sample}").isHomVar() ) || ((QD < 2.0 && vc.isSNP()) || (FS > 60.0 && vc.isSNP()) || (MQ < 40.0 && vc.isSNP()) || (SOR > 3.0 && vc.isSNP()) || (MQRankSum < -12.5 && vc.isSNP()) || (ReadPosRankSum < -8.0 && vc.isSNP()) || (QUAL < 30.0 && vc.isSNP()) ) || ((QD < 2.0  && vc.isIndel() ) || (FS > 200.0  && vc.isIndel() ) || (ReadPosRankSum < -20.0  && vc.isIndel()) || (QUAL < 30.0  && vc.isIndel()) )'
	"""


# 10)
rule create_indel_fasta:
	input:
		ref = "data/output/called/newref/{sample}.snp.ref.fasta",
		my_sample="data/output/called/vcf_extracted_SNP_INDELS/{sample}.snp.indel.vcf",
		my_mask="data/output/called/indelmasks/{sample}.mask.vcf"
	output:
		"data/output/called/newref/{sample}.snp.indel.ref.fasta"
	threads:
		config["n_cores"]
	shell: """
		gatk FastaAlternateReferenceMaker \
		-R {input.ref} \
		-V {input.my_sample} \
		-O {output} \
		--line-width 80 \
		--snp-mask {input.my_mask} \
		--snp-mask-priority \
		--use-iupac-sample "{wildcards.sample}" &&		
		printf '%s\n' '%s/>[0-9]./>/g' 'x' | ex {output}  &&
		printf '%s\n' '%s/:/ /g' 'x' | ex {output}
	"""

# 11) move second reference
rule move_snpindel_ref:
	input:
		"data/output/called/newref/{sample}.snp.indel.ref.fasta"
	output:
		protected("data/output/called/snp_indel_ref/{sample}.snp.indel.ref.fasta")
	run:
		os.renames(input[0], output[0])
	
# 7) index 2nd fasta
rule index_2nd_fasta:
	input:
		my_fasta="data/output/called/snp_indel_ref/{sample}.snp.indel.ref.fasta"
	output:
		"data/output/called/snp_indel_ref/{sample}.snp.indel.ref.fasta.fai"
	shell: """
		samtools faidx {input.my_fasta} -o {output}
	"""

# 8) create dictionary for 2nd fasta
rule create_dict_for_2nd_fasta:
	input:
		my_fasta="data/output/called/snp_indel_ref/{sample}.snp.indel.ref.fasta"
	output:
		"data/output/called/snp_indel_ref/{sample}.snp.indel.ref.dict"
	shell: """
		picard CreateSequenceDictionary -R {input.my_fasta} -O {output}
	"""