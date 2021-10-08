# Borrelia afzelii analysis workflow

Can be installed in a linux system with miniconda3 (https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) and snakemake (https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) by doing:

```
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda

conda env create -f envs/vc_borrelia_env.yml -n borrelia
```

And then activate the environment containing the appropriate programs with:

```
conda activate borrelia
```

run whole analysis with:

```
snakemake -s MY_SNAKEFILE --use-conda --cores 4 all
```

Where `MY_SNAKEFILE` can be either `snakefile_vc` for the concensus genomes or `snakefile_mlst` for the mlst analysis.

create a schematic representation of all jobs with:

```
snakemake -s MY_SNAKEFILE --forceall --dag | dot -Tpdf > graph_of_jobs.pdf
```


 - For the consensus genomes:

snakefile_vc: Pipeline for reference mapping and variant calling of _Borrelia_ samples. Two types of fasta consensus sequences are created with `bcftools consensus` (interesting recent alternative: `VCFCons`) and one plink bedfile for the stricter filtering.

The reference genome needs to be unzipped, indexed and placed inside `data/ref_genom/` and the input bam files are placed under `data/input/multi_refs`. All other files and folders are created by snakemake. NOTE: reference file name is searched exactly, if running verify that all file names are the same or change the names in the code.

The mapping uses `bwa mem` and the variant calling uses `gatk HaplotypeCaller` in gVCF mode.

Currently the code uses `GenomicsDBImport` to combine variants before jointly calling them which requires at least 100GB of free diskspace to run.


 - For the MLST analysis, all alleles for each housekeeping gene of interest should be downloaded and placed under `data/housekeeping_genes/{gene}_db/{gene}.fasta` (where "{gene}" the name of the gene).

The consensus genomes that will be used for the MLST analysis should be placed under `data/borrelia_samples/{sample}.fasta` (where "sample the name of the sample"), including the reference.

 - Installing other packages:

For other analyses install separate environments or install them in linux:

For phyml install with apt-get:
```
sudo apt-get update -y
sudo apt-get install -y phyml
```

trimAL (http://trimal.cgenomics.org/), mafft (https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html), beast (https://beast.community/), snp-sites (https://github.com/sanger-pathogens/snp-sites) and the MLST analysis run better in conda:
```
conda create --name trimal -c bioconda trimal

conda create --name mafft -c bioconda mafft

conda create --name snp-sites -c bioconda snp-sites

conda env create -f envs/mlst_env.yml -n mlst_env

conda create --name beast -c bioconda beast

conda create -n beast2 -c conda-forge -c bioconda beast2
```

Whenever a program needs to run, the appropriate environment should simply be activated.

Additional scripts:

`01_alignment_commands.txt: find samples with high genotyping rate and perform mafft alignment`

`02_alignment_partitions.txt : uses trimAL to remove sites with missing data`

`03_modeltest_and_phyml.sh: use jmodeltest to find substitution model and phyml to create ml tree`
