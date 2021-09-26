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
snakemake --use-conda --cores 4 all
```

create a schematic representation of all jobs with:

```
snakemake --forceall --dag | dot -Tpdf > graph_of_jobs.pdf
```

Pipeline for reference mapping and variant calling of _Borrelia_ samples. Two types of fasta consensus sequences are created with `bcftools consensus`. (interesting recent alternative: `VCFCons`)

The reference genome needs to be unzipped, indexed and placed inside "data/ref_genom/"

The mapping uses `bwa mem` and the variant calling uses `gatk HaplotypeCaller` in gVCF mode.

Currently the code uses `GenomicsDBImport` to combine variants before jointly calling them which requires at least 100GB of free diskspace to run.



For other analyses install separate environments or install them in linux:

For phyml install with apt-get:
```
sudo apt-get update -y
sudo apt-get install -y phyml
```

trimal, mafft, snp-sites and the MLST analysis run better in conda:
```
conda create --name phyml -c bioconda trimal

conda create --name phyml -c bioconda mafft

conda create --name phyml -c bioconda snp-sites

conda env create -f envs/mlst_env.yml -n mlst_env
```