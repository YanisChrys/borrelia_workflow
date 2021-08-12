# borrelia_workflow

Can be install in a linux with miniconda3 installed by doing:

```
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda install pysam

conda env create -f envs/environ_2.yml -n borrelia
```

Pipeline for mapping and calling Borrelia samples.

The reference genome needs to be unzipped, indexed and placed inside "data/ref_genom/"

The mapping uses bwa mem and the variant calling uses HaplotypeCaller in gVCF mode.
