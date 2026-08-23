
# Introduction

🚧 🛠️ 🏗️ 👷

This repo provides the bioinformatics pipeline, implemented in Nextflow, for our forthcoming preprint. It converts BAM files from scisoseq into tables of allele counts in single cells, to be used for downstream lineage tracing. 

# Prerequisites

Requirements:
* Slurm
* Nextflow
* Conda

Setup the conda environment:
```
conda create -n sclantern-nf
conda activate sclantern-nf

mamba install \
      python r \
      samtools minimap2 \
      vcftools bcftools tabix \
      bioconductor-rsamtools r-foreach r-doparallel \
      biopython numpy pandas pysam \
      gatk4 trgt=2.0
```

# Running

Run nextflow:
```
nextflow /PATH/TO/sclantern-nf/main.nf \
         -profile slurm \
         -resume \
         --sample_sheet /PATH/TO/SAMPLE/SHEET.csv \
         --ref_fa /PATH/TO/GRCh38_no_alt_analysis_set.fasta \
         --repeats_bed /PATH/TO/human_GRCh38_no_alt_analysis_set.platinumTRs-v1.0.trgt.bed \
         --outdir /PATH/TO/OUT/DIR
```

The sample sheet should be a csv file, here is an example:
```
sample_name,path
sample1,/PATH/TO/SAMPLE1/scisoseq.bam
sample2,/PATH/TO/SAMPLE2/scisoseq.bam
```

# Additional links

* [TRGT repeats file](https://zenodo.org/records/13178746) (used for `--repeats_bed`)
