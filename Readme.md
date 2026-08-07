
# Prerequisites

Setup and activate conda environment. Two alternative ways.
First, from the environment.yaml file:
```
conda create -f environment.yaml
```
Alternatively, create the conda env manually:
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

First, start a session to run nextflow
```
salloc --qos=desktop --partition=interactive_cpu -c 2 --mem 16G
```

Load modules and conda env
```
ml Nextflow
ml Micromamba
conda activate sclantern-nf
```

Run nextflow:
```
nextflow ~/devel/scrnallt-nf/main.nf \
         -profile shpc \
         -resume \
         --sample_sheet /PATH/TO/SAMPLE/SHEET.csv \
         --ref_fa /gstore/data/ctgbioinfo/taol9/mas/singularity/reference/GRCh38_no_alt_analysis_set.fasta \
         --repeats_bed /gstore/data/ctgbioinfo/taol9/mas/repeats/human_GRCh38_no_alt_analysis_set.platinumTRs-v1.0.trgt.bed \
         --outdir /PATH/TO/OUT/DIR
```

The sample sheet should be a csv file, here is an example:
```
sample_name,path
cells,/gstore/data/omni/biostat/lineage_tracing/hct116_cells_scisoseq_test.bam
nuclei,/gstore/data/omni/biostat/lineage_tracing/hct116_nuclei_scisoseq_test.bam
```
