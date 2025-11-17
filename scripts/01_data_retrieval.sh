#!/bin/bash
# Automated data retrieval from NCBI SRA using nf-core/fetchngs
# BioProject: PRJNA1097853

# Run nf-core fetchngs pipeline
nextflow run nf-core/fetchngs \
  -r 1.13.0 \
  -profile docker \
  --input sra_ids.csv \
  --outdir data/raw/

echo "Data retrieval completed. Raw FASTQ files are in data/raw/"
