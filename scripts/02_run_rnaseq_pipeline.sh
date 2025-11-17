#!/bin/bash
# RNA-seq Processing Pipeline using nf-core/rnaseq
# This script runs the complete nf-core/rnaseq pipeline for preprocessing, QC, and quantification

# Nextflow RNA-seq pipeline command
nextflow run nf-core/rnaseq \
  -r 3.18.0 \
  -name fares_300 \
  -profile docker \
  -work-dir /home/sinmo/fares_masters/outdir/work \
  -resume \
  -params-file nf-params.json

# After completion, check the MultiQC report for QC summary
echo "Pipeline completed. Check MultiQC report in results/multiqc/"
