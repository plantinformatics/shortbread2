#!/bin/bash
module load Nextflow
nextflow run /group/grains/git/shortbread2/main.nf -profile 'slurm' -resume -bg \
  --mode 'prod' \
  --aligner 'bwamem2' \
  --refgenome "/group/genomes/plant_species/Eucalyptus/Eglobulus-X46/EGLOB-X46.v1.0-CpMt-no-ChrUn.fa" \ #change this
  --indexdir "/group/genomes/plant_species/Eucalyptus/Eglobulus-X46/index/" \ #Change this
  --indexprefix "build" \
  --outdir "/group/trees/TBA-WGSanalysis/EGLOBULUS/240320-shorbread-analysis-of-referenceset/" \ #Change this
  -params-file /group/grains/git/shortbread2/params.yml \ #Leave this one as it is
  --keepbams true \ #Leave it
  --markduplicates true \ #This skips deduplication - set to false if you want deduplication to run
  --gatk true \ #This means only the alignment step will run
  --samplesheet "/group/trees/TBA-WGSanalysis/EGLOBULUS/240320-shorbread-analysis-of-referenceset/Samplesheet-reference-samples.csv" #Should be the full path


  #To run this script

  #nohup sh Run-shortbread.sh &