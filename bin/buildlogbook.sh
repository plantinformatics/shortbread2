#!/bin/bash
#Provides function to build logbook

# Get the path to the Nextflow pipeline
pipeline_path=$1

# Get the path to the RMD file
rmd_path=$2

# Copy the Nextflow pipeline to the RMD file
cp $pipeline_path $rmd_path

# Add a YAML header to the RMD file
echo "---" >> $rmd_path
echo "title: Nextflow Pipeline Documentation" >> $rmd_path
echo "author: Your Name" >> $rmd_path
echo "date: $(date)" >> $rmd_path
echo "---" >> $rmd_path

# Render the RMD file
rmarkdown::render($rmd_path)