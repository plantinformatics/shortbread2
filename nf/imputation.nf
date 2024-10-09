File imputelogpath=new File(params.runpath+"/Logs/GATK/") //Path to GATK log files
File mainimputepath=new File(params.runpath+"/Imputation/") //Path to output of GATK
if(!imputelogpath.exists())
    imputelogpath.mkdirs()
if(!mainimputepath.exists())
    mainimputepath.mkdirs()
    
process RUN_BEAGLE{
    module 'Beagle,BCFtools'
    input:
        path vcf
        val refgenome
        val otheroptions
        val knownvariants
    output:

    script:
    outputvcf="${mainimputepath}/${vcf.baseName}-BeagleImputed"
    """
    java -jar beagle.22Jul22.46e.jar \\
      ref = ${refgenome} \\
      vcf = ${vcf} \\
      out = ${outputvcf} \\
      known = ${knownvariants} \\
      nthreads = $task.cpus ${otheroptions}
    tabix -C Chr${i}-BeagleImputed.vcf.gz
    """
}

process PLINK_VCF {
    tag "$meta.id"
    label 'process_medium'
    input:
    tuple path(vcf)
    output:
    tuple path("*.bed"), emit: bed, optional: true
    tuple path("*.bim"), emit: bim, optional: true
    tuple path("*.fam"), emit: fam, optional: true
    script:
    def args =
    """
    plink \\
        --vcf ${vcf} \\
        $args \\
        --threads $task.cpus \\
        --out ${prefix}
    """
}

processs GENOMIC_DISTANCE{
    input:
        val SNPs_file
        val Genome_file
    script:
    """
    # Create a hash table to store the SNP positions.
    SNP_positions=()
    while IFS=, read -r chromosome position; do
      SNP_positions["$chromosome"]="$position"
    done < "$SNPs_file"

    # Calculate the genomic distance between each SNP.
    for chromosome in "${!SNP_positions[@]}"; do
      # Get the start and end positions of the chromosome.
      start_position=$(awk -F '\t' -v chromosome="$chromosome" '$1 == chromosome {print $2}' "$Genome_file")
      end_position=$(awk -F '\t' -v chromosome="$chromosome" '$1 == chromosome {print $3}' "$Genome_file")

      # Calculate the genomic distance for each SNP in the chromosome.
      for position in "${SNP_positions[$chromosome]}"
      do
        distance=$((end_position - position))
        echo "$chromosome $position $distance"
      done
    done
    """

}

