def tabixmaxsize = 530000000 //Current maximum chrom size used by Htslib tabix function during indexation

File mainlogpath=new File(params.runpath+"/Logs/04_GATK/") //Path to GATK log files
File maingatkpath=new File(params.runpath+"/GATK/") //Path to output of GATK
File veppath=new File(params.runpath+"/VEP/") //Path to output of Variant Effect Predictor
if(!mainlogpath.exists())
    mainlogpath.mkdirs()
if(!maingatkpath.exists())
    maingatkpath.mkdirs()

def bamoutspath="${launchDir}/Alignments/GATK/"
if(params.keepbams)
{
    bamoutspath=params.runpath+"/Alignments/GATK/"
}

/*
    Run  HaplotypeCaller
*/
process RUN_GATK_HAPLOTYPE_CALLER{
    tag "${sampleid}_${intervalname}"
    label 'usescratch'
    label "${params.gatkhap}"
    publishDir "${launchDir}/GVCFs/${sampleid}/", mode:'copy'
    label 'GATK'
    input:
        tuple val(intervals), val(sampleid), path(alignments),path(alignments_index)
        val refgenome
        val otheroptions
        val maxchromsize
        val allele
    output:
        tuple val("${intervals}"),path("${gvcf}"),path("${gvcf}.tbi"),env(index)
    script:
        outdir_GVCF="${maingatkpath}/GVCFs/${sampleid}/"
        intervalname="${intervals.replaceAll(':|-|__','_')}"
        gvcf="${sampleid}_SB2_${intervalname}.output.g.vcf.gz"
        logpath="${mainlogpath}/a_Haplotyping/${sampleid}/"
        bamsout=""
        if(params.GATKbamout)
        {
            bamsout="-bamout ${intervalname}_GATK.bam"
        }

    """
    #!/bin/bash
    set -euxo pipefail
    #create a logging directory for haplotype calls
    mkdir -p ${logpath}

    #Check if alignments actually exist
    echo Calling haplotypes on ${intervals}
    if [[ ! -e ${alignments} ]];
    then
        echo "${alignments} was not found!!!"
        exit 2
    fi

    index='true'
    #Check if the maximum chromosome size exceeds the limit in htslib and create temp index and set index to false
    if [[ ${maxchromsize} -ge ${tabixmaxsize} ]];
    then
        index='false'
        echo "--" > "${gvcf}.tbi"
    fi

    intervalvcf="${intervals}"

    #Check if a reference vcf has been passed
    if [[ -e "${allele}" && "${params.GATKreferenceVARSonly}" == 'true' ]];
    then
        #forcecall=true
        intervalvcf="${intervalname}.vcf"
        bcftools view --threads ${task.cpus} -r "${intervals}" \\
        -Ov -o \${intervalvcf} "${allele}.gz"
    fi

    allelevcf=${allele}
    if [[ -z "\${allelevcf}" || ! -f "${allele}" ]]; then
        allelevcf=null  # Explicitly set to 'null' if the condition is met
    fi

    #Run haplotype caller on each interval
    exec > Haplotype.log 2>&1
    gatk --java-options "-Xmx30g -Xms30g -XX:ParallelGCThreads=${task.cpus}" HaplotypeCaller  \\
                   -R "${refgenome}" \\
                   -I "${alignments}" \\
                   -O "${gvcf}" \\
                   -L "\${intervalvcf}" \\
                   --sample-name "${sampleid}" \\
                   -OVI \${index} \\
                   --pair-hmm-implementation "AVX_LOGLESS_CACHING_OMP" \\
                   --native-pair-hmm-threads ${task.cpus} ${bamsout} \\
                   -ERC GVCF ${otheroptions} --alleles  \${allelevcf}

    if [[ "${params.GATKbamout}" == "true" ]];
    then
        mkdir -p "${bamoutspath}/${sampleid}/"
        rsync -rvP "${intervalname}_GATK.bam" "${bamoutspath}/${sampleid}/"
    fi
    rsync -rvP Haplotype.log ${logpath}/${intervalname}.log
    """
}

/*
    Run GATK genomics DB import
*/
process BUILD_GENOMICSDBImport{
    tag "${interval}"
    label 'dbimport'
    label 'usescratch'
    label "${params.gatkmulti}"
    label 'GATK'
    input:
        tuple val(intervals), path(gvcfs), path(gvcfs_index),val(index)
        val otherGATKOptions
    output:
        tuple val("${chrom}"),val("${intervals}"), val("${outdir}"),val("${index}")
    script:
    interval="${intervals}".replaceAll(':|-','_')
    outdir="${maingatkpath}/GATKDBs/${interval}"
    chrom="${intervals}".split(":")[0]
    logpath="${mainlogpath}/b_GENOMICSDBImport/"
    """
    #!/bin/bash
    set -euxo pipefail
    echo "Running DB import on ${intervals} belonging to chromosome ${chrom}!"
    mkdir -p \$(dirname ${outdir})
    mkdir -p ${logpath}
    mapFilePath="mapfile.tsv"
    #Build a sample map file to be passed to dbimport
    for s in ${gvcfs};
    do
    {
        f=\${s}
        fname=\$(echo \$s|sed -e "s/_SB2_.*//")
        if [[ "${index}" == "false" ]];
        then
            bgzip -f -d --keep --threads $task.cpus \${f}
            ff=\${f/".gz"/""}
            gatk IndexFeatureFile --input \${ff}
            rm -f \${f}.tbi
            printf "\${fname}\\t\${ff}\\n">>\${mapFilePath}
        else
            printf "\${fname}\\t\${f}\\n">>\${mapFilePath}
        fi
    }
    done
    exec > dbimport.log 2>&1
    gatk --java-options "-Xmx30g -Xms30g -XX:ParallelGCThreads=${task.cpus}" GenomicsDBImport   \\
        --genomicsdb-workspace-path ${outdir} \\
        -L "${intervals}" \\
        --sample-name-map \${mapFilePath} \\
        --consolidate true \\
        -OVI "${index}" \\
        --overwrite-existing-genomicsdb-workspace true \\
        --reader-threads $task.cpus ${otherGATKOptions}

    if [ \$? -gt 0 ]
    then
        echo "ERROR; Failed to run GATK GenomicsDBImport"
        exit 1
    fi
    rsync -rvP dbimport.log ${logpath}/${interval}.log
    """
}

/*
Process to read gatk databases if only these have been provided
*/
process READ_GATKDBs{
    tag "${file(file(gatkdb).getParent()).baseName}"
    executor 'local'
    input:
        each gatkdb
    output:
         tuple val("${intervals}"),val("${chrom}"),val("${outdir}"),val("${index}")
    script:
    outdir=new File("${gatkdb}").getParent()
    intervals=new File(outdir).getName().toString().replaceFirst("_",":").replace("_","-")
    chrom=chrom="${intervals}".split(":")[0]
    index=false
    """
    #!/bin/bash
    set -euxo pipefail
    echo $gatkdb line
    #cat .command.log >> "${mainlogpath}/READ_GATKDBs.log"
    """
}

/*
    UPDATE GATK genomics DB import
*/
process UPDATE_GENOMICSDB{
    tag "${interval}"
    label 'dbimport'
    label 'usescratch'
    label "${params.gatkmulti}"
    label 'GATK'
    input:
        tuple val(intervals),path(gvcfs), path(gvcfs_index),val(index)
        val otherGATKOptions
        val GATKupdateexistingdb
        val GATKpathtodbs
    output:
        tuple val("${intervals}"),val("${chrom}"),val("${outdir}"),val("${index}")
    script:
    interval="${intervals}".replaceAll(':|-','_')
    outdir="${maingatkpath}/GATKDBs/${interval}"
    chrom="${intervals}".(":")[0]
    logpath="${mainlogpath}/b_UpdateGENOMICSDBImport/"
    """
    #!/bin/bash
    set -euxo pipefail
    mkdir -p ${logpath}
    if [[ ! -e ${outdir} ]];
    then
        mkdir -p \$(dirname ${outdir})
    fi
    currentdb="${GATKpathtodbs}/${interval}"
    if [[ ! -e \${currentdb} ]];
    then
        echo "\${currentdb}" does not exist
        exit 1
    fi

    mapFilePath="mapfile.tsv"
    #copy currentdb to temp directory to maintain consistence
    rsync -rvP \${currentdb}/ "${interval}_db"
    #Build a sample map file to be passed to dbimport
    for s in ${gvcfs};
    do
    {
        f=\${s}
        fname=\$(echo \$s|sed -e "s/_SB2_.*//")
        if [[ "${index}" == "false" ]];
        then
            bgzip -f -d --keep \${f}
            ff=\${f/".gz"/""}
            gatk IndexFeatureFile --input \${ff}
            rm -f \${f}.tbi
            printf "\${fname}\\t\${ff}\\n">>\${mapFilePath}
        else
            printf "\${fname}\\t\${f}\\n">>\${mapFilePath}
        fi
    }
    done
    exec > dbimport.log 2>&1
    gatk --java-options "-Xmx30g -Xms30g -XX:ParallelGCThreads=${task.cpus}" GenomicsDBImport   \\
        --genomicsdb-update-workspace-path "${interval}_db" \\
        --sample-name-map \${mapFilePath} \\
        --consolidate true \\
        -OVI "${index}" \\
        --reader-threads $task.cpus ${otherGATKOptions}
    if [ \$? -gt 0 ]
    then
        echo "ERROR; Failed to run GATK GenomicsDBImport"
        exit 1
    fi
    #Copy updated db to new location
    rsync -rvP "${interval}_db/" ${outdir}
    rsync -rvP dbimport.log ${logpath}/${interval}.log
    """
}

/*
    Perform joint Genotyping on GenomicsDB workspace created with GenomicsDBImport
*/
process RUN_GENOTYPEGVCFs {
    tag "${intervals}"
    label "${params.gatkmulti}"
    label 'usescratch'
    label 'GATK'
    input:
        tuple val(chrom),val(intervals),val(GATKdatabase),val(index)
        val otherGATKoptions
        val refgenome
        val allele
    output:
        tuple val("${chrom}"), path("${outputfile}"),val("${index}")
    script:
    logpath="${mainlogpath}/c_GENOTYPEGVCFs/"
    intervalname="${intervals.replaceAll(':|-|__','_')}"
    outputfile="${intervalname}-genotypes_output.vcf.gz"
    """
    #!/bin/bash
    set -euxo pipefail
    mkdir -p ${logpath}
    echo Perform Joint genotyping of the following intervals ${intervals} on chromosome ${chrom} !
    intervalvcf="${intervals}"
    #Check if a reference vcf has been passed
    if [[ -e "${allele}" && "${params.GATKreferenceVARSonly}" == 'true' ]];
    then
        #forcecall=true
        intervalvcf="intervals.vcf"
        bcftools view --threads ${task.cpus} -r "${intervals}" -Ov -o \${intervalvcf} "${allele}.gz"
    fi
    echo Starting genotyping
    exec > genotype.log 2>&1

    gatk --java-options "-Xmx40G -Xms20G -XX:ParallelGCThreads=${task.cpus}" GenotypeGVCFs \\
        -R ${refgenome} \\
        -V gendb://"${GATKdatabase}" \\
        -O genotyped.vcf.gz \\
        -L "\${intervalvcf}" ${otherGATKoptions} \\
        -OVI ${index}

    echo start adding more details
    #Add more details to interval VCF file
    bcftools view --threads ${task.cpus} genotyped.vcf.gz -Ou \\
    | bcftools +setGT -- -t q -n . -i 'FMT/DP=0' \\
    | bcftools +fill-tags -- -t AC_Hom,AC_Het,AC_Hemi,MAF,F_MISSING,NS,TYPE,CR:1=1-F_MISSING \\
    | bcftools +setGT -- -t q -n . -i 'FMT/DP=0' \\
    | bcftools +tag2tag -- -r --PL-to-GL \\
    | bcftools annotate --threads ${task.cpus} --set-id '%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' -Ou \\
    | bcftools sort -Oz9 -o ${outputfile}

    if [[ "${index}" == "false" ]];
    then
        bgzip -f -d --keep "${outputfile}"
        aa="${outputfile}"
        gatk IndexFeatureFile --input \${aa/".gz"/}
    fi
    echo "Running genotyping for ${intervals} done"
    rsync -rvP genotype.log ${logpath}/${intervalname}.log
    """
}

process RUN_VariantRecalibrator
{
    input:
        tuple val(chrom), path(vcffile),val(index)

    output:
    script:
        logpath="${mainlogpath}/VARIANTRECALIBRATOR/"
        intervalname="${intervals.replaceAll(':|-|__','_')}"
        outputfile="${intervalname}-genotypes_output.vcf.gz"
    """
    #!/bin/bash
    set -euxo pipefail
    echo Starting variant recalibration
    exec > genotype.log 2>&1
    cd data/; \
    gatk --java-options -Xms4G -Xmx4G -XX:ParallelGCThreads=${task.cpus} VariantRecalibrator \\
      -tranche 100.0 -tranche 99.95 -tranche 99.9 \\
      -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \\
      -tranche 95.0 -tranche 94.0 \\
      -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \\
      -R ${params.refgenome} \\
      -V  \\
      --resource:hapmap,known=false,training=true,truth=true,prior=15.0 \
      /fdb/GATK_resource_bundle/hg38/hapmap_3.3.hg38.vcf.gz  \
      --resource:omni,known=false,training=true,truth=false,prior=12.0 \
      /fdb/GATK_resource_bundle/hg38/1000G_omni2.5.hg38.vcf.gz \
      --resource:1000G,known=false,training=true,truth=false,prior=10.0 \
      /fdb/GATK_resource_bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
      -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR  \
      -mode SNP -O merged_SNP1.recal --tranches-file output_SNP1.tranches \
      --rscript-file output_SNP1.plots.R
    """
}


