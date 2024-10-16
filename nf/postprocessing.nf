def tabixmaxsize = 530000000 //Current maximum chrom size used by Htslib tabix function during indexation

File mainlogpath=new File(params.runpath+"/Logs/05_VCF_STATS/") //Path to VCF log files
File maingatkpath=new File(params.runpath+"/GATK/") //Path to output of GATK
File veppath=new File(params.runpath+"/VEP/") //Path to output of Variant Effect Predictor
if(!mainlogpath.exists())
    mainlogpath.mkdirs()
if(!maingatkpath.exists())
    maingatkpath.mkdirs()

/*
    Join interval level vcf files to chromosome level and perform filtering on VCF file
*/

process GATHER_VCF {
    tag "${chrom}"
    label 'usescratch'
    //label "${params.gathervcfs}"
    label "verylarge"
    input:
        tuple val(chrom),path(vcfs),val(indextype)
        val refgenome
        val type
    output:
        tuple val("${index}"),val("${genotypes}"),val("${chrom}"),val("${type}")
    script:
    outdir="${maingatkpath}/VCFs/"
    genotypes="${outdir}/${chrom}-${type}.vcf.gz"
    statsfile="${mainlogpath}/${chrom}-${type}.stats"
    index=indextype[0]
    """
    #!/bin/bash
    set -euxo pipefail
    echo "Gathering files belonging to chromosome ${chrom} together" !
    mkdir -p ${outdir}
    ls *.vcf.gz| split -l 100 - subset_vcfs
    for i in subset_vcfs*;
    do
    {
        vcfs_i=\$(cat \$i | tr '\\n' ' ')
        bcftools concat --threads $task.cpus --output-type v \${vcfs_i} \\
        |bcftools sort --output \${i}.vcf.gz
        bcftools index -c \${i}.vcf.gz
    }
    done

    #Combine the subsets and sort/index in one step
    bcftools concat --threads $task.cpus --output-type v subset_vcfs*.vcf.gz \\
    |bcftools sort -Oz9 -o "${genotypes}"
    bcftools stats --threads ${task.cpus} "${genotypes}" > "${statsfile}"


    # Conditional indexing (simplified)
    [[ "${index}" == "true" ]] && bcftools index -t "${genotypes}" || bcftools index -c "${genotypes}"

    vcftools --gzvcf ${genotypes} --missing-indv --out ${genotypes}
    #Plot missing genotypes

    #cat .command.log >> "${mainlogpath}/GATHER_VCF.log"
    """
}

/*
Process to filter called variants
*/
process FILTER_VAR {
    tag "${genotypes.baseName}"
    label "small"
    label 'usescratch'
    input:
        val CR
        val MAF
        val AC
        val MQ
        val otherfilters
        tuple val(chrom),path(genotypes),val(index)
        val samplecount
        val keepmultiallelic
        val gatkreferencevcf
        val keepindels
    output:
        tuple val("${chrom}"),path("${genotypefilt}"),val("${index}")
    script:
    genotypefilt="${genotypes.toString().replace('.vcf.gz','-filtered_with-')}MAF_${MAF}-CR_${CR}-AC_${AC}.vcf.gz"
    """
    #!/bin/bash
    set -euxo pipefail
    echo "Start filtering for ${genotypes} on ${CR}"
    bcftools index ${genotypes}
    if [[ ${samplecount} -gt 1 ]];
    then
        CR_comp=${CR}
        if [[ ${CR} -gt 1 ]]
        then
            CR_comp=\$(echo "scale=10; ${CR}/100"|bc)
        fi
        referencevcfintersect=""
        #Filter final VCF file using filtering options provided by the user
        filt_snps="MAF>=${MAF} && CR>=\${CR_comp} && AC>=${AC} && MQ > ${MQ}"
        if [[ ! -z "${otherfilters}" ]];
        then
            filt_snps="\${filt_snps} && ${otherfilters}"
        fi

        #Split vcf to SNPs and INDELs
        bcftools filter -Ou -i "TYPE=='SNP'" ${genotypes} \\
        |bcftools view --threads ${task.cpus} -Oz9 -i "\${filt_snps}" -o snps-filtered.vcf.gz
        bcftools index snps-filtered.vcf.gz
        if [[ "${keepindels}" == "yes" ]];
        then
            filt_indels="CR>=\${CR_comp}"
            bcftools view --threads ${task.cpus} -Ou -e "TYPE=='SNP'" ${genotypes} \\
            |bcftools view --threads ${task.cpus} -Oz -i "\${filt_indels}" -o indels-filtered.vcf.gz
            bcftools index indels-filtered.vcf.gz
        fi

        cohortfiltered="${chrom}-genotypes-filtered.vcf.gz"
        cohort="${chrom}-genotypes-multi.vcf.gz"

        bcftools concat -a -d exact --threads ${task.cpus} -Ou *-filtered.vcf.gz \\
        |bcftools sort -Oz9 -o \${cohort}

        bcftools index \${cohort}
        if [[ "${keepmultiallelic}" != "yes"  ]];
        then
            bcftools view --threads ${task.cpus} -m 2 -M 2 -O z9 \${cohort} -o \${cohortfiltered}
        else
            bcftools annotate --set-id +'%CHROM\\_%POS\\_%REF\\_%ALT' \${cohort} -O z9 -o \${cohortfiltered}
        fi
        bcftools index \${cohortfiltered}
        #Retain all positions in reference vcf file
        if [[ -n "${gatkreferencevcf}" && -f "${gatkreferencevcf}" ]];
        then
             # Reference VCF exists, filter and process
            bcftools view \\
                -R "${gatkreferencevcf}" \\
                --threads "${task.cpus}" \\
                -Oz9 \\
                -o Positions_in_reference.vcf.gz \\
                "${genotypes}"

            bcftools index Positions_in_reference.vcf.gz

            bcftools concat \\
                -a -d exact \\
                --threads "${task.cpus}" \\
                Positions_in_reference.vcf.gz \\
                "\${cohortfiltered}" \\
                | bcftools sort -Oz9 \\
                -o "${genotypefilt}"
        else
            #Generate stats for the filtered vcf file
            mv \${cohortfiltered} ${genotypefilt}
        fi

        if [[ "${index}" == "true" ]];
        then
            bcftools index -t ${genotypefilt}
        else
            bcftools index -c ${genotypefilt}
        fi
        echo "Filtering comple on ${chrom}"
    else
        echo "No filtering performed on the one sample"
    fi
    #cat .command.log >> "${mainlogpath}/FILTER_VAR.log"
    """
}

/*
Process to generate a list of variant positions without
Genotype information
*/
process GENERATE_VariantList{
    tag "${file(genotypes).baseName}"
    label "small"
    label 'usescratch'
    errorStrategy 'finish'
    input:
        val genotypes
    script:
    snplist="${genotypes.replaceAll('.vcf.gz','-List.vcf.gz')}"
    """
    #!/bin/bash
    set -euxo pipefail
    echo Generating variant list for "${genotypes}.baseName"
    echo "Start generating SNP list for ${genotypes} before final step"
    bcftools view --threads ${task.cpus} -G -O z9 -o SNPlist.vcf.gz ${genotypes}
    rsync -rvP SNPlist.vcf.gz ${snplist}
    bcftools index ${snplist}
    """
}

process GENERATE_VARIANTGRAPH
{
    tag "${chrom}"
    label 'usescratch'
    input:
        tuple val(index),val(genotypes),val(chrom),val(type)
    output:
        val "${output}"
    script:
    File output=new File(params.runpath+"/Variant_Graphs/${chrom}/${type}/")
    if(!output.exists())
        output.mkdirs()
    """
    #!/bin/bash
    set -euxo pipefail
    chromfasta="${chrom}.fasta"
    samtools faidx ${params.refgenome} ${chrom}>\${chromfasta}
    samtools faidx \${chromfasta}
    module load shifter
    shifter --image=donchbio/variantstore /opt/variantstore/variantstore construct \\
    -r \${chromfasta} \\
    -v ${genotypes} \\
    -p ${output}
    """
}


process VariantEffect
{
       input:
        path genotypes
        val species
       output:
       vep="${veppath}/${genotypes.replaceAll('.vcf.gz','_VEP.vcf')}"
       script:
       """
        vep -i input.txt -o output.txt --species triticum_aestivum --database --genomes
       """
}
