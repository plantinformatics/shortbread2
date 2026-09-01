def tabixmaxsize = 530000000 //Current maximum chrom size used by Htslib tabix function during indexation

File mainlogpath = new File(params.runpath + "/Logs/04_mpileup/") //Path to mpileup log files
File mainmpileuppath=new File(params.runpath+"/mpileup/") //Path to output of mpileup
File veppath=new File(params.runpath+"/VEP/") //Path to output of Variant Effect Predictor

def bamoutspath="${launchDir}/Alignments/mpileup/"
if(params.keepbams)
{
    bamoutspath=params.runpath+"/Alignments/mpileup/"
}

if(!mainlogpath.exists())
    mainlogpath.mkdirs()
if(!mainmpileuppath.exists())
    mainmpileuppath.mkdirs()

process RUN_GENOTYPE_MPILEUP {
    tag "${sampleid}_${intervalname}"
    label "${params.mpileupmulti}"
    label 'usescratch'
    label 'BCFTOOLS'
    input:
        tuple val(intervals), val(sampleid), path(alignments),path(alignments_index)
        val refgenome
        val otheroptions
        val maxchromsize
        val allele
    output:
        tuple val("${intervals}"), val("${chrom}"),path("${outputfile}"), path("${outputfileidx}"),env(index)
    script:
        mpileupLogPath = "${mainlogpath}/c_GENOTYPE_MPILEUP/${sampleid}/"
        bamsout=""
        if(params.GATKbamout)
        {
            bamsout="-bamout ${intervalname}_GATK.bam"
        }
    intervalname = "${intervals.replaceAll(':|-|__','_')}"
    outputfile = "${sampleid}_SB2_${intervalname}-genotypes_output.vcf.gz"
    outputfileidx = "${sampleid}_SB2_${intervalname}-genotypes_output.vcf.gz.csi"
    chrom="${intervals}".split(":")[0]
    """
    #!/bin/bash
    set -euxo pipefail
    mkdir -p ${mpileupLogPath}
    echo "bcftools calling for intervals ${intervals} on chromosome ${chrom}"

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
    fi

    intervalvcf="${intervals}"
    GVCF_BANDS="1,2,3,4,5,6,7,10,15,30"
    echo "Starting mpileup + call"
    exec > >(tee -i genotype.log) 2>&1
    bcftools mpileup \\
        -Ou \\
        -f ${refgenome} \\
        ${alignments} \\
        --annotate FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR \\
        -r "\${intervalvcf}" \\
    | bcftools call \\
        --gvcf \${GVCF_BANDS} -m -Ou \\
    | bcftools sort -Oz -o ${outputfile}

    pipeline_status=( "\${PIPESTATUS[@]}" )
    set -e

    if (( pipeline_status[0] != 0 ||
          pipeline_status[1] != 0 ||
          pipeline_status[2] != 0 )); then

        echo "ERROR: paired-end bcftools pipeline failed"
        echo "mpileup exit status: \${pipeline_status[0]}"
        echo "bcftools call exit status: \${pipeline_status[1]}"
        echo "bcftools sort exit status: \${pipeline_status[2]}"

        if (( pipeline_status[0] == 137 ||
              pipeline_status[1] == 137 ||
              pipeline_status[2] == 137 )); then
            echo "A pipeline component was killed, treating this as a probable OOM"
            exit 137
        fi

        exit 1
    fi



    bcftools index -f ${outputfile} -o ${outputfileidx}

    echo "bcftools calling for ${sampleid}_${intervals} complete"
    rsync -rvP genotype.log ${mpileupLogPath}/${sampleid}_${intervalname}.log

    """
}

process RUN_MERGE_MPILEUP {
    tag "${intervals}"
    label "${params.gatkmulti}"
    label 'usescratch'
    label 'BCFTOOLS'
    input:
        tuple val(intervals), val(chrom), path(gvcfs),path(gvcfsidx), val(index)
        val sample_order
    output:
        tuple val("${chrom}"), val("${outfile}"),val("${index}")
    script:
    interval="${intervals}".replaceAll(':|-','_')
    outdir="${mainmpileuppath}/mpileup/${interval}"
    outfile="${mainmpileuppath}/mpileup/${interval}.vcf.gz"
    chrom="${intervals}".split(":")[0]
    logpath="${mainlogpath}/d_RUN_MERGE_MPILEUP/"
    """
    #!/bin/bash
    set -euxo pipefail
    echo "Running merging on ${intervals} belonging to chromosome ${chrom}!"
    mkdir -p \$(dirname "${outdir}")
    mkdir -p "${logpath}"


    echo "Starting merge"
    exec > merge.log 2>&1
    merge_tmp_dir="merge_batches"
    mkdir -p "\${merge_tmp_dir}"
    trap 'rm -rf "\${merge_tmp_dir}"' EXIT

    # Create one input VCF path pr file
    echo "${gvcfs}" | tr ' ' '\n' > "\${merge_tmp_dir}/all_gvcfs.list"

    num_gvcfs=\$(wc -l < "\${merge_tmp_dir}/all_gvcfs.list")
    echo "Number of input VCFs: \${num_gvcfs}"

    if [[ "\${num_gvcfs}" -eq 0 ]]; then
        echo "No input VCFs were supplied"
        exit 1
    fi

    if [[ "\${num_gvcfs}" -le 200 ]]; then
        echo "Fewer than 200 vcfs, no issues with merging. Merging \${num_gvcfs} input VCFs directly"

        final_merge_inputs="\${merge_tmp_dir}/all_gvcfs.list"

    else
        echo "Too many vcfs for single pass, splitting \${num_gvcfs} VCFs into batches of 200"

        split \\
            -l 200 \\
            -d \\
            -a 4 \\
            "\${merge_tmp_dir}/all_gvcfs.list" \\
            "\${merge_tmp_dir}/batch_"

        intermediate_list="\${merge_tmp_dir}/intermediate_bcf.list"
        > "\${intermediate_list}"

        batch_number=0

        for batch_list in "\${merge_tmp_dir}"/batch_*; do
            batch_number=\$((batch_number + 1))

            intermediate_bcf="\${merge_tmp_dir}/intermediate_\$(printf '%04d' "\${batch_number}").bcf"
            batch_size=\$(wc -l < "\${batch_list}")

            echo "Merging batch \${batch_number} containing \${batch_size} VCFs"

            bcftools merge \\
                --file-list "\${batch_list}" \\
                --merge none \\
                --threads ${task.cpus} \\
                -Ob \\
                -o "\${intermediate_bcf}"

            # Intermediate BCFs need new index for merge.
            bcftools index \\
                --force \\
                --threads ${task.cpus} \\
                "\${intermediate_bcf}"

            printf "%s\\n" "\${intermediate_bcf}" >> "\${intermediate_list}"
        done

        num_intermediates=\$(wc -l < "\${intermediate_list}")
        echo "Created \${num_intermediates} intermediate BCFs"

        final_merge_inputs="\${intermediate_list}"
    fi

    echo "Running final merge and variant processing"

    bcftools merge \\
        --file-list "\${final_merge_inputs}" \\
        --merge none \\
        -Ou \\
    | bcftools +setGT -- -t q -n . -i 'FMT/DP=0' \\
    | bcftools view --threads ${task.cpus} -v snps,indels,mnps -Ou \\
    | bcftools +fill-tags -- -t AC_Hom,AC_Het,AC_Hemi,MAF,F_MISSING,NS,TYPE,CR:1=1-F_MISSING \\
    | bcftools +tag2tag -- -r --PL-to-GL \\
    | bcftools annotate --threads ${task.cpus} --set-id '%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' -Ou \\
    | bcftools sort -Oz -o "${outfile}"

    # In case sample order differs between chrom chunks we standardise now bsed on sample sheet
    printf "%s\n" ${sample_order.join(' ')} > desired_samples.txt
    bcftools query -l "${outfile}" > current_samples.txt
    grep -F -x -f current_samples.txt desired_samples.txt > reorder_samples.txt
    bcftools view --threads ${task.cpus} -S reorder_samples.txt -Oz -o "${outfile}.reordered.vcf.gz" "${outfile}"
    mv "${outfile}.reordered.vcf.gz" "${outfile}"
    touch "${outfile}"

    bcftools index -f "${outfile}"

    echo " BCFtools merging for ${intervals} done"
    rsync -rvP merge.log ${logpath}/${interval}.log
    """
}

process GATHER_VCF_MPILEUP {
    tag "${chrom}"
    label 'usescratch'
    //label "${params.gathervcfs}"
    label "verylarge"
    input:
        tuple val(chrom),path(vcfs),val(indextype)
        val refgenome
        val type
        path addedmetadata
    output:
        tuple val("${index}"),val("${genotypes}"),val("${chrom}"),val("${type}")
    script:
    outdir="${mainmpileuppath}/VCFs/"
    genotypes="${outdir}/${chrom}-${type}.vcf.gz"
    statsfile="${mainlogpath}/${chrom}-${type}.stats"
    index=indextype[0]
    """
    #!/bin/bash
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
    bcftools view --header-only ${genotypes} > fileheaderinfo.txt
    
    gawk '
    FNR == NR {

        # Store assembly line
        if (\$0 ~ /^##assembly=/) {
            assembly = \$0
        }

        # Store genome_url line
        if (\$0 ~ /^##genome_url/) {
            genome_url = \$0
        }

        # Start capturing metadata chunk
        if (\$0 ~ /##Shortbread2_analysis_start_date:/) {
            in_chunk = 1
        }

        # Store chunk lines, including start and end lines
        if (in_chunk) {
            chunk[++chunk_n] = \$0
        }

        # Stop capturing after Variant_call_method line
        if (\$0 ~ /##Shortbread2_variant_call_method:/) {
            in_chunk = 0
        }

        # Store contig-specific extra metadata from file2
        if (match(\$0, /^##contig=<ID=([^,>]+)/, m)) {
            id = m[1]

            if (match(\$0, /,species.*>/)) {
                contig_extra[id] = substr(\$0, RSTART, RLENGTH)
            }
        }

        next
    }


    # Insert metadata chunk below ##fileformat= as key word find
    /^##fileformat=/ {
        print

        for (i = 1; i <= chunk_n; i++) {
            print chunk[i]
        }

        next
    }

    # Insert assembly and genome_url below ##reference= as key word find
    /^##reference=/ {
        print

        if (assembly != "") {
            print assembly
        }

        if (genome_url != "") {
            print genome_url
        }

        next
    }

    # Update matching contig lines while preserving the length already reported so it remains the correct length in case of edge cases
    /^##contig=<ID=/ {

        if (match(\$0, /^##contig=<ID=([^,>]+)/, m)) {
            id = m[1]

            if (id in contig_extra) {
                line = \$0

                # Remove closing > from file1 contig line before appending metadata
                sub(/>\$/, "", line)

                print line contig_extra[id]
                next
            }
        }
    }

    {
        print
    }
    ' ${addedmetadata} fileheaderinfo.txt > updated_header.txt


    bcftools reheader -h updated_header.txt -o ${genotypes}.TMP ${genotypes} && mv "${genotypes}.TMP" "${genotypes}"




    """
}