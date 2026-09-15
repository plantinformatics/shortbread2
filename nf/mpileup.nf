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

    tag "${intervalname}"

    label "${params.mpileupmulti}"
    label 'usescratch'
    label 'BCFTOOLS'

    input:
        val(intervals)
        path(alignments)
        path(alignmentsIndexes)
        path refgenome
        val otheroptions
        val maxchromsize
        val allele
    output:
        tuple val("${chrom}"), path("${outputfile}"), path("${outputfileidx}"), val("${indexAvailable}")
    script:
        def maxchromsizeconverted = maxchromsize
            .toString()
            .replaceAll(',', '')
            .trim()
            .toLong()
        indexAvailable = maxchromsizeconverted < tabixmaxsize
        intervalname = intervals
            .toString()
            .replaceAll(':|-|__', '_')
        chrom = intervals
            .toString()
            .split(':')[0]
        outputfile = "${intervalname}.vcf.gz"
        outputfileidx = "${outputfile}.csi"
        mpileupLogPath = "${mainlogpath}/c_GENOTYPE_MPILEUP/"
        bamListArguments = alignments
            .collect { bam -> "\"${bam.toString()}\"" }
            .join(' ')
        indexAvailable = maxchromsizeconverted < tabixmaxsize

        if (alignments.size() != alignmentsIndexes.size()) {
            error """
            BAM/index count mismatch for interval ${intervals}:
              BAM files: ${alignments.size()}
              indexes:   ${alignmentsIndexes.size()}
            """.stripIndent()
        }

        intervalname = intervals
            .toString()
            .replaceAll(/:|-|__/, '_')

        chrom = intervals
            .toString()
            .split(':')[0]

        """
        #!/bin/bash
        set -uo pipefail

        mkdir -p "${mpileupLogPath}"

        exec > >(tee -i genotype.log) 2>&1

        echo "Starting joint bcftools calling"
        echo "Interval: ${intervals}"
        echo "Chromosome: ${chrom}"
        echo "Number of BAMs: ${alignments.size()}"

        intervalvcf="${intervals}"

        printf '%s\\n' ${bamListArguments} > bam_files.list

        num_bams=\$(grep -cve '^[[:space:]]*\$' bam_files.list)

        if [[ "\${num_bams}" -eq 0 ]]; then
            echo "ERROR: no BAM files were supplied"
            exit 1
        fi

        echo "Number of BAMs: \${num_bams}"
        echo "Open-file soft limit: \$(ulimit -Sn)"

        set +e

        bcftools mpileup \
            --threads ${task.cpus} \
            -Ou \
            -f "${refgenome}" \
            --bam-list bam_files.list \
            --annotate FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR \
            -r "\${intervalvcf}" \
        | bcftools call \
            --threads ${task.cpus} \
            -m \
            -v \
            -Ou \
        | bcftools +setGT \
            -Ou \
            -- -t q -n . -i 'FMT/DP=0' \
        | bcftools view \
            --threads ${task.cpus} \
            -v snps,indels,mnps \
            -Ou \
        | bcftools +fill-tags \
            -Ou \
            -- -t 'AC_Hom,AC_Het,AC_Hemi,MAF,F_MISSING,NS,TYPE,CR:1=1-F_MISSING' \
        | bcftools +tag2tag \
            -Ou \
            -- -r --PL-to-GL \
        | bcftools annotate \
            --threads ${task.cpus} \
            --set-id '%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' \
            -Ou \
        | bcftools sort \
            -Oz \
            -o "${outputfile}"

        pipeline_status=( "\${PIPESTATUS[@]}" )

        set -e

        component_names=(
            "bcftools mpileup"
            "bcftools call"
            "bcftools +setGT"
            "bcftools view"
            "bcftools +fill-tags"
            "bcftools +tag2tag"
            "bcftools annotate"
            "bcftools sort"
        )

        pipeline_failed=0

        for i in "\${!pipeline_status[@]}"; do
            status="\${pipeline_status[\${i}]}"

            echo "\${component_names[\${i}]} exit status: \${status}"

            if [[ "\${status}" -ne 0 ]]; then
                pipeline_failed=1

                if [[ "\${status}" -eq 137 ]]; then
                    echo "Pipeline component was killed; probable OOM"
                    exit 137
                fi
            fi
        done

        if [[ "\${pipeline_failed}" -ne 0 ]]; then
            echo "ERROR: joint bcftools pipeline failed"
            exit 1
        fi

        if [[ ! -s "${outputfile}" ]]; then
            echo "ERROR: output VCF was not created"
            exit 1
        fi

        bcftools index \\
            --force \\
            --csi \\
            --threads ${task.cpus} \\
            "${outputfile}"

        if [[ ! -s "${outputfileidx}" ]]; then
            echo "ERROR: output CSI index was not created"
            exit 1
        fi

        final_sample_count=\$(
            bcftools query -l "${outputfile}" |
            wc -l
        )

        if [[ "\${final_sample_count}" -ne "\${num_bams}" ]]; then
            echo "ERROR: VCF sample count does not match BAM count"
            echo "BAM count: \${num_bams}"
            echo "VCF sample count: \${final_sample_count}"
            exit 2
        fi

        echo "Final sample count: \${final_sample_count}"
        echo "Number of variant records:"
        bcftools index -n "${outputfile}"

        echo "Joint calling for ${intervals} complete"

        rsync \
            -rvP \
            genotype.log \
            "${mpileupLogPath}/${intervalname}.log"
        """
}

process GATHER_VCF_MPILEUP {
    tag "${chrom}"
    label 'usescratch'
    //label "${params.gathervcfs}"
    label "verylarge"
    publishDir "${mainmpileuppath}/VCFs/", mode: 'copy', overwrite: true
    input:
        tuple val(chrom),path(vcfs),val(indextype)
        val refgenome
        val type
        path addedmetadata
    output:
        tuple val("${index}"),path("${genotypes}"),val("${chrom}"),val("${type}")
    script:
    genotypes="${chrom}-${type}.vcf.gz"
    statsfile="${mainlogpath}/${chrom}-${type}.stats"
    index=indextype[0]
    """
    #!/bin/bash
    echo "Gathering files belonging to chromosome ${chrom} together" !
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

    # Recreate the index after reheadering so it matches the final VCF bytes.
    [[ "${index}" == "true" ]] && bcftools index -f -t "${genotypes}" || bcftools index -f -c "${genotypes}"
    """
}