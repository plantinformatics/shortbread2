/*
    Running trimming algorithm for each library
*/
//Check if merged bam files are to be retained
def alignments="${launchDir}/Alignments/"
if(params.keepbams)
{
    alignments=params.runpath+"/Alignments/"
}
alignments=alignments+params.aligner+"/"
File logpath= new File(params.runpath+"/Logs/02_Alignment/")
if(!logpath.exists())
{
    logpath.mkdirs()
}

process RUN_ALIGNMENT {
    tag "${samplename}"
    //publishDir "${logpath}", mode: 'copy',pattern: '*.log',saveAS: "${samplename}.log"
    label "${params.alignlabel}"
    label 'usescratch'
    label 'alignment'
    input:
        tuple val(sampleid),val(samplename)
        val trimmethod
        path samplesheet
        path results
        val  trimmeroptions
        val OFFSET
        val aligner
        val refindex
        val typeofseq
        val aligneroptions
        val keeptrimmedfqs
        val skiptrimming
    output:
         tuple val(sampleid), path("${out_samsorted}"),env(merge)
    script:
        out_samsorted="${samplename}_sorted.bam"
        logs="${logpath}/${aligner}/${sampleid}/"
        trimlog="${logpath}/${trimmethod}/"
    """
    #!/bin/bash
    
    samplesheetrow=\$(cat ${samplesheet}|grep -v "SampleName"|awk -v sample=${samplename} -F "," '{if(\$2==sample) print \$0}')
    echo \${samplesheetrow}
    sampleid=\$(echo \${samplesheetrow}|awk -F "," '{printf("%s",\$1)}')
    merged_bam="${alignments}/\${sampleid}_sorted.bam"
    mdcheckfile="\${merged_bam}.md5sum"
    #Create log directories
    mkdir -p ${logs}
    mkdir -p ${trimlog}
    merge="yes"
    exec > "${logs}/${samplename}.log" 2>&1
    if [[ -e "\${mdcheckfile}" && "\$(md5sum -c \${mdcheckfile}|grep OK)" ]];
    then
        ln -s "\${merged_bam}" "${out_samsorted}"
        merge="no"
    else
        echo -e "Running trimming with ${trimmethod} on ${samplename} fastq files"
        seqtype=\$(echo \${samplesheetrow}|awk -F "," '{printf("%s",\$NF)}')
        #Format of samplesheet
        #SampleID,SampleName,FCID,Lane,SampleIndex,RGID,RGLB,RGPL,RGPU,RGSM,Read1,Read2,SEQType
        #Get GROUP header details
        RGID=\$(echo \${samplesheetrow}|awk -F "," '{printf("%s",\$6)}')
        LB=\$(echo \${samplesheetrow}|awk -F "," '{printf("%s",\$7)}')
        PL=\$(echo \${samplesheetrow}|awk -F "," '{printf("%s",\$8)}')
        PU=\$(echo \${samplesheetrow}|awk -F "," '{printf("%s",\$9)}')
        SM=\$(echo \${samplesheetrow}|awk -F "," '{printf("%s",\$10)}')
        #ReadGroup=\"\${ReadGroup}\"
        Read1=\$(echo \${samplesheetrow}|awk -F "," '{printf("%s",\$(NF-2))}')
        Read2=\$(echo \${samplesheetrow}|awk -F "," '{printf("%s",\$(NF-1))}')
        fast1_paired="${samplename}_R1_paired.fq.gz"
        fast1_single="${samplename}_R1_single.fq.gz"
        fast2_paired=""
        fast2_single=""
        if [[ "\${seqtype}" == "PE" ]];
        then
            fast2_paired="${samplename}_R2_paired.fq.gz"
            fast2_single="${samplename}_R2_single.fq.gz"
        fi
        case "${trimmethod}" in
            fastp)
                if [ "\${seqtype}" == "PE" ];
                then
                    fastp -i \${Read1} -I \${Read2} \\
                    -o \${fast1_paired} \\
                    -O \${fast2_paired} \\
                    --unpaired1 \${fast1_single} \\
                    --unpaired2 \${fast2_single} -l 35 -g \\
                    --detect_adapter_for_pe -j "${trimlog}/${samplename}_fastp.json" \\
                    -h "${trimlog}/${samplename}_fastp.html ${trimmeroptions}" \\
                    --thread $task.cpus
                else
                    fastp -i \${Read1} -o \${fast1_paired} \\
                    --unpaired1 \${fast1_single} -l 35 -g \\
                    --detect_adapter_for_pe -j "${trimlog}/${samplename}_fastp.json" \\
                    -h "${trimlog}/${samplename}_fastp.html" ${trimmeroptions} \\
                    --thread $task.cpus
                fi
            ;;
            trimmomatic)
                if [ "\${seqtype}" == "PE" ];
                then
                    java -jar \$EBROOTTRIMMOMATIC/trimmomatic-0.38.jar \${seqtype} -phred${OFFSET} \\
                            -trimlog "${trimlog}/${samplename}.log" \\
                            \${Read1} \${Read2} \${fast1_paired} \${fast1_single} \${fast2_paired} \${fast2_single} \\
                            ILLUMINACLIP:\$EBROOTTRIMMOMATIC/adapters/NexteraPE-PE.fa:2:30:3:1:true LEADING:20 TRAILING:20  \\
                            SLIDINGWINDOW:3:15  AVGQUAL:20 MINLEN:35 ${trimmeroptions}
                else
                    java -jar \$EBROOTTRIMMOMATIC/trimmomatic-0.38.jar \${seqtype} -phred${OFFSET} \\
                    -trimlog "${trimlog}/${samplename}.log" \\
                    \${Read1} \${fast1_paired} \${fast1_single} \\
                    ILLUMINACLIP:\$EBROOTTRIMMOMATIC/adapters/TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 \\
                    SLIDINGWINDOW:4:15 MINLEN:36 ${trimmeroptions}

                fi
            ;;
            *)
                echo "No valid trimmining method selected!"
                exit 1
            ;;
        esac

        if [[ "${keeptrimmedfqs}" == "true" ]];
        then
            trimpath="${alignments}/${trimmethod}/${sampleid}/"
            mkdir -p \${trimpath}
            rsync -rvP *fq.gz \${trimpath}
        fi
        #Perform alignment
        echo "Running alignmentment for ${samplename} with ${aligner}"
        paired_align="${samplename}_paired.bam"
        R1align="${samplename}_R1_single.bam"
        R2align="${samplename}_R2_single.bam"
        rg=\$(cat ${samplesheet}|grep ${samplename}|awk -v sample=${samplename} -F "," '{if(\$2==sample) printf("%s",\$(NF-3))}')
        threads="--threads $task.cpus"
        threads_align=$task.cpus
        if [[ ! -z "${aligneroptions}" ]]
        then
            threads_align="\${threads_align} ${aligneroptions} "
        fi


        case "${aligner}" in
            bwamem)
            {
              bwa mem -M -t \${threads_align} ${refindex} \${fast1_paired} \${fast2_paired} \\
              |samtools view - \$threads -Su \\
              |samtools sort - -o \${paired_align} \$threads

              #Check if run is paired end
              if [[ "\${seqtype}" == "PE" ]];
              then
                bwa mem -M -t \${threads_align} ${refindex} \${fast1_single} \\
                |samtools view - \$threads -Su \\
                |samtools sort - -o \${R1align} \$threads

                bwa mem -M -t \${threads_align} ${refindex} \${fast2_single} \\
                |samtools view - \$threads -Su \\
                |samtools sort - -o \${R2align} \$threads
              fi
            }
            ;;
            bwamem2)
            {
               bwa-mem2 mem -t \${threads_align} ${refindex} \${fast1_paired} \${fast2_paired} \\
               |samtools view - \$threads -Su \\
               |samtools sort - -o \${paired_align} \$threads

               if [[ "\${seqtype}" == "PE" ]];
               then
                 bwa-mem2 mem -t \${threads_align} ${refindex} \${fast1_single} \\
                 | samtools view - \$threads -Su \\
                 | samtools sort - -o \${R1align} \$threads

                 bwa-mem2 mem -t \${threads_align} ${refindex} \${fast2_single} \\
                 | samtools view - \$threads -Su \\
                 | samtools sort - -o \${R2align} \$threads
               fi

            }
            ;;
            bowtie2)
            {
                  bowtie2 \${threads} -x ${refindex} -1 \${fast1_paired} -2 \${fast2_paired} "${aligneroptions}" \\
                  | samtools view - \${threads} -Su \\
                  | samtools sort - -o \${paired_align} \${threads} -O BAM

                  if [[ "\${seqtype}" == "PE" ]];
                  then
                    bowtie2 -x "${refindex}" -U "\${fast1_single}" \$threads "${aligneroptions}" \\
                    | samtools view - \$threads -Su \\
                    | samtools sort - -o \${R1align} \$threads -O BAM \\

                    bowtie2 -x "${refindex}" -U "\${fast2_single}" \$threads "${aligneroptions}" \\
                    | samtools view - \$threads -Su \\
                    | samtools sort - -o \${R2align} \$threads -O BAM
                  fi
            }
            ;;
            star)
            {

                STAR --genomeDir ${refindex} \\
                     --runThreadN $task.cpus \\
                     --readFilesIn \${fast1_paired} \${fast2_paired}  \\
                     --outFileNamePrefix \${paired_align} \\
                     --outSAMtype BAM SortedByCoordinate \\
                     --outSAMunmapped Within \\
                     --twopassMode Basic \\
                     --readFilesCommand zcat "${aligneroptions}"
                  mv *sortedByCoord.out.bam \${paired_align}

                  if [[ "\${seqtype}" == "PE" ]];
                  then
                    STAR --genomeDir ${refindex} \\
                         --runThreadN $task.cpus \\
                         --readFilesIn \${fast1_single}  \\
                         --outFileNamePrefix \${R1align} \\
                         --outSAMtype BAM SortedByCoordinate \\
                         --outSAMunmapped Within \\
                         --twopassMode Basic \\
                         --readFilesCommand zcat "${aligneroptions}"
                    mv *.out.bam \${R1align}
                    STAR --genomeDir ${refindex} \\
                         --runThreadN $task.cpus \\
                         --readFilesIn \${fast2_single}  \\
                         --outFileNamePrefix \${R2align} \\
                         --outSAMtype BAM SortedByCoordinate \\
                         --outSAMunmapped Within \\
                         --twopassMode Basic \\
                         --readFilesCommand zcat "${aligneroptions}"
                    mv *.out.bam \${R2align}
                  fi
            }
            ;;
            subread)
            {
              subread-align -i ${refindex} \\
                -r \${fast1_paired} \\
                -R \${fast2_paired} \\
                -t "${typeofseq}" \\
                -o \${paired_align} \\
                --sortReadsByCoordinates \\
                -T $task.cpus "${aligneroptions}"
              if [[ "\${seqtype}" == "PE" ]];
              then
                subread-align -i ${refindex} \\
                -r \${fast1_single} \\
                -t "${typeofseq}" -o \${R1align} \\
                --sortReadsByCoordinates \\
                -T $task.cpus "${aligneroptions}"
                subread-align -i ${refindex} \\
                -r \${fast2_single} \\
                -t "${typeofseq}" \\
                -o \${R2align} \\
                --sortReadsByCoordinates \\
                -T $task.cpus "${aligneroptions}"
              fi
            }
            ;;
            minimap2)
            {
                minimap2 -t $task.cpus -ax sr ${refindex} \${fast1_paired} \${fast2_paired} \\
                |samtools view - \$threads -Su \\
                |samtools sort - -o \${paired_align} \$threads

                if [[ "\${seqtype}" == "PE" ]];
                then
                minimap2 -t $task.cpus -ax sr ${refindex} \${fast1_single} \\
                |samtools view - \$threads -Su \\
                |samtools sort - -o \${R1align} \$threads

                minimap2 -t $task.cpus -ax sr ${refindex} \${fast2_single} \\
                |samtools view - \${threads} -Su \\
                |samtools sort - -o \${R2align} \${threads}
                fi
            }
           ;;
            *)
            echo Choose the correct aligner
            exit 0
            ;;
        esac

        #If paired end then merge bam files generated from the split trimmed fastq files
        if [[ "\${seqtype}" == "PE" ]];
        then
        {
            #Get list of bam files from the trimmed files
            BAMList=\$(find . -type f \\( -name "*paired.bam" -o -name "*single.bam" \\) \\
            | sort \\
            | awk '{for(i=1;i<=NF;i++)printf("I=%s%s", \$i,FS)}END{printf("\\n")}')

            #Merge bam files to a single sample - level file
            picard MergeSamFiles \${BAMList} \\
                                    O=merged.bam \\
                                    ASSUME_SORTED=true \\
                                    MERGE_SEQUENCE_DICTIONARIES=true \\
                                    USE_THREADING=true \\
                                    VALIDATION_STRINGENCY=LENIENT \\
                                    CREATE_INDEX=false
            rm -f \${BAMList}
        }
        else
            mv  \${paired_align} merged.bam
        fi
        #Add Read group details to the the bam files
        picard AddOrReplaceReadGroups \\
                       I=merged.bam \\
                       O=${out_samsorted} \\
                       RGID=\${RGID} \\
                       RGLB=\${LB} \\
                       RGPL=\${PL} \\
                       RGPU=\${PU} \\
                       RGSM=\${SM}
    fi
    echo "End alignmentment for ${samplename}"
    """
}


/*
Merge bam files belonging to the same sample ID
*/
process MERGE_BAMS_BYSAMPLEID{
    tag "${sampleid}"
    afterScript "cp .command.log ${logs}/${sampleid}.log"
    label "${params.mergelabel}"
    label 'usescratch'
    label 'samtools'
    if("${params.seqtype}"=="rna")
    label 'samtoolsgatk'
    input:
         tuple val(sampleid),path(bamlist),val(merge)
         val markduplicates
         path results
         val samtoolsoptions
         val mappingquality
         val aligner
         val seqtype
         val refgenome
         val max_chromsize
    output:
          tuple val(sampleid), val("${bamout}"), val("${bamout_index}")
    script:
        outputdir="${alignments}"
        bamout="${outputdir}${sampleid}_sorted.bam"
        bamout_index="${bamout}.csi"
        indextype="-c"
        /*if(Long.parseLong(max_chromsize)>params.tabixmaxsize)
        {
            bamout_index="${bamout}.bai"
            indextype="-b"
        }*/
        logs="${logpath}/Merging_and_deduplication/${sampleid}/"
    """
    #!/bin/bash
    
    mkdir -p "${outputdir}"
    mkdir -p ${logs}

    exec > ${logs}/${sampleid}.log 2>&1

    if [[ "${merge}" == "yes" && ! -e "${bamout}" ]]
    then
    {
        BAMList=\$(find . -name "*.bam")
        #Merge bam files
        if [[ \$(echo ${bamlist}|wc -w) -gt 1 ]]
        then
            BAMList=\$(echo ${bamlist}|sort -n|awk '{for(i=1;i<=NF;i++)printf("I=%s%s", \$i,FS)}END{printf("\\n")}')
            picard MergeSamFiles \${BAMList} \\
                O=${sampleid}.bam \\
                ASSUME_SORTED=true \\
                MERGE_SEQUENCE_DICTIONARIES=true \\
                USE_THREADING=true \\
                VALIDATION_STRINGENCY=LENIENT \\
                CREATE_INDEX=false
        else
            cp ${bamlist} ${sampleid}.bam
        fi

        #Sort the merged aligned bam file
        samtools sort -@ $task.cpus -o ${sampleid}.sorted.bam ${sampleid}.bam
        samtools index -@ $task.cpus -c ${sampleid}.sorted.bam ${sampleid}.sorted.bam.csi
        samtoolsview="-q ${mappingquality} -@ $task.cpus -b"

        if [[ ! -z "${samtoolsoptions}" ]];
        then
            samtoolsview="\${samtoolsview} ${samtoolsoptions}"
        fi

        samtools view \${samtoolsview} ${sampleid}.sorted.bam \\
        |samtools view -@ $task.cpus -o ${sampleid}.sorted-1.bam
        rm -f ${sampleid}.sorted.bam.csi
        mv ${sampleid}.sorted-1.bam ${sampleid}.sorted.bam

        samtools index -c -@ $task.cpus ${sampleid}.sorted.bam ${sampleid}.sorted.bam.csi
        if [[ "${markduplicates}" == "false" ]]
        then
            gatk --java-options "-Xms60G -Xmx60G -XX:+UseParallelGC -XX:ParallelGCThreads=$task.cpus" MarkDuplicatesSpark \\
                  -I ${sampleid}.sorted.bam \\
                  -O ${sampleid}.sorted.dedup.bam  \\
                  -M ${logs}/${sampleid}.metrics \\
                  -OBI false ${params.markduplicatesoptions}
            mv ${sampleid}.sorted.dedup.bam  ${sampleid}.sorted.bam
        fi

        if [[ "${seqtype}" == "rna" ]];
        then
            gatk --java-options "-Xmx30g -Xms30g -XX:ParallelGCThreads=${task.cpus}" SplitNCigarReads \\
                -R ${refgenome} \\
                -I ${sampleid}.sorted.bam \\
                -OBI false \\
                -O ${sampleid}.sorted.gatk.bam
            mv ${sampleid}.sorted.gatk.bam  ${sampleid}.sorted.bam
        fi

        #Perform base recalibration
        if [[ ! -z "${params.GATKknowsitesrecal}" ]];
        then
        {
            gatk --java-options "-Xms4G -Xmx4G -XX:ParallelGCThreads=${task.cpus}" BaseRecalibrator \\
              -I ${sampleid}.sorted.bam \\
              -R ${refgenome} \\
              -O ${sampleid}_bqsr.report \\
              --known-sites ${params.GATKknowsitesrecal} ${params.GATKbaserecaloptions}

            gatk --java-options "-Xms2G -Xmx2G -XX:ParallelGCThreads=${task.cpus}" ApplyBQSR \\
              -I ${sampleid}.sorted.bam \\
              -R ${refgenome} \\
              --bqsr-recal-file ${sampleid}_bqsr.report \\
              --create-output-bam-index false ${params.GATKApplyBQSRoptions} \\
              -O ${sampleid}_bqsr.bam
            rsync --remove-source-files -rvP ${sampleid}_bqsr.report ${logs}/${sampleid}_bqsr.report
            mv ${sampleid}_bqsr.bam ${sampleid}.sorted.bam
        }
        fi
        #Compute alignment stats
        samtools stats -@ $task.cpus ${sampleid}.sorted.bam > ${logs}/${sampleid}.stats
        #Copy files to working directory/resultsdir
        rsync --remove-source-files -rvP ${sampleid}.sorted.bam ${bamout}
        samtools index "${indextype}" -@ $task.cpus ${bamout}
        md5sum ${bamout} > "${bamout}.md5sum"
        #cat .command.log
        echo "End merging of alignmentments for ${sampleid}"
    }
    fi
    """
}
