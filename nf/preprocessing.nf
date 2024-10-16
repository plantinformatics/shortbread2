File logpath = new File(params.runpath+"/Logs/00_Preparation")
if(!logpath.exists())
    logpath.mkdirs()

def excludecontigs=params.excludesmallchrs
/*
    Process to get the maximum chromosome size and a list of chromosomes
*/
process GET_GENOMEINFO {
    tag "${file(refgenome).baseName}"
    label 'samtools'
    label 'varstore'
    executor 'local'
    input:
        val refgenome
        val GATKreferenceVCF
    output:
        tuple env(max_chromsize),
            env(chromosomes),
            env(min_chromsize),
            env(genomesize),
            env(dict),
            val(refvcf)
    script:
    resultsdir=params.runpath
    tabixmaxsize=params.tabixmaxsize
    refvcf=""
    if(new File(GATKreferenceVCF).exists())
        refvcf="${resultsdir}/${file(GATKreferenceVCF).baseName}.vcf"
    """
        #!/bin/bash
        set -euxo pipefail
        exec > ${logpath}/5_Get_chromsize.log 2>&1
        #Check if fasta file has an index
        fastaindex="${refgenome}.fai"
        if [[ ! -e \${fastaindex} ]];
        then
            echo indexing fasta file
            samtools faidx ${refgenome}
        fi
        #Check if reference dictionary exists, if not, build one
        dict=\$(echo ${refgenome}|sed -e "s/.[^.]*\$/.dict/")
        if [[ ! -e "\${dict}" ]]
        then
            code="CreateSequenceDictionary R=${refgenome} O=\${dict}"
            picard --java-options "-Xmx40G -Xms20G -XX:+UseParallelGC -XX:ParallelGCThreads=${task.cpus}" \${code}
        fi
        max_chromsize=\$(cat \${dict}|grep -v @HD |sort -k3 -r|head -n1|awk '{print \$3}'|sed -e "s/LN://g")
        min_chromsize=\$(cat \${dict}|grep -v @HD |sort -k3 -r|tail -n1|awk '{print \$3}'|sed -e "s/LN://g")
        chromosomes=\$(cat \${dict}|grep -v @HD|awk '{print \$2}'|sed -e "s/SN://")
        genomesize=\$(cat \${dict}|grep -v @HD |sed -e "s/LN://"|awk '{sum+=\$3;} END{print sum;}')
        echo "Checking if reference VCF has been provided"
        if [[ -e ${GATKreferenceVCF} ]];
        then
            bcftools view --threads ${task.cpus} -Ou ${GATKreferenceVCF} \\
            |bcftools sort -Ov -o ${refvcf}
            bgzip -kf ${refvcf}
            bcftools index -f ${refvcf}.gz
            gatk IndexFeatureFile -I ${refvcf}
        else
            refvcf="null"
        fi
    """
}

/*
   This process generates an interval list based on the size provided by the user
   The interval list is stored in memory rather than being written to file
   if the interval size provided by the user is greater than the maximum chromosome size, 
   then the interval size used is equivalent to the minimum chromosome size
*/
process PREPARE_INTERVALS{
    tag "${file(refgenome).baseName}"
    label 'samtools'
    executor 'local'
    input:
        val refgenome
        val intervalsize
        tuple val(max_chromsize),val(chromosomes),val(min_chromsize),val(genomsize),path(dict),val(refvcf)
        val mode
        val GATKpathtodbs
        val chromtoexclude
    output:
        env intervallist
    script:
    """
    #!/bin/bash
    set -euxo pipefail
    exec > ${logpath}/3_Prepare_intervals.log 2>&1
    intervallist=""
    if [[ "${GATKpathtodbs}" != "" ]];
    then
        #Get interval names used in building the database
            intervallist=\$(find ${GATKpathtodbs} -name "*\\\$*"|sed 's/\\(.*\\)\\\$/\\1-/'|sed 's/\\(.*\\)\\\$/\\1:/'|xargs echo|xargs -n1 basename)
    else
        excludecontigs="no"
        if [[ ${chromtoexclude} == *"*"* ]];
        then
            excludecontigs="yes"
        fi
        #Build interval list based on the dictionary from the reference genome
        #intervallist=\$(awk -v intsize=${intervalsize} -v exclude="\${excludecontigs}" -v chrexc=\$(echo ${chromtoexclude}|sed -e "s/\\*//g") '{split(\$2,chr,":"); if(exclude=="yes" && substr(chr[2],1,length(chrexc))==chrexc){next};if(chr[2]==chrexc){next};split(\$3,len,":");skip=exclude=="true"&&len[2]<intsize;if(!skip){for(i=1;i<=len[2];i+=intsize){j=i+intsize-1; if(j>len[2]){j=len[2]}; printf("%s:%d-%d ", chr[2], i, j )}}}' ${dict})
        intervallist=\$(grep "^@SQ" ${dict}|awk -v intsize=${intervalsize} -v exclude="\${excludecontigs}" -v chrexc=\$(echo ${chromtoexclude}|sed -e "s/\\*//g") '{split(\$2,chr,":"); if(exclude=="yes" && substr(chr[2],1,length(chrexc))==chrexc){next};if(chr[2]==chrexc){next};split(\$3,len,":");printf("%s:%d-%d ", chr[2], 1, len[2])}')
        if [[ -z "\${intervallist}" ]];
        then
            echo "\${intervallist}"
            echo "Interval list is empty, check interval size or set excludesmallchrs=false"
            exit 1
        fi
    fi
    if [[ "${mode}" == "test" ]];
    then
        intervallist=\$(echo \${intervallist}|awk '{for(i=1;i<=1;i++) printf("%s ",\$i)}')
    fi
    """
}



/*
Prepare samplesheet that includes sample IDs and sample names
- extract read group details from read 1
- if a samplesheet is provided then check if it has 2 columns [sample ID, read1]
*/
process PREPARE_SAMPLE_SHEET {
    publishDir "${results}",mode:'copy'
    tag "${file(samplesheetinput).baseName}"
    label 'samplesheet'
    executor 'local'
    input:
        val rawdata
        val samplesheetinput
        val results
        val RGPL
        val mode
    output:
        path "${outsheet}"
    script:
    outsheet="samplesheet_from_files.csv" //name of samplesheet saved in outdir
    md5samplecheck="${results}/.samplesheet.md5check" //md5 checksum used to check if to generate a new samplesheet
    samplesheet="copy_of_original_samplesheet.csv" //A copy of the samplesheet supplied by user is made to deal with write access issues
    numofcols=13 //Expected number of columns in a samplesheet generated by shortbread2
    """
    #!/bin/bash
    #set -euxo pipefail
    exec > ${logpath}/1_Prepare_samplesheet.log 2>&1
    numlibs=0 #Variable to hold the number of libraries generated in samplesheet
    #Make required folders
    mkdir -p "${results}" #Directory where the results from shortbread are saved
    SEQtype="PE"
    fastq1s=""
    sampleIds=""
    skip="no" #Skip generation of samplesheet
    #Build list of fastq1 files
    if [ -z ${samplesheetinput} ];
    then
        #Exit pipeline
        echo "No Samplesheet has been provided"
        exit 1
        #fastq1s=\$(find ${rawdata} -name "*[_][R]1[_][0][0][1].fastq*");
        #sampleIds=\$(for f in \${fastq1s};do basename \$f|sed -e "s/_[A-Z]*[0-9]*_L00[0-9]_R[12][_][0][0][1].fastq[.][g][z]//g"; done)
    else
         #check the number of columns
        num_columns=\$(head -n1 ${samplesheetinput}|awk -F ',' '{print NF}')
        if [[ \${num_columns} -lt 2 ]]; then
            echo "ERROR: Samplesheet requires at least 2 columns (SampleID, Read1)." >&2
            head "${samplesheetinput}" >&2
            exit 1
        elif [[ \${num_columns}  -eq ${numofcols} ]]; then
            echo "Copying existing samplesheet '${samplesheetinput}' to '${outsheet}'."
            cp "${samplesheetinput}" "${outsheet}"
            exit 0  # Exit successfully if the samplesheet is already complete
        fi


        #Make a copy of the original samplesheet and remove any windows encodings
        cat ${samplesheetinput}>${samplesheet}
        dos2unix ${samplesheet} #remove any dos related encodings

        #Check if samplesheet contains duplicates
        dups=\$(cat $samplesheet|cut -d"," -f2|sort|uniq -d|wc -l)

        if [[ \${dups} -gt 0 ]];
        then
            echo "The sample sheet provided contains duplicate Read 1 files as shown below"
            cat ${samplesheet}|cut -d"," -f2|sort|uniq -d
            exit 1
        fi

        #Build list of fastq1 files and sample IDs
        fastq1s=\$(cat ${samplesheet}|awk -F ',' '{printf("%s ",\$2)}')
        sampleIds=\$(cat ${samplesheet}|cut -d"," -f1)

        while read -r line;
        do
        {
            ID=\$(echo \$line|cut -d"," -f1)
            read1=\$(echo \$line|cut -d"," -f2)
            read1_md5=\$(echo -n \${read1}|md5sum|cut -d" " -f1)
            if [[ -e "\${read1}" ]];
            then
                echo "\${ID},\${read1},\${read1_md5}" >> samplesheet.tmp
            fi
        }
        done< ${samplesheet}
    fi
    sampleIds=\$(echo \$sampleIds|awk '{for (i=1;i<=NF;i++) if (!a[\$i]++) printf("%s%s",\$i,FS)}{printf("\\n")}')

    #If run mode is test, then sample 2 random samples
    if [[ "${mode}" == "test" ]];
    then
        arr=(\${sampleIds})
        list_length=\${#arr[@]}
        randoms=(\$(shuf --random-source='copy_of_original_samplesheet.csv' -i 0-\$((\${list_length} - 1)) -n 2))
        index1=\${randoms[0]}
        index2=\${randoms[1]}
        echo \$index1 \$index2
        # Extract the corresponding sampleIDs using the random indices
        r1="\${arr[\${index1}]}"
        r2="\${arr[\${index2}]}"
        fastq1s=\$(echo \${fastq1s}|awk -v ran1="\${r1}" -v ran2="\${r2}" '{for(i=1;i<=NF;i++){if(match(\$i,ran1)||match(\$i,ran2))printf("%s%s",\$i,FS)}}')
        if [ ! -z ${samplesheet} ];
        then
            fastq1s=\$(cat ${samplesheet}|grep "\${r1}\\|\${r2}"|cut -d"," -f2)
        fi
    fi


    #string to build hash value used to check if a samplesheet has already been created
    checkstring=\$(echo -n "\${fastq1s}"|sort|md5sum|cut -d" " -f1)

    if [[ -e "${md5samplecheck}" && -e "${results}/${outsheet}" ]]
    then
        if [[ "\${checkstring}" -eq "\$(cat ${md5samplecheck})" ]];
        then
            skip="yes"
            cat "${results}/${outsheet}" > ${outsheet}
        fi
    fi

    #Build sample sheet if one does not exist
    if [[ "\${skip}" == "no" ]];
    then
        #Write header of sample sheet
        
        echo "SampleID,SampleName,FCID,Lane,SampleIndex,RGID,RGLB,RGPL,RGPU,RGSM,Read1,Read2,SEQType">${outsheet}
        #RGID; RGLB; RGPL;RGPU; RGSM

        for read1 in \${fastq1s};  #Loop through each read1 fastq file and read FCID, lane, and sample index
        do
        {
            if [[ ! -f "\${read1}" ]]
            then
                echo "\${read1} not found"
                continue
            fi
            read2=""

            if echo \$read1|grep  "[_][R]1[_][0-9][0-9][0-9].fastq";
            then
                read2=\$(echo \${read1}|sed -e "s/_R1_/_R2_/g")
            elif echo \$read1|grep "1.fastq";
            then
                 read2=\$(echo \${read1}|sed -e "s/1.fastq/2.fastq/g")
            fi

            #Check if read2 exists
            if [[ ! -e "\${read2}" ]];
            then
              read2=""
              SEQtype="SE"
            else
              SEQtype="PE"
            fi
            readdata=\$(zcat \${read1}|head -n 1)
            FCID=\$(echo \${readdata}|cut -d ':' -f 3) #Read flow cell ID
            lane=\$(echo \${readdata}|cut -d ':' -f 4) #Read lane number from first read
            sampleIDX=\$(echo \${readdata}|awk -F ":" '{print \$NF}') #Read sample index
            sampleID=\$(basename \${read1}|sed -e "s/[_-][Ll][0-5]*_R1.*//g")
            read1md5=\$(echo -n \${read1}|md5sum|cut -d" " -f1)
            sample="\$(basename \${read1/.fastq*//})_\${read1md5}"
            sample="\$(echo \${sample}|sed -e 's/-/_/g')"

            if [[ ! -z "${samplesheet}" ]];
            then
                sampleID=\$(cat samplesheet.tmp|grep "\${read1md5}"|cut -d"," -f1|sed -e "s/ /-/g")
                if [[ -z \${sampleID} ]];
                then
                    echo "Sample ID is blank for \${read1}"
                    exit 1
                fi
            else
                #Extract sample ID if the naming of files is in illumina format
                if [[ "\$(basename \${read1}|grep [a-zA-Z0-9_\\-]*_[A-Z]*[0-9]*_L00[0-9]_R[12]_001.fastq.gz)" ]];
                then
                    sampleID=\$(basename \${read1}|sed "s/_[A-Z]*[0-9]*_L00[0-9]_R[12]_001.fastq.gz//g")
                fi
            fi

            #Format Group ID details that are needed by GATK
            RGPU="\${FCID}.\${lane}.\${sampleIDX}"
            #'@RG\tID:ID\tSM:SM\tLB:LB\tPU:PU\tPL:PL'
            RG="@RG\\\\tID:\\\\tSM:\${sampleID}\\\\tLB:\${sample}\\\\tPU:\${RGPU}\\\\tPL:${RGPL}"
            RG2="--rg-id \${FCID}.\${lane} --rg SM:\${sampleID} --rg LB:\${sample} --rg PU:\${RGPU} --rg PL:${RGPL}"
            RGID="\${FCID}.\${lane}"
            RGLB="\${sample}"
            RGSM="\${sampleID}"

            #Write to samplesheet
            #SampleID,SampleName,FCID,Lane,SampleIndex,RGID,RGLB,RGPL,RGPU,RGSM,Read1,Read2,SEQType
            echo "\${sampleID},\${sample},\${FCID},\${lane},\${sampleIDX},\${RGID},\${RGLB},${RGPL},\${RGPU},\${RGSM},\${read1},\${read2},\${SEQtype}" >>${outsheet}
            numlibs=\$((\$numlibs+1))
        }
        done
    else
        echo "\${checkstring}" > "${md5samplecheck}"
        numlibs=\$(wc -l "${outsheet}")
        numlibs=\$((\$numlibs -1))
    fi

    #Check that the samplesheet generated contains more than 1 file
    if [[ "\${numlibs}" -eq 0 ]];
    then
            echo "Files in samplesheet not found, please check that paths to read1 are valid"
            exit 1
    fi
    """
}

/*
 - Prepare reference genome depending on the aligner used
 -
*/
process PREPARE_GENOME{
    tag "${file(refgenome).baseName}"
    label 'alignment'
    if(params.genomeexecutorlocal)
        executor 'local'
    else
        label 'medium'
    input:
        val refgenome
        val refannotation
        val aligner
        val results
    output:
        env ref
    script:
    if("${refannotation}" == "")
        refannotation="-"
    """
    #!/bin/bash
    set -euxo pipefail
    exec > ${logpath}/2_Prepare_genome.log 2>&1
    #check if reference genome fasta file has been indexed
    if [ ! "\$(ls ${refgenome}.fai)" ];
    then
        echo "building fasta index file for reference genome"
        samtools faidx ${refgenome}
    fi

    #expected reference index
    ref="${params.refindex}"
    if ls \${ref}* >/dev/null 2>&1;
    then
        echo "Index already exists, skipping index building"
    else
        echo "Building index for ${aligner} in \${ref}"
        #Build index based on choice of aligner
        case "${aligner}" in
            bwamem)
                bwa index -p \${ref} ${refgenome}
            ;;
            bwamem2)
                bwa-mem2 index -p \${ref} ${refgenome}
            ;;
            bowtie2)
                bowtie2-build ${refgenome} \${ref} --threads $task.cpus --large-index
            ;;
            star)
                mkdir -p \${ref}
                STAR --runThreadN $task.cpus \\
                --runMode genomeGenerate \\
                --genomeDir \${ref} \\
                --genomeFastaFiles ${refgenome} \\
                --sjdbGTFfile ${refannotation} \\
                --limitGenomeGenerateRAM 38800807520
            ;;
            subread)
                subread-buildindex -o \${ref} ${refgenome}
            ;;
            minimap2)
                minimap2 -t $task.cpus -d \${ref} ${refgenome}
            ;;
            *)
            echo Choose the correct aligner
            exit 1
            ;;
        esac
    fi
    """
}

process SPLIT_INTERVALS()
{
    tag "${interval}"
    executor 'local'
    input:
        val interval
        val GATKpathtodbs
        val intervalsizegatk
    output:
         tuple val("${interval}"),val("${intevs}")
    script:
      chrom=interval.split(":")[0]
      start=interval.split(":")[1].split("-")[0].toLong()
      end=interval.split(":")[1].split("-")[1].toLong()
      l=end-start
      intevs=""
      def intervalsizegatk2=intervalsizegatk/params.splitintervalsgatk
      if(l>intervalsizegatk2)
      {
        for(i=start;i<=end;i+=intervalsizegatk2)
        {
            j=i+intervalsizegatk2-1
            if(j>end)
            {
                j=end
            }
            new_interval=chrom+":"+i+"-"+j
            if(intevs=="")
                intevs=new_interval
            else
                intevs+="="+new_interval
        }
      }
      else
          intevs=interval
    """
    #!/bin/bash
    set -euxo pipefail
    echo running split intervals
    cp .command.log ${logpath}/4_Split_intervals.log
    """
}


