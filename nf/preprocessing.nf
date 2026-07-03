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
            which picard
            echo Creating dictionary
            export _JAVA_OPTIONS="-Xmx40G -Xms20G -XX:+UseParallelGC -XX:ParallelGCThreads=${task.cpus}"
            picard CreateSequenceDictionary R=${refgenome} O=\${dict}
        fi
        max_chromsize=\$(cat \${dict}|grep -v @HD |awk '{print \$3}'|sed -e "s/LN://g"|sort -nr|head -n1)
        min_chromsize=\$(cat \${dict}|grep -v @HD |awk '{print \$3}'|sed -e "s/LN://g"|sort -nr|tail -n1)
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
        intervallist=\$(awk -v intsize=${intervalsize} -v exclude="\${excludecontigs}" -v chrexc=\$(echo ${chromtoexclude}|sed -e "s/\\*//g") '{split(\$2,chr,":"); if(exclude=="yes" && substr(chr[2],1,length(chrexc))==chrexc){next};if(chr[2]==chrexc){next};split(\$3,len,":");skip=exclude=="true"&&len[2]<intsize;if(!skip){for(i=1;i<=len[2];i+=intsize){j=i+intsize-1; if(j>len[2]){j=len[2]}; printf("%s:%d-%d ", chr[2], i, j )}}}' ${dict})
        #intervallist=\$(grep "^@SQ" ${dict}|awk -v intsize=${intervalsize} -v exclude="\${excludecontigs}" -v chrexc=\$(echo ${chromtoexclude}|sed -e "s/\\*//g") '{split(\$2,chr,":"); if(exclude=="yes" && substr(chr[2],1,length(chrexc))==chrexc){next};if(chr[2]==chrexc){next};split(\$3,len,":");printf("%s:%d-%d ", chr[2], 1, len[2])}')
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
    numofcols=14 //Expected number of columns in a normalized samplesheet generated by shortbread2
    oldnumofcols=13 //Legacy normalized samplesheet without RunID
    """
    #!/bin/bash
    exec > ${logpath}/1_Prepare_samplesheet.log 2>&1
    numlibs=0 #Variable to hold the number of libraries generated in samplesheet
    mkdir -p "${results}" #Directory where the results from shortbread are saved

    SEQtype="PE"
    sampleIds=""
    skip="no"
    read1_col=2
    declare -a fastq1s=()
    declare -A run_count

    if [ -z "${samplesheetinput}" ];
    then
        echo "No Samplesheet has been provided"
        exit 1
    else
        num_columns=\$(head -n1 "${samplesheetinput}" | awk -F ',' '{print NF}')

        if [[ \${num_columns} -lt 2 ]]; then
            echo "ERROR: Samplesheet requires either 2 columns (SampleID, Read1) or 3 columns (SampleID, RunID, Read1)." >&2
            head "${samplesheetinput}" >&2
            exit 1
        elif [[ \${num_columns} -eq ${numofcols} ]]; then
            echo "Copying existing normalized samplesheet '${samplesheetinput}' to '${outsheet}'."
            cp "${samplesheetinput}" "${outsheet}"
            exit 0
        elif [[ \${num_columns} -eq ${oldnumofcols} ]]; then
            echo "Converting legacy normalized samplesheet '${samplesheetinput}' to include RunID."
            awk -F',' 'BEGIN{OFS=","}
                NR==1{print "SampleID","RunID","SampleName","FCID","Lane","SampleIndex","RGID","RGLB","RGPL","RGPU","RGSM","Read1","Read2","SEQType"; next}
                {count[\$1]++; print \$1,"Run" count[\$1],\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13}' "${samplesheetinput}" > "${outsheet}"
            exit 0
        elif [[ \${num_columns} -gt 3 ]]; then
            echo "ERROR: Samplesheet requires either 2 columns (SampleID, Read1) or 3 columns (SampleID, RunID, Read1)." >&2
            head "${samplesheetinput}" >&2
            exit 1
        fi

        if [[ \${num_columns} -eq 3 ]]; then
            read1_col=3
        fi

        cat "${samplesheetinput}" > "${samplesheet}"
        dos2unix "${samplesheet}"

        dups=\$(tail -n +2 "${samplesheet}" | cut -d',' -f\${read1_col} | sort | uniq -d | wc -l)
        if [[ \${dups} -gt 0 ]]; then
            echo "The sample sheet provided contains duplicate Read 1 files as shown below"
            tail -n +2 "${samplesheet}" | cut -d',' -f\${read1_col} | sort | uniq -d
            exit 1
        fi

        if [[ \${num_columns} -eq 3 ]]; then
            run_dups=\$(tail -n +2 "${samplesheet}" | awk -F',' '{print \$1","\$2}' | sort | uniq -d | wc -l)
            if [[ \${run_dups} -gt 0 ]]; then
                echo "The sample sheet provided contains duplicate SampleID,RunID pairs as shown below"
                tail -n +2 "${samplesheet}" | awk -F',' '{print \$1","\$2}' | sort | uniq -d
                exit 1
            fi
        fi

        mapfile -t fastq1s < <(awk -F',' -v col=\${read1_col} 'NR>1 {print \$col}' "${samplesheet}")
        sampleIds=\$(awk -F',' 'NR>1 {print \$1}' "${samplesheet}")

        while IFS= read -r line || [[ -n "\$line" ]]; do
        {
          line=\$(printf '%s' "\$line" | tr -d '\r')
          [[ -z "\$line" ]] && continue
          ID=""
          runID=""
          read1=""

          if [[ \${num_columns} -eq 3 ]]; then
              IFS=',' read -r ID runID read1 _ <<< "\$line"
          else
              IFS=',' read -r ID read1 _ <<< "\$line"
          fi
          [[ -z "\$ID" ]] && continue

          if [[ -z "\${runID// }" ]]; then
              count=\${run_count["\$ID"]:-0}
              count=\$((count+1))
              run_count["\$ID"]=\$count
              runID="Run\${count}"
          fi

          read1_md5=\$(printf '%s' "\$read1" | md5sum | awk '{print \$1}')
          if [[ -e "\$read1" ]]; then
            echo "\$ID,\$runID,\$read1,\$read1_md5" >> samplesheet.tmp
          fi
        }
        done < <(tail -n +2 "${samplesheet}")
    fi

    sampleIds=\$(echo \$sampleIds | awk '{for (i=1;i<=NF;i++) if (!a[\$i]++) printf("%s%s",\$i,FS)} END{print ""}')

    if [[ "${mode}" == "test" ]]; then
        arr=(\${sampleIds})
        list_length=\${#arr[@]}
        randoms=(\$(shuf --random-source='copy_of_original_samplesheet.csv' -i 0-\$((\${list_length} - 1)) -n 2))
        index1=\${randoms[0]}
        index2=\${randoms[1]}
        echo \$index1 \$index2
        r1="\${arr[\${index1}]}"
        r2="\${arr[\${index2}]}"
        if [[ -n "${samplesheet}" ]]; then
            mapfile -t fastq1s < <(awk -F',' -v ran1="\${r1}" -v ran2="\${r2}" -v col=\${read1_col} 'NR>1 && (\$1==ran1 || \$1==ran2) {print \$col}' "${samplesheet}")
        fi
    fi

    checkstring=\$(for f in "\${fastq1s[@]}"; do echo "\$f"; done | sort | md5sum | cut -d' ' -f1)

    if [[ -e "${md5samplecheck}" && -e "${results}/${outsheet}" ]]; then
        if [[ "\${checkstring}" == "\$(cat "${md5samplecheck}")" ]]; then
            skip="yes"
            cat "${results}/${outsheet}" > "${outsheet}"
        fi
    fi

    if [[ "\${skip}" == "no" ]]; then
        echo "SampleID,RunID,SampleName,FCID,Lane,SampleIndex,RGID,RGLB,RGPL,RGPU,RGSM,Read1,Read2,SEQType" > "${outsheet}"

        for read1 in "\${fastq1s[@]}"; do
        {
            if [[ ! -f "\${read1}" ]]; then
                echo "\${read1} not found"
                continue
            fi

            read2=""
            if echo "\$read1" | grep "[_][R]1[_][0-9][0-9][0-9].fastq" >/dev/null; then
                read2=\$(echo "\${read1}" | sed -e 's/_R1_/_R2_/g')
            elif echo "\$read1" | grep '1.fastq' >/dev/null; then
                read2=\$(echo "\${read1}" | sed -e 's/1.fastq/2.fastq/g')
            elif echo "\$read1" | grep '1.fq' >/dev/null; then
                read2=\$(echo "\${read1}" | sed -e 's/1.fq/2.fq/g')
            fi

            if [[ ! -e "\${read2}" ]]; then
                read2=""
                SEQtype="SE"
            else
                SEQtype="PE"
            fi

            read1md5=\$(echo -n "\${read1}" | md5sum | cut -d' ' -f1)
            sampleID=\$(awk -F',' -v md5="\${read1md5}" '\$4==md5 {print \$1; exit}' samplesheet.tmp | sed -e 's/ /-/g')
            runID=\$(awk -F',' -v md5="\${read1md5}" '\$4==md5 {print \$2; exit}' samplesheet.tmp | sed -e 's/ /-/g')

            if [[ -z "\${sampleID}" ]]; then
                echo "Sample ID is blank for \${read1}"
                exit 1
            fi
            if [[ -z "\${runID}" ]]; then
                runID="Run1"
            fi

            readdata=\$(zcat "\${read1}" 2>/dev/null | head -n 1 || true)
            FCID=\$(echo "\${readdata}" | cut -d ':' -f 3)
            lane=\$(echo "\${readdata}" | cut -d ':' -f 4)
            sampleIDX=\$(echo "\${readdata}" | awk -F ':' '{print \$NF}')

            [[ -z "\${FCID}" ]] && FCID="unknown"
            [[ -z "\${lane}" ]] && lane="unknown"
            [[ -z "\${sampleIDX}" ]] && sampleIDX="unknown"

            sampleName="\${sampleID}__\${runID}"
            sampleName=\$(echo "\${sampleName}" | sed -E 's/[[:space:]]+/-/g; s/[^A-Za-z0-9_.-]/_/g')
            RGLB=\$(basename "\${read1/.fastq*//}")_\${read1md5}
            RGLB=\$(echo "\${RGLB}" | sed -e 's/-/_/g')

            rgSuffix=""
            if [[ "\${FCID}" == "unknown" || "\${lane}" == "unknown" || "\${sampleIDX}" == "unknown" ]]; then
                rgSuffix=".\${read1md5:0:8}"
            fi

            RGPU="\${FCID}.\${lane}.\${sampleIDX}\${rgSuffix}"
            RGID="\${FCID}.\${lane}\${rgSuffix}"
            RGSM="\${sampleID}"

            echo "\${sampleID},\${runID},\${sampleName},\${FCID},\${lane},\${sampleIDX},\${RGID},\${RGLB},${RGPL},\${RGPU},\${RGSM},\${read1},\${read2},\${SEQtype}" >> "${outsheet}"
            numlibs=\$((\$numlibs+1))
        }
        done
    else
        numlibs=\$(wc -l < "${outsheet}")
        numlibs=\$((\$numlibs -1))
    fi

    echo "\${checkstring}" > "${md5samplecheck}"

    if [[ "\${numlibs}" -eq 0 ]]; then
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
    
    exec > ${logpath}/2_Prepare_genome.log 2>&1


    # Safety cleanup: only rm -rf if ref is >1 level deeper than refgenome dir. Required for star which produces sub folders if it fails
    # Safe clean up will run a secondary check to confirm \${ref} is in fact a sub directory of where the reference genome is set. 
    # May be unnessary I just don't like having an recursive remove without a safety check in case of nonsense. 
    # If functionality for indexes to be built away from reference genomes, will need to revise this approach or remove the safety net. Seems that this is only required for star  
    safe_cleanup() {
        local ref_path="\$1"
        local rg_path="\$2"

        # Resolve absolute directories (works even if final files don't exist)
        local ref_dir
        local rg_dir
        ref_dir=\$(readlink -f "\$(dirname "\$ref_path")")
        rg_dir=\$(readlink -f "\$(dirname "\$rg_path")")

        if [[ "\$ref_dir" == "\$rg_dir" ]]; then
            echo "[WARN] Ref prefix is in the SAME directory as the reference genome."
            echo "[WARN] Will delete FILES only: \${ref_path}*"
            rm -f "\${ref_path}"* || true
        elif [[ "\$ref_dir" == "\$rg_dir"/* ]]; then
            echo "[INFO] Ref is at least one level deeper under the reference genome directory."
            echo "[INFO] Removing files AND directories: \${ref_path}*"
            rm -rf "\${ref_path}"* || true
        else
            # Different tree altogether; be conservative
            echo "[WARN] Ref is NOT under the reference genome directory."
            echo "[WARN] Will delete FILES only: \${ref_path}*"
            rm -f "\${ref_path}"* || true
        fi
    }

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
                # Added if to check for silent fail status and if failed to then remove the reference index partial so it can be rerun, finally exit in error code allowing for retry at increase RAM
                if ! bwa index -p \${ref} ${refgenome}; then
                    echo "Indexing failed, cleaning up for retry of \${ref}*"
                    rm -f \${ref}*
                    exit 137
                fi
            ;;
            bwamem2)
                if ! bwa-mem2 index -p \${ref} ${refgenome}; then
                    echo "Indexing failed, cleaning up for retry of \${ref}*"
                    rm -f \${ref}*
                    exit 137
                fi
            ;;
            bowtie2)
                if ! bowtie2-build ${refgenome} \${ref} --threads $task.cpus --large-index; then
                    echo "Indexing failed, cleaning up for retry of \${ref}*"
                    rm -f \${ref}*
                    exit 137
                fi
            ;;
            star)
                if [ ! -d \${ref} ]; then
                    mkdir -p \${ref}
                fi
                limit_ram=${ task.memory ? task.memory.toBytes() : 68800807520 }
                # Adding in a more dynamic memory limit for RAM in case genome is too large and so will increase with increasing task attempt.
                echo "[INFO] Using --limitGenomeGenerateRAM=\${limit_ram} bytes"
                if ! STAR --runThreadN $task.cpus \\
                --runMode genomeGenerate \\
                --genomeDir \${ref} \\
                --genomeFastaFiles ${refgenome} \\
                --sjdbGTFfile ${refannotation} \\
                --limitGenomeGenerateRAM \${limit_ram}; then 
                    echo "Indexing failed, cleaning up for retry of \${ref}*"
                    # safe cleanup required because of star producing sub directories in index not just file.
                    safe_cleanup "\${ref}" "${refgenome}"
                    exit 137
                fi
            ;;
            subread)
                if ! subread-buildindex -o \${ref} ${refgenome}; then
                    echo "Indexing failed, cleaning up for retry of \${ref}*"
                    rm -f \${ref}*
                    exit 137
                fi
            ;;
            minimap2)
                if ! minimap2 -t $task.cpus -d \${ref} ${refgenome}; then
                    echo "Indexing failed, cleaning up for retry of \${ref}*"
                    rm -f \${ref}*
                    exit 137
                fi
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
    
    echo running split intervals
    cp .command.log ${logpath}/4_Split_intervals.log
    """
}


