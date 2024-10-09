rawdata=$1
samplesheet=$2
results=$3
seqtype=$4
RGPL=$5
mode=$6

#Make required folders
if [[ ! -d ${results} ]]; then
  mkdir -p ${results}
fi
SEQtype=${seqtype}
fastq1s=""
sampleIds=""
if [ -z ${samplesheet} ]; then
  fastq1s=$(find ${rawdata} -name "*[_][R]1[_][0][0][1].fastq*")
  sampleIds=$(for f in ${fastq1s}; do basename $f | sed -e "s/_[A-Z]*[0-9]*_L00[0-9]_R[12][_][0][0][1].fastq[.][g][z]//g"; done)
else
  #check the number of columns
  num_columns=$(head -n1 ${samplesheet} | awk -F ',' '{print NF}')
  if [[ ${num_columns} -lt 2 ]]; then
    echo "The samplesheet requires at least 2 columns in the format - SampleID,Read1 "
    exit 5
  fi
  fastq1s=$(cat ${samplesheet} | cut -d"," -f2)
  sampleIds=$(cat ${samplesheet} | cut -d"," -f1)
fi
sample_details="${results}/Sample_number_of_reads.csv"
echo "SampleID,SampleName,FCID,Lane,RGPL,SampleIndex,PU,ReadGroup,ReadGroup2,Read1,Read2,SEQType" >samplesheet_from_files.csv
echo "SampleID,SampleName,R1.ReadCount,R2.ReadCount" >"${sample_details}"
sampleIds=$(echo $sampleIds | awk '{for (i=1;i<=NF;i++) if (!a[$i]++) printf("%s%s",$i,FS)}{printf("\\n")}')
if [[ "${mode}" == "test" ]]; then
  arr=(${sampleIds})
  list_length=${#arr[@]}
  idx=$((56824 % ${list_length}))
  idx=$(($idx + 1))
  r1=$(echo $sampleIds | awk -v idx=${idx} '{printf("%s",$idx)}')
  r2=""
  if [[ ${list_length} -gt 1 ]]; then
    idx2=$((568243 % ${list_length} + 1))
    if [[ ${idx2} -eq ${idx} || ${idx2} -gt ${list_length} ]]; then
      idx2=${list_length}
    fi
    r2=$(echo $sampleIds | awk -v idx=${idx2} '{printf("%s",$idx)}')
  fi
  fastq1s=$(echo ${fastq1s} | awk -v ran1="${r1}" -v ran2="${r2}" '{for(i=1;i<=NF;i++){if(match($i,ran1)||match($i,ran2))printf("%s%s",$i,FS)}}')
  if [ ! -z ${samplesheet} ]; then
    fastq1s=$(cat ${samplesheet} | grep "${r1}\\|${r2}" | cut -d"," -f2)
  fi
fi

for read1 in ${fastq1s}; do
  {
    echo ${read1}
    if [[ ! -e ${read1} ]]; then
      continue
    fi
    read2=""
    if [[ "${seqtype}" == "PE" || "${seqtype}" == "mixed" ]]; then
      if echo $read1 | grep "[_][R]1[_][0-9][0-9][0-9].fastq"; then
        read2=$(echo ${read1} | sed -e "s/_R1_/_R2_/g")
      elif echo $read1 | grep "1.fastq"; then
        read2=$(echo ${read1} | sed -e "s/1.fastq/2.fastq/g")
      fi
      #Check if read2 exists
      if [[ ! -e "${read2}" ]]; then
        read2=""
        SEQtype="SE"
      else
        SEQtype="PE"
      fi
    fi

    FCID=$(zcat ${read1} | head -n 1 | cut -d ':' -f 3)
    lane=$(zcat ${read1} | head -n 1 | cut -d ':' -f 4)
    sampleIDX=$(zcat ${read1} | head -n 1 | awk -F ":" '{print $NF}')
    sampleID=$(basename ${read1} | sed -e "s/[_-][Ll][0-5]*_R1.*//g")
    sample=$(echo -n ${read1} | md5sum | cut -d" " -f1)
    sample="$(basename ${read1/.fastq*//})_${sample}"
    sample="$(echo ${sample} | sed -e 's/-/_/g')"
    if [[ "$(basename ${read1} | grep [a-zA-Z0-9_\\-]*_[A-Z]*[0-9]*_L00[0-9]_R[12]_001.fastq.gz)" ]]; then
      echo "File name suggests Illumina sequencing, getting sample name"
      sampleID=$(basename ${read1} | sed "s/_[A-Z]*[0-9]*_L00[0-9]_R[12]_001.fastq.gz//g")
    fi
    f1_reads=$(echo $(zcat ${read1} | wc -l)/4 | bc)
    f2_reads=$(echo $(zcat ${read2} | wc -l)/4 | bc)
    echo "${sampleID},${sample},${f1_reads},${f2_reads}" >>"${sample_details}"
    if [[ ! -z "${samplesheet}" ]]; then
      sampleID=$(cat ${samplesheet} | grep ${read1} | cut -d"," -f1)
    fi
    pu="${FCID}.${lane}.${sampleIDX}"
    #'@RG\tID:ID\tSM:SM\tLB:LB\tPU:PU\tPL:PL'
    RG="@RG\\\\tID:${FCID}.${lane}\\\\tSM:${sampleID}\\\\tLB:${sample}\\\\tPU:${pu}\\\\tPL:${RGPL}"
    RG2="--rg-id ${FCID}.${lane} --rg SM:${sampleID} --rg LB:${sample} --rg PU:${pu} --rg PL:${RGPL}"
    echo "${sampleID},${sample},${FCID},${lane},${RGPL},${sampleIDX},${pu},${RG},${RG2},${read1},${read2},${SEQtype}" >>samplesheet_from_files.csv
  } &
  wait
done
