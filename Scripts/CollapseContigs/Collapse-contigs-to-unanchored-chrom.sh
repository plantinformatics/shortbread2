#!/bin/bash
function helppage() {
  echo
  echo
  echo -e "Usage: sh $(basename $0) -R Path to fasta file [-p Starting pattern of standard chromosomes] or [-l min length of contig to exclude from collapsing] "
  echo
  echo -e "\t\t[options]"
  echo -e "\t\t-R Absolute path to fasta file"
  echo -e "\t\t-p Starting pattern to chromosome names e.g chr in chr1"
  echo -e "\t\t-l minimum length of contig/scaffold to be excluded from unanchored chromosome"
  echo -e "\t\t-h for help"
  echo
}

if (! getopts "R:p:l:h" option); then
  helppage
  exit $E_OPTERROR
fi

fasta="" #path to fasta file
chrompattern="none" #$2 #Starting name for standard chromosomes
chromlength=0
HELP="no"

while getopts R:p:l:h option; do
  case "${option}" in
  R) fasta=${OPTARG} ;;
  p) chrompattern=${OPTARG} ;;
  l) chromlength=${OPTARG} ;;
  h) HELP=yes ;;
  esac
done

if [ "${HELP}" == "yes" ]; then

  ## Documentation on how to run script
  echo
  echo
  echo "====================================================================================================="
  echo "====================================================================================================="
  echo "========="
  echo "=== SCRIPT FOR Collapsing contigs in reference genomes"
  echo "========="
  echo "====================================================================================================="
  echo "====================================================================================================="
  echo
  helppage
  echo "===AATGAAATTGAATCAAGTCCAAAAGTGAATAAAAGAAAAATCAAATGAGAAATAAAATAAAGAATAAAATTTCAAGTGTACTAAAATACGTACGT==="
  echo "===TTACTTTAACTTAGTTCAGGTTTTCACTTATTTTCTTTTTAGTTTACTCTTTATTTTATTTCTTATTTTAAAGTTCACATGATTTTATGCATGCA==="
  echo "====================================================================================================="
  exit 0
fi
module purge
module load picard
module load SAMtools

#Check if fasta file has been indexed

if [[ ! -e "${fasta}.fai" ]];
then
  samtools faidx "${fasta}"
fi


if [[ "${chrompattern}" == "none" && ${chromlength} -eq 0 ]];
then
    echo "Provide a value to -p or -l !"
    helppage
    exit $E_OPTERROR
fi

tempdir=$(mktemp -d)
dict=${fasta%.*}".dict"
#Check if dictionary file exists
#Check if reference dictionary exists, if not, build one
if [[ ! -e "${dict}" ]]
then
      echo "Building dictionary for $(basename $fasta)"
      java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R="${fasta}" O=${dict}
fi

chroms_contigs=""
otherchroms=""
if [[ "${chrompattern}" == "none" ]];
then
  chroms_contigs=$(grep -v '@HD' ${dict}|sed -e "s/LN://"|sed -e "s/SN://"|sort -n -k3|awk -v ll=${chromlength} '{if($3<ll)printf $2" "}')
  otherchroms=$(grep -v '@HD' ${dict}|sed -e "s/LN://"|sed -e "s/SN://"|awk -v ll=${chromlength} '{if($3>=ll)printf $2" "}')
  chrompattern="chr"
else
  chroms_contigs=$(grep -v '@HD' ${dict}|awk '{print $2}'|grep -vi "${chrompattern}"|sed -e "s/LN://"|sort -n|awk '{split($0,a,":");printf a[2]" "}')
  otherchroms=$(grep -v '@HD' ${dict}|awk '{print $2}'|grep -i "${chrompattern}"|awk '{split($0,a,":");printf a[2]" "}')
fi

# Get the number of `N` characters to generate.
n=100
# Create a string with `n` `N` characters.
separate_Ns="N"
for (( i=1; i<$n; i++ )); do
  separate_Ns+="N"
done

#Create temp file for writing
allchroms="${tempdir}/allchroms_temp.fasta"

echo "Copying standard ${chrompattern} sequences to ${allchroms} temp file"
samtools faidx $fasta $otherchroms > ${allchroms}
echo ">${chrompattern}Un" >> ${allchroms}
i=1

#samtools faidx $fasta ${chroms_contigs}|sed -e "s/>.*/${separate_Ns}/">> ${allchroms}
for chrom in ${chroms_contigs};
do
    if [[ $i -gt 1 ]];
    then
      echo "${separate_Ns}" >> ${allchroms}
    fi
    samtools faidx $fasta $chrom|grep -v '>'>> ${allchroms}
    echo "Collapsing contig ${i}. $chrom"
    i=$(($i+1))
done

#index temp file


outfile=$(echo ${fasta}|awk '{split($0,a,".");sub("."a[length(a)],"_with_unanchored_chrom");print $0"."a[length(a)]}')
rm -fR ${outfile}

echo "Normalise length of fasta file ${allchroms} to $outfile"

awk '
BEGIN {
  RS = ">"
  FS = "\n"
  OFS = ""
  left_trim = 1
  width = 69
}
NR > 1 {
  print RS $1
  $1 = ""
  seq = $0
  for (i = left_trim; i < length(seq); i += width) {
    print substr(seq, i, width)
  }
}
' ${allchroms} > ${outfile}

#Index new
samtools faidx ${outfile}

module purge
module load picard
java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R=${outfile}

#Normalise final fasta file

#cp ${allchroms} ${outfile}
#
#java -jar $EBROOTPICARD/picard.jar NormalizeFasta \
#      I=${allchroms_ref} \
#      LINE_LENGTH=70 \
#     O=${outfile}

#delete temp directory
rm -fR ${tempdir}