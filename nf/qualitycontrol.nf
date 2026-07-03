/*
Process to run fastqc on the raw read files
*/
def output=params.runpath
process RUN_FASTQC {
    tag "${samplesheet.baseName}"
    label 'usescratch'
    label "${params.fastqclabel}"
    label 'qualitycontrol'
    errorStrategy 'finish'
    input:
        path samplesheet
        path results
    output:
        env logpath
    script:
    multiqc="${output}/Reports/"
    """
    #!/bin/bash
    logpath="${output}/Logs/01_FastQC/"
    mkdir -p \${logpath}
    mkdir -p "${multiqc}"
    # Build a robust list of FASTQ files from the normalized samplesheet.
    # Read1 and Read2 remain the last two read-path columns in the revised sheet:
    # SampleID,RunID,SampleName,FCID,Lane,SampleIndex,RGID,RGLB,RGPL,RGPU,RGSM,Read1,Read2,SEQType
    tail -n +2 ${samplesheet} \\
    | awk -F "," '{print \$(NF-2); if(\$(NF-1)!="") print \$(NF-1)}' \\
    | awk 'NF>0' > fastqc_reads.list

    if [[ ! -s fastqc_reads.list ]]; then
        echo "No FASTQ files found in normalized samplesheet for FastQC"
        exit 0
    fi

    # Split reads into groups of 100 to avoid overly long command lines
    split -l 100 fastqc_reads.list subset_reads_

    for ss in subset_reads_*;
    do
    {
        [[ ! -s "\${ss}" ]] && continue
        fastqc --outdir "\${logpath}" --threads $task.cpus \$(awk '{printf("%s ",\$0)}' "\${ss}")
    }
    done
    """
}


/*
Process to perform post processing of the log files
*/
process MULTIQC()
{
    tag "multiqc"
    label 'usescratch'
    errorStrategy 'finish'
    input:
        val genotypefilt
        path multiqc_files
        path multiqc_configfile
    script:
    multiqc="${output}/Reports/"
    """
    #!/bin/bash
    echo "Running multiqc step"
    multiqc -f --config ${multiqc_configfile} \\
            -n "Shortbread2-quality-report" \\
            --outdir "${multiqc}" "${multiqc_files}"
    """
}

process GENERATE_VCFSTATS {
    tag "${vcf.baseName}"
    label 'medium'
    label 'qualitycontrol'
    input:
        path vcf
    script:
    """
    #!/bin/bash
    hets=\$(bcftools view -H -i 'GT="0/1"' ${vcf}|wc -l)
    homs=\$(bcftools view -H -e 'GT="0/1"' ${vcf}|wc -l)
    multialleles=\$(bcftools view -m2 -v snps ${vcf}|wc -l)
    bcftools query -f "%CHROM:%POS\\t%INFO/ExcessHet\\t%INFO/MAF\\t%INFO/CR\\t%INFO/AF\\n" ${vcf} \\
    |awk 'BEGIN{print "POS\\tExcessHet\\tMAF\\tCR\\tAF"} {print \$0}' > /tmp/Stats.data
    Rscript -e "data <- read.table('/tmp/Stats.data', header = T,row.names=1);
                pdf('\$(pwd)/STATs-plot-RNAseq.pdf',height=10,width=10);
                par(mfrow=c(2,2));
                for(col in colnames(data))
                    hist(data[,col],main=col,xlab=col)
                dev.off();

                pdf('\$(pwd)/STATs-ExcessHets-vs-other-scatter.pdf',height=6,width=12);
                par(mfrow=c(1,2));
                for(stat in c('MAF','CR'))
                    plot(data[,'ExcessHet'],data[,stat],main=stat,pch=16,cex=0.5,xlab='ExcessHet',ylab=stat,cex.main=1);
                dev.off();"
    """
}

