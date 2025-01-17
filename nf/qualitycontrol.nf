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
    #Split reads into groups of 100
    cat ${samplesheet}|grep -vi 'Read1'|awk -F "," '{printf("%s\\n%s\\n",\$(NF-2),\$(NF-1))}'|split -l 100 - subset_reads
    for ss in subset_reads*;
    do
    {
        fastqc --outdir "\${logpath}" --threads $task.cpus \$(cat \${ss}|awk '{printf("%s ",\$0)}')
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

