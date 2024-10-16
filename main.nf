/*
 * Copyright (c) 2024
 */

/*
 * 'shortbread2' - Scaling PFG SNP calling capabilities with Nextflow
 *
 * This pipeline includes steps for Genotyping NGS data
 *
 */

/*
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

log.info """\

G'Day ${params.user}

Running shortbread2 version ${params.gittag}, git commit version [$params.gitver]

========================================================================================================================
$params.loginfo
========================================================================================================================
"""
System.out.println "Number of samples included in the analysis:"+String.format("%,d",params.numberofsamples)
/*
 * Include processes 
 */
include {
  PREPARE_SAMPLE_SHEET;
  PREPARE_GENOME;
  GET_GENOMEINFO;
  PREPARE_INTERVALS;
  SPLIT_INTERVALS;
} from './nf/preprocessing.nf'

include {
  RUN_ALIGNMENT;
  MERGE_BAMS_BYSAMPLEID;
} from './nf/alignment.nf'

include {
  READ_GATKDBs;
  RUN_GATK_HAPLOTYPE_CALLER;
  BUILD_GENOMICSDBImport;
  UPDATE_GENOMICSDB;
  RUN_GENOTYPEGVCFs;
} from './nf/GATK.nf'

include {
  FILTER_VAR;
  GATHER_VCF as GATHER_RAW_VCFs;
  GATHER_VCF as GATHER_FILTERED_VCFs;
  GENERATE_VariantList as GENERATE_RAW_VARLIST;
  GENERATE_VariantList as GENERATE_FILTERED_VARLIST;
} from './nf/postprocessing.nf'

include {
    RUN_FASTQC;
    GENERATE_VCFSTATS;
    MULTIQC;
} from './nf/qualitycontrol.nf'


/*
 * main pipeline logic
 */
println("Starting shortbread2 workflow")

workflow {
    // Step1 - Prepare data files
    if(!params.skipalignhapdb)
    {  
        //Build samplesheet
        PREPARE_SAMPLE_SHEET(
                params.rawdata,
                params.samplesheet,
                params.runpath,
                params.RGPL,
                params.mode
            )
        
        PREPARE_SAMPLE_SHEET.out.splitCsv(header: true).map{
                row ->row.SampleID
                }.unique().set{ samples }

        PREPARE_SAMPLE_SHEET.out.splitCsv(header: true)
        .map{ row ->
            tuple(row.SampleID,row.SampleName,row.Read1,row.Read2,row.SEQType)
        }.groupTuple().map{
                ID, samples,Read1,Read2,SEQType -> tuple( groupKey(ID, samples.size()), samples)
        }.transpose().set{ alignments }
        if(!params.fastqc||params.runfastqconly){
            //Check sequencing quality
            RUN_FASTQC(
            PREPARE_SAMPLE_SHEET.out,
            params.runpath
            )
        }
    }

    //Check whether to run fastqc only
    if(!params.runfastqconly)
    {
         GET_GENOMEINFO(
            params.refgenome,
            params.GATKreferenceVCF
            )

         GET_GENOMEINFO.out
        .set{chromsdata}

        //Prepare Reference genome details
        chroms=chromsdata.map{id->"${id[1]}".split(" ")}.flatten()
        maxchromsize=chromsdata.map{id->id[0]}
        genomesize=chromsdata.map{id->id[3]}
        gatkreferencevcf=chromsdata.map{id->id[5]}
        //Prepares index from the reference genome
        PREPARE_GENOME(
              params.refgenome,
              params.refannotation,
              params.aligner,
              params.runpath
        )
        if(!params.gatk){
            //Prepare intervals from the reference genome
           PREPARE_INTERVALS(
                params.refgenome,
                params.intervalsizegatk,
                GET_GENOMEINFO.out,
                params.mode,
                params.GATKpathtodbs,
                params.chromtoexclude
            )

            //Add group ID to intervals based on the number of samples
            PREPARE_INTERVALS.out.map{id->
                "${id}".split(" ")}
                .flatten()
                .map{interval->
                    return(groupKey(interval,params.numberofsamples))
                 }
                .set{intervals}

            //Split large intervals to smaller intervals before building DB
            SPLIT_INTERVALS(
                intervals,
                params.GATKpathtodbs,
                params.intervalsizegatk
            )
           SPLIT_INTERVALS.out.map{interval,subintervals->
                    subs=[]
                    subintervals.split( "=").each{subs.add(it)}
                    tuple(interval,subs)
                  }
                .transpose()
                .set{subintervals}
        }

       //Check whether to skip alignment, haplotyping step and database building step
       if(!params.skipalignhapdb)
       {    
            // Step2 - Alignment steps
            if(!params.alignment)
            {
                RUN_ALIGNMENT(
                    alignments,
                    params.trimmethod,
                    PREPARE_SAMPLE_SHEET.out,
                    params.runpath,
                    params.trimmeroptions,
                    params.OFFSET,
                    params.aligner,
                    PREPARE_GENOME.out,
                    params.seqtype,
                    params.aligneroptions,
                    params.keeptrimmedfqs,
                    params.trimming
                 )
                RUN_ALIGNMENT.out
                .groupTuple()
                .map{it->tuple(it.get(0),it.get(1).unique(),it.get(2)[0])}
                .set{mergealignments}

                // Call haplotypes for each sample and then build GATK db
                MERGE_BAMS_BYSAMPLEID(
                    mergealignments,
                    params.markduplicates,
                    params.runpath,
                    params.samtoolsoptions,
                    params.mappingquality,
                    params.aligner,
                    params.seqtype,
                    params.refgenome,
                    maxchromsize
                )
                bamqcwait=MERGE_BAMS_BYSAMPLEID.out
                if(!params.gatk)
                {
                    //Combine intervals with merged bam files
                    MERGE_BAMS_BYSAMPLEID.out.combine(intervals)
                    .map{
                        sampleid,bams,bam_index,intervals->
                        tuple(groupKey(intervals,sampleid.size()),sampleid,bams,bam_index)
                    }
                    .set{intervallist}
                }
            }
            else{
                //if bam files exist, read files and combine with interval list
                if(params.mode=="test")
                {
                    bamfiles=Channel.fromPath(params.bamdir+"/**.bam", checkIfExists:true).take(5)
                }
                else{
                    bamfiles=Channel.fromPath(params.bamdir+"/**.bam", checkIfExists:true)
                }

                PREPARE_INTERVALS.out.map{
                    id-> id.split(" ")
                    }
                    .flatten()
                    .combine(bamfiles)
                    .map{
                        intervals,bams ->tuple(intervals,bams.baseName.replaceAll("[_|.]sorted*",""),bams,bams+".csi")
                    }
                    .map{interval,sampleid,bam,bamindex->
                      tuple(groupKey(interval,sampleid.size()),sampleid,bam,bamindex)
                    }.set{intervallist}
                bamqcwait=intervallist
            }
       }

        //Check whether to skip all GATK steps
       if(!params.gatk)
        {
            //Check whether to skip haplotype and database step
            if(!params.skipalignhapdb)
            {
                // Step3 - GATK processing steps
                RUN_GATK_HAPLOTYPE_CALLER(
                    intervallist,
                    params.refgenome,
                    params.GATKHaplotypeoptions,
                    maxchromsize,
                    gatkreferencevcf
                    )
                subintervals.combine(RUN_GATK_HAPLOTYPE_CALLER.out
                .groupTuple(by:0,size:params.numberofsamples),by:0)
                .map{it->tuple(it.get(1),it.get(2),it.get(3),it.get(4)[0])}
                .set{gvcfs}
            }
            if(params.GATKupdateexistingdb)
            {
                UPDATE_GENOMICSDB(
                    gvcfs,
                    params.GATKDBImportoptions,
                    params.GATKupdateexistingdb,
                    params.GATKpathtodbs
                )
                UPDATE_GENOMICSDB.out.set{dbimport}
            }
            else if(params.skipalignhapdb) 
            {
                gatkdbs=Channel.fromPath(params.GATKpathtodbs+"/*/callset.json")
                READ_GATKDBs(gatkdbs)
                subintervals.combine(READ_GATKDBs.out,by:0)
                .map{it->
                    tuple(it[2],it[1],it[3],it[4])}
                .set{dbimport}
            }
            else
            {
                BUILD_GENOMICSDBImport(
                gvcfs,
                params.GATKDBImportoptions
                )
                BUILD_GENOMICSDBImport.out.set{dbimport}
            }
                   
            RUN_GENOTYPEGVCFs(
                dbimport,
                params.GATKGenotypingoptions,
                params.refgenome,
                gatkreferencevcf
            )

            //Group intervals by Chromosome in readiness for the gather step
            RUN_GENOTYPEGVCFs.out.groupTuple()
            .set{genotypes}
            FILTER_VAR(
                params.CRbcf,
                params.MAFbcf,
                params.ACbcf,
                params.MQbcf,
                params.otherfilters,
                RUN_GENOTYPEGVCFs.out,
                params.numberofsamples,
                params.keepmultiallelicbcf,
                gatkreferencevcf,
                params.keepindelsbcf)
            GATHER_RAW_VCFs(
                genotypes,
                params.refgenome,
                "raw"
            )
            GENERATE_RAW_VARLIST(GATHER_RAW_VCFs.out.map{it->it[1]})

            FILTER_VAR.out.groupTuple().set{filteredgenotypes}
            GATHER_FILTERED_VCFs(
                filteredgenotypes,
                params.refgenome,
                "filtered"
            )
            GENERATE_FILTERED_VARLIST(GATHER_FILTERED_VCFs.out.map{it->it[1]})
            GATHER_RAW_VCFs.out.collect()
            .combine(GATHER_FILTERED_VCFs.out.collect())
            .set{multiqcwait}
        }
        else{
            bamqcwait.collect().set{multiqcwait}
        }
        if(!params.skipalignhapdb&&!params.fastqc)
         multiqcwait.combine(RUN_FASTQC.out).set{multiqcwait}
    }
    else
    {
        RUN_FASTQC.out.set{multiqcwait}
    }
   MULTIQC(
        multiqcwait,
        "${params.runpath}/Logs",
        params.multiqcconfig
    )
}

