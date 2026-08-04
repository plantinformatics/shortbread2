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


// - Branch selection: convert "TRUE"/"FALSE" to boolean 
params.RunMPILEUP = (params.RunMPILEUP ?: 'FALSE')
def RUN_MPILEUP = ['TRUE','T','YES','Y','1'].contains(params.RunMPILEUP.toString().trim().toUpperCase())

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
  GENOME_METADATA;
  VCF_METADATA;
} from './nf/preprocessing.nf'

include {
  RUN_ALIGNMENT;
  MERGE_BAMS_BYSAMPLEID;
} from './nf/alignment.nf'

include {
  RUN_GENOTYPE_MPILEUP;
  RUN_MERGE_MPILEUP;
  GATHER_VCF_MPILEUP as GATHER_RAW_VCFs_mpileup;
  GATHER_VCF_MPILEUP as GATHER_FILTERED_VCFs_mpileup;
} from './nf/mpileup.nf'

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

// Git commit version, git url, added branch name, plus cast paramaters information for future incorporation into final vcfs
// Utilises First, Nextflow parameters and secondary failsafe using .git config file
// Note: I've tested this part and it seems robust but it's mostly ChatGPT, I'm not familiar enough with Nextflow and nextflow standard parameters to break down every component of it (JO)

def shortShaFromGit = { File projDir ->
    try {
        def head = new File(projDir, '.git/HEAD')
        if (!head.exists()) return null
        def txt = head.text.trim()
        if (txt.startsWith('ref:')) {
            def refPath = txt.split(':',2)[1].trim()
            def refFile = new File(projDir, ".git/${refPath}")
            return refFile.exists() ? refFile.text.trim().take(7) : null
        } else {
            return txt.take(7)
        }
    } catch (ignored) { null }
}

def headBranchFromGit = { File projDir ->
    try {
        def head = new File(projDir, '.git/HEAD')
        if (!head.exists()) return null
        def txt = head.text.trim()
        if (!txt.startsWith('ref:')) return null
        def refPath = txt.split(':',2)[1].trim()
        return refPath.tokenize('/').last()
    } catch (ignored) { null }
}

def readGitConfig = { File projDir ->
    def url = null
    def branches = []
    try {
        def cfg = new File(projDir, '.git/config')
        if (!cfg.exists()) return [url: null, branches: []]
        def curSection = ''
        cfg.eachLine { line ->
            def ln = line.trim()
            if (!ln) return
            def mSec = (ln =~ /^\[(.+?)\]\s*$/)
            if (mSec.matches()) {
                curSection = mSec[0][1]
                return
            }
            if (curSection.toLowerCase().startsWith('remote "origin"')) {
                def mUrl = (ln =~ /^url\s*=\s*(.+)$/)
                if (mUrl.matches()) url = mUrl[0][1].trim()
            }
            if (curSection.toLowerCase().startsWith('branch "')) {
                def mBr = (curSection =~ /^branch\s+"(.+)"$/)
                if (mBr.matches()) branches << mBr[0][1]
            }
        }
    } catch (ignored) { /* noop */ }
    [url: url, branches: branches]
}

// SAFE normalizer (no regex backrefs)
def normalizeRepoUrl = { String s ->
    if (!s) return ''
    String t = s.toString()
    if (t.startsWith('git@') && t.contains(':')) {
        int colon = t.indexOf(':')
        String host = t.substring('git@'.length(), colon)
        String path = t.substring(colon + 1)
        t = "https://${host}/${path}"
    } else if (t.startsWith('ssh://git@')) {
        String rest = t.substring('ssh://git@'.length())
        int slash = rest.indexOf('/')
        if (slash > 0) {
            String host = rest.substring(0, slash)
            String path = rest.substring(slash + 1)
            t = "https://${host}/${path}"
        }
    }
    if (t.endsWith('.git')) t = t.substring(0, t.length() - 4)
    return t
}

final File   _projDir   = new File( (workflow.projectDir ?: '.').toString() )
final String _repoRaw   = (workflow.repository ?: '').toString()
final String _commitSha = (workflow.commitId ?: shortShaFromGit(_projDir) ?: '').toString()
final String _revision  = (workflow.revision ?: '').toString()

final def    _cfg       = readGitConfig(_projDir)
final String _cfgUrl    = (_cfg.url ?: '')
final String _cfgBranch = headBranchFromGit(_projDir) ?: (_cfg.branches ? _cfg.branches[0] : '')

final String shortbread_repo_url = normalizeRepoUrl( _repoRaw ?: _cfgUrl )
final String shortbread_branch   = (_revision ?: _cfgBranch ?: '')

final String shortbread_version = _commitSha
    ? (shortbread_branch ? "${shortbread_branch}@${_commitSha.take(7)}" : "local@${_commitSha.take(7)}")
    : 'unknown@unknown'

log.info "[shortbread2] version: ${shortbread_version}"
log.info "[shortbread2] repo:    ${shortbread_repo_url ?: '(none)'}"
log.info "[shortbread2] branch:  ${shortbread_branch ?: '(none)'}"




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
                row -> row.SampleID
                }.unique().set{ samples }

        PREPARE_SAMPLE_SHEET.out.splitCsv(header: true)
        .map{ row ->
            tuple(
                row.SampleID,
                row.RunID,
                row.SampleName,
                row.Read1,
                row.Read2,
                row.SEQType,
                row.RGID,
                row.RGLB,
                row.RGPL,
                row.RGPU,
                row.RGSM
            )
        }.set{ alignments }
        if(!params.fastqc||params.runfastqconly){
            //Check sequencing quality
            RUN_FASTQC(
            PREPARE_SAMPLE_SHEET.out,
            params.runpath
            )
        }
    }



    GENOME_METADATA(
        params.accession
    )



    def variantcallmethod = RUN_MPILEUP ? "MPILEUP" : "GATK"


    VCF_METADATA(
        GENOME_METADATA.out,
        params.trimmethod,
        params.aligner,
        shortbread_version,
        variantcallmethod,
        workflow.start,
        (shortbread_repo_url   ?: ''),   
        (shortbread_branch ?: '')    
    )



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
                .map{ sampleid, bamlist, mergeflags ->
                    def mergeDecision = mergeflags.collect{ it.toString() }.contains('yes') ? 'yes' : 'no'
                    tuple(sampleid, bamlist.unique(), mergeDecision)
                }
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
            if(!RUN_MPILEUP)
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
                    "raw",
                    VCF_METADATA.out
                )
                GENERATE_RAW_VARLIST(GATHER_RAW_VCFs.out.map{it->it[1]})

                FILTER_VAR.out.groupTuple().set{filteredgenotypes}
                GATHER_FILTERED_VCFs(
                    filteredgenotypes,
                    params.refgenome,
                    "filtered",
                    VCF_METADATA.out
                )
                GENERATE_FILTERED_VARLIST(GATHER_FILTERED_VCFs.out.map{it->it[1]})
                GATHER_RAW_VCFs.out.collect()
                .combine(GATHER_FILTERED_VCFs.out.collect())
                .set{multiqcwait}
            }
        //End of GATK bits and bobs
            if(RUN_MPILEUP)
            {
                RUN_GENOTYPE_MPILEUP(
                    intervallist,
                    params.refgenome,
                    params.GATKHaplotypeoptions,
                    maxchromsize,
                    gatkreferencevcf
                )


                RUN_GENOTYPE_MPILEUP.out
                  .groupTuple(by: 0, size: params.numberofsamples)
                  .map { interval, chroms, files, fileindex, indexes ->
                      tuple(interval, chroms[0], files, fileindex, indexes[0])
                  }
                  .set { gvcfs }


                RUN_MERGE_MPILEUP(
                    gvcfs,
                    samples.collect()
                )

                RUN_MERGE_MPILEUP.out.groupTuple()
                .set{genotypes}
                FILTER_VAR(
                    params.CRbcf,
                    params.MAFbcf,
                    params.ACbcf,
                    params.MQbcf,
                    params.otherfilters,
                    RUN_MERGE_MPILEUP.out,
                    params.numberofsamples,
                    params.keepmultiallelicbcf,
                    gatkreferencevcf,
                    params.keepindelsbcf)
                GATHER_RAW_VCFs_mpileup(
                    genotypes,
                    params.refgenome,
                    "raw",
                    VCF_METADATA.out
                )
                GENERATE_RAW_VARLIST(GATHER_RAW_VCFs_mpileup.out.map{it->it[1]})

                FILTER_VAR.out.groupTuple().set{filteredgenotypes}
                GATHER_FILTERED_VCFs_mpileup(
                    filteredgenotypes,
                    params.refgenome,
                    "filtered",
                    VCF_METADATA.out
                )
                GENERATE_FILTERED_VARLIST(GATHER_FILTERED_VCFs_mpileup.out.map{it->it[1]})
                GATHER_RAW_VCFs_mpileup.out.collect()
                .combine(GATHER_FILTERED_VCFs_mpileup.out.collect())
                .set{multiqcwait}


            }
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

