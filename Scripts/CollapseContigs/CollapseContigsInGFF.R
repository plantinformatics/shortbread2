updateContigsGFF <- function (gff.file,fasta.dict,chrUn.name="chrUn")
{
  #check if gff format
  setwd("/group/grains/shortbread-work/test-run/")
  fasta_dict<-read.table(fasta.dict,sep="\t",header = F,skip = 1)
  fasta_dict <- fasta_dict[,-1]
  fasta_dict <- as.data.frame(apply(fasta_dict,2,function(x) gsub(".*:","",x)))
  fasta_dict$V3<-as.numeric(fasta_dict$V3)
  fasta_dict<-fasta_dict[!grepl("^chr",fasta_dict$V2,ignore.case = T),]
  #Sort data frame
  fasta_dict<-fasta_dict[order(fasta_dict$V3,decreasing = F),]
  starts<-cumsum(fasta_dict$V3)+101
  fasta_dict$Start<-c(1,starts[-length(starts)])
  gff.data<-read.table(gff.file,header = F,sep="\t")
  fasta_starts<-with(fasta_dict,split(Start,V2))
  contigs<-intersect(names(fasta_starts),gff.data$V1)
  match.contigs<-match(contigs,gff.data$V1)
  gff.data$V4[match.contigs]<-gff.data$V4[match.contigs]+as.numeric(fasta_starts[contigs])
  gff.data$V5[match.contigs]<-gff.data$V5[match.contigs]+as.numeric(fasta_starts[contigs])
  gff.data$V1[match.contigs]<-"chrUn"
  head(gff.data[grepl("chrUn",gff.data$V1),])
  ext<-strsplit(basename(gff.file),"\\.")[[1]]
  ext<-ext[length(ext)]
  output<-gsub(paste0(".",ext),paste0("_with_",chrUn.name,".",ext),gff.file)
  write.table(gff.data,output,quote = F,row.names = F,col.names = F)
}

