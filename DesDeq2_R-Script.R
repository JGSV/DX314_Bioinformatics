#load required libraries
library(DESeq2)
library(stats)
library(biomaRt)
library(dplyr)
library(gtools)
library(apeglm)

#Number of CPU threads to use
register(SnowParam(3))

#Working directory containing count table, treatment ID file, and output location
getwd()
#setwd("C:/")

#Define reference group
refGroup<-"Control"
#Define treatment groups
toComp<-c("1nM_atRA","100nM_atRA","1uM_NMP314","1uM_NMP314.1nM_atRA","1uM_NMP308","1uM_NMP308.1nM_atRA")

#Import Data and Organize
colData_fle<-read.csv("SampleTreatmentID.csv") #two columns. column 1 contains sample names from count table, column two contains repective treatment group defined above. First row is header
names(colData_fle[,2])<-c("Treatment")
colData<-data.frame(colData_fle[colData_fle$Treatment %in% c(refGroup,toComp),2])
rownames(colData)<-colData_fle[,1]
names(colData)<-c("Treatment")

#Imports count table
countData <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))

#remove only rows with no counts
countData<-countData[rowSums(countData[,-1]) >0,]

#creates DESeq dataset
dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData,tidy=F,design=~Treatment)

#defines and sorts factor levels
con <- as.character(colData(dds)$Treatment)
colData(dds)$Treatment <- factor(con, levels = c("Control","1nM_atRA","100nM_atRA","1uM_NMP314","1uM_NMP314.1nM_atRA","1uM_NMP308","1uM_NMP308.1nM_atRA"))
dds$Treatment<-relevel(dds$Treatment,refGroup)

#Run DESeq
dds<-DESeq(dds,parallel=T)

#import BioMART annotations and creates reference index
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
#attr<-listAttributes(ensembl)
bm <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "description","goslim_goa_description"),
            filters="ensembl_gene_id",
            values=rownames(dds),
            mart=ensembl)
bm <- arrange(bm, ensembl_gene_id)
idx <- match(rownames(dds), bm$ensembl_gene_id)

#Function to process each comparison and annotate
process<-function(toComp,refGroup=refGroup,dds=dds,bm=bm,idx=idx,c=c){
  res <- lfcShrink(dds, coef = c+2,type="apeglm",parallel=T)
  res$hgnc_symbol <- bm$hgnc_symbol[idx]
  res$description <- bm$description[idx]
  res$GO <-bm$goslim_goa_description[idx]
  res$FoldChange <- logratio2foldchange(res$log2FoldChange, base=2)
  idx2 <- c(grep(paste0("^",refGroup,"$"),colData(dds)$Treatment),grep(paste0("^",toComp,"$"),colData(dds)$Treatment))
  res$NormCounts <- counts(dds[,idx2],normalized=T)
  colnames(res$NormCounts) <- colData(dds[,idx2])$Treatment
  loc<-file.path(getwd(),paste(refGroup,"_vs_",toComp,sep=""))
  dir.create(loc)
  write.csv(res, file=paste(loc,"/",refGroup,"_vs_",toComp,"_analysis.csv",sep=""))
  write.csv(res, file=paste(getwd(),"/Combined/Raw/",c,"_",refGroup,"_vs_",toComp,"_analysis.csv",sep=""))
}

#Run analysis for each comparision
c<-0
for(i in toComp){
  process(i,refGroup,dds,bm,idx,c)
  c<-c+1
}
