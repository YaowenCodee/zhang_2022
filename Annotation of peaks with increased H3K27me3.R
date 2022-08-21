rm(list=ls())
options(stringsAsFactors=F)
library(ChIPseeker)
library(org.Dm.eg.db)
library(clusterProfiler)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene

#read narrowPeak files and insert title
narrowPeak<-dir(path="./",pattern="broadPeak$")
for(file in narrowPeak){
  
  nPeak<-read.table(paste("./",file,sep=""),head=F,)
  names(nPeak)<-c("chrom","chromStart","chromEnd")
                  #(,"name","score","strand","signalValue","pValue","qValue")
  write.table(nPeak,paste("./",file,".insert.title",sep=""),
              row.names = F,col.names = T,quote = F,sep="\t")
}

#read peak file
rkr<-readPeakFile("./H3K27me3_WTvs_401_DiffPeaks_MAvalues_1.broadPeak.insert.title")
#rwr<-readPeakFile("./H3K27me3_WT_vs_IP_401_all_MAvalues_Diffpeaks.broadPeak.insert.title")
#Chip peaks coverage plot, whole genome
covplot(rkr)
#Chip peaks coverage plot, specific chr
covplot(rkr,chrs=paste("chr",c('2L','2R','3L','3R',"X","Y"),sep=""))

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

tagMatrix <- getTagMatrix(rkr, windows=promoter)
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
peakHeatmap(rkr, TxDb=txdb, upstream=3000, downstream=3000, color="red")
plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
            xlab = "Genomic Region (5' -> 3')",
            ylab = "Read Count Frequency")
    
plotAvgProf(tagMatrix, xlim = c(-3000,3000),
            conf = 0.95, resample= 1000)
######
size<-5000
peakAnnoRKR <- annotatePeak(rkr, tssRegion=c(-size, size),TxDb=txdb, annoDb="org.Dm.eg.db",verbose=FALSE)
#peakAnnoRWR <- annotatePeak(rkr, tssRegion=c(-size, size),TxDb=txdb, annoDb="org.Dm.eg.db",verbose=FALSE)
plotAnnoBar(peakAnnoRKR)
plotAnnoPie(peakAnnoRKR)
#vennpie
vennpie(peakAnnoRKR)
#upsetplot(peakAnnoRKR, vennpie=TRUE)
#upsetplot(peakAnnoRKR) 
plotDistToTSS(peakAnnoRKR,title="Distribution of transcription factor-binding loci\nrelative to TSS")

#####
peakAnnoRKR <- as.data.frame(peakAnnoRKR) 
write.csv(peakAnnoRKR,file = "anno2.csv") 
