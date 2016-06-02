## http://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/
## http://www.homolog.us/blogs/blog/2013/11/22/cshl-keynote-talk-lior-pachter/
## https://www.biostars.org/p/133488/


##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
pathToCountFiles=args[1]
target=args[2]
trans_len=args[3]
gene_transcript_map=args[4]
#print(pathToCountFiles)
print(target)

library("edgeR")

files <- dir(path=pathToCountFiles,pattern=paste("^", target, ".*quant.counts$", sep = ""))
print(files)
data <- readDGE(files)
head(data$counts)
dataCounts_sorted=as.matrix(data$counts[order(rownames(data$counts)), ])
head(dataCounts_sorted)
rawCounts=rowSums(dataCounts_sorted)  ## raw counts mapping to transcripts

data2=read.table(trans_len, header=T)
data2_sorted=data2[with(data2, order(transcript)), ]
head(data2_sorted)
length=data2_sorted$length
stopifnot(all(rownames(dataCounts_sorted)==data2_sorted$transcript))
dataSummary=cbind(length,rawCounts)
head(dataSummary)

normCount=dataSummary[,2]/(dataSummary[,1]/1000) ## Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK)
normCountSum=sum(normCount) ## Count up all the RPK values in a sample
calcTPM=(dataSummary[,2] * 10^6) / ((dataSummary[,1]/1000)*normCountSum) ## Divide the RPK values by total RPK per 10^6 ## This formula can be written as: calcTPM=normCount/(normCountSum/10^6)
dataSummary2=cbind(dataSummary,calcTPM)
head(dataSummary2)

data3=read.table(gene_transcript_map)
data3_sorted= data3[with(data3, order(V2)), ]
genes=data3_sorted$V1
dataSummary2_sorted <- dataSummary2[ order(row.names(dataSummary2)), ]
stopifnot(all(rownames(dataSummary2_sorted)==data3_sorted$V2))
dataSummary3=data.frame(genes,dataSummary2_sorted)      ## Add gene names
head(dataSummary3)
#write.table(dataSummary3, file='dataSummary', sep='\t', quote=F, row.names=T, col.names=NA)

transcripts=rownames(dataSummary3)
dataSummary4=data.frame(transcripts,dataSummary3)

## calculate gene statistics
sumTPM=aggregate(dataSummary3$calcTPM, by=list(Category=dataSummary3$genes), FUN=sum)
colnames(sumTPM)=c("genes","sumTPM")
sumCounts=aggregate(dataSummary3$rawCounts, by=list(Category=dataSummary3$genes), FUN=sum)
colnames(sumCounts)=c("genes","sumCounts")
maxLength=aggregate(dataSummary3$length, by=list(Category=dataSummary3$genes), FUN=max)
colnames(maxLength)=c("genes","maxLength")

dataSummary5=merge(dataSummary4,sumTPM,by.x="genes",by.y="genes",sort=F)
dataSummary6=merge(dataSummary5,sumCounts,by.x="genes",by.y="genes",sort=F)
dataSummary7=merge(dataSummary6,maxLength,by.x="genes",by.y="genes",sort=F)

## calculate the TPM ratio per gene for each transcript
isoformTPM=(dataSummary7$calcTPM/dataSummary7$sumTPM)*100
isoformTPM[is.nan(isoformTPM)] <- 0
isoformTPM=format(round(isoformTPM, 3), nsmall = 3)
dataSummary8=data.frame(dataSummary7,isoformTPM)
#head(dataSummary8)

if (target == "") { 
   write.table(dataSummary8, file='dataSummary_comp', sep='\t', quote=F, row.names=T, col.names=NA)
} else { 
  write.table(dataSummary8, file=paste(target, 'dataSummary_comp', sep = "."), sep='\t', quote=F, row.names=T, col.names=NA)
}
################################
pdf("rawCount-hist-plot.pdf")
hist(rawCounts,xlim=c(0, 100), breaks = 10000000)
dev.off()

pdf("TPM-hist-plot.pdf")
hist(calcTPM,xlim=c(0, 20), breaks = 100000)
dev.off()

pdf("isoformTPM-hist-plot.pdf")
hist(as.numeric(isoformTPM),xlim=c(0, 100), breaks = 100)
dev.off()

pdf("isoformTPM-hist-plot-focus.pdf")
hist(as.numeric(isoformTPM),xlim=c(0, 10), breaks = 100)
dev.off()

pdf("isoformTPM-hist-plot-focus2.pdf")
hist(as.numeric(isoformTPM),xlim=c(0, 5), breaks = 1000)
dev.off()

#reducedData=extCounts[totCounts>20,]
#write.table(reducedData, file='reducedData', sep='\t', quote=F, row.names=T, col.names=T)

