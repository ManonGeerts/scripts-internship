#-----------------------------------------------------------------------------------------------------------
#### FEW CHECKS
#-----------------------------------------------------------------------------------------------------------


library(ggplot2)
library(reshape2)
library(plyr)
library(ggpubr)


#ReadCount <- read.table("map.stat1.txt")
#colnames(ReadCount) <- c("sample","mapped reads","unmapped reads","total number of reads")
#ReadCountSmall <- data.frame(sample = ReadCount$sample, mapped = ReadCount$`mapped reads`, unmapped = ReadCount$`unmapped reads`)

#MeltedReadCount = melt(ReadCountSmall, id=c('sample'))
#names(MeltedReadCount) <- c('sample', 'mapping', 'reads')


#-----------------------------------------------------------------------------------------------------------
#### MAPPING STATS: check how many reads mapped against reference genome
#-----------------------------------------------------------------------------------------------------------


setwd("C:/Users/mgeerts/Dropbox/Gambiense/Foefelare Manon/mapping")
ReadCount <- read.table("map.stat3.txt")
names(ReadCount) <- c('sample', 'mapping', 'reads')

# var = unmapped, mapped with a quality > 25, mappes with a quality < 25
numberofvar=3
numberofsamples=42
nrow(ReadCount) == numberofvar*numberofsamples

ReadsFraction <- ddply(ReadCount,.(sample),summarise,Count.Fraction = reads / sum(reads))
to_graph <- cbind(arrange(ReadCount, sample), fraction = ReadsFraction$Count.Fraction)

gp <- ggplot(data=to_graph, aes(x=sample, y=reads, fill=mapping)) +
  geom_bar(stat="identity",position = "stack",alpha=0.5) +
  geom_text(aes(label=paste(round(fraction*100),"%", sep="")),size = 3,vjust=1,position="stack",color="white") +
  theme(axis.text.x=element_text(face="bold",angle = 90))
gp

##### only ILLUMINA reads

ILLUMINA=c("1829_Alijo","Fontem_S10", "KP33_clone_16", "PTAG_129","TSW_187_78E")
ReadCount2 <- ReadCount[ReadCount$sample %in% ILLUMINA,]

ReadsFraction2 <- ddply(ReadCount2,.(sample),summarise,Count.Fraction = reads / sum(reads))
to_graph2 <- cbind(arrange(ReadCount2, sample), fraction = ReadsFraction2$Count.Fraction)

gp2 <- ggplot(data=to_graph2, aes(x=sample, y=reads, fill=mapping)) +
  geom_bar(stat="identity",position = "stack",alpha=0.5) +
  geom_text(aes(label=paste(round(fraction*100),"%", sep="")),size = 2,vjust=1,position="stack",color="white") +
  theme(axis.text.x=element_text(face="bold",angle=90,size=5))
gp2


##### only BGI reads

ReadCount3 <- ReadCount[!(ReadCount$sample %in% ILLUMINA),]

ReadsFraction3 <- ddply(ReadCount3,.(sample),summarise,Count.Fraction = reads / sum(reads))
to_graph3 <- cbind(arrange(ReadCount3, sample), fraction = ReadsFraction3$Count.Fraction)

gp3 <- ggplot(data=to_graph3, aes(x=sample, y=reads, fill=mapping)) +
  geom_bar(stat="identity",position = "stack",alpha=0.5) +
  geom_text(aes(label=paste(round(fraction*100),"%", sep="")),size = 2,vjust=1,position="stack",color="white") +
  theme(axis.text.x=element_text(face="bold",angle=90,size=5))
gp3

figure <- ggarrange(gp2, gp3, labels = c("Illumina?", "BGI"), ncol = 1, nrow = 2)
figure

#-----------------------------------------------------------------------------------------------------------
#### TRIMMING STATS ON UNMAPPED READS: check how many reads passed trimming (on all samples)
# ! Illumina reads are shorter than BGI reads --> other trimming requirements
#-----------------------------------------------------------------------------------------------------------


setwd("C:/Users/mgeerts/Dropbox/Gambiense/Foefelare Manon/komics")

###### BGI reads: -q 30 -l 100 -b 150

ReadCount <- read.table("trim.stat.BGI.txt")
names(ReadCount) <- c('sample', 'trimming', 'reads')

# var = passed trimming, not passed
numberofvar=2
numberofsamples=37
nrow(ReadCount) == numberofvar*numberofsamples

ReadsFraction <- ddply(ReadCount,.(sample),summarise,Count.Fraction = reads / sum(reads))
to_graph <- cbind(arrange(ReadCount, sample), fraction = ReadsFraction$Count.Fraction)

gp <- ggplot(data=to_graph, aes(x=sample, y=reads, fill=trimming)) +
  geom_bar(stat="identity",position = "stack",alpha=0.5) +
  geom_text(aes(label=paste(round(fraction*100),"%", sep="")),size = 2,vjust=1,position="stack",color="white") +
  theme(axis.text.x=element_text(face="bold",angle=90,size=5))
gp

###### Illumina reads: -q 30 -l 75 -b 125

ReadCount2 <- read.table("trim.stat.ILLUMINA.txt")
names(ReadCount2) <- c('sample', 'trimming', 'reads')

# var = passed trimming, not passed
numberofvar=2
numberofsamples=5
nrow(ReadCount2) == numberofvar*numberofsamples

ReadsFraction2 <- ddply(ReadCount2,.(sample),summarise,Count.Fraction = reads / sum(reads))
to_graph2 <- cbind(arrange(ReadCount2, sample), fraction = ReadsFraction2$Count.Fraction)

gp2 <- ggplot(data=to_graph2, aes(x=sample, y=reads, fill=trimming)) +
  geom_bar(stat="identity",position = "stack",alpha=0.5) +
  geom_text(aes(label=paste(round(fraction*100),"%", sep="")),size = 2,vjust=1,position="stack",color="white") +
  theme(axis.text.x=element_text(face="bold",angle=90,size=5))
gp2

figure <- ggarrange(gp2, gp, labels = c("Illumina?", "BGI"), ncol = 1, nrow = 2)
figure

###### all reads

ReadCount3 <- rbind(ReadCount,ReadCount2)
names(ReadCount3) <- c('sample', 'trimming', 'reads')

# var = passed trimming, not passed
numberofvar=2
numberofsamples=42
nrow(ReadCount3) == numberofvar*numberofsamples

ReadsFraction3 <- ddply(ReadCount3,.(sample),summarise,Count.Fraction = reads / sum(reads))
to_graph3 <- cbind(arrange(ReadCount3, sample), fraction = ReadsFraction3$Count.Fraction)

gp3 <- ggplot(data=to_graph3, aes(x=sample, y=reads, fill=trimming)) +
  geom_bar(stat="identity",position = "stack") +
  geom_text(aes(label=paste(round(fraction*100),"%", sep="")),size = 3,vjust=1,position="stack") +
  theme(axis.text.x=element_text(face="bold", color="#993333",angle = 90))
gp3


#--------------------------------------------------------------------------------------------------------------
### TRIMMING CHECK: check how many reads are lost after SamToFastq and trimming
# with different trimming requirements on a few samples
# samplename.MQ.MINLEN.MAXLEN
#--------------------------------------------------------------------------------------------------------------


setwd("C:/Users/mgeerts/Dropbox/Gambiense/Foefelare Manon/komics/trimming optim")

ReadCount <- read.table("trimming.optim.txt")
names(ReadCount) <- c('sample', 'loss', 'reads')

# var = passed trimming, not passed
numberofvar=3
numberofsamples=6
nrow(ReadCount) == numberofvar*numberofsamples

ReadsFraction <- ddply(ReadCount,.(sample),summarise,Count.Fraction = reads / sum(reads))
to_graph <- cbind(arrange(ReadCount, sample), fraction = ReadsFraction$Count.Fraction)

gp <- ggplot(data=to_graph, aes(x=sample, y=reads, fill=loss)) +
  geom_bar(stat="identity",position = "stack") +
  geom_text(aes(label=paste(round(fraction*100),"%", sep="")),size = 3,vjust=1,position="stack") +
  theme(axis.text.x=element_text(face="bold", color="#993333",angle = 90))
gp

### 108AT: 24% loss of reads after SamToFastq
### 45% after trimming with MQ=30, MINLEN=100, MAXLEN=150.
### MQ reduced to 27 --> 7% more reads


#--------------------------------------------------------------------------------------------------------------
### CSB3: check how many reads have the CSB3 region (minicircles)
# unmapped reads vs. mapped reads Q<60/Q>60 (all samples)
# MR = mapped reads; UNMR = unmapped reads 
#--------------------------------------------------------------------------------------------------------------


setwd("C:/Users/mgeerts/Dropbox/Gambiense/Foefelare Manon/mapping/")

ReadCount <- read.table("CSB.count.all.txt")
names(ReadCount) <- c('sample', 'loss', 'reads')

numberofvar=6
numberofsamples=42
nrow(ReadCount) == numberofvar*numberofsamples

ReadsFraction <- ddply(ReadCount,.(sample),summarise,Count.Fraction = reads / sum(reads))
to_graph <- cbind(arrange(ReadCount, sample), fraction = ReadsFraction$Count.Fraction)

gp <- ggplot(data=to_graph, aes(x=sample, y=reads, fill=loss)) +
  geom_bar(stat="identity",position = "stack", alpha=0.5) +
  theme(axis.text.x=element_text(face="bold", angle=90))
  #geom_text(aes(label=paste(round(fraction*100),"%", sep="")),size = 1,vjust=1,position="stack")
gp

###### compare reads havings CSB3

ReadCount2 <- ReadCount[ReadCount$loss!="UNMR_without_CSB3",]
ReadCount2 <- ReadCount2[ReadCount2$loss!="MR_without_CSB3",]

numberofvar=4
numberofsamples=42
nrow(ReadCount2) == numberofvar*numberofsamples

ReadsFraction2 <- ddply(ReadCount2,.(sample),summarise,Count.Fraction = reads / sum(reads))
to_graph2 <- cbind(arrange(ReadCount2, sample), fraction = ReadsFraction2$Count.Fraction)

gp2 <- ggplot(data=to_graph2, aes(x=sample, y=reads, fill=loss)) +
  geom_bar(stat="identity",position = "stack", alpha=0.5) +
  geom_text(aes(label=paste(round(fraction*100),"%", sep="")),size = 2,vjust=1,position="stack",color="white") +
  theme(axis.text.x=element_text(face="bold",angle=90, size=6))
gp2

## on subset 
subset <- c("108AT", "Nabe","Fontem_S10","KP33_clone_16")

ReadCount4 <- ReadCount2[(ReadCount2$sample %in% subset),]

ReadsFraction4 <- ddply(ReadCount4,.(sample),summarise,Count.Fraction = reads / sum(reads))
to_graph4 <- cbind(arrange(ReadCount4, sample), fraction = ReadsFraction4$Count.Fraction)

gp4 <- ggplot(data=to_graph4, aes(x=sample, y=reads, fill=loss)) +
  geom_bar(stat="identity",position = "stack", alpha=0.5) +
  geom_text(aes(label=paste(round(fraction*100),"%", sep="")),size = 2,vjust=1,position="stack",color="white") +
  theme(axis.text.x=element_text(face="bold"))
gp4

figure <- ggarrange(gp2, gp4, labels = c("Reads without CSB3", "on subset"), ncol = 1, nrow = 2)
figure

###### compare reads without CSB3

ReadCount3 <- rbind(ReadCount[ReadCount$loss=="UNMR_without_CSB3",],ReadCount[ReadCount$loss=="MR_without_CSB3",])

ReadsFraction3 <- ddply(ReadCount3,.(sample),summarise,Count.Fraction = reads / sum(reads))
to_graph3 <- cbind(arrange(ReadCount3, sample), fraction = ReadsFraction3$Count.Fraction)

gp3 <- ggplot(data=to_graph3, aes(x=sample, y=reads, fill=loss)) +
  geom_bar(stat="identity",position = "stack",alpha=0.5) +
  theme(axis.text.x=element_text(face="bold",angle=90))
  #geom_text(aes(label=paste(round(fraction*100),"%", sep="")),size = 3,vjust=1,position="stack") +
gp3


#--------------------------------------------------------------------------------------------------------------
### CSB3: check how many reads with CSB3 region are left after SamTo FastQ and trimming
# on a few samples
#--------------------------------------------------------------------------------------------------------------


setwd("C:/Users/mgeerts/Dropbox/Gambiense/Foefelare Manon/komics/trimming optim")

ReadCount <- read.table("CSB.count.aftershit.txt")
names(ReadCount) <- c('sample', 'loss', 'reads')

numberofvar=3
numberofsamples=4
nrow(ReadCount) == numberofvar*numberofsamples

ReadsFraction <- ddply(ReadCount,.(sample),summarise,Count.Fraction = reads / sum(reads))
to_graph <- cbind(arrange(ReadCount, sample), fraction = ReadsFraction$Count.Fraction)

gp <- ggplot(data=to_graph, aes(x=sample, y=reads, fill=loss)) +
  geom_bar(stat="identity",position = "stack") +
  geom_text(aes(label=paste(round(fraction*100),"%", sep="")),size = 3,vjust=1,position="stack") +
  theme(axis.text.x=element_text(face="bold", color="#993333"))
gp


#--------------------------------------------------------------------------------------------------------------
### CSB3: coverage reads with CBS3 region
# on a few samples
#--------------------------------------------------------------------------------------------------------------

setwd("C:/Users/mgeerts/Dropbox/Gambiense/Foefelare Manon/komics/trimming optim")

#samples <- c("108AT","Nabe")
#for (s in 1:length(samples)) {
#  table[s] <- read.table(paste(samples[s],".MR_CBS3.freq-chrom.txt", sep=""))
#  }


AT108 <- read.table("108AT.MR_CSB3.freq-chrom.txt")
Nabe <- read.table("Nabe.MR_CSB3.freq-chrom.txt")
KP33_clone_16 <- read.table("KP33_clone_16.MR_CSB3.freq-chrom.txt")
Fontem_S10 <- read.table("Fontem_S10.MR_CSB3.freq-chrom.txt")

allmerged <- Reduce(function(x, y) merge(x, y, by="V1", all=TRUE), list(AT108, Nabe, KP33_clone_16,Fontem_S10))
#rownames(merged) <- merged[,1]
#merged <- merged[,-1]
colnames(allmerged) <- c("chrom","108AT","Nabe","KP33_clone_16","Fontem_S10")
allmerged[is.na(allmerged)] <- 0
melted = melt(allmerged, id="chrom")
colnames(melted) <- c("chrom","sample","reads")

gp <- ggplot(data=melted, aes(x=chrom, y=reads, fill=sample)) +
  geom_bar(stat="identity",position = "stack",alpha=0.5) +
  theme(axis.text.x=element_text(face="bold", color="#993333", angle=90))
  #geom_text(aes(label=value), vjust=1, color="white", size=2) +
gp

melted2 <- melted[melted$reads>2,]
gp2 <- ggplot(data=melted2, aes(x=chrom, y=reads, fill=sample)) +
  geom_bar(stat="identity",position = "stack", alpha=0.5) +
  #geom_text(aes(label=reads), vjust=1, color="white", size=2) +
  theme(axis.text.x=element_text(face="bold", size=6,angle=90))
gp2


#--------------------------------------------------------------------------------------------------------------
### Minicircles: extracted circ CSB3 contigs, different kmers
# 108AT/Nabe mapped vs unmapped reads
# UNMR.R1: kmin=79, kmax=139, kstep=10
# UNMR.R2: kmin=49, kmax=139, kstep=10
# UNMR.R3: kmin=49, kmax=69, kstep=10
# UNMR.R4: kmin=99, kmax=119, kstep=10
# UNMR.R5: kmin=89, kmax=129, kstep=10
# ALLR R1: kmin=79, kmax=139, kstep=10
# ALLR R2: kmin=99, kmax=119, kstep=10
#--------------------------------------------------------------------------------------------------------------

####### CHECK KOMICS DIFFERENT SETTINGS


### DATA PREP

setwd("C:/Users/mgeerts/Dropbox/Gambiense/Foefelare Manon/komics/trimming optim")

circ <- read.table(file = "Nabe.length.circ.minicircles.diffR.txt", sep="_")
circ$V3 <- as.numeric(gsub("len", "", circ$V3))
circ$V1 <- gsub(">","",circ$V1)
circ <- circ[,c(1,3)]
colnames(circ) <- c("parameter","length contig")
parameters <- unique(circ$parameter)

#gp1 <- ggplot(circ[which(circ$parameter=="108AT.UNMR.R1"),], aes(x=`length MC`)) + geom_histogram(binwidth = 1,color="yellow",fill="yellow",alpha=0.5)
#gp2 <- ggplot(circ[which(circ$parameter=="108AT.UNMR.R2"),], aes(x=`length MC`)) + geom_histogram(binwidth = 1,color="red",fill="red",alpha=0.5)
#gp3 <- ggplot(circ[which(circ$parameter=="108AT.UNMR.R3"),], aes(x=`length MC`)) + geom_histogram(binwidth = 1,color="blue",fill="blue",alpha=0.5)
#gp4 <- ggplot(circ[which(circ$parameter=="108AT.UNMR.R4"),], aes(x=`length MC`)) + geom_histogram(binwidth = 1,color="green",fill="green",alpha=0.5)
#gp5 <- ggplot(circ[which(circ$parameter=="108AT.UNMR.R5"),], aes(x=`length MC`)) + geom_histogram(binwidth = 1,color="orange",fill="orange",alpha=0.5)
#gp6 <- ggplot(circ[which(circ$parameter=="108AT.MR"),], aes(x=`length MC`)) + geom_histogram(binwidth = 1,color="grey",fill="grey",alpha=0.5)
#figure <- ggarrange(gp1, gp2, gp3, labels = c("R1", "R2", "R3","R4","R5","MR"), ncol = 3, nrow = 2)
#figure
#av_length <- ddply(circ,.(parameters), summarise, grp.mean=mean(`length contig`))
#head(av_length)
#geom_vline(data=av_length,aes(xintercept=grp.mean,color=parameter), linetype="dashed")
#ggplot(circ, aes(x=`length MC`, fill=sample, color=sample)

pdf("Nabe - different runs.pdf")

## HISTOGRAM freq/read length of 6 different runs 

gp1 <- ggplot(circ, aes(x=`length contig`, fill=parameter, color=parameter)) + 
  geom_histogram(position="identity", binwidth = 1,alpha=0.5) + 
  facet_grid(parameter ~ .) + 
  ggtitle("Histrogram: Frequency of circularized minicircle length per run") + 
  theme(plot.title = element_text(hjust=0.5, margin=margin(10,0,40,0)))
gp1

table(circ)

setwd("C:/Users/mgeerts/Desktop")

fasta <- read.dna("108AT.minicircles.diffR.fasta",'fasta')
### 774 sequences

mini <- as.numeric(gsub('.*_len|_cir.*','',attr(fasta, 'names')))
names(mini) <- gsub('_.*','',attr(fasta, 'names'))
df <- as.data.frame(matrix(nrow=length(mini), ncol=2))
df[,1] <- names(mini)
df[,2] <- mini

allruns <- unique(names(mini))

p.108AT <- ggplot(df, aes(x=V2, fill=V1, color=V1)) + xlab("minicircle length") + ggtitle("108AT") +
  geom_boxplot() 
p.108AT 

## BARPLOT number of circ contigs per run 

circ_summary <- matrix(nrow=length(parameters),ncol=3)
colnames(circ_summary) <- c("parameters" ,"circ contig len < 1100", "circ contigs len >= 1100")
for (i in 1:length(parameters)) {
  circ_summary[i,1] <- parameters[i]
  circ_summary[i,2] <- nrow(subset(circ, circ$parameter==parameters[i] & circ$length <= 1100))
  circ_summary[i,3] <- nrow(subset(circ, circ$parameter==parameters[i] & circ$length > 1100))
}
circ_summary <- as.data.frame(circ_summary)

circ_melted = melt(circ_summary, id=c('parameters'))
names(circ_melted) <- c('parameters', 'length', 'contigs')
circ_melted$contigs <- as.numeric(circ_melted$contigs)

ReadsFraction <- ddply(circ_melted,.(parameters),summarise,Count.Fraction = contigs / sum(contigs))
to_graph <- cbind(arrange(circ_melted, parameters), fraction = ReadsFraction$Count.Fraction)

gp2 <- ggplot(data=to_graph, aes(x=parameters, y=contigs, fill=length)) +
  geom_bar(stat="identity",position = "stack",alpha=0.5) +
  geom_text(aes(label=contigs),size = 3,vjust=5,color="white") +
  theme(axis.text.x=element_text(size=4)) +
  ggtitle("Barplot: Count of circularized minicircles per run") + 
  theme(plot.title = element_text(hjust=0.5, margin=margin(10,0,40,0)))
gp2


####### ALL MINICIRCLES

### DATA PREP

file <- read.table(file = "Nabe.length.minicircles.diffR.txt")
allmc <- matrix(ncol = 2,nrow=nrow(file))
colnames(allmc) <- c("parameter","length contig")
allmc <- as.data.frame(allmc)
allmc[2] <- as.numeric(gsub('.*_len|_cir.*','', file$V1))
allmc[1] <- gsub('_c.*','',file$V1)
allmc$parameter <- gsub('>','',allmc$parameter)
parameters <- unique(allmc$parameter)

## HISTOGRAM freq/read length of 6 different runs

gp3 <- ggplot(allmc, aes(x=`length contig`, fill=parameter, color=parameter)) +
  geom_histogram(position="identity", binwidth = 1,alpha=0.5) + 
  facet_grid(parameter ~ .) +
  ggtitle("Histrogram: Frequency of all minicircles length per run") + 
  theme(plot.title = element_text(hjust=0.5, margin=margin(10,0,40,0)))
gp3

## BARPLOT number of all minicircles per run 

table(allmc)

allmc_summary <- matrix(nrow=length(parameters),ncol=3)
colnames(allmc_summary) <- c("parameters" ,"all minicircles len < 1100", "all minicircles len >= 1100")
for (i in 1:length(parameters)) {
  allmc_summary[i,1] <- parameters[i]
  allmc_summary[i,2] <- nrow(subset(allmc, allmc$parameter==parameters[i] & allmc$length <= 1100))
  allmc_summary[i,3] <- nrow(subset(allmc, allmc$parameter==parameters[i] & allmc$length > 1100))
}
allmc_summary <- as.data.frame(allmc_summary)

allmc_melted = melt(allmc_summary, id=c('parameters'))
names(allmc_melted) <- c('parameters', 'length', 'contigs')
allmc_melted$contigs <- as.numeric(allmc_melted$contigs)

ReadsFraction <- ddply(allmc_melted,.(parameters),summarise,Count.Fraction = contigs / sum(contigs))
to_graph <- cbind(arrange(allmc_melted, parameters), fraction = ReadsFraction$Count.Fraction)

gp4 <- ggplot(data=to_graph, aes(x=parameters, y=contigs, fill=length)) +
  geom_bar(stat="identity",position = "stack",alpha=0.5) +
  geom_text(aes(label=contigs),size = 3,vjust=1,color="white") +
  theme(axis.text.x=element_text(size=4)) +
  ggtitle("Barplot: Count of all minicircles per run") + 
  theme(plot.title = element_text(hjust=0.5, margin=margin(10,0,40,0)))
#geom_text(aes(label=paste(round(fraction*100),"%", sep="")),size = 1,vjust=1,position="stack")
gp4

#figure <- ggarrange(gp1, gp2, gp3, gp4, labels = c("A", "B", "C","D"), ncol = 2, nrow = 2)
#figure

dev.off()



