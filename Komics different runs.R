library(ggplot2)
library(reshape2)
library(plyr)
library(ggpubr)
library(readxl)
library(ape)

#################################### MINICIRCLES 

setwd("C:/Users/mgeerts/Dropbox/Gambiense/Foefelare Manon/komics/trimming optim")

pdf("Different runs - 108AT.pdf", width=15)

##------------------------------------------------------------------------------------------------------------
## CALCULATIONS
##------------------------------------------------------------------------------------------------------------


#### FASTA 
dnamini <- read.dna('108AT.minicircles.diffR.fasta', 'fasta')
mini <- as.numeric(gsub('.*_len|_cir.*','',attr(dnamini, 'names')))
names(mini) <- gsub('_.*','',attr(dnamini, 'names'))
runs <- unique(names(mini))


#### UC
total_MSCs <- GAPs <- vector()
INS <- DEL <- list()
uc <- read.table("108AT.minicircles.diffR.uc")
uc_C <- uc[which(uc$V1=='C'),]
uc_C2 <- uc_C[,2]
uc_H <- uc[which(uc$V1 =='H'),]
  
total_MSCs <- length(uc_C2)
  
uc_H2 <- uc_H[-c(which(as.numeric(gsub('.*_len|_circ.*','',uc_H$V9))>1200), which(as.numeric(gsub('.*_len|_circ.*','',uc_H$V10))>1200)),]
perfectalignments <- sort(c(which(apply(uc_H, 1, function(x) (x[8] == paste(x[3], 'M', sep='')))), which(uc_H$V8 == "=")))
GAPs <- 1-length(perfectalignments)/length(uc_H$V8)
uc_H3 <- uc_H[-perfectalignments,]
insertions <- gsub('[0-9]*D','',gsub('[0-9]*M','',uc_H3$V8))
INS <- table(insertions)
names(INS) <- gsub('I','',names(INS))
deletions <- gsub('[0-9]*I','',gsub('[0-9]*M','',uc_H3$V8))
DEL <- table(deletions)
names(DEL) <- gsub('D','',names(DEL))

clusters <- vector(mode = 'list', length = length(uc_C2))
for (i in 1:length(uc_C2)) {
  if (uc_C[uc_C$V2==uc_C2[i], 3] == 1) {
    clusters[[i]] <- as.character(uc_C[uc_C$V2==uc_C2[i], 9])
  } else {
    clusters[[i]] <- unique(c(as.character(uc_H[uc_H$V2==uc_C2[i], 9]), as.character(uc_H[uc_H$V2==uc_C2[i],10])))
  }
}

clust_mat <- matrix(nrow=length(clusters), ncol = length(runs))
colnames(clust_mat) <- runs
rownames(clust_mat) <- paste('C', uc_C[,2], sep = '')

for (t in 1:length(runs)) {
  temp <- lapply(lapply(clusters, function(x) table(gsub('_contig.*','',x))), function(x) x[which(names(x) == runs[t])])
  
  for (ct in 1:length(temp)) {
    if (length(temp[[ct]]) > 0) {
      clust_mat[ct,t] <- temp[[ct]]
    } else if (length(temp[[ct]]) == 0) {
      clust_mat[ct,t] <- 0 #temp[[ct]]
    }
  }
}


##------------------------------------------------------------------------------------------------------------
## BOXPLOT ALL MC
##------------------------------------------------------------------------------------------------------------

df <- as.data.frame(matrix(nrow=length(mini), ncol=2))
df[,1] <- names(mini)
df[,2] <- mini
colnames(df) <- c("run","length MC")

box <- ggplot(df, aes(x=`length MC`, fill=run, color=run)) + xlab("minicircle length") +
  ggtitle("108AT: Frequency minicircle length") + geom_boxplot(alpha=0.5) 
box


##------------------------------------------------------------------------------------------------------------
## BARPLOT ALL MC
##------------------------------------------------------------------------------------------------------------

df_summary <- as.data.frame(matrix(nrow=length(runs),ncol=3))
colnames(df_summary) <- c("run" ,"MC len < 1200", "MC len >= 1200")
for (i in 1:length(runs)) {
  df_summary[i,1] <- runs[i]
  df_summary[i,2] <- nrow(subset(df, df$run==runs[i] & df$`length MC` <= 1200))
  df_summary[i,3] <- nrow(subset(df, df$run==runs[i] & df$`length MC` > 1200))
}
melted <- melt(df_summary, id='run')
names(melted) <- c('run', 'length', 'number of MC')
#allmc_melted$contigs <- as.numeric(allmc_melted$contigs)

bar <- ggplot(melted, aes(x=run, y=`number of MC`, fill=length)) + 
  geom_bar(stat="identity",position = "stack",alpha=0.5) + xlab('') +
  geom_text(aes(label=`number of MC`),size = 3,vjust=1,color="white") +
  theme(axis.text.x=element_text(size=8, angle = 45)) + ggtitle("108AT: Number of minicircles per run") 
  #theme(plot.title = element_text(hjust=0.5, margin=margin(10,0,40,0)))

N_between <- as.data.frame(matrix(ncol=2, nrow=length(runs)))
N_between[,1] <- runs
N_between[,2] <- colSums(clust_mat)

point <- ggplot(N_between, aes(x = V1, y=V2, label=V2,color=V1)) + geom_point() +
  xlab('') + ylim(0,160) + ylab('') + geom_text(vjust = 0, nudge_y = 1) + 
  theme(axis.text.x=element_text(size=8,angle=45)) + ggtitle("idem")

figure <- ggarrange(bar,point)
figure


##------------------------------------------------------------------------------------------------------------
## BOXPLOT CIRC MC
##------------------------------------------------------------------------------------------------------------

circ <- read.table(file = "108AT.length.circ.minicircles.diffR.txt", sep="_")
circ$V3 <- as.numeric(gsub("len", "", circ$V3))
circ$V1 <- gsub(">","",circ$V1)
circ <- circ[,c(1,3)]
colnames(circ) <- c("run","length MC")
runs <- unique(circ$run)

box <- ggplot(circ, aes(x=`length MC`, fill=run, color=run)) + xlab("minicircle length") +
  ggtitle("108AT: Frequency circularized minicircle length") + geom_boxplot(alpha=0.5) 
box


##------------------------------------------------------------------------------------------------------------
## BARPLOT CIRC MC
##------------------------------------------------------------------------------------------------------------


circ_summary <- as.data.frame(matrix(nrow=length(runs),ncol=3))
colnames(circ_summary) <- c("run" ,"circ MC len < 1200", "circ MC len >= 1200")
for (i in 1:length(runs)) {
  circ_summary[i,1] <- runs[i]
  circ_summary[i,2] <- nrow(subset(circ, circ$run==runs[i] & circ$`length MC` <= 1200))
  circ_summary[i,3] <- nrow(subset(circ, circ$run==runs[i] & circ$`length MC` > 1200))
}

circ_melted = melt(circ_summary, id='run')
names(circ_melted) <- c('run', 'length', 'number of MC')

bar <- ggplot(circ_melted, aes(x=run, y=`number of MC`, fill=length)) + 
  geom_bar(stat="identity",position = "stack",alpha=0.5) + xlab('') +
  geom_text(aes(label=`number of MC`),size = 3,vjust=1,color="white") +
  theme(axis.text.x=element_text(size=8, angle=45)) + ggtitle("108AT: Number of circ minicircles per run") 

point <- ggplot(circ_melted, aes(x=run, y=`number of MC`, fill=length, color=run, label=`number of MC`)) + 
  geom_point(stat="identity",alpha=0.5) + xlab('') + geom_text(vjust = 0, nudge_y = 1) + 
  theme(axis.text.x=element_text(size=8, angle = 45)) + ggtitle("idem") 

figure <- ggarrange(bar,point)
figure


##------------------------------------------------------------------------------------------------------------
## HEATMAP
##------------------------------------------------------------------------------------------------------------


melted <- melt(clust_mat)

heatmap <- ggplot(melted, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + ylab('') + xlab('') +
  scale_fill_gradient(low="dodgerblue3", high="firebrick") + ggtitle("108AT: Presence of clusters in different runs") +
  theme(axis.text.x = element_text(angle = 90, size=6))
heatmap

dev.off()
