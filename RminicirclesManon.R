##############################################################################################################
#  MINICIRCLES 
##############################################################################################################

library(ggplot2)
library(reshape2)
library(plyr)
library(ggpubr)
library(readxl)
library(ape)
library(openxlsx)
library(network)
library(pegas)
library(ggrepel)
library(gridExtra)
library(RColorBrewer)

##------------------------------------------------------------------------------------------------------------
## species2
##------------------------------------------------------------------------------------------------------------

setwd('C:/Users/mgeerts/Dropbox/Gambiense/Foefelare Manon/SNP_calling')

species2 <- read.csv("list sub species.csv", sep=',', header=T); species2 <- species2[,-1]
species2 <- species2[order(species2$strain),]
colnames(species2)[6] <- "DNA quantity (µg)" 
colnames(species2)[7] <- "test result BGI"
colnames(species2)[2] <- "sub species"

Tbb <- as.character(species2[which(species2$`sub species`=='T.b. brucei'),1])
Tbg <- as.character(species2[which(species2$`sub species`=="T.b. gambiense"),1])
TbgII <- as.character(species2[which(species2$`sub species`=="T.b. gambiense II"),1])
Tbr <- as.character(species2[which(species2$`sub species`=="T.b. rhodesiense"),1])
uncertain <- as.character(species2[which(species2$`sub species`=='uncertain'),1])

strains <- c(Tbb,Tbg,TbgII,Tbr,uncertain)

subset <- c("PTAG_129", "1829_Alijo", "Fontem_S10", "KP33_clone_16","TSW_187_78E")
species2 <- species2[-which(species2$strain %in% subset),] ## 38
Tbb <- as.character(species2[which(species2$`sub species`=='T.b. brucei'),1])
Tbg <- as.character(species2[which(species2$`sub species`=="T.b. gambiense"),1])
TbgII <- as.character(species2[which(species2$`sub species`=="T.b. gambiense II"),1])
Tbr <- as.character(species2[which(species2$`sub species`=="T.b. rhodesiense"),1])
uncertain <- as.character(species2[which(species2$`sub species`=='uncertain'),1])

strains <- c(Tbb,Tbg,TbgII,Tbr,uncertain)

write.csv(table(species2$`sub species`), 
          'C:/Users/mgeerts/Dropbox/Gambiense/Foefelare Manon/komics/minicircles/species2subset.csv')

##------------------------------------------------------------------------------------------------------------
## select MC: circ seq and seq with 800>len>1200 (based on fasta files)
##------------------------------------------------------------------------------------------------------------

setwd('C:/Users/mgeerts/Desktop/minicircles/')

files <- list.files(pattern ='fasta')
MC_sel <- as.data.frame(matrix(nrow=length(files), ncol=4))
colnames(MC_sel) <- c("total n of MC", "n of MC with len<1200", "n of MC with len>800", "n of circ MC") 
rownames(MC_sel) <- gsub("\\..*", '', files)
tmp.dna <- list()
for (i in 1:length(files)) {
  tmp.dna <- read.dna(files[i], 'fasta')
  MC_sel[i,1] <- length(tmp.dna)
  tmp.dna <- tmp.dna[which(as.numeric(gsub('.*_len|_cir.*','',attr(tmp.dna, 'names')))<1200)]
  MC_sel[i,2] <- length(tmp.dna)
  tmp.dna <- tmp.dna[which(as.numeric(gsub('.*_len|_cir.*','',attr(tmp.dna, 'names')))>800)]
  MC_sel[i,3] <- length(tmp.dna)
  tmp.dna <- tmp.dna[grep("circ", attr(tmp.dna, 'names'))]
  MC_sel[i,4] <- length(tmp.dna)
  write.dna(tmp.dna, paste("selected MC/",gsub("\\..*", '', files[i]),'.minicircles.circ.fasta',sep=''), "fasta")
}
MC_sel <- cbind(MC_sel, apply(MC_sel, 1, function(x) x[4]/x[1]))
colnames(MC_sel)[5] <- "% of selected MC"; MC_sel <- MC_sel[order(rownames(MC_sel)),]
MC_sel <- rbind(MC_sel, colSums(MC_sel)); rownames(MC_sel)[39] <- "total" 

write.csv(MC_sel, 'C:/Users/mgeerts/Dropbox/Gambiense/Foefelare Manon/komics/MC_sel.csv')

###### perform vsearch

##------------------------------------------------------------------------------------------------------------
## stats (based on fasta files)
##------------------------------------------------------------------------------------------------------------
## ----
### MINI STATS
## ----
setwd('C:/Users/mgeerts/Desktop/minicircles/selected MC')
files <- list.files(pattern ='fasta')
ministats <- vector()
for (n in 1:length(files)) {
  ministats[n] <- length(read.dna(files[n], 'fasta'))
}
names(ministats) <- gsub('.mini.*','',files)
ministats <- as.data.frame(ministats)
colnames(ministats) <- 'n of MC'  #rownames(ministats) == species2$strain
#write.csv(ministats, 'C:/Users/mgeerts/Dropbox/Gambiense/Foefelare Manon/komics/Nmini.csv', col.names=F)

setwd("C:/Users/mgeerts/Dropbox/Gambiense/")
pdf("Foefelare Manon/komics/figures/A.N of MC.pdf", useDingbats = F, width = 15)
g <- ggplot(ministats, aes(x=rownames(ministats), y=`n of MC`, color=species2$`sub species`, fill=species2$`sub species`)) + geom_bar(alpha=0.5, show.legend = F, stat="identity") + 
  xlab('') + ylab('n of MC') + ggtitle("Number of MC") + facet_grid(. ~ species2$`sub species`, scale='free', space='free') + 
  theme(axis.text.x=element_text(angle=90, hjust=1), strip.text.x = element_text(size=8))
g; dev.off()

## ----
# ALL STATS CSB3
## ----
allstats <- read.xlsx("RESULTS.TABLES.xlsx", sheet = "Sheet1", rowNames=T)
colnames(allstats) <- gsub('\\.', ' ',colnames(allstats))

g <- ggplot(allstats, aes(x=`total number of unmapped reads (UNMR)`, y=`UNMR with CSB3`, color=species2$`sub species`)) + 
  geom_point(size=2, alpha=0.5) + geom_text_repel(label=rownames(allstats), size=3) + ggtitle("Relation between n of UNMR and n of UNMR with CBS3")
g
g <- ggplot(allstats, aes(x=`total number of unmapped reads (UNMR)`, y=`n of MC`, color=species2$`sub species`)) + 
  geom_point(size=2, alpha=0.5) + 
  geom_text_repel(label=rownames(allstats), size=3) + ggtitle("Relation between n of UNMR and n of MC")
g; dev.off()

## ----
### DNA QUANTITY
## ----

DNAquantity <- as.data.frame(matrix(nrow=38, ncol=2))
DNAquantity[,1] <- MC_sel[1:38,"n of circ MC"]; rownames(DNAquantity) <- rownames(MC_sel)[-39]  
DNAquantity[,2] <- species2[,"DNA quantity (µg)"]; #DNAquantity[,3] <- rownames(MC_sel)[-39]

# outliers from previous barplot
subset <- c("MBO-NG-74-R10","Bosendja","LiTat-1-5-P9","MBOT-GM-77-GB2","MSUS-CI-78-TSW382")

mean(DNAquantity[which(rownames(DNAquantity)%in%Tbg),2]); mean(DNAquantity[which(rownames(DNAquantity)%in%Tbg),1])
mean(DNAquantity[which(rownames(DNAquantity)%in%TbgII),2]);mean(DNAquantity[which(rownames(DNAquantity)%in%TbgII),1])
mean(DNAquantity[which(rownames(DNAquantity)%in%Tbb),2]);mean(DNAquantity[which(rownames(DNAquantity)%in%Tbb),1])
mean(DNAquantity[which(rownames(DNAquantity)%in%Tbr),2]);mean(DNAquantity[which(rownames(DNAquantity)%in%Tbr),1])
mean(DNAquantity[which(rownames(DNAquantity)%in%uncertain),2]);mean(DNAquantity[which(rownames(DNAquantity)%in%uncertain),1])

library(scales)
show_col(hue_pal()(5))
mycolors <- hue_pal()(5)

pdf("Foefelare Manon/komics/figures/X.DNA quantity.pdf", useDingbats = F, width = 15)
g <- ggplot(DNAquantity, aes(x=V1, y=V2, color=species2$`sub species`, shape=species2$`sub species`)) + ggtitle("DNA quantity (µg)") +
  geom_point(size=4, alpha=0.5) + geom_text_repel(label=ifelse(rownames(DNAquantity)%in%subset, rownames(DNAquantity),'')) + 
  theme(legend.title = element_blank()) + xlab('n of MC') + ylab('DNA quantity (µg)') +
  
g; dev.off()





##------------------------------------------------------------------------------------------------------------
## minicircle length (based on fasta)
##------------------------------------------------------------------------------------------------------------
setwd("C:/Users/mgeerts/Dropbox/Gambiense/Foefelare Manon/komics/minicircles")

dnamini <- read.dna('all.minicircles.circ.BGI.fasta', 'fasta')
mini <- as.numeric(gsub('.*_len|_cir.*','',attr(dnamini, 'names')))
mini <- as.data.frame(mini) 
mini[,2] <- gsub('_co.*','',attr(dnamini, 'names'))
#sort(unique(gsub('_c.*','',attr(dnamini, 'names')))) == sort(strains)

mini[which(mini$V2%in%Tbg),3] <- "T.b. gambiense (n=18)"; mini[which(mini$V2%in%TbgII),3] <- "T.b. gambiense II (n=3)"
mini[which(mini$V2%in%Tbr),3] <- "T.b. rhodesiense (n=4)"; mini[which(mini$V2%in%Tbb),3] <- "T.b. brucei (n=7)"
mini[which(mini$V2%in%uncertain),3] <- "uncertain (n=6)"
colnames(mini) <- c("MC length", "strain", "sub species")

pdf('figures/B.Frequency MC length.pdf', useDingbats = F, width = 15)
g <- ggplot(mini, aes(x=`MC length`, fill=`sub species`, color=`sub species`)) +  ggtitle('Frequency MC length') + 
  geom_histogram(binwidth=10, alpha=0.5, show.legend = F) + facet_grid(. ~ `sub species`) 
g; dev.off() 

##------------------------------------------------------------------------------------------------------------
## calculate some stuff (based on uc)
##------------------------------------------------------------------------------------------------------------
setwd('C:/Users/mgeerts/Dropbox/Gambiense/Foefelare Manon/komics/minicircles/')


id <- c(70,80,90,95,96,97,98,99,100)
MSCs <- GAPs <- MSCs.tbg <- gaps.tbg <- vector()
INS <- DEL <- INS.tbg <- DEL.tbg <- clust_mat_ids <- list()
for (n in 1:length(id)) {
  uc <- read.table(paste('all.minicircles.circ.BGI.id',id[n], '.uc', sep=''))
  uc_C <- uc[which(uc$V1=='C'),]
  uc_C2 <- uc_C[,2]
  uc_C.tbg <- na.omit(uc_C[which(gsub('_co.*','',uc$V9)%in% Tbg),])
  uc_C2.tbg <- uc_C.tbg[,2]
  uc_H <- uc[which(uc$V1 =='H'),]
  uc_H.tbg <- na.omit(uc_H[which(gsub('_co.*','',uc$V10)%in% Tbg),])
  
  MSCs[n] <- length(uc_C2)
  MSCs.tbg[n] <- length(uc_C2.tbg)
  
  perfectalignments <- sort(c(which(apply(uc_H, 1, function(x) (x[8] == paste(x[3], 'M', sep='')))), which(uc_H$V8 == "=")))
  perfectalignments.tbg <- sort(c(which(apply(uc_H.tbg, 1, function(x) (x[8] == paste(x[3], 'M', sep='')))), which(uc_H.tbg$V8 == "=")))
  GAPs[n] <- 1-length(perfectalignments)/length(uc_H$V8)
  gaps.tbg[n] <- 1-length(perfectalignments.tbg)/length(uc_H.tbg$V8)
  if (id[n] != 100) {
    uc_H3 <- uc_H[-perfectalignments,]
    insertions <- gsub('[0-9]*D','',gsub('[0-9]*M','',uc_H3$V8))
    INS[[n]] <- table(insertions)
    names(INS[[n]]) <- gsub('I','',names(INS[[n]]))
    deletions <- gsub('[0-9]*I','',gsub('[0-9]*M','',uc_H3$V8))
    DEL[[n]] <- table(deletions)
    names(DEL[[n]]) <- gsub('D','',names(DEL[[n]]))

    uc_H3.tbg <- uc_H.tbg[-perfectalignments.tbg,]
    insertions.tbg <- gsub('[0-9]*D','',gsub('[0-9]*M','',uc_H3.tbg$V8))
    INS.tbg[[n]] <- table(insertions.tbg)
    names(INS.tbg[[n]]) <- gsub('I','',names(INS.tbg[[n]]))
    deletions.tbg <- gsub('[0-9]*I','',gsub('[0-9]*M','',uc_H3.tbg$V8))
    DEL.tbg[[n]] <- table(deletions.tbg)
    names(DEL.tbg[[n]]) <- gsub('D','',names(DEL.tbg[[n]]))
  }
  
  clusters <- vector(mode = 'list', length = length(uc_C2))
  for (i in 1:length(uc_C2)) {
    if (uc_C[uc_C$V2==uc_C2[i], 3] == 1) {
      clusters[[i]] <- as.character(uc_C[uc_C$V2==uc_C2[i], 9])
    } else {
      clusters[[i]] <- unique(c(as.character(uc_H[uc_H$V2==uc_C2[i], 9]), as.character(uc_H[uc_H$V2==uc_C2[i],10])))
    }
  }
  
  clust_mat <- matrix(nrow=length(clusters), ncol = length(strains))
  colnames(clust_mat) <- strains
  rownames(clust_mat) <- paste('C', uc_C[,2], sep = '')
  for (t in 1:length(strains)) {
    temp <- lapply(lapply(clusters, function(x) table(gsub('_contig.*','',x))), function(x) x[which(names(x) == strains[t])])
    
    for (ct in 1:length(temp)) {
      if (length(temp[[ct]]) > 0) {
        clust_mat[ct,t] <- temp[[ct]]
      } else if (length(temp[[ct]]) == 0) {
        clust_mat[ct,t] <- 0 #temp[[ct]]
      }
    }
  }
  clust_mat_ids[[n]] <- clust_mat
}


##------------------------------------------------------------------------------------------------------------
## Number of minicircle sequence classes per % identity
##------------------------------------------------------------------------------------------------------------

#for (i in 1:length(id)) {
#  tmp <- apply(clust_mat_ids[[i]][,Tbg], 1, function(x) sum(x!=0))
#  MSCs_tbg[i,1] <- length(tmp)
#  MSCs_tbg[i,2] <- sum(tmp!=0)
#}

pdf('figures/C.MSC and gaps per id.pdf', useDingbats = F, width = 15)
MSCs <- cbind(MSCs,MSCs.tbg); rownames(MSCs) <- id
melted <- melt(MSCs); colnames(melted) <- c("% identity", "MSC", "n of MSC")
g <- ggplot(melted, aes(x=`% identity`, y=`n of MSC`, col=MSC, fill=MSC)) + geom_point(alpha=0.5, size=3) + geom_line() + 
  theme(legend.title = element_blank()) + ggtitle('N of MCSs per % identity')
g
##------------------------------------------------------------------------------------------------------------
## Number of gaps per % identity
##------------------------------------------------------------------------------------------------------------

namesINS <- sort(unique(do.call(c,lapply(INS, names))))[-1]
matrixINS <- matrix(ncol = length(namesINS), nrow = length(INS))
colnames(matrixINS) <- namesINS
rownames(matrixINS) <- id[-9]
for (i in 1:length(INS)) {
  for (t in 1:length(names(INS[[i]]))) {
    if (names(INS[[i]])[t] %in% namesINS) {
      matrixINS[i,names(INS[[i]])[t]] <- INS[[i]][t]
    }
  }
}

namesDEL <- sort(unique(do.call(c,lapply(DEL, names))))[-1]
matrixDEL <- matrix(ncol = length(namesDEL), nrow = length(DEL))
colnames(matrixDEL) <- namesDEL
rownames(matrixDEL) <- id[-9]
for (i in 1:length(DEL)) {
  for (t in 1:length(names(DEL[[i]]))) {
    if (names(DEL[[i]])[t] %in% namesDEL) {
      matrixDEL[i,names(DEL[[i]])[t]] <- DEL[[i]][t]
    }
  }
}

matrixINS[is.na(matrixINS)] <- 0
matrixDEL[is.na(matrixDEL)] <- 0

#write.csv2(matrixINS, 'INSERTIONS.csv')
#write.csv2(matrixDEL, 'DELETIONS.csv')

namesINS.tbg <- sort(unique(do.call(c,lapply(INS.tbg, names))))[-1]
matrixINS.tbg <- matrix(ncol = length(namesINS.tbg), nrow = length(INS.tbg))
colnames(matrixINS.tbg) <- namesINS.tbg
rownames(matrixINS.tbg) <- id[-9]
for (i in 1:length(INS.tbg)) {
  for (t in 1:length(names(INS.tbg[[i]]))) {
    if (names(INS.tbg[[i]])[t] %in% namesINS.tbg) {
      matrixINS.tbg[i,names(INS.tbg[[i]])[t]] <- INS.tbg[[i]][t]
    }
  }
}

namesDEL.tbg <- sort(unique(do.call(c,lapply(DEL.tbg, names))))[-1]
matrixDEL.tbg <- matrix(ncol = length(namesDEL.tbg), nrow = length(DEL.tbg))
colnames(matrixDEL.tbg) <- namesDEL.tbg
rownames(matrixDEL.tbg) <- id[-9]
for (i in 1:length(DEL.tbg)) {
  for (t in 1:length(names(DEL.tbg[[i]]))) {
    if (names(DEL.tbg[[i]])[t] %in% namesDEL.tbg) {
      matrixDEL.tbg[i,names(DEL.tbg[[i]])[t]] <- DEL.tbg[[i]][t]
    }
  }
}

matrixINS.tbg[is.na(matrixINS.tbg)] <- 0
matrixDEL.tbg[is.na(matrixDEL.tbg)] <- 0


##------------------------------------------------------------------------------------------------------------
## Length of gaps per % identity
##------------------------------------------------------------------------------------------------------------

INSlength <- matrix(ncol = 10, nrow = length(INS))
for (n in 1:10) {
  INSlength[,n] <- unlist(lapply(INS, function(x) sum(x[which(names(x) == n)])))
}

DELlength <- matrix(ncol = 10, nrow = length(DEL))
for (n in 1:10) {
  DELlength[,n] <- unlist(lapply(DEL, function(x) sum(x[which(names(x) == n)])))
}

gaps <- INSlength+DELlength
G1 <- apply(gaps, 1, function(x) x[2]/sum(x))
G2 <- apply(gaps, 1, function(x) x[3]/sum(x))

INS.length.tbg <- matrix(ncol = 10, nrow = length(INS.tbg))
for (n in 1:10) {
  INS.length.tbg[,n] <- unlist(lapply(INS.tbg, function(x) sum(x[which(names(x) == n)])))
}

DEL.length.tbg <- matrix(ncol = 10, nrow = length(DEL.tbg))
for (n in 1:10) {
  DEL.length.tbg[,n] <- unlist(lapply(DEL.tbg, function(x) sum(x[which(names(x) == n)])))
}

gaps.tbg <- INS.length.tbg+DEL.length.tbg
G1.tbg <- apply(gaps.tbg, 1, function(x) x[2]/sum(x))
G2.tbg <- apply(gaps.tbg, 1, function(x) x[3]/sum(x))

rownames(gaps.tbg) <- rownames(gaps) <- id[-9]
colnames(gaps.tbg) <- colnames(gaps) <- c("1st gap", "2nd gap", "3rd gap", "4th gap", "5th", "6th", "7th", "8th","9th","10th")
melted <- melt(gaps)
melted <- cbind(melted, "gaps all species")
melted.tbg <- melt(gaps.tbg)
melted.tbg <- cbind(melted.tbg, "gaps tbg")
colnames(melted)[4] <- colnames(melted.tbg)[4] <- "gaps"
melted <- rbind(melted, melted.tbg); colnames(melted)[1:3] <- c("% identity", "nth gap", "n of alignments") 

g2 <- ggplot(melted[melted$`nth gap`%in%c("2nd gap","3rd gap"),], aes(x=`% identity`, y=`n of alignments`, color=gaps, fill=`nth gap`, shape=`nth gap`)) + 
  geom_point(size=3, alpha=0.5) + geom_line() + theme(legend.title = element_blank()) +  ggtitle("N of gaps")
g2
figure <- ggarrange(g,g2)
figure; dev.off()


##------------------------------------------------------------------------------------------------------------
## Number of minicircles averaged across % identity
##------------------------------------------------------------------------------------------------------------

N_between <- matrix(ncol = length(strains), nrow = length(id))
colnames(N_between) <- strains
rownames(N_between) <- id
for (i in 1:length(id)) {
  N_between[i,] <- apply(clust_mat_ids[[i]], 2, function(x) sum(x>0))
}

melted <- melt(N_between)
melted[which(melted$Var2%in%Tbg),4] <- "T.b. gambiense (n=18)"; melted[which(melted$Var2%in%TbgII),4] <- "T.b. gambiense II (n=3)"
melted[which(melted$Var2%in%Tbr),4] <- "T.b. rhodesiense (n=4)"; melted[which(melted$Var2%in%Tbb),4] <- "T.b. brucei (n=7)"
melted[which(melted$Var2%in%uncertain),4] <- "uncertain (n=6)"
colnames(melted) <- c("% identity", "strain", "Number of MSC", "sub species")

g <- ggplot(melted, aes(x=strain, y=`Number of MSC`, color=`sub species`, fill=`sub species`)) + geom_boxplot(alpha=0.5, show.legend = F) + 
  theme(axis.text.x = element_text(angle=90, hjust=1)) + facet_grid(. ~ `sub species`, scale='free', space='free') + xlab('') + 
  ggtitle('Number of MSC per sub species') + theme(strip.text.x = element_text(size=6))
g

range <- as.data.frame(apply(N_between, 2, function(x) x[9]-x[1]))
range <- cbind(range, colSums(N_between)/9)
range <- range[order(rownames(range)), ]
g3 <- ggplot(range, aes(x=range[,2], y=range[,1], color=species2$`sub species`, shape=species2$`sub species`)) + 
  geom_point(size=3, alpha=0.5, show.legend = F) + ylab('range max min') + xlab('average n of MSC') + 
  ggtitle('Average n of MSC vs range max-min') +  theme(plot.margin = margin(0.2, 1, 5, 1, "cm"))
  #facet_wrap(. ~ species2$`sub species`) +geom_smooth(se = FALSE, alpha=0.5, size=0.2)
g3

figure <- ggarrange(g,g3, ncol=2, nrow=1, widths = c(2, 1))
pdf('figures/D.n of MSC per sub species across ids.pdf', useDingbats = F, width = 15)
figure

# line plots 
g4 <- ggplot(melted, aes(x = `% identity`, y = `Number of MSC`, fill = strain, color=`sub species`)) + 
  geom_line() +  ggtitle('Number of MSC per sub species')
g4; dev.off()



##------------------------------------------------------------------------------------------------------------
## Number of shared and unique minicircles across % id
##------------------------------------------------------------------------------------------------------------

N_shared <- matrix(ncol = length(strains), nrow = length(id))
colnames(N_shared) <- strains
Tbgcirc <- allnonTbgcirc <- circles <- list()
for (i in 1:length(id)) {
  Tbgcirc[[i]] <- apply(clust_mat_ids[[i]][,Tbg],1, function(x) sum(x!=0))/length(Tbg)
  #TbgIIcirc[[i]] <- apply(clust_mat_ids[[i]][,TbgII],1, function(x) sum(x!=0))/length(TbgII)
  allnonTbgcirc[[i]] <- apply(clust_mat_ids[[i]][,c(TbgII,Tbr,Tbb,uncertain)],1, function(x) sum(x!=0))/length(c(TbgII,Tbr,Tbb,uncertain))
  #unknowncirc[[i]] <- apply(clust_mat_ids[[i]][,unknown],1, function(x) sum(x!=0))/length(unknown)

  tmp <- rbind(Tbgcirc[[i]],allnonTbgcirc[[i]]) #,allnonTbgcirc[[i]],unknowncirc[[i]])
  circles[[i]] <- table(apply(tmp, 2, function(x) paste(which(x != 0), collapse='')))
}

circles <- lapply(circles, function(x) x/sum(x))
groups <- c("1","12","2")
shared <- matrix(ncol = length(groups),nrow=length(circles))
colnames(shared) <- groups
for (i in 1:length(circles)) {
  for (n in 1:length(groups)) {
    res <- circles[[i]][which(names(circles[[i]]) == groups[n])]
    if (length(res) > 0) shared[i,groups[n]] <- res
  }
}
apply(shared, 1, function(x) sum(x, na.rm = T))
shared[is.na(shared)] <- 0

shared <- as.data.frame(shared); colnames(shared) <- c("Tbg", "shared", "all non-Tbg")
melted <- melt(shared); melted[,3] <- seq(1,9,by=1); colnames(melted) <- c("shared","frequency", "% identity")

pdf('figures/F.MSC shared vs unique per id.pdf', useDingbats = F)
g <- ggplot(melted, aes(x=`% identity`,y=frequency, fill=shared)) + geom_bar(stat="identity",alpha=0.8) + 
  theme(legend.title = element_blank(), axis.text.x = element_blank()) + ylab('relative frequency') + 
  ggtitle('N of shared vs unique MSC') 
g; dev.off()



##------------------------------------------------------------------------------------------------------------
## Number of MSCs present in every strain of sub specie per % id
##------------------------------------------------------------------------------------------------------------

t <- matrix(nrow=length(id), ncol=5) 
rownames(t) <- id
for (i in 1:9) {
  t[i,1] <- length(which(apply(clust_mat_ids[[i]][,Tbb], 1, function(x) sum(x!=0)) == length(Tbb))) 
  t[i,2] <- length(which(apply(clust_mat_ids[[i]][,Tbg], 1, function(x) sum(x!=0)) == length(Tbg)))
  t[i,3] <- length(which(apply(clust_mat_ids[[i]][,TbgII], 1, function(x) sum(x!=0)) == length(TbgII)))
  t[i,4] <- length(which(apply(clust_mat_ids[[i]][,Tbr], 1, function(x) sum(x!=0)) == length(Tbr)))
  t[i,5] <- length(which(apply(clust_mat_ids[[i]][,uncertain], 1, function(x) sum(x!=0)) == length(uncertain)))
}

colnames(t) <- c("Tbb", 'Tbg', "TbgII", "Tbr", "uncertain")
write.csv(t, "n of specific MSC with full coverage per sub specie.csv")
melted <- melt(t)

g <- ggplot(melted, aes(Var1, value, fill=Var1, color=Var2)) + geom_point(size=3, alpha=0.5, show.legend = F) + geom_line() + 
  xlab('% identity') + ylab('n of MSCs') + ggtitle('n of MCSs - full coverage sub species')
g
g2 <- ggplot(melted[melted$Var2 != "TbgII",], aes(Var1, value, fill=Var1, color=Var2)) + geom_point(size=3, alpha=0.5, show.legend = F) + geom_line() + 
  xlab('% identity') + ylab('n of MSCs') + ggtitle('n of MCSs - full coverage sub species')
g2

setwd("C:/Users/mgeerts/Dropbox/Gambiense/Foefelare Manon/komics/")
pdf("figures/H.N of MCS present in every strain of specie per id.pdf", useDingbats = F, width = 15)
figure <- ggarrange(g,g2)
figure;dev.off()

##------------------------------------------------------------------------------------------------------------
## select those clusters with full coverage tbg in other sub species
##------------------------------------------------------------------------------------------------------------

setwd('C:/Users/mgeerts/Dropbox/Gambiense/Foefelare Manon/komics/minicircles')

tbgC <- list()
clust_tbg <- tmp_tbg <- matrix(ncol=38)
for (n in 1:length(id)) {
  tbgC[[n]] <- which(apply(clust_mat_ids[[n]][,Tbg], 1, function(x) sum(x!=0)) == length(Tbg))
  tmp_tbg <- as.matrix(clust_mat_ids[[n]][tbgC[[n]],])
  clust_tbg <- rbind(clust_tbg, tmp_tbg)
}

clust_tbg <- clust_tbg[-1,]
clust_tbg <- rbind(clust_tbg, clust_mat_ids[[8]][tbgC[[8]],]) 
rownames(clust_tbg)[38] <- "C2095"
colnames(clust_tbg) <- colnames(clust_mat_ids[[8]])

C_contigs <- H_contigs <- list()
for (n in 1:8) { #length(tbgC)) {
  uc <- read.table(paste('all.minicircles.circ.BGI.id',id[n], '.uc', sep=''))
  uc_C <- uc[which(uc$V1=='C'),]
  uc_H <- uc[which(uc$V1=='H'),]
  dat <- uc_C[which(paste('C', uc_C$V2, sep='') %in% names(tbgC[[n]])),]
  dat2 <- uc_H[which(paste('C', uc_H$V2, sep='') %in% names(tbgC[[n]])),]
  C_contigs[[n]] <- as.character(dat$V9)
  H_contigs[[n]] <- unique(c(as.character(dat2$V9),as.character(dat2$V10)))
  #dna2 <- dnamini[C_contigs]
  #names(dna2) <- paste(gsub('_contig.*','',names(dna2)),id[n],names(tbgC[[n]]), sep='.')
  #write.dna(dna2, paste('tbgclusters.id',id[n],'.fasta',sep=''), format = 'fasta')
}


dna2 <- dnamini[unique(C_contigs_all)]
write.dna(dna2,"all.tbgclusters.id70-99.fasta","fasta")

## align & splitstree

clust_tbg <- as.data.frame(clust_tbg)
clust_tbg[,39] <- unlist(C_contigs)
y <- list()
for (n in 1:8) {y[[n]] <- rep(rownames(t)[n], t[n,"Tbg"])}
clust_tbg[,40] <- unlist(y)
clust_tbg[,41] <- rownames(clust_tbg)


uc <- read.table(paste('all.minicircles.circ.BGI.id95.uc', sep=''))
uc_C <- uc[which(uc$V1=='C'),]
uc_H <- uc[which(uc$V1=='H'),]
dat2 <- uc_H[which(paste('C', uc_H$V2, sep='') == "C1107"),]
HIT <- unique(c(as.character(dat2$V9),as.character(dat2$V10)))
dnahit <- dnamini[HIT]
names(dnahit) <- gsub('_con.*','',names(dnahit))
write.dna(dnahit, "C1107.fasta", "fasta")


##------------------------------------------------------------------------------------------------------------
## heatmap
##------------------------------------------------------------------------------------------------------------

melted <- melt(clust_tbg[,c(1:38,41)])
melted[which(melted$strain%in%Tbg),4] <- "Tbg"; melted[which(melted$strain%in%TbgII),4] <- "TbgII"
melted[which(melted$strain%in%Tbr),4] <- "Tbr"; melted[which(melted$strain%in%Tbb),4] <- "Tbb"
melted[which(melted$strain%in%uncertain),4] <- "unc"
for (t in 1:nrow(melted)) {
  for (x in 1:8) {if (melted[t,1]%in%names(tbgC[[x]])) {melted[t,6] <- id[x]}}
}
melted[,6] <- melt(clust_tbg[,c(1:38,40)])[,1]
melted[,5] <- melt(clust_tbg[,1:39])[,1]
melted[,5] <- gsub("_len.*","",melted[,5])

colnames(melted) <- c("C", "strain", "value", "sub species", "contig","id")

groupI <- c("MHOM-ZM-80-TRPZ-23_contig10061", "MCAP-CI-91-BALEA-2_contig14961","15BT-relapse_contig1310")
groupII <- c("MBOT-GM-77-GB2_contig10794","GPAP-CI-82-KP10-29_contig344", "Jua_contig395")
groupIII <- c("GPAP-CI-82-KP10-29_contig16723", "Jua_contig1860")
groupIV <- c("MHOM-SD-82-BIYAMINA_contig456","GPAP-CI-82-KP10-29_contig17489","108AT_contig721")
groupV <- c("GPAP-CI-82-KP10-29_contig7664","NDMI_contig3299") 
groupVI <- c("108AT_contig635")

melted[which(melted$contig %in% groupI),7]  <- "I"  
melted[which(melted$contig %in% groupII),7]  <- "II" 
melted[which(melted$contig %in% groupIII),7]  <- "III"
melted[which(melted$contig %in% groupIV),7]  <- "IV"
melted[which(melted$contig %in% groupV),7]  <- "V"
melted[which(melted$contig %in% groupVI),7]  <- "VI"
colnames(melted)[7] <- "group"

library(RColorBrewer)
display.brewer.pal(6, "Set2")
display.brewer.pal(9, "YlOrRd")
brewer.pal(6, "Set2")


pdf("figures/I. Tbg MCS present in other sub species per id.pdf", useDingbats = F, width = 15, height=10)
mycolors <- c("I"="#66C2A5", "II"="#FC8D62", "III"="#8DA0CB","IV"="#E78AC3", "V"="#A6D854","VI"="#FFD92F")
g <- ggplot(melted, aes(x=strain,y=contig, fill=value, color=value)) + geom_tile() +
  theme(axis.text.x=element_text(angle=90, hjust=1), strip.text.x = element_text(size=8)) + 
  facet_grid(id ~ `sub species`, space='free', scale='free') + xlab('') + ylab('MSC') + 
  ggtitle('Presence of specific (?) Tbg clusters in other sub species') +
  theme(axis.text.y=element_text(color=mycolors))
g

mycolors <- c( "70"="#FFEDA0", "80"="#FED976", "90"="#FEB24C", 
               "95"="#FD8D3C", "96"="#FC4E2A","97"="#E31A1C", "98"="#BD0026","99"="#800026")
g2 <- ggplot(melted, aes(x=strain,y=C, fill=value, color=value)) + geom_tile(show.legend = F) +
  theme(axis.text.x=element_text(angle=90, hjust=1), strip.text.x = element_text(size=8)) + 
  facet_grid(group ~ `sub species`, space='free', scale='free') + xlab('') + ylab('MSC') + 
  ggtitle('Presence of specific (?) Tbg clusters in other sub species') +
  theme(axis.text.y=element_text(color=mycolors))

g3 <- ggplot(melted, aes(x=strain,y=contig, fill=value, color=value)) + geom_tile() +
  theme(axis.text.x=element_text(angle=90, hjust=1), strip.text.x = element_text(size=8),
        axis.text.y=element_text(size=7, angle=45)) + 
  facet_grid(group ~ `sub species`, space='free', scale='free') + xlab('') + ylab('') + 
  ggtitle('Presence of specific (?) Tbg clusters in other sub species') 

figure <- ggarrange(g2,g3)
figure

g4 <- ggplot(melted, aes(x=strain,y=contig, fill=value, color=value)) + geom_tile() +
  theme(axis.text.x=element_text(angle=90, hjust=1),strip.text.x = element_text(size=8),
        strip.text.y = element_text(size=8, angle=0)) + 
  facet_grid(id+group ~ `sub species`, space='free', scale='free') + xlab('') + ylab('MSC') + 
  ggtitle('Presence of specific (?) Tbg clusters in other sub species') +
  theme(panel.spacing=unit(0, "lines"))
g4;dev.off()


clust_tbg[which(gsub('_len.*','',clust_tbg$V39) %in% groupI),42]  <- "I"  
clust_tbg[which(gsub('_len.*','',clust_tbg$V39) %in% groupII),42]  <- "II" 
clust_tbg[which(gsub('_len.*','',clust_tbg$V39) %in% groupIII),42]  <- "III"
clust_tbg[which(gsub('_len.*','',clust_tbg$V39) %in% groupIV),42]  <- "IV"
clust_tbg[which(gsub('_len.*','',clust_tbg$V39) %in% groupV),42] <- "V"
clust_tbg[which(gsub('_len.*','',clust_tbg$V39) %in% groupVI),42]  <- "VI"


####################################################################################################
####################################################################################################

##------------------------------------------------------------------------------------------------------------
## clust mat id 97
##------------------------------------------------------------------------------------------------------------

setwd('C:/Users/mgeerts/Dropbox/Gambiense/Foefelare Manon/komics/minicircles/')
uc <- read.table('all.minicircles.circ.BGI.id97.uc')
uc_C <- uc[which(uc$V1=='C'),]
uc_C2 <- uc_C[,2]
uc_H <- uc[which(uc$V1 =='H'),]

write.csv(clust_mat_ids[[6]], "clustmat.id97.csv")
clust_mat <- clust_mat_ids[[6]] 

##------------------------------------------------------------------------------------------------------------
## clusts present in every Tbg strain id 97
##------------------------------------------------------------------------------------------------------------


list <- which(apply(clust_mat_ids[[6]][,Tbg], 1, function(x) sum(x==1)) == length(Tbg))

clust_tbg <- as.data.frame(clust_mat_ids[[6]][list,])
clust_tbg[,39] <- rownames(clust_tbg)
melted <- melt(clust_tbg)

melted[which(melted$variable%in%Tbg),4] <- "T.b. gambiense (n=18)"; melted[which(melted$variable%in%TbgII),4] <- "T.b. gambiense II (n=3)"
melted[which(melted$variable%in%Tbr),4] <- "T.b. rhodesiense (n=4)"; melted[which(melted$variable%in%Tbb),4] <- "T.b. brucei (n=7)"
melted[which(melted$variable%in%uncertain),4] <- "uncertain (n=6)"

pdf("figures/G.Presence of specific Tbg clusters in other sub specie.pdf", width=15, useDingbats = F)
g <- ggplot(melted, aes(x=variable,y=V39, fill=value, color=value)) + geom_point(show.legend = F) +
  theme(axis.text.x=element_text(angle=90, hjust=1), strip.text.x = element_text(size=8)) + 
  facet_grid(. ~ V4, space='free', scale='free') + xlab('') + ylab('MSC') + 
  ggtitle('Presence of specific (?) Tbg clusters in other sub species')
g; dev.off()


##------------------------------------------------------------------------------------------------------------
## Number of minicircles averaged per species with id of 97%
##------------------------------------------------------------------------------------------------------------

N_between <- N_between[,order(colnames(N_between))]
melted <- melt(N_between)
melted[which(melted$Var2%in%Tbg),4] <- "T.b. gambiense (n=18)"; melted[which(melted$Var2%in%TbgII),4] <- "T.b. gambiense II (n=3)"
melted[which(melted$Var2%in%Tbr),4] <- "T. b. rhodesiense (n=4)"; melted[which(melted$Var2%in%Tbb),4] <- "T.b. brucei (n=7)"
melted[which(melted$Var2%in%uncertain),4] <- "uncertain (n=6)"
colnames(melted) <- c("% identity", "strain", "Number of MSC", "sub species")
melted <- melted[order(melted$strain),]
c <- melted[melted$`% identity`=="97",]; c <- as.data.frame(c[,-1])

pdf('figures/DII.MSCav id97 in boxplot.pdf', useDingbats = F, width = 15)
g <- ggplot(c, aes(x=`sub species`, y=`Number of MSC`, color=`sub species`, fill=`sub species`)) + 
  geom_boxplot(alpha=0.5, show.legend = F) + xlab('') + ggtitle('Number of MSC per sub species') 
  #geom_text_repel(data = c$`Number of MSC`<60, label = c$strain)
g;dev.off()

##------------------------------------------------------------------------------------------------------------
## how many MSC are specific for Tbg with id 97%
##------------------------------------------------------------------------------------------------------------

setwd('C:/Users/mgeerts/Dropbox/Gambiense/Foefelare Manon/komics/')
## which MSCs are not shared?
clustTbb <-  apply(clust_mat[,Tbb], 1, function(x) sum(x!=0))/length(Tbb)
clustTbg <- apply(clust_mat[,Tbg], 1, function(x) sum(x!=0))/length(Tbg)
clustTbgII <- apply(clust_mat[,TbgII], 1, function(x) sum(x!=0))/length(TbgII)
clustTbr <- apply(clust_mat[,Tbr], 1, function(x) sum(x!=0))/length(Tbr)
clustunc <- apply(clust_mat[,uncertain], 1, function(x) sum(x!=0))/length(uncertain)
#clustTbgrest <- apply(clust_mat[,c(TbgII,Tbr,Tbb,uncertain)], 1, function(x) sum(x!=0))/length(c(TbgII,Tbr,Tbb,uncertain))

freqs <- rbind(clustTbb,clustTbg,clustTbgII,clustTbr,clustunc) 
freqslist <- apply(freqs, 2, function(x) paste(which(x != 0), collapse=''))

freqslist2 <- as.data.frame(table(freqslist)); colnames(freqslist2) <- c("x", "frequency")
#freqslist2[,1] <- c("unique for Tbg","shared", "unique for the rest")
pdf("figures/FII.MSC Tbg spec id97 per sub specie.pdf", width=15, useDingbats = F)

g <- ggplot(freqslist2, aes(x=x, y=frequency, fill=x)) + geom_bar(stat='identity', show.legend = F, alpha=0.5) + 
  geom_text(aes(label=frequency), vjust=-1,color="gray50", size=5) + xlab('') +
  annotate("text", x=2.5, y=900, label="T.b. brucei", fontface = 'italic') + 
  annotate("text", x=10, y=300, label="T.b. gambiense", fontface = 'italic') +
  annotate("text", x=12, y=450, label="T.b. gambiense II", fontface = 'italic') + 
  annotate("text", x=15, y=700, label="T.b. rhodesiense", fontface = 'italic') +
  annotate("text", x=17, y=625, label="uncertain", fontface = 'italic') + 
  ggtitle("N of unique MSC per sub specie") 
g; dev.off()

##------------------------------------------------------------------------------------------------------------
## select Tbg specific clusters id 97
##------------------------------------------------------------------------------------------------------------

## MC_select specific for Tbg
clustnames <- names(which(freqslist=='2'))  # 187

## select the ones with len<1200 and in H
setwd('C:/Users/mgeerts/Dropbox/Gambiense/Foefelare Manon/komics/minicircles/clusters Tbg')

for (n in 1:length(clustnames)) {
  dat <- uc_H[which(paste('C', uc_H$V2, sep='') == clustnames[n]),]
  datcontigs <- unique(c(as.character(dat$V9), as.character(dat$V10)))
  dna2 <- dnamini[datcontigs]
  #if(length(dna2) != 0) {dna2 <- dna2[-which(sapply(dna2, is.null))]}
  names(dna2) <- gsub('_contig.*','',names(dna2))
  if(length(dna2) != 0) {write.dna(x = dna2, file = paste(clustnames[n],'.fasta',sep=''), format = 'fasta')} 
}

## create alignments for each MSC

  
  
##------------------------------------------------------------------------------------------------------------
## heatmap
##------------------------------------------------------------------------------------------------------------

## calculate some stuff
files <- list.files(pattern ='.fasta')
clusters <- gsub('.fasta','',files)
Tbgclusters <- Tbgtemp <- list()
for (i in 1:length(clusters)) {
  Tbgtemp <- read.dna(files[i], 'fasta')
  Tbgclusters[[i]] <- Tbgtemp
}
names(Tbgclusters) <- clusters

MCincluster <- matrix(nrow = length(clusters), ncol = 1)
colnames(MCincluster) <- "clusters"
rownames(MCincluster) <- clusters
for (i in 1:length(clusters)) {
  if (length(Tbgclusters[[i]])>1) {MCincluster[i,] <- length(Tbgclusters[[i]])}
  else {MCincluster[i,] <- 1}
}
for(i in 1:length(MCincluster)) {if (MCincluster[i,]>20) {MCincluster[i,] <- 1}}

x <- clust_mat[clusters,]
Cs <- c("C1318","C1518","C1790","C2042","C2503","C3021")
melted <- melt(x)
melted[which(melted$Var2%in%Tbg),4] <- "Tbg (n=18)"; melted[which(melted$Var2%in%TbgII),4] <- "Tbg II (n=3)"
melted[which(melted$Var2%in%Tbr),4] <- "Tbr (n=4)"; melted[which(melted$Var2%in%Tbb),4] <- "Tbb (n=7)"
melted[which(melted$Var2%in%uncertain),4] <- "uncertain (n=6)"

setwd('C:/Users/mgeerts/Dropbox/Gambiense/Foefelare Manon/komics/')
pdf("figures/heatmap Tbgcluster.pdf", height=15)
g <- ggplot(melted, aes(y=Var1, x=Var2, fill=value)) + coord_equal() +
  geom_tile(show.legend = F, alpha=0.5) + xlab('') + ylab('') + scale_fill_gradientn(colours = c("white", "darkorchid")) + 
  theme(axis.text.y = element_text(size=8, hjust=1, color=ifelse(melted$Var1%in%Cs, "red", "gray")), 
        axis.text.x = element_text(size=4, angle = 45, hjust=1)) 
g
dev.off()

D <- table(apply(x[,colnames(x) %in% Tbg], 1, function(x) sum(x!=0)))


subset <- melted[melted$Var1%in%names(cl18),]
melted$value <- as.factor(melted$value)
pdf("figures/heatmap Tbgcluster subsets freq18.pdf", width=10)
g18 <- ggplot(subset, aes(y=Var1, x=Var2, color=as.factor(value))) + #, fill=value, color=value)) + #coord_equal() +
  geom_point(show.legend = F, size=3, shape=18) + xlab('') + ylab('') +#+ scale_fill_gradientn(colours = c("white", "darkorchid")) + 
  theme(axis.text.y = element_text(size=14, hjust=1,face = "bold"), 
        axis.text.x = element_text(size=8, angle = 45, hjust=1), strip.text.x = element_text(size=10)) +
  facet_grid(.~subset$V4, space='free', scale='free') + ylab('MSC') + scale_color_manual(values=c("white", "darkorchid"))
g18
dev.off()

pdf("figures/heatmap Tbgcluster subsets.pdf", width=15)
figure <- ggarrange(g18,g17,g16,g15, labels=c("Cl freq=18", "Cl freq=17", "Cl freq=16", "Cl freq=15"))
figure; dev.off()
which(melted$Var2%in%subset)
names(cl18)
dim(x)
cl18 <- which(apply(x[,colnames(x) %in% Tbg], 1, function(x) sum(x!=0))==18)
cl17 <- which(apply(x[,colnames(x) %in% Tbg], 1, function(x) sum(x!=0))==17)
cl16 <- which(apply(x[,colnames(x) %in% Tbg], 1, function(x) sum(x!=0))==16)
cl15 <- which(apply(x[,colnames(x) %in% Tbg], 1, function(x) sum(x!=0))==15)
D

for (n in 1:length(Cs)) {
  dat <- uc_H[which(paste('C', uc_H$V2, sep='') == Cs[n]),]
  datcontigs <- unique(c(as.character(dat$V9), as.character(dat$V10)))
  dna2 <- dnamini[datcontigs]
  #if(length(dna2) != 0) {dna2 <- dna2[-which(sapply(dna2, is.null))]}
  names(dna2) <- gsub('_contig.*','',names(dna2))
  names(dna2) <- paste(names(dna2), Cs[n],sep='_')
  write.dna(x = dna2, file = paste(Cs[n],'.fasta',sep=''), format = 'fasta') 
}


dnamini[Cs]


##------------------------------------------------------------------------------------------------------------
## PCA minicircles with id 97
##------------------------------------------------------------------------------------------------------------

pca <- prcomp(x = clust_mat[,c(TbgII,Tbr,Tbb,Tbg)])
pca2 <- prcomp(x = clust_mat[,c(TbgII,Tbr,Tbb,Tbg,uncertain)])

pca.r <- as.data.frame(pca["rotation"])
pca2.r <- as.data.frame(pca2["rotation"])

subset=c("LiTat-1-5-P9", "Bosendja")

for (n in 1:nrow(pca2.r)) {
  if(rownames(pca2.r)[n] %in% Tbg) {pca2.r[n,39] <- "T. b. gambiense"}
  if(rownames(pca2.r)[n] %in% TbgII) {pca2.r[n,39] <- "T. b. gambiense II"}
  if(rownames(pca2.r)[n] %in% Tbb) {pca2.r[n,39] <- "T. b. brucei"}
  if(rownames(pca2.r)[n] %in% Tbr) {pca2.r[n,39] <- "T. b. rhodesiense"}
  if(rownames(pca2.r)[n] %in% uncertain) {pca2.r[n,39] <- "uncertain"}
}
setwd('C:/Users/mgeerts/Dropbox/Gambiense/Foefelare Manon/komics/')
pdf("figures/J.PCA.pdf", useDingbats = F, width=15)
g <- ggplot(pca.r, aes(x=rotation.PC1, y=rotation.PC2, col=V33, fill=V33, shape=V33)) + geom_point(size=4, alpha=0.5) +
  stat_ellipse() + geom_text_repel(label=ifelse(rownames(pca.r)%in%subset, rownames(pca.r),''),vjust=2)
g
g <- ggplot(pca2.r, aes(x=rotation.PC1, y=rotation.PC2, col=V39, fill=V39, shape=V39)) + geom_point(size=4, alpha=0.5) +
  stat_ellipse() + geom_text_repel(label=ifelse(rownames(pca2.r)%in%subset, rownames(pca2.r),''),vjust=2)
g;dev.off()










setwd('C:/Users/mgeerts/Dropbox/Gambiense/Foefelare Manon/komics/')

dim(clust_mat)

clust_mat_melted <- melt(clust_mat) 
clust_mat_melted[which(clust_mat_melted$Var2%in%Tbg),4] <- "Tbg (n=18)"; clust_mat_melted[which(clust_mat_melted$Var2%in%TbgII),4] <- "Tbg II (n=3)"
clust_mat_melted[which(clust_mat_melted$Var2%in%Tbr),4] <- "Tbr (n=4)"; clust_mat_melted[which(clust_mat_melted$Var2%in%Tbb),4] <- "Tbb (n=7)"
clust_mat_melted[which(clust_mat_melted$Var2%in%uncertain),4] <- "uncertain (n=6)"


h <- ggplot(clust_mat_melted, aes(y=Var1, x=Var2, fill=value,color=value)) + 
  geom_tile(show.legend = F) + xlab('') +
  theme(axis.text.y = element_blank(), 
        axis.text.x = element_text(size=8, angle = 90, hjust=1), strip.text.x = element_text(size=10)) +
  facet_grid(. ~ clust_mat_melted$V4, space='free', scale='free') + ylab('MSC (n=3053)') 
h
pdf("figures/heatmap clust_mat id97.pdf", useDingbats = F, width = 15);h; dev.off()


CA <- unique(as.character(species2$region))[1]
Ca <- as.character(species2[which(species2$region==CA),1])
WA <- unique(as.character(species2$region))[2]
Wa <- as.character(species2[which(species2$region==WA),1])
SH <- unique(as.character(species2$region))[5]
Sh <- as.character(species2[which(species2$region==SH),1])

length(clustnames)
clust_mat_melted_tbg <- melt(clust_mat[clustnames,Tbg])

species2[,c(1,4)]
species2[which(species2$region==CA),1]
species2[which(species2$region==WA),1]


h <- ggplot(clust_mat_melted_tbg, aes(y=Var1, x=Var2, fill=value,color=value)) + 
  geom_tile(show.legend = F) + xlab('') + ylab('MSC (n=184)') +
  theme(axis.text.y = element_blank(), 
        axis.text.x = element_text(size=8, angle = 90, hjust=1), strip.text.x = element_text(size=10)) +
  facet_grid(. ~ clust_mat_melted_tbg$V5, space='free', scale='free')
h
qspdf("figures/heatmap clust_mat id97 per Tbg cluster.pdf", useDingbats = F, width = 15);h; dev.off()

cluster1 <- c("OUSOU","ROUPO-VAVOUA--80-MURAZ-14","LiTat-1-5-P9")
cluster2 <- c("NKOUA", "NDMI")
cluster3 <- c("MHOM-SD-82-MUSIKIA-cloneA")
cluster4 <- c("Nabe","Pakwa","Bosendja")
cluster5 <- c("MBA","LOGRA")
cluster6 <- c("108AT","108BT","57AT","348BT","15BT-relapse")
cluster7 <- "Jua"
cluster8 <- "BIM-AnTat-8-1-P8"
cluster9 <- c(cluster1, cluster2,cluster3)
clust_mat_melted_tbg[which(clust_mat_melted_tbg$Var2%in%Sh),4] <- "Sahel (n=1)" 
clust_mat_melted_tbg[which(clust_mat_melted_tbg$Var2%in%Ca),4] <- "Central Africa (n=14)"
clust_mat_melted_tbg[which(clust_mat_melted_tbg$Var2%in%Wa),4] <- "West Africa (n=3)"

clust_mat_melted_tbg[which(clust_mat_melted_tbg$Var2%in%cluster1),5] <- "cluster 1 (WA)"
clust_mat_melted_tbg[which(clust_mat_melted_tbg$Var2%in%cluster2),5] <- "cluster 2 (CB)"
clust_mat_melted_tbg[which(clust_mat_melted_tbg$Var2%in%cluster3),5] <- "cluster 3 (Sudan)"
clust_mat_melted_tbg[which(clust_mat_melted_tbg$Var2%in%cluster4),5] <- "cluster 4 (DRC)"
clust_mat_melted_tbg[which(clust_mat_melted_tbg$Var2%in%cluster5),5] <- "cluster 5 (DRC)"
clust_mat_melted_tbg[which(clust_mat_melted_tbg$Var2%in%cluster6),5] <- "cluster 6 (DRC)"
clust_mat_melted_tbg[which(clust_mat_melted_tbg$Var2%in%cluster7),5] <- "cluster 7 (DRC)"
clust_mat_melted_tbg[which(clust_mat_melted_tbg$Var2%in%cluster8),5] <- "cluster 8 (DRC)"
clust_mat_melted_tbg[which(clust_mat_melted_tbg$Var2%in%cluster9),5] <- "cluster 9 (DRC/Sudan/CB)"
