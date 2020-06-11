##------------------------------------------------------------------------------------------------------------
## FUNCTIONS
##------------------------------------------------------------------------------------------------------------
read.geno <- function(file) {
  geno <- fread(file, data.table = F, header = F)[,-1]
  rownames(geno) <- as.character(read.table(paste(file, 'indv', sep='.'))[,1])
  genopos <- read.table(paste(file, 'pos', sep='.'))
  colnames(geno) <- as.character(paste(genopos[,1], genopos[,2], sep='.'))
  return(as.data.frame(geno))
}

geno2gl <- function(geno) {
  list <- as.list(as.data.frame(t(geno)))
  loci <- colnames(geno)
  positions <- as.character(lapply(str_split(loci,'_'), function(x) x[2]))
  chromosomes <- as.character(lapply(str_split(loci,'_'), function(x) x[1]))
  gl <- new('genlight', as.list(as.data.frame(t(geno))))
  gl@chromosome <- as.factor(chromosomes)
  gl@position <- as.factor(positions)
  gl@loc.names <- loci
  return(gl)
}

rmfix <- function(X) {
  fix.ref <- which(apply(X, 2, function(x) sum(x==0)) == nrow(X))
  fix.het <- which(apply(X, 2, function(x) sum(x==1)) == nrow(X))
  fix.alt <- which(apply(X, 2, function(x) sum(x==2)) == nrow(X))
  return(X[, colnames(X)[-c(fix.ref, fix.alt, fix.het)]])
}

rmfix2 <- function(X) {
  fix.ref <- which(apply(X, 2, function(x) sum(x==0)) == nrow(X))
  fix.alt <- which(apply(X, 2, function(x) sum(x==2)) == nrow(X))
  return(X[, colnames(X)[-c(fix.ref, fix.alt)]])
}

rmfix3 <- function(X) {
  fix.ref <- which(apply(X, 2, function(x) sum(x==0)) == nrow(X))
  return(X[, colnames(X)[-fix.ref]])
}

## --------------------------------------------------------------------------------
## packages
##------------------------------------------------------------------------------------------------------------
library(data.table)
library(adegenet)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(ggplot2)
library(readxl)
library(stringr)
library(reshape2)

##------------------------------------------------------------------------------------------------------------
## subspecies
##------------------------------------------------------------------------------------------------------------
setwd('C:/Users/mgeerts/Dropbox/Gambiense/Foefelare Manon/')

#sub species <- read_excel("CHARHAT - lijst pellets DNA extractie fenol chloroform.xlsx", sheet=3)
# excel sorteert anders dan vcftools (capital sensitive)
subspecies <- read_excel("CHARHAT - lijst pellets DNA extractie fenol chloroform.xlsx", sheet=4)

##------------------------------------------------------------------------------------------------------------
## read geno
##------------------------------------------------------------------------------------------------------------
setwd('C:/Users/mgeerts/Dropbox/Gambiense/MANON')

geno <- read.geno('combined.genotyped.SNP.filtered.pASS.vcf.012') # n of SNP loci = 378.081

stats <- as.data.frame(matrix(nrow=43,ncol=5))
stats[,1] <- rownames(geno) #rownames(geno)==subspecies$strain
stats[,2] <- missing <- apply(geno, 1, function(x) sum(x==-1))
stats[,3] <- hetero <- apply(geno, 1, function(x) sum(x==1))
stats[,4] <- homo <- apply(geno, 1, function(x) sum(x==2))
stats[,5] <- subspecies$`sub species`
colnames(stats) <- c("strain", "missing", "hetero", "homo", "sub species")
melted <- reshape2::melt(stats)

setwd('C:/Users/mgeerts/Dropbox/Gambiense/Foefelare Manon/SNP_calling')
pdf("SNPstats.pdf", useDingbats = F, width = 15)
g <- ggplot(melted, aes(x=strain, y=value, fill=variable)) + xlab('') + ylab('n of SNPs') + 
  geom_bar(position="dodge", stat="identity", alpha=0.5, show.legend = F) + ggtitle("SNP stats") + 
  theme(axis.text.x=element_text(angle=90, size=8, hjust=1)) + ylim(0,130000) +
  facet_grid(variable ~ `sub species`, scale='free', space='free') 
g; dev.off()
#plot(homo, hetero, pch = 16, type = 'n')
#text(homo, hetero, labels = names(homo), cex = .5)


pdf('homo.vs.hetero.pdf', useDingbats = F, width = 15)
g <- ggplot(stats, aes(x=homo, y=hetero,fill=`sub species`, color=`sub species`)) + geom_point(size=2, alpha=0.5) +
  #geom_text_repel(data = subset(stats, hetero > 75000), label = stats$strain)
  geom_text(aes(label=ifelse(hetero>75000,as.character(strain),'')),vjust=3,hjust=0,size=3) + 
  ggtitle("plot homo SNPs versus hetero SNPs")
g; dev.off()

##------------------------------------------------------------------------------------------------------------
## pca plot 
##------------------------------------------------------------------------------------------------------------

missing.cols <- apply(geno, 2, function(x) sum(x==-1))
geno.no.missing <- geno[,names(which(missing.cols==0))] # n of SNP loci without missing data = 315.963

## all strains
genogl <- geno2gl(geno.no.missing)
genogl.pca <- glPca(genogl)
pdf('PCA all strains.pdf', useDingbats = F, width = 15)
g <- ggplot(as.data.frame(genogl.pca$scores), aes(x=PC1, y=PC2, fill=subspecies$`sub species`, color=subspecies$`sub species`)) + 
  geom_point(show.legend = F) + geom_text_repel(label=subspecies$strain, vjust=1, size=3) + 
  guides(color=guide_legend(title="sub species")) + ggtitle('PCA all strains (sub species A)')
g

#plot(genogl.pca$scores[,1], genogl.pca$scores[,2], pch = 16, cex = 2, xlab = 'PC1', ylab = 'PC2')
#text(genogl.pca$scores[,1], genogl.pca$scores[,2], labels = rownames(genogl.pca$scores), cex = .5)

##------------------------------------------------------------------------------------------------------------
## subspecies2
##------------------------------------------------------------------------------------------------------------

### subsets based on pCA and splitstree
subspecies2 <- subspecies
subspecies2[subspecies2$strain=="MSUS-CI-78-TSW-157",2] <- "T.b. gambiense II"
#sub species2[subspecies2$strain==names(which(genogl.pca$scores[,1]<(-50))),2]
subspecies2[subspecies2$strain=="MHOM-SD-82-MUSIKIA-cloneA",2] <- "T.b. gambiense"
tmp.Tbb <- names(which(genogl.pca$scores[,2]>50)) 
Tbb <- tmp.Tbb[c(-1,-5,-6,-9)] # 10
subspecies2[which(subspecies2$strain%in%Tbb),2] <- "T.b. brucei"
subspecies2[which(subspecies2$strain%in%names(which(genogl.pca$scores[,1]>150))),2] <- "T.b. rhodesiense"
uncertain <- c("STIB-851", "FEO", "FEO-AnTat-16-1", "MHOM-SD-82-BIYAMINA", "MBO-NG-74-R10","GMOM-ZM-83-TRPZ-317") # 6
subspecies2[which(subspecies2$strain%in%uncertain),2] <- "uncertain"

table(subspecies$`sub species`); table(subspecies2$`sub species`)

g <- ggplot(as.data.frame(genogl.pca$scores), aes(x=PC1, y=PC2, fill=subspecies2$`sub species`, color=subspecies2$`sub species`)) + 
  geom_point(show.legend = F) + geom_text_repel(label=subspecies$strain, vjust=0, size=3) + 
  guides(color=guide_legend(title="sub species")) + ggtitle('PCA all strains (sub species B)')
g; dev.off()

write.csv(subspecies2, "list sub species.csv")

Tbb <- unlist(as.list(subspecies2[which(subspecies2$`sub species`=="T.b. brucei"),1]))
Tbg <- unlist(as.list(subspecies2[which(subspecies2$`sub species`=="T.b. gambiense"),1]))
TbgII <- unlist(as.list(subspecies2[which(subspecies2$`sub species`=="T.b. gambiense II"),1]))
Tbr <- unlist(as.list(subspecies2[which(subspecies2$`sub species` =="T.b. rhodesiense"),1]))
uncertain <- unlist(as.list(subspecies2[which(subspecies2$`sub species`=="uncertain"),1]))

##------------------------------------------------------------------------------------------------------------
## pca plot per group
##------------------------------------------------------------------------------------------------------------

## subset 1: Tbb + TbgII + MBO
subset1 <- names(which(genogl.pca$scores[,2]>50))
geno.subset1 <- geno.no.missing[subset1,]
genogl.subset1 <- geno2gl(geno.subset1)
genogl.pca.subset1 <- glPca(genogl.subset1)
subspecies_subset1 <- subspecies2[subspecies2$strain%in%subset1,]
pdf('pCA.subset1.pdf', useDingbats = F, width = 15)
q <- ggplot(as.data.frame(genogl.pca.subset1$scores), aes(x=PC1, y=PC2, color=subspecies_subset1$`sub species`)) + 
  geom_point(size=2, alpha=0.5) + geom_text_repel(label=rownames(genogl.pca.subset1$scores), vjust=1, size=4) + 
  ggtitle("PCA - Subset 1") + theme(legend.title = element_blank())
  #color = ifelse(subset1%in%Tbb, "grey50", ifelse(subset1%in%TbgII, "forestgreen", "red"))) 
q ; dev.off()

## subset 2: Tbg minus FEO strains
subset2 <- names(which(genogl.pca$scores[,1]<(-50)))
geno.subset2 <- geno.no.missing[subset2,]
genogl.subset2 <- geno2gl(geno.subset2)
genogl.pca.subset2 <- glPca(genogl.subset2)
subspecies_subset2 <- subspecies2[subspecies2$strain%in%subset2,]
pdf('pCA.subset2.Tbg.pdf', useDingbats = F, width = 15)
q <- ggplot(as.data.frame(genogl.pca.subset2$scores), aes(x=PC1, y=PC2, color=subspecies_subset2$host)) + 
  geom_point(size=2, alpha=0.5, show.legend = F) + geom_text_repel(label=rownames(genogl.pca.subset2$scores), vjust=1, size=4) +  
  #color=sub species2[which(sub species2$strain%in%subset2),]$country) + #color = ifelse(subset2%in%Tbg, "grey50", "red")) + 
  ggtitle("PCA - Subset 2 (Tbg) by host") + guides(color=guide_legend(title="host"))
q
q2 <- ggplot(as.data.frame(genogl.pca.subset2$scores), aes(x=PC1, y=PC2, color=subspecies_subset2$country)) + 
  geom_point(size=2, alpha=0.5, show.legend = F) + geom_text_repel(label=rownames(genogl.pca.subset2$scores), vjust=1, size=4) +  
  ggtitle("PCA - Subset 2 (Tbg) by country") + guides(color=guide_legend(title="country"))
q2

#glMDSplot(genogl.pca.subset2$scores, groups=Tbg$sub species, labels = rownames(genogl.pca.subset2$scores))

q3 <- ggplot(as.data.frame(genogl.pca.subset2$scores), aes(x=PC1, y=PC2, color=subspecies_subset2$region)) + 
  geom_point(size=2, alpha=0.5, show.legend = F) + geom_text_repel(label=rownames(genogl.pca.subset2$scores), vjust=1, size=4) +  
  ggtitle("PCA - Subset 2 (Tbg) by region") + guides(color=guide_legend(title="region"))
q3; dev.off()

#plot(genogl.pca.tbg$scores[,1], genogl.pca.tbg$scores[,2], pch = 16, cex = 2, xlab = 'PC1', ylab = 'PC2')
#text(genogl.pca.tbg$scores[,1], genogl.pca.tbg$scores[,2], labels = rownames(genogl.pca.tbg$scores), cex = .5)

## subset 3: Tbg minus strains from Cameroon
subset3 <- names(which(genogl.pca$scores[,1]<(-50)))
subset3 <- subset3[c(-7,-9,-10)]
geno.subset3 <- geno.no.missing[subset3,]
genogl.subset3 <- geno2gl(geno.subset3)
genogl.pca.subset3 <- glPca(genogl.subset3)
subspecies_subset3 <- subspecies2[subspecies2$strain%in%subset3,]
pdf('pCA.subset3.pdf', useDingbats = F, width = 15)
q <- ggplot(as.data.frame(genogl.pca.subset3$scores), aes(x=PC1, y=PC2, color=subspecies_subset3$country)) + 
  geom_point(size=2, alpha=0.5, show.legend = F) + geom_text_repel(label=rownames(genogl.pca.subset3$scores), vjust=1, size=4) + 
  ggtitle("PCA - Subset 3") + guides(color=guide_legend(title="country"))
q 
q2 <- ggplot(as.data.frame(genogl.pca.subset3$scores), aes(x=PC1, y=PC2, color=subspecies_subset3$region)) + 
  geom_point(size=2, alpha=0.5, show.legend = F) + geom_text_repel(label=rownames(genogl.pca.subset3$scores), vjust=1, size=4) + 
  ggtitle("PCA - Subset 3") + guides(color=guide_legend(title="region"))
q2 ; dev.off()

##------------------------------------------------------------------------------------------------------------
## overal SNP count
##------------------------------------------------------------------------------------------------------------

### n SNPs tegenover ref genome vs variation between sub groups 
## rmfix removes all SNPs where variation = 0
N_SNP <- as.data.frame(matrix(ncol=6))
colnames(N_SNP) <- c(names(table(subspecies2$`sub species`)),"Total")
N_SNP[1,1:5] <- table(subspecies2$`sub species`)
N_SNP[1,6] <- sum(N_SNP[1,1:5])
N_SNP[3,1] <- as.numeric(ncol(rmfix(geno.no.missing[Tbb,])))
N_SNP[3,2] <- as.numeric(ncol(rmfix(geno.no.missing[Tbg,])))
N_SNP[3,3] <- as.numeric(ncol(rmfix(geno.no.missing[TbgII,])))
N_SNP[3,4] <- as.numeric(ncol(rmfix(geno.no.missing[Tbr,])))
N_SNP[3,5] <- as.numeric(ncol(rmfix(geno.no.missing[uncertain,])))
N_SNP[3,6] <- as.numeric(ncol(rmfix(geno.no.missing)))

#rmfix3(geno.no.missing[Tbb,]) zou hetzelfde moeten zijn als derde rij N_SNP
SNPcount <- t(read.table("SNP.count.txt"))
N_SNP[2,] <- SNPcount[2,]
rownames(N_SNP) <- c("n of strains", "n of SNP loci vs Tb927","n of SNP loci within each sub species")

write.csv(N_SNP, 'N_SNP.csv')

tN_SNP <- as.data.frame(t(N_SNP[,1:5]))

pdf("variation vs SNP.pdf", useDingbats = F, width = 10)
g <- ggplot(tN_SNP, aes(`n of SNP loci within each sub species`,`n of SNP loci vs Tb927`)) + 
  geom_point(alpha=0.5,size=3) + geom_text_repel(label=rownames(tN_SNP), vjust=1, size=4)
g; dev.off()

##------------------------------------------------------------------------------------------------------------
## BARPLOT: distribution of SNPs with a 150 SNP window per sub species
##------------------------------------------------------------------------------------------------------------

sp <- c("Tbb", "Tbg", "TbgII", "Tbr","uncertain")
tmp.count <- as.data.frame(matrix(nrow=5, ncol=150)); rownames(tmp.count) <- sp
breaks <- seq(1,ncol(geno.no.missing), by=150)
N_SNP.150w <- matrix(ncol=5, nrow=length(breaks)); colnames(N_SNP.150w) <- sp
for (i in 1:length(breaks)) {
  colnames(tmp.count) <- colnames(geno.no.missing[,breaks[i]:(breaks[i]+149)])
  tmp.count[1,] <- apply(geno.no.missing[Tbb,breaks[i]:(breaks[i]+149)], 2, function(x) sum(x!=0))/length(Tbb)
  tmp.count[2,] <- apply(geno.no.missing[Tbg,breaks[i]:(breaks[i]+149)], 2, function(x) sum(x!=0))/length(Tbg)
  tmp.count[3,] <- apply(geno.no.missing[TbgII,breaks[i]:(breaks[i]+149)], 2, function(x) sum(x!=0))/length(TbgII)
  tmp.count[4,] <- apply(geno.no.missing[Tbr,breaks[i]:(breaks[i]+149)], 2, function(x) sum(x!=0))/length(Tbr)
  tmp.count[5,] <- apply(geno.no.missing[uncertain,breaks[i]:(breaks[i]+149)], 2, function(x) sum(x!=0))/length(uncertain)
  
  N_SNP.150w[i,1:5] <- rowSums(tmp.count)
}

N_SNP.150w <- N_SNP.150w[-2107,]
melted <- reshape2::melt(N_SNP.150w)

pdf("BARPLOT - distribution of SNPs with a 150 SNP window per sub species.pdf", useDingbats = F, width = 15)
g <- ggplot(melted, aes(x=Var1,y=value)) + geom_bar(stat='identity', alpha=0.5, color="darkblue") + 
  facet_grid(Var2 ~ .) + theme(axis.text.x = element_blank()) + #theme(axis.text.x = element_text(angle=90))  
  xlab('SNP loci') + ylab('n of SNPs per 150 SNP loci') + ggtitle("SNP distribution per sub species") + ylim(0,150)
g; dev.off()


##------------------------------------------------------------------------------------------------------------
## BARPLOT: distribution of SNP variation with a 150 SNP window between sub groups
##------------------------------------------------------------------------------------------------------------

tmp.count <- as.data.frame(matrix(nrow=5, ncol=150)); rownames(tmp.count) <- sp
N_SNPvar.150w <- matrix(ncol=5, nrow=length(breaks)); colnames(N_SNPvar.150w) <- sp
for (i in 1:length(breaks)) {
  colnames(tmp.count) <- colnames(geno.no.missing[,breaks[i]:(breaks[i]+149)])
  tmp.count["Tbb",] <- as.numeric(colnames(tmp.count) %in% colnames(rmfix(geno.no.missing[Tbb,breaks[i]:(breaks[i]+149)])))                 
  tmp.count["Tbg",] <- as.numeric(colnames(tmp.count) %in% colnames(rmfix(geno.no.missing[Tbg,breaks[i]:(breaks[i]+149)]))) 
  tmp.count["TbgII",] <- as.numeric(colnames(tmp.count) %in% colnames(rmfix(geno.no.missing[TbgII,breaks[i]:(breaks[i]+149)]))) 
  tmp.count["Tbr",] <- as.numeric(colnames(tmp.count) %in% colnames(rmfix(geno.no.missing[Tbr,breaks[i]:(breaks[i]+149)]))) 
  tmp.count["uncertain",] <- as.numeric(colnames(tmp.count) %in% colnames(rmfix(geno.no.missing[uncertain,breaks[i]:(breaks[i]+149)]))) 
  
  N_SNPvar.150w[i,1:5] <- rowSums(tmp.count)
}  

N_SNPvar.150w <- N_SNPvar.150w[-2107,]
melted <- reshape2::melt(N_SNPvar.150w)

pdf("BARPLOT - distribution of SNP variation with a 150 SNP window per sub species.pdf", useDingbats = F, width = 15)
g <- ggplot(melted, aes(x = Var1, y=value)) + geom_bar(stat = "identity", alpha=0.5, color="darkblue") + 
  facet_grid(Var2 ~ .) + theme(axis.text.x = element_blank(), legend.title=element_blank()) + ylim(0,150) +
  xlab('SNP variation') + ylab('count per 150 SNP loci') + ggtitle("SNP variation between sub species")
g; dev.off()


##------------------------------------------------------------------------------------------------------------
## BARPLOT: count SNP variation per chromosome
##------------------------------------------------------------------------------------------------------------

# \. matches a dot
# [^\.] matches everything but a dot
# * specifies that the previous expression (everything but a dot) may occur between 0 and unlimited times
# $ marks the end of the string
chr <- unique(gsub("\\.[^\\.]*$",'',colnames(geno.no.missing))) # 45
chr_subset <- chr[1:11]
#chromosomes <- as.list(read.table("chr.Tb927.txt")) # 127

tmp.count <- as.data.frame(matrix(nrow=5, ncol=150)); rownames(tmp.count) <- sp
window150 <- matrix(ncol=5, nrow=length(breaks)); colnames(window150) <- sp
chr <- list()
for (i in 1:length(breaks)) {
  colnames(tmp.count) <- colnames(geno.no.missing[,breaks[i]:(breaks[i]+149)])
  tmp.count["Tbb",] <- as.numeric(colnames(tmp.count) %in% colnames(rmfix(geno.no.missing[Tbb,breaks[i]:(breaks[i]+149)])))                 
  tmp.count["Tbg",] <- as.numeric(colnames(tmp.count) %in% colnames(rmfix(geno.no.missing[Tbg,breaks[i]:(breaks[i]+149)]))) 
  tmp.count["TbgII",] <- as.numeric(colnames(tmp.count) %in% colnames(rmfix(geno.no.missing[TbgII,breaks[i]:(breaks[i]+149)]))) 
  tmp.count["Tbr",] <- as.numeric(colnames(tmp.count) %in% colnames(rmfix(geno.no.missing[Tbr,breaks[i]:(breaks[i]+149)]))) 
  tmp.count["uncertain",] <- as.numeric(colnames(tmp.count) %in% colnames(rmfix(geno.no.missing[uncertain,breaks[i]:(breaks[i]+149)]))) 
  
  window150[i,1:5] <- rowSums(tmp.count)
  chr[i] <- unique(gsub("\\.[^\\.]*$",'',colnames(tmp.count)))[1]
}  

window150 <- window150[-2107,]

melted2 <- reshape2::melt(window150)
melted2[,4] <- unlist(chr)


g <- ggplot(melted2, aes(x = Var1, y=value)) + geom_bar(stat = "identity", alpha=0.5, color="darkblue") + 
  facet_grid(Var2 ~ V4, space='free', scale='free') + 
  theme(axis.text.x = element_blank(), legend.title=element_blank(), strip.text.x = element_text(angle = 90)) + 
  xlab('SNP variation') + ylab('count per 150 SNP loci') + ggtitle("SNP variation between sub sub species") + ylim(0,150)
g
pdf("HUPLLAAAAAAA BARPLOT - distribution of SNP variation with a 150 SNP window per sub species EN PER CHR.pdf", useDingbats = F, width = 15)
g2 <- ggplot(melted2[which(melted2$V4%in%chr_subset),], aes(x = Var1, y=value)) + geom_bar(stat = "identity", alpha=0.5, color="darkblue") + 
  facet_grid(Var2 ~ V4, space='free', scale='free') + theme(axis.text.x = element_blank(), legend.title=element_blank(), strip.text.x = element_text(angle = 90)) + 
  xlab('SNP variation') + ylab('count per 150 SNP loci') + ggtitle("SNP variation between sub sub species") + ylim(0,150)
g2; dev.off()




##------------------------------------------------------------------------------------------------------------
## rommel
##------------------------------------------------------------------------------------------------------------


ncol(geno.no.missing)/10000

x <- as.data.frame(matrix(nrow=5, ncol=10000))
rownames(x) <- c("Tbg", "Tbb", "TbgII", "Tbr","uncertain")
breaks <- seq(1,ncol(geno.no.missing), by=10000)
bla <- matrix(nrow=1,ncol=5)
for (i in 1:length(breaks)) {
  colnames(x) <- colnames(geno.no.missing[,breaks[i]:(breaks[i]+9999)])
  x[1,] <- apply(geno.no.missing[Tbg,breaks[i]:(breaks[i]+9999)], 2, function(x) sum(x!=0))/length(Tbg)
  x[2,] <- apply(geno.no.missing[Tbb,breaks[i]:(breaks[i]+9999)], 2, function(x) sum(x!=0))/length(Tbb)
  x[3,] <- apply(geno.no.missing[TbgII,breaks[i]:(breaks[i]+9999)], 2, function(x) sum(x!=0))/length(TbgII)
  x[4,] <- apply(geno.no.missing[Tbr,breaks[i]:(breaks[i]+9999)], 2, function(x) sum(x!=0))/length(Tbr)
  x[5,] <- apply(geno.no.missing[uncertain,breaks[i]:(breaks[i]+9999)], 2, function(x) sum(x!=0))/length(uncertain)
  
  bla <- rbind(bla, rowSums(x))
}

bla <- bla[-1,]
melted <- reshape2::melt(bla)
pdf("10k window uitprobeersel.pdf", useDingbats = F, width = 15)
g <- ggplot(melted, aes(x=Var1,y=value, color=Var2, fill=Var2)) + geom_bar(stat='identity', alpha=0.5) + 
  facet_grid(Var2 ~ .) + theme(axis.text.x = element_blank()) #theme(axis.text.x = element_text(angle=90))
g; dev.off()



