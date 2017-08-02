
#Project Drug Repositioning
#Renal Analysis
#Joshua Ellis
#CSC7300 Bioinformatics 1
source("https://bioconductor.org/biocLite.R")

############
#GSE1563
#Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE1563", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL8300", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "11111111000000111111111XXXXXXXX11111111100000001111111111XXXXX"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT1 <- topTable(fit2, adjust="fdr", sort.by="B", number=5000)

tT1 <- subset(tT1, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT1, file=stdout(), row.names=F, sep="\t")

save(tT1, file="tT1.Rdata")

############
#GSE21374
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE21374", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0("11111001110111000111111111111100110001011100111111",
               "11110000101111111000111101111101110100001110101110",
               "11001110111001111101110111111111110111111000111110",
               "11111111111111110111101001111111111010001011011110",
               "01111111101111111111111111111110101110111011111110",
               "01100011101000110110110100110110")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=5000)

tT2 <- subset(tT2, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT2, file=stdout(), row.names=F, sep="\t")
save(tT2, file="tT2.Rdata")
#############
#GSE36059

library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE36059", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0("11111111010100011111111100011011001011110001000011",
               "11110000110111110010100111111111000010110110101101",
               "00011110100110011111111111111000111111111111111100",
               "11111111111111110111011111111111111110111010111111",
               "10100111111111110111011100011111001111011110101110",
               "11111011111101111111111011111111111111111111110100",
               "10011111011110101011111110101110111111111111111111",
               "11111111111011111110111111111111111110111101111001",
               "10011110011")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT3 <- topTable(fit2, adjust="fdr", sort.by="B", number=5000)

tT3 <- subset(tT3, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT3, file=stdout(), row.names=F, sep="\t")
save(tT3, file="tT3.Rdata")
############
#GSE25092


library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE25092", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6246", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "000111111"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT4 <- topTable(fit2, adjust="fdr", sort.by="B", number=5000)

tT4 <- subset(tT4, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT4, file=stdout(), row.names=F, sep="\t")
save(tT4, file="tT4.Rdata")
##############
#GSE50058
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE50058", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0("00000000000000000000001111111100011111110111111010",
               "11111100111100000111101110000000000000010011001000",
               "0")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=1000)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")
save(tT,file="tT.Rdata")
##########
# load series and platform data from GEO

gset <- getGEO("GSE50058", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0("00000000000000000000001111111100011111110111111010",
               "11111100111100000111101110000000000000010011001000",
               "0")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=1000)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")
save(tT, file="tT.Rdata")
#########
#GROUP and Fisher's
Z=rbind(tT,tT1,tT2,tT3,tT4)
save(Z,file="Z.Rdata")
View(Z)


biocLite("metap")
library(metap)
biocLite("dplyr")
library(dplyr)
Zp=Z %>% group_by(Gene.symbol) %>% filter(n()>4)

Zp=Zp[order(Zp$Gene.symbol),]
View(Zp)
L=length(Zp$ID)
#EFFECT SIZE
for (i in 1:L){
  if (abs(Zp$logFC[i])<=1){
    Zp$Gene.symbol[i]=""
    
  }
}

Zp=Zp[!(is.na(Zp$Gene.symbol) | Zp$Gene.symbol==""), ]
save(Zp,file="Zp.Rdata")

L=length(Zp$Gene.symbol)

d=0
gname=0
j=0
e=1
C=0
for (i in 1:L){
  if (Zp$Gene.symbol[i]==Zp$Gene.symbol[i+1]){
    counts=sum(Zp$Gene.symbol==Zp$Gene.symbol[i])
    e=i
    for (k in 1:counts){
      C[k]=c(Zp$adj.P.Val[e])
      e=e+1
   
    }
    Fis=sumlog(C)
    d[i]=Fis$p
    gname[i]=as.character(Zp$Gene.symbol[i])
    if (as.numeric(as.character(d[i]))>0.002){
      gname[i]="NA"
    }
  }

  
  else{
    d[i]="NA"
    gname[i]="NA"
  }
  j=j+1
  
}
require(reshape2)
d= melt(data.frame(d,gname))
colnames(d)= c("P","Gene.Symbol")

#Remove NA
d=d[!(is.na(d$Gene.Symbol) | d$Gene.Symbol=="NA"), ]
dr=d[!duplicated(d$Gene.Symbol),]
save(dr,file="dr.Rdata")
View(dr)


##########
#Pathway
biocLite("hgu133a2.db")
library(hgu133a2.db)
ls("package:hgu133a2.db")
"org.Hs.egGENENAME"
## Bimap interface:
x <- hgu133a2ENTREZID
# Get the gene names that are mapped to an entrez gene identifier
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
l=length(dr$Gene.Symbol)
if(length(xx) > 0) {
  # Get the Ensembl gene IDs for the first five genes
  E=xx[1:l]
  # Get the first one
  #xx[[1]]
  
}
dr$ENTREZ=E
View(dr)
View(E)
