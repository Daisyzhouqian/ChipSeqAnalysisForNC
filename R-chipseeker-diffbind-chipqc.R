#source ("https://bioconductor.org/biocLite.R")
#biocLite("ChIPseeker")
#biocLite("TxDb.Mmusculus.UCSC.mm10.knownGene")
#biocLite("ReactomePA")
#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("ChIPseeker")
#BiocManager::install("org.Mm.eg.db")
#BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
#BiocManager::install("clusterProfiler")
#BiocManager::install("ReactomePA")
library("ChIPseeker")
library("org.Mm.eg.db")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
library("clusterProfiler")
library(DOSE)
library("ReactomePA")
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

#setwd("F:/Data/chip_seq/aligned")
#suz12<-readPeakFile("M2EZH2-200-12-5.filter.bed")
#查看peak在全基因组的位置
covplot(suz12,weightCol=5)
#covplot(peak.gr,weightCol=5)
#自定义染色体位置
covplot(suz12, weightCol=5, chrs=c("chr4", "chr5"), xlim=c(4.5e7, 5e7))

#ChIP peaks结合TSS 区域的情况
#TSS:transcription start site
#首先，计算ChIP peaks结合在TSS区域的情况。这就需要准备TSS区域，这一般定义在TSS位点的侧翼序列（默认-3000~+3000）。然后比对 map到这些区域的peak，并生成tagMatrix
filelist=list.files()
filelist
diffpeak1<-readPeakFile(filelist[15])
diffpeak2<-readPeakFile(filelist[16])
diffpeak1$"Diff"="decrease"
diffpeak2$"Diff"="increase"
gene1 <- seq2gene(diffpeak1, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
gene2 <- seq2gene(diffpeak2, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
gene<-c(gene1,gene2)

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(diffpeak, windows=promoter)
#Heatmap of ChIP binding to TSS regions
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
#peakHeatmap(suz12, TxDb=txdb, upstream=3000, downstream=3000, color="blue") 一步生成
#Average Profile of ChIP peaks binding to TSS region
plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
            conf=0.95,resample = 1000,
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

peakAnno <- annotatePeak(diffpeak1, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Mm.eg.db")
#peakAnno <- annotatePeak(diffpeak2, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Mm.eg.db")
#可视化基因组注释
plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)#注释会重叠
upsetplot(peakAnno) #重叠
upsetplot(peakAnno, vennpie=TRUE)

write.csv(peakAnno %>% data.frame(),file = "../ZygoteGenome1-M2WT-control-increase.csv")

library(ReactomePA)

#gene <- seq2gene(diffpeak, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
pathway2 <- enrichPathway(gene,pvalueCutoff = 0.2,organism = "mouse")
dotplot(pathway2,showCategory=15)
#write.csv(pathway2,"../Zygote-Con-genome12-diff-pathway.csv")
write.csv(pathway2,"Zygote-EZH2-Con-genome2-D$I-reactome.csv")


ego <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = 'mmu',
  pvalueCutoff      = 0.2,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.5
)
write.csv(ego,"Zygote-EZH2-Con-genome2-D$I-kegg.csv")
dotplot(ego,showCategory=15)

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                ont = "BP", #BP,CC,MF
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 1,
                readable = TRUE)
write.csv(ego,"Zygote-EZH2-Con-genome2-D$I-GO-BP.csv")
write.csv(ego,"Zygote-EZH2-Con-genome2-D$I-GO-CC.csv")
write.csv(ego,"Zygote-EZH2-Con-genome2-D$I-GO-MF.csv")
dotplot(ego,showCategory=15)


symbol=as.character(gene[,1])
eg = bitr(symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
id = as.character(eg[,2])
ids <- bitr(symbol, fromType="SYMBOL", toType=c("UNIPROT", "ENSEMBL"), OrgDb="org.Hs.eg.db")
head(id)


##################### multiple peaks ###################
#读取进去后重新命名
filelist=list.files()
H3K27me2_WT_1<-readPeakFile(filelist[8])
H3K27me2_WT_2<-readPeakFile(filelist[1])
H3K27me2_EZH2<-readPeakFile(filelist[5])
H3K27me3_WT<-readPeakFile(filelist[2])
#H3K27me3_WT_1<-readPeakFile(filelist[3])
#H3K36me3_WT<-readPeakFile(filelist[4])

H3K27me2_ZygoteWT_Maternal<-readPeakFile(filelist[10])
H3K27me2_ZygoteWT_Paternal<-readPeakFile(filelist[11])
H3K27me2_ZygoteEZH2_Maternal<-readPeakFile(filelist[18])
H3K27me2_ZygoteEZH2_Paternal<-readPeakFile(filelist[19])

#promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
#tagMatrix <- getTagMatrix(H3K27me2_ZygoteEZH2_Paternal, windows=promoter)
#Heatmap of ChIP binding to TSS regions
#tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
#plotAvgProf(tagMatrix, xlim=c(-3000, 3000),conf=0.95,resample = 1000,xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")


peaks<-list(H3K27me2_WT_1=H3K27me2_WT_1,H3K27me2_WT_2=H3K27me2_WT_2,
            H3K27me2_EZH2=H3K27me2_EZH2,H3K27me3_WT=H3K27me3_WT)

peaks<-list(H3K27me2_ZygoteWT_Maternal=H3K27me2_ZygoteWT_Maternal,H3K27me2_ZygoteWT_Paternal=H3K27me2_ZygoteWT_Paternal,
            H3K27me2_ZygoteEZH2_Maternal=H3K27me2_ZygoteEZH2_Maternal,H3K27me2_ZygoteEZH2_Paternal=H3K27me2_ZygoteEZH2_Paternal)


#peaks <- list(cbx7=cbx7,ring1B=ring1B,suz12=suz12)
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(peaks, getTagMatrix, windows=promoter)

#plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=c("blue","red2","seagreen","purple"))

peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,tssRegion=c(-3000, 3000), verbose=FALSE)
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)

#upsetplot(peakAnnoList[[1]], vennpie=TRUE)
plotAnnoPie(peakAnnoList[[1]])
plotAnnoPie(peakAnnoList[[2]])
plotAnnoPie(peakAnnoList[[3]])
plotAnnoPie(peakAnnoList[[4]])

genes= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
vennplot(genes)

write.csv(peakAnnoList[[1]] %>% data.frame(),file = "../H3K27me2_ZygoteWT_Maternal.csv")
write.csv(peakAnnoList[[2]] %>% data.frame(),file = "../H3K27me2_ZygoteWT_Paternal.csv")
write.csv(peakAnnoList[[3]] %>% data.frame(),file = "../H3K27me2_ZygoteEZH2_Maternal.csv")
write.csv(peakAnnoList[[4]] %>% data.frame(),file = "../H3K27me2_ZygoteEZH2_Paternal.csv")

genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
#names(genes) = sub("_", "\n", names(genes))
compKEGG <- compareCluster(geneCluster   = genes,organism="mmu",
                           fun = "enrichKEGG", #enrichGO OrgDb='org.Mm.eg.db',
                           pvalueCutoff  = 0.05, 
                           pAdjustMethod = "BH")
dotplot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")

compKEGG <- compareCluster(geneCluster   = genes,OrgDb='org.Mm.eg.db',
                           fun = "enrichGO", #enrichGO 
                           pvalueCutoff  = 0.05, 
                           pAdjustMethod = "BH")
dotplot(compKEGG, showCategory = 15, title = "GO Enrichment Analysis")







########################## CHIPQC 包学习 ############################
#source("http://bioconductor.org/biocLite.R")
#biocLite("ChIPQC")
## Load libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("ChIPQC")
library("DiffBind")
library(ChIPQC)
## Load sample data
samples <- read.csv('meta/samplesheet_chr12.csv')
View(samples)
## Create ChIPQC object
chipObj <- ChIPQC(samples, annotation="hg19") 
## Create ChIPQC report
ChIPQCreport(chipObj, reportName="ChIP QC report: Nanog and Pou5f1", reportFolder="ChIPQCreport")

#### example 1
##Experiment sample sheet
samples = read.csv(file.path(system.file("extdata", package="ChIPQC"),"example_QCexperiment.csv"))
exampleExp = ChIPQC(samples,annotaiton="hg19")
exampleExp

data(example_QCexperiment)
exampleExp
ChIPQCreport(exampleExp)

#### example 2
library(ChIPQC)
samples = read.csv(file.path(system.file("extdata", package="ChIPQC"),"tamoxifenQC.csv"))
data(blacklist_hg19)
tamoxifen = ChIPQC(samples, consensus=TRUE, bCount=TRUE, summits=250,
                   annotation="hg19", chromosomes="chr18",
                   blacklist = blacklist.hg19)
tamoxifen
data(tamoxifen_QC)
tamoxifen
ChIPQCreport(tamoxifen,facetBy=c("Tissue","Condition"))
#Plotting Coverage Histograms
plotCoverageHist(tamoxifen,facetBy=c("Tissue","Condition"))
#Plotting Cross-Coverage
plotCC(tamoxifen,facetBy=c("Tissue","Condition"))

#Plotting Relative Enrichment of reads in Genomic Intervals
plotRegi(tamoxifen,facetBy=c("Tissue","Condition"))

#Plotting Peak Profiles
plotPeakProfile(tamoxifen,facetBy=c("Tissue","Condition"))

#Plotting Reads overlapping Peaks and the Blacklist
plotRap(tamoxifen,facetBy=c("Tissue","Condition"))
plotFribl(tamoxifen,facetBy=c("Tissue","Condition"))

#Plotting Sample Clustering
plotCorHeatmap(tamoxifen,attributes=c("Tissue","Factor","Condition","Replicate"))

#####Single sample assessment
CTCF1 = ChIPQCsample("reads/SRR568129.bam", peaks="peaks/SRR568129_chr22_peaks.bed")
CTCF1 = QCsample(exampleExp,1)
QCmetrics(CTCF1)
plotCC(CTCF1)

########################## DiffBind ###############################
#browseVignettes("DiffBind")
library(DiffBind)
##过程如下
tamoxifen <- dba(sampleSheet="tamoxifen.csv")
tamoxifen <- dba.count(tamoxifen)
tamoxifen <- dba.contrast(tamoxifen)
tamoxifen <- dba.analyze(tamoxifen)
tamoxifen.DB <- dba.report(tamoxifen)
##example 1
samples <- read.csv(file.path(system.file("extra", package="DiffBind"),"tamoxifen.csv"))
samples

tamoxifen <- dba(sampleSheet="tamoxifen.csv",dir=system.file("extra", package="DiffBind"))
tamoxifen

tamoxifen <- dba.count(tamoxifen, summits=250)
tamoxifen <- dba.contrast(tamoxifen, categories=DBA_CONDITION)
tamoxifen <- dba.analyze(tamoxifen)
plot(tamoxifen, contrast=1)
tamoxifen.DB <- dba.report(tamoxifen)

##########################
filelist<-list.files()

M2EZH2Chip<-read.table("M2EZH2-ChIP-genome1.depth.txt",sep='\t',header = FALSE)
M2EZH2Input<-read.table("M2EZH2-input-genome1.depth.txt",sep='\t',header = FALSE)
M2EZH2_Depth<-M2EZH2Chip$V7-M2EZH2Input$V7

M2EZH2_Depth[M2EZH2_Depth<=0]=NA


M2WTChip<-read.table(filelist[3],sep='\t',header = FALSE)
M2WTInput<-read.table(filelist[4],sep='\t',header = FALSE)
M2WT_Depth<-M2WTChip$V7-M2WTInput$V7

M2WT_Depth[M2WT_Depth<=0]=NA

data<-data.frame(Depth=c(M2EZH2_Depth,M2WT_Depth),Group=rep(c("M2EZH2","M2WT"),c(length(M2EZH2_Depth),length(M2WT_Depth))))
data<-na.omit(data)

library(ggpubr)

my_comparisons<-c("M2EZH2","M2WT")
ggviolin(data,x="Group",y="Depth",fill="Group",palette = c("#00AFBB","#FC4E07"),add = "boxplot",
         add.params = list(fill="white"))+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  stat_compare_means(label.y = 1)

########################

filelist<-list.files()
#filelist<-filelist[grep("depth",filelist)]

M2EZH2Chip<-read.table("M2EZH2-ChIP-genome1.depth.txt",sep='\t',header = FALSE)
M2EZH2Input<-read.table("M2EZH2-input-genome1.depth.txt",sep='\t',header = FALSE)
M2EZH2_Depth<-M2EZH2Chip$V7-M2EZH2Input$V7
M2EZH2_Depth[M2EZH2_Depth<=0]=NA
M2WTChip<-read.table("2-WT-Chip-depth.txt",sep='\t',header = FALSE)
M2WTInput<-read.table("",sep='\t',header = FALSE)
M2WT_Depth<-M2WTChip$V7-M2WTInput$V7

M2WT_Depth[M2WT_Depth<=0]=NA

data<-data.frame(Depth=c(M2EZH2_Depth,M2WT_Depth),Group=rep(c("M2EZH2","M2WT"),c(length(M2EZH2_Depth),length(M2WT_Depth))))
data<-na.omit(data)

library(ggpubr)
my_comparisons<-c("M2EZH2","M2WT")
ggviolin(data,x="Group",y="Depth",fill="Group",palette = c("#00AFBB","#FC4E07"),add = "boxplot",
         add.params = list(fill="white"))+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  stat_compare_means(label.y = 1)

M2EZH2_Depth<-M2EZH2Chip$V7
M2WT_Depth<-M2WTChip$V7
data<-data.frame(Depth=c(M2EZH2_Depth,M2WT_Depth),Group=rep(c("M2EZH2","M2WT"),c(length(M2EZH2_Depth),length(M2WT_Depth))))
data<-na.omit(data)
my_comparisons<-c("M2EZH2","M2WT")
ggviolin(data,x="Group",y="Depth",fill="Group",palette = c("#00AFBB","#FC4E07"),add = "boxplot",
         add.params = list(fill="white"))+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  stat_compare_means(label.y = 1)
##############################################################################

Zygote_CChip_1<-read.table(filelist[5],sep='\t',header = FALSE)
Zygote_CInput_1<-read.table(filelist[7],sep='\t',header = FALSE)
Zygote_CChip_2<-read.table(filelist[6],sep='\t',header = FALSE)
Zygote_CInput_2<-read.table(filelist[8],sep='\t',header = FALSE)

Zygote_C_1_depth<-Zygote_CChip_1$V7-Zygote_CInput_1$V7
Zygote_C_1_depth[Zygote_C_1_depth<=0]=NA

Zygote_C_2_depth<-Zygote_CChip_2$V7-Zygote_CInput_2$V7
Zygote_C_2_depth[Zygote_C_2_depth<=0]=NA


Zygote_EChip_1<-read.table(filelist[9],sep='\t',header = FALSE)
Zygote_EInput_1<-read.table(filelist[11],sep='\t',header = FALSE)
Zygote_EChip_2<-read.table(filelist[10],sep='\t',header = FALSE)
Zygote_EInput_2<-read.table(filelist[12],sep='\t',header = FALSE)

Zygote_E_1_depth<-Zygote_EChip_1$V7-Zygote_EInput_1$V7
Zygote_E_1_depth[Zygote_E_1_depth<=0]=NA

Zygote_E_2_depth<-Zygote_EChip_2$V7-Zygote_EInput_2$V7
Zygote_E_2_depth[Zygote_E_2_depth<=0]=NA

data<-data.frame(Depth=c(Zygote_C_1_depth,Zygote_C_2_depth,Zygote_E_1_depth,Zygote_E_2_depth),
                 Zygote=rep(c("WT_M","WT_P","EZH2_M","EZH2_P"),c(length(Zygote_C_1_depth),length(Zygote_C_2_depth),length(Zygote_E_1_depth),length(Zygote_E_2_depth))))
data<-na.omit(data)
data$Zygote<-factor(data$Zygote,levels=c("WT_M","WT_P","EZH2_M","EZH2_P"))
my_comparisons<-list(c("WT_M","WT_P"),c("EZH2_M","EZH2_P"),c("WT_M","EZH2_M"),c("WT_P","EZH2_P"))
ggviolin(data,x="Zygote",y="Depth",fill="Zygote",palette = c("#00AFBB","#FC4E07","#E7B800","purple"),add = "boxplot",
         add.params = list(fill="white"))+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  theme(legend.position = "none")+
  scale_x_discrete(label=c("WT_M","WT_P","KO_M","KO_P"))
  #stat_compare_means(label.y = 1.5)

ggboxplot(data,x="Zygote",y="Depth",color="Zygote",palette = c("#00AFBB","#FC4E07","#E7B800","purple"),
          fill = "white",add = "dotplot", add.params = list(dotsize=0.1))+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+#stat_compare_means(label.y = 0.75)+
  theme(legend.position = "none",axis.title.x = element_blank())+ylab("Coverage")+
  theme(axis.text.x = element_text(size=10))

ggdotplot(data,x="Zygote",y="Depth",color="Zygote", 
          palette = c("#00AFBB","#FC4E07","#E7B800","purple"), binwidth = 0.001, ylab = "Zygote",
          add = "median_se", add.params = list(size=0.9,color="white"))

#################### coverage #############################
library(tidyverse)

fileslist=list.files()
index<-grep("depth",fileslist)
fileslist[index]

chr<-c(paste("chr",1:19,sep = ''),'chrX','chrY')

M2_EZH2<-read.table(fileslist[index][1],sep='\t',header = TRUE)
M2_EZH2<-M2_EZH2[M2_EZH2$Chr_name %in% chr,]
M2_EZH2$"Group"=rep("M2_EZH2",21)

M2_WT<-read.table(fileslist[index][2],sep='\t',header = TRUE)
M2_WT<-M2_WT[M2_WT$Chr_name %in% chr,]
M2_WT$"Group"=rep("M2_WT",21)

Zygote_WT_Maternal<-read.table(fileslist[index][3],sep='\t',header = TRUE)
Zygote_WT_Maternal<-Zygote_WT_Maternal[Zygote_WT_Maternal$Chr_name %in% chr,]
Zygote_WT_Maternal$"Group"=rep("WT_M",21)

Zygote_WT_Paternal<-read.table(fileslist[index][4],sep='\t',header = TRUE)
Zygote_WT_Paternal<-Zygote_WT_Paternal[Zygote_WT_Paternal$Chr_name %in% chr,]
Zygote_WT_Paternal$"Group"=rep("WT_P",21)

Zygote_EZH2_Maternal<-read.table(fileslist[index][5],sep='\t',header = TRUE)
Zygote_EZH2_Maternal<-Zygote_EZH2_Maternal[Zygote_EZH2_Maternal$Chr_name %in% chr,]
Zygote_EZH2_Maternal$"Group"=rep("KO_M",21)

Zygote_EZH2_Paternal<-read.table(fileslist[index][6],sep='\t',header = TRUE)
Zygote_EZH2_Paternal<-Zygote_EZH2_Paternal[Zygote_EZH2_Paternal$Chr_name %in% chr,]
Zygote_EZH2_Paternal$"Group"=rep("KO_P",21)

data<- rbind(Zygote_WT_Maternal,Zygote_WT_Paternal) %>% rbind(Zygote_EZH2_Maternal) %>% rbind(Zygote_EZH2_Paternal) %>%
  data.frame #%>% select(Chr_name,Coverage,Group)
data$Chr_name<-factor(data$Chr_name,levels = chr)
data$Coverage_new <- round(data$Coverage,4)
data$label<-paste(data$Coverage_new*100,'%',sep = '')

library(ggpubr)
ggplot(data,aes(x=Chr_name,y=Coverage_new,fill=Chr_name))+ylim(0,0.6)+
  geom_col(position = position_dodge2(preserve = "single"), width = 0.4)+
  facet_grid(Group~., scales = "free_x", space = "free_x")+theme_classic()+
  geom_text(aes(label=data$label),size=2.5,hjust=0.5, vjust=0.5,y=0.4)+
  theme(legend.position = 'none',axis.title.x = element_blank())+
  ylab("Coverage")


my_comparisons<-list(c("WT_M","WT_P"),
                     c("KO_M","KO_P"),
                     c("KO_M","WT_M"),
                     c("KO_P","WT_P"))
                     #c("M2_EZH2","M2_WT"),
                     #c("Zygote_WT_M","M2_WT"),
                     #c("Zygote_EZH2_M","M2_EZH2")


ggdotplot(data,x="Group",y="Coverage_new",color="Group",palette = c("blue","red","purple","green"),
          fill = "white", bins=30,add = "median_iqr", add.params = list(size=0.9),shape="Group")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+#stat_compare_means(label.y = 0.75)+
  theme(legend.position = "none",axis.title.x = element_blank())+ylab("Coverage")+
  theme(axis.text.x = element_text(size=10))


ggboxplot(data,x="Group",y="Coverage_new",color="Group",palette = c("#00AFBB","#FC4E07","#E7B800","purple"),
          fill = "white",add = "dotplot", add.params = list(dotsize=0.8))+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+#stat_compare_means(label.y = 0.75)+
  theme(legend.position = "none")+ylab("Coverage")+
  theme(axis.text.x = element_text(size=10))+xlab("Zygote")

write.csv(data,"../depth.csv")
#########################################

filelist<-list.files()

Zygote_CChip_1<-read.table(filelist[5],sep='\t',header = FALSE)
Zygote_CInput_1<-read.table(filelist[7],sep='\t',header = FALSE)
Zygote_CChip_2<-read.table(filelist[6],sep='\t',header = FALSE)
Zygote_CInput_2<-read.table(filelist[8],sep='\t',header = FALSE)

Zygote_C_1_depth<-Zygote_CChip_1$V7-Zygote_CInput_1$V7
Zygote_C_1_depth[Zygote_C_1_depth<=0]=NA

Zygote_C_2_depth<-Zygote_CChip_2$V7-Zygote_CInput_2$V7
Zygote_C_2_depth[Zygote_C_2_depth<=0]=NA


Zygote_EChip_1<-read.table(filelist[9],sep='\t',header = FALSE)
Zygote_EInput_1<-read.table(filelist[11],sep='\t',header = FALSE)
Zygote_EChip_2<-read.table(filelist[10],sep='\t',header = FALSE)
Zygote_EInput_2<-read.table(filelist[12],sep='\t',header = FALSE)

Zygote_E_1_depth<-Zygote_EChip_1$V7-Zygote_EInput_1$V7
Zygote_E_1_depth[Zygote_E_1_depth<=0]=NA

Zygote_E_2_depth<-Zygote_EChip_2$V7-Zygote_EInput_2$V7
Zygote_E_2_depth[Zygote_E_2_depth<=0]=NA

data<-data.frame(Depth=c(Zygote_C_1_depth,Zygote_C_2_depth,Zygote_E_1_depth,Zygote_E_2_depth),
                 Zygote=rep(c("WT_M","WT_P","EZH2_M","EZH2_P"),c(length(Zygote_C_1_depth),length(Zygote_C_2_depth),length(Zygote_E_1_depth),length(Zygote_E_2_depth))))
data<-na.omit(data)
data$Zygote<-factor(data$Zygote,levels=c("WT_M","WT_P","EZH2_M","EZH2_P"))
my_comparisons<-list(c("WT_M","WT_P"),c("EZH2_M","EZH2_P"),c("WT_M","EZH2_M"),c("WT_P","EZH2_P"))
ggviolin(data,x="Zygote",y="Depth",fill="Zygote",palette = c("#00AFBB","#FC4E07","#E7B800","purple"),add = "dotpoint",add.params = list(size=0.1))+
  stat_compare_means(comparisons = my_comparisons,label = "p.format")+
  theme(legend.position = "none")+
  scale_x_discrete(label=c("WT_M","WT_P","KO_M","KO_P"))


ggboxplot(data,x="Zygote",y="Depth",fill="Zygote",palette = c("#00AFBB","#FC4E07","#E7B800","purple"),add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,label = "p.format")+
  theme(legend.position = "none")+
  scale_x_discrete(label=c("WT_M","WT_P","KO_M","KO_P"))


library(PMCMR)
library(PMCMRplus)

library(pgirmess)
library(coin)
library(multcomp)

kruskal.test(Depth~Zygote,data=data)

kruskalmc(Depth~Zygote,data=data,probs=0.05)
mult <- oneway_test(Depth~Zygote,
                    data = data,
                    ytrafo = function(data) trafo(data, numeric_trafo = rank),
                    xtrafo = function(data) trafo(data, factor_trafo = function(x)
                      model.matrix(~x - 1) %*% t(contrMat(table(x), "Tukey"))),
                    teststat = "max",distribution = approximate(B = 90000))

pvalue(mult, method = "single-step")

kruskal_test(Depth~Zygote,data=data)
kruskalTest(Depth~Zygote,data=data)

###############################################






















