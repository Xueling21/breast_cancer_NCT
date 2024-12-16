#Due to file upload restrictions on the GitHub website, we have removed input and output files for code analysis processes that 
# are over 10MB in size. We have tried to save files for key steps as much as possible. If further analysis is needed,
# please download the data and analysis according to the comments. Please ensure that your computer server has 300GB 
# of running memory
####GSE163882_limma(Organize data）####
setwd("BRCA")
###
library(tidyverse)
library(GEOquery)
chooseBioCmirror()
gset = getGEO('GSE163882', destdir=".", AnnotGPL = F, getGPL = F)
class(gset)
gset[[1]]
pdata <- pData(gset[[1]])
exprSet <- read.table("GSE163882_all.data.tpms_222Samples.csv",
                      comment.char="!",
                      stringsAsFactors=F,
                      header=T,sep = ",")
class(exprSet)
exprSet <- exprSet[!duplicated(exprSet$annotation),]
rownames(exprSet) <- exprSet[,224]
exprSet <- exprSet[,-224]
exprSet <- exprSet[,-1]
exp <- exprSet

exp <- exp[apply(exp,1,mean)>0,]
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=2,varwidth=T)
# Select and filter out samples with large outliers in clinical information BA00154
del <- which(pdata$description%in%c("BA00154"))
pdata <- pdata[-del,]
which(colnames(exp)%in%c("BA00154"))
colnames(exp[,c(217)])
exp=as.data.frame(exp[,-c(217)])
library(limma)
exp = data.frame(exp)
exp=normalizeBetweenArrays(exp)
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=2)
range(exp)
#exp <- ceiling(2^(exp)-1)
exp <- log2(exp+1)
range(exp)
dev.off()

exp1 <- as.data.frame(exp)
x<-data.frame(GeneSymbol=rownames(exp1))
rownames(exp1)<-NULL
exp1 <- cbind(x,exp1)
write.table(exp1, file = "exp-GSE882.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(pdata, file = "pdata221.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

####exp221####
exp <- read.table("exp-GSE882.txt",sep = "\t")
exp <- t(exp)
colnames(exp) <- exp[1,]
exp <- exp[-1,]
exp <- as.data.frame(exp)
exp <- exp[order(exp$V1),]
rownames(exp) <- exp[,1]
exp <- exp[,-1]
exp <- t(exp)
exp <- as.data.frame(exp)

saveRDS(exp,"exp221,rds")
write.table(exp,"exp221.txt",sep = "\t",row.names=T,col.names=T)
exp <- read.table("exp221.txt",sep = "\t")

####meta_data221####
pdata <- read.table("pdata221.txt",sep = "\t")
colnames(pdata)
table(pdata$'response.to.nac.ch1')
meta_data2<- pdata$'response.to.nac.ch1'
meta_data1<- pdata$'description'
meta_data = cbind(meta_data1,meta_data2)
colnames(meta_data)=c("Sample","Group")
meta_data <- as.data.frame(meta_data)
meta_data <- meta_data[order(meta_data$Sample),]
rownames(meta_data) <- meta_data$Sample
saveRDS(meta_data,"meta_data221.rds")

####GSE163882_limma####
library(limma)
#1 no
#meta_data221
feature <- 'Group'
stage="RD"
test_idx <- which(meta_data[[feature]] == stage)
control_idx <- which(meta_data[[feature]] != stage)
test_num <- length(test_idx) 
control_num <- length(control_idx) 
class <- c(rep("2",test_num),rep("1",control_num))#Determine the comparison group, 2-1 (in alphabetical order)
class
design <- model.matrix(~class)#Determine the comparison group, where 2 equals 1
test_idx <- as.data.frame(test_idx)
control_idx <- as.data.frame(control_idx)
fit <- lmFit(rt[,c(test_idx,control_idx)],design)
fit <- eBayes(fit)
allDiff <- topTable(fit,adjust='BH',number=200000)
deg <- allDiff
#2
#meta_data221
group_list <- ifelse(str_detect(meta_data$Group, "RD"), "RD",
                     "pCR")
group_list = factor(group_list,
                    levels = c("pCR","RD"))
group_list
exp <- read.table("exp221.txt",sep = "\t")
exp <- log2(exp+1)
range(exp)
design=model.matrix(~group_list)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg2=topTable(fit,coef=2,number = Inf)
#end
write.table(deg2, file = "deg-GSE163882_all.txt",sep = "\t",row.names = T,col.names = T)



####heatmap####
#inpute deg
deg <- read.table("deg-GSE163882_all.txt",sep = "\t",header = T,row.name=1)
logFC=0.5
P.Value = 0.05
k1 = (deg$P.Value < P.Value)&(deg$logFC < -logFC)
k2 = (deg$P.Value < P.Value)&(deg$logFC > logFC)
deg$change = ifelse(k1,"down",ifelse(k2,"up","stable"))
table(deg$change)
#  down  stable   up 
#   178  46797   121 
rownames(deg) <- deg$gene_sample
cg = rownames(deg)[deg$change !="stable"]
diff=exp[cg,]
install.packages("pheatmap")
library(pheatmap)
annotation_col=data.frame(group=group_list)
rownames(annotation_col)=colnames(diff) 
pheatmap(diff,
         annotation_col=annotation_col,
         scale = "row",
         show_rownames = F,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         fontsize = 10,
         fontsize_row=3,
         fontsize_col=3)
dev.off()
####volcano plot####
install.packages('readxl')
BiocManager::install('ggrepel',force = TRUE)
library(readxl)
library(ggplot2)
library(ggrepel)


deg <- read.table("deg-GSE163882_all.txt",sep = "\t",header = T)
exprSet <- deg
#exprSet$'Gene Symbol'<- rownames(exprSet)
#colnames(exprSet) <- c("log2FC", "AveExpr", "t", "Pvalue", "adj.P.Val", "B", "Gene Symbol")
#saveRDS(exprSet,"huoshandeg.rds")
#exprSet <- huoshandeg
cut_off_Pvalue = 0.05
cut_off_logFC = 0.5

cut_off_Pvalue1 = 0.000001
cut_off_logFC1 = 0.5

exprSet$Sig <- ifelse(exprSet$P.Value < cut_off_Pvalue & 
                        abs(exprSet$logFC) >= cut_off_logFC, 
                      ifelse(exprSet$logFC > cut_off_logFC ,'Up','Down'),'no-DEGs')

exprSet <- data.frame(exprSet)
table(exprSet$Sig)

#1 gene
ggplot(exprSet, aes(x = logFC, y = -log10(P.Value), colour=Sig)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757")) + xlim(c(-2, 2)) + 
  geom_vline(xintercept=c(-cut_off_logFC,cut_off_logFC),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_Pvalue),
             lty=4,col="black",lwd=0.8) +
  labs(x="log(fold change)", y="-log10 (P value)") +
  theme_bw() +
  ggtitle("pCR vs RD") +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank() 
  ) +  
  geom_text_repel(
    data = subset(exprSet, exprSet$gene_sample==c("CCNA2","CCND1","KIAA1467","ESR1","CD52")),
    aes(label = gene_sample), size = 3,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
#2 line
ggplot(exprSet, aes(x = logFC, y = -log10(P.Value), colour=Sig)) +
  geom_point(alpha=0.4, size=2) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757")) + xlim(c(-2, 2)) + 
  geom_vline(xintercept=c(-cut_off_logFC,cut_off_logFC),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_Pvalue),
             lty=4,col="black",lwd=0.8) +
  labs(x="log(fold change)", y="-log10 (P value)") +
  theme_bw() +
  ggtitle("pCR vs RD") +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank() 
  ) +  
  geom_text_repel(
    data = subset(exprSet, exprSet$P.Value < cut_off_Pvalue1 & abs(exprSet$logFC) >= cut_off_logFC1),
    aes(label = gene_sample), size = 2.5,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE ,
    fontface = "bold" ,max.overlaps = 20)
#end
####GEO_enrichment analysis####
setwd("")
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
deg <- read.table("deg-GSE163882_all.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
logFC=0.3
P.Value = 0.05
k1 = (deg$P.Value < P.Value)&(deg$logFC < -logFC)
k2 = (deg$P.Value < P.Value)&(deg$logFC > logFC)
deg$change = ifelse(k1,"down",ifelse(k2,"up","stable"))
table(deg$change)
deg <- deg %>% filter(change!="stable")

DEG <- deg
#DEG <- DEG %>% rownames_to_column("Gene")

genelist <- bitr(DEG$gene_sample, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by=c("gene_sample"="SYMBOL"))


#GO分析
ego <- enrichGO(gene = DEG$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)

ego_res <- ego@result
write.table(ego_res,"brca1294-GO.txt",sep = "\t")
save(ego,ego_res,file = "GO_GSE163882_deg-0.3-0.05.Rdata")

##3.1 
barplot(ego, showCategory = 20,color = "pvalue")
##3.2 
dotplot(ego, showCategory = 20,label_format=100)
##3.3 
barplot(ego, drop = TRUE, showCategory =10,split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')
dotplot(ego,showCategory = 10,split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')

#KEGG
kk <- enrichKEGG(gene         = DEG$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05,
                 qvalueCutoff =0.05)
kk_res <- kk@result
kk_res1 <- kk_res[which(kk_res$p.adjust < 0.05),]
write.table(kk_res1,"brca1294-kegg.txt",sep = "\t")
save(kk,kk_res,file = "KEGG_GSE163882_DEG-0.3-0.05.Rdata")

#
barplot(kk, showCategory = 20,color = "pvalue")
#
dotplot(kk, showCategory = 20,label_format=100)

dev.off()

#1
test_name ="./result/GSE119409_300DEG_ttest.txt"
data_test=read.table(test_name,sep="\t",header=T)
#2
data <- 
  deg <- read.table("deg-GSE163882_all.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#end
data_test <- deg
#Select the top ranked genes by descending the absolute value of t.
data_test = arrange(data_test,-abs(data_test$p.value))
data_test = data_test[data_test$p.value<0.05,]

data=data[rownames(data) %in% deg$id,]
data=data[match(data_test$id,rownames(data)),]
data=data[1:30,]
data=t(data)

pdf(outpdf,height = 8,width = 15)
par(las=1,mar=c(8,5,3,3))                                                       
x=c(1:ncol(data))                                                               
y=c(1:ncol(data))
# 
plot(x,y,                                          
     xlim = c(0,3*ncol(data)),ylim = c(min(data)-0.02,max(data)+0.02),
     main="Top 30 DEGs of GSE119409",
     xlab = "",
     ylab = "Gene Expression",
     cex.main=1.8,
     cex.lab=1.8,
     cex.axis=1.8,
     pch=21,
     col="white",
     xaxt="n"
)
plot(legend("topright", fill = c("blue", "red"), legend = c("R", "NR"), title = "Sensitivity"))


library(vioplot)
#Draw a violin diagram for each immune cell cycle. Resistance is indicated in blue, sensitivity is indicated in red.

ncex=3
for(i in 1:ncol(data)){
  normalData=data[1:response,i]
  tumorData=data[(response+1):(response+non),i]
  vioplot(normalData,at=ncex*(i-1),lty=1,add=T,col="blue")
  vioplot(tumorData,at=ncex*(i-1)+1,lty=1,add=T,col="red")
  p=data_test$p.value[i]
  mx=max(c(normalData,tumorData))
  lines(c(x=ncex*(i-1)+0.1,x=ncex*(i-1)+0.8),c(mx,mx))                          
  text(x=ncex*(i-1)+0.5,y=mx+0.15,labels = convertP(p),cex=1.2)
  text(seq(1,ncex*ncol(data),ncex),1.5,xpd=NA,labels = colnames(data),cex = 1,srt=45,pos = 2) 
}

dev.off()

av <- available.packages(filters=list())
####Batch survival analysis(rnaseq)####
library(RTCGA)
infoTCGA <- infoTCGA() 
library(RTCGA.clinical)
library(dplyr)
clin <- survivalTCGA(BRCA.clinical) 
head(clin)
library(RTCGA.rnaseq)
dim(BRCA.rnaseq)

esprSet <- BRCA.rnaseq
table(substr(esprSet$bcr_patient_barcode,14,16))
esprSet <- esprSet[substr(esprSet$bcr_patient_barcode,14,16)%in% c("01A"),]
#1080

esprSet <- esprSet %>%
  # then make it a tibble (nice printing while debugging)
  as_tibble() %>%
  # then trim the barcode (see head(clin), and ?substr)
  mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) %>%
  # then join back to clinical data
  inner_join(clin, by="bcr_patient_barcode")
str(esprSet)
#colnames(esprSet[1:40])
#esprSet <- esprSet[,-(2:30)]
esprSet1 <- inner_join(clin,esprSet,by="bcr_patient_barcode")
#esprSet <- as.data.frame(esprSet)
#esprSet <- esprSet[apply(esprSet[,4:20534],1,mean)>0,]
esprSet1 <- esprSet[,c(1,20533,20534)]
esprSet <- esprSet[,-c(1,20533,20534)]
for (col in names(esprSet)) {
  if (!is.factor(esprSet[[col]])) {
    esprSet[[col]] <- as.numeric(esprSet[[col]])
  }
}
esprSet <- cbind(esprSet1,esprSet)
str(esprSet)
library(survival)
my.surv <- Surv(esprSet$times, esprSet$patient.vital_status)

#1 no
library(parallel)
cl.cores <- detectCores()
cl <- makeCluster(30)
clusterExport(cl,c("esprSet","my.surv"))
log_rank_p <- parApply(cl,esprSet[,25:length(names(esprSet))],2,function(values){
  
  group=ifelse(values>median(na.omit(values)),'high','low')
  
  if (TRUE & FALSE %in% group) {kmfit2 <- survival::survfit(my.surv~group,data=esprSet)
  
  #plot(kmfit2)
  
  data.survdiff=survival::survdiff(my.surv~group)
  
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  
}})
stopCluster(cl)
outTab <- log_rank_p

#2
outTab <- NULL
for(gene in colnames(esprSet[,4:ncol(esprSet)])){
  
  a=esprSet[,gene]<=median(esprSet[,gene])
  
  if (TRUE & FALSE %in% a){
    
    diff=survdiff(my.surv ~a,data = esprSet)
    
    pValue=1-pchisq(diff$chisq,length(diff$n) - 1)
    
    HR = (diff$obs[2]/diff$exp[2])/(diff$obs[1]/diff$exp[1])
    
    outTab=rbind(outTab,cbind(gene=gene,pvalue=pValue,HR=HR))
    
    fit <- survfit(my.surv ~ a, data = esprSet)
    
    summary(fit)
    
    pValue[is.na(pValue)] <- 1
    
    
}}
#end
outTab1 <- as.data.frame(outTab)
outTab1<-outTab1[order(as.numeric(as.vector(outTab1$pvalue))),]
saveRDS(outTab1,"gene_df_c.rds")

gene_diff <- outTab1$gene[outTab1$pvalue<0.05]
gene_diff_po <- outTab1$gene[(outTab1$pvalue<0.05)&(outTab1$HR > 1)]
gene_diff_na <- outTab1$gene[(outTab1$pvalue<0.05)&(outTab1$HR < 1)]

library(survminer)
library(survival)
library(tidyverse)
gene_df <- readRDS("gene_df_c.rds")
gene_df1<-as.data.frame(str_split(gene_df$gene,'[|]',simplify = T)[,1])
gene_df1 <- cbind(gene_df,gene_df1)
gene_df1[,1] <- gene_df1[,4]
gene_df1 <- gene_df1[,-4]
colnames(gene_df1) <- c("gene_sample","P.Value","HR")
which(gene_df1$gene_sample%in%c("KIAA1467","TYROBP","CD5","CD2","RAC2","KRT17","IGJ","CHI3L1","MMP7","ZIC1",
                                "TRAT1","CXCL1","CD52","KRT5","KRT16","EVL","SERPINB5","APOBEC3G","HLA-DOB","MIA"))
gene_df2=as.data.frame(gene_df1$gene_sample[c(141,505,602,652,719,827,1030,1186,1623,1697,1987,1999,2249,2313,2512,2862,2969,3210,3255,
                                              3474)])
colnames(gene_df2) <- "gene_sample"
gene_df1 <- merge(gene_df1,gene_df2,by='gene_sample')
write.table(gene_df1,"gene_df20",sep = "\t")


custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5)
    )
}
esprSet1 <- esprSet[,-(4:32)]
esprSet2 <- esprSet1[,c(1:3)]
esprSet1 <- esprSet1[,-c(1:3)]
esprSet3<-as.character(str_split(colnames(esprSet1),'[|]',simplify = T)[,1])
colnames(esprSet1) <- esprSet3
esprSet1 <- cbind(esprSet2,esprSet1)
fig <- function(a,b,c,d){
  group <- ifelse(a>median(a),'high','low')
  sfit <- survfit(Surv(times, patient.vital_status)~group, data=esprSet1)
  p <- ggsurvplot(sfit,palette=c("#E7B800","#2E9FDF"),
                   surv.median.line = "hv",
                   pval =F,
                   conf.int =TRUE,
                   xlab ="Times(days)",
                   title=b,
                   ggtheme=custom_theme())
  
  p$plot = p$plot + ggplot2::annotate("text",x = 1300, y = 0.2,
                                        label = paste(c))+
    ggplot2::annotate("text",x = 1300, y = 0.1,
                      label = paste(d))
  
}

gene_df20 <- read.table("gene_df20",sep = "\t")
round(gene_df20$P.Value, digits=3)
round(gene_df20$HR, digits=2)
p1 <- fig(esprSet1$KIAA1467,"KIAA1467","p= 0.008","Hazard Ratio= 0.64")#5x5
p1
p2 <- fig(esprSet1$APOBEC3G,"APOBEC3G","p= 0.027","Hazard Ratio= 1.69")
p2
p3 <- fig(esprSet1$CD2,"CD2","p= 0.014","Hazard Ratio= 1.63")
p4 <- fig(esprSet1$CD5,"CD5","p= 0.041","Hazard Ratio= 1.50")
p5 <- fig(esprSet1$CD52,"CD52","p= 0.003","Hazard Ratio= 1.83")
p6 <- fig(esprSet1$CHI3L1,"CHI3L1","p= 0.033","Hazard Ratio= 1.52")
p7 <- fig(esprSet1$CXCL1,"CXCL1","p= 0.023","Hazard Ratio= 1.57")
p8 <- fig(esprSet1$EVL,"EVL","p= 0.004","Hazard Ratio= 1.77")
p9 <- fig(esprSet1$`HLA-DOB`,"HLA-DOB","p= 0.024","Hazard Ratio= 1.56")
p9
p10 <- fig(esprSet1$IGJ,"IGJ","p= 0.018","Hazard Ratio= 1.59")
p10
p11 <- fig(esprSet1$KRT16,"KRT16","p= 0.035","Hazard Ratio= 1.51")
p12 <- fig(esprSet1$KRT17,"KRT17","p= 0.003","Hazard Ratio= 1.80")
p13 <- fig(esprSet1$KRT5,"KRT5","p= 0.006","Hazard Ratio= 1.71")
p14 <- fig(esprSet1$MIA,"MIA","p= 0.004","Hazard Ratio= 1.78")
p15 <- fig(esprSet1$MMP7,"MMP7","p= 0.013","Hazard Ratio= 1.63")
p16 <- fig(esprSet1$RAC2,"RAC2","p= 0.0003","Hazard Ratio= 2.08")
p17 <- fig(esprSet1$SERPINB5,"SERPINB5","p= 0.002","Hazard Ratio= 1.83")
p18 <- fig(esprSet1$TRAT1,"TRAT1","p= 0.018","Hazard Ratio= 1.61")
p19 <- fig(esprSet1$TYROBP,"TYROBP","p= 0.041","Hazard Ratio= 1.51")
p20 <- fig(esprSet1$ZIC1,"ZIC1","p= 0.046","Hazard Ratio= 1.48")


pz<- ggarrange(p1, p2, p3, p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15, ncol = 3, nrow = 5,
          labels = c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O"), 
          label.x = 0, label.y = 1,
          font.label = list(size = 20, face = "bold")) 
pz
pz2 <- ggarrange(p16,p17,p18,p19,p20, ncol = 3, nrow = 5,
          labels = c("P","Q","R","S","T"), 
          label.x = 0, label.y = 1, 
          font.label = list(size = 20, face = "bold")) 
pz2
#A4 paper size, 21x29.7, with 2 blank spaces on top, bottom, left, and right, 17x25
# end

####vnn####
library(tidyverse)
exprSet <- read.table("GSE163882_all.data.tpms_222Samples.csv",
                      comment.char="!",      
                      stringsAsFactors=F,
                      header=T,sep = ",")
class(exprSet)
exprSet <- exprSet[!duplicated(exprSet$annotation),]

exprSet[,2] <- exprSet1[,224]
exprSet$'Gene ID' <- exprSet[,1]
exprSet$gene_sample <- exprSet[,2]
exprSet <- exprSet[,-1]
deg300 <- read.table("deg300_GSE163882.txt",sep = "\t")
deg300$gene_sample <- rownames(deg300)

deg300 <- inner_join(exprSet,deg300,by="gene_sample")
#deg300 <- merge(exprSet,deg300,by='gene_sample',
#                          all = T, sort = T)
write.table(deg300,"deg300-GSE163882.txt",sep = "\t")
deg300 <- read.table("deg300-GSE163882.txt",header=T)
#sur <- read.table("brca_survival_all.txt",sep = "\t")
#rownames(sur) <-sur[,1] 
#colnames(sur) <-sur[1,] 
#sur <- sur[-1,-1]
#sur$`Gene ID`<- substr(sur$`Gene ID`,1,15)

venn_dat2 <- as.data.frame(gene_diff)
venn_dat1 <- as.data.frame(deg300$`gene_sample`)
venn_dat2<-as.data.frame(str_split(venn_dat2$gene_diff,'[|]',simplify = T)[,1])
colnames(venn_dat2) <- 'gene_sample'



library(stats)
phyper(77-1,3569,47096-3569,299,lower.tail=F)
#7.016496e-22
#Among 47096 sequenced genes, first search for TCGA genes, prognostic related gene sets, then search for GEO differential gene sets, and finally take the intersection; 77 chip genes, prognosis related, GEO differential genes; 3569 TCGA genes are prognostic related, with a total of 47096 sequenced genes,
#GEO differentially expressed genes among 299 sequenced genes；

####VNN（select up genes）####
library(tidyverse)
write.table(deg300,"deg300-GSE163882.txt",sep = "\t")
deg300 <- read.table("deg300-GSE163882.txt",header=T)
deg_up <- deg300 %>% filter(change=="up")
#sur <- read.table("brca_survival_all.txt",sep = "\t")
#rownames(sur) <-sur[,1] 
#colnames(sur) <-sur[1,] 
#sur <- sur[-1,-1]
#sur$`Gene ID`<- substr(sur$`Gene ID`,1,15)
deg300 <- deg_up
gene_diff <- gene_diff_na
venn_dat2 <- as.data.frame(gene_diff)
venn_dat1 <- as.data.frame(deg300$`gene_sample`)
venn_dat2<-as.data.frame(str_split(venn_dat2$gene_diff,'[|]',simplify = T)[,1])
colnames(venn_dat2) <- 'gene_sample'


library(VennDiagram)
venn_list <- list(deg_gene = venn_dat1$'deg300$gene_sample', sur_gene = venn_dat2$'gene_sample')
venn.diagram(venn_list, 
             
             filename = 'deg_sur_neg.png', 
             
             imagetype = 'png', 
             
             euler.d = TRUE,
             
             scaled = F, 
             
             fill=c("red","blue"), 
             
             alpha=c(0.65,0.3), 
             
             cex=2,
             
             cat.cex = 1.5, 
             category.names = c("sur_neg","deg_up"),
             cat.pos = c(-5,5), 
             
             margin = 0.1); #Number giving the amount of whitespace around the diagram in grid units

inter <- get.venn.partitions(venn_list) 
interset_mRNA<-as.data.frame(inter$..values..[1]) 
interset_mRNA$'gene_sample' <- interset_mRNA[,1]
deg77vnn <- inner_join(interset_mRNA,deg300,by=c("gene_sample"))
write.table(deg77vnn, 'deg_sur_neg.txt', col.names = T,row.names = T, sep = '\t', quote = FALSE)

####VNN（select down genes）####
library(tidyverse)
deg300 <- read.table("deg300-GSE163882.txt",header=T)
deg_down <- deg300 %>% filter(change=="down")
#sur <- read.table("brca_survival_all.txt",sep = "\t")
#rownames(sur) <-sur[,1] 
#colnames(sur) <-sur[1,] 
#sur <- sur[-1,-1]
#sur$`Gene ID`<- substr(sur$`Gene ID`,1,15)
deg300 <- deg_down
gene_diff <- gene_diff_po

venn_dat2 <- as.data.frame(gene_diff)
venn_dat1 <- as.data.frame(deg300$`gene_sample`)
venn_dat2<-as.data.frame(str_split(venn_dat2$gene_diff,'[|]',simplify = T)[,1])
colnames(venn_dat2) <- 'gene_sample'


library(VennDiagram) 

venn_list <- list(deg_gene = venn_dat1$'deg300$gene_sample', sur_gene = venn_dat2$'gene_sample')

venn.diagram(venn_list, 
             
             filename = 'deg_sur_pos.png', 
             
             imagetype = 'png', 
             
             euler.d = TRUE,
             
             scaled = F, 
             
             fill=c("red","blue"), 
             
             alpha=c(0.65,0.3), 
             
             cex=2, 
             
             cat.cex = 1.5, 
             category.names = c("sur_pos","deg_down"),
             cat.pos = c(-5,5), 
             
             margin = 0.1); #Number giving the amount of whitespace around the diagram in grid units
inter <- get.venn.partitions(venn_list) ##get.venn.partitions
interset_mRNA<-as.data.frame(inter$..values..[1]) 
interset_mRNA$'gene_sample' <- interset_mRNA[,1]
deg77vnn <- inner_join(interset_mRNA,deg300,by=c("gene_sample"))
write.table(deg77vnn, 'deg_sur_pos.txt', col.names = T,row.names = T, sep = '\t', quote = FALSE)

interset_mRNA$gene_sample

####VNN（Integrate）####
outTab1 <- readRDS("gene_df_c.rds")
outTab1 <- outTab1[-which(str_split(outTab1$gene,'[|]',simplify = T)[,1]=="?"),]

gene_diff_po <- outTab1$gene[(outTab1$pvalue<0.05)&(outTab1$HR > 1)]
gene_diff_na <- outTab1$gene[(outTab1$pvalue<0.05)&(outTab1$HR < 1)]
gene_diff_po <-as.data.frame(str_split(gene_diff_po,'[|]',simplify = T)[,1])#HR<1
gene_diff_na <- as.data.frame(str_split(gene_diff_na,'[|]',simplify = T)[,1])#HR>1

library(tidyverse)
deg300 <- read.table("deg300-GSE163882.txt",header=T)
deg_down <- deg300 %>% filter(change=="down")
deg_up <- deg300 %>% filter(change=="up")


venn_dat1 <- as.data.frame(deg_up$`gene_sample`)
venn_dat2 <- as.data.frame(deg_down$`gene_sample`)
venn_dat3 <- as.data.frame(gene_diff_na)
venn_dat4 <- as.data.frame(gene_diff_po)
head(venn_dat4)
#venn_dat1-4
colnames(venn_dat1) <- 'gene_sample'
colnames(venn_dat2) <- 'gene_sample'
colnames(venn_dat3) <- 'gene_sample'
colnames(venn_dat4) <- 'gene_sample'

library(VennDiagram) 
venn_list <- list(venn_dat1$gene_sample,venn_dat2$gene_sample,venn_dat3$gene_sample,venn_dat4$gene_sample)
venn.diagram(venn_list, 
             
             filename = 'deg_sur_2.png', 
             
             imagetype = 'png', 
             
             euler.d = TRUE,
             
             scaled = F, 
             
             fill=c("#D693BE","#8EC8ED","#F5B3A5","#AED594"), 
             
             alpha=c(0.6,0.5,0.4,0.3), 
             
             cex=1.5, 
             
             cat.cex = 1, 
             category.names = c("Up-regulated \n DEGs","Down-regulated \n DEGs","HR>1","HR<1"),
             cat.pos = c(0,0,0,0), 
             cat.dist = c(0.25,0.25,0.1,0.1),
             margin = 0.1)
inter <- get.venn.partitions(venn_list) 
interset_mRNA<-as.data.frame(inter$..values..[1]) 
interset_mRNA$'gene_sample' <- interset_mRNA[,1]
deg77vnn <- inner_join(interset_mRNA,deg300,by=c("gene_sample"))
write.table(deg77vnn, 'deg_sur_pos.txt', col.names = T,row.names = T, sep = '\t', quote = FALSE)

interset_mRNA$gene_sample


####Machine Learning####
####RF####
#exp221
exp<- readRDS("exp221.rds")
which(rownames(exp)%in%c("KIAA1467","TYROBP","CD5","CD2","RAC2","KRT17","IGJ","CHI3L1","MMP7","ZIC1","TRAT1",
                          "CXCL1","CD52","KRT5","KRT16","EVL","SERPINB5","APOBEC3G","HLA-DOB","MIA",
                         "CCND1","ESR1","CCL2","PTPRC","MMP9","STAT1","CXCL8","CD274","GATA3","CCNA2"))
rownames(exp[c(307,1531,1627,1943,2416,3416,3569,3742,3777,4418,4624,5015,6008,6014,6579,6688,7492,8679,9587,
               11017,11101,12575,12580,16101,16281,17028,19145,30477,30955,39858),])
exp30=as.data.frame(exp[c(307,1531,1627,1943,2416,3416,3569,3742,3777,4418,4624,5015,6008,6014,6579,6688,7492,8679,9587,
                          11017,11101,12575,12580,16101,16281,17028,19145,30477,30955,39858),])
data <-  as.data.frame(exp30)
str(data)
for (col in names(data)) {
  if (!is.factor(data[[col]])) {
    data[[col]] <- as.numeric(data[[col]])
  }
}
str(data)
data <- log2(data+1)
range(data)
#meta_data221
meta_data<- readRDS("meta_data221.rds")
feature <- 'Group'
stage="RD"
test_idx <- which(meta_data[[feature]] == stage)
control_idx <- which(meta_data[[feature]] != stage)
test_num <- length(test_idx) 
control_num <- length(control_idx) 
class <- c(rep("1",control_num),rep("2",test_num))


#install.packages('randomForest')
library('randomForest')
#randomForest

control_idx <- as.data.frame(control_idx)
test_idx <- as.data.frame(test_idx)
datal = data[,control_idx[,1]]
data2 = data[,test_idx[,1]]
data = cbind(datal,data2)
#RF
x=as.matrix(t(data))
col_names <- gsub("-","", colnames(x))
colnames(x) <- col_names
y=class
set.seed(78650)
rf=randomForest(as.factor(y)~.,data=x,ntree=500)
plot(rf, main="Random forest", lwd=2)
dev.off()
optionTrees=which.min(rf$err.rate[,1])
optionTrees
#293
rf2=randomForest(as.factor(y)~., data=x,ntree=optionTrees)
importance=importance(x=rf2)
rfGenes=importance[order(importance[,"MeanDecreaseGini"],decreasing = TRUE),]
rfGenes=names(rfGenes[rfGenes>3])
randomForest::varImpPlot(rf2,main="")

rfGenes=names(rfGenes[1:14])
#"CCNA2"    "ESR1"     "GATA3"    "CXCL1"    "CD2"      "CHI3L1"   "CD52"     "CCL2"     "APOBEC3G"
#"MMP9"     "TRAT1"    "KIAA1467" "CCND1"    "CD274"  
saveRDS(rfGenes,"rfGenes1.rds")
####SVM####
data <-  as.data.frame(exp30)
for (col in names(data)) {
if (!is.factor(data[[col]])) {
data[[col]] <- as.numeric(data[[col]])
}
}
data <- log2(data+1)
range(data)
#meta_data221
data <- t(data)
data <- as.data.frame(data)
data$Sample<- rownames(data)
data <- merge(meta_data,data,by="Sample")
rownames(data) <- data[,1]
data <- data[,-1]
str(data)
data$Group <- as.factor(data$Group)

set.seed(346351)
library(tidyverse)
library(e1071)
library(caret)
source("msvmRFE.R")
input <- data

#1 no
nfold = 10 
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, input, k=10, halve.above=100)
top.features = WriteFeatures(results, input, save=F)
featsweep = lapply(1:30, FeatSweep.wrap, results, input)
no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))
#pdf("BRCA19_svm.pdf", height = 6, width = 8)
PlotErrors(errors, no.info=no.info)
top.features22 <- top.features[1:22,]
saveRDS(top.features8,"svmGenes.rds")

#2
library(parallel)
nfold = 10
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
#make a cluster
#cl <- makeMPIcluster(mpi.universe.size())
cl <- makeCluster(50)
clusterExport(cl, list("input","svmRFE","getWeights","svm"))
results <-parLapply(cl,folds, svmRFE.wrap, input, k=10, halve.above=100)
top.features = WriteFeatures(results, input, save=F)
clusterExport(cl, list("top.features","results", "tune","tune.control"))
featsweep = parLapply(cl,1:30, FeatSweep.wrap, results, input)
stopCluster(cl)
no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))
PlotErrors(errors, no.info=no.info)
plot(top.features)
top.features16 <- top.features[1:16,]
saveRDS(top.features16$FeatureName,"svmGenes.rds")
#end

####lasso####
library(glmnet)
set.seed(734568)
#60 65 68
train <- data#SVMzheng li data
x <- as.matrix(train[,-1])
y <- ifelse(train$Group == "pCR", 0,1) 
fit = glmnet(x, y, family = "binomial", alpha = 1, lambda = NULL)
plot(fit, xvar = "dev", label = TRUE)
cvfit = cv.glmnet(x, y,
                  nfold=10,
                  family = "binomial", type.measure = "class")
plot(cvfit)
cvfit$lambda.min
myCoefs <- coef(cvfit, s="lambda.min")
lasso_fea <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
(lasso_fea <- lasso_fea[-1])
saveRDS(lasso_fea,"lassoGenes.rds")
####VNN####
library(VennDiagram)
rfGenes <- readRDS("rfGenes1.rds")
svmGenes <- readRDS("svmGenes.rds")
lassoGenes <- readRDS("lassoGenes.rds")
venn_list <- list(   SVM_RFE=svmGenes,     RandomForest=rfGenes,   LASSO=lassoGenes)
venn.plot<-venn.diagram(venn_list,
                        filename = 'machinelearn12.png', 
                        imagetype = 'png', 
                        euler.d = TRUE,
                        scaled = F, 
                        fill=c("#7EA6D9","#D4E6BC","#F8F2A4"), 
                        alpha=c(0.6,0.45,0.3), 
                        cex=2, 
                        cat.cex = 1.3, 
                        cat.dist = c(0.1,0.15,0.2),
                        margin = 0.2, #Number giving the amount of whitespace around the diagram in grid units
                        category.names = c("SVM_RFE","LASSO","RandomForest"),
                        cat.pos = c(-30,-22,60))

dev.new()
inter <- get.venn.partitions(venn_list) 
interset_mRNA<-as.data.frame(inter$..values..[1]) 
interset_mRNA$'gene_sample' <- interset_mRNA[,1]
saveRDS(interset_mRNA,"vnn5.rds")
interset_mRNA$gene_sample

library(VennDiagram)
rfGenes <- readRDS("rfGenes.rds")
svmGenes <- readRDS("svmGenes.rds")
lassoGenes <- readRDS("lassoGenes.rds")
venn_list <- list(   SVM_RFE=svmGenes,     LASSO=lassoGenes,RandomForest=rfGenes)
venn.plot<-venn.diagram(venn_list,
                        filename = 'machinelearn12.png', 
                        imagetype = 'png', 
                        euler.d = TRUE,
                        scaled = F, 
                        fill=c("red","blue","yellow"),
                        alpha=c(0.65,0.3,0.475), 
                        cex=2, 
                        cat.cex = 1.5, 
                        cat.dist = c(0.08,0.08,0.05),
                        margin = 0.1, 
                        category.names = c("SVM_RFE","LASSO","RandomForest"),
                        force.unique=T)
dev.new()
inter <- get.venn.partitions(venn_list) 
interset_mRNA<-as.data.frame(inter$..values..[1]) 
interset_mRNA$'gene_sample' <- interset_mRNA[,1]
saveRDS(interset_mRNA,"vnn5.rds")
interset_mRNA$gene_sample

####nomogram####
vnn5 <- readRDS("vnn5.rds")
exp <- readRDS("exp178.rds")
exp <- readRDS("exp221.rds")
which(rownames(exp)%in%c("CCNA2","CCND1","KIAA1467","ESR1","CD52"))
rownames(exp[c(1627,1943,3742,8679,12580),])
exp5=as.data.frame(exp[c(1627,1943,3742,8679,12580),])
data <-  as.data.frame(exp5)
for (col in names(data)) {
  if (!is.factor(data[[col]])) {
    data[[col]] <- as.numeric(data[[col]])
  }
}
data <- log2(data+1)
range(data)
#meta_data221
data <- t(data)
data <- as.data.frame(data)
data$Sample<- rownames(data)
data <- merge(meta_data,data,by="Sample")
rownames(data) <- data[,1]
data <- data[,-1]
str(data)
#data$Group <- as.factor(data$Group)
Mydata <- data
library(Hmisc) 
library(rms)
dd<-datadist(Mydata)  
options(datadist="dd") 
f_lrm <-lrm(Group~KIAA1467+CCND1+CCNA2+ESR1+CD52
            , data=Mydata)     
summary(f_lrm)    
par(mgp=c(1.6,0.6,0),mar=c(5,5,3,1))     
nomogram <- nomogram(f_lrm,fun=function(x)1/(1+exp(-x)),  
                     fun.at = c(0.01,0.05,0.2,0.5,0.9,0.99),
                     funlabel = "Prob of RD",    
                     conf.int = F,  
                     abbrev = F)  
plot(nomogram)    
saveRDS(nomogram,"BRCA30_nomogram.RDS")
####DCA####
#par(mgp=c(2.4,0.3,0),mar=c(6,5,3,3))     
fit2<-lrm(Group~KIAA1467+CCND1+CCNA2+ESR1+CD52
          , data=Mydata,x=TRUE,y=TRUE) 
cal1<-calibrate(fit2,method = "boot",B=1000)
P1 <- predict(fit2, type = "fitted")
Group_numeric <- as.factor(Mydata$Group)
val.prob(P1, Group_numeric)

par(mar=c(5,5,2,1))
plot(cal1,
     xlim = c(0,1.0),
     ylim = c(0,1.0),
     xlab = "Nomogram Predicted RD",
     ylab="Actual RD",
     xaxs = "i",
     yaxs = "i",
     legend =FALSE,
     subtitles = FALSE,
     cex=1.2,
     cex.axis=1.2,
     cex.lab=1.2)
abline(0,1,col="#4370B4",lty=2,lwd=2)
lines(cal1[,c("predy","calibrated.orig")],type="l",lwd=2,col="#C30078")
lines(cal1[,c("predy","calibrated.corrected")],type="l",lwd=2,col="#549F9A")
legend(x=0.58,y=0.35,
legend=c("Ideal","Apparent","Bias-corrected"),
lty = c(2,1,1),
lwd = c(2,1,1),
col = c("#4370B4","#C30078","#549F9A"),
bty="n",
cex=1.2)

####DCA####
library(rmda) 
dim(Mydata) 
head(Mydata) 
Mydata$Group <- dplyr::recode(Mydata$Group, "pCR" = "0")
Mydata$Group <- dplyr::recode(Mydata$Group, "RD" = "1")
str(Mydata) 
Mydata$Group <- as.numeric(Mydata$Group)
fit3<- decision_curve(Group~KIAA1467+CCND1+CCNA2+ESR1+CD52
                      , data=Mydata,study.design = "cohort",
                      bootstraps =500) 
plot_decision_curve(fit3, curve.names = "DCA",
                    cost.benefit.axis = F, 
                    confidence.intervals = "none") 
####pROC####
nom <- nomogram
library(nomogramFormula)
results<-formula_rd(nomogram=nom)
Mydata$points<-points_cal(formula = results$formula,rd=Mydata)
pre<-Mydata$points
library(pROC)
plot.roc(Mydata$Group, pre,
         main="ROC Curve", percent=TRUE,
         print.auc=TRUE,
         ci=TRUE,)
rocplot1 <- roc(Mydata$Group,pre)
ci.auc(rocplot1)
rocobj1 <- plot.roc(Mydata$Group, Mydata$KIAA1467,percent=TRUE, col="#95A8D2")
rocobj2 <- lines.roc(Mydata$Group, Mydata$CCND1, percent=TRUE, col="#C3EAB5")
rocobj3 <- lines.roc(Mydata$Group, Mydata$CCNA2, percent=TRUE, col="#C5767B")
rocobj4 <- lines.roc(Mydata$Group, Mydata$ESR1, percent=TRUE, col="#870A4C")
rocobj5 <- lines.roc(Mydata$Group, Mydata$CD52, percent=TRUE, col="#000084")
legend("bottomright", legend=c("KIAA1467", "CCND1","CCNA2","ESR1","CD52"), col=c("#95A8D2", "#C3EAB5","#C5767B","#870A4C","#000084"),
lwd=1.5)
rocobj5 <- lines.roc(Mydata$Group, pre, percent=TRUE, col="black")
legend("bottomright", legend=c("KIAA1467","CCND1","CCNA2","ESR1","CD52","5 Genes"), col=c("#95A8D2", "#C3EAB5","#C5767B","#870A4C","#000084","black"),
lwd=1.5)
####PRROC####
library(PRROC)
library(tidyverse)
pre <- as.numeric(pre)
weight <- ifelse(Mydata$Group == "pCR", 0,1)
# 
pr <- pr.curve(scores.class0 = pre, weights.class0 = weight,
               curve=TRUE,rand.compute = TRUE)
plot(pr, rand.plot = TRUE, auc.main=T,legend=F,lwd = 2,cex.main=1.2,col="black")


####GSE20271_limma（Organize data）####
library(tidyverse)
library(GEOquery)
chooseBioCmirror()
options('download.file.method.GEOquery'='auto')
options('GEOquery.inmemory.gpl'=FALSE)
gset = getGEO('GSE20271', destdir=".", AnnotGPL = T)
gset[[1]]
class(gset)
pdata <- pData(gset[[1]])
exprSet = exprs(gset[[1]])
class(exprSet)
exprSet <- as.data.frame(exprSet)
exp <- exprSet
exp <- exp[apply(exp,1,mean)>0,]
colnames(pdata)
group_list<-str_split(pdata$characteristics_ch1,' ',simplify = T)[,4]
group_list = factor(group_list,levels = c("pCR","RD"))  
table(group_list)
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=2,varwidth=T)
library(limma)
exp=normalizeBetweenArrays(exp)
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=2)
range(exp)
index = gset[[1]]@annotation

#1
gse_gp<-getGEO('GPL96',destdir =".")
annotation <- gse_gp@dataTable@table
platform_file_set=annotation[,c(1,11)]
exp$ID <- rownames(exp)
exp <- merge(x = exp, y = platform_file_set, by.x = "ID")
#2 no
BiocManager::install("hgu133a.db")
library(hgu133a.db)
ls("package:hgu133a.db")
ids <- toTable(hgu133aSYMBOL)
head(ids)
which(ids$symbol%in%c("KIAA1467"))
#length(unique(ids$symbol))
#table(sort(table(ids$symbol)))
#id转换
library(tidyverse)
exp <- as.data.frame(exp)
exp <- exp %>% mutate(probe_id=rownames(exp))
exp <- exp %>% inner_join(ids,by="probe_id") 
#end

exp <- exp[!duplicated(exp$`Gene Symbol`),]
rownames(exp) <- exp$`Gene Symbol`
exp <- exp[,-c(1,180)]
write.table(pdata, file = "pdata178.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
exp1 <- exp

####exp178####
saveRDS(exp,"exp178.rds")

####meta_data178####
pdata <- read.table("pdata178.txt",sep = "\t",header = T,row.names = 1)
colnames(pdata)
group_list<-str_split(pdata$characteristics_ch1,' ',simplify = T)[,4]
table(group_list)
meta_data2<- group_list
meta_data1<- pdata$geo_accession
meta_data = as.data.frame(cbind(meta_data1,meta_data2))
colnames(meta_data)=c("Sample","Group")
rownames(meta_data) <- meta_data$Sample
saveRDS(meta_data,"meta_data178.rds")


####pROC####
exp1 <- readRDS("exp178.rds")
which(rownames(exp1)%in%c("KIAA1467","CCND1","CCNA2","ESR1","CD52"))
rownames(exp1[c(2070,3041,3501,6403,8372),])
exp5=as.data.frame(exp1[c(2070,3041,3501,6403,8372),])
data1 <-  as.data.frame(exp5)
data1 <- log2(data1+1)
range(data1)
#meta_data221
meta_data1 <- readRDS("meta_data178.rds")
data1 <- t(data1)
data1 <- as.data.frame(data1)
data1$Sample<- rownames(data1)
data1 <- merge(meta_data1,data1,by="Sample")
rownames(data1) <- data1[,1]
data1 <- data1[,-1]
str(data1)
#data$Group <- as.factor(data$Group)
Mydata1 <- data1
nomogram <- readRDS("BRCA30_nomogram.RDS")
nom <- nomogram
library(nomogramFormula)
results<-formula_rd(nomogram=nom)
Mydata1$points<-points_cal(formula = results$formula,rd=Mydata1)
pre<-Mydata1$points
library(pROC)
plot.roc(Mydata1$Group, pre,
         main="ROC Curve", percent=TRUE,
         print.auc=TRUE,
         ci=TRUE,)
library(PRROC)
library(tidyverse)
pre <- as.numeric(pre)
weight <- ifelse(Mydata1$Group == "pCR", 0,1)
# 
pr <- pr.curve(scores.class0 = pre, weights.class0 = weight,
               curve=TRUE,rand.compute = TRUE)
plot(pr, rand.plot = TRUE, auc.main=T,legend=F,lwd = 2,cex.main=1.2,col="black")


####oncopredict####
options(stringsAsFactors = F)
BiocManager::install("oncoPredict")
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggplot2)
library(ggpubr)
th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
GDSC2_Expr = readRDS(file='GDSC2_Expr (RMA Normalized and Log Transformed).rds')
GDSC2_Res = readRDS("GDSC2_Res.rds")

GDSC2_Expr <- log2(GDSC2_Expr+1)
range(GDSC2_Expr)
testExpr <- readRDS("exp221.rds")
testExpr <- t(testExpr)
testExpr <-  as.data.frame(testExpr)
for (col in names(testExpr)) {
  if (!is.factor(testExpr[[col]])) {
    testExpr[[col]] <- as.numeric(testExpr[[col]])
  }
}
testExpr <- log2(testExpr+1)
range(testExpr)
nomogram <- readRDS("BRCA30_nomogram.RDS")
nom <- nomogram
library(nomogramFormula)
results<-formula_rd(nomogram=nom)
testExpr$points<-points_cal(formula = results$formula,rd=testExpr)
pre<-testExpr$points
###meta_data221_pre###
meta_data<- readRDS("meta_data221.rds")
table(meta_data$Group)
#pCR  RD
#80 141
round(80/(80+141),2)
#0.36
testExpr$Group<- ifelse(testExpr$points>quantile(testExpr$points,.36),"High","Low")
table(testExpr$Group)
#High  Low 
#141   80 
meta_data$Grouppr <- testExpr$Group
saveRDS(meta_data,"meta_data221_predict.rds")
write.table(meta_data,file = "meta_data221_predict.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
testExpr1 <- testExpr[,-47097]
testExpr1 <- testExpr1[,-47097]
testExpr1 <- t(testExpr1)

gc()
calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr1,
              batchCorrect = 'eb',  #   "eb" for ComBat
              powerTransformPhenotype = F,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 20,
              printOutput = TRUE,
              removeLowVaringGenesFrom = 'homogenizeData' )


testPtype <- read.csv('./calcPhenotype_Output/DrugPredictions.csv', 
                      row.names = 1,check.names = F)
testPtype[1:4, 1:4]
dim(testPtype)
identical(colnames(testPtype),colnames(GDSC2_Res))
identical(rownames(testExpr),rownames(testPtype))

res <- read.csv("./calcPhenotype_Output/DrugPredictions.csv")
res <- as.data.frame(res)
gene <- as.data.frame(testPtype)
group <- as.data.frame(res[,-c(2:199)])
group$group <- as.factor(testExpr[,47098])
#gene$group <- as.factor(testExpr[,47098])
colnames(group) <- c("sample","group")
#gene$sample <- group$sample
gene <- t(gene)
gene <- as.data.frame(gene)
result <- NULL
for (n in 1:nrow(gene)) {
    gene_n <- data.frame(t(gene[n,]))
    gene_id <- names(gene_n)[1]
    names(gene_n)[1] <- 'gene'
    
    gene_n$sample <- rownames(gene_n)
    gene_n <- merge(gene_n, group,by= 'sample'  , all.x = TRUE)
    
    gene_n$group <- factor(gene_n$group)
    p_value <- wilcox.test(gene~group, gene_n)$p.value
    if (!is.na(p_value) & p_value < 0.05) {
      stat <- summaryBy(gene~group, gene_n, FUN = c(mean, median))
      result <- rbind(result, c(gene_id, as.character(stat[1,1]), stat[1,2], stat[1,3], 
                                as.character(stat[2,1]), stat[2,2], stat[2,3], p_value))
    }
}
result <- data.frame(result)
names(result) <- c('Drug', 'Group1', 'mean1', 'median1', 'Group2', 'mean2', 'median2', 'p_value')
result$p_adjust <- p.adjust(result$p_value, method = 'BH') 

#durg3
result <- result[order(result$p_adjust),]
result$Drug[grepl("Docet", result$Drug)]
#"Docetaxel_1007" "Docetaxel_1819"
result$Drug[grepl("Pacli", result$Drug)]
#"Paclitaxel_1080"
which(result$Drug%in%c("Docetaxel_1007","Docetaxel_1819","Paclitaxel_1080"))
#62,75,98
result3=as.data.frame(result[c(62,75,98),])
result3$p_adjust
#5.567698e-11    3.747709e-10    5.671554e-09

res <- read.csv("./calcPhenotype_Output/DrugPredictions.csv")
which(colnames(res)%in%c("Docetaxel_1007","Docetaxel_1819","Paclitaxel_1080"))
#6  37 141
res3=as.data.frame(res[,c(1,6,37,141)])
sample_group <- group$group
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggsci)
res3 %>% 
  select(1:4) %>% 
  bind_cols(sample_group = sample_group) %>% 
  pivot_longer(2:4,names_to = "drugs",values_to = "logIC50") %>% 
  ggplot(., aes(sample_group,logIC50))+
  geom_boxplot(aes(fill=sample_group))+
  scale_fill_jama()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none")+
  facet_wrap(vars(drugs),scales = "free_y",nrow = 1)+
  stat_compare_means()

#drug30
log(50)
#3.912
str(result1)
result$mean1 <- exp(as.numeric(result$mean1))
result$mean2 <- exp(as.numeric(result$mean2))
result$median1 <- exp(as.numeric(result$median1))
result$median2 <- exp(as.numeric(result$median2))
str(result)
result1<- subset(result, p_adjust < 0.01)
result1$"fold change" <- result1[, "mean2"]/result1[, "mean1"]
write.table(result1, 'brca221_5_durg.wilcox.txt', sep = '\t', row.names = FALSE, quote = FALSE)
drugre <- read.table('brca221_5_durg.wilcox.txt', sep = '\t')
#result2<- subset(result, mean1 < 3.5 & mean1 < mean2)
#result2$Drug
#"SB505124_1194"   "AZD2014_1441"    "Uprosertib_1553"
result1 <- result1[order(result1$"fold change",decreasing = T),][1:30,]
result1 <- result1[order(result1$`fold change`), ]
original_order <- unique(result1$Drug)
ggplot(result1, aes(x = result1$`fold change`, y = factor(result1$Drug, levels = original_order), size = result1$mean1, color = result1$p_adjust)) +
  geom_point(alpha = 0.8) +
  scale_size(range = c(2, 8)) +
  scale_color_gradient(low = "#2171b5", high = "#D95319") + 
  labs(x = "Fold Change", y = NULL, size = "IC50", color = "p.adgust") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

#end

####Estimate####
library(utils)
#install.packages("estimate", repos="http://r-forge.r-project.org", dependencies=TRUE)
library(estimate)
## start analysis
exp.file <- readRDS("exp221.rds")
str(exp.file)
exp.file <-  as.data.frame(exp.file)
for (col in names(exp.file)) {
  if (!is.factor(exp.file[[col]])) {
    exp.file[[col]] <- as.numeric(exp.file[[col]])
  }
}
exp.file <- log2(exp.file+1)
range(exp.file)
in.gct.file = "ESTIMATE_input.gct"
# format input file into GCT format
outputGCT(exp.file, in.gct.file)
write.table(exp.file,"expEST221.txt",sep="\t",quote = F)
filterCommonGenes(input.f = "expEST221.txt",
                  output.f = in.gct.file,
                  id = "GeneSymbol")
# "Merged dataset includes 10209 genes (203 mismatched)."

out.score.file = "ESTIMATE_score.gct"
estimateScore(in.gct.file, out.score.file, platform = "illumina")
scores <- read.table(out.score.file,skip = 2,header = T)
rownames(scores) <- scores[,1]
scores <- t(scores[,3:ncol(scores)])
ESTIMATE_score <- as.data.frame(scores)
str(ESTIMATE_score)
ESTIMATE_score$TumorPurity = cos(0.6049872018+0.0001467884 * ESTIMATE_score[,3])

exp.file <- t(exp.file)
exp.file <-  as.data.frame(exp.file)
str(exp.file)
nomogram <- readRDS("BRCA30_nomogram.RDS")
nom <- nomogram
library(nomogramFormula)
results<-formula_rd(nomogram=nom)
exp.file$points<-points_cal(formula = results$formula,rd=exp.file)
pre<-exp.file$points
#meta_data221
meta_data<- readRDS("meta_data221.rds")
table(meta_data$Group)
#pCR  RD
#80 141
round(80/(80+141),2)
#0.36
exp.file$Group<- ifelse(exp.file$points>quantile(exp.file$points,.36),"High","Low")
table(exp.file$Group)
identical(rownames(exp.file),rownames(ESTIMATE_score))
ESTIMATE_score$group = exp.file$Group

# boxplot
library(ggplot2)
library(ggpubr)
library(ggsci)
ggplot(ESTIMATE_score, aes(x = group, y = TumorPurity, fill = group)) +
geom_boxplot(position = position_dodge(0.8)) +
geom_point(position = position_jitterdodge()) +
scale_fill_nejm() +
labs(x = "", y = 'ESTIMATE Tumor Purity') +
stat_compare_means() +
stat_compare_means(comparisons = combn(unique(ESTIMATE_score$group), 2, simplify =FALSE))+
theme_bw(base_size = 11)+theme(legend.position = "none")

ggplot(ESTIMATE_score, aes(x = group, y = ImmuneScore, fill = group)) +
geom_boxplot(position = position_dodge(0.8)) +
geom_point(position = position_jitterdodge()) +
scale_fill_nejm() +
labs(x = "", 
     y = 'ESTIMATE Immune Score') +
  stat_compare_means() +
  stat_compare_means(comparisons = combn(unique(ESTIMATE_score$group), 2, 
                                       simplify = FALSE)) +
  theme_bw(base_size = 11)+theme(legend.position = "none")

####IPS####
library(IOBR)
library(tidyverse)
library(ggpubr)

exp <- readRDS("exp221.rds")
str(exp)
exp <-  as.data.frame(exp)
for (col in names(exp)) {
  if (!is.factor(exp[[col]])) {
    exp[[col]] <- as.numeric(exp[[col]])
  }
}
exp <- log2(exp+1)
range(exp)
exp_in <- as.matrix(exp)
ips<-deconvo_tme(eset = exp_in, method = "ips", plot= FALSE)
meta_data <- readRDS("meta_data221_predict.rds")
table(meta_data$Grouppr)
identical(ips$ID,rownames(meta_data))
ips$group <- meta_data$Grouppr
#ESTIMATE_score$group = exp.file$Group
p1<- ggplot(ips,aes(group,MHC_IPS))+
  geom_violin(aes(fill=group),cex=1.2)+  
  scale_fill_manual(values = c('#FB5554','#579ABB'))+
  geom_boxplot(width=0.1,cex=1.2)+
  theme_classic(base_size = 11)+
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'none',axis.title.x = element_blank(),)+
  stat_compare_means(comparisons = combn(unique(ips$group), 
                                         2, simplify =FALSE))
p2<- ggplot(ips,aes(group,EC_IPS))+
  geom_violin(aes(fill=group),cex=1.2)+  
  scale_fill_manual(values = c('#FB5554','#579ABB'))+
  geom_boxplot(width=0.1,cex=1.2)+
  theme_classic(base_size = 11)+
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'none',axis.title.x = element_blank(),)+
  stat_compare_means(comparisons = combn(unique(ips$group), 
                                         2, simplify =FALSE))
p3<- ggplot(ips,aes(group,SC_IPS))+
  geom_violin(aes(fill=group),cex=1.2)+  
  scale_fill_manual(values = c('#FB5554','#579ABB'))+
  geom_boxplot(width=0.1,cex=1.2)+
  theme_classic(base_size = 11)+
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'none',axis.title.x = element_blank(),)+
  stat_compare_means(comparisons = combn(unique(ips$group), 
                                         2, simplify =FALSE))
p4<- ggplot(ips,aes(group,CP_IPS))+
  geom_violin(aes(fill=group),cex=1.2)+  
  scale_fill_manual(values = c('#FB5554','#579ABB'))+
  geom_boxplot(width=0.1,cex=1.2)+
  theme_classic(base_size = 11)+
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'none',axis.title.x = element_blank(),)+
  stat_compare_means(comparisons = combn(unique(ips$group), 
                                         2, simplify =FALSE))
p5<- ggplot(ips,aes(group,IPS_IPS))+
  geom_violin(aes(fill=group),cex=1.2)+  
  scale_fill_manual(values = c('#FB5554','#579ABB'))+
  geom_boxplot(width=0.1,cex=1.2)+
  theme_classic(base_size = 11)+
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'none',axis.title.x = element_blank(),)+
  stat_compare_means(comparisons = combn(unique(ips$group), 
                                         2, simplify =FALSE))
p5
ggsave("brca221_IPS.png",width = 8,height = 6)
combined_plot <- p1 + p2 + p3 + p4
combined_plot
ggsave("brca221_IPS_细分.png",width = 9.33,height = 7)
####TIDE_log2####
rt <- readRDS("exp221.rds")
str(rt)
rt <-  as.data.frame(rt)
for (col in names(rt)) {
  if (!is.factor(rt[[col]])) {
    rt[[col]] <- as.numeric(rt[[col]])
  }
}
rt <- log2(rt+1)
range(rt)
Expr <- rt
Expr <-t(apply(Expr,1,function(x)x-(mean(x))))
write.table(data.frame("gene symbol"=rownames(Expr),Expr),file="brca221(TIDEmatix_log2).txt",
            sep = "\t",quote = F,row.names = F)
#http://tide.dfci.harvard.edu/
#read "brca221_TIDE_results.csv" file
output <- read.csv("brca221_TIDE_resultslog2.csv")
tide <- output[,c(1,4)]
#rfGenes=tide[order(tide[,"MeanDecreaseGini"],decreasing = TRUE),]
tide1 <- tide[order(tide$Patient),]
meta_data <- readRDS("meta_data221_predict.rds")
identical(tide1$Patient,rownames(meta_data))
tide1$group <- meta_data$Grouppr
library(ggplot2)
library(ggpubr)
ggplot(tide1,aes(group,TIDE))+
  geom_violin(aes(fill=group),cex=1.2)+  
  scale_fill_manual(values = c('#FB5554','#579ABB'))+
  geom_boxplot(width=0.1,cex=1.2)+
  theme_classic(base_size = 11)+
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'none',axis.title.x = element_blank())+
  labs(y="TIDE score")+
  stat_compare_means(comparisons = combn(unique(tide1$group), 
                                         2, simplify =FALSE))
ggsave("brca221_TIDElog2.png",width = 4,height = 6)

###Predict immunotherapy###
TIDE <- output
meta <- meta_data
TIDE <- TIDE[order(TIDE$Patient),]
row.names(TIDE) <- TIDE$Patient
identical(rownames(meta),rownames(TIDE))
s <- intersect(rownames(meta),rownames(TIDE))
meta <- meta[s,]
TIDE <- TIDE[s,]
meta_TIDE <- cbind(meta,TIDE)
#write.csv(meta_TIDE,"meta_TIDE.csv")
res <- meta_TIDE
table(res$Responder,res$Grouppr)
f = fisher.test(table(res$Grouppr,res$Responder))

library(dplyr)
dat=count(res,Grouppr,Responder)
dat=dat%>%group_by(Grouppr)%>%
  summarise(Responder=Responder,n=n/sum(n))
dat$Responder=factor(dat$Responder,levels=c("False","True"))
dat

library(ggplot2)
p2=ggplot(data=dat)+
  geom_bar(aes(x=Grouppr,y=n,
               fill=Responder),
           stat="identity")+
  scale_fill_manual(values=c("#e04030","#6cb8d2"))+
  geom_label(aes(x=Grouppr,y=n,
                 label=scales::percent(n),
                 fill=Responder),
             color="white",
             size=4,label.size=0,
             show.legend = FALSE,
             position=position_fill(vjust=0.5))+
  ylab("Percentage")+
  theme_minimal()+
  guides(fill = guide_legend(title = "Responder"))  # 仅保留一个图例
p2
library(patchwork)
p2+plot_layout(widths=c(3,2),guides="collect")&theme(axis.text = element_text(color = 'black'),
                                                     axis.title.x = element_blank())
ggsave('brca221_immunpredictlog2.png',width = 5.5,height = 6)
str(res)


install.packages("tinyarray")
library(tinyarray)
res$Grouppr <- factor(res$Grouppr,levels = c("High","Low"))
colnames(res)
dat <- t(res[,c(17,8,10:12,14:15,18)])
head(dat)[1:4,1:4]
draw_boxplot(dat,res$Grouppr)+
  facet_wrap(~rows,scales ="free",nrow = 2) +
  scale_fill_manual(values = c("High" = "#e04030", "Low" = "#6cb8d2")) 
ggsave("brca221_TIDE_xifen.png",width = 14,height = 9.33)


####correlation analysis####
td <- readRDS("exp221.rds")
td <- t(td)
td <-  as.data.frame(td)
for (col in names(td)) {
  if (!is.factor(td[[col]])) {
    td[[col]] <- as.numeric(td[[col]])
  }
}
td <- log2(td+1)
range(td)

nomogram <- readRDS("BRCA30_nomogram.RDS")
nom <- nomogram
library(nomogramFormula)
results<-formula_rd(nomogram=nom)
td$points<-points_cal(formula = results$formula,rd=td)
pre<-td$points
which(colnames(td)%in%c("PDCD1","CD274","CTLA4","GZMB","LAG3"))
#which(colnames(td)%in%c("NCOA4","FTH1"))
colnames(td[,c(1848,2315,5015,11044,16643)])
#colnames(td[,c(12234,41079)])
#"LAG3"   "GZMB"  "CD274"-PDL1   "CTLA4"     "PDCD1"-PD1

scores <- td[,"points"]
gene1 <- td[,"PDCD1"]
gene2 <- td[,"CD274"]
gene3 <- td[,"CTLA4"]
gene4 <- td[,"GZMB"]
gene5 <- td[,"LAG3"]
#gene6 <- td[,"NCOA4"]
#gene7 <- td[,"FTH1"]
range(gene3)
cor.test (scores, gene1, method="spearman")
#cor.test (scores, gene6, method="spearman")
#cor.test (scores, gene7, method="spearman")
library(ggpubr)
p1<- ggscatter(td, x = "points", y = "PDCD1", 
          color = "blue1",fill = "lightgray",
          add = "reg.line", conf.int = TRUE, 
          add.params = list( color = "black",fill = "lightgray",fill = "lightgray"),
          cor.coef = T,
          cor.coeff.args = list(),
          cor.method = "spearman",
          cor.coef.coord = c(93, 4),
          cor.coef.size = 5,xlab = "Risk score",
          ylab = "PD-1")+
  coord_cartesian(ylim = c(0, 4))
p1
p2<- ggscatter(td, x = "points", y = "CD274", 
               color = "blue1",fill = "lightgray",
               add = "reg.line", conf.int = TRUE, 
               add.params = list( color = "black",fill = "lightgray",fill = "lightgray"),
               cor.coef = T,
               cor.coeff.args = list(),
               cor.method = "spearman",
               cor.coef.coord = c(93, 5),
               cor.coef.size = 5,xlab = "Risk score",
               ylab = "PD-L1")+
  coord_cartesian(ylim = c(0, 5))
p2
p3<- ggscatter(td, x = "points", y = "CTLA4", 
               color = "blue1",fill = "lightgray",
               add = "reg.line", conf.int = TRUE, 
               add.params = list( color = "black",fill = "lightgray",fill = "lightgray"),
               cor.coef = T,
               cor.coeff.args = list(),
               cor.method = "spearman",
               cor.coef.coord = c(93, 3),
               cor.coef.size = 5,xlab = "Risk score",
               ylab = "CTLA4")+
  coord_cartesian(ylim = c(0, 3))
p3
p4<- ggscatter(td, x = "points", y = "GZMB", 
               color = "blue1",fill = "lightgray",
               add = "reg.line", conf.int = TRUE, 
               add.params = list( color = "black",fill = "lightgray",fill = "lightgray"),
               cor.coef = T,
               cor.coeff.args = list(),
               cor.method = "spearman",
               cor.coef.coord = c(93, 8),
               cor.coef.size = 5,xlab = "Risk score",
               ylab = "GZMB")+
  coord_cartesian(ylim = c(0, 8))
p4
p5<- ggscatter(td, x = "points", y = "LAG3", 
               color = "blue1",fill = "lightgray",
               add = "reg.line", conf.int = TRUE, 
               add.params = list( color = "black",fill = "lightgray",fill = "lightgray"),
               cor.coef = T,
               cor.coeff.args = list(),
               cor.method = "spearman",
               cor.coef.coord = c(93, 4),
               cor.coef.size = 5,xlab = "Risk score",
               ylab = "LAG3")+
  coord_cartesian(ylim = c(0, 4))
p5
combined_plot <- ggarrange(p1, p2, p3, p4, p5, nrow = 1)
combined_plot


####tumor-subtypes_pROC_GSE163882####
###meta_data221###
pdata <- read.table("pdata221.txt",sep = "\t")
colnames(pdata)
table(pdata$'estrogen.receptor.status.ch1')
table(pdata$'her2.receptor.status.ch1')
table(pdata$'progesterone.receptor.status.ch1')
specific_char <- "P"

f_ER_P_HER_N <- pdata[!grepl(specific_char, pdata$her2.receptor.status.ch1) 
                       &grepl(specific_char, pdata$estrogen.receptor.status.ch1) , ]
f_HER_P <- pdata[grepl(specific_char, pdata$her2.receptor.status.ch1), ]
f_ER_N_HER_N <- pdata[!grepl(specific_char, pdata$her2.receptor.status.ch1) 
                      &!grepl(specific_char, pdata$estrogen.receptor.status.ch1) , ]

ft_TNBC <- pdata[!grepl(specific_char, pdata$her2.receptor.status.ch1) #TNBC
                 &!grepl(specific_char, pdata$estrogen.receptor.status.ch1) 
                 &!grepl(specific_char, pdata$progesterone.receptor.status.ch1), ]
f_ER_N_HER_N == ft_TNBC

#meta_data68ERP_HER2N
meta_data2<- f_ER_P_HER_N$'response.to.nac.ch1'
meta_data1<- f_ER_P_HER_N$'description'
meta_data = cbind(meta_data1,meta_data2)
colnames(meta_data)=c("Sample","Group")
meta_data <- as.data.frame(meta_data)
meta_data <- meta_data[order(meta_data$Sample),]
rownames(meta_data) <- meta_data$Sample
meta_dataER <- meta_data
saveRDS(meta_data,"meta_data68ERP_HER2N.rds")
#meta_data63HERP
meta_data2<- f_HER_P$'response.to.nac.ch1'
meta_data1<- f_HER_P$'description'
meta_data = cbind(meta_data1,meta_data2)
colnames(meta_data)=c("Sample","Group")
meta_data <- as.data.frame(meta_data)
meta_data <- meta_data[order(meta_data$Sample),]
rownames(meta_data) <- meta_data$Sample
meta_dataHE <- meta_data
saveRDS(meta_data,"meta_data63HERP.rds")
#meta_data90TN
meta_data2<- f_ER_N_HER_N$'response.to.nac.ch1'
meta_data1<- f_ER_N_HER_N$'description'
meta_data = cbind(meta_data1,meta_data2)
colnames(meta_data)=c("Sample","Group")
meta_data <- as.data.frame(meta_data)
meta_data <- meta_data[order(meta_data$Sample),]
rownames(meta_data) <- meta_data$Sample
meta_dataTN <- meta_data
saveRDS(meta_data,"meta_data90TN.rds")

###exp221###
exp <- read.table("exp-GSE882.txt",sep = "\t")
exp <- t(exp)
colnames(exp) <- exp[1,]
exp <- exp[-1,]
exp <- as.data.frame(exp)
exp <- exp[order(exp$V1),]
rownames(exp) <- exp[,1]
exp <- exp[,-1]
exp <- t(exp)
exp <- as.data.frame(exp)

exp <- t(exp)
exp <- as.data.frame(exp)
exp$sample <- rownames(exp)
#exp68ER
common_characters <- intersect(exp$sample, meta_dataER$Sample)
expER <- exp[exp$sample %in% common_characters,]
expER <- expER[,-47097]
expER <- as.data.frame(t(expER))
saveRDS(expER,"exp68ER.rds")
#exp63HE
common_characters <- intersect(exp$sample, meta_dataHE$Sample)
expHE <- exp[exp$sample %in% common_characters,]
expHE <- expHE[,-47097]
expHE <- as.data.frame(t(expHE))
saveRDS(expHE,"exp63HE.rds")
#exp90TN
common_characters <- intersect(exp$sample, meta_dataTN$Sample)
expTN <- exp[exp$sample %in% common_characters,]
expTN <- expTN[,-47097]
expTN <- as.data.frame(t(expTN))
saveRDS(expTN,"exp90TN.rds")

###pROC_ER###
#exp68
exp1 <- readRDS("exp68ER.rds")
which(rownames(exp1)%in%c("KIAA1467","CCND1","CCNA2","ESR1","CD52"))
rownames(exp1[c(1627,1943,3742,8679,12580),])
exp5=as.data.frame(exp1[c(1627,1943,3742,8679,12580),])
str(exp5)
for (col in names(exp5)) {
  if (!is.factor(exp5[[col]])) {
    exp5[[col]] <- as.numeric(exp5[[col]])
  }
}
data1 <-  as.data.frame(exp5)
data1 <- log2(data1+1)
range(data1)
#meta_data68
meta_data1 <- readRDS("meta_data68ERP_HER2N.rds")
table(meta_data1$Group)
data1 <- t(data1)
data1 <- as.data.frame(data1)
data1$Sample<- rownames(data1)
data1 <- merge(meta_data1,data1,by="Sample")
rownames(data1) <- data1[,1]
data1 <- data1[,-1]
str(data1)
#data$Group <- as.factor(data$Group)
Mydata1 <- data1
nomogram <- readRDS("BRCA30_nomogram.RDS")
nom <- nomogram
library(nomogramFormula)
results<-formula_rd(nomogram=nom)
Mydata1$points<-points_cal(formula = results$formula,rd=Mydata1)
pre<-Mydata1$points
library(pROC)
plot.roc(Mydata1$Group, pre,
         main="ROC Curve", percent=TRUE,
         print.auc=TRUE,
         ci=TRUE,)
library(PRROC)
library(tidyverse)
pre <- as.numeric(pre)
weight <- ifelse(Mydata1$Group == "pCR", 0,1)
pr <- pr.curve(scores.class0 = pre, weights.class0 = weight,
               curve=TRUE,rand.compute = TRUE)
plot(pr, rand.plot = TRUE, auc.main=T,legend=F,lwd = 2,cex.main=1.2,col="black")

###pROC_HER###
#exp63
exp1 <- readRDS("exp63HE.rds")
which(rownames(exp1)%in%c("KIAA1467","CCND1","CCNA2","ESR1","CD52"))
rownames(exp1[c(1627,1943,3742,8679,12580),])
exp5=as.data.frame(exp1[c(1627,1943,3742,8679,12580),])
str(exp5)
for (col in names(exp5)) {
  if (!is.factor(exp5[[col]])) {
    exp5[[col]] <- as.numeric(exp5[[col]])
  }
}
data1 <-  as.data.frame(exp5)
data1 <- log2(data1+1)
range(data1)
#meta_data63
meta_data1 <- readRDS("meta_data63HERP.rds")
table(meta_data1$Group)
data1 <- t(data1)
data1 <- as.data.frame(data1)
data1$Sample<- rownames(data1)
data1 <- merge(meta_data1,data1,by="Sample")
rownames(data1) <- data1[,1]
data1 <- data1[,-1]
str(data1)
#data$Group <- as.factor(data$Group)
Mydata1 <- data1
nomogram <- readRDS("BRCA30_nomogram.RDS")
nom <- nomogram
library(nomogramFormula)
results<-formula_rd(nomogram=nom)
Mydata1$points<-points_cal(formula = results$formula,rd=Mydata1)
pre<-Mydata1$points
library(pROC)
plot.roc(Mydata1$Group, pre,
         main="ROC Curve", percent=TRUE,
         print.auc=TRUE,
         ci=TRUE,)
library(PRROC)
library(tidyverse)
pre <- as.numeric(pre)
weight <- ifelse(Mydata1$Group == "pCR", 0,1)
pr <- pr.curve(scores.class0 = pre, weights.class0 = weight,
               curve=TRUE,rand.compute = TRUE)
plot(pr, rand.plot = TRUE, auc.main=T,legend=F,lwd = 2,cex.main=1.2,col="black")

###pROC_TNBC###
#exp90
exp1 <- readRDS("exp90TN.rds")
which(rownames(exp1)%in%c("KIAA1467","CCND1","CCNA2","ESR1","CD52"))
rownames(exp1[c(1627,1943,3742,8679,12580),])
exp5=as.data.frame(exp1[c(1627,1943,3742,8679,12580),])
str(exp5)
for (col in names(exp5)) {
  if (!is.factor(exp5[[col]])) {
    exp5[[col]] <- as.numeric(exp5[[col]])
  }
}
data1 <-  as.data.frame(exp5)
data1 <- log2(data1+1)
range(data1)
#meta_data90
meta_data1 <- readRDS("meta_data90TN.rds")
table(meta_data1$Group)
data1 <- t(data1)
data1 <- as.data.frame(data1)
data1$Sample<- rownames(data1)
data1 <- merge(meta_data1,data1,by="Sample")
rownames(data1) <- data1[,1]
data1 <- data1[,-1]
str(data1)
#data$Group <- as.factor(data$Group)
Mydata1 <- data1
nomogram <- readRDS("BRCA30_nomogram.RDS")
nom <- nomogram
library(nomogramFormula)
results<-formula_rd(nomogram=nom)
Mydata1$points<-points_cal(formula = results$formula,rd=Mydata1)
pre<-Mydata1$points
library(pROC)
plot.roc(Mydata1$Group, pre,
         main="ROC Curve", percent=TRUE,
         print.auc=TRUE,
         ci=TRUE,)
library(PRROC)
library(tidyverse)
pre <- as.numeric(pre)
weight <- ifelse(Mydata1$Group == "pCR", 0,1)
pr <- pr.curve(scores.class0 = pre, weights.class0 = weight,
               curve=TRUE,rand.compute = TRUE)
plot(pr, rand.plot = TRUE, auc.main=T,legend=F,lwd = 2,cex.main=1.2,col="black")

####tumor-subtypes_pROC_GSE20271####
###meta_data178###
pdata <- read.table("pdata178.txt",sep = "\t",header = T,row.names = 1)
colnames(pdata)
#characteristics_ch1.12
#GSM508058 don't have "ER% positive",use "NA"replace
#characteristics_ch1.11
which(colnames(pdata)%in%c("characteristics_ch1.11"))
which(colnames(pdata)%in%c("characteristics_ch1.24"))
which(rownames(pdata)%in%c("GSM508058"))
row <- 47
pdata[row, 22:34] <- pdata[row, 21:33]
pdata[row, 21] <- NA

table(pdata$'characteristics_ch1.12')#ER
table(pdata$'characteristics_ch1.14')#PR
table(pdata$'characteristics_ch1.15')#HER
specific_char <- "P"

f_ER_P_HER_N <- pdata[!grepl(specific_char, pdata$characteristics_ch1.15) 
                      &grepl(specific_char, pdata$characteristics_ch1.12) , ]
f_HER_P <- pdata[grepl(specific_char, pdata$characteristics_ch1.15), ]
f_ER_N_HER_N <- pdata[!grepl(specific_char, pdata$characteristics_ch1.12) 
                      &!grepl(specific_char, pdata$characteristics_ch1.15) , ]

ft_TNBC <- pdata[!grepl(specific_char, pdata$characteristics_ch1.12) #TNBC
                 &!grepl(specific_char, pdata$characteristics_ch1.14) 
                 &!grepl(specific_char, pdata$characteristics_ch1.15), ]
#meta_data89ERP_HER2N
library(tidyverse)
group_list<-str_split(f_ER_P_HER_N$characteristics_ch1,' ',simplify = T)[,4]
table(group_list)
meta_data2 <- group_list
meta_data1<- f_ER_P_HER_N$geo_accession
meta_data = as.data.frame(cbind(meta_data1,meta_data2))
colnames(meta_data)=c("Sample","Group")
rownames(meta_data) <- meta_data$Sample
meta_dataER <- meta_data
saveRDS(meta_data,"meta_data89_ERP_HERN_178.rds")
#meta_data26HERP
group_list<-str_split(f_HER_P$characteristics_ch1,' ',simplify = T)[,4]
table(group_list)
meta_data2 <- group_list
meta_data1<- f_HER_P$geo_accession
meta_data = as.data.frame(cbind(meta_data1,meta_data2))
colnames(meta_data)=c("Sample","Group")
rownames(meta_data) <- meta_data$Sample
meta_dataHE <- meta_data
saveRDS(meta_data,"meta_data26_HERP_178.rds")
#meta_data63ERN_HERN
group_list<-str_split(f_ER_N_HER_N$characteristics_ch1,' ',simplify = T)[,4]
table(group_list)
meta_data2 <- group_list
meta_data1<- f_ER_N_HER_N$geo_accession
meta_data = as.data.frame(cbind(meta_data1,meta_data2))
colnames(meta_data)=c("Sample","Group")
rownames(meta_data) <- meta_data$Sample
meta_dataTN <- meta_data
saveRDS(meta_data,"meta_data63_ERN_HERN_178.rds")
###exp178###
exp <- readRDS("exp178.rds")
exp <- t(exp)
exp <- as.data.frame(exp)
exp$sample <- rownames(exp)
#exp89ER
common_characters <- intersect(exp$sample, meta_dataER$Sample)
expER <- exp[exp$sample %in% common_characters,]
expER <- expER[,-13517]
expER <- as.data.frame(t(expER))
saveRDS(expER,"exp89ER_178.rds")
#exp26HE
common_characters <- intersect(exp$sample, meta_dataHE$Sample)
expHE <- exp[exp$sample %in% common_characters,]
expHE <- expHE[,-13517]
expHE <- as.data.frame(t(expHE))
saveRDS(expHE,"exp26HE_178.rds")
#exp63TN
common_characters <- intersect(exp$sample, meta_dataTN$Sample)
expTN <- exp[exp$sample %in% common_characters,]
expTN <- expTN[,-13517]
expTN <- as.data.frame(t(expTN))
saveRDS(expTN,"exp63TN_178.rds")

###pROC_ER###
#exp89
exp1 <- readRDS("exp89ER_178.rds")
which(rownames(exp1)%in%c("KIAA1467","CCND1","CCNA2","ESR1","CD52"))
rownames(exp1[c(2070,3041,3501,6403,8372),])
exp5=as.data.frame(exp1[c(2070,3041,3501,6403,8372),])
str(exp5)
data1 <-  as.data.frame(exp5)
data1 <- log2(data1+1)
range(data1)
#meta_data89
meta_data1 <- readRDS("meta_data89_ERP_HERN_178.rds")
data1 <- t(data1)
data1 <- as.data.frame(data1)
data1$Sample<- rownames(data1)
data1 <- merge(meta_data1,data1,by="Sample")
rownames(data1) <- data1[,1]
data1 <- data1[,-1]
str(data1)
#data$Group <- as.factor(data$Group)
Mydata1 <- data1
nomogram <- readRDS("BRCA30_nomogram.RDS")
nom <- nomogram
library(nomogramFormula)
results<-formula_rd(nomogram=nom)
Mydata1$points<-points_cal(formula = results$formula,rd=Mydata1)
pre<-Mydata1$points
library(pROC)
plot.roc(Mydata1$Group, pre,
         main="ROC Curve", percent=TRUE,
         print.auc=TRUE,
         ci=TRUE,)
library(PRROC)
library(tidyverse)
pre <- as.numeric(pre)
weight <- ifelse(Mydata1$Group == "pCR", 0,1)
pr <- pr.curve(scores.class0 = pre, weights.class0 = weight,
               curve=TRUE,rand.compute = TRUE)
plot(pr, rand.plot = TRUE, auc.main=T,legend=F,lwd = 2,cex.main=1.2,col="black")

###pROC_HER###
#exp26
exp1 <- readRDS("exp26HE_178.rds")
which(rownames(exp1)%in%c("KIAA1467","CCND1","CCNA2","ESR1","CD52"))
rownames(exp1[c(2070,3041,3501,6403,8372),])
exp5=as.data.frame(exp1[c(2070,3041,3501,6403,8372),])
str(exp5)
data1 <-  as.data.frame(exp5)
data1 <- log2(data1+1)
range(data1)
#meta_data26
meta_data1 <- readRDS("meta_data26_HERP_178.rds")
data1 <- t(data1)
data1 <- as.data.frame(data1)
data1$Sample<- rownames(data1)
data1 <- merge(meta_data1,data1,by="Sample")
rownames(data1) <- data1[,1]
data1 <- data1[,-1]
str(data1)
#data$Group <- as.factor(data$Group)
Mydata1 <- data1
nomogram <- readRDS("BRCA30_nomogram.RDS")
nom <- nomogram
library(nomogramFormula)
results<-formula_rd(nomogram=nom)
Mydata1$points<-points_cal(formula = results$formula,rd=Mydata1)
pre<-Mydata1$points
library(pROC)
plot.roc(Mydata1$Group, pre,
         main="ROC Curve", percent=TRUE,
         print.auc=TRUE,
         ci=TRUE,)
library(PRROC)
library(tidyverse)
pre <- as.numeric(pre)
weight <- ifelse(Mydata1$Group == "pCR", 0,1)
pr <- pr.curve(scores.class0 = pre, weights.class0 = weight,
               curve=TRUE,rand.compute = TRUE)
plot(pr, rand.plot = TRUE, auc.main=T,legend=F,lwd = 2,cex.main=1.2,col="black")

###pROC_TNBC###
#exp63
exp1 <- readRDS("exp63TN_178.rds")
which(rownames(exp1)%in%c("KIAA1467","CCND1","CCNA2","ESR1","CD52"))
rownames(exp1[c(2070,3041,3501,6403,8372),])
exp5=as.data.frame(exp1[c(2070,3041,3501,6403,8372),])
str(exp5)
data1 <-  as.data.frame(exp5)
data1 <- log2(data1+1)
range(data1)
#meta_data63
meta_data1 <- readRDS("meta_data63_ERN_HERN_178.rds")
data1 <- t(data1)
data1 <- as.data.frame(data1)
data1$Sample<- rownames(data1)
data1 <- merge(meta_data1,data1,by="Sample")
rownames(data1) <- data1[,1]
data1 <- data1[,-1]
str(data1)
#data$Group <- as.factor(data$Group)
Mydata1 <- data1
nomogram <- readRDS("BRCA30_nomogram.RDS")
nom <- nomogram
library(nomogramFormula)
results<-formula_rd(nomogram=nom)
Mydata1$points<-points_cal(formula = results$formula,rd=Mydata1)
pre<-Mydata1$points
library(pROC)
plot.roc(Mydata1$Group, pre,
         main="ROC Curve", percent=TRUE,
         print.auc=TRUE,
         ci=TRUE,)
library(PRROC)
library(tidyverse)
pre <- as.numeric(pre)
weight <- ifelse(Mydata1$Group == "pCR", 0,1)
pr <- pr.curve(scores.class0 = pre, weights.class0 = weight,
               curve=TRUE,rand.compute = TRUE)
plot(pr, rand.plot = TRUE, auc.main=T,legend=F,lwd = 2,cex.main=1.2,col="black")


####subtypes_pcr_information####
#GSE163882#
pdata <- read.table("pdata221.txt",sep = "\t")
estrogen.receptor.status<- pdata$'estrogen.receptor.status.ch1'
her2.receptor.status<- pdata$'her2.receptor.status.ch1'
progesterone.receptor.status<- pdata$'progesterone.receptor.status.ch1'
response.to.nac<- pdata$'response.to.nac.ch1'
description<- pdata$'description'
meta_data = cbind(description,response.to.nac,estrogen.receptor.status,progesterone.receptor.status
                  ,her2.receptor.status)
write.table(meta_data,file = "subtype_PCR_882.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

#20271#
pdata <- read.table("pdata178.txt",sep = "\t",header = T,row.names = 1)
#characteristics_ch1.12
#GSM508058 don't have "ER% positive",use "NA"replace
#characteristics_ch1.11
which(colnames(pdata)%in%c("characteristics_ch1.11"))
which(colnames(pdata)%in%c("characteristics_ch1.24"))
which(rownames(pdata)%in%c("GSM508058"))
row <- 47
pdata[row, 22:34] <- pdata[row, 21:33]
pdata[row, 21] <- NA
table(pdata$'characteristics_ch1.12')#ER
table(pdata$'characteristics_ch1.14')#PR
table(pdata$'characteristics_ch1.15')#HER
estrogen.receptor.status<- pdata$'characteristics_ch1.12'
her2.receptor.status<- pdata$'characteristics_ch1.15'
progesterone.receptor.status<- pdata$'characteristics_ch1.14'
library(tidyverse)
group_list<-str_split(pdata$characteristics_ch1,' ',simplify = T)[,4]
table(group_list)
response.to.nac <- group_list
description<- pdata$geo_accession

meta_data = cbind(description,response.to.nac,estrogen.receptor.status,progesterone.receptor.status
                  ,her2.receptor.status)
write.table(meta_data,file = "subtype_PCR_271.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

