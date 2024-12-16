####seurat####
setwd("BRCA")
library(Seurat)
library(dplyr)
library(patchwork)
#library(mindr)
library(Matrix)
#library(DT)

dir_name=c('CID3921','CID3946','CID4461','CID3586','CID3838','CID3941','CID3948','CID3963','CID4040','CID4066','CID4067','CID4290A','CID4398','CID4463','CID4465','CID4471',
           'CID4495','CID4513','CID4515','CID4523','CID4530N','CID4535','CID44041','CID44971','CID44991','CID45171')
datalist=list()
for (i in 1:length(dir_name)){
  dir.10x = paste0("GSE176078/",dir_name[i])
  my.data <- Read10X(data.dir = dir.10x,gene.column = 1) 
  datalist[[i]]=CreateSeuratObject(counts = my.data, project = dir_name[i], 
                                   min.cells = 3, min.features = 250)
}

folders=list.files('./GSE176078')
folders
sce.big <- merge(datalist[[1]],y = datalist[-1],
                 add.cell.ids = folders, 
                 project = "BRCA3")
sce.big
table(sce.big$orig.ident)
save(sce.big,file = 'sce.big.brca26-pbmc.Rdata')
readrds(file = "sce.big.brca26-pbmc.Rdata")
pbmc <- sce.big


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

#1
sce.big<- pbmc
sce.big <- ScaleData(object = sce.big, 
                     vars.to.regress = c('nCount_RNA'), 
                     model.use = 'linear', 
                     use.umi = FALSE)
sce.big <- FindVariableFeatures(object = sce.big, 
                                mean.function = ExpMean, 
                                dispersion.function = LogVMR, 
                                x.low.cutoff = 0.0125, 
                                x.high.cutoff = 4, 
                                y.cutoff = 0.5)
length(VariableFeatures(sce.big)) 
sce.big <- RunPCA(object = sce.big, pc.genes = VariableFeatures(sce.big))

sce.big <- JackStraw(sce.big, num.replicate = 100)
sce.big <- ScoreJackStraw(sce.big, dims = 1:20)
JackStrawPlot(sce.big, dims = 1:15)
#sce.big <- RunICA(sce.big )
sce.big <- RunTSNE(sce.big )
#sce.big <- RunUMAP(sce.big,dims = 1:10)
#VizPCA( sce.big, pcs.use = 1:2)
#DimPlot(object = sce.big, reduction = "pca") 
#DimPlot(object = sce.big, reduction = "ica")
DimPlot(object = sce.big, reduction = "tsne")

sce.big <- FindNeighbors(object = sce.big, dims = 1:20, verbose = FALSE) 
sce.big <- FindClusters(object = sce.big, resolution = 0.3,verbose = FALSE)
DimPlot(object = sce.big, reduction = "tsne",group.by = 'RNA_snn_res.0.3')
DimPlot(object = sce.big, reduction = "tsne",
        group.by = 'orig.ident')
table(sce.big$orig.ident,sce.big@meta.data$RNA_snn_res.0.3)

#2
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = VariableFeatures(object = pbmc),vars.to.regress = "percent.mt")
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Examine and visualize PCA results a few different ways
#print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
#DimPlot(pbmc, reduction = "pca") + NoLegend()
#ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:30, reduction = "pca")
pbmc <- FindClusters(pbmc, resolution = 0.5, cluster.name = "unintegrated_clusters")

pbmc <- RunUMAP(pbmc, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(pbmc, reduction = "umap.unintegrated", group.by = c("orig.ident", "seurat_clusters"))

#1
library(foreach)
library(doParallel)
detectCores()
cl.cores <- detectCores()
cl<- makeCluster(ceiling(cl.cores*2/3))
registerDoParallel(cl)

stopCluster(cl)
#2
library(future)
library(future.apply)
availableCores()
plan("multisession",workers=8)
nbrOfWorkers()
options(future.globals.maxSize=3000*1024^2)
#end

pbmc <- IntegrateLayers(
  object = pbmc, 
  method = RPCAIntegration,
  orig.reduction = "pca", 
  new.reduction = "integrated.rpca",
  verbose = FALSE)

saveRDS(pbmc, file = "pbmc-brca26.rds")
#load(file = "pbmc-brca26.rds")
#pbmc<-readRDS('pbmc-brca26.rds')


pbmc <- FindNeighbors(pbmc, reduction = "integrated.rpca", dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 1, cluster.name = "rpca_clusters")
pbmc <- RunUMAP(pbmc, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
DimPlot(pbmc,reduction = "umap.rpca",group.by = c("orig.ident", "rpca_clusters"),combine = FALSE,label.size = 2)
pbmc <- RunTSNE(pbmc, reduction = "integrated.rpca", dims = 1:30, reduction.name = "tsne.rpca")
DimPlot(pbmc,reduction = "tsne.rpca",group.by = c("orig.ident", "rpca_clusters"),combine = FALSE,label.size = 2)
#end

pbmc <- JoinLayers(pbmc)
cluster2.markers <- FindAllMarkers(pbmc, only.pos = TRUE, 
                                   min.pct = 0.25, 
                                   test.use = "wilcox", # wlicox,bimod,roc,Student`s t test,poisson,negbinom,LR,MAST,DEseq2
                                   logfc.threshold = 0.25)

saveRDS(cluster2.markers, "cluster2.markers.rds")
table(cluster2.markers$cluster)
#------->
sce.markers <- cluster2.markers



top10 <- sce.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

write.table(top10, file = "top10_brca26.csv", 
            quote = TRUE, sep = ",",row.names = F,
            col.names = TRUE)

pbmc1 <- pbmc
top10$gene
DoHeatmap(pbmc1, features = top10$gene) + NoLegend()
#ggsave("DoHeatmap1.png",path="/home/shuangshuang/R/Rstudio/03_1.MethyICIBERSORT/picture/picture_Transcriptome",width = 16, height = 12)

VlnPlot(pbmc1, features = top10$gene[1:10],pt.size=0)

DotPlot(pbmc1, features = unique(top10$gene)[1:20])+RotatedAxis()
ggsave("DotPlot1.pdf",,width = 10, height = 7)


#------------------------------

#-------------------------Annotation of cell taxa using SingleR--------------------
options(repos="http://mirrors.tuna.tsinghua.edu.cn/CRAN/")
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
BiocManager::install("pheatmap") 
install.packages('ExperimentHub')

library(SingleR)
library(celldex)
library(matrixStats)
#rm(list=ls())
#options(stringsAsFactors = F)
library(Seurat)
library(pheatmap)


sce_for_SingleR <- GetAssayData(pbmc1, layer="data")
#rownames(sce_for_SingleR)[1:10]
clusters=pbmc1@meta.data$seurat_clusters

#1
ref.data <- HumanPrimaryCellAtlasData(ensembl=F)
ref.data@colData
rownames(ref.data)
pred.humanImmu <- SingleR(test = sce_for_SingleR, ref = ref.data, labels = ref.data$label.main,
                          clusters = clusters, 
                          assay.type.test = "logcounts", assay.type.ref = "logcounts")
meta=pbmc1@meta.data 
table(meta$seurat_clusters)
table(pred.humanImmu$labels)
#MSC 
#14 
saveRDS(ref.data, "ref.data.rds")
pred.humanImmu$labels
#pred.humanImmu$first.labels
head(pred.humanImmu$scores)
SingleR::plotScoreHeatmap(pred.humanImmu)


summary(is.na(pred.humanImmu$pruned.labels))
#   Mode   FALSE 
#logical      14 
to.remove <- is.na(pred.humanImmu$pruned.labels)
table(Label=pred.humanImmu$labels, Removed=to.remove)



library(ExperimentHub)
tools::R_user_dir("ExperimentHub", which="cache")
list.files("/home/zz/.cache/R/ExperimentHub")
library(AnnotationHub)
snapshotDate()
ExperimentHub::snapshotDate()
#linux
#rm -rf dir/home/zz/.cache/R/ExperimentHub/
#------------
ref.data1 <- celldex::DatabaseImmuneCellExpressionData()
ref.data1 <- ref_Monaco
ref.data1@colData
#saveRDS(ref.data1,"ref_data1.rds")
load(file = "ref_data1.rds")
rownames(ref.data1)
pred.humanImmu1 <- SingleR(test = sce_for_SingleR, ref = ref.data1, labels = ref.data1$label.main,
                           clusters = clusters, 
                           assay.type.test = "logcounts", assay.type.ref = "logcounts")
table(pred.humanImmu1$labels)

ref.data2 <- ref_Hematopoietic
rownames(ref.data2)
pred.humanImmu2 <- SingleR(test = sce_for_SingleR, ref = ref.data2, labels = ref.data2$label.main,
                           clusters = clusters, 
                           assay.type.test = "logcounts", assay.type.ref = "logcounts")
table(pred.humanImmu2$labels)
SingleR::plotScoreHeatmap(pred.humanImmu2)
#Monocytes 
#14 
saveRDS(ref.data1, "ref.data1.rds")

pred.humanImmu1
pred.humanImmu1$labels
#pred.humanImmu1$first.labels
pred.humanImmu1$scores
SingleR::plotScoreHeatmap(pred.humanImmu1)


#-----------Cell type------------
cellType=data.frame(ClusterID=levels(pbmc1@meta.data$seurat_clusters),
                    ref.data=pred.humanImmu$labels,
                    ref.data1=pred.humanImmu1$labels,
                    ref.data2=pred.humanImmu2$labels)
head(cellType)
pbmc1@meta.data$singleR=cellType[match(clusters,cellType$ClusterID),'ref.data']
pbmc1@meta.data$singleR
DimPlot(pbmc1,reduction = "tsne.rpca",label=T, label.size = 5,group.by = 'singleR')
DimPlot(pbmc1,reduction = "umap.rpca",label=T, label.size = 5,group.by = 'singleR')
#DimPlot(pbmc1,reduction = "umap.rpca",label=T,split.by ='orig.ident',group.by = 'singleR')

#-------------------------------

#----------------------Cell annotation------------------------
new.cluster.ids <- c("T cells",   
                     "T cells",    
                     "Epithelial",  
                     "T cells",             
                     "Epithelial",         
                     "T cells",  
                     "Macrophage",           
                     "Epithelial",            
                     "Endothelial",
                     "B cells",    
                     "Epithelial",  
                     "Fibroblasts",             
                     "Pericytes",         
                     "Fibroblasts",  
                     "B cells",
                     "Epithelial",
                     "Monocyte",
                     "T cells",
                     "Pericytes",
                     "Endothelial",
                     "Epithelial",
                     "T cells",
                     "Endothelial",
                     "Basal cells",
                     "Macrophage",
                     "T cells",
                     "pDC",
                     "Endothelial")

names(new.cluster.ids) <- levels(pbmc1)
#new.cluster.ids
FeaturePlot(pbmc, features = c("KRT17", "EPCAM", "PECAM1", "COL1A1", "RGS5", "FCN1", "CD68", "CD3D",
                               "CD79A","CLEC4C"),reduction = "umap.rpca")
FeaturePlot(pbmc, features = c("CD19","CD79A","CD79B","MS4A1","CD37"),reduction = "umap.rpca")


Idents(pbmc1)
table(Idents(pbmc1))

pbmc1 <- RenameIdents(pbmc1, new.cluster.ids)
Idents(pbmc1)
table(Idents(pbmc1))
head(pbmc1@meta.data)
colnames(pbmc1@meta.data)
table(pbmc1@meta.data$singleR)

DimPlot(pbmc1, reduction = "umap.rpca", label = TRUE, pt.size = 1.2) + NoLegend()
DimPlot(pbmc1, reduction = "tsne", label = TRUE, pt.size = 1.2) + NoLegend()

saveRDS(pbmc1, file = "pbmc1_brca26_zhushi.rds")
#load(file = "pbmc1_brca26_zhushi.rds")
#pbmc1<-readRDS('pbmc1_brca26_zhushi.rds')

pbmc1@meta.data$cell_anno <- Idents(pbmc1)
table(pbmc1@meta.data$cell_anno)

rownames(pbmc1@assays$RNA@layers$data) <- pbmc1[["RNA"]]@features %>% rownames()

colnames(pbmc1@assays$RNA@layers$data) <- pbmc1[["RNA"]]@cells %>% rownames()

data <- as.data.frame(pbmc1@assays$RNA@layers$data)
dim(data)
range(data)
#[1] 18379  5189

phe <- pbmc1@meta.data
dim(phe)
#[1] 5189    7

phe1 <- t(phe)
dim(phe1)
#[1]    7   5189

all(colnames(phe1) %in% colnames(data))
#[1] TRUE

all(colnames(data) %in% colnames(phe1))
#[1] TRUE


loc = match(colnames(data),colnames(phe1))
loc

a1 <- rbind(phe1,data)
dim(a1)
#[1]  [1] 18386  5189

saveRDS(a1, file = "a1.rds")
#load(file = "a1.rds")
#a1<-readRDS('a1.rds')
a1 <- a1[-c(1:6),]
dim(a1)


a2 <- a1

head(colnames(a2))
write.csv(pbmc1@meta.data,file = "pbmc1_metadata.csv")



#-----------------------------------------------------

colnames(a2) <- phe$cell_anno
a2 <- a2[-1,]


#----------Number of cells for each cell type----------
table(pbmc1@meta.data$cell_anno)
length(which(colnames(a2)=="T_cells"))
#[1] 1741

length(which(colnames(a2)=="Fibroblasts"))
#[1] 147


#---------------------------------------

#-----------------------------Organize data-----------------
rm(list = ls())
gc()
setwd("BRCA")
library(dplyr)
library(stringr)



clusters <- as.data.frame(table(colnames(a2)))
#clusters <- clusters[!(clusters$Var1 == "Gene"), ]

dim(a2)
#1
a2[,99316] <- rownames(a2)
a2 <- a2[,c(99316,1:99315)]
colnames(a2)[1] <- "gene"

geneid = a2[, 1]
#rownames(a2) = a2[, 1]
a2 <- a2[,-1]
#2
a2[,99316] <- rownames(a2)
colnames(a2)[99316] <- "gene"
geneid = a2[, 99316]
a2 <- a2[,-99316]
geneid <- rownames(a2)
#end
saveRDS(a2, file = "a2.rds")
#write.csv(a2,file = "a2.csv",row.names = TRUE)


newdata = as.data.frame(NULL)
datai = select(a2, starts_with(as.character(clusters[1, 1]))) 
newdata = rbind(newdata, datai)
for (i in 2:length(clusters$Var1)) {
  datai = select(a2, starts_with(as.character(clusters[i, 1])))
  newdata = cbind(newdata, datai)
}
rownames(newdata) = geneid
data_sample_weight = data.frame()
weight = 3
num = floor(clusters$Freq / weight)
datai = select(a2, starts_with(as.character(clusters[1, 1])))
datai = sample(datai, num[1], replace = FALSE) 
data_sample_weight = rbind(data_sample_weight, datai)
for (i in 2:length(clusters$Var1)) {
  datai = select(a2, starts_with(as.character(clusters[i, 1])))
  datai = sample(datai, num[i], replace = FALSE)
  data_sample_weight = cbind.data.frame(data_sample_weight, datai)
}


data <- data_sample_weight
col1 <- data[, 1]
rownames(data) <- col1

col2 <- data[1,]
colnames(data) <- col2
data <- data[-1,-1]

#colnames(data_sample_weight) = geneid
#colnames(data)=str_split(colnames(data),'[.]',simplify=T)[,1]
#data_sample_weight <- cbind(geneid, data_sample_weight)

library(future)
library(future.apply)
availableCores()
plan("multisession",workers=8)
nbrOfWorkers()
options(future.globals.maxSize=3000*1024^2)

data <-  as.data.frame(a2)
str(data)
for (col in names(data)) {
  if (!is.factor(data[[col]])) {
    data[[col]] <- as.numeric(data[[col]])
  }
}
str(data)


#col_names <- gsub("[ ]","_", colnames(data))
#colnames(data) <- col_names

PID = c("Basal cells","Epithelial","Endothelial","Fibroblasts","Pericytes",
        "Monocyte","Macrophage","T cells","B cells","pDC")
new_data = matrix(NA, nrow(data), length(PID))
for(i in 1:length(PID)){
  rep_columns = select(data, contains(PID[i]))
  new_data[, i] = future_apply(rep_columns, 1, mean,na.rm = TRUE)
}

#class(data1)
#table(is.na(data1))
#mean(data1$`data$Endothelial`)

new_data = as.data.frame(new_data)
colnames(new_data) = PID
rownames(new_data) <- rownames(data)

saveRDS(new_data, file = "brca26_cellsgene.rds")

brca <- brca26_cellsgene
x<-data.frame(GeneSymbol=rownames(brca))
rownames(brca)<- NULL
brca <- cbind(x,brca)
write.table(brca,'brca26_cellsgene.txt',sep='\t',col.names=T,row.names = F)
brca <- read.table('brca26_cellsgene.txt',sep='\t')

#brca26_cellsgene is the final file uploaded to the official website, named brca26_cellsgene file

setwd()
####Draw single-cell proportion analysis####
install.packages("ggplot2")
library(dplyr)
library(reshape2)
library(ggpubr)
## 2.
#Organize brca,BRCA

#1 BRCA10
BRCA <- read.table("CIBERSORTx_Job1_Results.txt",sep="\t",header=T)
BRCA <- BRCA[order(BRCA$Mixture),]
rownames(BRCA) <- BRCA[,1]
saveRDS(BRCA,"CIBERSORTx-Results,rds")
#meta_data221.rds
table(meta_data$Group)
group_list <- meta_data$Group %>% 
  factor(.,levels = c("pCR","RD"))
table(group_list) 
TME_data <- BRCA
TME_data$group <- group_list
TME_data <- TME_data[,-1]
TME_data$sample <- row.names(TME_data)

# 2.2 
TME_New = melt(TME_data, id.vars = c("group","sample"))
## Using group, sample as id variables

colnames(TME_New)=c("Group","Sample","Celltype","Composition")
head(TME_New)
TME_New <- TME_New[1:2210,]

#2 brca22
#Upload the bulk file to the CIBERSORTx website, use the LM22 file provided by the website, and run it to obtain CIBERSORT-Results. txt
brca <- read.table("CIBERSORT-Results.txt",sep="\t",header=T)
brca <- brca[order(brca$Mixture),]
rownames(brca) <- brca[,1]
saveRDS(brca,"CIBERSORT-Results,rds")
table(group_list) 
TME_data <- brca
TME_data$group <- group_list
TME_data <- TME_data[,-1]
TME_data$sample <- row.names(TME_data)
# 2.2 
TME_New = melt(TME_data, id.vars = c("group","sample"))
## Using group, sample as id variables

colnames(TME_New)=c("Group","Sample","Celltype","Composition")
head(TME_New)
TME_New <- TME_New[1:4862,]
#end

# 3.3 
plot_order = TME_New[TME_New$Group=="RD",] %>% 
  group_by(Celltype) %>% 
  summarise(m = median(Composition)) %>% 
  arrange(desc(m)) %>% 
  pull(Celltype)
TME_New$Celltype = factor(TME_New$Celltype,levels = plot_order)
# 3.3 
if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"), 
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  ) }



box_TME <- ggplot(TME_New, aes(x = Celltype, y = Composition))+ 
  labs(y="Cell composition",x= NULL,title = "Cell composition")+  
  geom_boxplot(aes(fill = Group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
  scale_fill_manual(values = c("#1CB4B8", "#EB7369"))+
  theme_classic() + mytheme + 
  stat_compare_means(aes(group =  Group),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = T)
#box_TME + scale_y_continuous(breaks=seq(0, 1, 0.1))
box_TME;ggsave("BRCA-10cellcomposition.pdf",box_TME,height=12,width=32,unit="cm")
box_TME;ggsave("brca-22cellcomposition.pdf",box_TME,height=15,width=40,unit="cm")

####Single cell level correlation analysis####
###Organize data###
a2 <- readRDS("a2.rds")
a2 <- as.data.frame(a2)
str(a2)
for (col in names(a2)) {
  if (!is.factor(a2[[col]])) {
    a2[[col]] <- as.numeric(a2[[col]])
  }
}
str(a2)
range(a2)
library(dplyr)
Tcell <- a2 %>% select(contains("T cells"))
Bcell <- a2 %>% select(contains("B cells"))
pDC <- a2 %>% select(contains("pDC"))
monocyte <- a2 %>% select(contains("monocyte"))
macrophage <- a2 %>% select(contains("macrophage"))
Myeloid <- cbind(pDC, monocyte, macrophage)
str(Myeloid)
which(rownames(a2)%in%c("CD52"))
#338

####Tcell####
str(Tcell)
split_value_row <- Tcell[338, ]
split_value_row <- as.numeric(split_value_row)
median_value <- median(split_value_row)
group1 <- Tcell[, split_value_row <= median_value]
group2 <- Tcell[, split_value_row > median_value]

library(limma)
combined_matrix <- cbind(group2, group1)
group <- c(rep(2, ncol(group2)), rep(1, ncol(group1)))
design <- model.matrix(~ factor(group))
fit=lmFit(combined_matrix,design)
fit=eBayes(fit)
deg2=topTable(fit,coef=2,number = Inf)
write.table(deg2, file = "deg-Tcells.txt",sep = "\t",row.names = T,col.names = T)

#Tcell enrichment analysis#
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
deg <- read.table("deg-Tcells.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
deg$gene_sample <- rownames(deg)
P.value = 0.0000000001
deg1 <- subset(deg, (deg$P.Value < P.value)&(deg$logFC > 0))
deg2 <- subset(deg, (deg$P.Value < P.value)&(deg$logFC < 0))
#deg1  500  deg2  444
#Tcell_up genes#
DEG <- deg1
genelist <- bitr(DEG$gene_sample, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by=c("gene_sample"="SYMBOL"))
#GO
ego <- enrichGO(gene = DEG$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)
ego_res <- ego@result
saveRDS(ego,"go_ego_t_up.rds")
#t_up####
ego <- readRDS("go_ego_t_up.rds")
# BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(ggplot2)
ego1 <- ego
ego1@result$Description[1:18]
ego1@result <- ego1@result[-9, ]
ego1@result <- ego1@result[-10, ]
dotplot_obj <- dotplot(ego1, showCategory = 10, label_format = 200)+
  theme(axis.text.y = element_text(size = 23))+theme(legend.position = "top",legend.justification = c(1, 1),legend.key.width = unit(13, "mm"))+
  theme(plot.margin = unit(c(0, 0.5, 0, 0.5), "cm"))+
  theme(axis.text.x = element_text(size = 20))+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 48)) +
  theme(axis.title = element_text(size = 20),  
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20))+
  theme(legend.title = element_text(vjust = 1),
        legend.text = element_text(vjust = 0.5))
dotplot_obj
#+theme(legend.key.height = unit(50, "pt"))
ggsave("go_t_up.pdf", plot = dotplot_obj, device = "pdf", width = 9, height = 7.6)
#5.5 9.1
write.table(ego_res,"GO_Tcell-1*10-10_up.txt",sep = "\t")
save(ego,ego_res,file = "GO_Tcell-1*10-10_up.Rdata")
#Tcell_down_genes#
DEG <- deg2
genelist <- bitr(DEG$gene_sample, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by=c("gene_sample"="SYMBOL"))
#GO
ego <- enrichGO(gene = DEG$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)
ego_res <- ego@result
saveRDS(ego,"go_ego_t_down.rds")
#t_down####
ego <- readRDS("go_ego_t_down.rds")
dotplot_obj <- dotplot(ego, showCategory = 10, label_format = 200)+
  theme(axis.text.y = element_text(size = 23))+theme(legend.position = "top",legend.justification = c(1, 1),legend.key.width = unit(16.3, "mm"))+
  theme(plot.margin = unit(c(0, 0.5, 0, 0.5), "cm"))+
  theme(axis.text.x = element_text(size = 20))+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 45)) +
  theme(axis.title = element_text(size = 20),  
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 20))+
  theme(legend.title = element_text(vjust = 1),
        legend.text = element_text(vjust = 0.5))
dotplot_obj
#+theme(legend.key.height = unit(50, "pt"))
ggsave("go_t_down.pdf", plot = dotplot_obj, device = "pdf", width = 9, height = 7.6)
#5  9
write.table(ego_res,"GO_Tcell-1*10-10_down.txt",sep = "\t")
save(ego,ego_res,file = "GO_Tcell-1*10-10_down.Rdata")

####Bcell####
str(Bcell)
split_value_row <- Bcell[338, ]
split_value_row <- as.numeric(split_value_row)
median_value <- median(split_value_row)
group1 <- Bcell[, split_value_row <= median_value]
group2 <- Bcell[, split_value_row > median_value]

library(limma)
combined_matrix <- cbind(group2, group1)
group <- c(rep(2, ncol(group2)), rep(1, ncol(group1)))
design <- model.matrix(~ factor(group))
fit=lmFit(combined_matrix,design)
fit=eBayes(fit)
deg2=topTable(fit,coef=2,number = Inf)
write.table(deg2, file = "deg-Bcells.txt",sep = "\t",row.names = T,col.names = T)

#Bcell enrichment analysis#
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
deg <- read.table("deg-Bcells.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
deg$gene_sample <- rownames(deg)
P.value = 0.0000000001
deg1 <- subset(deg, (deg$P.Value < P.value)&(deg$logFC > 0))
deg2 <- subset(deg, (deg$P.Value < P.value)&(deg$logFC < 0))
#2541  496
#Bcell_up genes#
DEG <- deg1
genelist <- bitr(DEG$gene_sample, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by=c("gene_sample"="SYMBOL"))
#GO
ego <- enrichGO(gene = DEG$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)
ego_res <- ego@result
write.table(ego_res,"GO_Bcell-1*10-10_up.txt",sep = "\t")
save(ego,ego_res,file = "GO_Bcell-1*10-10_up.Rdata")
str(ego@result[["p.adjust"]])
saveRDS(ego,"go_ego_b_up.rds")
#b_up####
ego <- readRDS("go_ego_b_up.rds")
ego@result[["p.adjust"]] <- as.numeric(formatC(ego@result[["p.adjust"]], format = "e", digits = 0))
dotplot_obj <- dotplot(ego, showCategory = 10, label_format = 200)+
  theme(axis.text.y = element_text(size = 23))+theme(legend.position = "top",legend.justification = c(1, 1),legend.key.width = unit(10, "mm"))+
  theme(plot.margin = unit(c(0, 0.5, 0, 0.5), "cm"))+
  theme(axis.text.x = element_text(size = 20))+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 47)) +
  theme(axis.title = element_text(size = 20),
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 20))+
  theme(legend.title = element_text(vjust = 1),
        legend.text = element_text(vjust = 0.5))
dotplot_obj
#+theme(legend.key.height = unit(50, "pt"))
ggsave("go_b_up.pdf", plot = dotplot_obj, device = "pdf", width = 9, height = 7.6)
ggsave("go_b_up.pdf", plot = dotplot_obj, device = "pdf", width = 9.1, height = 5.5)
#Bcell_down genes#
DEG <- deg2
genelist <- bitr(DEG$gene_sample, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by=c("gene_sample"="SYMBOL"))
#GO
ego <- enrichGO(gene = DEG$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)
ego_res <- ego@result
saveRDS(ego,"go_ego_b_down.rds")
#b_down####
ego <- readRDS("go_ego_b_down.rds")
ego1 <- ego
library(stringr)
ego1@result$Description[1:10]
ego1@result$Description[10]<- "adaptive immune response based on somatic recombination of immune receptors"
dotplot_obj <- dotplot(ego1, showCategory = 10, label_format = 200)+
  theme(axis.text.y = element_text(size = 23))+theme(legend.position = "top",legend.justification = c(1, 1),legend.key.width = unit(13.2, "mm"))+
  theme(plot.margin = unit(c(0, 0.5, 0, 0.5), "cm"))+
  theme(axis.text.x = element_text(size = 20))+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 47)) +
  theme(axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))+
  theme(legend.title = element_text(vjust = 1),
        legend.text = element_text(vjust = 0.5))
dotplot_obj
#+theme(legend.key.height = unit(50, "pt"))
ggsave("go_b_down.pdf", plot = dotplot_obj, device = "pdf", width = 9, height = 7.6)
ggsave("go_b_down.pdf", plot = dotplot_obj, device = "pdf", width = 9.1, height = 5.5)
# dotplot(ego, showCategory = 15,label_format=200)+ theme(axis.text.y = element_text(size = 12))
write.table(ego_res,"GO_Bcell-1*10-10_down.txt",sep = "\t")
save(ego,ego_res,file = "GO_Bcell-1*10-10_down.Rdata")

####Myeloid####
str(Myeloid)
split_value_row <- Myeloid[338, ]
split_value_row <- as.numeric(split_value_row)
median_value <- median(split_value_row)
group1 <- Myeloid[, split_value_row <= median_value]
group2 <- Myeloid[, split_value_row > median_value]

library(limma)
combined_matrix <- cbind(group2, group1)
group <- c(rep(2, ncol(group2)), rep(1, ncol(group1)))
design <- model.matrix(~ factor(group))
fit=lmFit(combined_matrix,design)
fit=eBayes(fit)
deg2=topTable(fit,coef=2,number = Inf)
write.table(deg2, file = "deg-Myeloid.txt",sep = "\t",row.names = T,col.names = T)

#Myeloid enrichment analysis#

library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
deg <- read.table("deg-Myeloid.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
deg$gene_sample <- rownames(deg)
P.value = 0.0000000001
deg1 <- subset(deg, (deg$P.Value < P.value)&(deg$logFC > 0))
deg2 <- subset(deg, (deg$P.Value < P.value)&(deg$logFC < 0))
#956  376
#Myeloid_up genes#
DEG <- deg1
genelist <- bitr(DEG$gene_sample, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by=c("gene_sample"="SYMBOL"))
#GO
ego <- enrichGO(gene = DEG$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)
ego_res <- ego@result
saveRDS(ego,"go_ego_myeloid_up.rds")
#m_up####
ego <- readRDS("go_ego_myeloid_up.rds")
ego@result[["p.adjust"]] <- as.numeric(formatC(ego@result[["p.adjust"]], format = "e", digits = 0))
dotplot_obj <- dotplot(ego, showCategory = 10, label_format = 200)+
  theme(axis.text.y = element_text(size = 23))+theme(legend.position = "top",legend.justification = c(1, 1),legend.key.width = unit(15, "mm"))+
  theme(plot.margin = unit(c(0, 0.5, 0, 0.5), "cm"))+
  theme(axis.text.x = element_text(size = 20))+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 42)) +
  theme(axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))+
  theme(legend.title = element_text(vjust = 1),
        legend.text = element_text(vjust = 0.5))
dotplot_obj
#+theme(legend.key.height = unit(50, "pt"))
ggsave("go_myeloid_up.pdf", plot = dotplot_obj, device = "pdf", width = 9, height = 7.6)
ggsave("go_myeloid_up.pdf", plot = dotplot_obj, device = "pdf", width = 9.1, height = 5.5)
write.table(ego_res,"GO_Myeloid-1*10-10_up.txt",sep = "\t")
save(ego,ego_res,file = "GO_Myeloid-1*10-10_up.Rdata")
#Myeloid_down genes#
DEG <- deg2
genelist <- bitr(DEG$gene_sample, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by=c("gene_sample"="SYMBOL"))
#GO
ego <- enrichGO(gene = DEG$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)
ego_res <- ego@result
saveRDS(ego,"go_ego_myeloid_down.rds")
#m_down####
ego <- readRDS("go_ego_myeloid_down.rds")
dotplot_obj <- dotplot(ego, showCategory = 10, label_format = 200)+
  theme(axis.text.y = element_text(size = 23))+theme(legend.position = "top",legend.justification = c(1, 1),legend.key.width = unit(10, "mm"))+
  theme(plot.margin = unit(c(0, 0.5, 0, 0.5), "cm"))+
  theme(axis.text.x = element_text(size = 20))+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 42)) +
  theme(axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))+
  theme(legend.title = element_text(vjust = 1),
        legend.text = element_text(vjust = 0.5))
dotplot_obj
#+theme(legend.key.height = unit(50, "pt"))
ggsave("go_myeloid_down.pdf", plot = dotplot_obj, device = "pdf", width = 9, height = 7.6)
ggsave("go_myeloid_down.pdf", plot = dotplot_obj, device = "pdf", width = 9.1, height = 5.5)
write.table(ego_res,"GO_Myeloid-1*10-10_down.txt",sep = "\t")
save(ego,ego_res,file = "GO_Myeloid-1*10-10_down.Rdata")
#5.8 9.6
