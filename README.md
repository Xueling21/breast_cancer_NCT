# breast_cancer_NCT
Input:
1, exp26HE_178.rds, Her2-positive patients' gene expression profiles in the test set, GSE20271, MAS5 data as downloaded and then normalize by "normalizeBetweenArrays" funcation of "limma" R package.

2, exp63HE.rds, Her2-positive patients' gene expression profiles in the training set, GSE163882, TPM data normalize by "normalizeBetweenArrays" funcation of "limma" R package.

3, exp63TN_178.rds, ER-negative and Her2-negative patients' gene expression profiles in the test set, GSE20271, MAS5 data as downloaded and then normalize by "normalizeBetweenArrays" funcation of "limma" R package.

4, exp68ER.rds, ER+Her2- patients' gene expression profiles in the training set, GSE163882, TPM data normalize by "normalizeBetweenArrays" funcation of "limma" R package.

5, exp89ER_178.rds, ER+Her2- patients' gene expression profiles in the test set, GSE20271, MAS5 data as downloaded and then normalize by "normalizeBetweenArrays" funcation of "limma" R package.

6, exp90TN.rds, ER-Her2- patients' gene expression profiles in the training set, GSE163882, TPM data normalize by "normalizeBetweenArrays" funcation of "limma" R package.

7, exp178.rds. all gene expression profiles in the test set, GSE20271, MAS5 data as downloaded and then normalize by "normalizeBetweenArrays" funcation of "limma" R package.

8, GDSC2_Res.rds, Drug IC50 Reference Document input in "oncoPredict", where each row is a cell line and each column is a drug.

9, pdata178.txt and pdata221.txt are information of clinical phenotypes in test and training dataset, respectively, obtained by pData funciton of the GEOquery package. 





Output:
1,BRCA30_nomogram.RDS,    output by logistic regression model,    f_lrm <-lrm(Group~KIAA1467+CCND1+CCNA2+ESR1+CD52
            , data=Mydata)

2,brca221_5_durg.wilcox.txt,   output by "calcPhenotype" funcation of"oncoPredict" R package,   Drug,   Group1,high-risk;	mean1and median1, mean and median of IC50 of high-risk group in the training data;    Group2, low-risk; mean2 and median2, mean and median of IC50 of low-risk group in the training data; p_value, p_adjust, fold change of High-risk group over low-risk group mean

3,brca221_TIDE_resultslog2.csv, output by TIDE,  http://tide.dfci.harvard.edu/;

4,brca1294-GO.txt, DEGs between pCR and RD group from the training data set by logFC=0.3, P.Value = 0.05, output by "enrichGO" funcation of "clusterProfiler" R package

5,brca1294-kegg.txt, DEGs between pCR and RD group from the training data set by logFC=0.3, P.Value = 0.05, output by "enrichKEGG" funcation of "clusterProfiler" R package

6,deg_sur_neg.txt, The intersect of up genes in RD breat cancer and OS poor prognosis-related genes, and the corresponding logFC and p-values between pCR and RD, output by "venn.diagram" funcation of "VennDiagram" R package

7,deg_sur_pos.txt, The intersect of down genes in RD breat cancer and OS better prognosis-related genes, and the corresponding logFC and p-values between pCR and RD, output by "venn.diagram" funcation of "VennDiagram" R package

8, deg300-GSE163882.txt, the top 300 differentially expressed genes between pCR and RD groups according to logFC with p-value less than 0.05. output by "lmFit" and "eBayes" funcation of "limma" R package

9, deg-GSE163882_all.txt,  differential analysis of the whole genome genes between pCR and RD groups in the training data, output by "lmFit" and "eBayes" funcation of "limma" R package

10, ESTIMATE_score.gct, ESTIMATE_score and ImmuneScore in each patient in the training set output by "outputGCT" funcation of "estimate" R package

11, gene_df_c.rds, genome-wide overall survival analysis results, in which each row is a gene, HR is hazard ratio and its p-value, with the median of each gene expression as the threshold of high- and low- groups. 

12, gene_df20, the 20 candidate genes with consistent prognosis, i.e., up in RD and poor prognosis or down in RD and better prognosis, and there p-values and HR resulted from overall survival analysis.

13, GO_GSE163882_deg-0.3-0.05.Rdata, Gene Ontology Enrichment results from DEGs in the training dataset between pCR and RD groups.

14, KEGG_GSE163882_DEG-0.3-0.05.Rdata, KEGG Enrichment results from DEGs in the training dataset between pCR and RD groups.

15, lassoGenes.rds, 17 genes selected by Lasso selection.

16, meta_data178.rds, meta_data221.rds, meta_data221_predict.rds, meta_data26_HERP_178.rds, meta_data63_ERN_HERN_178.rds, meta_data63HERP.rds, meta_data68ERP_HER2N.rds, meta_data89_ERP_HERN_178.rds, meta_data90TN.rds are subtype clinical information of  PCR and RD in test dataset, pCR and RD information in training dataset, pCR and RD and asssigned high- or low-group by the prediction model in training dataset, same information in ER-HER+ test dataset, ER-HER- test dataset, HER+ training dataset, ER+HER- training dataset, ER+HER- test dataset, ER-HER- training dataset, respectively.

17, rfGenes1.rds, svmGenes.rds, and vnn5.rds are genes selected from random forest, support vector machine and the overlapping genes selected by the three machine learning classifiers, random forest, SVM and lasso, respectively.


18, brca26_cellsgene.txt, 10-cell-type gene expression signature matrix of whole genome, where each row is a gene, and each column is a celltype

19,  CIBERSORT-Results,rds and CIBERSORTx_Job1_Results.txt are deconvoluted results from training dataset with home-made 10-cell-type gene expression signature matrix and LM22, respectively.

20, deg-Bcells.txt, deg-Myeloid.txt, deg-Tcells.txt, are DEGs from the CD52 high- and low- groups of T, B, and myeloid cells, respectively.

21, GO_Bcell-1e-10_down.txt, GO_Bcell-1e-10_up.txt, GO_Myeloid-1e-10_down.txt, GO_Myeloid-1e-10_up.txt, GO_Tcell-1e-10_down.txt, and GO_Tcell-1e-10_up.txt are Gene Ontology enrichment analysis results of the DEGs from the CD52 high- and low- groups of T, B, and myeloid cells, respectively.
       
22, pbmc1_metadata.csv and top10_brca26.csv are cell type annotated results of the scRNA-Seq data and the top 10 marker genes for each cluster, respectively.  


Codes:
Bulk_analysis.R, codes for bulk data analysis for training dataset, GSE163882, and test datset, GSE202712.
scRNA_analysis.R, codes for scRNA-Seq data analysis, GSE176078.





