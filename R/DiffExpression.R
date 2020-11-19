#DESeq2 performs for each gene a hypothesis test to see whether evidence is sufficient to decide against the null hypothesis that there is no effect of the treatment on the gene and that the observed difference between treatment and control was merely caused by experimental variability (i. e., the type of variability that you can just as well expect between different samples in the same treatment group)

#results(dds, alpha=alpha)

#--------------------------------------------------------------------------------------------------------
#Load Libraries
install.packages(c("yaml","gridExtra","checkmate","RSQLite","gplots","RColorBrewer","ashr"), dependencies = TRUE)
library(gridExtra)
library(yaml)
library(checkmate)
library(RSQLite)
library(gplots)
library(RColorBrewer)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DESeq2","data.table","vctrs","VennDiagram","pcaExplorer","org.Hs.eg.db","topGO","BiocParallel","pheatmap","genefilter","limma","edgeR","topconfects","NBPSeq","clusterProfiler","enrichplot","ggupset","pathview","EnhancedVolcano", "STRINGdb"))

library(topGO)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(vctrs)
library(DESeq2)
library(data.table)
library(pcaExplorer)
library(BiocParallel)
library(pheatmap)
library(genefilter)
library(edgeR)
library(topconfects)
library(limma)
library(NBPSeq)
library(matrixStats)
library(clusterProfiler)
library(VennDiagram)
library(enrichplot)
library(ggupset)
library(pathview)
library(EnhancedVolcano)
library(STRINGdb)
library(dplyr)

register(MulticoreParam(4))

#-------------------------------------------------------------------------------------------------------------------------------------------
#Sample Info
sampleinfo=read.delim('Sampletable_noM57RV.txt',header=TRUE)
sampleinfo$tissuebmi<-paste0(sampleinfo$tissue,sampleinfo$bmi)
coldata=sampleinfo
countdata = read.csv("raw_counts_noM57RV.csv",row.names=1,check.names=TRUE, header = TRUE)

#geneannotation
#columns(org.Hs.eg.db)
anno_df <- data.frame(gene_id = rownames(countdata),
                      stringsAsFactors=FALSE)

anno_df$gene_name<- mapIds(org.Hs.eg.db,
                           keys=anno_df$gene_id,
                           column="SYMBOL",
                           keytype="ENSEMBL",
                           multiVals="first")

anno_df$genedescript <- mapIds(org.Hs.eg.db,
                               keys=anno_df$gene_id,
                               column="GENENAME",
                               keytype="ENSEMBL",
                               multiVals="first")
anno_df$ont <- mapIds(org.Hs.eg.db,
                      keys=anno_df$gene_id,
                      column="ONTOLOGYALL",
                      keytype="ENSEMBL",
                      multiVals="first")

anno_df$transcipt <- mapIds(org.Hs.eg.db,
                            keys=anno_df$gene_id,
                            column="ENSEMBLTRANS",
                            keytype="ENSEMBL",
                            multiVals="first")

anno_df$Entrez <- mapIds(org.Hs.eg.db,
                         keys=anno_df$gene_id,
                         column="ENTREZID",
                         keytype="ENSEMBL",
                         multiVals="first")

anno_df$KEGG<-mapIds(org.Hs.eg.db,
                     keys=anno_df$gene_id,
                     column="PATH",
                     keytype="ENSEMBL",
                     multiVals="first")

anno_df$Uni<-mapIds(org.Hs.eg.db,
                    keys=anno_df$gene_id,
                    column="UNIPROT",
                    keytype="ENSEMBL",
                    multiVals="first")

anno_df$Uni<-mapIds(org.Hs.eg.db,
                    keys=anno_df$gene_id,
                    column="UNIPROT",
                    keytype="ENSEMBL",
                    multiVals="first")


rownames(anno_df) <- anno_df$gene_id
#-------------------------------------------------------------------------------------------------------------------------------------------
#DeSeq Read in Sample Info
#run DeSeq Pipeline
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData = coldata,
  rowData=rownames(countdata),
  design = ~ tissuebmi)

#run DeSeq Pipeline
dds<-DESeq(ddsFullCountTable, parallel=TRUE)


bg_ids <- rownames(dds)[rowSums(counts(dds)) > 0]
pca2go <- pca2go(dds,
                 annotation = anno_df,
                 organism = "Hs",
                 ensToGeneSymbol = TRUE,
                 background_genes = bg_ids)

write.table(pca2go, file = "GO AnnotationAdipose.txt", row.names = T, sep = "\t")


#Confects
dconfects <- deseq2_confects(dds, contrast=c("tissuebmi","AdipNorm","AdipOb"), fdr=.05,step=0.05, cooksCutoff=Inf,independentFiltering=FALSE)

dconfects <- deseq2_confects(dds, contrast=c("tissuebmi","MyoNorm","MyoOW"), fdr=.1,step=0.05)

shrunk <- lfcShrink(dds, contrast=c("tissuebmi","AdipNorm","AdipOb"), type="normal")
dconfects$table$shrunk <- shrunk$log2FoldChange[dconfects$table$index]
t<-dconfects$table
nu<-dconfects$table[,3==0]

#Here the usual logFC values estimated by limma are shown as dots, with lines to the confect value.
confects_plot(dconfects) + 
  geom_point(aes(x=shrunk, size=baseMean, color="lfcShrink"), alpha=0.75)


confects_plot(dconfects)

confects_plot_me(confects)


#Create All results
dir.create("DESeq2_Final", showWarnings = FALSE)
contras=levels(dds$tissuebmi)
contrcmb<-combn(contras, 2, simplify = TRUE)
cat(contras,"\t",length(contras),"elements","\n",ncol(contrcmb),"possible pairwise comparisons","\n",print(contrcmb),file="DESeq2_readcontrast.txt")

for(i in seq(1, ncol(contrcmb))){
  res <- results(dds,contrast=c("tissuebmi",as.character(contrcmb[1,i]),as.character(contrcmb[2,i])))
  res <- res[order(res$pvalue),]
  res1 <- as.data.frame(res)
  res1<-cbind(res1, anno_df$gene_name[match(rownames(res1), rownames(anno_df))],
              anno_df$genedescript[match(rownames(res1), rownames(anno_df))],
              anno_df$ont[match(rownames(res1), rownames(anno_df))],
              anno_df$transcipt[match(rownames(res1), rownames(anno_df))],
              anno_df$Entrez[match(rownames(res1), rownames(anno_df))],
              anno_df$KEGG[match(rownames(res1), rownames(anno_df))],
              anno_df$Uni[match(rownames(res1), rownames(anno_df))]
              )
  
  colnames(res1)[7]<-"genename"
  colnames(res1)[8]<-"GeneDescription"
  colnames(res1)[9]<-"GOAnnoation"
  colnames(res1)[10]<-"TrasncriptID"
  colnames(res1)[11]<-"EntrezID"
  colnames(res1)[12]<-"KEGG"
  colnames(res1)[13]<-"Uniprot"
  
  
  res1$regulation<-factor(ifelse(res1$log2FoldChange > 0, "up", "down"))
  write.table(res1,file=paste("DESeq2_Final/DE_",contrcmb[1,i],"VS",contrcmb[2,i],".txt",sep=""),sep="\t",col.names=TRUE, row.names = T)
  
  #save only diffExpresed where FDR is less than .1 and absolute value of the fold change is greater than 1
  selectedET <- res1$pvalue < 0.05 & abs(res1$log2FoldChange) > 3
  selectedET <- res1[selectedET, ]
  selectedET<-selectedET[complete.cases(selectedET),]
  write.table(selectedET,file=paste("DESeq2_Final/DiffExpress_",contrcmb[1,i],"VS",contrcmb[2,i],".txt",sep=""),sep="\t",col.names=TRUE, row.names = T)
  
  #png(paste("DESeq2_Final/DE_",contrcmb[1,i],"VS",contrcmb[2,i],"_Volcano Plot.png", sep=""))
  pic<-ggplot(data=res1, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
    geom_point() +
    theme_minimal() +
    geom_text_repel() +
    scale_color_manual(values=c("blue", "black", "red")) +
    geom_vline(xintercept=c(-0.6, 0.6), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red")+
    ggtitle("Differntial Gene Expression Volcano Plot")
    
     EnhancedVolcano(res1,
                   lab=paste(res1$genename,
                   rownames(res1),sep=","),
                   lab=res1$genename,
                   x="log2FoldChange",
                   y="pvalue",
                   xlim=c(-10,10),
                   ylim=c(0,8),
                   pCutoff = .05,
                   FCcutoff = 3,
                   hline = c(10e-2,10e-3,10e-4),
                   hlineCol = c('grey0', 'grey25','grey50'),
                   pointSize=c(ifelse((res1$log2FoldChange>1||res1$log2FoldChange< -1), 4,.5)),
                   hlineType = "dash",
                   hlineWidth = .8,
                   gridlines.major = FALSE,
                   gridlines.minor=FALSE,
                   labSize = 3,
                   title="Differntial Gene Expression Volcano Plot",
                   caption="FC cuttoff =|3|, alpha=.05",
                   legendPosition = "right",
                   subtitle=paste(contrcmb[1,i],"VS",contrcmb[2,i]))
  print(pic)
  dev.off()

  png(paste("DESeq2_Final/DE_",contrcmb[1,i],"VS",contrcmb[2,i],"_MA Plot.png", sep=""))
   picb<-DESeq2::plotMA(res,
                        main=c(paste(contrcmb[1,i]," VS ",contrcmb[2,i]," MA Plot",sep="")),
                        ylim=c(-2,2))
   print(picb)
   dev.off()
  
  
}

#baseMean: the average of the normalized count values, dividing by size factors, taken over all samples. 
#log2FoldChange: is the effect size estimate. It tells us how much the gene's expression seems to have changed due to treatment with DPN in comparison to control. This value is reported on a logarithmic scale to base 2: for example, a log2 fold change of 1.5 means that the gene's expression is increased by a multiplicative factor of 2^1.5 ??? 2.82.
#lfcSE: the standard error estimate for the log2 fold change estimate.
#p value indicates the probability that a fold change as strong as the observed one, or even stronger, would be seen under the situation described by the null hypothesis.)
#padj: The DESeq2 software automatically performs independent filtering which maximizes the number of genes which will have adjusted p value less than a critical value (by default, alpha is set to 0.1)
#---------------------------------------------------------------------------------------------------------------------------------------------
#EdgeR


Edge.List<- DGEList(counts=countdata, group=factor(sampleinfo$tissuebmi))
e.list<-Edge.List

#normalize samples by keeping samples with row counts above 100 for at least 11 samples
apply(e.list$counts,2,sum)
keep<-rowSums(cpm(e.list)>100)<=11
d<- e.list[keep,]
d$samples$lib.size<-colSums(d$counts)
d<-calcNormFactors(d)

#Estimate Dispersions
d1<-estimateCommonDisp(d, verbose = T)
d1<-estimateTagwiseDisp(d1)

design.mat<-model.matrix(~0 + d$samples$group)
colnames(design.mat)<-levels(d$samples$group)
d2<-estimateGLMCommonDisp(d,design.mat)
d2<- estimateGLMTrendedDisp(d2, design.mat, method="power")
d2<- estimateGLMTagwiseDisp(d2, design.mat)

#DataVisualization
plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("topright", as.character(unique(d$samples$group)), col=1:6, cex=.7, pch=1)

#plots tagwise biological coefficient of variation (square root of dispersions) against log2-CPM
plotBCV(d1)

#plot GLM Tagwaise dispersion
plotBCV(d2)

#DEG EdgeR
dir.create("EdgeR_Final", showWarnings = FALSE)
contras=unique(sampleinfo$tissuebmi)
contrcmb<-combn(contras, 2, simplify = TRUE)
cat(contras,"\t",length(contras),"elements","\n",ncol(contrcmb),"possible pairwise comparisons","\n",print(contrcmb),file="EdgeR_readcontrast.txt")

for(i in seq(1, ncol(contrcmb))){
  res <- exactTest(d1,pair=c(as.character(contrcmb[1,i]),as.character(contrcmb[2,i])))
  resa <- res[order(res$table$PValue),]
  res1 <- as.data.frame(resa)
  res1<-cbind(res1, anno_df$gene_name[match(rownames(res1), rownames(anno_df))],
              anno_df$genedescript[match(rownames(res1), rownames(anno_df))],
              anno_df$ont[match(rownames(res1), rownames(anno_df))],
              anno_df$transcipt[match(rownames(res1), rownames(anno_df))],
              anno_df$Entrez[match(rownames(res1), rownames(anno_df))],
              anno_df$KEGG[match(rownames(res1), rownames(anno_df))],
              anno_df$Uni[match(rownames(res1), rownames(anno_df))]
  )
  
  colnames(res1)[4:10]<-c("genename","GeneDescription","GOAnnoation","TrasncriptID","EntrezID","KEGG","Uniprot")
  
  
  res1$regulation<-factor(ifelse(res1$logFC > 0, "up", "down"))
  write.table(res1,file=paste("EdgeR_Final/EdgeRDE_",contrcmb[1,i],"VS",contrcmb[2,i],".txt",sep=""),sep="\t",col.names=TRUE, row.names = T)
  
  #save only diffExpresed where FDR is less than .1 and absolute value of the fold change is greater than 1
  selectedET <- res1$Pvalue < 0.05 & abs(res1$logFC) > 3
  selectedET <- res1[selectedET, ]
  selectedET<-selectedET[complete.cases(selectedET),]
  write.table(selectedET,file=paste("EdgeR_Final/EdgeRDiffExpress_",contrcmb[1,i],"VS",contrcmb[2,i],".txt",sep=""),sep="\t",col.names=TRUE, row.names = T)
  
   png(paste("EdgeR_Final/DE_",contrcmb[1,i],"VS",contrcmb[2,i],"_Volcano Plot.png", sep=""))
   pic<-EnhancedVolcano(res1,
                        #lab=paste(res1$genename, rownames(res1),sep=","),
                        lab=res1$genename,
                        x="logFC",
                        y="PValue",
                        xlim=c(-10,10),
                        ylim=c(0,10),
                        pCutoff = .05,
                        FCcutoff = 3,
                        #hline = c(10e-2,10e-3,10e-4),
                        hlineCol = c('grey0', 'grey25','grey50'),
                        #pointSize=c(ifelse((res1$logFC>1||res1$log2FC< -1), 20,10)),
                        hlineType = "dash",
                        hlineWidth = .8,
                        gridlines.major = FALSE,
                        gridlines.minor=FALSE,
                        labSize = 3,
                        title="EdgeR: Differntial Gene Expression Volcano Plot",
                        caption=paste(contrcmb[1,i],"VS",contrcmb[2,i]),
                        legendPosition = "right",
                        subtitle="FC cuttoff =|3|, alpha=.05")
   print(pic)
   dev.off()
}
#--------------------------------------------------------------------------------------------------------------------
#Limma Voom




#---------------------------------------------------------------------------------------------------------------------------------------------
#read in Deseqs

  fn<-list.files("DESeq2_Final/")
  fn.names<-gsub(".txt","",fn)
  fn.names.t<-gsub("DE_","",fn.names)
  fn.top<-c("AdipNormVSAdipOb","AdipNormVSAdipOW","AdipObVSAdipOW","MyoNormVSMyoOb","MyoNormVSMyoOW","MyoObVSMyoOW")
  for(i in seq(1, length(fn.names))){
    assign(fn.names.t[i], read.table(paste("DESeq2_presentation/",fn[i], sep=""), row.names = 1, header= TRUE)) 
  }
  
  fn<-list.files("EdgeR/")
  fn.names<-gsub(c(".txt"),"",fn)
  #fn.names.t<-gsub("DE_","",fn.names)
  #fn.top<-c("AdipNormVSAdipOb","AdipNormVSAdipOW","AdipObVSAdipOW","MyoNormVSMyoOb","MyoNormVSMyoOW","MyoObVSMyoOW")
  for(i in seq(1, length(fn.names))){
    assign(fn.names[i], read.table(paste("EdgeR/",fn[i], sep=""), row.names = 1, header= TRUE)) 
  }
  
  #find matching samples
  a<-EdgeRDE_AdipNormVSAdipOb
  b<-AdipNormVSAdipOb
  n<-nrow(EdgeRDE_AdipNormVSAdipOb)
  m<-nrow(AdipNormVSAdipOb)
  q<-0
  small<-0
  big<-0
  if(n>m){
  q<-n
  big<-a
  small<-b
  }else{
    q<-m
    big<-b
    small<-a
  }
  A.NormvOB<-merge(big,small, by="row.names",all.y=TRUE)
  write.table(selectedET,file=paste("CombineDiffExpress_AdiposeNormVOB",".txt",sep=""),sep="\t",col.names=TRUE, row.names = T)

  
-----------------------------------------------------------------------------------------------------  
  venn.diagram(list(rownames(MyoObVSMyoOW),rownames(EdgeRDE_MyoOWVSMyoOb)),
               category.names = c("DeSeq","EdgeR"),
               fill=c("red","green"),
               alpha=c(.05,.05),
               cex=2,
               cat.fontface=4,
               filename = "test.png",
               output=TRUE)
  
  #VennDiagrams of overlapp

   png("venn_diagram_NormaltoObese.png")
  venn.plot.normob <- venn.diagram(list(rownames(DiffExpress_AdipNormVSAdipOb), rownames(DiffExpress_MyoNormVSMyoOb)), NULL, fill=c("red", "green"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("Adipose", "Myocardium"), 
                                   main="Normal to Obese Transition",
                                   sub="p-value < .05 & Fold Change > |3|",
                                   main.cex = 2,
                                   sub.cex = 1,
                                   sub.fontface = 4)
  grid.draw(venn.plot.normob)
  dev.off()
  
  png("venn_diagram_NormaltoOverweight.png")
  venn.plot.normow <- venn.diagram(list(rownames(DiffExpress_AdipNormVSAdipOW), rownames(DiffExpress_MyoNormVSMyoOW)), NULL, fill=c("red", "green"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("Adipose", "Myocardium"), main="Normal to Overweight Transition",
                                   sub="p-value < .05 & Fold Change > |3|",
                                   main.cex = 2,
                                   sub.cex = 1,
                                   sub.fontface = 4)
  grid.draw(venn.plot.normow)
  dev.off()
  
  png("venn_diagram_ObesetoOverweight.png")
  venn.plot.owob <- venn.diagram(list(rownames(DiffExpress_AdipObVSAdipOW),rownames(DiffExpress_MyoObVSMyoOW)), NULL, fill=c("red", "green"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("Adipose", "Myocardium"), main="Overweight to Obese Transition",
                                 sub="p-value < .05 & Fold Change > |3|",
                                 main.cex = 2,
                                 sub.cex = 1,
                                 sub.fontface = 4)
  grid.draw(venn.plot.owob)
  dev.off()
  
  
  venn.plot.myo<- venn.diagram(list
                                 (rownames(DiffExpress_MyoObVSMyoOW), 
                                   rownames(DiffExpress_MyoNormVSMyoOb),
                                   rownames(DiffExpress_MyoNormVSMyoOW)), 
                                 NULL, fill=c("red", "green","blue"), 
                                 alpha=c(0.5,0.5,0.5), cex = 2, cat.fontface=4, 
                                 category.names=c("Obese to Overweight",
                                                  "Normal to Obese",
                                                  "Normal to Overweight"
                                                  ), 
                                 main="Myocardium BMI Transitions",
                               sub="p-value < .05 & Fold Change > |3|",
                               main.cex = 2,
                               sub.cex = 1,
                               sub.fontface = 4)
  png("venn_diagram_Myocardium.png")
  grid.draw(venn.plot.myo)
  dev.off()
  
  venn.plot.adip<- venn.diagram(list
                               (rownames(DiffExpress_AdipObVSAdipOW), 
                                 rownames(DiffExpress_AdipNormVSAdipOb),
                                 rownames(DiffExpress_AdipNormVSAdipOW)), 
                               NULL, fill=c("red", "green","blue"), 
                               alpha=c(0.5,0.5,0.5), cex = 2, cat.fontface=4, 
                               category.names=c("Obese to Overweight",
                                 "Normal to Obese",
                                 "Normal to Overweight"
                               ), 
                               main="Adipose BMI Transitions",
                               sub="p-value < .05 & Fold Change > |3|",
                               main.cex = 2,
                               sub.cex = 1,
                               sub.fontface = 4)
  png("venn_diagram_Adipose.png")
  grid.draw(venn.plot.adip)
  dev.off()
  
  
#GGO: Identify purpose test Adipose normal V Obese
  OrgDb <- org.Hs.eg.db
n<-1
  
  for(i in list(AdipNormVSAdipOb,AdipNormVSAdipOW,AdipObVSAdipOW,MyoNormVSMyoOb,MyoNormVSMyoOW,MyoObVSMyoOW)){
    
    dir.create(fn.top[n])

    print(getwd())
    
    #go
    geneList <- as.vector(i$log2FoldChange)
    print("gene list entered")
    names(geneList) <- as.character(i$EntrezID)
    geneList<- sort(geneList, decreasing = TRUE)
    print("names entered")
    gene <- as.character(na.omit(i$EntrezID))
    print("gene  entered")
    
    #Kegg
    geneListk <- as.vector(i$log2FoldChange)
    print("gene list entered")
    names(geneListk) <- as.character(i$KEGG)
    geneListk<- sort(geneListk, decreasing = TRUE)
    print("names entered")
    genek <- as.character(na.omit(i$KEGG))
    print("gene  entered")
    
    
    
    #Biological Process
    ggo <- clusterProfiler::groupGO(gene     = gene,
                                    OrgDb    = OrgDb,
                                    ont      = "BP",
                                    level    = 3,
                                    readable = TRUE)
    #head(summary(ggo)[,-5])
    png(paste0(fn.top[n],"_BiologicalProcess.png"), width=1000, height=1000)
    print(barplot(ggo, drop=TRUE, showCategory=25))
    dev.off()
    
    #Cellular Component
    ggo <- clusterProfiler::groupGO(gene     = gene,
                                    OrgDb    = OrgDb,
                                    ont      = "CC",
                                    level    = 3,
                                    readable = TRUE)
    #head(summary(ggo)[,-5])
    png(paste0(fn.top[n],"_CellularComponent.png"), width=1000, height=1000)
  print(barplot(ggo, drop=TRUE, showCategory=25))
    dev.off()
    
    #Molecular Function
    ggo <- clusterProfiler::groupGO(gene= gene,
                                    OrgDb    = OrgDb,
                                    ont      = "MF",
                                    level    = 3,
                                    readable = TRUE)
    #head(summary(ggo)[,-5])
    png(paste0(fn.top[n],"_MolecularFunction.png"), width=1000, height=1000)
    print(barplot(ggo, drop=TRUE, showCategory=25))
    dev.off()

    #Gene Disease Enrichment
    b<-enrichDGN(gene=genek,
                 pvalueCutoff = .05,
                 minGSSize = 5,
                 readable = TRUE)
    edox <- setReadable(b, 'org.Hs.eg.db', 'ENTREZID')
    # p1 <- heatplot(edox)
    p2 <- heatplot(edox, foldChange=geneList)
    
    png(paste0(fn.top[n],"_Enrichment HeatMap Functional Classification.png"), width=1000, height=1000)
    
    print(cowplot::plot_grid(p2, ncol=1, labels=paste(fn.top[n],"Enrichment")))
    dev.off()
    
    # p1 <- emapplot(b)
    # p2 <- emapplot(b, pie_scale=1.5)
    # p3 <- emapplot(b,layout="kk")
    p4 <- emapplot(b, pie_scale=1.5,layout="kk") 
    png(paste0(fn.top[n],"_Enrichment NetworkMap.png"), width=1000, height=1000)
    print(cowplot::plot_grid(p4, ncol=1, labels=paste0(fn.top[n],"_Enrichment HeatMap Network Map Classification.png")))
    dev.off()
    
    

      #Overrepresentation test: representation is just whether it is there or not. over-representation is the same as enrichment which means that, say you are working with genes, several genes from the same pathway are highly expressed or even several genes from the same pathway are diffenrentially expressed comapring two samples. In contrast, representation would than mean that genes are expressed but in no certain / specific manner.
    
    #used total as background universe
 
  #Biological Procces 
  over.go <- clusterProfiler::enrichGO(gene= gene,
                                   OrgDb= OrgDb,
                                   ont = "BP",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.05,
                                   qvalueCutoff  = 0.05, 
                                   readable      = TRUE)
  png(paste0(fn.top[n],"_BiologicalProcces_Overrepresentation.png"), width=1000, height=1000)
  print(clusterProfiler::dotplot(over.go, showCategory=25))
  dev.off()
  
 
  
  #Cellular Component 
  over.go <- clusterProfiler::enrichGO(gene= gene,
                                       OrgDb= OrgDb,
                                       ont = "CC",
                                       pAdjustMethod = "BH",
                                       pvalueCutoff  = 0.05,
                                       qvalueCutoff  = 0.05, 
                                       readable      = TRUE)
  png(paste0(fn.top[n],"_CellularComponent_Overrepresentation.png"), width=1000, height=1000)
  print(clusterProfiler::dotplot(over.go, showCategory=25))
  dev.off()
  
  #MolecularFunction 
  over.go <- clusterProfiler::enrichGO(gene= gene,
                                       OrgDb= OrgDb,
                                       ont = "MF",
                                       pAdjustMethod = "BH",
                                       pvalueCutoff  = 0.05,
                                       qvalueCutoff  = 0.05, 
                                       readable      = TRUE)
  png(paste0(fn.top[n],"_MolecularFunction_Overrepresentation.png"), width=1000, height=1000)
  print(clusterProfiler::dotplot(over.go, showCategory=25))
  dev.off()
  
  
  #Kegg
  hs<-search_kegg_organism('hsa', by='kegg_code')
  kk <- enrichKEGG(gene = genek,
                   organism = 'hsa',
                   pvalueCutoff = 0.05)
  
  write.table(kk@result,file=paste(fn[n],"_Kegg Overrepresentation",".txt",sep=""),sep="\t",col.names=TRUE, row.names = T)
  
  png(paste0(fn.top[n],"_EnrichedConcepts Barplot.png"), width=1000, height=1000)
  print(barplot(kk,
                Title="Enriched Kegg Pathways",
                font.size = 50,
                showCategory = 25))
  dev.off()
  
  png(paste0(fn.top[n],"_EnrichedConcepts Network.png"), width=1000, height=1000)
  print(cnetplot(kk, categorySize="pvalue", foldChange=geneListk))
  dev.off()
  
  for(i in 1:3){
  dme <- pathview(gene.data=geneListk, pathway.id=kk@result$ID[i], species = hs)}
  
  
  
  kk2 <- gseKEGG(geneList     = geneListk,
                 organism     = 'hsa',
                 nPerm        = 1000,
                 minGSSize    = 10,
                 pvalueCutoff = .1,
                 verbose      = FALSE)
  
  write.table(kk2@result,file=paste(fn[n],"_Kegg Gene Set Enrichment",".txt",sep=""),sep="\t",col.names=TRUE, row.names = T)
  

  
  png(paste0(fn.top[n],"_EnrichedConcepts.png"), width=1000, height=1000)
  print(upsetplot(kk2))
  dev.off()
  
  
  n<-n+1
  
  print(paste("WD Reset n is ",n))
  print(fn.top[n])
  }
  
  
#Sequencing Depth Per Library
depth<-colSums(counts(dds))
barplot(depth, las=2, main="Counts Bar Plot")


#histogram of pvalue
par(mfrow=c(1,2))
#png("Histogram of pValues and FDR.png")
hist(AdipObVSAdipOW$pvalue, main="Obese vs Normal PValue", xlim = c(0,.06),freq=FALSE)


hist(AdipObVSMyoNorm$padj, main="Obese vs Normal FDR")

#dev.off()

par(mfrow=c(1,1))

#count which pass:
fdr <- as.numeric('0.1')
length(which(AdipObVSMyoNorm$padj < fdr)) #12082 pass

alpha <- as.numeric('0.05')
#ob.alpha<-summary(na.omit(res.obese.df$pvalue))[2]
length(which(AdipNormVSAdipOb$pvalue < .05)) #1271 pass


#MAPlot, mean of normalized counts
#MA-plot shows the log2 fold changes from the condition over the mean of normalized counts, i.e. the average of counts normalized by size factor.Items in red are genes with adjusted pvalues below the threshhold of .1. Points which fall out of the window are plotted as open triangles pointing either up or down.

#dispersion plot
par(mfrow=c(1,1))
png("Dispersion Plot.png", width=1000, height=1000)
plotDispEsts( dds, ylim = c(1e-6, 1e2),
              main="Gene Dispersal Plot")
dev.off()


#ratio of pvals

#PCA

rld <- rlog(dds)
#head( assay(rld))
#Heatmap plot
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists )

rownames(sampleDistMatrix) <- paste(rld$tissuebmi,
                                    rld$sampleName,
                                    sep="-")
colnames(sampleDistMatrix) <- paste(rld$tissuebmi,
                                    rld$sampleName,
                                    sep="-")
sampleDistMatrix <-sampleDistMatrix[,sort(colnames(sampleDistMatrix))]
colours = colorRampPalette(rev(brewer.pal(9, "Reds")))(255)
png("HeatMap_dropM57RVM.png", width=1000, height=1000)
heatmap.2(sampleDistMatrix, 
          trace="none", 
          Rowv=T,
          Colv = T,
          #dendrogram="row",
          col=colours,
          cexRow = .7,
          cexCol = .7)
dev.off()

#PCA
par(mfrow=c(1,1))
pc<-plotPCA(rld, intgroup="tissuebmi", returnData=TRUE)
pc<-cbind(pc, sampleinfo$tissue, sampleinfo$bmi)

# bc<-plotPCA(rld, intgroup="eatcondition", returnData=TRUE)

png("PCA by BMI Symbols.png", width=1000, height=1000)

plot(pc[,1],pc[,2], 
     col=as.numeric(pc$`sampleinfo$bmi`),
     cex=3,
     pch=as.numeric(pc$`sampleinfo$tissue`),
     xlab="PC1",
     ylab="PC2",
     main="All Tissues PCA by BMI")
#text(pc,labels=rownames(pc), cex=1, font=2,col=as.numeric(pc$`sampleinfo$bmi`))
#legend("top", legend=c("Myocardium Normal BMI","Adipose Normal BMI","Myocardium Overweight BMI","Adipose Overweight BMI","Myocardium Obese BMI","Adipose Obese BMI"),col=(1:6), pch=19, title=c("BMI"))

legend("bottom", legend=c("Normal BMI","Overweight BMI","Obese BMI"),col=(c(1,3,2)), pch=19, title=c("Body Mass Index"))
dev.off()

#GeneClustering


topVarGenes <- head(order(rowVars(assay(rld) ), decreasing=TRUE ), 35)

png("Top35 Genes HM.png", width=1000, height=1000)
heatmap.2( assay(rld)[ topVarGenes, ], scale="row",
           trace="none", dendrogram="column",
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
dev.off()


#Pc Loadings
rlog<-rlogTransformation(dds)
pcaobj<-prcomp(t(assay(rlog)))

rest<-correlatePCs(pcaobj,colData(dds))
# plotPCcorrs(rest, pc=1, logp=TRUE)
# plotPCcorrs(rest, pc=2, logp=TRUE)
# plotPCcorrs(rest, pc=3, logp=TRUE)
# plotPCcorrs(rest, pc=4, logp=TRUE)
par(mfrow=c(2,2))
barplot(rest[1,2:5],
        main=paste("Relationships between PC1 and Covariates"),
        ylab = "-log10 of pValue",
        col = "green")
barplot(rest[2,2:5],
        main=paste("Relationships between PC2 and Covariates"),
        ylab = "-log10 of pValue",
        col = "green")
barplot(rest[3,2:5],
        main=paste("Relationships between PC3 and Covariates"),
        ylab = "-log10 of pValue",
        col = "red")
barplot(rest[4,2:5],
        main=paste("Relationships between PC4 and Covariates"),
        ylab = "-log10 of pValue",
        col = "blue")

#StringDB
string_db <- STRINGdb$new( version="10", species=9606,
                           score_threshold=200, input_directory="")

  #Read in data as string object
adipobnorm_mapped <- string_db$map( AdipNormVSAdipOb, "genename", removeUnmappedRows = TRUE)
hits <- adipobnorm_mapped$STRING_id[1:100]
string_db$plot_network(hits)

#---------------------------------------------------
#Gene Set Enrichment dpylr

gene.list=AdipNormVSAdipOb$log2FoldChange
names(gene.list)=AdipNormVSAdipOb$genename
gene.list = sort(gene.list, decreasing = TRUE)
gene.list = gene.list[!duplicated(names(gene.list))]
head(gene.list)


#---------------------------------------
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
