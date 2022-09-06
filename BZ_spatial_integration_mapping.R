library(dplyr)
library(dbplyr)
library(ggplot2)
library(Seurat)
library(RColorBrewer)
library(matrixStats)


#load all visium objects
V39 <- Load10X_Spatial("/Users/kevinking/Downloads/10X_2022_08_01//V39/outs/") #sham1 + NP
V39$orig.ident = "sham1+np"

V40 <- Load10X_Spatial("/Users/kevinking/Downloads/10X_2022_08_01//V40/outs/") #sham2 + NP
V40$orig.ident = "sham2+np"

load(file = "V3.robj")#30 min, I/R, D3
V3$orig.ident <- "I/R"
load(file = "V5.robj") # 1 hr
V5$orig.ident <- "1hr"

load(file = "V6.robj") # 4hr
V6$orig.ident <- "4hr"

V36 <- Load10X_Spatial("/Users/kevinking/Downloads/10X_2022_08_01//V36/outs/") #ISO
V36$orig.ident <- "ISO"

load(file = "V1.robj")#D3
V1$orig.ident <- "D3_1"
V37 <- Load10X_Spatial("/Users/kevinking/Downloads/10X_2022_08_01//V37/outs/") #D3
V37$orig.ident <- "D3_2"
V38 <- Load10X_Spatial("/Users/kevinking/Downloads/10X_2022_08_01//V38/outs/") #D3
V38$orig.ident <- "D3_3"

load(file = "V20_NP_D3.robj") #NP_1
V20_NP_D3$orig.ident

V7 <- Load10X_Spatial("/Users/kevinking/Border_zone_project/V7/") #TAC + NP
V7$orig.ident  <- "TAC+np"

load(file = "V20_NP_D3.robj")
V20_NP_D3$orig.ident <- "NP"

#integration
object_list <- c(V39,V5,V6, V40,V1, V37, V38,V3,V36,V20_NP_D3)
object_list <- lapply(X = object_list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = object_list)
object_list <- lapply(X = object_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = object_list, dims = 1:20, reference = 5)
object_integrated <- IntegrateData(anchorset = anchors, dims = 1:20)
object_integrated <- ScaleData(object_integrated, verbose = FALSE)
object_integrated <- RunPCA(object_integrated, verbose = FALSE)
object_integrated <- RunUMAP(object_integrated, dims =1:20)
object_integrated <- FindNeighbors(object_integrated, reduction = "pca", dims = 1:20)
object_integrated <- FindClusters(object_integrated, verbose = FALSE, resolution = .5)

save(object_integrated, file = "object_integrated_visium.robj")

#reordering clusters
object_integrated <- FindClusters(object_integrated, verbose = FALSE, resolution = .8)
UMAPPlot(object_integrated, label = T)
Idents(object_integrated) <- "integrated_snn_res.0.8"
object_integrated<-RenameIdents(object_integrated, '13'=1,'0' = 1,'5' = 2,'3' = 3, 
                                '9' = 4,'7' = 5, '11' =6, '12'=7,
                                '6' = 8,'2' =9,'10'=10,'4'=11,'8'=12,'14'=12)
UMAPPlot(object_integrated, label = T)
object_integrated <- StashIdent(object_integrated, save.name = "high_resolution")
data.frame(object_integrated$RZ_genes, object_integrated$BZ1_genes, object_integrated$BZ2_genes, object_integrated$high_resolution)

#sn/scRNA-seq derived gene-set scores 

load(file ="/Users/home/hrt.integrated_sc_nuclei_filtered.robj") #load integrated dataset
subclass_list <- unique(hrt.integrated_2$subclass)
centroids <- matrix(,nrow = length(subclass_list), ncol = 2)

for (i in 1:length(subclass_list)) {
  input.data <- data.frame(hrt.integrated_2@reductions$umap@cell.embeddings[hrt.integrated_2$subclass == subclass_list[i],"UMAP_1"],hrt.integrated_2@reductions$umap@cell.embeddings[hrt.integrated_2$subclass == subclass_list[i],"UMAP_2"])
  centroids[i,] <- colMeans(input.data)
}

dist_matrix <- as.matrix(dist(centroids, diag = F))
rownames(dist_matrix) <- subclass_list
colnames(dist_matrix) <- subclass_list
dist_matrix[dist_matrix == 0] = colMeans(dist_matrix)
min_values <- colMins(dist_matrix)
max_values <- colMaxs(dist_matrix)

max_markers <- c()
min_markers <- c()
subset.markers_spatial_class <- list()
subset.markers_spatial_subclass <- list()
Idents(hrt.integrated_2) <- "subclass"
DefaultAssay(hrt.integrated_2) <- "RNA"
for (i in 1:length(subclass_list)) {
  max_cluster <- subclass_list[which( dist_matrix[,i] == max_values[i])]
  min_cluster <- subclass_list[which( dist_matrix[,i] == min_values[i])]
  subset.markers_spatial_class[[i]] <- FindMarkers(hrt.integrated_2,ident.1 = subclass_list[i],ident.2 =as.character(max_cluster), only.pos = TRUE, min.pct = 0.25,logfc.threshold = 2)
  subset.markers_spatial_subclass[[i]] <- FindMarkers(hrt.integrated_2,ident.1 = subclass_list[i], only.pos = TRUE, min.pct = 0.25,logfc.threshold = 2)
}


lapply(subset.markers_spatial_subclass, function(x) length(rownames(x)))
#6,13, 16 had fewer than 10 DEGs with logfc > 2... so lower logfc for those closer to 1 
subclass_list[6]
i = 16 #iteratively adjusting 
max_cluster <- subclass_list[which( dist_matrix[,i] == max_values[i])]
min_cluster <- subclass_list[which( dist_matrix[,i] == min_values[i])]
subset.markers_spatial_subclass[[i]] <- FindMarkers(hrt.integrated_2,ident.1 = subclass_list[i], only.pos = TRUE, min.pct = 0.25,logfc.threshold = 1)
lapply(subset.markers_spatial_subclass, function(x) length(rownames(x)))

#sorting by logfc
subset.markers_spatial_class <- lapply(subset.markers_spatial_class, function(x) x[order(x$avg_log2FC,decreasing = T),])
subset.markers_spatial_subclass <- lapply(subset.markers_spatial_subclass, function(x) x[order(x$avg_log2FC,decreasing = T),])

#checking for high ^Rp 
lapply(subset.markers_spatial_subclass, function(x) length(grep(x, pattern = "^Rp")))

#removing Rp genes
x <- subset.markers_spatial_subclass[[15]]
grep(rownames(x), pattern = "^Rp")
subset.markers_spatial_subclass[[15]] <- subset.markers_spatial_subclass[[15]][-grep(rownames(x), pattern = "^Rp"),]  

save(subset.markers_spatial_class, file = "subset.markers_spatial_class.robj")
save(subset.markers_spatial_subclass, file = "subset.markers_spatial_subclass.robj")
load("/Users/home/Desktop/subset.markers_spatial_class.robj")
load("/Users/home/Desktop/subset.markers_spatial_subclass.robj")

#Mapping of CM gene-set scorers to visium clusters

class_threshold = 5
for (i in 1:length(subclass_list)){
  DefaultAssay(object_integrated) <- "Spatial"
  markers_spatial_subclass <- rownames(subset.markers_spatial_subclass[[i]])
  markers_spatial_class <- rownames(subset.markers_spatial_class[[i]])
  
  genes_in_spatial_and_cluster <- markers_spatial_subclass[markers_spatial_subclass %in% rownames(object_integrated)]
  print(i)
  print(length(genes_in_spatial_and_cluster))
  gene_counts <- GetAssayData(object_integrated, assay = "Spatial", slot = "counts")[genes_in_spatial_and_cluster[1:10],]
  
  genes_in_spatial_and_cluster_class <- markers_spatial_class[markers_spatial_class %in% rownames(object_integrated)]
  gene_counts_class <- GetAssayData(object_integrated, assay = "Spatial", slot = "counts")[genes_in_spatial_and_cluster_class[1:10],]
  
  subclass_score_temp <- Matrix::colSums(gene_counts)/object_integrated@meta.data$nCount_Spatial*10000
  class_score <- Matrix::colSums(gene_counts_class)/object_integrated@meta.data$nCount_Spatial*10000
  subclass_score_temp[class_score < class_threshold] <- 0
  
  object_integrated <- AddMetaData(object_integrated,subclass_score_temp, col.name = as.character(subclass_list[[i]]))
}

#renaming CM scores for convenience
colnames(object_integrated@meta.data)[colnames(object_integrated@meta.data) =="X38_CM_Myh6"] = "RZ_genes"
colnames(object_integrated@meta.data)[colnames(object_integrated@meta.data) =="X39_CM_Ankrd1"] = "BZ1_genes"
colnames(object_integrated@meta.data)[colnames(object_integrated@meta.data) =="X41_CM_Xirp2"] = "BZ2_genes"


#classification

zone_markers <- FindAllMarkers(object_integrated, only.pos = T, logfc.threshold = .25)
top_zone_markers <- zone_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
high_resolution_markers <- top_zone_markers
DoHeatmap(object_integrated, features = c(top_zone_markers$gene), 
          group.colors =  c("#7570B3",  "#BA3859"  , "#FF0000","darkgrey","darkgrey")) +  scale_fill_gradientn(colors = c("grey","white",brewer.pal(8,"Dark2")[3]) )

DoHeatmap(object_integrated, features = c(top_zone_markers$gene), 
          group.colors = c(brewer.pal(8,"Set2")[1:8],brewer.pal(8,"Accent")[1:8]) ) +  scale_fill_gradientn(colors = c("grey","white",brewer.pal(8,"Dark2")[3]) )

low_resolution_markers <- top_zone_markers

#CM Score
DefaultAssay(object_integrated) <- "Spatial"
i = 26
markers_spatial_subclass <- rownames(subset.markers_spatial_subclass[[i]])
markers_spatial_class <- rownames(subset.markers_spatial_class[[i]])

genes_in_spatial_and_cluster_class <- markers_spatial_class[markers_spatial_class %in% rownames(object_integrated)]
gene_counts_class <- GetAssayData(object_integrated, assay = "Spatial", slot = "counts")[genes_in_spatial_and_cluster_class[1:10],]

subclass_score_temp <- Matrix::colSums(gene_counts)/object_integrated@meta.data$nCount_Spatial*10000
class_score <- Matrix::colSums(gene_counts_class)/object_integrated@meta.data$nCount_Spatial*10000
subclass_score_temp[class_score < class_threshold] <- 0

object_integrated <- AddMetaData(object_integrated,class_score, col.name = "CM_Score")

#create new assay with gene set scores 
gene_scores <- CreateAssayObject(t(object_integrated@meta.data[,c(58:84,93)]))
object_integrated[["genescores"]] <- gene_scores

###CM or IZ?
Idents(object_integrated) <- "orig.ident"
object_integrated_sham <- subset(object_integrated, idents = "sham1")
object_integrated_D3_1 <- subset(object_integrated, idents =  "D3_1")
hist(object_integrated$CM_Score, breaks = 30)

data <- data.frame(type = c(rep("sham",length(object_integrated_sham$CM_Score)),rep("D3",length(object_integrated_D3_1$CM_Score))),
                   value = c(object_integrated_sham$CM_Score,object_integrated_D3_1$CM_Score))

#histograms of CM score to determine negative values/clusters for ROC analysis
p <- data %>%
  ggplot( aes(x=value, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c( "#FF3333","#00CCCC")) +
  labs(fill="")+xlim(0,200)+theme_light()+NoLegend()+theme( axis.title.y = element_blank())+ theme( axis.title.x = element_blank())

p
### ROC analysis
FindMarkers(object_integrated,ident.1 = 12, ident.2 = c(10,11), features = rownames(object_integrated)[28], test.use = "roc")
VlnPlot(object_integrated, "CM_Score", pt.size = -1,cols =c(brewer.pal(8,"Set2")[1:8],brewer.pal(8,"Accent")[1:8]))+theme_classic()+NoLegend()+ theme(axis.text.y = element_blank())+ theme(axis.text.x = element_blank())+ theme( axis.title.y = element_blank())+ theme( axis.title.x = element_blank())+ theme( tit   = element_text(size = 0)) +NoLegend()
Idents(object_integrated) <- "high_resolution"
object_integrated_CMs <- subset(object_integrated, idents = c(1:9,13)) #subseting CMs in new object


#RZ or BZ?
DefaultAssay(object_integrated)<-"genescores"
object_integrated_sham <- subset(object_integrated_CMs, idents = "sham1")
object_integrated_D3_1 <- subset(object_integrated_CMs, idents =  "D3_3")

data <- data.frame(type = c(rep("sham",length(object_integrated_sham$BZ1_genes)),rep("D3",length(object_integrated_D3_1$BZ1_genes))),
                   value = c(object_integrated_sham$BZ1_genes,object_integrated_D3_1$BZ1_genes))

#histogram of BZ1 scores (sham v D3) to define negative for ROC analysis
p <- data %>%
  ggplot( aes(x=value, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c( "#FF3333","#00CCCC")) +
  labs(fill="")+xlim(0,75)+theme_light()+NoLegend()+theme( axis.title.y = element_blank())+ theme( axis.title.x = element_blank())

p

VlnPlot(object_integrated_CMs, "BZ1_genes")

FindMarkers(object_integrated,ident.1 = 12, ident.2 = c(5), features = rownames(object_integrated)[25], test.use = "roc")

#subsetting BZ clusters in new object
object_integrated_CMs_BZ <- subset(object_integrated_CMs, idents = c(5,6,7,8,9,13))

Idents(object_integrated_CMs_BZ) <- "orig.ident"
bject_integrated_sham <- subset(object_integrated_CMs_BZ, idents = "sham1")
object_integrated_D3_1 <- subset(object_integrated_CMs_BZ, idents =  "D3_3")

data <- data.frame(type = c(rep("D3",length(object_integrated_D3_1$BZ2_genes))),
                   value = c(object_integrated_D3_1$BZ2_genes))

# Histogram
p <- data %>%
  ggplot( aes(x=value, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c( "#FF3333","#00CCCC")) +
  labs(fill="")+xlim(25,100)+theme_light()+NoLegend()+theme( axis.title.y = element_blank())+ theme( axis.title.x = element_blank())

p

#produce scatter plot comparing BZ1 and BZ2 scores in D3 sample
data_2 <- object_integrated_D3_1@meta.data
ggplot(data_2, aes(x=BZ1_genes, y=BZ2_genes))+ylim(25,100) + geom_point(alpha = .7, size = 1)+theme_light()+NoLegend()+theme( axis.title.y = element_blank())+ theme( axis.title.x = element_blank())

#ROC analysis
FindMarkers(object_integrated,ident.1 = 9, ident.2 = 5, features = rownames(object_integrated)[28], test.use = "roc")

###annotating spatial regions
object_integrated<-RenameIdents(object_integrated, '9'=1,'' = 1,'1' = 1,'2' = 1, 
                                '5' = 2,'7' = 2, '11' =2, '4'=2,
                                '6' = 3, '8' =3,
                                '3' = 4,'10' = 4, '9' = 4)

object_integrated <- StashIdent(object_integrated, save.name = "annotated")


###AUC as function of gene-set length
#creating a matrix of RZ scores.  Each row representing number of genes (1-20)
RZ_scores_gene_length <- matrix(,nrow = 20, ncol = length(colnames(object_integrated)))

for (x in 2:20){
  DefaultAssay(object_integrated) <- "Spatial"
  markers_spatial_subclass <- rownames(subset.markers_spatial_subclass[[i]])
  markers_spatial_class <- rownames(subset.markers_spatial_class[[i]])
  
  genes_in_spatial_and_cluster <- markers_spatial_subclass[markers_spatial_subclass %in% rownames(object_integrated)]
  print(i)
  print(length(genes_in_spatial_and_cluster))
  print(x)
  print(length(genes_in_spatial_and_cluster[1:x]))
  gene_counts <- GetAssayData(object_integrated, assay = "Spatial", slot = "counts")[genes_in_spatial_and_cluster[1:x],]
  subclass_score_temp <- Matrix::colSums(gene_counts)/object_integrated@meta.data$nCount_Spatial*10000
  
  genes_in_spatial_and_cluster_class <- markers_spatial_class[markers_spatial_class %in% rownames(object_integrated)]
  gene_counts_class <- GetAssayData(object_integrated, assay = "Spatial", slot = "counts")[genes_in_spatial_and_cluster_class[1:10],]
  
  class_score <- Matrix::colSums(gene_counts_class)/object_integrated@meta.data$nCount_Spatial*10000
  subclass_score_temp[class_score < class_threshold] <- 0
  RZ_scores_gene_length[x-1,] <- subclass_score_temp
  
  #object_integrated <- AddMetaData(object_integrated,subclass_score_temp, col.name = paste("RZ_score_genes",x,sep = ""))
  
}

head(RZ_scores_gene_length)[1:5,1:5]
colnames(RZ_scores_gene_length) <- colnames(object_integrated)
rownames(RZ_scores_gene_length)  <- 1:20
dim(object_integrated)
dim(RZ_scores_gene_length)
RZ_scores_by_genes <- CreateAssayObject(RZ_scores_gene_length)
object_integrated[["RZ_scores_gene_length"]] <- RZ_scores_by_genes
DefaultAssay(object_integrated) <- "RZ_scores_gene_length"

RZ_ROC_results <- list()
for (x in 1:19) {
  ROC_results_temp <- FindMarkers(object_integrated, ident.1 = 1, ident.2 = 10:11,features = rownames(object_integrated)[x] ,test.use = "roc")
  RZ_ROC_results[x] <- ROC_results_temp$myAUC
}

RZ_ROC_results <- unlist(RZ_ROC_results)
RZ_ROC_results


#BZ1 now
i <- 27
BZ1_scores_gene_length <- matrix(,nrow = 20, ncol = length(colnames(object_integrated)))

for (x in 2:20){
  DefaultAssay(object_integrated) <- "Spatial"
  markers_spatial_subclass <- rownames(subset.markers_spatial_subclass[[i]])
  markers_spatial_class <- rownames(subset.markers_spatial_class[[i]])
  
  genes_in_spatial_and_cluster <- markers_spatial_subclass[markers_spatial_subclass %in% rownames(object_integrated)]
  print(i)
  print(length(genes_in_spatial_and_cluster))
  print(x)
  print(length(genes_in_spatial_and_cluster[1:x]))
  gene_counts <- GetAssayData(object_integrated, assay = "Spatial", slot = "counts")[genes_in_spatial_and_cluster[1:x],]
  subclass_score_temp <- Matrix::colSums(gene_counts)/object_integrated@meta.data$nCount_Spatial*10000
  
  genes_in_spatial_and_cluster_class <- markers_spatial_class[markers_spatial_class %in% rownames(object_integrated)]
  gene_counts_class <- GetAssayData(object_integrated, assay = "Spatial", slot = "counts")[genes_in_spatial_and_cluster_class[1:10],]
  
  class_score <- Matrix::colSums(gene_counts_class)/object_integrated@meta.data$nCount_Spatial*10000
  subclass_score_temp[class_score < class_threshold] <- 0
  BZ1_scores_gene_length[x-1,] <- subclass_score_temp
  
  #object_integrated <- AddMetaData(object_integrated,subclass_score_temp, col.name = paste("BZ1_score_genes",x,sep = ""))
  
}

head(BZ1_scores_gene_length)[1:5,1:5]
colnames(BZ1_scores_gene_length) <- colnames(object_integrated)
rownames(BZ1_scores_gene_length)  <- 1:20
dim(object_integrated)
dim(BZ1_scores_gene_length)
BZ1_scores_by_genes <- CreateAssayObject(BZ1_scores_gene_length)
object_integrated[["BZ1_scores_gene_length"]] <- BZ1_scores_by_genes
DefaultAssay(object_integrated) <- "BZ1_scores_gene_length"


FeaturePlot(object_integrated,features = rownames(object_integrated)[1:4]  )


BZ1_ROC_results <- list()
for (x in 1:19) {
  ROC_results_temp <- FindMarkers(object_integrated, ident.1 = 5, ident.2 = 1,features = rownames(object_integrated)[x] ,test.use = "roc")
  BZ1_ROC_results[x] <- ROC_results_temp$myAUC
}

BZ1_ROC_results <- unlist(BZ1_ROC_results)
BZ1_ROC_results

#BZ2 now
i <- 25
BZ2_scores_gene_length <- matrix(,nrow = 20, ncol = length(colnames(object_integrated)))

for (x in 2:20){
  DefaultAssay(object_integrated) <- "Spatial"
  markers_spatial_subclass <- rownames(subset.markers_spatial_subclass[[i]])
  markers_spatial_class <- rownames(subset.markers_spatial_class[[i]])
  
  genes_in_spatial_and_cluster <- markers_spatial_subclass[markers_spatial_subclass %in% rownames(object_integrated)]
  print(i)
  print(length(genes_in_spatial_and_cluster))
  print(x)
  print(length(genes_in_spatial_and_cluster[1:x]))
  gene_counts <- GetAssayData(object_integrated, assay = "Spatial", slot = "counts")[genes_in_spatial_and_cluster[1:x],]
  subclass_score_temp <- Matrix::colSums(gene_counts)/object_integrated@meta.data$nCount_Spatial*10000
  
  genes_in_spatial_and_cluster_class <- markers_spatial_class[markers_spatial_class %in% rownames(object_integrated)]
  gene_counts_class <- GetAssayData(object_integrated, assay = "Spatial", slot = "counts")[genes_in_spatial_and_cluster_class[1:10],]
  
  class_score <- Matrix::colSums(gene_counts_class)/object_integrated@meta.data$nCount_Spatial*10000
  subclass_score_temp[class_score < class_threshold] <- 0
  BZ2_scores_gene_length[x-1,] <- subclass_score_temp
  
  #object_integrated <- AddMetaData(object_integrated,subclass_score_temp, col.name = paste("BZ2_score_genes",x,sep = ""))
  
}

head(BZ2_scores_gene_length)[1:5,1:5]
colnames(BZ2_scores_gene_length) <- colnames(object_integrated)
rownames(BZ2_scores_gene_length)  <- 1:20
dim(object_integrated)
dim(BZ2_scores_gene_length)
BZ2_scores_by_genes <- CreateAssayObject(BZ2_scores_gene_length)
object_integrated[["BZ2_scores_gene_length"]] <- BZ2_scores_by_genes
DefaultAssay(object_integrated) <- "BZ2_scores_gene_length"



BZ2_ROC_results <- list()
for (x in 1:19) {
  ROC_results_temp <- FindMarkers(object_integrated, ident.1 = 9, ident.2 = 5,features = rownames(object_integrated)[x] ,test.use = "roc")
  BZ2_ROC_results[x] <- ROC_results_temp$myAUC
}

BZ2_ROC_results <- unlist(BZ2_ROC_results)
BZ2_ROC_results
object_integrated$high_resolution
Idents(object_integrated) <-"high_resolution"
UMAPPlot(object_integrated, label = T)
object_integrated<-RenameIdents(object_integrated, '1'=1,'2' = 1,'3' = 1,'4' = 1, 
                                '5' = 2,'6' = 2, '7' =2, '8'=2,
                                '9' = 3,
                                '10' =4,'11'=4,'12'=4)
object_integrated <- StashIdent(object_integrated, save.name = "annotated")

save(object_integrated, file = "object_integrated_final.robj")
load("object_integrated_rev2.robj")

