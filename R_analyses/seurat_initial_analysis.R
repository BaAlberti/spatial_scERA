# Libraries
library(Matrix)
library(Seurat)
library(ggplot2)
library(viridisLite)
library(scCustomize)
library(dplyr)


# Preprocessing of all three replicates 
##Creating seurat objects and computing ribosomal and mitochondrial gene expression percentages
rawdata_61 <-Read10X(data.dir = "../data/scRNA-seq_matrix/filtered_feature_bc_matrix_61/filtered_feature_bc_matrix_61/")
data_61 <-CreateSeuratObject(counts= rawdata_61, min.cells =3, min.genes =200, project ="rep1") 
data_61 # 24886 cells
data_61[["percent.rb"]] <- PercentageFeatureSet(data_61, pattern="m?Rp|28SrRNA|18SrRNA")
data_61[["percent.mt"]] <- PercentageFeatureSet(data_61, pattern="mt:")

rawdata_62 <-Read10X(data.dir = "../data/scRNA-seq_matrix/filtered_feature_bc_matrix_62/filtered_feature_bc_matrix_62/")
data_62 <-CreateSeuratObject(counts= rawdata_62, min.cells =3, min.genes =200, project ="rep2") 
data_62 # 23125 cells
data_62[["percent.rb"]] <- PercentageFeatureSet(data_62, pattern="m?Rp|28SrRNA|18SrRNA")
data_62[["percent.mt"]] <- PercentageFeatureSet(data_62, pattern="mt:")

rawdata_63 <-Read10X(data.dir = "../data/scRNA-seq_matrix/filtered_feature_bc_matrix_63/filtered_feature_bc_matrix_63/")
data_63 <-CreateSeuratObject(counts= rawdata_63, min.cells =3, min.genes =200, project ="rep3") 
data_63 # 23437 cells
data_63[["percent.rb"]] <- PercentageFeatureSet(data_63, pattern="m?Rp|28SrRNA|18SrRNA")
data_63[["percent.mt"]] <- PercentageFeatureSet(data_63, pattern="mt:")

# Quality control values to select only what I consider a real cell
minGene = 500
maxGene = 6000
minUmi = 2000
maxUmi = 50000
maxpct_rb = 36 # Going lower will remove amnioserosa cells from the analysis
maxpct_mt = 15

# Subsetting to keep only the "real" cells
data_filtered_61 <- subset(data_61, subset=nFeature_RNA > minGene & nFeature_RNA < maxGene & nCount_RNA > minUmi & nCount_RNA < maxUmi & percent.mt < maxpct_mt & percent.rb < maxpct_rb)
data_filtered_61 #1934 cells
data_filtered_62 <- subset(data_62, subset=nFeature_RNA > minGene & nFeature_RNA < maxGene & nCount_RNA > minUmi & nCount_RNA < maxUmi & percent.mt < maxpct_mt & percent.rb < maxpct_rb)
data_filtered_62 #1963 cells
data_filtered_63 <- subset(data_63, subset=nFeature_RNA > minGene & nFeature_RNA < maxGene & nCount_RNA > minUmi & nCount_RNA < maxUmi & percent.mt < maxpct_mt & percent.rb < maxpct_rb)
data_filtered_63 #2691 cells

# Adding a "replicate" meta data column
data_filtered_61$replicate=1
data_filtered_62$replicate=2
data_filtered_63$replicate=3

# Merging all three replicates together and adding a replicate ID in front of each cell to keep track of their replicate of origin
data6_filtre_list <- merge(x = data_filtered_61, y =c(data_filtered_62,data_filtered_63), add.cell.ids=c("rep1","rep2","rep3"), project="stage6")
data6_filtre_list <- SplitObject(data6_filtre_list, split.by="replicate")

# Free some space in the memory
rm(rawdata_61,rawdata_62,rawdata_63)
rm(data_61,data_62,data_63)

# After merging the first processing steps are done on each replicate separatly 
for (i in 1:length(data6_filtre_list)){
  data6_filtre_list[[i]] <- NormalizeData(data6_filtre_list[[i]])
  data6_filtre_list[[i]] <- FindVariableFeatures(data6_filtre_list[[i]], selection.method = "mvp")
}

# Integrating the replicates together / removing batch effect 
data6_filtre_list.anchors <- FindIntegrationAnchors(object.list = data6_filtre_list,dims = 1:30)
data6_filtre_list.integrated <- IntegrateData(anchorset = data6_filtre_list.anchors, dims = 1:30)

# scRNA-seq analysis pipeline is done on the integrated assay of the seurat object.
DefaultAssay(data6_filtre_list.integrated) <- "integrated"
data6_filtre_list.integrated <- ScaleData(data6_filtre_list.integrated)
data6_filtre_list.integrated <- RunPCA(data6_filtre_list.integrated, npcs=30)
ElbowPlot(data6_filtre_list.integrated, ndims=30) # I choose 30 as the curve is realy flat after that
data6_filtre_list.integrated <- FindNeighbors(data6_filtre_list.integrated, dims = 1:30)
data6_filtre_list.integrated <- FindClusters(data6_filtre_list.integrated, resolution = 0.5)
data6_filtre_list.integrated <- RunUMAP(data6_filtre_list.integrated, dims=1:30)

# Changing cells name and saving them and the matrix for PCR and novosparc analyses
id_cells_6_int <- Cells(data6_filtre_list.integrated)
id_cells_6_int <- gsub("-1","",as.character(id_cells_6_int))
data6_filtre_list.integrated <- RenameCells(data6_filtre_list.integrated,new.names=id_cells_6_int)
write.csv(id_cells_6_int, file="../data/seurat_outputs/id_cells_6_int_full_clusters.csv", row.names = FALSE)
write.table(data6_filtre_list.integrated@assays$RNA@counts,file = "../data/seurat_outputs/matrice_6_int_full_clusters.csv",quote=FALSE,sep=',',col.names=TRUE)

# Saving the seurat object before annotation
saveRDS(data6_filtre_list.integrated,file="../data/seurat_objects/data6_fullclustering_no_label.RDS")

# Identifying the different cell types in the clusters
cluster.markers_6 <- FindAllMarkers(data6_filtre_list.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
list_10_gene=c()
for (i in 0:15){
  clusters=subset(cluster.markers_6,cluster.markers_6$cluster==i)
  top10=head(clusters,10)
  list_10_gene=c(list_10_gene,top10[,7])
}
write.table(list_10_gene, file="../data/seurat_outputs/list_markers_6_full_version.csv",row.names=F,col.names=F)

# #Annotation for each cluster (do not run if you want to label the cluster as in the paper)
# data6_filtre_list.integrated <- RenameIdents(data6_filtre_list.integrated, 
#                                              `0` = "ecto1", 
#                                              `1` = "meso1", 
#                                              `2` = "ecto2",
#                                              `3` = "unknown1", 
#                                              `4` = "amnio1", 
#                                              `5` = "pro ecto1", 
#                                              `6` = "pro ecto2", 
#                                              `7` = "mesecto1", 
#                                              `8` = "dors ecto1",
#                                              `9` = "yolk",
#                                              `10` = "endo1",
#                                              `11` = "dors ecto2",
#                                              `12` = "unknown2",
#                                              `13` = "meso2",
#                                              `14` = "pole cells",
#                                              `15` = "hemocyte") 
# DimPlot(data6_filtre_list.integrated, reduction = "umap", label=TRUE, label.size = 5, pt.size = 1.15,repel = TRUE,label.box = TRUE)

# # Collecting all cells ID for each cluster to plot them on novosparc later
# Id_cells_per_clusters <- data.frame(matrix(ncol=2,nrow=0))
# names(Id_cells_per_clusters) <- c("cluster","sc_cells")
# for (clust in levels(data6_filtre_list.integrated@active.ident)) {
#   temp_data <- subset(data6_filtre_list.integrated,idents=clust)
#   list_cell <- paste(Cells(temp_data),collapse=",")
#   Id_cells_per_clusters[nrow(Id_cells_per_clusters) +1,] = c(clust,list_cell)
# }
# write.csv(x = Id_cells_per_clusters,file = "Id_cells_per_clusters_full_no_label.csv",row.names=FALSE)

# Paper' annotation
data6_filtre_list.integrated <- RenameIdents(data6_filtre_list.integrated, 
                                             `0` = "ectoderm", 
                                             `2` = "ectoderm", 
                                             `3` = "unknown",
                                             `12` = "unknown",
                                             `5` = "procephalic region",
                                             `6` = "procephalic region",
                                             `8` = "dorsal ectoderm", 
                                             `11` = "dorsal ectoderm",
                                             `1` = "mesoderm",
                                             `13` = "mesoderm",
                                             `4` = "amnioserosa",
                                             `10` = "endoderm",
                                             `7` = "mesectoderm",
                                             `9` = "yolk",
                                             `14` = "pole cells",
                                             `15` = "hemocytes")

# Saving UMAPs
DimPlot_scCustom(data6_filtre_list.integrated, reduction = "umap", label=TRUE, label.size = 5, pt.size = 1.15,repel = TRUE,label.box = TRUE,ggplot_default_colors = T)


# Collecting all cells ID for each cluster to plot them on novosparc later
Id_cells_per_clusters <- data.frame(matrix(ncol=2,nrow=0))
names(Id_cells_per_clusters) <- c("cluster","sc_cells")
for (clust in levels(data6_filtre_list.integrated@active.ident)) {
  temp_data <- subset(data6_filtre_list.integrated,idents=clust)
  list_cell <- paste(Cells(temp_data),collapse=",")
  Id_cells_per_clusters[nrow(Id_cells_per_clusters) +1,] = c(clust,list_cell)
}
write.csv(x = Id_cells_per_clusters,file = "../data/seurat_outputs/Id_cells_per_clusters_full.csv",row.names=FALSE)

# Saving the seurat object before annotation
saveRDS(data6_filtre_list.integrated,file="../data/seurat_objects/data6_fullclustering_with_label.RDS")


### Post processing of the full version of the clustering
# Number of cell per enhancer
DefaultAssay(data6_filtre_list.integrated) <- "RNA"
data6_filtre_list.integrated[["CellName"]] <- colnames(data6_filtre_list.integrated)
list_enhancer=c("eve-late-variant","h-stripe1","salm-blastoderm-early-enhancer","twi-CHIP-42","Vnd-743","CHIP-27","prd01","pscE14","shg-A",
                "SoxN-5830","GMR77A12","GMR83E01","CRM1","CRM2","CRM3","CRM4","CRM5","CRM6","CRM7","CRM8","CRM9","CRM10","CRM11","CRM12","CRM13")

Id_cells_with_enhancers_sc <- data.frame(matrix(ncol=3,nrow=0))
names(Id_cells_with_enhancers_sc) <- c("enhancer","sc_cells_sc_only","nb_cells_sc_only")
for (elem in list_enhancer){
  print(elem)
  res <- try(FetchData(data6_filtre_list.integrated,vars = elem))
  if (inherits(res, "try-error")){Id_cells_with_enhancers_sc[nrow(Id_cells_with_enhancers_sc)+1,] = c(elem,'',0)}
  else{
    expr <- FetchData(data6_filtre_list.integrated,vars = elem)
    res2 <- try(Cells(data6_filtre_list.integrated[,which(x=expr > 0)]))
    if (inherits(res2,"try-error")){Id_cells_with_enhancers_sc[nrow(Id_cells_with_enhancers_sc)+1,] = c(elem,'',0)}
    else{
      list_cell <- Cells(data6_filtre_list.integrated[,which(x=expr > 0)])
      Id_cells_with_enhancers_sc[nrow(Id_cells_with_enhancers_sc) +1,] = c(elem,paste(list_cell,collapse = ","),length(list_cell))
    }
  }
}
write.csv(x = Id_cells_with_enhancers_sc,file = "../data/seurat_outputs/ID_cells_with_enhancers_sc_fullclustering.csv")

# Computing the number of cell per enhancer per sc cluster from sc only
for (elem in list_enhancer){
  print(elem)
  if (elem %in% c("CRM5","CRM8")){next}
  expr <- FetchData(data6_filtre_list.integrated,vars = elem)
  list_sc <- Cells(data6_filtre_list.integrated[,which(x=expr > 0)])
  test_subset <- subset(data6_filtre_list.integrated, subset = CellName %in% list_sc)
  print(table(Idents(test_subset)))
}

## Switching to the reduced version of the dataset where three clusters are removed as they are not at the surface of the embryo

# Subsetting to remove the three "internal" clusters
data6_wo_rb <- subset(data6_filtre_list.integrated,idents=c("unknown","yolk","hemocytes"), invert=TRUE)

# redo the analysis
DefaultAssay(data6_wo_rb) <- "integrated"
data6_wo_rb <- ScaleData(data6_wo_rb)
data6_wo_rb <- RunPCA(data6_wo_rb, features = VariableFeatures(object = data6_wo_rb))
ElbowPlot(data6_wo_rb, ndims = 30)+xlab("Principal Components")+theme(axis.title=element_text(size=16))
data6_wo_rb <- FindNeighbors(data6_wo_rb, dims = 1:30)
data6_wo_rb <- FindClusters(data6_wo_rb, resolution = 0.5)
data6_wo_rb <- RunUMAP(data6_wo_rb, dims = 1:30,min.dist=0.5)

# Changing cells name and saving them and the matrix for PCR and novosparc analyses
id_cells_6_int <- Cells(data6_wo_rb)
id_cells_6_int <- gsub("-1","",as.character(id_cells_6_int))
data6_wo_rb <- RenameCells(data6_wo_rb,new.names=id_cells_6_int)
write.csv(id_cells_6_int, file="../data/seurat_outputs/id_cells_6_int_reduced_clusters.csv", row.names = FALSE)
write.table(data6_wo_rb@assays$RNA@counts,file = "../data/seurat_outputs/matrice_6_int_reduced_clusters.csv",quote=FALSE,sep=',',col.names=TRUE)

# Saving the seurat object before annotation
saveRDS(data6_wo_rb,file="../data/seurat_objects/data6_reducedclustering_no_label.RDS")

# Identifying the different cell types in this new dataset
cluster.markers_6 <- FindAllMarkers(data6_wo_rb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
list_10_gene=c()
for (i in 0:12){
  clusters=subset(cluster.markers_6,cluster.markers_6$cluster==i)
  top10=head(clusters,10)
  list_10_gene=c(list_10_gene,top10[,7])
}
write.table(list_10_gene, file="list_markers_6_reduced version.csv",row.names=F,col.names=F)

# #Annotating all clusters individually (do not run if you want to obtain the same cell-type as in the paper)
# data6_wo_rb <- RenameIdents(data6_wo_rb,     
#                             `0` = "ecto1", 
#                             `1` = "meso1", 
#                             `2` = "ecto2",
#                             `3` = "amnio", 
#                             `4` = "pro ecto1", 
#                             `5` = "pro ecto2", 
#                             `6` = "endo1", 
#                             `7` = "dors ecto1", 
#                             `8` = "dors ecto2",
#                             `9` = "meso2",
#                             `10` = "head meso",
#                             `11` = "mesecto",
#                             `12` = "pole cells")
# DimPlot(data6_wo_rb, reduction = "umap", label=TRUE, label.size = 5, pt.size = 1.15,repel = TRUE,label.box = TRUE) +NoLegend()

# # Collecting all cells ID for each cluster to plot them on novosparc later
# Id_cells_per_clusters <- data.frame(matrix(ncol=2,nrow=0))
# names(Id_cells_per_clusters) <- c("cluster","sc_cells")
# for (clust in levels(data6_wo_rb@active.ident)) {
#   temp_data <- subset(data6_wo_rb,idents=clust)
#   list_cell <- paste(Cells(temp_data),collapse=",")
#   Id_cells_per_clusters[nrow(Id_cells_per_clusters) +1,] = c(clust,list_cell)
# }
# write.csv(x = Id_cells_per_clusters,file = "Id_cells_per_clusters_reduced_no_label.csv",row.names=FALSE)

# Paper' annotation
data6_wo_rb <- RenameIdents(data6_wo_rb,     
                            `0` = "ectoderm", 
                            `2` = "ectoderm", 
                            `1` = "mesoderm",
                            `9` = "mesoderm",
                            `3` = "amnioserosa",
                            `4` = "procephalic region",
                            `5` = "procephalic region",
                            `7` = "dorsal ectoderm", 
                            `8` = "dorsal ectoderm",
                            `6` = "endoderm",
                            `10` = "head mesoderm",
                            `11` = "mesectoderm",
                            `12` = "pole cells")

# Saving UMAPs
DimPlot_scCustom(data6_wo_rb, reduction = "umap", label=TRUE, label.size = 5, pt.size = 1.15,repel = TRUE,label.box = TRUE,ggplot_default_colors = T)

saveRDS(object = data6_wo_rb, file="../data/seurat_objects/data6_reducedclustering_with_label.RDS")

# Collecting all cells ID for each cluster to plot them on novosparc later
Id_cells_per_clusters <- data.frame(matrix(ncol=2,nrow=0))
names(Id_cells_per_clusters) <- c("cluster","sc_cells")
for (clust in levels(data6_wo_rb@active.ident)) {
  temp_data <- subset(data6_wo_rb,idents=clust)
  list_cell <- paste(Cells(temp_data),collapse=",")
  Id_cells_per_clusters[nrow(Id_cells_per_clusters) +1,] = c(clust,list_cell)
}
write.csv(x = Id_cells_per_clusters,file = "../data/seurat_outputs/Id_cells_per_clusters_reduced.csv",row.names=FALSE)

### Post processing
# Number of cell per enhancer
DefaultAssay(data6_wo_rb) <- "RNA"
list_enhancer=c("eve-late-variant","h-stripe1","salm-blastoderm-early-enhancer","twi-CHIP-42","Vnd-743","CHIP-27","prd01","pscE14","shg-A",
                "SoxN-5830","GMR77A12","GMR83E01","CRM1","CRM2","CRM3","CRM4","CRM5","CRM6","CRM7","CRM8","CRM9","CRM10","CRM11","CRM12","CRM13")

Id_cells_with_enhancers_sc <- data.frame(matrix(ncol=3,nrow=0))
names(Id_cells_with_enhancers_sc) <- c("enhancer","sc_cells_sc_only","nb_cells_sc_only")
for (elem in list_enhancer){
  print(elem)
  res <- try(FetchData(data6_wo_rb,vars = elem))
  if (inherits(res, "try-error")){Id_cells_with_enhancers_sc[nrow(Id_cells_with_enhancers_sc)+1,] = c(elem,'',0)}
  else{
    expr <- FetchData(data6_wo_rb,vars = elem)
    res2 <- try(Cells(data6_wo_rb[,which(x=expr > 0)]))
    if (inherits(res2,"try-error")){Id_cells_with_enhancers_sc[nrow(Id_cells_with_enhancers_sc)+1,] = c(elem,'',0)}
    else{
      list_cell <- Cells(data6_wo_rb[,which(x=expr > 0)])
      Id_cells_with_enhancers_sc[nrow(Id_cells_with_enhancers_sc) +1,] = c(elem,paste(list_cell,collapse = ","),length(list_cell))
    }
  }
}
write.csv(x = Id_cells_with_enhancers_sc,file = "../data/seurat_outputs/ID_cells_with_enhancers_sc_reduced_clustering.csv")

for (elem in list_enhancer){
  print(elem)
  if (elem %in% c("CRM5","CRM8","pscE14")){next}
  expr <- FetchData(data6_filtre_list.integrated,vars = elem)
  list_sc <- Cells(data6_filtre_list.integrated[,which(x=expr > 0)])
  test_subset <- subset(data6_filtre_list.integrated, subset = CellName %in% list_sc)
  print(table(Idents(test_subset)))
}

