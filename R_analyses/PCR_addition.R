# Libraries
library(Matrix)
library(Seurat)
library(ggplot2)
library(viridisLite)
library(scCustomize)
library(dplyr)

# Loading Seurat Objects
data6_full <- readRDS("../data/seurat_objects/data6_fullclustering_with_label.RDS")
data6_reduced <- readRDS("../data/seurat_objects/data6_reducedclustering_with_label.RDS")
data6_full[["CellName"]] <- colnames(data6_full)
data6_reduced[["CellName"]] <- colnames(data6_reduced)
DefaultAssay(data6_full) <- "RNA"
DefaultAssay(data6_reduced) <- "RNA"

list_enhancer=c("eve-late-variant","h-stripe1","salm-blastoderm-early-enhancer","twi-CHIP-42","Vnd-743","CHIP-27","prd01","pscE14","shg-A",
                "SoxN-5830","GMR77A12","GMR83E01","CRM1","CRM2","CRM3","CRM4","CRM5","CRM6","CRM7","CRM8","CRM9","CRM10","CRM11","CRM12","CRM13")


#### Comparing PCR discoveries against scRNA-seq results
## In every for loop forward you can change both the PCR list of cells but depending on the one you choose, you have to change the scdataset in the loop
## ( full -> data6_filtre_list.intgrated / reduced -> data6_reduced)

df_all_list_PCR_full <- read.csv("../python_analysis/trimming_PCR/merged_full/lists_cells_for_each_enhancer_full_version_merged.tsv",sep="\t", header=T)
#df_61_list_PCR_full <- read.csv("../python_analysis/trimming_PCR/individual_full/lists_cells_for_each_enhancer_full_version_6_1.tsv",sep="\t", header=T)
#df_62_list_PCR_full <- read.csv("../python_analysis/trimming_PCR/individual_full/lists_cells_for_each_enhancer_full_version_6_2.tsv",sep="\t", header=T)
#df_63_list_PCR_full <- read.csv("../python_analysis/trimming_PCR/individual_full/lists_cells_for_each_enhancer_full_version_6_3.tsv",sep="\t", header=T)
df_all_list_PCR_reduced <- read.csv("../python_analysis/trimming_PCR/merged_reduced/lists_cells_for_each_enhancer_reduced_version_merged.tsv",sep="\t",header=T)
#df_61_list_PCR_reduced <- read.csv("../python_analysis/trimming_PCR/individual_full/lists_cells_for_each_enhancer_reduced_version_6_1.tsv",sep="\t", header=T)
#df_62_list_PCR_reduced <- read.csv("../python_analysis/trimming_PCR/individual_full/lists_cells_for_each_enhancer_reduced_version_6_2.tsv",sep="\t", header=T)
#df_63_list_PCR_reduced <- read.csv("../python_analysis/trimming_PCR/individual_full/lists_cells_for_each_enhancer_reduced_version_6_3.tsv",sep="\t", header=T)

# computing the number of cell discovered by the PCR alone (full)
nb_cell_PCR_alone_full <- data.frame(matrix(ncol=2,nrow=0))
names(nb_cell_PCR_alone_full) <- c("enhancer","nb_cell_PCR_only")
for (elem in df_all_list_PCR_full$enhancer){
  print(elem)
  nb_cell <- length(unlist(strsplit(df_all_list_PCR_full$liste[df_all_list_PCR_full$enhancer==elem],",")))
  nb_cell_PCR_alone_full[nrow(nb_cell_PCR_alone_full) +1,] = c(elem,length(nb_cell))
}
write.csv(x = nb_cell_PCR_alone_full,file = "nb_cell_PCR_alone_full.csv")

# computing the number of cell discovered by the PCR alone (reduced)
nb_cell_PCR_alone_reduced <- data.frame(matrix(ncol=2,nrow=0))
names(nb_cell_PCR_alone_reduced) <- c("enhancer","nb_cell_PCR_only")
for (elem in df_all_list_PCR_reduced$enhancer){
  print(elem)
  nb_cell <- length(unlist(strsplit(df_all_list_PCR_reduced$liste[df_all_list_PCR_reduced$enhancer==elem],",")))
  nb_cell_PCR_alone_reduced[nrow(nb_cell_PCR_alone_reduced) +1,] = c(elem,length(nb_cell))
}
write.csv(x = nb_cell_PCR_alone_reduced,file = "nb_cell_PCR_alone_reduced.csv")


for (elem in df_all_list_PCR_full$enhancer){
  print(elem)
  print(length(unlist(strsplit(df_all_list_PCR_full$liste[df_all_list_PCR_full$enhancer==elem],","))))
}

# Computing the number of cell per enhancer per sc cluster from the PCR only (full)
for (elem in df_all_list_PCR_full$enhancer){
  print(elem)
  spe_cell <- unlist(strsplit(df_all_list_PCR_full$liste[df_all_list_PCR_full$enhancer==elem],","))
  test_subset <- subset(data6_full, subset = CellName %in% spe_cell)
  print(table(Idents(test_subset)))
}

# Computing the number of cell per enhancer per sc cluster from the PCR only (reduced)
for (elem in df_all_list_PCR_reduced$enhancer){
  print(elem)
  spe_cell <- unlist(strsplit(df_all_list_PCR_reduced$liste[df_all_list_PCR_reduced$enhancer==elem],","))
  test_subset <- subset(data6_reduced, subset = CellName %in% spe_cell)
  print(table(Idents(test_subset)))
}

# Computing the number of cell per enhancer per sc cluster after merging unique cell barcodes from both PCR and sc analyses
# Also make a table with the common cells between scRNA-seq and PCR (full)
nb_cell_common_full <- data.frame(matrix(ncol=2,nrow=0))
names(nb_cell_common_full) <- c("enhancer","nb_cell_common")
for (elem in list_enhancer){
  print(elem)
  res <- try(FetchData(data6_full,vars = elem))
  if (inherits(res, "try-error")){list_sc <- c()}
  else{
    expr <- FetchData(data6_full,vars = elem)
    res <- try(Cells(data6_full[,which(x=expr > 0)]))
    if (inherits(res, "try-error")){list_sc <- c()}
    else{list_sc <- Cells(data6_full[,which(x=expr > 0)])}
  }
  list_pcr <- df_all_list_PCR_full$liste[df_all_list_PCR_full$enhancer==elem]
  nb_cell <- length(intersect(list_sc,unlist(strsplit(df_all_list_PCR_full$liste[df_all_list_PCR_full$enhancer==elem],","))))
  nb_cell_common_full[nrow(nb_cell_common_full) +1,] = c(elem,nb_cell)
  double_list <- paste(list_sc,list_pcr,sep=",")
  unique_list <- unique(unlist(strsplit(double_list,",")))
  if (is.null(unique_list)){next}
  else{
    test_subset <- subset(data6_full, subset = CellName %in% unique_list)
    print(table(Idents(test_subset)))
  }
}
write.csv(x=nb_cell_common_full,file="nb_cell_common_full.csv")

# Computing the number of cell per enhancer per sc cluster after merging unique cell barcodes from both PCR and sc analyses
# Also make a table with the common cells between scRNA-seq and PCR (full)
nb_cell_common_reduced <- data.frame(matrix(ncol=2,nrow=0))
names(nb_cell_common_reduced) <- c("enhancer","nb_cell_common")
for (elem in list_enhancer){
  print(elem)
  res <- try(FetchData(data6_reduced,vars = elem))
  if (inherits(res, "try-error")){list_sc <- c()}
  else{
    expr <- FetchData(data6_reduced,vars = elem)
    res <- try(Cells(data6_reduced[,which(x=expr > 0)]))
    if (inherits(res, "try-error")){list_sc <- c()}
    else{list_sc <- Cells(data6_reduced[,which(x=expr > 0)])}
  }
  list_pcr <- df_all_list_PCR_reduced$liste[df_all_list_PCR_reduced$enhancer==elem]
  nb_cell <- length(intersect(list_sc,unlist(strsplit(df_all_list_PCR_reduced$liste[df_all_list_PCR_reduced$enhancer==elem],","))))
  nb_cell_common_reduced[nrow(nb_cell_common_reduced) +1,] = c(elem,nb_cell)
  double_list <- paste(list_sc,list_pcr,sep=",")
  unique_list <- unique(unlist(strsplit(double_list,",")))
  if (is.null(unique_list)){next}
  else{
    test_subset <- subset(data6_reduced, subset = CellName %in% unique_list)
    print(table(Idents(test_subset)))
  }
}
write.csv(x=nb_cell_common_reduced,file="nb_cell_common_reduced.csv")


