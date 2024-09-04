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

df_all_list_PCR_full <- read.csv("../python_analyses/merged_full/lists_cells_for_each_enhancer_full_version_merged.tsv",sep="\t", header=T)
df_all_list_PCR_reduced <- read.csv("../python_analyses/merged_reduced/lists_cells_for_each_enhancer_reduced_version_merged.tsv",sep="\t",header=T)



DimPlot(data6_full,reduction="umap", group.by = "replicate")
ggsave(filename = "batch_effect_integration.png",dpi=300,width=8,height=7.5)

DimPlot_scCustom(data6_full, reduction = "umap", group.by = 'seurat_clusters', label=F, pt.size = 1.15,figure_plot = TRUE,ggplot_default_colors = T) +NoLegend()
ggsave(filename = "DimPlot_full_without_label.png",dpi=300,width=9,height=7)
DimPlot_scCustom(data6_full, reduction = "umap", label=F, pt.size = 1.15,figure_plot = TRUE,ggplot_default_colors = T) +NoLegend()
ggsave(filename = "DimPlot_full_with_label.png",dpi=300,width=9,height=7)
DimPlot_scCustom(data6_reduced, reduction = "umap", group.by = 'seurat_clusters', label=F, pt.size = 1.15,figure_plot = TRUE,ggplot_default_colors = T) +NoLegend()
ggsave(filename = "DimPlot_reduced_without_label.png",dpi=300,width=9,height=7)
DimPlot_scCustom(data6_reduced, reduction = "umap", label=F, pt.size = 1.15,figure_plot = TRUE,ggplot_default_colors = T) +NoLegend()
ggsave(filename = "DimPlot_reduced_with_label.png",dpi=300,width=9,height=7)

FindVariableFeatures(data6_full, selection.method = "mvp")
data6_full <- ScaleData(data6_full,features = rownames(data6_full))
cluster.markers_6 <- FindAllMarkers(data6_full, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
list_3_gene=c()
for (i in c("ectoderm","unknown","procephalic ectoderm","dorsal ectoderm","mesoderm","amnioserosa","endoderm","mesectoderm","yolk","pole cells","hemocytes")){
  clusters=subset(cluster.markers_6,cluster.markers_6$cluster==i)
  top3=head(clusters,3)
  list_3_gene=c(list_3_gene,top3[,7])
}
Clustered_DotPlot(data6_full,features=list_3_gene,k=11,colors_use_idents=scCustomize_Palette(11,ggplot_default_colors = T),colors_use_exp=c("lightgrey","#9972ec","green","blue"),column_label_size = 18,legend_label_size=16,row_label_size = 18,plot_km_elbow = FALSE) #colors_use_idents=scCustomize_Palette(11,ggplot_default_colors = T)
ggsave(filename="clustered_dotplot_full.png",dpi=300,width=12000,height=10000,units="px")

FindVariableFeatures(data6_reduced, selection.method = "mvp")
data6_reduced <- ScaleData(data6_reduced,features = rownames(data6_reduced))
cluster.markers_6 <- FindAllMarkers(data6_reduced, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
list_3_gene=c()
for (i in c("ectoderm","unknown","procephalic ectoderm","dorsal ectoderm","mesoderm","amnioserosa","endoderm","mesectoderm","yolk","pole cells","hemocytes")){
  clusters=subset(cluster.markers_6,cluster.markers_6$cluster==i)
  top3=head(clusters,3)
  list_3_gene=c(list_3_gene,top3[,7])
}
DoHeatmap(data6_reduced,features=c(list_3_gene,"sog"),size=4)
list_col <- show_col(hue_pal()(9))
Clustered_DotPlot(data6_reduced,features=list_3_gene,k=11,colors_use_idents=scCustomize_Palette(11,ggplot_default_colors = T),colors_use_exp=c("lightgrey","#9972ec","green","blue"),column_label_size = 18,legend_label_size=16,row_label_size = 18,plot_km_elbow = FALSE) #colors_use_idents=scCustomize_Palette(11,ggplot_default_colors = T)
ggsave(filename="clustered_dotplot_full.png",dpi=300,width=12000,height=10000,units="px")

example_genes <- c("Trim9","tup","sog","oc","ken","stumps","sisA","Ptx1","stai","gcm","Tom")
for (elem in example_genes) {
  FeaturePlot_scCustom(data6_full,features = elem ,pt.size = 1.25,figure_plot = T,colors_use =c('lightgrey','blue') ) + theme(aspect.ratio=1)
  ggsave(filename = paste(elem,"FeaturePlot_full.png",sep=""),dpi=300,width = 9,height = 7)
  FeaturePlot_scCustom(data6_reduced,features = elem ,pt.size = 1.25,figure_plot = T,colors_use =c('lightgrey','blue') ) + theme(aspect.ratio=1)
  ggsave(filename = paste(elem,"FeaturePlot_reduced.png",sep=""),dpi=300,width = 9,height = 7)
}

for (elem in list_enhancer) {
  res <- try(FeaturePlot_scCustom(data6_full,features = elem ,colors_use=viridis_plasma_dark_high,pt.size = 1.25,figure_plot = T))
  if (inherits(res, "try-error")){next}
  else{
    FeaturePlot_scCustom(data6_full,features = elem ,colors_use=viridis_plasma_dark_high,pt.size = 1.25,figure_plot = T)
    ggsave(filename = paste(elem,"FeaturePlot_full.png",sep=""),dpi=300,width = 9,height = 7)
    VlnPlot(data6_full,features = elem)
    ggsave(filename= paste(elem,"VlnPlot_full.png",sep=""),dpi=300)
  }
}
for (elem in list_enhancer) {
  res <- try(FeaturePlot_scCustom(data6_reduced,features = elem ,colors_use=viridis_plasma_dark_high,pt.size = 1.25,figure_plot = T))
  if (inherits(res, "try-error")){next}
  else{
    FeaturePlot_scCustom(data6_reduced,features = elem ,colors_use=viridis_plasma_dark_high,pt.size = 1.25,figure_plot = T)
    ggsave(filename = paste(elem,"FeaturePlot_reduced.png",sep=""),dpi=300,width = 9,height = 7)
    VlnPlot(data6_reduced,features = elem)
    ggsave(filename= paste(elem,"VlnPlot_reduced.png",sep=""),dpi=300)
  }
}

for (elem in list_enhancer){
  res <- try(FetchData(data6_full,vars = elem))
  if (inherits(res, "try-error")){next}
  else{
    expr <- FetchData(data6_full,vars = elem)
    res2 <- try(Cells(data6_full[,which(x=expr > 0)]))
    if (inherits(res2, "try-error")){next}
    else{
      list_sc <- Cells(data6_full[,which(x=expr > 0)])
      Cell_Highlight_Plot(data6_full,cells_highlight = list(elem=list_sc),figure_plot=T,
                          highlight_color = "red",background_color = "lightblue") + theme(legend.position = "none")
      ggsave(filename= paste(elem,"_sc_dimplot_full.png"),dpi=300,width=9,height=7)
    }
  }
}

for (elem in list_enhancer){
  res <- try(FetchData(data6_reduced,vars = elem))
  if (inherits(res, "try-error")){next}
  else{
    expr <- FetchData(data6_reduced,vars = elem)
    res2 <- try(Cells(data6_reduced[,which(x=expr > 0)]))
    if (inherits(res2, "try-error")){next}
    else{
      list_sc <- Cells(data6_reduced[,which(x=expr > 0)])
      Cell_Highlight_Plot(data6_reduced,cells_highlight = list(elem=list_sc),figure_plot=T,
                          highlight_color = "red",background_color = "lightblue") + theme(legend.position = "none")
      ggsave(filename= paste(elem,"_sc_dimplot_reduced.png"),dpi=300,width=9,height=7)
    }
  }
}

# Plotting where the cells for each enhancer for PCR only are in the UMAP
for (elem in df_all_list_PCR_full$enhancer){
  Cell_Highlight_Plot(data6_full,figure_plot = T,cells_highlight = list("PCR only"=unlist(strsplit(df_all_list_PCR_full$liste[df_all_list_PCR_full$enhancer==elem],","))),
                      highlight_color = "red",background_color = "lightblue")
  ggsave(filename= paste(elem,"_PCR_dimplot_full.png"),dpi=300,width=9,height=7)
}

# Plotting where the cells for each enhancer for PCR only are in the UMAP
for (elem in df_all_list_PCR_reduced$enhancer){
  Cell_Highlight_Plot(data6_reduced,figure_plot = T,cells_highlight = list("PCR only"=unlist(strsplit(df_all_list_PCR_reduced$liste[df_all_list_PCR_reduced$enhancer==elem],","))),
                      highlight_color = "red",background_color = "lightblue")
  ggsave(filename= paste(elem,"_PCR_dimplot_reduced.png"),dpi=300,width=9,height=7)

# Plotting where the cells for each enhancer for sc and PCR together are in the UMAP
for (elem in list_enhancer){
  print(elem)
  res <- try(FetchData(data6_full,vars = elem))
  if (inherits(res, "try-error")){list_sc <- c()}
  else{
    expr <- FetchData(data6_full,vars = elem)
    res2 <- try(Cells(data6_full[,which(x=expr > 0)]))
    if (inherits(res2, "try-error")){list_sc <- c()}
    else{list_sc <- Cells(data6_full[,which(x=expr > 0)])}
  }
  list_pcr <- df_all_list_PCR_full$liste[df_all_list_PCR_full$enhancer==elem]
  double_list <- paste(list_sc,list_pcr,sep=",")
  unique_list <- unique(unlist(strsplit(double_list,",")))
  Cell_Highlight_Plot(data6_full,figure_plot = TRUE,cells_highlight = list(elem=unique_list),highlight_color = "red",background_color = "lightblue")+NoLegend() +theme(legend.position = "none")
  ggsave(filename= paste(elem,"_merged_sc_PCR_dimplot_full.png"),dpi=300,width=9,height=7)
}  

for (elem in list_enhancer){
  print(elem)
  res <- try(FetchData(data6_reduced,vars = elem))
  if (inherits(res, "try-error")){list_sc <- c()}
  else{
    expr <- FetchData(data6_reduced,vars = elem)
    res2 <- try(Cells(data6_reduced[,which(x=expr > 0)]))
    if (inherits(res2, "try-error")){list_sc <- c()}
    else{list_sc <- Cells(data6_reduced[,which(x=expr > 0)])}
  }
  list_pcr <- df_all_list_PCR_reduced$liste[df_all_list_PCR_reduced$enhancer==elem]
  double_list <- paste(list_sc,list_pcr,sep=",")
  unique_list <- unique(unlist(strsplit(double_list,",")))
  Cell_Highlight_Plot(data6_reduced,figure_plot = TRUE,cells_highlight = list(elem=unique_list),highlight_color = "red",background_color = "lightblue")+NoLegend() +theme(legend.position = "none")
  ggsave(filename= paste(elem,"_merged_sc_PCR_dimplot_reduced.png"),dpi=300,width=9,height=7)
}  

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
  list_pcr <- unlist(strsplit(list_pcr,","))
  both <- intersect(list_sc,list_pcr)
  unique_sc <- setdiff(list_sc,both)
  unique_pcr <- setdiff(list_pcr,both)
  cells <- list("SC only" = unique_sc,
                "PCR only"=unique_pcr,
                "common"=both)
  Cell_Highlight_Plot(seurat_object = data6_full, cells_highlight = cells,highlight_color = c("#4ed419","#3B4CC0","#e8a70c"),background_color = "lightblue",pt.size=1.5,figure_plot=T) + NoLegend() + theme(legend.position="none")
  ggsave(filename= paste(elem,"_addition_PCR_dimplot_full.png"),dpi=300,width=9,height=7)
}

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
  list_pcr <- unlist(strsplit(list_pcr,","))
  both <- intersect(list_sc,list_pcr)
  unique_sc <- setdiff(list_sc,both)
  unique_pcr <- setdiff(list_pcr,both)
  cells <- list("SC only" = unique_sc,
                "PCR only"=unique_pcr,
                "common"=both)
  Cell_Highlight_Plot(seurat_object = data6_reduced, cells_highlight = cells,highlight_color = c("#4ed419","#3B4CC0","#e8a70c"),background_color = "lightblue",pt.size=1.5,figure_plot=T) + NoLegend() + theme(legend.position="none")
  ggsave(filename= paste(elem,"_addition_PCR_dimplot_reduced.png"),dpi=300,width=9,height=7)
}

for (celltype in levels(data6_full@active.ident)){
  Cluster_Highlight_Plot(data6_full,cluster_name=celltype,highlight_color = "darkblue",background_color = "lightgray",figure_plot = TRUE) +NoLegend()
  ggsave(filename= paste(celltype,"_dimplot_full.png"),dpi=300,width=9,height=7)
}

for (celltype in levels(data6_reduced@active.ident)){
  Cluster_Highlight_Plot(data6_reduced,cluster_name=celltype,highlight_color = "darkblue",background_color = "lightgray",figure_plot = TRUE)
  ggsave(filename= paste(celltype,"_dimplot_reduced.png"),dpi=300,width=9,height=7)
}

#extracting cells where four genes are expressed
spe_gene=c("Mef2","Lim1","vnd","Scr")
Id_cells_with_genes_sc <- data.frame(matrix(ncol=2,nrow=0))
names(Id_cells_with_genes_sc) <- c("gene","sc_cells")
for (elem in spe_gene) {
  expr <- FetchData(data6_reduced,vars = elem)
  list_cell <- paste(Cells(data6_reduced[,which(x=expr > 1.5)]),collapse = ",")
  Id_cells_with_genes_sc[nrow(Id_cells_with_genes_sc) +1,] = c(elem,list_cell)
}
write.csv(x = Id_cells_with_genes_sc,file = "Id_cells_with_four_genes_fig3C.csv")