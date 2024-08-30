#extracting cells where four genes are expressed
FeaturePlot(data6_reduced,c("Mef2","Lim1","vnd","Scr"))

spe_gene=c("Mef2","Lim1","vnd","Scr")
Id_cells_with_genes_sc <- data.frame(matrix(ncol=2,nrow=0))
names(Id_cells_with_genes_sc) <- c("gene","sc_cells")
for (elem in spe_gene) {
  expr <- FetchData(data6_reduced,vars = elem)
  list_cell <- paste(Cells(data6_reduced[,which(x=expr > 1.5)]),collapse = ",")
  Id_cells_with_genes_sc[nrow(Id_cells_with_genes_sc) +1,] = c(elem,list_cell)
}
write.csv(x = Id_cells_with_genes_sc,file = "Id_cells_with_four_genes_fig3C.csv")





clusters=c("endoderm")
for (celltype in clusters){
  Cluster_Highlight_Plot(data6_reduced,cluster_name=celltype,highlight_color = "darkblue",background_color = "lightgray",figure_plot = TRUE)
  ggsave(filename= paste(celltype,"_dimplot_reduced.png"),dpi=300,width=9,height=7)
}
