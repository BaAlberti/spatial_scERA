# Packages
library(Matrix)
library(Seurat)
library(ggplot2)

# Loading dataset
# the normalized dataset can be downloaded at : http://bimsbstatic.mdc-berlin.de/rajewsky/DVEX/dge_normalized.txt.gz
data_karaiskos <- read.table("data/dge_normalized.txt",header=TRUE,sep="\t")
data_karaiskos <- CreateSeuratObject(data_karaiskos, project="zinzen")

# Seurat classic pipeline
data_karaiskos <- FindVariableFeatures(object = data_karaiskos,selection.method = "mvp",verbose = FALSE)
data_karaiskos <- ScaleData(object = data_karaiskos)
data_karaiskos <- RunPCA(object = data_karaiskos,features = VariableFeatures(object = data_karaiskos))

#Selecting the number of component
ElbowPlot(data_karaiskos, ndims = 50)
nPC = 16

#Performing clustering and UMAP representation
data_karaiskos <- FindNeighbors(data_karaiskos, dims = 1:nPC)
data_karaiskos <- FindClusters(data_karaiskos, resolution = 0.5)
data_karaiskos <- RunUMAP(data_karaiskos, dims = 1:nPC)

DimPlot(data_karaiskos,reduction="umap",pt.size=2)

# Find differentially expressed genes per cluster
cluster.markers <- FindAllMarkers(data_karaiskos, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Gathering the top10 differentially expressed genes per cluster 
list_10_gene=c()
for (i in 0:5){
  clusters=subset(cluster.markers,cluster.markers$cluster==i)
  top10=head(clusters,10)
  list_10_gene=c(list_10_gene,top10[,7])
  print(top10)
}

# Extracting the list of 60 differentially expressed genes
write.table(list_10_gene, file="top_differentially_expressed_genes_karaiskos.csv",row.names=F,col.names=F)


