#####  method2:seurat   #####
#Rscript $0 outdir inrds pro pc res
#pc=15,res=1

#Transfer parameters
args <- commandArgs(T)
outdir <- args[1]
inrds <- args[2] 
pro <- args[3]
pc <- as.numeric(args[4])
res <- as.numeric(args[5])

set.seed(6)
library(ggplot2)
library(Seurat)

pal <- c("#1f77b4","#ff7f0e","#279e68","#d62728","#aa40fc","#8c564b",
         "#e377c2","#b5bd61","#17becf","#aec7e8","#ffbb78","#98df8a",
         "#ff9896","#c5b0d5","#c49c94","#FFFF00","#1CE6FF","#FF34FF",
         "#FC0D05","#008941","#006FA6","#A30059","#9400D3","#4CC700B6",
         "#3B5DFF","#FF2F80","#FF7C75FA","#00FFCC","#612E21","#FF4A46",
         "#0000A6","#63FFAC","#B79762","#004D43","#8FB0FF","#997D87",
         "#5A0007","#809693","#1B4400","#4FC601")

if (!dir.exists(outdir)) {
  dir.create(outdir)
}
setwd(outdir)

obj <- readRDS(inrds)
obj <- SCTransform(obj, assay = "Spatial", ncells = 3000, verbose = FALSE)
#obj <- NormalizeData(obj)
#obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
#obj <- ScaleData(obj, features = rownames(obj))
obj <- RunPCA(obj,assay = "SCT",features = VariableFeatures(object = obj))
p <- ElbowPlot(obj, ndims = 50)
ggsave(paste0(pro,"_Elbow.png"),plot = p,width = 8,height = 8)
obj <- RunUMAP(obj,reduction = "pca", dims = 1:pc)
obj <- FindNeighbors(obj,reduction = "pca", dims = 1:pc)
obj <- FindClusters(obj,resolution = res, verbose = FALSE)
saveRDS(obj,file = paste0(pro,"_dim",pc,"res",res,"_Cluster.rds"))
write.table(obj@meta.data,file=paste0(pro,"_dim",pc,"res",res,"_Cluster.meta.txt"),sep="\t",quote=F)

#plot
col <- pal[1:length(unique(obj$seurat_clusters))]
names(col) <- unique(obj@meta.data$seurat_cluster)
plot1 <- DimPlot(obj, reduction = "umap", label = TRUE,cols = col)
plot2 <- SpatialDimPlot(obj, stroke = NA,
                        cols = col)+
  ggtitle(pro)
pdf(paste0(pro,"_Cluster.pdf"),width = 8,height = 8)
print(plot1 + plot2)
dev.off()
plot3 <- SpatialDimPlot(obj, stroke = NA,
                        cols = col)+
  NoLegend()
ggsave(paste0(pro,"_Cluster.png"),plot = plot3,width = 5,height = 5)





