#####  method3:BayesSpace   #####
#Rscript $0 outdir inrds pro n_cluster

args <- commandArgs(T)
outdir <- args[1]
inrds <- args[2] 
pro <- args[3]
n <- as.numeric(args[4])  #number of cluster

library(SingleCellExperiment)
library(ggplot2)
library(DBI)
library(BayesSpace)
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

# create sce obj
counts <-obj@assays$Spatial@counts
colData = coord <- obj@images$image@coordinates
colData$row = colData$imagerow <- -coord$col
colData$col = colData$imagecol <- coord$row
rowData <-data.frame(gene_name = rownames(obj@assays$Spatial))
rownames(rowData)<- rowData$gene_name
sce <- SingleCellExperiment(assays=list(counts=as(counts, "dgCMatrix")),
                            rowData=rowData,
                            colData=colData)
sce
# log+PCA
set.seed(2000)
sce <- spatialPreprocess(sce, platform="Visium",  
                         n.PCs=15, n.HVGs=2000, log.normalize=TRUE)
# model:Selecting the number of clusters
sce <- qTune(sce, qs=seq(2,20), platform="Visium", d=15) 
pdf(file=paste0(pro,".bys.qplot.pdf"), 8,8)
qPlot(sce)
dev.off()

# spatial cluster with Bayesspace
set.seed(2000)    
sce <- spatialCluster(sce, q=n, platform="Visium", d=15,
                      init.method="mclust", model="t", gamma=2, 
                      nrep=10000, burn.in=100, 
                      save.chain=TRUE)
sce
head(colData(sce))
saveRDS(sce,paste0(pro,"_BayesSpace.rds"))

# visualize
col <- pal[1:n]
pdf(file=paste0(pro,".bayes.cluster.pdf"), 8,8)
clusterPlot(sce,label = "spatial.cluster",platform = "Visium",palette = col,color = NA) +
  labs(fill="BayesSpace\ncluster",title = pro)
dev.off()
p <- clusterPlot(sce,label = "spatial.cluster",platform = "Visium",palette = col,color = NA)+
  NoLegend()
ggsave(paste0(pro,"_BayesSpace.png"),plot = p,width = 5,height = 5)



