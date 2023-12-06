library(ggplot2)
library(Seurat)
library(DAIS)
dir <- getwd()
setwd(paste0(dir,'//Output'))


############################## 00.Input data ##############################
mat <- read.table(paste0(dir,"//Input//Stereo-seq//LC00-M_FJ2_bin50_CelltypeTrans_MaligComb_Spotlight_First.txt.gz"),header=T)
mat$predict_CellType <- gsub("Cholangiocyte|Malignant","Cholang.Malig",mat$predict_CellType)
mat$predict_CellType <- gsub("T.cell|NK","T.NK",mat$predict_CellType)
cellType <- c("B.cell","Plasma","T.NK","Macrophage","DC","Cholang.Malig","Hepatocyte","Endothelial","Fibroblast")
mat$predict_CellType <- factor(mat$predict_CellType,levels = cellType)
colnames(mat) <- c('bin_name','row','col','celltype')

# Cell annotation
cls <- c("B.cell"="#808ad0","Plasma"="#0e61ba","T.NK"="#dcbeff","Macrophage"="pink",
         "DC"="#f58231","Cholang.Malig"="#C0C0C0","Hepatocyte"="#f9dca7","Endothelial"="#97d1b4",
         "Fibroblast"="#0ea387")

p0 <- ggplot(mat,aes(x=row,y=col)) +
  geom_tile(aes(fill=celltype)) +
  theme_test() +
  coord_fixed(ratio = 1)  +
  scale_fill_manual(values=cls)  +
  guides(fill=guide_legend(title="Cell type"))
print(p0)


###----- 1.Fig1C- stroma region -----
############################## 01.Select cells ##############################
pre <- "LC0-B_FJ2.stroma"
max_immnum <- 99999
max_area <- 99999
celltype_select <- c('Fibroblast','Endothelial')

metaall <- mat
meta <- subset(metaall,celltype %in% celltype_select)

p1 <- ggplot(meta,aes(x=row,y=col)) +
  geom_tile(aes(fill = celltype),alpha = 1)+
  theme_test()  +
  coord_fixed(ratio = 1)   +
  guides(fill=guide_legend(title="Cell type"))
print(p1)

############################## 02.Density clustering ##############################
type.pos <- meta[,c('row','col')]  #row,col
k1 <- 15
eps1 <- 3
min_immnum <- 50
fx <- Dbscan.cluster(type.pos = type.pos, eps = eps1, MinPts = k1)
fx <- deletelines(fx,min_immnum)  # result of clustering
num_TLS <- length(unique(fx$cluster))
num_TLS
table(fx$cluster)

#Clustering result
p2 <- ggplot(fx,aes(x=row,y=col)) +
  geom_tile(aes(fill = as.factor(cluster)),alpha = 1)+
  theme_test() +
  ggtitle(paste0("k=",k1,",eps=",eps1," min=",min_immnum)) +
  coord_fixed(ratio = 1)  +
  guides(fill=guide_legend(title="Cluster"))
print(p2)

write.table(fx,file=paste0(pre,".","k",k1,"eps",eps1,".tsv",sep="."),quote=F,sep="\t")
stat1 <- as.data.frame(table(fx$cluster))
stat1

############################## 03.Contour recognition ##############################
p <- 1
mer1 <- list()
mer2 <- list()
colnames(fx)[3] <- "Region"
clusterX <- unique(fx$Region)
alpha <- 6
for(i in clusterX){
  C1 <- extractCluster2(fx,i,0.5,alpha = alpha)  #r=0.5,select a-shape method
  mer1[[p]] <- C1[[1]]  ## polygen list
  mer2[[p]] <- C1[[2]]  ## dataframe list
  p <- p +1
}

############################## 04.Calculate Area ##############################
area1 <- as.data.frame(matrix(nrow=length(clusterX),ncol=2,dimnames =list(clusterX,c("cluster","area"))))
for (i in 1:length(mer1)){
  area1[i,2] <- st_area(mer1[[i]])
  area1[i,1] <- clusterX[i]
}
area1$merge <- paste(area1[,1],round(area1[,2],1),sep=":")

dataframes_with_label <- lapply(seq_along(mer2), function(i) {
  df <- mer2[[i]]
  df$source <- as.factor(paste0("Dataframe ", i))
  return(df)
})
combined_data <- bind_rows(dataframes_with_label)

cls2 <- c("B.cell"="#FAFAD2","T.NK"="#FAFAD2", "Plasma"="#FAFAD2","DC"="#FFE4E1","Macrophage"="#FFE4E1","Endothelial"="#AFEEEE","Fibroblast"="#AFEEEE",
         "Cholang.Malig"="#808080","Hepatocyte"="#DCDCDC")
p3 <- ggplot(metaall,aes(x=row,y=col)) +
  geom_tile(aes(fill = celltype),alpha = 1)+
  theme_test() +
  guides(fill=guide_legend(title="Cell type")) +
  geom_path(data=combined_data,aes(x=x,y=y,group=cluster),color="black") +
  scale_fill_manual(values = cls2) +
  labs(x="x",y="y") + ggtitle(paste(pre,",k=",k1,",eps=",eps1,",alpha=",alpha,sep=" "))  +
  coord_fixed(ratio = 1)
print(p3)

############################## 05.Determine if the spot is in the region ##############################
dex <- matrix(ncol = 5)
colnames(dex) <- c(colnames(metaall),"cluster")
for(n in 1:length(mer2)){
  mat3 <- as.data.frame(matrix(unlist(mer2[[n]]),ncol=3))
  colnames(mat3) <- c("row","col","cluster") # 更新 x- row
  yy <- rangexy(mat3,metaall,n, p="IM" ) # p没用
  dex <- rbind(dex,yy)
  print(n)
}
dex <- dex[complete.cases(dex), ]

#number of spots in every region
stat <- as.data.frame(table(dex$cluster))
stat

#plot the outlines of different clusters with different colors
p4 <- ggplot(metaall,aes(x=row,y=col)) +  geom_tile(aes(fill = celltype),alpha = 1) +
  theme_test() +
  geom_path(data=combined_data,aes(x=x,y=y,color=as.factor(cluster)),linewidth=1.5) +
  scale_fill_manual(values = cls) +
  labs(x="x",y="y") +  coord_fixed(ratio = 1) +
  ggtitle(paste(pre,",k=",k1,",eps=",eps1,",alpha=",alpha,sep=" "))  +
  scale_linetype_manual(values = 3)  + guides(color=guide_legend(title="Cluster"))
print(p4)

write.table(metaall,file=paste0(pre,".k",k1,".eps",eps1,".alpha",alpha,".cluster.tsv"),quote=F,sep="\t")
write.table(combined_data,file=paste0(pre,".k",k1,".eps",eps1,".alpha",alpha,".lines.tsv"),quote=F,sep="\t")
write.table(area1,file=paste0(pre,".k",k1,".eps",eps1,".alpha",alpha,".area.tsv"),quote=F,sep="\t")



###----- 2.Fig1C- tumor region -----
############################## 01.Select cells ##############################
pre <- "LC0-B_FJ2.tumor"
celltype_select <- c('Cholang.Malig')
metaall <- mat
meta <- subset(metaall,celltype %in% celltype_select)

p5 <- ggplot(meta,aes(x=row,y=col)) +
  geom_tile(aes(fill = celltype),alpha = 1)+
  theme_test()  +
  coord_fixed(ratio = 1)  +
  guides(fill=guide_legend(title="Cell type"))
print(p5)

type.pos <- meta[,c('row','col')]  #row,col
k1 <- 21
eps1 <- 3
min_immnum <- 100
fx <- Dbscan.cluster(type.pos = type.pos, eps = eps1, MinPts = k1)
length(unique(fx$cluster))

############################## 02.Density clustering ##############################
fx <- deletelines(fx,min_immnum)  # fy记录第一次聚类结果：row,col,clusterid
num_TLS <- length(unique(fx$cluster))
num_TLS
table(fx$cluster)

p6 <- ggplot(fx,aes(x=row,y=col)) +
  geom_tile(aes(fill = as.factor(cluster)),alpha = 1)+
  theme_test() +
  ggtitle(paste0("k=",k1,",eps=",eps1," min=",min_immnum)) +
  coord_fixed(ratio = 1)  +
  guides(fill=guide_legend(title="Cluster"))
print(p6)

dim(fx)
fx <- unique(fx)
write.table(fx,file=paste0(pre,".","k",k1,"eps",eps1,".tsv",sep="."),quote=F,sep="\t")
stat1 <- as.data.frame(table(fx$cluster))
stat1

############################## 03.Contour recognition ##############################
p <- 1
mer1 <- list()
mer2 <- list()
colnames(fx)[3] <- "Region"
clusterX <- unique(fx$Region)
alpha <- 8

for(i in clusterX){
  C1 <- extractCluster2(fx,i,0.5,alpha = alpha)
  mer1[[p]] <- C1[[1]]  ## polygen list
  mer2[[p]] <- C1[[2]]  ## dataframe list
  p <- p +1
}

############################## 04.Calculate Area ##############################
area1 <- as.data.frame(matrix(nrow=length(clusterX),ncol=2,dimnames =list(clusterX,c("cluster","area"))))
for (i in 1:length(mer1)){
  area1[i,2] <- st_area(mer1[[i]])
  area1[i,1] <- clusterX[i]
}
area1$merge <- paste(area1[,1],round(area1[,2],1),sep=":")

dataframes_with_label <- lapply(seq_along(mer2), function(i) {
  df <- mer2[[i]]
  df$source <- as.factor(paste0("Dataframe ", i))
  return(df)
})
combined_data <- bind_rows(dataframes_with_label)

p7 <- ggplot(metaall,aes(x=row,y=col)) +
  geom_tile(aes(fill = celltype),alpha = 1)+
  theme_test() +
  geom_path(data=combined_data,aes(x=x,y=y,group=cluster),color="black") +
  scale_fill_manual(values = cls) +
  labs(x="x",y="y") +
  ggtitle(paste(pre,",k=",k1,",eps=",eps1,",alpha=",alpha,sep=" "))  +
  coord_fixed(ratio = 1) +
  guides(fill=guide_legend(title="Cell type"))
print(p7)

############################## 05.Determine if the spot is in TLS region ##############################
dex2 <- matrix(ncol = 5)
colnames(dex2) <- c(colnames(metaall),"cluster")
for(n in 1:length(mer2)){
  mat3 <- as.data.frame(matrix(unlist(mer2[[n]]),ncol=3))
  colnames(mat3) <- c("row","col","cluster")
  yy <- rangexy(mat3,metaall,n, p="IM" )
  dex2 <- rbind(dex2,yy)
  print(n)
}
dex2 <- dex2[complete.cases(dex2), ]

stat <- as.data.frame(table(dex2$cluster))
stat

p8 <- ggplot(metaall,aes(x=row,y=col)) +
  geom_tile(aes(fill = celltype),alpha = 1) +
  theme_test() +
  geom_path(data=combined_data,aes(x=x,y=y,color=as.factor(cluster)),linewidth=1.5) +
  scale_fill_manual(values = cls) +
  labs(x="x",y="y") +  coord_fixed(ratio = 1) +
  ggtitle(paste(pre,",k=",k1,",eps=",eps1,",alpha=",alpha,sep=" "))  +
  scale_linetype_manual(values = 3)  +
  guides(color=guide_legend(title="Cluster"))
print(p8)

write.table(dex2,file=paste0(pre,".k",k1,".eps",eps1,".alpha",alpha,".cluster.tsv"),quote=F,sep="\t")
write.table(combined_data,file=paste0(pre,".k",k1,".eps",eps1,".alpha",alpha,".lines.tsv"),quote=F,sep="\t")
write.table(area1,file=paste0(pre,".k",k1,".eps",eps1,".alpha",alpha,".area.tsv"),quote=F,sep="\t")


# line chart
path.strom <- read.table("LC0-B_FJ2.stroma.k15.eps3.alpha6.lines.tsv",sep="\t")
path.tumor <- read.table("LC0-B_FJ2.tumor.k21.eps3.alpha8.lines.tsv",sep="\t")
path1 <- subset(path.tumor, cluster==1)
path6 <- subset(path.tumor, cluster==6)

p9 <- ggplot(mat,aes(x=row,y=col)) +
  geom_tile(aes(fill = celltype),alpha = 1)+  theme_test() +
  geom_path(data=path.strom,aes(x=x,y=y,group=cluster),color="black") +
  scale_fill_manual(values = cls2) +
  geom_path(data=path1,aes(x=x,y=y),color="#e6194B") +
  geom_path(data=path6,aes(x=x,y=y),color="#e6194B") +
  labs(x="x",y="y") + ggtitle("LC00-M_raw")  +
  coord_fixed(ratio = 1)
print(p9)


pdf(file="Stereo-seq.Border.pdf",9,9)
print(p0)
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
print(p8)
print(p9)
dev.off()



sessionInfo()
# R version 4.2.2 (2022-10-31 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)
#
# Matrix products: default

# locale:
#   [1] LC_COLLATE=Chinese (Simplified)_China.utf8
# [2] LC_CTYPE=Chinese (Simplified)_China.utf8
# [3] LC_MONETARY=Chinese (Simplified)_China.utf8
# [4] LC_NUMERIC=C
# [5] LC_TIME=Chinese (Simplified)_China.utf8
#
# attached base packages:
#   [1] stats4    grid      stats     graphics  grDevices utils     datasets  methods
# [9] base
#
# other attached packages:
#   [1] DAIS_1.0.1            VAM_1.0.0             MASS_7.3-58.2
# [4] UCell_2.2.0           tidyr_1.3.0           singscore_1.18.0
# [7] sf_1.0-13             scMetabolism_0.2.1    rsvd_1.0.5
# [10] rjson_0.2.21          RColorBrewer_1.1-3    pheatmap_1.0.12
# [13] patchwork_1.1.2       pagoda2_1.0.10        igraph_1.4.1
# [16] Matrix_1.5-3          magrittr_2.0.3        GSVA_1.46.0
# [19] GSEABase_1.60.0       graph_1.76.0          annotate_1.76.0
# [22] XML_3.99-0.14         AnnotationDbi_1.60.2  IRanges_2.32.0
# [25] S4Vectors_0.36.2      Biobase_2.58.0        BiocGenerics_0.44.0
# [28] ggsci_3.0.0           ggpubr_0.6.0          fpc_2.2-10
# [31] factoextra_1.0.7.999  dplyr_1.1.0           dbscan_1.1-11
# [34] data.table_1.14.8     cowplot_1.1.1         ComplexHeatmap_2.14.0
# [37] clusterProfiler_4.6.2 AUCell_1.20.2         alphahull_2.5
# [40] ggplot2_3.4.2         SeuratObject_4.1.3    Seurat_4.3.0.1
#
# loaded via a namespace (and not attached):
#   [1] ica_1.0-3                   class_7.3-22
# [3] ps_1.7.2                    foreach_1.5.2
# [5] lmtest_0.9-40               crayon_1.5.2
# [7] rhdf5filters_1.10.1         nlme_3.1-162
# [9] backports_1.4.1             GOSemSim_2.24.0
# [11] rlang_1.1.0.9000            XVector_0.38.0
# [13] HDO.db_0.99.1               ROCR_1.0-11
# [15] irlba_2.3.5.1               callr_3.7.3
# [17] limma_3.54.2                BiocParallel_1.32.6
# [19] bit64_4.0.5                 glue_1.6.2
# [21] splancs_2.01-44             sctransform_0.3.5
# [23] parallel_4.2.2              processx_3.8.0
# [25] spatstat.sparse_3.0-2       classInt_0.4-9
# [27] DOSE_3.24.2                 spatstat.geom_3.2-1
# [29] tidyselect_1.2.0            SummarizedExperiment_1.28.0
# [31] usethis_2.2.1               fitdistrplus_1.1-11
# [33] zoo_1.8-12                  xtable_1.8-4
# [35] cli_3.6.0                   zlibbioc_1.44.0
# [37] rstudioapi_0.14             miniUI_0.1.1.1
# [39] sp_2.0-0                    fastmatch_1.1-3
# [41] treeio_1.22.0               shiny_1.7.4.1
# [43] BiocSingular_1.14.0         clue_0.3-64
# [45] pkgbuild_1.4.2              gson_0.1.0
# [47] cluster_2.1.4               tidygraph_1.2.3
# [49] sgeostat_1.0-27             KEGGREST_1.38.0
# [51] tibble_3.2.0                ggrepel_0.9.3
# [53] ape_5.7-1                   listenv_0.9.0
# [55] Biostrings_2.66.0           png_0.1-8
# [57] future_1.33.0               withr_2.5.0
# [59] bitops_1.0-7                ggforce_0.4.1
# [61] plyr_1.8.8                  e1071_1.7-13
# [63] pillar_1.9.0                GlobalOptions_0.1.2
# [65] cachem_1.0.7                fs_1.6.1
# [67] flexmix_2.3-19              kernlab_0.9-32
# [69] GetoptLong_1.0.5            DelayedMatrixStats_1.20.0
# [71] vctrs_0.6.2                 ellipsis_0.3.2
# [73] generics_0.1.3              urltools_1.7.3
# [75] devtools_2.4.5              tools_4.2.2
# [77] munsell_0.5.0               tweenr_2.0.2
# [79] proxy_0.4-27                fgsea_1.24.0
# [81] DelayedArray_0.23.2         fastmap_1.1.1
# [83] compiler_4.2.2              pkgload_1.3.2
# [85] abind_1.4-5                 httpuv_1.6.9
# [87] sessioninfo_1.2.2           plotly_4.10.2
# [89] GenomeInfoDbData_1.2.9      gridExtra_2.3
# [91] edgeR_3.40.2                lattice_0.20-45
# [93] deldir_1.0-9                utf8_1.2.3
# [95] later_1.3.0                 jsonlite_1.8.4
# [97] scales_1.2.1                ScaledMatrix_1.6.0
# [99] dendsort_0.3.4              tidytree_0.4.2
# [101] pbapply_1.7-2               carData_3.0-5
# [103] sparseMatrixStats_1.10.0    lazyeval_0.2.2
# [105] promises_1.2.0.1            car_3.1-2
# [107] doParallel_1.0.17           R.utils_2.12.2
# [109] goftest_1.2-3               spatstat.utils_3.0-3
# [111] reticulate_1.30             brew_1.0-8
# [113] Rtsne_0.16                  downloader_0.4
# [115] uwot_0.1.14                 HDF5Array_1.26.0
# [117] Rook_1.2                    survival_3.5-5
# [119] prabclus_2.3-2              htmltools_0.5.4
# [121] memoise_2.0.1               profvis_0.3.7
# [123] modeltools_0.2-23           locfit_1.5-9.8
# [125] graphlayouts_0.8.4          viridisLite_0.4.2
# [127] digest_0.6.31               mime_0.12
# [129] N2R_1.0.1                   units_0.8-2
# [131] RSQLite_2.3.1               yulab.utils_0.0.6
# [133] future.apply_1.11.0         remotes_2.4.2
# [135] urlchecker_1.0.1            blob_1.2.4
# [137] R.oo_1.25.0                 drat_0.2.3
# [139] labeling_0.4.2              splines_4.2.2
# [141] Rhdf5lib_1.20.0             RCurl_1.98-1.12
# [143] broom_1.0.5                 rhdf5_2.42.1
# [145] colorspace_2.1-0            mnormt_2.1.1
# [147] GenomicRanges_1.50.2        shape_1.4.6
# [149] aplot_0.1.10                nnet_7.3-19
# [151] Rcpp_1.0.10                 mclust_6.0.0
# [153] RANN_2.6.1                  circlize_0.4.15
# [155] enrichplot_1.18.4           fansi_1.0.4
# [157] parallelly_1.36.0           R6_2.5.1
# [159] ggridges_0.5.4              lifecycle_1.0.3
# [161] ggsignif_0.6.4              leiden_0.4.3
# [163] robustbase_0.99-0           qvalue_2.30.0
# [165] RcppAnnoy_0.0.20            iterators_1.0.14
# [167] spatstat.explore_3.2-1      stringr_1.5.0
# [169] htmlwidgets_1.6.2           beachmat_2.14.2
# [171] polyclip_1.10-4             triebeard_0.4.1
# [173] purrr_1.0.1                 RMTstat_0.3.1
# [175] shadowtext_0.1.2            gridGraphics_0.5-1
# [177] mgcv_1.8-42                 globals_0.16.2
# [179] spatstat.random_3.1-5       progressr_0.13.0
# [181] codetools_0.2-19            matrixStats_1.0.0
# [183] GO.db_3.16.0                prettyunits_1.1.1
# [185] psych_2.3.9                 SingleCellExperiment_1.20.1
# [187] R.methodsS3_1.8.2           GenomeInfoDb_1.34.9
# [189] gtable_0.3.3                DBI_1.1.3
# [191] ggfun_0.1.1                 tensor_1.5
# [193] httr_1.4.6                  KernSmooth_2.23-21
# [195] stringi_1.7.12              reshape2_1.4.4
# [197] farver_2.1.1                diptest_0.76-0
# [199] viridis_0.6.3               ggtree_3.6.2
# [201] sccore_1.0.3                BiocNeighbors_1.16.0
# [203] interp_1.1-4                ggplotify_0.1.1
# [205] scattermore_1.2             DEoptimR_1.1-2
# [207] bit_4.0.5                   scatterpie_0.2.1
# [209] MatrixGenerics_1.10.0       spatstat.data_3.0-1
# [211] ggraph_2.1.0                pkgconfig_2.0.3
# [213] rstatix_0.7.2
