library(ggplot2)
library(Seurat)
library(DAIS)

dir <- getwd()
outdir <- paste0(dir,"//Output")
setwd(outdir)
indir <- paste0(dir,"//Input//10X")

############################## 00.Seurat object data readin   ##############################
###if infile is a SeuratObject
#countexp.Seurat <- readRDS(infile)

###if infile is a dataset generated with the Visium technology from 10x Genomics
#countexp.Seurat <- Seurat::Load10X_Spatial(infile)

###if infile is a dataset such as gem.gz generated with the Visium technology from Stereo-seq,please execute the following command
#infile2 <- 'LC00-M_FJ2_bin50.rds'  #exanple
#Transform.Mtx2SeuratObj(infile = infile,bin.size=50,outdir=outdir,outfile=infile2)
#countexp.Seurat <- readRDS(infile2)

###if infile is a dataset generated with the Visium technology from FFPE
exp.mat <- Read10X_h5(filename = paste0(indir,"//GSM5924035_ffpe_c_20_filtered_feature_bc_matrix.h5"))
countexp.Seurat <- Seurat::CreateSeuratObject(counts = exp.mat,
                                              project = 'TLS', assay = 'Spatial',
                                              min.cells=5, min.features=5)
img <- Seurat::Read10X_Image(image.dir = indir)
Seurat::DefaultAssay(object = img) <- 'Spatial'
img <- img[colnames(x = countexp.Seurat)]
countexp.Seurat[['image']] <- img


############################## 01.compute geneset score   ##############################
signaturesFile <- as.character(paste0(indir,"//TLS.signature.v2.gmt"))
method <- "AUCell"
signature_exp <- sc.Meta(object = countexp.Seurat,assays = "Spatial",
                         signaturesFile = signaturesFile,method = method)
write.table(signature_exp, file=paste(method,'score.tsv',sep="."),quote=F,sep="\t")


# --------------part2 plot  real TLS in  the original article-----------
cells <- colnames(countexp.Seurat)
mat1 <- matrix(ncol=2,nrow=length(cells))
mat1[,1] <- cells
mat1[,2] <- 1
colnames(mat1) <- c("Barcode","score")

annofile <- paste0(indir,"//TLS_annotation.mod.csv") #position of real TLS in the original article
anno <- read.csv(annofile,header=T)
mg1 <- merge(mat1,anno,by="Barcode",all=T)
all(mg1$Barcode == colnames(countexp.Seurat))
mg1[is.na(mg1)] <- "NO"
mg1$score <- ifelse(mg1$TLS_2_cat =="NO", 0, ifelse(mg1$TLS_2_cat=="NO_TLS",1,2))  # 0 = not annotated, 1= not TLS, 2=TLS
countexp.Seurat@meta.data$TLS <- mg1$score

#position of real TLS in the original article
realTLS <- subset(countexp.Seurat@meta.data,TLS==2)
p1 <- SpatialPlot(countexp.Seurat,stroke = NA, alpha = 0.8, pt.size.factor = 0.8,  cells.highlight = as.vector(rownames(realTLS)), cols.highlight = c("gold","darkblue"))   
print(p1)



############################## 02.cell select ##############################
# part3 1E
xy <- countexp.Seurat@images$image@coordinates
merg <- cbind(xy,signature_exp)
colnames(merg)[8] <- "Score"
quantile(merg[,8],seq(0.9,1,0.01)) # 99%
i <- round(quantile(merg[,8],0.99),3)

p2 <- ggplot(merg,aes(x=row,y=col)) +
  geom_point(aes(fill = Score),stroke=NA,shape=23, alpha = 0.5, size=2) +
  scale_fill_gradientn(colours = c("#0776B3","#faf6bf","red"))  +
  ggtitle("Geneset=Overlap, method=AUCell") +
  theme_test() + coord_fixed(ratio = 1) #  guides(fill=guide_legend(title="Score"))
print(p2)


f1 <- subset(merg,Score >= i)
p3 <- ggplot(f1,aes(x=row,y=col)) +  geom_point(aes(fill = Score),stroke=NA,shape=23,alpha = 1, size=2.4)+  theme_test() +
  scale_fill_gradientn(colours = c("#faf6bf","red")) + ggtitle(paste("TLS cutoff is ",i,sep=" ")) +
  coord_fixed(ratio = 1)
print(p3)


############################## 03.Density clustering ##############################
type.pos <- f1[,c(2,3)] #spots used to cluster
k1 <- 4
dbscan::kNNdistplot(type.pos, k =  k1)
eps1 <- 6
fx <- Dbscan.cluster(type.pos,eps = eps1,MinPts = k1)
table(fx$cluster)

###Clustering result graph after filter
p4 <- ggplot(fx,aes(x=row,y=col)) +
  geom_point(aes(fill = as.factor(cluster)), shape=23,stroke=NA,alpha = 1,size=3)+
  scale_fill_manual(values =  c("#f8766d")) +
  coord_fixed(ratio = 1) +
  guides(fill=guide_legend(title="Cluster")) +
  theme_test() +
  ggtitle("Cluster for TLS prediction")
print(p4)
write.table(fx,file=paste0("k",k1,"eps",eps1,".cutoff_",i,".tsv",sep="."),quote=F,sep="\t")
stat1 <- as.data.frame(table(fx$Region))


############################## 04.Contour recognition ##############################
fy <- fx
p <- 1
mer1 <- list()
mer2 <- list()

clusterX <- unique(fy$cluster)
for(i in clusterX){
  C1 <- extractCluster(fy, i,0.5)  #Spread radius = 0.5
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
area1

############################## 05.Get the outline of TLS region ##############################
###outline of TLS region
p5 <- ggplot() + geom_sf(data=mer1[[1]],alpha=0.2)  + theme_test()  +
  ggtitle(paste0("TLS area =",area1$merge))
print(p5)

#mat2 <- merg

###outline of TLS region in all cell map 
p6 <- ggplot(merg,aes(x=row,y=col)) +
  geom_point(aes(fill = Score), shape=23,alpha = 1, size=2,stroke=NA)+
  theme_test() +
  scale_fill_gradientn(colours = c("#0776B3","#faf6bf","red")) +
  geom_path(data=mer2[[1]],aes(x=row,y=col,group=cluster)) +
  labs(x="x",y="y") +
  ggtitle("TLS candidate region")  +
  coord_fixed(ratio = 1)
print(p6)

p7 <- ggplot(merg,aes(x=row,y=col)) +
  geom_point(aes(fill = Score), shape=23,alpha = 1, size=2,stroke=NA)+
  theme_test() +
  scale_fill_gradientn(colours = c("#0776B3","#faf6bf","red")) +
  geom_path(data=mer2[[1]],aes(x=row,y=col,group=cluster)) +
  labs(x="x",y="y") +
  ggtitle("TLS candidate region")
print(p7)


############################## 06.Determine if the spot is in TLS region ##############################
dex2 <- matrix(ncol = 5)
mat2 <- merg[,c(1,2,3,8)]
colnames(dex2) <- c(colnames(mat2),"cluster")
for(n in 1:length(mer2)){
  mat3 <- as.data.frame(matrix(unlist(mer2[[n]]),ncol=3))
  colnames(mat3) <- c("row","col","cluster")
  yy <- rangexy(mat3,mat2,n, p="IM" )
  dex2 <- rbind(dex2,yy)
  print(n)
}
dex2 <- dex2[complete.cases(dex2), ]
stat <- as.data.frame(table(dex2$cluster))
stat

write.table(dex2,file=paste0("k",k1,".eps",eps1,".cluster.tsv"),quote=F,sep="\t")

p8 <- ggplot(mat2,aes(x=row,y=col)) +
  geom_point(aes(fill = Score), alpha = 1, size=2,shape=23,stroke=NA)+
  theme_test() +
  scale_fill_gradientn(colours = c("#0776B3","#faf6bf","red")) +
  ggtitle(paste0("k=",k1,",eps=",eps1,sep=" ")) +
  geom_path(data=mer2[[1]],aes(x=row,y=col,group=as.factor(cluster),color=as.factor(cluster)),linewidth=1) +
  scale_linetype_manual(values = 3) + #指定linetype顺序
  guides(color=guide_legend(title="Cluster")) +
  coord_fixed(ratio = 1)
print(p8)


pdf("10X.enrich.pdf",6,6)
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
print(p8)


############################## Fig1G ##############################
library(VennDiagram)
pred1 <- dex2
pred1 <- subset(pred1,cluster !=0)
DAIS <- rownames(pred1)
tls1 <- subset(mg1,score >1)
TLS <- tls1[,1]


# create Venn plot
vennplot <- venn.diagram(
  x = list(DAIS = DAIS, TLS = TLS),
  category.names = c("DAIS predicted spots", "verified TLS"),
  filename = NULL
)
plot.new()
grid.draw(vennplot)
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

