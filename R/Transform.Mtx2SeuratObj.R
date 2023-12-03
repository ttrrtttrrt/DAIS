#' @title Create a Spatial Seurat object from file generated with the Visium technology such as Stereo-seq
#'
#' @param infile file contains at least 4 columns such as the coordinates of spots(named x,y),geneID and UMICount
#' @param bin.size number of bins a spot contains
#' @param outdir output directory
#' @param outfile name of output file
#' @param ncores number of cores used
#'
#' @return a Spatial Seurat object
#' @export
#'
#' @examples
Transform.Mtx2SeuratObj <- function (infile,bin.size=NULL,outdir,outfile,ncores = NULL) {
  if (is.null(ncores)) {
    ncores=5
  }
  if (is.null(ncores)) {
    bin.size=50
  }

  library(Seurat)
  library(clusterProfiler)
  library(dplyr)
  library(data.table)
  library(Matrix)
  library(rjson)
  library(ggplot2)
  library(ggsci)
  library(patchwork)
  library(RColorBrewer)
  library(pheatmap)
  library(cowplot)
  library(ggpubr)


  bs <- bin.size
  if (!dir.exists(outdir)) {
    dir.create(outdir)
  }

  pro <- tail(unlist(strsplit(infile,"/")),1)
  pro <- gsub(".txt|.tsv|.gem|.gz|.lasso.gem.gz|_lasso.gem.gz","",pro)

  ############################## 1. bin data  ##############################
  dat <- fread(file = infile)
  if(length(grep("MIDCounts|MIDCount",colnames(dat))>0)){
    colnames(dat) <- gsub("MIDCounts|MIDCount","UMICount",colnames(dat))}
  out <- as.data.frame(dat)

  dat$x <- trunc((dat$x - min(dat$x)) / bs + 1)
  dat$y <- trunc((dat$y - min(dat$y)) / bs + 1)

  out <- cbind(dat$y,dat$x,out)
  colnames(out)[1:2] <- c(paste0("bin",bs,".y"),paste0("bin",bs,".x"))


  dat <- dat[, sum(UMICount), by = .(geneID, x, y)]
  dat$bin_ID <- max(dat$x) * (dat$y - 1) + dat$x
  bin.coor <- dat[, sum(V1), by = .(x, y)]


  ##
  geneID <- seq(length(unique(dat$geneID)))
  hash.G <- data.frame(row.names = unique(dat$geneID), values = geneID)
  gen <- hash.G[dat$geneID, 'values']


  ##
  bin_ID <- unique(dat$bin_ID)
  hash.B <- data.frame(row.names = sprintf('%d', bin_ID), values = bin_ID)
  bin <- hash.B[sprintf('%d', dat$bin_ID), 'values']


  ##
  cnt <- dat$V1


  ##
  rm(dat)
  gc()


  ##
  tissue_lowres_image <- matrix(1, max(bin.coor$y), max(bin.coor$x))

  tissue_positions_list <- data.frame(row.names = paste('BIN', rownames(hash.B), sep = '.'),
                                      tissue = 1,
                                      row = bin.coor$y, col = bin.coor$x,
                                      imagerow = bin.coor$y, imagecol = bin.coor$x)

  scalefactors_json <- toJSON(list(fiducial_diameter_fullres = 1,
                                   tissue_hires_scalef = 1,
                                   tissue_lowres_scalef = 1))


  ##
  mat <- sparseMatrix(i = gen, j = bin, x = cnt)

  rownames(mat) <- rownames(hash.G)
  colnames(mat) <- paste('BIN', sprintf('%d', seq(max(hash.B[, 'values']))), sep = '.')

  #seurat_spatialObj <- CreateSeuratObject(mat, project = 'Spatial', assay = 'Spatial',min.cells=5, min.features=5)
  seurat_spatialObj <- CreateSeuratObject(mat, project = 'Spatial', assay = 'Spatial',min.cells=1, min.features=1)

  ##
  generate_spatialObj <- function(image, scale.factors, tissue.positions, filter.matrix = TRUE)
  {
    if (filter.matrix) {
      tissue.positions <- tissue.positions[which(tissue.positions$tissue == 1), , drop = FALSE]
    }

    unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef

    spot.radius <- unnormalized.radius / max(dim(image))

    return(new(Class = 'VisiumV1',
               image = image,
               scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef,
                                            fiducial = scale.factors$fiducial_diameter_fullres,
                                            hires = scale.factors$tissue_hires_scalef,
                                            lowres = scale.factors$tissue_lowres_scalef),
               coordinates = tissue.positions,
               spot.radius = spot.radius))
  }

  spatialObj <- generate_spatialObj(image = tissue_lowres_image,
                                    scale.factors = fromJSON(scalefactors_json),
                                    tissue.positions = tissue_positions_list)


  ##
  spatialObj <- spatialObj[Cells(seurat_spatialObj)]
  DefaultAssay(spatialObj) <- 'Spatial'

  seurat_spatialObj[['slice1']] <- spatialObj


  rm(mat)
  rm(bin.coor)
  rm(hash.G)
  rm(hash.B)
  rm(bin)
  rm(gen)
  rm(cnt)

  #saveRDS(seurat_spatialObj,file=paste0(pro,"_bin",bs,".rds"))
  saveRDS(seurat_spatialObj,file=outfile)

}
