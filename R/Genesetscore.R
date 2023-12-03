#' @title Compute the score of specific geneset on single cell level with 9 optional methods
#' @description Calculate the expression levels of each program (geneset) on single cell level with 9 optional methods such as AUCell, ssGSEA, gsva, AddModuleScore, mean, pagoda2, VAM, singscore, UCell.
#'
#'
#' @param object Seurat object
#' @param assays name of assays to store gene matrix
#' @param signaturesFile gmt file of geneset
#' @param method name of calculation method
#' @param imputation whether to use normalized data,default = FALSE
#' @param ncores number of cores used
#'
#' @return a data frame with module scores; each module is stored as name# for each module program present in features; 1 cell per row,1 module per column
#' @export
#'
#' @examples
sc.Meta <- function (object,
                     assays = NULL,
                     signaturesFile,
                     method,
                     imputation = FALSE, ncores = NULL) {
  if (is.null(assays)) {
    assays <- 'RNA'
  }
  if (is.null(ncores)) {
    ncores <- 5
  }
  set.seed(666)

  DefaultAssay(object) <- assays
  object <- NormalizeData(object)
  countexp <- as.data.frame(as.matrix(object[[assays]]@data),stringsAsFactors = F)

  if (imputation == FALSE) {
    countexp2 <- countexp
  } else {
    cat("Start imputation...\n")
    cat("Citation: George C. Linderman, Jun Zhao, Yuval Kluger. Zero-preserving imputation of scRNA-seq data using low-rank approximation. bioRxiv. doi: https://doi.org/10.1101/397588 \n")
    result.completed <- alra(as.matrix(countexp))
    countexp2 <- result.completed[[3]]
    row.names(countexp2) <- row.names(countexp)
  }


  cat("Start quantify the geneset score\n")


  ############################## Method1:AUCell  ##############################
  if (method == "AUCell") {
    library(AUCell)
    library(GSEABase)
    cells_rankings <- AUCell_buildRankings(as.matrix(countexp2),
                                           nCores = ncores,
                                           plotStats = F)
    geneSets <- getGmt(signaturesFile)
    cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
    signature_exp <- data.frame(t(getAUC(cells_AUC)))
  }


  ############################## Method2:ssGSEA  ##############################
  if (method == "ssGSEA") {
    library(GSVA)
    library(GSEABase)
    geneSets <- getGmt(signaturesFile)
    gsva_es <- gsva(as.matrix(countexp2), geneSets, method = c("ssgsea"),
                    kcdf = c("Poisson"), parallel.sz = ncores)
    signature_exp <- data.frame(t(gsva_es))
  }


  ############################## Method3:gsva  ##############################
  if (method == "gsva") {
    library(GSVA)
    library(GSEABase)
    geneSets <- getGmt(signaturesFile)
    gsva_es <- gsva(as.matrix(countexp2), geneSets, method = c("gsva"),
                    kcdf = c("Poisson"), parallel.sz = ncores)
    signature_exp <- data.frame(t(gsva_es))
  }


  ############################## Method4:AddModuleScore(Seurat)  ##############################
  if (method == "AddModuleScore"){
    library(GSEABase)
    object1 <- object
    gene <- read.gmt(signaturesFile)
    geneSet <- split(gene$gene,gene$term)
    object1 <- AddModuleScore(object1,features =  geneSet,name = paste(names(geneSet),"score",sep="_"),
                              seed =666,ctrl=100,nbin = 24)
    colnames(object1@meta.data) <- gsub("_score[0-9]*", "", colnames(object1@meta.data))
    signature_exp <- data.frame(object1@meta.data[,intersect(colnames(object1@meta.data),names(geneSet)),drop=F])
  }


  ############################## Method5:mean  ##############################
  if (method == "mean"){
    library(pagoda2)
    library(GSEABase)
    geneSets <- getGmt(signaturesFile)
    paths <- names(geneSets)
    mat <- t(countexp2)
    enrich <- matrix(nrow=dim(mat)[1], ncol=length(paths),dimnames = list(rownames(mat),paths))
    count <- 0
    for(i in 1:length(paths)){
      geneset1 <- names(geneSets[i])
      count <- count + 1
      enrich[,count] <-  score.cells.nb0(mat,geneSets[[i]]@geneIds) # 基因集是character
    }
    signature_exp <- data.frame(enrich)
  }


  ############################## Method6:pagoda2  ##############################
  if (method == "pagoda2"){
    library(pagoda2)
    library(GSEABase)
    geneSets <- getGmt(signaturesFile)
    paths <- names(geneSets)
    mat <- t(countexp2)
    enrich <- matrix(nrow=dim(mat)[1], ncol=length(paths),dimnames = list(rownames(mat),paths))
    count <- 0
    for(i in 1:length(paths)){
      geneset1 <- names(geneSets[i])
      count <- count + 1
      enrich[,count] <-  score.cells.puram(mat,geneSets[[i]]@geneIds,correct=TRUE,show.plot=FALSE) # 基因集是character
    }
    signature_exp <- data.frame(enrich)
  }


  ############################## Method7:VAM  ##############################
  if (method == "VAM"){
    library(VAM)
    object1 <- object
    #transfer geneSymbol to geneID
    feature.data = read.delim("HG.tsv",header = T, stringsAsFactors = FALSE)
    ensembl.ids = feature.data[,1]
    gene.names = feature.data[,2]
    genes.after.QC = rownames(countexp2)
    indices.to.keep = unlist(sapply(genes.after.QC, function(x) {which(gene.names == x)[1]}))
    ensembl.ids = ensembl.ids[indices.to.keep]
    #readin geneset
    library(GSEABase)
    geneSets <- getGmt(signaturesFile)
    paths <- names(geneSets)
    enrich <- matrix(nrow=dim(countexp2)[2], ncol=length(paths),dimnames = list(colnames(countexp2),paths))
    count <- 0
    for(i in 1:length(paths)){
      gene.set.name = paths[i]
      blymphocyte.gene.ids = geneSets[[i]]@geneIds
      # Create a collection list for this gene set based on the Ensembl IDs
      gene.set.id.list = list()
      gene.set.id.list[[1]] = blymphocyte.gene.ids
      names(gene.set.id.list)[1] = gene.set.name
      gene.set.id.list
      gene.set.collection = createGeneSetCollection(gene.ids=ensembl.ids, gene.set.collection=gene.set.id.list)
      object1 <- vamForSeurat(seurat.data=object1, gene.set.collection= gene.set.collection,  center=F, gamma=T, sample.cov=F, return.dist=T)
      count <- count +1
      enrich[,count] <- as.matrix(object1@assays$VAMcdf@data)
    }
    signature_exp <- data.frame(enrich)
  }


  ############################## Method8:singscore  ##############################
  if (method == "singscore") {
    library(singscore)
    library(GSEABase)
    geneSets <- getGmt(signaturesFile)
    paths <- names(geneSets)
    enrich <- matrix(nrow=dim(countexp2)[2], ncol=length(paths),dimnames = list(colnames(countexp2),paths))  # 1 cell per row,1 pathway per column
    count <- 0
    for(i in 1:length(paths)){
      geneset1 <- names(geneSets[i])
      count <- count + 1
      Mat <-  rankGenes(as.data.frame(countexp2))
      Mat2 <- simpleScore(Mat, upSet =geneSets[[i]]) #可设置downSet,默认NULL
      enrich[,count] <- Mat2[,1]
    }
    signature_exp <- data.frame(enrich)
  }


  ############################## Method9:UCell  ##############################
  if (method == "UCell") {
    library(UCell)
    library(clusterProfiler)
    object1 <- object
    # slot(obj$SCT@SCTModel.list[[1]], 'median_umi') = median(obj$SCT@SCTModel.list[[1]]@cell.attributes$umi)
    gene <- read.gmt(signaturesFile)
    geneSet <- split(gene$gene,gene$term)
    input.pathway <- unique(gene$term)
    object1 <- AddModuleScore_UCell(object1, features = geneSet)
    colnames(object1@meta.data) <- gsub("_UCell", "", colnames(object1@meta.data))
    signature_exp <- data.frame(object1@meta.data[,(ncol(object1@meta.data)-length(geneSet)+1):ncol(object1@meta.data),drop=F])
  }

  return(signature_exp)
}
