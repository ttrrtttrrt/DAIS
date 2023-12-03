# Density clustering with dscan method
set.seed(123)

#' @title Select cell according to gene set score
#'
#' @param metaall dataframe contain gene set score
#' @param genename name of gene set stored in the column names of the metaall
#' @param cutoff value of gene set score to cutoff
#' @param absolute use absolute values(TRUE) or percentages(FALSE),default = FALSE
#'
#' @return meta data frame of filtered spots
#' @export
#'
#' @examples
score.subset <- function (metaall, genename, cutoff, absolute = TRUE) {
  if (absolute == TRUE) {
    meta <- metaall[metaall[,genename]>cutoff,]
  } else {
    cutoff2 <- quantile(metaall[,genename],(1-cutoff))
    meta <- metaall[metaall[,genename]>cutoff2,]
  }
  return(meta)
}


#' @title Density clustering
#' @description Density clustering using the coordinates of spots using dbscan method.
#'
#' @param type.pos data frame contain corordinate of spots
#' @param eps distance measure that will be used to locate the points in the neighborhood of any point
#' @param MinPts the minimum number of points (a threshold) clustered together for a region to be considered dense
#'
#' @return a data frame containing the cell coordinates and cluster id of clustering
#' @export
#'
#' @examples fx <- Dbscan.cluster(type.pos = as.data.frame(matrix(rnorm(100), ncol = 2)), eps = 3, MinPts = 5)
Dbscan.cluster <- function (type.pos, eps = NULL, MinPts) {
  if (is.null(eps)) {
    epsall <- dbscan::kNNdist(type.pos, k = MinPts)
    eps <- quantile(epsall,0.25)
  }
  db <- fpc::dbscan(type.pos, eps = eps, MinPts = MinPts)

  #
  Mer <- cbind(type.pos,db$cluster)
  colnames(Mer)[3] <- c("cluster")
  #filter
  fx <- subset(Mer, cluster != 0)
  return(fx)
}


#' @title Clustering results filtering
#' @description Filter cluster accoring to the number of cell in clusters.
#'
#' @param x data frame contain clustering results
#' @param cutoff the minimum number of spots contained in the cluster to be filtered
#' @param colname the column name of clustering result ID in x
#'
#' @return Clustering results after filtering
#' @export
#'
#' @examples
deletelines <- function(x,cutoff,colname = NULL) {
  if (is.null(colname)) { colname <-  'cluster'}
  stand.col <- x[, colname] # 设根据指定列如cluster列进行删除操作,即根据cell数过滤cluster
  count <- table(stand.col)
  if (all(count < 1)){ stop("no repeated records")}
  else {
    ind <- sapply(stand.col, function(t) ifelse(count[as.character(t)] > as.numeric(cutoff), TRUE, FALSE))   }
  return(x[ind, ])
}


#' @title Specific cluster selecting and density clustering
#' @description Select a specific cluster and perform density clustering.
#'
#' @param fx data frame contain clustering results
#' @param clusterID ID of cluster to select
#' @param k the minimum number of points (a threshold) clustered together for a region to be considered dense
#' @param eps distance measure that will be used to locate the points in the neighborhood of any point
#'
#' @return data frame contains the result of a specific cluster reclustering
#' @export
#'
#' @examples
subsetCluster <- function(fx,clusterID,k,eps){
  mat1 <- subset(fx,cluster==clusterID)
  type.pos <- mat1[,c('row','col')]
  fz <- Dbscan.cluster(type.pos,MinPts = k,eps = eps)
  fz <- deletelines(fz,min_immnum)
  cluster.new <- paste(rep(clusterID,nrow(fz)),fz$cluster,sep=".") #cluster ID replace
  fz$cluster <- cluster.new
  table(fz$cluster)
  return(fz)
}


#' @title Draw the outline with convex hull method
#' @description Get coordinates of polygon point for specific cluster with convex hull method,return a list contains two formats such as polygen and data frame.
#'
#' @param x data frame contain the cluster id and coordinates of cells named cluster,row,col
#' @param subcluster cluster id
#' @param e the expanding radius
#'
#' @return a list outer contour points contains two formats such as polygen and data frame
#' @export
#'
#' @examples
extractCluster <- function(x, subcluster,e=NULL){
  if (is.null(e)) {
    e <- 2
  }
  C0 <- as.matrix(subset(x,cluster==subcluster)[,c('row','col')])
  Cpt0 <- st_multipoint(C0,dim="XY")    #create obj
  Cpg0 <- st_convex_hull(Cpt0) #polygon
  Cb0 <- st_buffer(Cpg0, e) #extend out
  Cor0 <-  as.matrix(Cb0)  #transform polygon to matrix
  mg0 <- cbind(Cor0,rep(subcluster,nrow(Cor0)))  #Outer contour coordinates：x y cluster
  mg0 <- as.data.frame(mg0)
  colnames(mg0) <- c("row","col","cluster")
  mg0$row <- as.numeric(mg0$row)
  mg0$col <- as.numeric(mg0$col)
  #mg0$cluster <- as.factor(mg0$cluster)
  mer0 <- list(Cb0,mg0)  #Two different forms of the outer contour：polygen, data frame
  return(mer0)
}


#' @title Draw the outline with alpha hull (a-shape) method
#' @description Get coordinates of polygon point for specific cluster with alpha hull (a-shape) method,return a list contains two formats such as polygen and data frame.
#'
#' @param x data frame contain the cluster id and coordinates of cells named cluster,row,col
#' @param cluster cluster id
#' @param e the expanding radius
#' @param alpha Value of α
#'
#' @return a list of outer contour points contains two formats such as polygen and data frame
#' @export
#'
#' @examples
extractCluster2 <- function(x, cluster,e, alpha){
  C0 <- as.matrix(subset(x,Region==cluster)[,1:2])
  dat.ashape <- ashape(C0, alpha= alpha)
  #plot(dat.ashape)
  ###get the point order
  edges <- dat.ashape$edges
  mer <- matrix(nrow = dim(edges)+1,ncol = 3) %>% `colnames<-`(c('node','x','y'))
  mer[1:2,1] <- dat.ashape$edges[1,1:2]
  edges2 <- rbind(edges[,c(1,2)],edges[,c(2,1)])
  xlst <- edges2[,'ind1']
  ylst <- edges2[,'ind2']
  for (i in 2:(nrow(mer)-2)) {
    index <- which(xlst==mer[i,'node'])
    index <- index[-which(ylst[index] %in% mer[,'node'])]
    mer[i+1,'node'] <- edges2[index,'ind2']
  }
  mer[nrow(mer),'node'] <- mer[1,'node']

  ###get the point coordinates
  pointdf <- dplyr::distinct(as.data.frame(rbind(edges[,c(1,3,4)],edges[,c(2,5,6)]))) %>%
    `colnames<-`(c('node','x','y'))
  mer[,'x'] <- pointdf[match(mer[,'node'],pointdf[,'node']),'x']
  mer[,'y'] <- pointdf[match(mer[,'node'],pointdf[,'node']),'y']
  #write.table(mer,'ashape_border.txt',sep = "\t",quote = F)

  # ggplot()+
  #    geom_point(as.data.frame(C0),mapping = aes(x=C0[,1],y=C0[,2]),size=1,shape=15)+
  #    theme_test() + geom_path(data = as.data.frame(mer),aes(x=x,y=y)) +
  #    ggrepel::geom_text_repel(data = as.data.frame(mer[-nrow(mer),]),
  #      aes(x=x, y=y, label = node), size=5, fontface="bold", lineheight=1)

  lincC <- st_linestring(as.matrix(mer[,2:3]),dim="XY")  # convert point to line
  # ggplot() + geom_sf(data=lincC)
  newp1 <- st_polygonize(lincC)  # convert line to polygon
  #  ggplot() + geom_sf(data=newp1,color="red")

  Cb0 <- st_buffer(newp1, e) #extend out,polygon
  Cor0 <-  as.matrix(Cb0)
  mg0 <- cbind(Cor0,rep(cluster,nrow(Cor0)))
  mg0 <- as.data.frame(mg0)
  colnames(mg0) <- c("x","y","cluster")
  mg0$x <- as.numeric(mg0$x)
  mg0$y <- as.numeric(mg0$y)
  mg0$cluster <- as.factor(mg0$cluster)
  mer0 <- list(Cb0,mg0)
  return(mer0)
}


#' @title Get coordinates of spots within specific cluster
#'
#' @param frame a data frame that contains the point coordinates for each cluster polygon. The row and col columns will contain the point coordinates as rows and columns, respectively, while the cluster column will indicate the cluster ID to which each point belongs
#' @param annofile a data frame of density clustering result,contains at least 3 columns:coordinates(row,col) and cluster id(cluster)
#' @param n specific cluster id
#' @param p index of specific cluster in polygon list
#'
#' @return a data frame contains coordinates of spots within specific cluster
#' @export
#'
#' @examples
rangexy <- function(frame,annofile,n,p){
  minx <- min(frame$row)
  maxx <- max(frame$row)
  miny <- min(frame$col)
  maxy <- max(frame$col)
  for(i in 1:nrow(annofile)){
    if(annofile[i,2] >= minx && annofile[i,2] <= maxx && annofile[i,3] >= miny && annofile[i,3] <= maxy){
      if(point_in_polygon(frame,c(annofile[i,2],annofile[i,3]))){
        annofile[i,5] <- n
      }
    }
  }
  colnames(annofile)[5] <- "cluster"
  annofile <- subset(annofile,cluster == n )
  return(annofile)
}
point_in_polygon <- function(polygon, point){
  odd = FALSE
  i = 0
  j = nrow(polygon) - 1
  while(i < nrow(polygon) - 1){
    i = i + 1
    if (((polygon[i,2] > point[2]) != (polygon[j,2] > point[2]))
        && (point[1] < ((polygon[j,1] - polygon[i,1]) * (point[2] - polygon[i,2]) / (polygon[j,2] - polygon[i,2])) + polygon[i,1]))
    {
      odd = !odd
    }
    j = i
  }
  return (odd)
}

