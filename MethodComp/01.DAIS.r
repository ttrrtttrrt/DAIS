#####  method1:DAIS   #####
#Rscript $0 pro inrds pc res

args <- commandArgs(T)
outdir <- args[1]
inrds <- args[2]  
pre <- args[3]
topprop <- as.numeric(args[4])  
k1 <- as.numeric(args[5]) 
eps1 <- as.numeric(args[6])  
r <- as.numeric(args[7])  
min_immnum <- as.numeric(args[8])  
scorecut <- as.numeric(args[9])  

library(DAIS)
library(ggplot2)
library(Seurat)

if (!dir.exists(outdir)) {
  dir.create(outdir)
}
setwd(outdir)

obj <- readRDS(inrds)
obj[["RNA"]] <- CreateAssayObject(data = obj@assays$Spatial@counts) 

#-----score compute-----
method <- 'AUCell'
signaturesFile <- '/hwfssz1/ST_HEALTH/P21Z10200N0134/USER/tianru/01.TLS/01.DAIS/04.FFPE_TLS/01.DAIS/00.marker/TLS.signature.v2.gmt'
score <- sc.Meta(object = obj,assays = "RNA",signaturesFile = signaturesFile,method = method,imputation = F)
coord <- obj@images$image@coordinates
signature_exp <- cbind(coord,score[rownames(coord),])
colnames(signature_exp)[1] <- 'Row.names'
signature_exp$Row.names <- rownames(signature_exp)
write.table(signature_exp, file=paste(pre,method,'rds.score.tsv',sep="."),quote=F,sep="\t")


#-----clustering-----
pal <- c("#1f77b4","#ff7f0e","#279e68","#d62728","#aa40fc","#8c564b",
         "#e377c2","#b5bd61","#17becf","#aec7e8","#ffbb78","#98df8a",
         "#ff9896","#c5b0d5","#c49c94","#FFFF00","#1CE6FF","#FF34FF",
         "#FC0D05","#008941","#006FA6","#A30059","#9400D3","#4CC700B6",
         "#3B5DFF","#FF2F80","#FF7C75FA","#00FFCC","#612E21","#FF4A46",
         "#0000A6","#63FFAC","#B79762","#004D43","#8FB0FF","#997D87",
         "#5A0007","#809693","#1B4400","#4FC601")

############################## 00.Seurat object data readin   ##############################
#signature_exp <- read.table(scorefile,sep = "\t",header = T)
merg <- signature_exp

############################## 02.cell select ##############################
# part3 1E
colnames(merg)[8] <- "Score"
quantile(merg[,8],seq(0.9,1,0.01)) # 99%
if (is.na(scorecut)) {
  scorecut <- round(quantile(merg[,8],(1-topprop)),3)
}


p2 <- ggplot(merg,aes(x=row,y=col)) +
  geom_point(aes(fill = Score),stroke=NA,shape=23, alpha = 0.5, size=2) +
  scale_fill_gradientn(colours = c("#0776B3","#faf6bf","red"))  +
  ggtitle("Geneset=Overlap, method=AUCell") +
  theme_test()
#print(p2)


f1 <- subset(merg,Score >= scorecut)
p3 <- ggplot(f1,aes(x=row,y=col)) +  geom_point(aes(fill = Score),stroke=NA,shape=23,alpha = 1, size=2.4)+  theme_test() +
  scale_fill_gradientn(colours = c("#faf6bf","red")) + ggtitle(paste("TLS cutoff is ",scorecut,sep=" "))
#print(p3)


############################## 03.Density clustering ##############################
type.pos <- f1[,c(2,3)] #spots used to cluster
#dbscan::kNNdistplot(type.pos, k =  k1)


fx <- Dbscan.cluster(type.pos,eps = eps1,MinPts = k1)
table(fx$cluster)
fx <- deletelines(fx,min_immnum)  # result of clustering
table(fx$cluster)


###Clustering result graph after filter
p4 <- ggplot(fx,aes(x=row,y=col)) +
  geom_point(aes(fill = as.factor(cluster)), shape=23,stroke=NA,alpha = 1,size=3)+
  scale_fill_manual(values =  pal) +
  #coord_fixed(ratio = 1) +
  guides(fill=guide_legend(title="Cluster")) +
  theme_test() +
  ggtitle("Cluster for TLS prediction")
#print(p4)
stat1 <- as.data.frame(table(fx$Region))


############################## 04.Contour recognition ##############################
fy <- fx
p <- 1
mer1 <- list()
mer2 <- list()

clusterX <- unique(fy$cluster)
for(i in clusterX){
  C1 <- extractCluster(fy, i,r)  #Spread radius = 0.5
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
dataframes_with_label <- lapply(seq_along(mer2), function(i) {
  df <- mer2[[i]]
  df$source <- as.factor(paste0("Dataframe ", i))
  return(df)
})
combined_data <- bind_rows(dataframes_with_label)

###outline of TLS region in all cell map 
p7 <- ggplot(merg,aes(x=row,y=col)) +
  geom_point(aes(fill = Score), shape=23,alpha = 1, size=2,stroke=NA)+
  theme_test() +
  scale_fill_gradientn(colours = c("#0776B3","#faf6bf","red")) +
  geom_path(data=combined_data,aes(x=row,y=col,group=cluster),color="black") +
  labs(x="x",y="y") +
  ggtitle("TLS candidate region")
#print(p7)


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

write.table(fx,file=paste0("k",k1,"eps",eps1,".cutoff_",scorecut,".tsv",sep="."),quote=F,sep="\t")
write.table(dex2,file=paste0("k",k1,".eps",eps1,".cluster.tsv"),quote=F,sep="\t")

pdf(paste0(pre,"_10X.enrich.pdf"),6,6)
print(p2)
print(p3)
print(p4)
print(p7)
dev.off()

p8 <- ggplot(merg,aes(x=row,y=col)) +
  geom_point(aes(fill = Score), shape=23,alpha = 1, size=2,stroke=NA)+
  theme_test() +
  scale_fill_gradientn(colours = c("#0776B3","#faf6bf","red")) +
  geom_path(data=combined_data,aes(x=row,y=col,group=cluster),color="black") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none')
png(filename = paste0(pre,"_DAIS.enrich.png"),width = 400,height = 400)
print(p8)
dev.off()

#save.image(paste0(pre,".cluster.RData"))



