# DAIS

## Identifying Specific Area based on Density Clustering of Spatial Omics Data

## Introduction
DAIS is an R package for identifying density-based area using spatial omics data. DAIS takes the spot coordinates, spot annotation or signature gene sets as input and outputs the coordinate of each spot and the contour of the identified zones for downstream quantitative analysis.

Several methods of calculating the gene set score within this package were developed based on the originals implemented within [scMetabolism](https://github.com/wu-yc/scMetabolism)

## Installation

### Requirements

```R
install.packages(c("Seurat","VAM","pagoda2"))
BiocManager::install(c("AUCell","GSVA","singscore","UCell")) devtools::install_github("YosefLab/VISION")
devtools::install_github("wu-yc/scMetabolism")
```

### Install
```R
devtools::install_github('ttrrtttrrt/DAIS') 
```

## Tutorials
- Identify stromal and tumor areas using the results of celltype annotation :(https://github.com/ttrrtttrrt/DAIS/blob/main/Tutorials/Fig1BC.anno.R)
- Identify TLS regions using the results of enrichment scores for specific gene sets:(https://github.com/ttrrtttrrt/DAIS/blob/main/Tutorials/Fig.1F.enrich.R)





