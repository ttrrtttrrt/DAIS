#####  method4: SpaGCN #####
#python $0 outdir h5file posfile pro n_clusters

#1.library
import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
from scipy.sparse import issparse
import random, torch
import warnings
warnings.filterwarnings("ignore")
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import SpaGCN as spg
import cv2

import sys
outdir = sys.argv[1]
h5file = sys.argv[2]
posfile = sys.argv[3]
pro = sys.argv[4]
n_clusters = sys.argv[5]


if not os.path.isdir(outdir):
  os.makedirs(outdir)

os.chdir(outdir)  

#3. Read in data：gene matrix，position.txt，tif(optional)
#Read original data and save it to h5ad
from scanpy import read_10x_h5
adata = read_10x_h5(h5file)
spatial=pd.read_csv(posfile,sep=",",header=None,na_filter=False,index_col=0) 
cells=list(set(adata.obs.index).intersection(set(spatial.index)))
spatial=spatial.loc[cells]
adata.obs["x1"]=spatial[1]
adata.obs["x2"]=spatial[2]
adata.obs["x3"]=spatial[3]
adata.obs["x4"]=spatial[4]
adata.obs["x5"]=spatial[5]
adata.obs["x_array"]=adata.obs["x2"]
adata.obs["y_array"]=adata.obs["x3"]
adata.obs["x_pixel"]=adata.obs["x4"]
adata.obs["y_pixel"]=adata.obs["x5"]

#Select captured samples 
#adata=adata[adata.obs["x1"]==1]
adata.var_names=[i.upper() for i in list(adata.var_names)]
adata.var["genename"]=adata.var.index.astype("str")
adata.write_h5ad(str(pro)+".data.h5ad")


#4. Integrate gene expression and histology into a Graph
#Set coordinates
x_array=adata.obs["x_array"].tolist()
y_array=adata.obs["y_array"].tolist()
x_pixel=adata.obs["x_pixel"].tolist()
y_pixel=adata.obs["y_pixel"].tolist()

#If histlogy image is not available, SpaGCN can calculate the adjacent matrix using the fnction below
adj=spg.calculate_adj_matrix(x=x_pixel,y=y_pixel, histology=False)
np.savetxt(str(pro)+'.adj.csv', adj, delimiter=',')


#5. Spatial domain detection using SpaGCN
#5.1 Expression data preprocessing
#adata=sc.read(str(pro)+".data.h5ad")
#adj=np.loadtxt(str(pro)+'.adj.csv', delimiter=',')
adata.var_names_make_unique()
spg.prefilter_genes(adata,min_cells=3) # avoiding all genes are zeros
spg.prefilter_specialgenes(adata)
#Normalize and take log for UMI
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

#5.2 Set hyper-parameters
p=0.5   
#Find the l value given p
l=spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)

#n_clusters
n_clusters= int(n_clusters)
#Set seed
r_seed=t_seed=n_seed=100
#Search for suitable resolution
res=spg.search_res(adata, adj, l, n_clusters, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)


#5.3 Run SpaGCN
clf=spg.SpaGCN()
clf.set_l(l)
#Set seed
random.seed(r_seed)
torch.manual_seed(t_seed)
np.random.seed(n_seed)
#Run:
clf.train(adata,adj,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)
y_pred, prob=clf.predict()
adata.obs["pred"]= y_pred
adata.obs["pred"]=adata.obs["pred"].astype('category')
#Do cluster refinement(optional)
#shape="hexagon" for Visium data, "square" for ST data.
adj_2d=spg.calculate_adj_matrix(x=x_array,y=y_array, histology=False)
refined_pred=spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), dis=adj_2d, shape="hexagon")
adata.obs["refined_pred"]=refined_pred
adata.obs["refined_pred"]=adata.obs["refined_pred"].astype('category')
#Save results
adata.write_h5ad(str(pro)+".results.h5ad")


#5.4  Plot spatial domains
#adata=sc.read(str(pro)+".results.h5ad")
#adata.obs should contain two columns for x_pixel and y_pixel
#Set colors used
plot_color=["#F56867","#FEB915","#C798EE","#59BE86","#7495D3","#D1D1D1","#6D1A9C","#15821E","#3A84E6","#997273","#787878","#DB4C6C","#9E7A7A","#554236","#AF5F3C","#93796C","#F9BD3F","#DAB370","#877F6C","#268785",
"#1f77b4","#ff7f0e","#279e68","#d62728","#aa40fc","#8c564b",  
"#e377c2","#b5bd61","#17becf","#aec7e8","#ffbb78","#98df8a",  
"#ff9896","#c5b0d5","#c49c94","#FFFF00","#1CE6FF","#FF34FF",  
"#FC0D05","#008941","#006FA6","#A30059","#9400D3","#4CC700B6",
"#3B5DFF","#FF2F80","#FF7C75FA","#00FFCC","#612E21","#FF4A46",  
"#0000A6","#63FFAC","#B79762","#004D43","#8FB0FF","#997D87",  
"#5A0007","#809693","#1B4400","#4FC601"]

#1.Plot spatial domains
domains="pred"
num_celltype=len(adata.obs[domains].unique())
adata.uns[domains+"_colors"]=list(plot_color[:num_celltype])
ax=sc.pl.scatter(adata,alpha=1,x="y_pixel",y="x_pixel",color=domains,title=domains,color_map=plot_color,show=False,size=300000/adata.shape[0])
ax.set_aspect('equal', 'box')
ax.axes.invert_yaxis()
plt.savefig(str(pro)+"_spaGCN.pred.png", dpi=600)
#plt.close()

#2.Plot refined spatial domains
domains="refined_pred"
num_celltype=len(adata.obs[domains].unique())
adata.uns[domains+"_colors"]=list(plot_color[:num_celltype])
ax=sc.pl.scatter(adata,alpha=1,x="y_pixel",y="x_pixel",color=domains,title=domains,color_map=plot_color,show=False,size=300000/adata.shape[0])
ax.set_aspect('equal', 'box')
ax.axes.invert_yaxis()
plt.savefig(str(pro)+"_spaGCN.refined_pred.png", dpi=600)
plt.savefig(str(pro)+"_spaGCN.refined_pred.pdf", dpi=600)
#plt.close()

#plot2
ax=sc.pl.scatter(adata,alpha=1,x="y_pixel",y="x_pixel",color=domains,title="",color_map=plot_color,show=False,size=300000/adata.shape[0])
ax.set_aspect('equal', 'box')
ax.axes.invert_yaxis()
# remove axis
ax.axes.xaxis.set_visible(False)
ax.axes.yaxis.set_visible(False)
# remove legend
plt.legend('', frameon=False)
plt.savefig(str(pro)+"_spaGCN.refined_pred2.png", dpi=600)
#plt.close()
