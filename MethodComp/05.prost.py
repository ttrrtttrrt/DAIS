#####  method5: PROST #####
#python $0 outdir inputdir h5filename pro n_clusters

#1.library
import sys

import pandas as pd 
import numpy as np 
import scanpy as sc 
import os 
import warnings 
warnings.filterwarnings("ignore") 
import matplotlib as mpl 
import matplotlib.pyplot as plt 
import PROST 

#2.Set up the working environment and import data 
# the location of R (used for the mclust clustering)
ENVpath = "/hwfssz5/ST_HEALTH/P21Z10200N0134/07.TLS/2.tools/Miniconda3/envs/PROST_ENV"      # refer to 'How to use PROST' section
os.environ['R_HOME'] = f'{ENVpath}/lib/R'
os.environ['R_USER'] = f'{ENVpath}/lib/python3.7/site-packages/rpy2'

# Set seed
SEED = 818
PROST.setup_seed(SEED)


import sys
output_dir = sys.argv[1]
input_dir = sys.argv[2]   
count_file = sys.argv[3]  #filename
pro = sys.argv[4]
n_clusters = sys.argv[5]


# Set directory
if not os.path.isdir(output_dir):
  os.makedirs(output_dir)

os.chdir(output_dir)  

#%% Read in data
# Read data from input_dir
adata = sc.read_visium(path=input_dir, count_file=count_file)
adata.var_names_make_unique()

#-----SVG-----
#3.Calculate and save PI
# Calculate PI
adata = PROST.prepare_for_PI(adata, platform="visium") 
adata = PROST.cal_PI(adata, platform="visium")  

# Calculate spatial autocorrelation statistics and do hypothesis test
PROST.spatial_autocorrelation(adata, k = 10, permutations = None)

# Save PI result
adata.write_h5ad(str(pro)+"_PI_result.h5")


#4.Draw SVGs detected by PI:
# Draw SVGs detected by PI
#PROST.plot_gene(adata, platform="visium",size = 2, sorted_by = "PI", top_n = 25,save_path = output_dir)


#-----clustering-----
# Set the number of clusters
n_clusters = int(n_clusters)

#1.Read PI result and Expression data preprocessing:
#adata = sc.read(str(pro)+"_PI_result.h5")
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata = PROST.feature_selection(adata, by = "prost", n_top_genes = 3000)

#2.Run PROST clustering
PROST.run_PNN(adata,
              platform="visium",
              key_added = "PROST",
              init="mclust",                         
              n_clusters = n_clusters,                        
              lap_filter = 2,                                  
              lr = 0.1,                         
              SEED=SEED,                          
              max_epochs = 500,                        
              tol = 5e-3,                        
              post_processing = True,                        
              pp_run_times = 3,
              cuda = False)

#3.Save result
adata.write_h5ad(str(pro)+"_PNN_result.h5")   
clustering = adata.obs
clustering.to_csv(str(pro)+"_clusters.csv",header = True)
embedding = adata.obsm["PROST"]
np.savetxt(str(pro)+"_embedding.txt",embedding) 

#4.Plot clustering results
# Read data
#adata = sc.read(str(pro)+"_PNN_result.h5")

# Set colors
plot_color=["#F56867","#FEB915","#C798EE","#59BE86","#7495D3","#D1D1D1","#6D1A9C","#15821E","#3A84E6","#997273","#787878","#DB4C6C","#9E7A7A","#554236","#AF5F3C","#93796C","#F9BD3F","#DAB370","#877F6C","#268785",
"#1f77b4","#ff7f0e","#279e68","#d62728","#aa40fc","#8c564b",  
"#e377c2","#b5bd61","#17becf","#aec7e8","#ffbb78","#98df8a",  
"#ff9896","#c5b0d5","#c49c94","#FFFF00","#1CE6FF","#FF34FF",  
"#FC0D05","#008941","#006FA6","#A30059","#9400D3","#4CC700B6",
"#3B5DFF","#FF2F80","#FF7C75FA","#00FFCC","#612E21","#FF4A46",  
"#0000A6","#63FFAC","#B79762","#004D43","#8FB0FF","#997D87",  
"#5A0007","#809693","#1B4400","#4FC601"]
cmp = mpl.colors.ListedColormap(plot_color)

# Plot clustering results
plt.rcParams["figure.figsize"] = (4,4)
sc.pl.spatial(adata, 
              img_key = "hires", 
              color = ["clustering","pp_clustering"],
              title = ['clustering','post-processedclustering'],                
              na_in_legend = False,
              ncols = 2,
              size = 1)
plt.savefig(str(pro)+"_PROST_pred.pdf", dpi=600)

#plot result of pp_clustering(post-processedclustering)
plt.rcParams["figure.figsize"] = (4,4)
sc.pl.spatial(adata, 
              img_key = "hires", 
              color = ["pp_clustering"],
              title = [''],                
              na_in_legend = False,
              ncols = 1,
              size = 1,
              frameon=False)  
#remove legend
plt.legend('', frameon=False)
plt.savefig(str(pro)+"_PROST_pred2.png", dpi=600)

