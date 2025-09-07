#细胞分割

#配准align
#质控
#分群聚类
#自动化注释

#L1-early_b_new
import os,sys,getopt,re
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal
def_path = '/hwfssz5/ST_SUPERCELLS/P21Z10200N0206/jiazihan/'
sys.path.append(def_path)
from supplement import bin1tobinx,space_multi,lasso_pre
from supplement import basic_stats, basic_stats_multi, map_counts
from supplement import symbol2fbgn,gene2tissue,enrichment_sf,anno_fullstage

import warnings
import pandas as pd
import spateo as st
warnings.filterwarnings("ignore")
import math
import operator
from typing import Optional, Union
import matplotlib as mpl
import numpy as np
from anndata import AnnData
from scipy.sparse import isspmatrix
from sklearn.preprocessing import minmax_scale
import anndata as ad

import matplotlib.pyplot as plt
import seaborn as sns

import dynamo as dyn
import scanpy as sc

import logging as logg
import time

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal 

from scipy.sparse import issparse

#############################
#name = sys.argv[1]
name = 'L3-early_new'
project = name
'''
#############################
########### align ###########
#############################

binning_folder = '/hwfssz5/ST_SUPERCELLS/P21Z10200N0206/project/drosophila/03.ST/03.basic_anlysis/second_batch/'+name+'/01_segmentation/bin20_data_2'
cellbin_folder = '/hwfssz5/ST_SUPERCELLS/P21Z10200N0206/project/drosophila/03.ST/03.basic_anlysis/second_batch/'+name+'/01_segmentation/final_seg_data_2'
save_folder = '/hwfssz5/ST_SUPERCELLS/P21Z10200N0206/project/drosophila/03.ST/03.basic_anlysis/second_batch/'+name+'/02_alignment'

##############################
save_binning_folder = os.path.join(save_folder, "bin20_align_h5ad")
save_binning_image_folder = os.path.join(save_folder, "bin20_align_image")

save_cellbin_folder = os.path.join(save_folder, "cellbin_align_h5ad")
save_cellbin_image_folder = os.path.join(save_folder, "cellbin_align_image")

for folder in [binning_folder, save_binning_folder, save_binning_image_folder, save_cellbin_folder, save_cellbin_image_folder]:
    if not os.path.exists(folder):
        os.mkdir(folder)

# Some parameters in the visualization part need to be adjusted
########################################################################################################################

#读取
binning_files = [filename for root1, dirs1, files in os.walk(binning_folder) for filename in files]
binning_files.sort()
print(binning_files)
binning_slices_raw = [ad.read(os.path.join(binning_folder, bin_file)) for bin_file in binning_files]

cellbin_files = [filename for root2, dirs2, files in os.walk(cellbin_folder) for filename in files]
cellbin_files.sort()
print(cellbin_files)
cellbin_slices_raw = [ad.read(os.path.join(cellbin_folder, bin_file)) for bin_file in cellbin_files]

#配准
cellbin_slices, binning_slices = st.tl.slices_align_ref(
    slices=cellbin_slices_raw,
    slices_ref=binning_slices_raw,
    spatial_key="spatial",
    key_added="align_spatial",
    dissimilarity="euc",
    numItermax=500,
    numItermaxEmd=500000,
     device="cpu"
)


for adata in binning_slices:
    adata.write_h5ad(filename=os.path.join(save_binning_folder, f'{adata.obs["slices"][0]}.h5ad'), compression="gzip")

for adata in cellbin_slices:
    adata.write_h5ad(filename=os.path.join(save_cellbin_folder, f'{adata.obs["slices"][0]}.h5ad'), compression="gzip")


#整合配准后的bin和seg
binning_adata_raw = st.tl.integrate(
    adatas=[i.copy() for i in binning_slices_raw], batch_key="slices"
)
cellbin_adata_raw = st.tl.integrate(
    adatas=[i.copy() for i in cellbin_slices_raw], batch_key="slices"
)

space_multi(adata=binning_adata_raw, groupby=None, spatial_key="spatial", slice_key="slices", filename=os.path.join(save_binning_image_folder, "binning_alignment_all_raw.pdf"))
space_multi(adata=cellbin_adata_raw, groupby=None, spatial_key="spatial", slice_key="slices", filename=os.path.join(save_cellbin_image_folder, "cellbin_alignment_all_raw.pdf"))

bin20_adata = st.tl.integrate(adatas=binning_slices, batch_key="slices")
cellbin_adata = st.tl.integrate(adatas=cellbin_slices, batch_key="slices")

bin20_adata.write_h5ad(filename=os.path.join(save_folder, f'{bin20_adata.obs["slices"][0][:-4]}_bin20.h5ad'), compression="gzip")
cellbin_adata.write_h5ad(filename=os.path.join(save_folder, f'{cellbin_adata.obs["slices"][0][:-4]}_cellbin.h5ad'), compression="gzip")

space_multi(adata=bin20_adata, groupby=None, spatial_key="align_spatial", slice_key="slices", filename=os.path.join(save_binning_image_folder, "binning_alignment_all.pdf"))

space_multi(adata=cellbin_adata, groupby=None, spatial_key="align_spatial", slice_key="slices", filename=os.path.join(save_cellbin_image_folder, "cellbin_alignment_all.pdf"))


#############################
############ qc #############
#############################

project = name
align_cellbin_folder = '/hwfssz5/ST_SUPERCELLS/P21Z10200N0206/project/drosophila_1212/03.ST/03.basic_anlysis/second_batch/'+name+'/02_alignment/cellbin_align_h5ad'
save_folder = '/hwfssz5/ST_SUPERCELLS/P21Z10200N0206/project/drosophila_1212/03.ST/03.basic_anlysis/second_batch/'+name+'/03_qc_normalization'

#############################

print('===== qc & normalization start! =====')

save_h5ad_folder = os.path.join(save_folder, "cellbin_norm_h5ad")
save_integrate_folder = os.path.join(save_folder, "cellbin_integrate_h5ad")
save_bin_norm_file = save_integrate_folder+'/'+project+"_cellbin_norm.h5ad"
save_basic_stats_img1_file = save_integrate_folder+'/'+project+"_cellbin_norm_beforeqc"
save_basic_stats_img2_file = save_integrate_folder+'/'+project+"_cellbin_norm_afterqc"
save_area_img_file = save_integrate_folder+'/'+project+"_cellbin_norm_afterqc_area.pdf"

for folder in [save_folder, save_h5ad_folder, save_integrate_folder]:
    if not os.path.exists(folder):
        os.mkdir(folder)

##################################################################
bin_files = [filename for root1, dirs1, files in os.walk(align_cellbin_folder) for filename in files]
bin_files.sort()

# preprocessing
raw_slices = []
norm_slices = []

for bin_file in bin_files:
    adata = ad.read(os.path.join(align_cellbin_folder, bin_file))
    z = adata.obs["slices"].map(lambda x: int(x[-16:-15]) * 14)
    adata.obsm["align_spatial"] = np.c_[adata.obsm["align_spatial"], z]
    basic_stats(adata=adata, mito_label="mt:")
    raw_slices.append(adata)
    #计算线粒体
    adata = adata[adata.obs.pMito < 0.1, :]
    st.pp.filter.filter_cells(adata=adata, min_area=20, min_expr_genes=20, inplace=True)
    st.pp.filter.filter_genes(adata=adata, min_cells=3, min_counts=1, inplace=True)
    basic_stats(adata=adata, mito_label="mt:")
    #皮尔逊残差
    st.tl.pearson_residuals(adata=adata, n_top_genes=None)
    adata.write_h5ad(filename=os.path.join(save_h5ad_folder, bin_file), compression="gzip")
    adata.layers["raw_count"] = adata.X.copy()
    
    norm_adata = adata.copy()
    norm_adata.X = norm_adata.obsm["pearson_residuals"]
    del norm_adata.layers, norm_adata.obsm["pearson_residuals"]
    norm_slices.append(norm_adata)

# integration
raw_adata = st.tl.integrate(adatas=raw_slices, batch_key="slices")
norm_adata = st.tl.integrate(adatas=norm_slices, batch_key="slices")
new_norm_adata = map_counts(norm_adata=norm_adata, counts_adata=raw_adata, key_added="counts_X")
new_norm_adata.write_h5ad(filename=save_bin_norm_file, compression="gzip")

space_multi(adata=new_norm_adata, groupby=None, slice_key="slices", spatial_key="align_spatial", point_size=0.2,filename=save_area_img_file, colormap="black")

# basic_stats before qc
basic_stats_multi(adatas=raw_slices, slice_key="slices", inner="box", save_show_or_return="save",save_kwargs={"path": save_basic_stats_img1_file, "prefix": 'basic_stats', "dpi": 300, "ext": 'pdf'})

# basic_stats after qc
basic_stats_multi(adatas=norm_slices, slice_key="slices", inner="box", save_show_or_return="save",save_kwargs={"path": save_basic_stats_img2_file, "prefix": 'basic_stats', "dpi": 300, "ext": 'pdf'})

print('===== cellbin qc & normalization finished! =====')

'''
#############################
############ clu ############
#############################

cellbin_h5ad = '/hwfssz5/ST_SUPERCELLS/P21Z10200N0206/project/drosophila_1212/03.ST/03.basic_anlysis/second_batch/'+name+'/03_qc_normalization/cellbin_integrate_h5ad/'+name+'_cellbin_norm.h5ad'
outdir = '/hwfssz5/ST_SUPERCELLS/P21Z10200N0206/project/drosophila_1212/03.ST/03.basic_anlysis/second_batch/'+name+'/04_clustering'
#############################

def get_pca_components(adata):
    from sklearn.decomposition import PCA
    target = 0.9
    sc.tl.pca(adata, svd_solver='arpack')
    res = PCA(n_components=target).fit_transform(adata.varm['PCs'])
    print("original shape: ", adata.varm['PCs'].shape)
    print("transformed shape: ", res.shape)
    return res.shape[1]

def cluster(norm_adata,save_folder,resolution,sample_name = "S"):
    key = "scc_"+str(resolution)
    save_image_folder = os.path.join(save_folder,key)
    if not os.path.exists(save_image_folder):
        os.makedirs(save_image_folder)
    st.tl.scc(
        adata=norm_adata, 
        spatial_key="align_spatial", key_added=key, pca_key="X_pca", 
        e_neigh=30, s_neigh=6,cluster_method="leiden", 
        resolution=resolution, 
        copy=False
    )
    space_multi(
        adata=norm_adata, 
        groupby=key, spatial_key="align_spatial", slice_key="slices",
        point_size=0.2, 
        filename= f'{save_image_folder}/{sample_name}_{key}_spatial_cluster.pdf')
    dyn.pl.umap(
            norm_adata, 
            color=[key, "slices"], 
            pointsize=0.1, alpha=0.8, theme='blue', figsize=(10, 5),show_legend="upper left", save_show_or_return='save',
            save_kwargs={"path": f'{save_image_folder}/', "prefix": f'{sample_name}_{key}_umap', "dpi": 300, "ext": 'pdf'}
        )
    
    sc.tl.rank_genes_groups(adata=norm_adata, groupby=key, n_genes=50, method="t-test")
    dedf = sc.get.rank_genes_groups_df(norm_adata, key="rank_genes_groups",group =None)
    dedf.to_csv(f'{save_image_folder}/{sample_name}_{key}_marker.csv', sep=",", index=False)
    sc.tl.dendrogram(norm_adata,groupby = key)
    del norm_adata.uns[key+"_colors"], norm_adata.uns["slices_colors"]
    sc.pl.rank_genes_groups_heatmap(adata=norm_adata, key="rank_genes_groups", 
                                    n_genes=10,vmin=0.3,vmax = 3,
                                    groupby = key,
                                    show_gene_labels=True,
                                    cmap="hot_r",
                                    var_group_labels=False, swap_axes=True, show=False)
    plt.tight_layout()
    plt.savefig(f'{save_image_folder}/{sample_name}_{key}_heatmap.pdf', dpi=300)
    return(norm_adata,set(dedf.groupby('group').head(20).names))

norm_h5ad = cellbin_h5ad
save_folder = os.path.join(outdir,'cellbin')
save_bin_norm_file = os.path.join(save_folder,f'{project}_cellbin.h5ad')
res = [1.2,1.4,1.6,1.8,2.0,2.2]
norm_adata = ad.read(norm_h5ad)
norm_adata.X = norm_adata.X.astype(np.float64)

print(norm_adata.X[np.isnan(norm_adata.X)])

st.tl.pca_spateo(adata=norm_adata, n_pca_components=50, pca_key="X_pca")
# umap
dyn.tl.reduceDimension(adata=norm_adata, basis="umap", n_components=2)
del norm_adata.uns["umap_fit"]  

## clustering
marker_genes = []

for resolution in res:
    norm_adata,genes_1 = cluster(norm_adata,save_folder,resolution = resolution,sample_name = project)
    marker_genes.append(genes_1)
    print('===== resolution =',resolution,' finished ! =====')

norm_adata.write_h5ad(filename=save_bin_norm_file, compression="gzip")

#pic
genelist = marker_genes[0]|marker_genes[1]|marker_genes[2]|marker_genes[3]|marker_genes[4]|marker_genes[5]

save_bin_genes_pattern_folder = os.path.join(save_folder, f"makrer")
if not os.path.exists(save_bin_genes_pattern_folder):
    os.makedirs(save_bin_genes_pattern_folder)
    
for gene in genelist:
    gene = str(gene)
    gene_filename = gene.replace(":", "_") if ":" in gene else gene
    filename = os.path.join(save_bin_genes_pattern_folder, f"{gene_filename}.png")
    space_multi(adata=norm_adata, groupby=gene, spatial_key="align_spatial", slice_key="slices", point_size=0.2, filename=filename)


print('===== cellbin clustering finished ! =====')

#############################
######### auto-anno #########
#############################

adata_path = '/hwfssz5/ST_SUPERCELLS/P21Z10200N0206/project/drosophila_1212/03.ST/03.basic_anlysis/second_batch/'+name+'/04_clustering/cellbin/'+name+'_cellbin.h5ad'

outdir = '/hwfssz5/ST_SUPERCELLS/P21Z10200N0206/project/drosophila_1212/03.ST/03.basic_anlysis/second_batch/'+name+'/04_clustering/cellbin/'

#############################
def_path = '/hwfssz5/ST_SUPERCELLS/P21Z10200N0206/public/drosophila_pipe_v1/'
sys.path.append(def_path)

fbgn_path = '/hwfssz5/ST_SUPERCELLS/P21Z10200N0206/public/drosophila_pipe_v1/supplement/autoanno_file/deml_fbgn.tsv.gz'
fullstage_path = '/hwfssz5/ST_SUPERCELLS/P21Z10200N0206/public/drosophila_pipe_v1/supplement/autoanno_file/full_stage_anno.csv'

adata = ad.read(adata_path)
adata_autosave = os.path.join(outdir,name+'_cellbin_autoanno.h5ad')
time = re.split('_',adata.obs['slices'][0])[0]
res_list = adata.obs.columns[-6:].tolist()

for i in res_list:
    out = os.path.join(outdir,i)
    save_umap_autoanno =  os.path.join(out,i+'_umap_autoanno.pdf')
    save_autoanno_cluster = os.path.join(out,i+'_cluster_autoanno.pdf')
    save_autoanno_heatmap = os.path.join(out,i+'_heatmap_autoanno.pdf')
    save_autoanno_csv = os.path.join(out,i+'_autoanno.csv')
    #rank_genes_groups
    sc.tl.rank_genes_groups(adata, i, method='t-test')
    anno_fullstage(adata=adata,times=[time],gene_nametype="symbol",
                       groupby=i,groups=None,
                       n_genes=50,
                       obs_key=i+"_auto_anno",
                       key="rank_genes_groups",
                       key_added=i+"_marker_genes_anno",
                       copy=False,
                       fbgn_path=fbgn_path,
                       fullstage_path=fullstage_path)
    #dotplt
    space_multi(
                adata=adata, 
                groupby=i+"_auto_anno", spatial_key="align_spatial", slice_key="slices",
                point_size=0.2, 
                filename=save_autoanno_cluster
            )
    #umap
    sc.pl.umap(adata,color=i+"_auto_anno")
    plt.tight_layout()
    plt.savefig(save_umap_autoanno, dpi=300,bbox_inches = 'tight')
    #heatmap
    sc.pl.rank_genes_groups_heatmap(
                adata=adata, 
                key="rank_genes_groups", groupby=i+"_auto_anno", groups=None,
                n_genes=10, vmin=0, vmax=5, 
                min_logfoldchange=None, show_gene_labels=True,
                cmap="hot_r", var_group_labels=False,swap_axes=True, 
                show=False
            )
    plt.tight_layout()
    plt.savefig(save_autoanno_heatmap, dpi=300,bbox_inches = 'tight')
    #csv
    df = pd.DataFrame.from_dict(adata.uns[i+"_marker_genes_anno"],orient='index')
    df.columns = [i+'_auto_anno']
    df['cluster'] = df.index
    df.to_csv(save_autoanno_csv,sep=",", index = False)
    print(f'------ Auto annotation {time} {i} finished!')

adata.write_h5ad(adata_autosave,compression='gzip')
print(f'------ Auto annotation result saved in {adata_autosave}!')
