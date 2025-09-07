#E14-16h_d_new
import os,sys,getopt,re
def_path = '/hwfssz5/ST_SUPERCELLS/P21Z10200N0206/public/drosophila_pipe_v1/'
sys.path.append(def_path)
#from supplement import bin1tobinx,space_multi,lasso_pre
#from supplement import basic_stats, basic_stats_multi, map_counts
#from supplement import symbol2fbgn,gene2tissue,enrichment_sf,anno_fullstage

import warnings
import pandas as pd
import dynamo as dyn
import spateo as st
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

warnings.filterwarnings("ignore")
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
name = 'L3-early_all'
project = name

#############################
########### align ###########
#############################

binning_folder = '/hwfssz5/ST_SUPERCELLS/P21Z10200N0206/project/drosophila_1212/03.ST/03.basic_anlysis/second_batch/'+name+'/01_segmentation/bin20_data'
cellbin_folder = '/hwfssz5/ST_SUPERCELLS/P21Z10200N0206/project/drosophila_1212/03.ST/03.basic_anlysis/second_batch/'+name+'/01_segmentation/final_seg_data'
save_folder = '/hwfssz5/ST_SUPERCELLS/P21Z10200N0206/project/drosophila_1212/03.ST/03.basic_anlysis/second_batch/'+name+'/02_alignment'

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



binning_slices, binning_slices_ref, pis = st.tl.models_align_ref(
    models=binning_slices_raw,
    models_ref=None,
    n_sampling=2000,
    sampling_method="trn",
    spatial_key='spatial',
    key_added="align_spatial",
    mapping_key_added= "models_align",
    device="0",
)


cellbin_slices = []
for i, (binning_slice, raw_cellbin_slice) in enumerate(zip(binning_slices, cellbin_slices_raw)):
    cellbin_slice = raw_cellbin_slice.copy()
   
    if i == 0:
        cellbin_slice.obsm["align_spatial"] = cellbin_slice.obsm["spatial"]
    else:
        cellbin_slice = st.tl.alignment.transform.paste_transform(
            adata=cellbin_slice,
            adata_ref=binning_slice,
            spatial_key="spatial",
            key_added="align_spatial",
            mapping_key="models_align",
        )
        cellbin_slice.uns["models_align"] = binning_slice.uns["models_align"]
    cellbin_slices.append(cellbin_slice)


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

'''
space_multi(adata=binning_adata_raw, groupby=None, spatial_key="spatial", slice_key="slices", filename=os.path.join(save_binning_image_folder, "binning_alignment_all_raw.pdf"))

space_multi(adata=cellbin_adata_raw, groupby=None, spatial_key="spatial", slice_key="slices", filename=os.path.join(save_cellbin_image_folder, "cellbin_alignment_all_raw.pdf"))
'''
bin20_adata = st.tl.integrate(adatas=binning_slices, batch_key="slices")

cellbin_adata = st.tl.integrate(adatas=cellbin_slices, batch_key="slices")

bin20_adata.write_h5ad(filename=os.path.join(save_folder, f'{bin20_adata.obs["slices"][0][:-4]}_bin20.h5ad'), compression="gzip")

cellbin_adata.write_h5ad(filename=os.path.join(save_folder, f'{cellbin_adata.obs["slices"][0][:-4]}_cellbin.h5ad'), compression="gzip")
'''
space_multi(adata=bin20_adata, groupby=None, spatial_key="before_3d_align_spatial", slice_key="slices", filename=os.path.join(save_binning_image_folder, "binning_alignment_all.pdf"))

space_multi(adata=cellbin_adata, groupby=None, spatial_key="align_spatial", slice_key="slices", filename=os.path.join(save_cellbin_image_folder, "cellbin_alignment_all.pdf"))
'''
