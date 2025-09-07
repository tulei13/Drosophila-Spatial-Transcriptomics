#cell_segmentation

import os,re,sys,getopt
import cv2
import numpy as np
import warnings
from pathlib import Path
from typing import Literal, Optional, Union,Tuple
import spateo as st
from anndata import AnnData
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
import pandas as pd
import anndata as ad
import torch
from skimage.filters import unsharp_mask
from matplotlib.backends.backend_pdf import PdfPages

#
def_path = '/hwfssz5/ST_SUPERCELLS/P21Z10200N0206/public/drosophila_pipe_v1/'
sys.path.append(def_path)
from supplement import crop_by_ssdna,bin1tobinx,custom_segmentation,cellpose,distance_to_center
from supplement import qc_seg,qc_cellbin,qc_bin20,equalhist_transfer,output_img
from supplement import check_nan_gene
from supplement import deHaze, output_img, segmentation_spateo

file = '/hwfssz5/ST_SUPERCELLS/P21Z10200N0206/jiazihan/L3-early/L3-early.csv'
seg_para = pd.read_csv(file)

stage = 'L3-early'

indir = '/hwfssz5/ST_SUPERCELLS/P21Z10200N0206/project/drosophila/03.ST/03.basic_anlysis/second_batch'
outdir = '/hwfssz5/ST_SUPERCELLS/P21Z10200N0206/jiazihan/L3-early/'
aout = outdir + '/' + stage
#makedirs
name_list = ['seg_data','qc']
for a in ['seg_01']:
  out = os.path.join(aout,a)
  if not os.path.exists(out):
      os.makedirs(out)
  for i in name_list:
      outdir = os.path.join(out,i)
      if not os.path.exists(outdir):
          os.makedirs(outdir)
            
out1 = os.path.join(aout,'seg_01')
qc_txt_save_01 = os.path.join(out1,'qc','all_slices_qc.txt')
seg_out_qc = os.path.join(out1,'qc')

for i in range(len(seg_para)):
  new_bin1_save = indir + '/' + seg_para['stage'][i]  + '/01_segmentation/bin1_data/' + seg_para['stage'][i] + '_' + seg_para['slice'][i] + '_SS200000174BL_bin1_final.gem.gz'
  cropped_img_save = indir + '/' + seg_para['stage'][i]  + '/01_segmentation/cropped_img/' + seg_para['stage'][i] + '_' + seg_para['slice'][i]+'_SS200000174BL_cropped_img.tif'
  if __name__ == '__main__':   
    project = re.split('/',new_bin1_save)[-1][:-18]
    seg_data_save = os.path.join(out1,project+'_seg.h5ad')
    seg_gem_save_01 = os.path.join(out1,project+'_seg.gem.gz')
		# 03_segmentation_xr
    seg_adata = custom_segmentation(
        path_bin1 = new_bin1_save,
        path_ssdna = cropped_img_save,
        dilation_distance = seg_para['dilation_distance'][i],
        amount = seg_para['amount'][i],
        block_size = seg_para['block_size'][i],
        radius = seg_para['radius'][i],
        otsu_class = seg_para['otsu_class'][i],
        otsu_index = seg_para['otsu_index'][i]
        )
    st.cs.expand_labels(seg_adata, "cellpose_labels", distance = 5, max_area=np.inf, out_layer = "cell_labels_expanded")
    print('===== ',project,'cellbin segmentation finished! =====')
		# qc_seg
    qc_seg(project = project, seg_adata = seg_adata, cropped_img = cropped_img_save, save = seg_out_qc)
		#####################################################
		###################### cellbin ######################
    #####################################################
    cell_layer = "cell_labels_expanded"
    # cell labels
    cells = pd.DataFrame(seg_adata.layers[cell_layer])
    cells["x"] = cells.index
    cells = pd.melt(cells, id_vars=["x"])
    cells = cells[cells["value"] != 0]
    cells.columns = ["x", "y", cell_layer]
    cells.sort_values(by=[cell_layer, "x", "y"], inplace=True)
    cells.to_csv(
        seg_gem_save_01,
        sep="\t",
        index=False,
        )    
    adata = st.io.read_bgi(
        path=new_bin1_save,
        segmentation_adata=seg_adata,
        labels_layer='cell_labels_expanded',
        seg_binsize=1,
        )
    adata.obs["slices"] = str(project)
		#adata.uns['raw_min'] = {project+'_x':a,project+'_y':b}
    adata.write(seg_data_save, compression="gzip")
		# qc_cellbin
    qc_cellbin(
        project = project,
        adata = adata,
        save = seg_out_qc,
        )
    print('===== seg by xr :',project,' cellbin finished! =====')
    print('===== all slices segmentation finished ! =====')

