#before PS
import spateo as st
import matplotlib.pyplot as plt
import sys
import os
from PIL import Image

project = sys.argv[1]
path_bin1 = sys.argv[2]
outdir = sys.argv[3]

if not os.path.exists(outdir):
    os.makedirs(outdir)

path_img_bin1 = outdir+'/'+project+'_bin1_beforePS.tiff'
print(path_img_bin1)

adata = st.io.read_bgi_agg(path_bin1)
img = adata.X.todense()
img[img>=10]=10
img = img/10*255
img = Image.fromarray(img)
img = img.convert("L")

img.save(path_img_bin1)
