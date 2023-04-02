import os
import glob
from astropy.table import Table

#check whether galaxies in the mucho-galfit output directories on the virtual machine contain mask fits images for r-band and w3

vf_cat_path = '/mnt/astrophysics/muchogalfit-input-cats/v2_snrcoadd.fits'
vf_cat = Table.read(vf_cat_path)

subsample_flag = vf_cat['sgacut_flag']
vf_subsample = vf_cat[subsample_flag]
objnames = vf_subsample['objnames']

path_to_dir = '/mnt/astrophysics/muchogalfit-output/'    

mask_w3_suffix = '-custom-image-wise-mask.fits'
mask_r_suffix = '-custom-image-r-mask.fits'

group_count=0
no_mask_list_r=[]
no_mask_list_w3=[]

for galname in objnames:
    galpath = path_to_dir+galname+'/'
    if os.path.exists(galpath):
        rmask = galpath+galname+mask_r_suffix
        w3mask = galpath+galname+mask_w3_suffix
        if not os.path.exists(rmask):
            no_mask_list_r.append(galname)
        if not os.path.exists(w3mask):
            no_mask_list_w3.append(galname)  
    else:
        group_count+=1
    
print(f'# group galaxies: {group_count}')
print(f'galaxies without rmask: {no_mask_list_r}')
print(f'galaxies without w3mask: {no_mask_list_w3}')
    