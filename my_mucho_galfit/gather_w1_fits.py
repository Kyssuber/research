'''
AIM: cp the W1 galaxy cutouts and out1, out2 GALFIT parameters from each individual directory in the virtual machine into a separate directory, which I will then scp to my local machine. This amounts to ~496*3 FITS files. 

I do not intend to be fancy schmancy here.
'''

import glob
import numpy as np
import os
from astropy.table import Table
homedir = os.getenv("HOME")

path_to_dir = '/mnt/astrophysics/muchogalfit-output/'
cutout_ext = '-custom-image-W1.fits'   #generic image extension for the images
destination_dir = '/mnt/astrophysics/muchogalfit-output/W1_fits/'

vf = Table.read(homedir+'/vf_snrcoadd.fits')
vf_cut = vf[vf['subsample_flag']]  #reduces number of rows to 496, which corresponds with the number of subsample galaxies (pre-GALFIT)

dirnames = vf_cut['VFID']   #all directory names are simply the VFIDs

for n in range(len(dirnames)):
  os.chdir(path_to_dir+dirnames[n])   #cd to correct directory
  
  cutout_fits = glob.glob('*'+cutout_ext)[0]
  out1_fits = glob.glob('*-W1-out1.fits')[0]   #unconvolved model parameters
  out2_fits = glob.glob('*-W1-out2.fits')[0]   #convolved model parameters

  os.system('cp cutout_fits out1_fits out2_fits destination_dir')

  print(f'Galaxy {} completed.')
