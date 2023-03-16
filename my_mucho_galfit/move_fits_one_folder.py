'''
Aim: collect all relevant FITS files from /mnt/astrophysics/muchogalfit-output/ and port them to the fits_folder located in the home directory. These files will be necessary when running the build_html_website.py script (in my_mucho_galfit) either locally (will have to scp the folder) or on the virtual machine.
'''

import os
import sys
import glob
import numpy as np
from astropy.table import Table

#should be two per galaxy - rband and w3
def grab_input_cutouts(catalog, host_folder_path, target_folder):
  VFIDs = catalog['VFID']
  objnames = catalog['objname']
  
  for i in range(len(catalog)):
    galaxy_folder = host_folder_path+VFIDs[i]+'/'
    input_cutouts = glob.glob(galaxy_folder+objnames[i]+'-custom-image-*')
    
    for im in input_cutouts:
      os.system('cp '+im+' '+target_folder)
  
  #should be four per galaxy - rband (nopsf, psf) and w3 (nopsf, psf)
  def grab_output_cutouts(catalog, host_folder_path, target_folder):
    VFIDs = catalog['VFID']
    objnames = catalog['objname']

    for i in range(len(catalog)):
      galaxy_folder = host_folder_path+VFIDs[i]+'/'
      output_mosaics = glob.glob(galaxy_folder+objnames[i]+'*out*.fits')

      for im in output_mosaics:
        os.system('cp '+im+' '+target_folder)

if __name__ == '__main__':
  homedir=os.getenv("HOME")
  vf = Table.read(homedir+'/sgacut_coadd.fits')
  onefolder_path = homedir+'/all_input_fits'
  host_folder_path = '/mnt/astrophysics/muchogalfit-output/'

  'Moving postage stamp cutouts for rband and W3...'
  grab_input_cutouts(vf, host_folder_path, onefolder_path)
  'Moving GALFIT output mosaics for rband and w3...'
  grab_output_cutouts(vf, host_folder_path, onefolder_path)
