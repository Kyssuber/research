'''
Aim: collect all relevant FITS files from /mnt/astrophysics/muchogalfit-output/ and port them to the all_input_fits located in /mnt/astrophysics/kconger_wisesize/. These files will be necessary when running the build_html_website.py script (in my_mucho_galfit) either locally (will have to scp the folder) or on the virtual machine.
'''

import os
import sys
import glob
import numpy as np
from astropy.table import Table

#convert cutouts from .fz to .fits
def fz_to_fits(path_to_im, galaxy_name):
    galaxy_w3 = path_to_im+galaxy_name+'-custom-image-W3.fits.fz'
    galaxy_r = path_to_im+galaxy_name+'-custom-image-r.fits.fz'
    galaxies = [galaxy_w3,galaxy_r]
    for galaxy in galaxies:
        fz = fits.getdata(galaxy)
        fz_header = (fits.open(galaxy)[1]).header
        fits.writeto(galaxy[:-3],fz,header=fz_header,overwrite=True)    #[:-3] removes the .fz from location string

#should be two per galaxy - rband and w3
def grab_input_cutouts(catalog, input_cutouts_path, target_folder):
    VFIDs = catalog['VFID']
    objnames = catalog['objname']
    RAs = catalog['RA_1']
  
    #have to retrieve these images from the RA directories.
    for i in range(len(catalog)):
        print(VFIDs[i])
        ra_int = int(np.floor(RAs[i]))

        ra_folder = input_cutouts_path+str(ra_int)+'/'
        galaxy_folder = ra_folder+objnames[i]+'/'
        galaxy_folder_group = ra_folder+objnames[i]+'_GROUP/'

        print(galaxy_folder)

        if os.path.isdir(galaxy_folder):
            fz_to_fits(galaxy_folder,objnames[i])
            input_cutouts = glob.glob(galaxy_folder+objnames[i]+'-custom-image-*.fits')
        else:
            fz_to_fits(galaxy_folder_group,objnames[i])
            input_cutouts = glob.glob(galaxy_folder_group+objnames[i]+'-custom-image-*.fits')

        for im in input_cutouts:
            print('Moving '+im)
            os.system('cp '+im+' '+target_folder)

#should be four per galaxy - rband (nopsf, psf) and w3 (nopsf, psf)
#if galfit 'failed', then out* images will not appear in the folder. 
def grab_output_cutouts(catalog, host_folder_path, target_folder):
    VFIDs = catalog['VFID']
    objnames = catalog['objname']

    for i in range(len(catalog)):
        print('Moving '+VFIDs[i])
        galaxy_folder = host_folder_path+VFIDs[i]+'/'
        output_mosaics = glob.glob(galaxy_folder+objnames[i]+'*out*.fits')

        for im in output_mosaics:
            print('Moving '+im)
            os.system('cp '+im+' '+target_folder)

if __name__ == '__main__':
  
  homedir=os.getenv("HOME")
  vf = Table.read(homedir+'/sgacut_coadd.fits')
  
  host_folder_path = '/mnt/astrophysics/muchogalfit-output/'
  input_cutouts_path = '/mnt/virgofilaments-data/'

  onefolder_path = '/mnt/astrophysics/kconger_wisesize/all_input_fits/'
  
  print('Creating target directory '+onefolder_path)
  os.system('mkdir '+onefolder_path)
  
  'Moving postage stamp cutouts for rband and W3...'
  grab_input_cutouts(vf, input_cutouts_path, onefolder_path)
  'Moving GALFIT output mosaics for rband and w3...'
  grab_output_cutouts(vf, host_folder_path, onefolder_path)
