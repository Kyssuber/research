'''
Aim: collect all relevant FITS files from /mnt/astrophysics/muchogalfit-output/ and port them to the all_input_fits located in /mnt/astrophysics/kconger_wisesize/. These files will be necessary when running the build_html_website.py script (in my_mucho_galfit) either locally (will have to scp the folder) or on the virtual machine.
'''

import os
import glob
import numpy as np
from astropy.table import Table
from astropy.io import fits

#convert cutouts from .fz to .fits, then save .fits to target_folder
def fz_to_fits(path_to_im, galaxy_name, target_folder, group_name=None, group_flag=False):
    galaxy_w3 = galaxy_name+'-custom-image-W3.fits.fz'
    galaxy_r = galaxy_name+'-custom-image-r.fits.fz'
    galaxies = [galaxy_w3,galaxy_r]
    if group_name is not None:
        group_w3 = group_name+'-custom-image-W3.fits.fz'
        group_r = group_name+'-custom-image-r.fits.fz'
        groups = [group_w3,group_r]
    for n in range(2):
        galaxy_path = path_to_im+galaxies[n]
        if group_flag==True:
            galaxy_path = path_to_im+groups[n]     
        try:
            fz = fits.getdata(galaxy_path)
            fz_header = (fits.open(galaxy_path)[1]).header
            print('Adding '+target_folder+galaxies[n][:-3])
            fits.writeto(target_folder+galaxies[n][:-3],fz,header=fz_header,overwrite=True)    #[:-3] removes the .fz from filename string
        except:
            print('Unable to pull galaxy cutout file.')

#should be two per galaxy - rband and w3
def grab_input_cutouts(catalog, sgaparams, input_cutouts_path, target_folder):
    VFIDs = catalog['VFID']
    VFID_V1s = catalog['VFID_V1']
    objnames = catalog['objname']
    RAs = catalog['RA_1']
    group_names = sgaparams['GROUP_NAME']
  
    #have to retrieve these images from the RA directories.
    for i in range(len(catalog)):
        
        print(VFIDs[i]+' input time')
        group_name = group_names[sgaparams['VFID']==VFID_V1s[i]]
        group_name = group_name[0]
        #print(group_name)

        ra_int = int(np.floor(RAs[i]))
        ra_folder = input_cutouts_path+str(ra_int)+'/'
        
        galaxy_folder = ra_folder+objnames[i]+'/'

        if os.path.isdir(galaxy_folder):
            print(galaxy_folder)
            fz_to_fits(galaxy_folder,objnames[i],target_folder)
        
        else:
            galaxy_folder_group = ra_folder+group_name+'/'
            print(galaxy_folder_group)
            fz_to_fits(galaxy_folder_group,group_name,target_folder,group_name,group_flag=True)

#should be four per galaxy - rband (nopsf, psf) and w3 (nopsf, psf)
#if galfit 'failed', then out* images will not appear in the folder. 
def grab_output_cutouts(catalog, host_folder_path, target_folder):
    VFIDs = catalog['VFID']
    objnames = catalog['objname']

    for i in range(len(catalog)):
        print('Moving '+VFIDs[i]+' output file, if any.')
        galaxy_folder = host_folder_path+VFIDs[i]+'/'
        output_mosaics_r = glob.glob(galaxy_folder+'*r-out*')
        output_mosaics_w3 = glob.glob(galaxy_folder+'*W3-out*')
        output_mosaics = np.concatenate([output_mosaics_r,output_mosaics_w3])

        for im in output_mosaics:   #if no images in output_mosaics, then none will be cp'd. if only one, then only one will be cp'd. usw.
            print('Moving '+im)
            os.system('cp '+im+' '+target_folder)

def grab_mask_images(catalog, host_folder_path, target_folder):
    VFIDs = catalog['VFID']
    objnames = catalog['objname']
    
    for i in range(len(catalog)):
        print('Moving '+VFIDs[i]+' mask(s), if any.')
        galaxy_folder = host_folder_path+VFIDs[i]+'/'
        rmask = glob.glob(galaxy_folder+'*r-mask.fits')
        w3mask = glob.glob(galaxy_folder+'*wise-mask.fits')
        masks = np.concatenate([rmask,w3mask])
        
        for im in masks:   #if no images in masks, then none will be cp'd. if only one, then only one will be cp'd. usw.
            print('Moving '+im)
            os.system('cp '+im+' '+target_folder)
    
if __name__ == '__main__':
  
  homedir=os.getenv("HOME")
  vf = Table.read(homedir+'/sgacut_coadd.fits')   #contains objnames, RAs, and VFIDs
  sga_params = Table.read(homedir+'/sga_params.fits')   #contains list of group names
  
  host_folder_path = '/mnt/astrophysics/muchogalfit-output/'
  input_cutouts_path = '/mnt/virgofilaments-data/'

  onefolder_path = '/mnt/astrophysics/kconger_wisesize/all_input_fits/'
  
  print('Creating target directory '+onefolder_path)
  os.system('mkdir '+onefolder_path)
  
  print('Moving postage stamp cutouts for rband and W3...')
  grab_input_cutouts(vf, sga_params, input_cutouts_path, onefolder_path)
  print('Moving GALFIT output mosaics for rband and w3...')
  grab_output_cutouts(vf, host_folder_path, onefolder_path)
  print('Moving r-band and W3 mask images...')
  grab_mask_images(vf, host_folder_path, onefolder_path)

  
