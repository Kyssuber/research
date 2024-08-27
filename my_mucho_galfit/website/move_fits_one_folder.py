'''
Aim: collect all relevant FITS files from /mnt/astrophysics/muchogalfit-output/ and port them to the all_input_fits located in /mnt/astrophysics/kconger_wisesize/. These files will be necessary when running the build_html_website.py script (in my_mucho_galfit) either locally (will have to scp the folder) or on the virtual machine.
'''

import os
import glob
import numpy as np
from astropy.table import Table
from astropy.io import fits

#convert cutouts from .fz to .fits, then save .fits to target_folder. not currently needed (all input files are already in the VFID output directories).
def fz_to_fits(path_to_im, galaxy_name, target_folder, group_name=None, group_flag=False):
    galaxy_w3 = galaxy_name+'-custom-image-W3.fits.fz'
    #galaxy_r = galaxy_name+'-custom-image-r.fits.fz'
    #galaxy_w1 = galaxy_name+'custom-image-W1.fits.fz'
    #galaxies = [galaxy_w3,galaxy_w1,galaxy_r]
    galaxiew=[galaxy_w3]
    if group_name is not None:
        group_w3 = group_name+'-custom-image-W3.fits.fz'
        #group_r = group_name+'-custom-image-r.fits.fz'
        #group_w1 = group_name+'-custom-image-W1.fits.fz'
        #groups = [group_w3,group_w1, group_r]
        groups = [group_w3]

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

#should give 2 images per galaxy -- one input cutout fits, and one output (out2.fits) fits.
def gather_w3_fits(catalog, cutouts_path, target_folder, fix=False):

    nonprimaries = []
    
    VFIDs = catalog['VFID']
    objnames = catalog['objname']
  
    #have to retrieve these images from the RA directories (LOLJK the input files are also in the output directories if galfit ran successfully!)
    for i in range(len(catalog)):
        
        print(VFIDs[i]+' input time')
        galaxy_folder = cutouts_path+VFIDs[i]+'/'
        try:
            input_w3 = glob.glob(galaxy_folder+'*-custom-image-W3.fits')

            for input_im in input_w3:  
                print('Moving '+input_im)
                os.system('cp '+input_im+' '+target_folder)
    
            print('Moving '+VFIDs[i]+' output file, if any.')
            output_mosaics_w3 = glob.glob(galaxy_folder+'*W3-out2*')   #currently only want convolution case!
            if fix:
                output_mosaics_w3 = glob.glob(galaxy_folder+'*W3*out2.fits')   #fixed BA, PA is *-W1-fixBA-out2.fits

            for im in output_mosaics_w3:   #if no images in output_mosaics, then none will be cp'd. if only one, then only one will be cp'd. usw.
                print(im)
                print('Moving '+im)
                os.system('cp '+im+' '+target_folder)
        except:
            nonprimaries.append(VFIDs[i])
            print(f'{VFIDs[i]} is a non-primary group galaxy.')

    print('list of non-primary group galaxies skipped:')
    print(nonprimaries)
#if galfit 'failed', then out* images will not appear in the folder. 

def grab_mask_images(catalog, host_folder_path, target_folder):
    VFIDs = catalog['VFID']
    objnames = catalog['objname']
    
    for i in range(len(catalog)):
        print('Moving '+VFIDs[i]+' mask(s), if any.')
        galaxy_folder = host_folder_path+VFIDs[i]+'/'
        #rmask = glob.glob(galaxy_folder+'*r-mask.fits')
        try:
            w3mask = glob.glob(galaxy_folder+'*wise-mask.fits')
            #masks = np.concatenate([rmask,w3mask])
            masks = w3mask
    
            for im in masks:   #if no images in masks, then none will be cp'd. if only one, then only one will be cp'd. usw.
                print('Moving '+im)
                os.system('cp '+im+' '+target_folder)
        except:
            print(f'{VFIDs[i]} is a non-primary group galaxy.')
            
def gather_w1_fits(catalog, host_folder_path_w1, target_folder, fix=False):
    
    dirnames = catalog['VFID']   #all directory names are simply the VFIDs
    for n in range(len(dirnames)):

      try:
          os.chdir(host_folder_path+dirnames[n])   #cd to correct directory
          print(dirnames[n])
          cutout_fits = glob.glob('*-custom-image-W1.fits')   #FITS cutout of galaxy
          for im in cutout_fits:
              print(im)
              os.system(f'cp {im} {target_folder}')
          #out1_fits = glob.glob('*-W1-out1.fits')   #unconvolved model parameters
          out2_fits = glob.glob('*-W1-out2.fits')   #convolved model parameters
          if fix:
            out2_fits = glob.glob('*-W1*out2.fits')   #fixed BA, PA is *-W1-fixBA-out2.fits
         
          #for imout1 in out1_fits:
          #  os.system(f'cp {imout1} {target_folder}')
          
          for imout2 in out2_fits:
              os.system(f'cp {imout2} {target_folder}')
      except:
          print(f'{dirnames[n]} is a non-primary group galaxy.')

def gather_r_fits(catalog, host_folder_path_r, target_folder):

    dirnames = catalog['VFID']   #all directory names are simply the VFIDs
    RAs = catalog["RA"]
    objname = catalog["group_name"]
    for n in range(len(dirnames)):
      try:
          os.chdir(host_folder_path+dirnames[n])   #cd to correct directory
          print(dirnames[n])
          cutout_fits = glob.glob('*-custom-image-r.fits')   #FITS cutout of galaxy

          for im in cutout_fits:
              print(im)
              os.system(f'cp {im} {target_folder}')
          
          #out1_fits = glob.glob('*-r-out1.fits')   #unconvolved model parameters
          out2_fits = glob.glob('*-r-out2.fits')   #convolved model parameters
         
          #for imout1 in out1_fits:
          #  os.system(f'cp {imout1} {target_folder}')
          
          for imout2 in out2_fits:
              os.system(f'cp {imout2} {target_folder}')
      except:
          print(f'{dirnames[n]} is a non-primary group galaxy or r-band cutout is not in the indicated directory.')

if __name__ == '__main__':
  
  homedir=os.getenv("HOME")
  vf = Table.read(homedir+'/VF_WISESIZE_photSNR.fits')   #contains objnames, RAs, and VFIDs; subsample_flag
  vf = vf[vf['subsample_flag']]
  host_folder_path = '/mnt/astrophysics/muchogalfit-output/'

  onefolder_path = '/mnt/astrophysics/kconger_wisesize/vf_html_w1_v3/all_input_fits/'
  
  print('Creating target directory '+onefolder_path)
  try:
    os.system('mkdir '+onefolder_path)
  except:
      print('Error: target directory already exists.')

  #gather_w3_fits(vf, host_folder_path, onefolder_path, fix=True)
  #grab_mask_images(vf, host_folder_path, onefolder_path)
  #gather_w1_fits(vf, host_folder_path, onefolder_path, fix=True)
  gather_r_fits(vf, host_folder_path, onefolder_path)
  
    
  #print('''
  #Move W3 images --> cutouts and out2 parameters...if fix=True, will ALSO pull fixed PA, BA results
  #gather_w3_fits(vf, host_folder_path, onefolder_path, fix=False)
  
  #Move W3 mask images (uncomment relevant lines for rband as well -- only need W3 masks if using W1 in place of rband)...
  #grab_mask_images(vf, host_folder_path, onefolder_path)
  
  #Move W1 images --> cutouts and out2 parameters...if fix=True, will ALSO pull fixed PA, BA results
  #gather_w1_fits(vf, host_folder_path, onefolder_path, fix=False)

  #Move r-band images --> cutouts and out2 parameters...
  #gather_r_fits(vf, host_folder_path, onefolder_path)

  #Note that vf, host_folder_path, and onefolder_path should already be defined. Type variables into terminal to print the results; if not desired, you may change these variables.
  #''')

  os.system('cd /mnt/astrophysics/kconger_wisesize/github/research/my_mucho_galfit/website/')
