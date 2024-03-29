
import os
import numpy as np
import argparse
from astropy.io import fits
from astropy.table import Table
from astropy.io import ascii
from matplotlib import pyplot as plt
import glob


homedir = os.getenv('HOME')

class galaxy():
    def __init__(self,name='galname',vfid='VFID',band='3'):
        self.band = band
        self.galname = name
        self.vfid = vfid
        self.image_rootname = self.galname+'-unwise-w'+str(self.band)
        self.image = self.image_rootname+'-img-m.fits'

    #MASK
    def mask_gui(self,output_dir):
        os.chdir('/mnt/astrophysics/wisesize/'+str(self.vfid))
        im = glob.glob('*w3-img-m.fits')[0]
        print(im)
        #the --auto arg permits the running of maskwrapper.py without necessitating the manual closing of the gui to proceed to the next galaxy
        string = "python ~/github/halphagui/maskwrapper.py --image {}".format(im)+" --auto"
        os.system(string)
        im_mask = glob.glob('*mask.fits')[0]
        im_mask_inv = glob.glob('*mask.fits')[1]
        self.mask_image = im_mask
        #copy mask image to output directory 'output_dir'
        os.system('cp '+ im_mask + ' ' + str(output_dir))
        os.system('cp '+ im_mask_inv + ' ' + str(output_dir))

    def mask_gui_one(self,output_dir):
        os.chdir('/mnt/astrophysics/wisesize/'+str(self.vfid))
        im = glob.glob('*w3-img-m.fits')[0]
        print(im)
        string = "python ~/github/halphagui/maskwrapper.py --image {}".format(im)
        os.system(string)
        im_mask = glob.glob('*mask.fits')[0]
        im_mask_inv = glob.glob('*mask.fits')[1]
        self.mask_image = im_mask
        os.system('cp '+ im_mask + ' ' + str(output_dir))
        os.system('cp '+ im_mask_inv + ' ' + str(output_dir))

    
    def mask_haim(self,output_dir):
        os.chdir('/mnt/astrophysics/wisesize/'+str(self.vfid))
        im = glob.glob('*w3-img-m.fits')[0]
        print(im)
        get_ipython().run_line_magic('run','/mnt/astrophysics/kconger_wisesize/github/HalphaImaging/python3/uat_mask.py --R '+ im +' --nods9')
        im_mask = glob.glob('*mask.fits')[0]
        im_mask_inv = glob.glob('*inv-mask.fits')[0]
        self.mask_image = im_mask
        #copy mask image to output directory 'output_dir'
        os.system('cp '+ im_mask + ' ' + str(output_dir))
        os.system('cp '+ im_mask_inv + ' ' + str(output_dir))

        #loop through galaxies in input fits file, generate masks and cp results to '/mnt/astrophysics/kconger_wisesize/github/masks' directory

def run_mask_one(galaxy_sample,type):
    output_dir = '/mnt/astrophysics/kconger_wisesize/github/masks'
    vfid = galaxy_sample['VFID']
    name = galaxy_sample['prefix']
    g = galaxy(name=name,vfid=vfid,band='3')
    print(name)
    if str(type) == 'haim':
        g.mask_haim(output_dir == output_dir)
    if str(type) == 'gui':
        g.mask_gui_one(output_dir == output_dir)
    else:
        print('types: gui (for halphagui) or haim (for HalphaImaging)')

    
    
def run_mask_all(galaxy_sample,type):
    output_dir = '/mnt/astrophysics/kconger_wisesize/github/masks'
    for n in range(0,len(galaxy_sample)):
        vfid = galaxy_sample['VFID'][n]
        name = galaxy_sample['prefix'][n]
        g = galaxy(name=name,vfid=vfid,band='3')
        print(name)
        if str(type) == 'haim':
            g.mask_haim(output_dir == output_dir)
        if str(type) == 'gui':
            g.mask_gui(output_dir == output_dir)
        else:
            print('types: gui (for halphagui) or haim (for HalphaImaging)')



if __name__ == '__main__':
    homedir = os.getenv("HOME")
    vf = Table.read(homedir+'/vfcut.fits',format='ascii')
    print('run_mask_one(galaxy,type), run_mask_all(galaxy_sample,type)')
    print('types: haim for HalphaImaging, gui for halphagui')
    
