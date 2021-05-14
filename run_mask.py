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
    def mask(self,output_dir):
        os.chdir('/mnt/astrophysics/wisesize/'+str(self.vfid))
        im = glob.glob('*w3-img-m.fits')[0]
        get_ipython().run_line_magic('run','~/github/halphagui/maskwrapper.py'+' --image '+im)
        im_mask = glob.glob('*mask.fits')[0]
        im_mask_inv = glob.glob('*mask.fits')[1]
        self.mask_image = im_mask
        #copy mask image to output directory 'output_dir'
        os.system('cp '+ im_mask + ' ' + str(output_dir))
        os.system('cp '+ im_mask_inv + ' ' + str(output_dir))


#loop through galaxies in input fits file, generate masks and cp results to '~/github/masks' directory
def run_mask_all(galaxy_sample):
    output_dir = homedir+'/github/masks'
    for n in range(0,len(galaxy_sample)):
        vfid = galaxy_sample['VFID'][n]
        name = galaxy_sample['prefix'][n]
        g = galaxy(name=name,vfid=vfid,band='3')
        print(name)
        g.mask(output_dir=output_dir)
