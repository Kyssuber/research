'''
GOAL:
- run galfit on all scripts generated using inputcomp_notrim.py
OUTPUT:
- .FITS output files containing mosaic and output parameters (via header information); output_params.py and output_mosaic.py will then attend to these clunky buggers.
'''

import numpy as np
import warnings
from astropy.io import fits
import os
import sys
import argparse
from astropy import units as u
from astropy.wcs import WCS
import astropy.wcs as wcs
from astropy.table import Table
from astropy.io import ascii
import glob

homedir = os.getenv("HOME")
os.sys.path.append(homedir+'/github/virgowise/')
#import rungalfit as rg #has needed functions already defined


def run_galfit(cat):
    
    '''
    
    GOAL:
        run GALFIT on all sample galaxies.
        note --> all VFID repositores are labeled with VFID_V1 names
        
        first iteration of this function included two optional parameters, convflag and sersic_start, which are now embedded in a different script (inputcomp_notrim.py)
       
    '''
    
    for i in range(0,len(cat)):
        
        print(cat['VFID'][i])
        os.chdir('/mnt/astrophysics/kconger_wisesize/github/gal_output')
        input_script = glob.glob(cat['VFID'][i]+'*'+'galfit.input.'+'*')[0]
        psf_img = glob.glob(cat['VFID'][i]+'*'+'psf.fits')[0]
        
        os.system('cp '+ input_script + ' /mnt/astrophysics/wisesize/'+cat['VFID_V1'][i])
        os.system('cp '+ psf_img + ' /mnt/astrophysics/wisesize/'+cat['VFID_V1'][i])
        
        os.chdir('/mnt/astrophysics/wisesize/'+cat['VFID_V1'][i])
        
        s = 'galfit '+input_script
        print('run the following: ',s)
        os.system(s)
        

def organize_output(cat):
    
    '''
    
    GOAL:
        transfer galfit output fits files into gal_output directory
        note --> all VFID repositories are labeled with VFID_V1 names
    '''

    for i in cat['VFID_V1']:
        os.chdir('/mnt/astrophysics/wisesize/'+str(i))
        try:
            output_script = glob.glob('*'+'out.fits')[0]
            os.system('cp '+output_script+' /mnt/astrophysics/kconger_wisesize/github/gal_output')
        except:
            print(i+' failed to cp due to GALFIT error.')
        
if __name__ == '__main__':
  

    homedir = os.getenv("HOME")

    #catalogs
    dummycat = Table.read(homedir+'/dummycat.fits')
    sgacut = Table.read(homedir+'/sgacut_coadd.fits')

    
    #run scripts external to ipython environment!
    #helpful particularly when running with psf --> can use nohup to pipe output into a .txt file, then can terminate the window without terminating the program
    
    if '-h' in sys.argv or '--help' in sys.argv:
        print("Usage: %s [-organize (0 for False, 1 for True; will cp all output fits files into gal_output directory)] [-cat sga (currently sga only)] [-range_min integer] [-range_max integer]" % sys.argv[0])
        sys.exit(1)
    
    if '-organize' in sys.argv:
        p = sys.argv.index('-organize')
        organize = int(sys.argv[p+1])
        if (organize != 0) & (organize != 1):
            print('organize=0 or organize=1 only')
            sys.exit(1)
        if organize == 1:
            organize = True
        if organize == 0:
            organize = False
    
    if '-cat' in sys.argv:
        p = sys.argv.index('-cat')
        cat = sys.argv[p+1]
        if cat == 'sga':
            vf = sgacut
    
    if '-range_min' in sys.argv:
        p = sys.argv.index('-range_min')
        range_min = int(sys.argv[p+1])

    if '-range_max' in sys.argv:
        p = sys.argv.index('-range_max')
        range_max = int(sys.argv[p+1])
        try:
            vf = sgacut[range_min:range_max]
            print('d')
        except:
            print('enter correct min and max ranges!')
            
    run_galfit(vf)
    if organize==True:
        organize_output(vf)
    
    
