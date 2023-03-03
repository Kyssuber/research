'''
GOAL:
- run galfit on all scripts generated using make_inputcomp.py
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

def run_galfit(cat,gal_output_path):
    '''run GALFIT on all sample galaxies.
        note --> the first iteration of this function included two optional parameters, convflag and sersic_start, which are now embedded in a different script (make_inputcomp.py)
    '''
    
    for i in range(0,len(cat)):
        
        print(cat['VFID'][i])
        try:
            os.chdir(gal_output_path)
        except:
            os.chdir('/mnt/astrophysics/kconger_wisesize/github/gal_output/')
        
        s = 'galfit '+input_script
        print('Executing: '+s)
        os.system(s)
        
if __name__ == '__main__':
    
    #run scripts external to ipython environment!
    #helpful particularly when running with psf --> can use nohup to pipe output into a .txt file, then can terminate the window without terminating the program.
    
    homedir = os.getenv("HOME")
    
    if '-h' in sys.argv or '--help' in sys.argv:
        print("Usage: %s [-param_file /path/to/paramfile.txt] [-range_min integer] [-range_max integer]" % sys.argv[0])
        sys.exit(1)
    
    if '-param_file' in sys.argv:
        p = sys.argv.index('-param_file')
        param_file = str(sys.argv[p+1])
    
    if '-range_min' in sys.argv:
        p = sys.argv.index('-range_min')
        range_min = int(sys.argv[p+1])

    if '-range_max' in sys.argv:
        p = sys.argv.index('-range_max')
        range_max = int(sys.argv[p+1])
    
    #create dictionary with keywords and values, from parameter .txt file

    param_dict = {}
    with open(homedir+'/'+param_file) as f:
        for line in f:
            try:
                key = line.split()[0]
                val = line.split()[1]
                param_dict[key] = val
            except:
                continue    
    
    #now...extract parameters and assign to relevantly-named variables
    cat_path = param_dict['vf_catalog']
    cat = Table.read(cat_path)
    dummycat_path = param_dict['dummycat']
    dummycat = Table.read(dummycat_path)
    gal_output_path = param_dict['gal_output_path']   
    
    try:
       cat=cat[range_min:range_max]
    except:
        print('if you intended to "cut the cat," please enter correct min and max ranges!')
    
    run_galfit(cat,gal_output_path)
    
    
