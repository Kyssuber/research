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


def run_galfit(cat,gal_output_path,directory):
    
    '''
    GOAL:
        run GALFIT on all sample galaxies.
        note --> all VFID repositores are labeled with VFID_V1 names
        
    --- first iteration of this function included two optional parameters, convflag and sersic_start, which are now embedded in a different script (make_inputcomp.py)
        
    --- make_inputcomp.py has already cp'd .psf files into gal_output directory   
    '''
    
    for i in range(0,len(cat)):
        #need all files to be in directory in which galfit input script is run
        #I could write the full pathname in the input scripts (via make_inputcomp.py), but...I did not.
        print(cat['VFID'][i])
        os.chdir(gal_output_path)
        input_script = glob.glob(cat['VFID'][i]+'*'+'galfit.input.'+'*')[0]
        psf_img = glob.glob(cat['VFID'][i]+'*'+'psf.fits')[0]
        
        os.system('cp '+ input_script + ' '+directory+cat['VFID_V1'][i])
        os.system('cp '+ psf_img + ' '+directory+cat['VFID_V1'][i])
        
        #
        os.chdir(directory+cat['VFID_V1'][i])
        
        s = 'galfit '+input_script
        print('run the following: ',s)
        os.system(s)
        
        #tidy the directory a bit
        os.system('rm '+psf_img)   #remove psf image
        os.system('rm galfit.*')   #remove all galfit iteration logs
        

def organize_output(cat,gal_output_path,directory):
    
    '''
    GOAL:
        transfer galfit output fits files into gal_output directory
        note --> all VFID repositories are labeled with VFID_V1 names
    '''

    for i in cat['VFID_V1']:
        os.chdir(directory+str(i))
        try:
            output_script = glob.glob('*'+'out.fits')[0]
            os.system('cp '+output_script+' '+gal_output_path)
        except:
            print(i+' failed to cp due to GALFIT error.')
        
if __name__ == '__main__':
  
    #run scripts external to ipython environment!
    #helpful particularly when running with psf --> can use nohup to pipe output into a .txt file, then can terminate the window without terminating the program
    
    if '-h' in sys.argv or '--help' in sys.argv:
        print("Usage: %s [-param_file (name of parameter file, no single or double quotation marks)] [-range_min integer] [-range_max integer]" % sys.argv[0])
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
        try:
            vf = sgacut[range_min:range_max]
        except:
            print('enter correct min and max ranges!')
  
       
    #create dictionary with keywords and values, from parameter .txt file
    
    homedir = os.getenv("HOME")
    
    param_dict = {}
    with open(homedir+'/'+param_file) as f:
        for line in f:
            try:
                key = line.split()[0]
                val = line.split()[1]
                param_dict[key] = val
            except:
                continue
    
    ############
    
    cat_path = param_dict['vf_sample']
    vf = Table.read(cat_path)
    
    dummycat_path = param_dict['dummycat']
    dummycat = Table.read(dummycat_path)
    
    directory = param_dict['directory']
    gal_output_path = param_dict['gal_output_path']
    
    organize = int(param_dict['organize'])
    
    ############
    
    run_galfit(vf,gal_output_path=gal_output_path,directory=directory)
    
    if organize==1:
        organize_output(vf,gal_output_path=gal_output_path,directory=directory)
    
    
