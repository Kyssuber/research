'''
GOAL:
- generate .html files that comprise a website interface for viewing vf galaxy images and models
INPUT:
- 
OUTPUT:
- 
General instructions for opening .html locally:
- 
'''

import numpy as np
from matplotlib import pyplot as plt
import os
from astropy.table import Table
from astropy.io import ascii
import glob


#NOTE: once psf is properly established, may no longer be necessary to include no_psf row. Until then, do include, and be sure both are the correct length (which will not be so if GALFIT fails with ncomp>1 for psf and not non_psf, or vice versa)


def build_gal_html(i, vf_sample, galpath, htmlpath, galfit_params_nopsf, galfit_params_psf, dummycat = None, fixed_galfit_params_nopsf = None, fixed_galfit_params_psf = None, ):
    
    #VFID (v2) html name of webpage
    htmlname = str(vf_sample[i]['VFID']) + '.html'
    
    #full local path to where this .html is stored
    full_htmlpath = htmlpath + htmlname
    
    #begin to write text in full_htmlpath file
    with open(full_htmlpath, 'w') as html:
        
        #I think much of this is style preparation, similar to the initial lines in a LaTeX document
        html.write('<html><body>\n')
        html.write('<title>' + str(vf_sample['prefix'][i]) + '<\title>\n')  #title that appears on the browser tab
        html.write('<style type="text/css">\n')
        html.write('.img-container{text-align: left;}')
        html.write('table, td, th {padding: 5px; text-align: center; border: 1px solid black;}\n')
        html.write('p {display: inline-block;;}\n')
        html.write('</style>\n')
        
        #name of central galaxy
        html.write('<fontsize="40">Central Galaxy: ' + str(vf_sample['prefix'][i]) + '</font><br />\n')
        
        #add hyperlinks that return user to either homepage, next galaxy (if i != 0), or prev galaxy (if i != len(sample) - 1)
        html.write('<a href=main.html>Return to Homepage</a></br />\n')
        if i != len(vf_sample) - 1:
            html.write('<a href='+str(vf_sample['VFID'][i+1])+'.html>Next Galaxy</a></br />\n')
        if i != 0:
            html.write('<a href='+str(vf_sample['VFID'][i-1])+'.html>Previous Galaxy</a></br />\n')
        
        '''
        unclear how we will treat postage stamps with 2+ Sersic obj, so this and following relevant sections will be commented out until further notice.
        
        #if zero, then ncomp = 1
        ncomp = len(np.where(dummycat['central galaxy'] == vf_sample['VFID'][i])[0]) + 1
        '''
        
        ncomp = 1
        
    



def build_html_all(sample,paba_comparison=False):
    for i in range(0,len(sample)):
        build_html_one(sample,i,ku_or_siena,paba_comparison)
        print(sample['VFID'][i])
    build_htmlhome_galfit(sample,ku_or_siena)




if __name__ == '__main__':    
    
    if '-h' in sys.argv or '--help' in sys.argv:
        print("Usage: %s [-param_file (name of parameter file, no single or double quotation marks)]")
        sys.exit(1)
    
    if '-param_file' in sys.argv:
        p = sys.argv.index('-param_file')
        param_file = str(sys.argv[p+1])
    
    
    homedir = os.getenv("HOME")
       
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
    
    #catalogs and pathnames
    cat_path = param_dict['vf_catalog']
    cat = Table.read(cat_path)
    dummycat_path = param_dict['dummycat']
    dummycat = Table.read(dummycat_path)
    directory = param_dict['directory']
    psf_filepath = param_dict['psf_filepath']
    psf_filepath = homedir+psf_filepath
    gal_output_path = param_dict['gal_output_path']
