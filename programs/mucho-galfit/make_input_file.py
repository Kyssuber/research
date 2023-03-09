
'''
Goal: generate .txt input file for running galfit in parallel on all sample galaxies.
-accepts vf subsample. already defined as a variable below, so only need to call write_file() to run.
'''

import os
import numpy as np
from astropy.table import Table

homedir=os.getenv("HOME")
    
def write_file():
    #table contains all v2 galaxies; 'sgacut_flag' is the flag column which trims vf such that only galaxies are those in sga_coadd subsample
    vf = Table.read(homedir+'/v2_snrcoadd.fits')
    vf = vf[vf['sgacut_flag']]
    file_obj = open('galaxy_file.txt','w')
    for i in range(0,len(vf)):
        ra = vf['RA_1'][i]
        ra_int = int(np.floor(ra))
        objname = vf['objname'][i]
        file_obj.write(str(ra_int)+'/'+objname)
        file_obj.write('\n')
    file_obj.close()
