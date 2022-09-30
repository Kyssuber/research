
import os
import numpy as np
from astropy.table import Table

homedir=os.getenv("HOME")
    
def write_file():    
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
