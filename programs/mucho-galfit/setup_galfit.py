#!/usr/bin/env python
'''

call as:
setup_galfit.py W3
where W3 can be any wavelength band [g,r,z,W1 through W4]
-will create repositories for sample galaxies, each containing a sourcelist with which to run galfit

'''
import os
import sys

wavelength = sys.argv[1]

from astropy.table import Table

# filepath which contains VF tables 
filepath = '/mnt/astrophysics/muchogalfit-input-cats/'


# read in the main table
main_name = 'vf_v2_main.fits'
main_table = Table.read(filepath+main_name)


# read in table that is line-matched to v2 main table that says galaxy is in sample
line_matched_name = 'v2_snrcoadd.fits'
matched_table = Table.read(filepath+line_matched_name)


# cut main catalog based on sample flag
sampleflag = matched_table['sgacut_flag']
finalcat = matched_table[sampleflag]

# for each galaxy, create a directory and write sourcelist
for i in range(len(finalcat)):
    
    galpath = '/mnt/astrophysics/muchogalfit-output/'+finalcat['VFID'][i]
    os.mkdir(galpath)

    sourcefile = galpath+'/{}sourcelist'.format(finalcat['VFID'][i])
    sourcelist = open(sourcefile,'w')
    
    # write out one line with VFID, objname, RA, DEC, wavelength
    output_string = finalcat['VFID'][i] + ' ' + finalcat['objname'][i] + ' ' + str(finalcat['RA_1'][i]) + ' ' + str(finalcat['DEC_1'][i]) + ' ' + str(wavelength) + ' \n'.format()
    sourcelist.write(output_string)
    sourcelist.close()
