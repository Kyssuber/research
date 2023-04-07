#!/usr/bin/env python

"""
GOAL: 
* figure out which galaxies are in John's groups
* get list of groups

USAGE:
* run from muchogalfit-output directory



"""
import numpy as np
import os
from astropy.table import Table
import glob


topdir = '/mnt/astrophysics/rfinn/muchogalfit-output/'
try:
    os.chdir(topdir)
except FileNotFoundError: # assuming that we are running on virgo vms
    topdir = '/mnt/astrophysics/muchogalfit-output/'
    os.chdir(topdir)

# look through output directories
dirlist = []
t = open('Dirs.txt','r')
for line in t:
    dirlist.append(line.rstrip())
t.close()

# set up arrays to track group galaxies and central galaxies
groupflag = np.zeros(len(dirlist),'bool')
primaryflag = np.zeros(len(dirlist),'bool')

groupdirs_list = []
primaryvfids=[]
for i,dir in enumerate(dirlist):
    #print(i,os.getcwd(),dir)
    # search for galaxies with no fits files
    os.chdir(dir)
    fitsfiles = glob.glob("*.fits")

    if len(fitsfiles) > 0:
        os.chdir(topdir)
        continue
    elif len(fitsfiles) == 0:
        groupflag[i] = True

        # read in input file to get galname, objname, ra, dec, and bandpass
        sourcefile = open(dir+'sourcelist','r')
        galaxies = sourcefile.readlines()

        # parse information from file
        vfid, objname, ra, dec, bandpass = galaxies[0].rstrip().split()
        ra = float(ra)
        dec = float(dec)

        # check to see if this is the primary galaxy

        group_dir = f'/mnt/astrophysics/virgofilaments-data/{int(ra)}/{objname}_GROUP'
        #print(group_dir)
        if os.path.exists(group_dir):
            primaryflag[i] = True
            #print('group and found primary')
            groupdirs_list.append(group_dir)
            primaryvfids.append(dir)
            os.chdir(topdir)
        else:
            #print('group and but not primary')            
            os.chdir(topdir)
            continue

        
# write out list of groupdirs
outfile = open('groupDirs.txt','w')
for g in groupdirs_list:
    outfile.write(g+'\n')
outfile.close()

# write out list of groupdirs
outfile = open('groupPrimaryVFID.txt','w')
for p in primaryvfids:
    outfile.write(p+'\n')
outfile.close()

# save information about group
tab = Table([dirlist,groupflag,primaryflag],names=['VFID','groupflag','primary'])
tab.write('groupinfo.fits',format='fits',overwrite=True)
