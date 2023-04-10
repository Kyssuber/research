import os
import numpy as np
import glob
from astropy.table import Table
'''
AIM: append ncomp, group membership, and primary galaxy flag columns to my subsample catalog
-located at /mnt/astrophysics/muchogalfit-input-cats/v2_snrcoadd.fits
-I write this scipt in somewhat of a frenzy. no need for exotic formatting with function and class wizardry.
'''

#read table
path_to_table = '/mnt/astrophysics/muchogalfit-input-cats/'
table_name = 'v2_snrcoadd.fits'
path_to_galfit = '/mnt/astrophysics/muchogalfit-output/'
vf = Table.read(path_to_table+table_name).copy()   #creating a copy for testing purposes. we would not want to overwrite the original table with false information, would we? WOULD WE?

#set up empty arrays with which we will fashion our new columns.
mask = np.zeros(len(vf),dtype=bool)
groupGalaxy = np.zeros(len(vf),dtype=bool)
primaryGroup = np.zeros(len(vf),dtype=bool)
ncomp = np.zeros(len(vf),dtype=int)
column_names = ['group_flag','primaryGroup_flag','ncomp']

for n in range(len(vf)):
    if vf['sgacut_flag'][n]:  #if the galaxy is part of the subsample, then proceed
        if os.path.exists(path_to_galfit+VFID+'/galsFOV.txt'):
            fovtab = ascii.read(path_to_galfit+VFID+'/galsFOV.txt')
            group_vfids = self.fovtab['col1']    #the length of this textfile (i.e., number of entries)
            
            for num in range(len(fovtab)):
                if num==0:
                    ncomp[n] = len(fovtab)
                    primaryGroup[n] = True
                    groupGalaxy[n] = True
                    mask[n] = False
                else:
                    ncomp[vf['VFID']==group_vfids[num]] = len(fovtab)
                    primaryGroup[vf['VFID']==group_vfids[num]] = False
                    groupGalaxy[vf['VFID']==group_vfids[num]] = True
                    mask[vf['VFID']==group_vfids[num]] = False

        else:
            mask[n] = False   #galaxy is part of subsample. will not be masked out.
            #if the VFID is already part of a group, do NOT overwrite the flag with "False." This line is a fail-safe in case the primary galaxy containing this group galaxy runs first (in which case the groupGalaxy flag will be True here). If not, then the flag will be False by default (and will ultimately be overwritten with True above if the galaxy is, in fact, part of a group).
            groupGalaxy[n] = False if groupGalaxy[n] != True else True  
            primaryGroup[n] = False   #galaxy will not be a primary galaxy, as evidenced by the lack of a .txt file in the directory
            ncomp[n] = 1 if ncomp[n]==0 else ncomp[n]  #similar problem-o here. if ncomp[n] is already populated by a nonzero integer, lEaVe It AlOnE.
            
    else:
        groupGalaxy[n] = False   #placeholder value; will be masked out
        primaryGalaxy[n] = False   #placeholder value; will be masked out
        ncomp[n] = 1   #placeholder value; will be masked out
        mask[n] = True  #don't include a value at this index, since the galaxy is not part of the subsample

vf.write(path_to_table+'vf_snrcoadd_test.fits',overwrite=True)
