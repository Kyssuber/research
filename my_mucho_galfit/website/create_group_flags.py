import os
import numpy as np
import glob
from astropy.table import Table
from astropy.io import ascii
from astropy.table import MaskedColumn

'''
AIM: append ncomp, group membership, and primary galaxy flag columns to my subsample catalog
-located at /mnt/astrophysics/muchogalfit-input-cats/v2_snrcoadd.fits
-I write this scipt in somewhat of a frenzy. no need for exotic formatting with function and class wizardry.
'''

#read table
path_to_table = '/mnt/astrophysics/muchogalfit-input-cats/'
table_name = 'v2_snrcoadd.fits'
path_to_galfit = '/mnt/astrophysics/muchogalfit-output/'
vf = Table.read(path_to_table+table_name)

#set up empty arrays with which we will fashion our new columns.
mask = np.zeros(len(vf),dtype=bool)
groupGalaxy = np.zeros(len(vf),dtype=bool)
primaryGroup = np.zeros(len(vf),dtype=bool)
ncomp = np.zeros(len(vf),dtype=int)
group_name = np.zeros(len(vf),dtype=str)

columns = [groupGalaxy,primaryGroup,ncomp,group_name]
column_names = ['group_flag','primaryGroup_flag','ncomp','group_name']

#read list of primary galaxy pathnames
primary_pathnames = ascii.read(path_to_galfit+'groupDirs.txt')
primary_counter = 0   #will monitor the number of primary galaxies the loop below encounters; this number will serve as the index for the primary pathname list (e.g., if the loop find its first primary galaxy, then navigate to index 0 of list and pull the pathname)

for n in range(len(vf)):
    VFID = vf['VFID'][n]
    if vf['sgacut_flag'][n]:  #if the galaxy is part of the subsample, then proceed
        if os.path.exists(path_to_galfit+VFID+'/galsFOV.txt'):
            fovtab = ascii.read(path_to_galfit+VFID+'/galsFOV.txt')
            group_vfids = fovtab['col1']    #the length of this textfile (i.e., number of entries)

            for num in range(len(fovtab)):
                if num==0:
                    ncomp[n] = len(fovtab)
                    primaryGroup[n] = True
                    mask[n] = False
                    primary_pathname = primary_pathnames[primary_counter]
                    print(primary_pathname)
                    name = primary_pathname.split('/')[-1] #divide pathname into strings that were separated by '/', then isolate group name
                    group_name[n] = name
                    primary_counter+=1   #set up primary counter for the next primary galaxy!
                else:
                    ncomp[vf['VFID']==group_vfids[num]] = len(fovtab)
                    primaryGroup[vf['VFID']==group_vfids[num]] = False
                    mask[vf['VFID']==group_vfids[num]] = False
                    group_name[vf['VFID']==group_vfids[num]] = name  #will have same group name as primary galaxy

        else:
            mask[n] = False   #galaxy is part of subsample. will not be masked out.
            primaryGroup[n] = False   #galaxy will not be a primary galaxy, as evidenced by the lack of a .txt file in the directory
            ncomp[n] = 1 if ncomp[n]==0 else ncomp[n]  #if ncomp[n] is already populated by a nonzero integer, lEaVe It AlOnE.
            group_name[n] = '' if len(group_name[n])==0 else vf['objname'][n]  #if group galaxy, then entry will either not be empty or will be replaced with group name later on. Similar idea as above --> if already populated, LEAVE IT ALONE.

        #if ncomp>1, then galaxy is part of a group. full stop.
        groupGalaxy[n] = True if ncomp[n]>1 else False

    else:
        groupGalaxy[n] = False   #placeholder value; will be masked out
        primaryGroup[n] = False   #placeholder value; will be masked out
        ncomp[n] = 1   #placeholder value; will be masked out
        mask[n] = True  #don't include a value at this index, since the galaxy is not part of the subsample
        group_name[n] = vf['objname'][n]   #let the group name be the objname of the galaxy

for col in range(len(columns)):
    c=MaskedColumn(columns[col],name=column_names[col],mask=mask)
    vf.add_column(c)

vf.write(path_to_table+'v2_snrcoadd.fits',overwrite=True)
