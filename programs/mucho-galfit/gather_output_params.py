#!/usr/bin/env python
# coding: utf-8

'''
GOAL:
- gather all GALFIT output parameters, as well as the galnames and 'central' designations, and compile into an Astropy table.
- only includes galaxies with VFIDs. entries in 'dummy catalog' which do not have a VFID are not added to the table.
- be sure directories contain all GALFIT output .fits files before running. 
- output will be multiple .fits tables, one for each central galaxy and one that combines parameters of all galaxies.
'''

import numpy as np
import os
homedir = os.getenv("HOME")
import glob
from astropy.table import Table
from astropy.wcs import WCS
from astropy.io import ascii
from astropy.io import fits

#have to call sample
homedir = os.getenv("HOME")
mnt_cat_path = '/mnt/astrophysics/muchogalfit-input-cats/'
cat = Table.read(mnt_cat_path+'sgacut_coadd.fits')
dummycat = Table.read(mnt_cat_path+'dummycat.fits')
vf = Table.read(mnt_cat_path+'vf_v2_main.fits')


header = ['VFID','xc','xc_err','yc','yc_err','mag','mag_err','re','re_err','nsersic','nsersic_err','BA','BA_err','PA','PA_er','sky','sky_err','err_flag','chi2nu','central_flag']
    

class output_galaxy:
    def __init__(self, vfid=None, objname=None, outimage=None, convflag=None, fovtab=None, band='W3', ncomp=1):
        
        self.vfid=vfid
        self.objname=objname
        self.fovtab=fovtab
        self.ncomp=ncomp
        
        #if band is an integer, than prepend with w to indicate that the 'band' is a WISE channel
        self.band = band

        outimage = glob.glob(str(self.objname)+'*'+str(self.band)+'-out1.fits')[0]
        if convflag==str(1):
            outimage = glob.glob(str(self.objname)+'*'+str(self.band)+'out2.fits')[0]
        self.outimage=outimage
        print(self.outimage)
        
    #takes GALFIT .fits output header and extracts the parameter/error values, outputs list of fit parameters
    def parse_galfit(self,printflag=False):
        
        #the following assumes that the maximum number of VFID sersic objects in a given cutout is 4.
        if self.ncomp == 1:
            header_keywords1=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_SKY','CHI2NU']
            header_keywords=[header_keywords1]
        if self.ncomp == 2:
            header_keywords1=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','3_SKY','CHI2NU']
            header_keywords2=['2_XC','2_YC','2_MAG','2_RE','2_N','2_AR','2_PA','3_SKY','CHI2NU']
            header_keywords=[header_keywords1,header_keywords2]
        if self.ncomp == 3:
            header_keywords1=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','4_SKY','CHI2NU']
            header_keywords2=['2_XC','2_YC','2_MAG','2_RE','2_N','2_AR','2_PA','4_SKY','CHI2NU']
            header_keywords3=['3_XC','3_YC','3_MAG','3_RE','3_N','3_AR','3_AR','4_SKY','CHI2NU']
            
        if self.ncomp == 3:
            header_keywords1=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','4_SKY','CHI2NU']
            header_keywords2=['2_XC','2_YC','2_MAG','2_RE','2_N','2_AR','2_PA','4_SKY','CHI2NU']
            header_keywords3=['3_XC','3_YC','3_MAG','3_RE','3_N','3_AR','3_PA','4_SKY','CHI2NU']
            header_keywords4=['4_XC','4_YC','4_MAG','4_RE','4_N','4_AR','4_PA','5_SKY','CHI2NU']
            header_keywords=[header_keywords1,header_keywords2,header_keywords3,header_keywords4]

        fit_parameters=[]
        
        #if galaxy is the 'primary' of a group, isolate all VFIDs in the group
        if self.fovtab is None:
            group_vfids = self.fovtab['col1']
        
        #parse parameters for each galaxy in the cutout/model/etc.
        for n in range(self.ncomp):
            
            numerical_error_flag=0   #default assumption is that galfit ran successfully on the galaxy. hooray.
            
            temp=[]   #create temporary list that will hold the parameters for galaxy n
        
            if self.ncomp>1:
                temp.append(group_vfids[n])   #if ncomp>1, use the list of VFIDs to create first list item
            else:
                temp.append(self.VFID)  #if ncomp=1, simply use the central galaxy's VFID
            
            image_header = fits.getheader(self.outimage,2)   #grab header information from the model image
            for hkey in header_keywords[n]:
                s=str(image_header[hkey])
                if s.find('[') > -1:
                    s=s.replace('[','')
                    s=s.replace(']','')
                    t=s.split('+/-')
                    values=(float(t[0]),0.)   #fit and error
                else:
                    t=s.split('+/-')
                    try:
                        values=(float(t[0]),float(t[1]))   #fit and error
                        temp.append(values[0])
                        temp.append(values[1])
                    except ValueError:
                        # look for * in the string, which indicates numerical problem
                        if t[0].find('*') > -1:
                            numerical_error_flag=1   #numerical error is now a 1. not hooray.
                            t[0]=t[0].replace('*','')
                            t[1]=t[1].replace('*','')
                            values=(float(t[0]),float(t[1]))   #fit and error
                            temp.append(values[0])
                            temp.append(values[1])
                    except IndexError:   #for CHI2NU
                        chi2nu=float(t[0])
                        continue
                if printflag:
                    print('{:6s}: {:s}'.format(hkey,s))
            temp.append(numerical_error_flag)
            temp.append(chi2nu)
            if n == 0:   #if galaxy is "primary" or central, then add a '1' flag
                temp.append(1)
            else:
                temp.append(0)

            fit_parameters.append(temp)
        return fit_parameters

if __name__ == '__main__':
    
    #create empty table
    dtypes=[str,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float]
    full_sample_table = Table(names=header,dtype=dtypes)
    
    #for every galaxy in the vf subsample, add an empty row to the table
    for i in range(len(cat)):
        full_sample_table.add_row()
    
    #populate the VFID column with all subsample VFIDs
    full_sample_table['VFID'] = cat['VFID']
    
    #try automating? naw.
    convflag = input('conv? enter 0 (n) or 1 (y): ')
    band = input('band? enter W1-4, r for r-band: ')

    #for every galaxy in the VF subsample catalog
    for i in range(len(cat)):        
        
        #define output_galaxy class; sets up the vf subsample catalog, the objname, the VFID, and the conflag/band.
        g = output_galaxy(objname=cat['objname'][i], vfid=cat['VFID'][i], fovtab=None, convflag=convflag, band=band) 

        #define the 'parent' galfit path as well as the directory for the current galaxy 
        galfit_dir = '/mnt/astrophysics/muchogalfit-output/'
        galfit_dir_one = galfit_dir+g.vfid
        
        #check whether the galaxy is a member of a group; if so, then the directory will contain galsFOV.txt file
        if os.path.exists(galfit_dir_one+'/'+'galsFOV.txt'):
            g.fovtab=ascii.read(galfit_dir_one+'/'+'galsFOV.txt')   #define as class variable; will need for g.parse_galfit()
            g.ncomp=len(g.fovtab)   #the length of this textfile (i.e., number of entries)
        
        #next, check whether this galaxy is 'primary' or part of a group (if the latter, then there will be no galfit output in the directory). if not, then the the main loop will proceed onward to the next VFID
        if os.path.isfile(galfit_dir_one+'/'+g.outimage):
            os.chdir(galfit_dir_one)
            
            #create a list of param lists (if ncomp=1, then one param list [[]]; if ncomp=2, then two param lists [[],[]])
            param_rows = g.parse_galfit()

            for row in range(len(param_rows)):
                row_vfid = row[0]  #VFID is the first entry in the row list
                full_sample_table[full_sample_table['VFID']==row_vfid] = row   #change zeros row to parameter row
                            
            index=len(full_sample_table)-1

    print(full_sample_table)
            
    if int(convflag) == 1:
        full_sample_table.write(galfit_dir+'output_params_'+g.band+'_psf.fits', format='fits', overwrite=True)
    if int(convflag) == 0:
        full_sample_table.write(galfit_dir+'output_params_'+g.band+'_nopsf.fits', format='fits', overwrite=True)