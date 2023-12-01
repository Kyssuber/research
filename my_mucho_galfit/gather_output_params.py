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
    def __init__(self, vfid=None, objname=None, outimage=None, convflag=None, fovtab=None, band='W3', ncomp=1, galfit_dir=None):
        
        self.vfid=vfid
        self.objname=objname
        self.fovtab=fovtab
        self.ncomp=ncomp
        self.galfit_dir=galfit_dir
        self.vfid_dir=self.galfit_dir+self.vfid+'/'

        self.band = band

        #I guess 1+ sersic obj fitting is not set up for W1 yet, so this is a failsafe.
        if self.ncomp>1:
            print(self.ncomp)
            print(self.ncomp)
            print(self.ncomp)

        outimage_path = self.vfid_dir+str(self.objname)+'*'+str(self.band)+'-out1.fits'
        if convflag==str(1):
            outimage_path = self.vfid_dir+str(self.objname)+'*'+str(self.band)+'-out2.fits'
        try:
            self.outimage = glob.glob(outimage_path)[0]
            print(self.outimage)
        except IndexError:
            print(self.vfid+' is not a primary galaxy. skipping.')
            self.outimage='nay'  #I pass some dummy definition here since self.outimage requires a definition in order for the script to run correctly.
        
    #takes GALFIT .fits output header and extracts the parameter/error values, outputs list of fit parameters
    def parse_galfit(self,printflag=False):
        #the following assumes that the maximum number of VFID sersic objects in a given cutout is 5.
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
            header_keywords=[header_keywords1,header_keywords2,header_keywords3]            
        if self.ncomp == 4:
            header_keywords1=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','5_SKY','CHI2NU']
            header_keywords2=['2_XC','2_YC','2_MAG','2_RE','2_N','2_AR','2_PA','5_SKY','CHI2NU']
            header_keywords3=['3_XC','3_YC','3_MAG','3_RE','3_N','3_AR','3_PA','5_SKY','CHI2NU']
            header_keywords4=['4_XC','4_YC','4_MAG','4_RE','4_N','4_AR','4_PA','5_SKY','CHI2NU']
            header_keywords=[header_keywords1,header_keywords2,header_keywords3,header_keywords4]
        if self.ncomp == 5:
            header_keywords1=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','6_SKY','CHI2NU']
            header_keywords2=['2_XC','2_YC','2_MAG','2_RE','2_N','2_AR','2_PA','6_SKY','CHI2NU']
            header_keywords3=['3_XC','3_YC','3_MAG','3_RE','3_N','3_AR','3_PA','6_SKY','CHI2NU']
            header_keywords4=['4_XC','4_YC','4_MAG','4_RE','4_N','4_AR','4_PA','6_SKY','CHI2NU']
            header_keywords5=['5_XC','5_YC','5_MAG','5_RE','5_N','5_AR','5_PA','6_SKY','CHI2NU']
            header_keywords=[header_keywords1,header_keywords2,header_keywords3,header_keywords4,header_keywords5]
        if self.ncomp == 9:
            header_keywords1=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','10_SKY','CHI2NU']
            header_keywords2=['2_XC','2_YC','2_MAG','2_RE','2_N','2_AR','2_PA','10_SKY','CHI2NU']
            header_keywords3=['3_XC','3_YC','3_MAG','3_RE','3_N','3_AR','3_PA','10_SKY','CHI2NU']
            header_keywords4=['4_XC','4_YC','4_MAG','4_RE','4_N','4_AR','4_PA','10_SKY','CHI2NU']
            header_keywords5=['5_XC','5_YC','5_MAG','5_RE','5_N','5_AR','5_PA','10_SKY','CHI2NU']
            header_keywords6=['6_XC','6_YC','6_MAG','6_RE','6_N','6_AR','6_PA','10_SKY','CHI2NU']
            header_keywords7=['7_XC','7_YC','7_MAG','7_RE','7_N','7_AR','7_PA','10_SKY','CHI2NU']
            header_keywords8=['8_XC','8_YC','8_MAG','8_RE','8_N','8_AR','8_PA','10_SKY','CHI2NU']
            header_keywords9=['9_XC','9_YC','9_MAG','9_RE','9_N','9_AR','9_PA','10_SKY','CHI2NU']
            header_keywords=[header_keywords1,header_keywords2,header_keywords3,header_keywords4,header_keywords5,header_keywords6,
                header_keywords7,header_keywords8,header_keywords9]

        fit_parameters=[]
        
        #if galaxy is the 'primary' of a group, isolate all VFIDs in the group
        if self.fovtab is not None:
            group_vfids = self.fovtab['col1']
        
        #parse parameters for each galaxy in the cutout/model/etc.
        for n in range(self.ncomp):
            
            numerical_error_flag=0   #default assumption is that galfit ran successfully on the galaxy. hooray.
            
            temp=[]   #create temporary list that will hold the parameters for galaxy n
        
            if self.ncomp>1:
                temp.append(group_vfids[n])   #if ncomp>1, use the list of VFIDs to create first list item
            else:
                temp.append(self.vfid)  #if ncomp=1, simply use the central galaxy's VFID
            
            image_header = fits.getheader(self.outimage,2)   #grab header information from the model image
            
            #if there are more than 1 sersic objects in the cutout, check that GALFIT was set up such that it could successfully run on these group galaxies (if the outskirt galax(y/ies) is/are disregarded, then the sky header will only be 2_SKY and not [self.ncomp+1]_SKY).
            if (self.ncomp>1):
                sky_num = str(self.ncomp+1)
                sky_head = sky_num+'_SKY'
                print(sky_head)
                if (sky_head not in image_header):
                    print('keywords not in image_header')
                    header_keywords = [header_keywords1]
                    break
            
            for hkey in header_keywords[n]:

                s=str(image_header[hkey])
            
                if s.find('[')+1 > 0:   #if a parameter is held fixed (meaning the parameter is enclosed with brackets and does not contain an uncertainty)
                    s=s.replace('[','')
                    s=s.replace(']','')
                    s=f'{s} +/- {0.}'   #fit and error
                
                t=s.split('+/-')
                
                try:
                    values=(float(t[0]),float(t[1]))   #fit and error
                    temp.append(values[0])
                    temp.append(values[1])
                except ValueError:   #for numerical errors
                    # look for * in the string, which indicates numerical problem
                    if t[0].find('*')+1 > 0:
                        numerical_error_flag=1   #numerical error is now a 1. not hooray.
                        t[0]=t[0].replace('*','')
                        t[1]=t[1].replace('*','')
                        values=(float(t[0]),float(t[1]))   #fit and error
                        temp.append(values[0])
                        temp.append(values[1])
                except IndexError:   #for CHI2NU
                        chi2nu=float(t[0])
                        continue
        
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
    
    #define the 'parent' galfit path
    galfit_dir = '/mnt/astrophysics/muchogalfit-output/'

    #for every galaxy in the VF subsample catalog
    for i in range(len(cat)):        

        #define output_galaxy class; sets up the vf subsample catalog, the objname, the VFID, and the conflag/band.
        g = output_galaxy(objname=cat['objname'][i], vfid=cat['VFID'][i], fovtab=None, convflag=convflag, band=band, galfit_dir=galfit_dir) 
        
        galfit_dir_one = g.vfid_dir

        #check whether the galaxy is a member of a group; if so, then the directory will contain galsFOV.txt file
        if os.path.exists(galfit_dir_one+'/'+'galsFOV.txt'):
            g.fovtab=ascii.read(galfit_dir_one+'/'+'galsFOV.txt')   #define as class variable; will need for g.parse_galfit()
            g.ncomp=len(g.fovtab)   #the length of this textfile (i.e., number of entries)
        #print(g.ncomp)
        #if g.band == 'W1':
        #    g.ncomp=1          #1+ sersic obj fitting not yet set up for W1. this line is a failsafe.

        #next, check whether this galaxy is 'primary' or part of a group (if the latter, then there will be no galfit output in the directory). if not, then the the main loop will proceed onward to the next VFID
        #recall that g.outimage is the FULL pathname to the galfit output file(s) for one galaxy
        if os.path.isfile(g.outimage):
            os.chdir(galfit_dir_one)
            
            #create a list of param lists (if ncomp=1, then one param list [[]]; if ncomp=2, then two param lists [[],[]])
            param_rows = g.parse_galfit()
            for row in param_rows:
                row_vfid = row[0]  #VFID is the first entry in the row list
                #if len(row)==len(dtypes):   #if no parameters are set fixed, then add row to table (in a few group cases, Rose set parameters to the fixed no-psf values in order to avoid galfit crashing)
                full_sample_table[full_sample_table['VFID']==row_vfid] = row   #change zeros row to parameter row
                            
            index=len(full_sample_table)-1

    print(full_sample_table)
            
    if int(convflag) == 1:
        full_sample_table.write(galfit_dir+'output_params_'+g.band+'_psf.fits', format='fits', overwrite=True)
    if int(convflag) == 0:
        full_sample_table.write(galfit_dir+'output_params_'+g.band+'_nopsf.fits', format='fits', overwrite=True)
