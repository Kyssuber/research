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

header = ['galname','xc','xc_err','yc','yc_err','mag','mag_err','re','re_err','nsersic','nsersic_err','BA','BA_err','PA','PA_er','sky','sky_err','err_flag','chi2nu','central_flag']
    

class output_galaxy:
    def __init__(self,galname=None,vfid=None,objname=None,vfid_v1=None,outimage=None,convflag=None,band='r',ncomp=1):
        
        self.galname=galname
        self.vfid=vfid
        self.objname=objname
        self.vfid_v1=vfid_v1
        
        if vfid in dummycat['central galaxy']:
            self.ncomp = len(np.where(dummycat['central galaxy'] == vfid)[0]) + 1
        else:
            self.ncomp=ncomp
        
        #if band is an integer, than prepend with w to indicate that the 'band' is a WISE channel
        try:
            int(band)
            self.band = 'w'+str(band)
        except:
            self.band = band

        #outimage = str(self.galname)+'-'+str(self.band)+'-'+str(self.ncomp)+'Comp-galfit-out.fits'
        #self.outimage=outimage
        outimage = str(self.objname)+'-'+str(self.band)+'-out2.fits'
        self.outimage=outimage
        print(self.outimage)
        
        
    #takes GALFIT .fits output header and extracts the parameter/error values, outputs list of fit parameters
    def parse_galfit(self,printflag=False):
        
        numerical_error_flag=0
        #the following assumes that the maximum number of sersic objects in a given cutout is 3. number may differ for WISESize.
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

        fit_parameters=[]
        
        for n in range(0,self.ncomp):
            temp=[]
            if n == 0:
                temp.append(self.galname)   #if no external galaxies, only include central galaxy; if external galaxies, then first index will represent the central galaxy...then proceed to next index.
            else:
                indices = np.where(dummycat['central galaxy'] == self.vfid)[0]  #find where external(s) have self.vfid as 'host'
                index = indices[n-1]   #only want one index; if n=1, then we want the first external galaxy, meaning the 0th element in the indices list
                temp.append(dummycat['ID'][index])   #append the name (ID) of the external galaxy

            working_dir=os.getcwd()+'/'
            
            image_header = fits.getheader(self.outimage,2)
            for hkey in header_keywords[n]:
                s=str(image_header[hkey])
                if s.find('[') > -1:
                    s=s.replace('[','')
                    s=s.replace(']','')
                    t=s.split('+/-')
                    values=(float(t[0]),0.)# fit and error
                else:
                    t=s.split('+/-')
                    try:
                        values=(float(t[0]),float(t[1]))# fit and error
                        temp.append(values[0])
                        temp.append(values[1])
                    except ValueError:
                        # look for * in the string, which indicates numerical problem
                        if t[0].find('*') > -1:
                            numerical_error_flag=1
                            t[0]=t[0].replace('*','')
                            t[1]=t[1].replace('*','')
                            values=(float(t[0]),float(t[1]))# fit and error
                            temp.append(values[0])
                            temp.append(values[1])
                    except IndexError: # for CHI2NU
                        chi2nu=float(t[0])
                        continue
                if printflag:
                    print('{:6s}: {:s}'.format(hkey,s))
            temp.append(numerical_error_flag)
            temp.append(chi2nu)
            if n == 0:   #if galaxy is "central", then add a '1' flag
                temp.append(1)
            else:
                temp.append(0)
            print(temp)
            fit_parameters.append(temp)
        #print(len(fit_parameters))
        return fit_parameters


if __name__ == '__main__':
    
    dtype=[str,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float]
    full_sample_table = Table(names=header,dtype=dtype)
    
    convflag = input('conv? enter 0 (n) or 1 (y): ')
    #band = input('band? enter 1-4 for w, r for r-band: ')
    
    for i in range(0,len(cat)):
           
        #add row of zeros for ith central galaxy and each "off-centered" sersic object
        g = output_galaxy(galname=cat['prefix'][i],objname=cat['objname'][i], vfid=cat['VFID'][i], vfid_v1=cat['VFID_V1'][i], band='r') 
        num_rows = int(g.ncomp)   
        for i in range(0,num_rows):
            zero_row = np.zeros(len(dtype))
            zero_row = np.ndarray.tolist(zero_row) #convert to list, as the np array is all floats and won't allow the zeroth index element to be replaced with a string.
            zero_row[0] = '        ' #number of spaces corresponding to number of characters for VFID####. this formatting is necessary, for whatever reason. If I only inserted one space, then the resulting cell would read 'V'. Unhelpful.
            full_sample_table.add_row(zero_row)
        
        #if path exists (os.path.exists() returns True), then galaxy directory exists...
        #helpful for testing purposes
        galfit_dir = '/mnt/astrophysics/muchogalfit-output/'
        galfit_dir_one = galfit_dir+g.objname
        if os.path.exists(galfit_dir_one):
            os.chdir(galfit_dir_one)

            one_gal_table = Table(names=header,dtype=[str,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float])

            if convflag == 1:
                g.outimage = str(g.galname)+str(g.band)+'-'+str(g.ncomp)+'Comp-galfit-out-conv.fits'

            param_rows = g.parse_galfit()

            for n in range(0,len(param_rows)):
                
                '''FUNCTIONALITY EXAMPLE: if ncomp=2, then there should be 2 additional zero rows for ith galaxy, appended at the bottom of the table. 
                Beginning with n=1 (central galaxy), its corresponding zeroth row will be at index tab_length - (n) - 1. the -1 is to accommodate the length being some number N but the final index of the table being N-1.
                We replace that row of zeros with the output parameters, if the galaxy output directory exists, then proceed to the n=2 galaxy. repeat.'''
                
                current_table_length = len(full_sample_table)
                zeroth_row_index = current_table_length-n-1
                zeroth_row_index = int(zeroth_row_index)
                #repopulate row at this index with n
                full_sample_table[zeroth_row_index] = param_rows[n]
#no need for such trickery with the central galaxy table...simply add the row
                one_gal_table.add_row(param_rows[n])
                
                
                    
            #write table of central (+external, if any) galaxy parameters in directory
            if int(convflag) == 1:
                one_gal_table.write(galfit_dir+g.objname+'/output_params_'+g.band+'_psf.fits',format='fits',overwrite=True)
            if int(convflag) == 0:
                one_gal_table.write(galfit_dir+g.objname+'/output_params_'+g.band+'_nopsf.fits',format='fits',overwrite=True)
            
            
            index=len(full_sample_table)-1
            
            
    band=g.band  
        
    print(full_sample_table)
    for i in full_sample_table['galname']:
        if 'index' in i:
            print('external ID not in vf catalog, removed')
            full_sample_table.remove_row(np.where(full_sample_table['galname']==i)[0])
            
    if int(convflag) == 1:
        full_sample_table.write(galfit_dir+'output_params_'+band+'_psf.fits', format='fits', overwrite=True)
    if int(convflag) == 0:
        full_sample_table.write(galfit_dir+'output_params_'+band+'_nopsf.fits',format='fits',overwrite=True)

