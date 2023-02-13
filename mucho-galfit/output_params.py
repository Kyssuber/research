#!/usr/bin/env python
# coding: utf-8

'''
GOAL:
- gather all GALFIT output parameters, as well as the galnames and 'central' designations, and compile into an Astropy table.
- only includes galaxies with VFIDs. entries in 'dummy catalog' which do not have a VFID are not added to the table.
- be sure directories contain all GALFIT output .fits files before running. 
- output will be multiple .fits tables, one for each central galaxy and one that combines parameters of all galaxies.
'''

import sys
import numpy as np
import os
homedir = os.getenv("HOME")
import glob
from astropy.table import Table
from astropy.wcs import WCS
from astropy.io import ascii
from astropy.io import fits
    

class output_galaxy:
    def __init__(self,galname=None,vfid=None,objname=None,vfid_v1=None,outimage=None,convflag=None,band=None,ncomp=1):
        
        self.galname=galname
        self.vfid=vfid
        self.objname=objname
        self.vfid_v1=vfid_v1

        '''
        if vfid in dummycat['central galaxy']:
            self.ncomp = len(np.where(dummycat['central galaxy'] == vfid)[0]) + 1
        else:
            self.ncomp=ncomp
        '''
        self.band = band
        self.ncomp = ncomp
        
        pathname = '/mnt/astrophysics/muchogalfit-output/'+str(self.vfid)+'/'
        outimage = glob.glob(pathname+self.objname+'-'+str(self.band)+'-out*')[0]
        if convflag==1:
            outimage = glob.glob(pathname+self.objname+'-'+str(self.band)+'-out*')[1]
        
        self.outimage=outimage
        print(self.outimage)
        
    #takes GALFIT .fits output header and extracts the parameter/error values, outputs list of fit parameters
    def parse_galfit_1comp(self,printflag=False):
        
        numerical_error_flag=0
        #the following assumes that the maximum number of sersic objects in a given cutout is 3. number may differ in WISESize.
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
            #print(temp)
            fit_parameters.append(temp)
        #print(len(fit_parameters))
        return fit_parameters


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
    with open('/mnt/astrophysics/kconger_wisesize/github/research/mucho-galfit/'+param_file) as f:
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
    gal_output_dir = param_dict['gal_output_path']
    band = param_dict['band']
    convflag = int(param_dict['convflag'])
    test = bool(param_dict['test'])
    

    header = ['galname_v2','xc','xc_err','yc','yc_err','mag','mag_err','re','re_err','nsersic','nsersic_err','BA','BA_err','PA','PA_er','sky','sky_err','err_flag','chi2nu','central_flag']
    
    dtype=[np.object,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float]
    
    full_sample_table = Table(names=header,dtype=dtype)
    
    for i in range(0,len(cat)):
        
        g = output_galaxy(galname=cat['prefix'][i], vfid=cat['VFID'][i], vfid_v1 = cat['VFID_V1'][i],convflag=convflag, band=band)
        
        if test == True:
            
            num_rows = int(g.ncomp)
            for i in range(num_rows):
                zero_row = np.zeros(len(dtype))
                zero_row = np.ndarray.tolist(zero_row) 
                zero_row[0] = '        '
                full_sample_table.add_row(zero_row)
            
            os.chdir(gal_output_dir+self.vfid)
            
            one_gal_table = Table(names=header,dtype=dtype)
            
            #if convflag == 1:
            #    g.outimage = str(g.galname)+str(g.band)+'-'+str(g.ncomp)+'Comp-galfit-out-conv.fits'
                
            param_rows = g.parse_galfit_1comp()
            
            for n in range(len(param_rows)):
                '''FUNCTIONALITY EXAMPLE: if ncomp=2, then there should be 2 additional zero rows for ith galaxy, appended at the bottom of the table. 
                Beginning with n=1 (central galaxy), its corresponding zeroth row will be at index tab_length - (n) - 1. the -1 is to accommodate the length being some number N but the final index of the table being N-1.
                We replace that row of zeros with the output parameters, if the galaxy output directory exists, then proceed to the n=2 galaxy. repeat.'''
                
                current_table_length = len(full_sample_table)
                zero_row_index = int(current_table_length - n - 1)
                #repopulate this row with nth param_row
                full_sample_table[zero_row_index] = param_rows[n]
                one_gal_table.add_row(param_rows[n])  
        
        
        if test == False:

            #if convflag == 1:
            #    g.outimage = str(g.galname)+str(g.band)+'-'+str(g.ncomp)+'Comp-galfit-out-conv.fits'

            param_rows = g.parse_galfit_1comp()
            for n in param_rows:
                if 'index' in n[0]:
                    print('external ID not in vf catalog')
                else:
                    t.add_row(n)
            band=g.band  
        
    print(full_sample_table)            
    if int(convflag) == 1:
        full_sample_table.write(homedir+'/output_params_'+band+'_psf.txt', format='ascii', overwrite=True)
    if int(convflag) == 0:
        full_sample_table.write(homedir+'/output_params_'+band+'_nopsf.txt',format='ascii',overwrite=True)
      
    

