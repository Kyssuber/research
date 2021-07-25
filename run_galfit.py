#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#galaxy class is adapted from Rose's code.
#only added vfid arg to galaxy class and edited get_wise_image to grab cutouts with corrected oversubtraction halos

###########################################################################################################################

import os
#import pyds9
import numpy as np
import argparse
from astropy.io import fits
from astropy.visualization import simple_norm
from astropy import units as u
from astropy.wcs import WCS
from astropy.table import Table
from astropy.io import ascii
from matplotlib import pyplot as plt
from scipy.stats import scoreatpercentile
import wget
import tarfile
import glob
import gzip

import csv
import shutil

import astropy.wcs as wcs

homedir = os.getenv("HOME")

os.sys.path.append(homedir+'/github/virgowise/')
import rungalfit as rg #This code has all the defined functions that I can use
os.sys.path.append(homedir+'/github/HalphaImaging/python3/')
import plot_cutouts_ha as cutouts #This code has all the defined functions that I can use




       
class galaxy():
    def __init__(self,ra,dec,size,name='galname',vfid='VFID',band='3'):
        '''
        galaxy for wise analysis

        params:
        -------
        * ra = ra of gal center in deg
        * dec = dec of gal center in deg
        * size = D25/2, or some other estimate of radius in arcsec

        optional params:
        ---------------
        * name = name of galaxy to use when saving images, default is galname
        * band = WISE band, default is '3'
          - 1 = W1
          - 2 = W2
          - 3 = W3
          - 4 = W4
        
        '''
        self.ra = ra#*u.deg
        self.dec = dec#*u.deg
        self.radius = size#*u.arcsec
        self.band = band
        self.galname = name
        self.vfid = vfid
        self.image_rootname = self.galname+'-unwise-w'+str(self.band)
        self.image = self.image_rootname+'-img-m.fits'

        self.mask_image = self.galname+'-unwise-mask.fits'
        self.sigma_image = self.image_rootname+'-std-m.fits'
        self.invvar_image = self.image_rootname+'-invvar-m.fits'

        # remove the log file if it exists
        self.logfilename = self.galname+'-unwise-'+'w'+str(self.band)+'-log.txt'
        if os.path.exists(self.logfilename):
            os.remove(self.logfilename)

        # write the header line into the log file
        output=open(self.logfilename,'w')
        output.write('# xc xc_err yc yc_err mag mag_err re re_err nsersic nsrsic_err BA BA_err PA PA_err sky sky_err error chi2nu galname \n')
        # close log file
        output.close()    
        
        
#obsolete, but I will keep for documentation purposes
   
    def get_wise_image_59(self,makeplots=False):
        base_dir = homedir+'/github/unwise_fixed/'
        
        self.image = base_dir + str(self.vfid) + '/unwise-' + str(self.vfid) + '-w3-img-m.fits'
        self.sigma_image = base_dir + str(self.vfid) + '/unwise-' + str(self.vfid) + '-w3-std-m.fits'
        temp = fits.getdata(self.image)
        print(temp.shape)
        self.ximagesize,self.yimagesize=temp.shape



# Definition modified such that if galaxy VFID is in list containing the 59 subsample, then grab corrected fits image from 
# directories created by Dustin.
# Else...use Rose's routine. :-)

    def get_wise_image(self,makeplots=False):
        '''
        GOAL: Get the unWISE image from the unWISE catalog

        INPUT: nsaid used to grab unwise image information

        OUTPUT: Name of file to retrieve from

        '''
        vf = Table.read(homedir+'/sample_main.fits')
        #vfmain = Table.read(homedir+'/github/research/sample_main.fits')
        base_dir = homedir+'/github/unwise_100px_dstn'
        vfmain_list = []
        for i in range(0,len(vf)):
            vfmain_list.append(vf['VFID'][i])
        
        if str(self.vfid) in vfmain_list: 
            self.image = base_dir + '/' + str(self.vfid) + '/unwise-' + str(self.vfid) + '-w3-img-m.fits'
            self.sigma_image = base_dir +'/' + str(self.vfid) + '/unwise-' + str(self.vfid) + '-w3-std-m.fits'
            temp = fits.getdata(self.image)
            print(temp.shape)
            self.ximagesize,self.yimagesize=temp.shape
        else:
            baseurl = 'http://unwise.me/cutout_fits?version=allwise'
            imsize = self.radius*2

            imagenames,weightnames,multiframe = cutouts.get_unwise_image(self.ra,self.dec,galid=self.galname,pixscale=1,imsize=self.radius*2,bands=self.band,makeplots=makeplots,subfolder=None)
        
            self.image = imagenames[0]
        
            self.sigma_image = weightnames[0]

            # read in image and get image size

            temp = fits.getdata(self.image)
            print(temp.shape)
            self.ximagesize,self.yimagesize = temp.shape
        
            #print(self.sigma_image)
            #print(imagenames)

    
    
    #this function is my contribution...
    
    #grab cutout from allwise catalog
    
    def get_allwise_image(self):
        ra_ex = self.ra
        dec_ex = self.dec
        size_ex = self.radius
        name_ex = self.galname  

        #the following yields a .csv file containing allwise four (w1-w4) data rows for the given coords.
        irsa_call = '\"https://irsa.ipac.caltech.edu/ibe/search/wise/allwise/p3am_cdd?POS=' + str(ra_ex) + ',' + str(dec_ex) + '&ct=csv&mcen\"'
        filename = homedir+'/github' + '/allwise.csv'
        if os.path.exists(filename):
            os.remove(filename)

        #copies .csv to pathname (i.e., filename)
        os.system("curl -o %s %s" % (filename, irsa_call))

        #this .csv information is then extractable to acquire the desired .fits cutouts, like so:
        outdir = homedir+'/github/allwise_cutouts/'

        with open(filename) as csvfile:
            reader = csv.DictReader(csvfile)
            for line in reader:
                band = line['band']
                if int(band) != 3:
                    pass
                else:
                    coadd_id = line['coadd_id']
                    coaddgrp = coadd_id[:2]
                    coadd_ra = coadd_id[:4]
                    #int for intensity, unc for uncertainty; only pulls w3 files
                    filename_int = coadd_id + '-w3'+'-int-3.fits.gz'
                    filename_unc = coadd_id + '-w3'+'-unc-3.fits.gz'
                    #url dongxi
                    irsa_url = "%s/%s/%s/" % (coaddgrp, coadd_ra, coadd_id)
                    irsa_url_int = irsa_url + filename_int
                    irsa_url_unc = irsa_url + filename_unc
                    irsa_url_int = '\"https://irsa.ipac.caltech.edu/ibe/data/wise/allwise/p3am_cdd/' + irsa_url_int +'?center=' + str(ra_ex) + ',' + str(dec_ex) +'&size=' + str(size_ex*2) + 'arcsec' + '\"'
                    irsa_url_unc = '\"https://irsa.ipac.caltech.edu/ibe/data/wise/allwise/p3am_cdd/' + irsa_url_unc +'?center=' + str(ra_ex) + ',' + str(dec_ex) +'&size=' + str(size_ex*2) + 'arcsec' + '\"'
                    filename_int = outdir + filename_int
                    filename_unc = outdir + filename_unc
                    #grab int .gz
                    os.system("curl -o %s %s" % (filename_int, irsa_url_int))
                    #rename file
                    name_desired_int = outdir + str(name_ex) + '-int.fits.gz'
                    os.rename(filename_int, name_desired_int)
                    #grab unc .gz
                    os.system("curl -o %s %s" % (filename_unc, irsa_url_unc))
                    #rename file
                    name_desired_unc = outdir + str(name_ex) + '-unc.fits.gz'
                    os.rename(filename_unc, name_desired_unc)                       
        
        
        #decompress resulting gzip files (from Stackoverflow)
        
        #int, name_desired_int, name_desired_int_nogz
        name_desired_int_nogz = outdir + str(name_ex) + '-int.fits'
        with gzip.open(name_desired_int,'rb') as f_in:
            with open(name_desired_int_nogz,'wb') as f_out:
                shutil.copyfileobj(f_in,f_out)
        
        #unc, name_desired_unc, name_desired_unc_nogz
        name_desired_unc_nogz = outdir + str(name_ex) + '-unc.fits'
        with gzip.open(name_desired_unc,'rb') as f_in:
            with open(name_desired_unc_nogz,'wb') as f_out:
                shutil.copyfileobj(f_in,f_out)                  
                          
                          
        self.image = name_desired_int_nogz
        
        self.sigma_image = name_desired_unc_nogz

        # read in image and get image size

        temp = fits.getdata(self.image)
        print(temp.shape)
        self.ximagesize,self.yimagesize = temp.shape                  
                          
                          
                          
                          
                          
        
    def set_image_names(self):
        '''
        GOAL:
        * Set the psf image name and some additional parameters for galfit 

        PARAMS:
        * none

        OUTPUT: PSF image and optional diffusion kernel, other parameters  

        '''
        #just using center til, doesn't matter usually
        os.system('cp '+homedir+'/github/virgowise/wise_psfs/wise-w3-psf-wpro-09x09-05x05.fits .') 
        self.psf_image = 'wise-w3-psf-wpro-09x09-05x05.fits' 
        self.psf_oversampling = 8
        #mask_image = 'testimage_mask.fits' no mask image 
        self.xminfit=0
        self.yminfit=0
        self.xmaxfit=self.ximagesize
        self.ymaxfit=self.yimagesize
        self.convolution_size=self.ximagesize

        ### NEED TO UPDATE MAGZP AND PSCALE FOR UNWISE
        self.magzp=22.5 #For WISE
        self.pscale=0.344 #For WISE
        
        #convflag=1 # apply psf convolution
        #constraintflag=1 # add a constraint file?
        #self.fitallflag=0
        self.ncomp=2

    def getpix(self):
        '''
        GOAL: Get pixel values that correspond to ra and dec

        INPUT: nsaid
        
        OUTPUT: xc, yc in pixels

        '''

        w=wcs.WCS(self.image)
        self.xc, self.yc = w.wcs_world2pix(self.ra, self.dec,1)


    def set_sersic_params_random(self):
        '''
        GOAL: Set random parameters for galfit

        INPUT: nsaid

        OUTPUT: 5 random parameters for nsersic, magnitude, effective radius, axis ratio, and position angle

        '''
        self.nsersic = 5.5*np.random.random()+.5
        self.mag =14*np.random.random()+2
        self.re = 60*np.random.random()
        self.BA = np.random.random()
        self.PA =181*np.random.random()-89.0

    def set_sersic_manual(self,n=2,m=7,re=5,BA=1,PA=0):
        '''
        GOAL: Set random parameters for galfit

        INPUT: nsaid

        OUTPUT: 5 random parameters for nsersic, magnitude, effective radius, axis ratio, and position angle

        '''
        self.nsersic = n
        self.mag =m
        self.re = re
        self.BA = BA
        self.PA = PA
                
    def initialize_galfit(self,convflag=True):
        '''
        GOAL: Preparing file to be run in galfit. Initialize galfit image parameters 

        INPUT: nsaid 

        OUTPUT: A definition of everything from galname to convflag, necessary for running galfit

        '''
        print('self.psfimage = ',self.psf_image)
        
        self.gal1 = rg.galfit(galname=self.image_rootname,image=self.image, mask_image = self.mask_image, sigma_image=self.sigma_image,psf_image=self.psf_image,psf_oversampling=self.psf_oversampling,xminfit=self.xminfit,yminfit=self.yminfit,xmaxfit=self.xmaxfit,ymaxfit=self.ymaxfit,convolution_size=self.convolution_size,magzp=self.magzp,pscale=self.pscale,ncomp=self.ncomp,convflag=convflag)
        
        
        
    def run_galfit_wise(self,fitBA=1,fitPA=1):
        '''
        GOAL: 
        * run galfit on one image

        optional params:
        ----------------
        * fitBA = set to 1 to let galfit fit the axis ratio BA of the galaxy
        * fitPA = set to 1 to let galfit fit the PA position angle of the galaxy

        OUTPUT: several output files

        '''
        #os.system('cp '+self.psf_image+' .')
        self.gal1.set_sersic_params(xobj=self.xc,yobj=self.yc,mag=self.mag,rad=self.re,nsersic=self.nsersic,BA=self.BA,PA=self.PA,fitmag=1,fitcenter=1,fitrad=1,fitBA=fitBA,fitPA=fitPA,fitn=1,first_time=0)
     
        
        self.gal1.set_sky(0)
        self.gal1.run_galfit()
    def get_galfit_results(self,printflag = False):
        '''
        GOAL: 
        -----
        * Grab results from galfit (xc, yc, mag, re, nsersic, BA, PA, sky, error, chi2nu) 
          and parse them into self.filename

        PARAMS:
        -------
        * none

        OPTIONAL PARAMS:
        ----------------
        * printflag = print fit results, default is False

        OUTPUT:
        -------
        * stores fit results in variables 

        '''

        self.filename = self.gal1.output_image
        #self.galname+'-unwise-'+'w'+str(self.band)+'-1Comp-galfit-out.fits'
        t = rg.parse_galfit_1comp(self.filename)
        if printflag:
            self.gal1.print_galfit_results(self.filename)
        
        if int(self.ncomp) == 1:
            header_keywords=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_SKY','ERROR','CHI2NU','galname']
            self.xc, self.xc_err = t[0]
            self.yc, self.yc_err = t[1]
            self.mag, self.mag_err = t[2]
            self.re, self.re_err = t[3]
            self.nsersic, self.nsersic_err = t[4]
            self.BA, self.BA_err = t[5]
            self.PA, self.PA_err = t[6]
            self.sky, self.sky_err = t[7]
            self.error = t[8]
            self.chi2nu = t[9]        
        
        if int(self.ncomp) == 2:
            header_keywords= ['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_XC','2_YC','2_MAG','2_RE','2_N','2_AR','2_PA','3_SKY','CHI2NU','galname']
            self.xc, self.xc_err = t[0]
            self.yc, self.yc_err = t[1]
            self.mag, self.mag_err = t[2]
            self.re, self.re_err = t[3]
            self.nsersic, self.nsersic_err = t[4]
            self.BA, self.BA_err = t[5]
            self.PA, self.PA_err = t[6]
            self.xc2, self.xc2_err = t[7]
            self.yc2, self.yc2_err = t[8]
            self.mag2, self.mag2_err = t[9]
            self.re2, self.re2_err = t[10]
            self.nsersic2, self.nsersic2_err = t[11]
            self.PA2, self.PA2_err = t[12]
            self.BA2, self.BA2_err = t[13]            
            self.sky, self.sky_err = t[14]
            self.error = t[15]
            self.chi2nu = t[16]

    
        if int(self.ncomp) == 3:
            header_keywords= ['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_XC','2_YC','2_MAG','2_RE','2_N','2_AR','2_PA','3_XC','3_YC','3_MAG','3_RE','3_N','3_AR','3_PA','4_SKY','CHI2NU','galname']
            self.xc, self.xc_err = t[0]
            self.yc, self.yc_err = t[1]
            self.mag, self.mag_err = t[2]
            self.re, self.re_err = t[3]
            self.nsersic, self.nsersic_err = t[4]
            self.BA, self.BA_err = t[5]
            self.PA, self.PA_err = t[6]
            self.xc2, self.xc2_err = t[7]
            self.yc2, self.yc2_err = t[8]
            self.mag2, self.mag2_err = t[9]
            self.re2, self.re2_err = t[10]
            self.nsersic2, self.nsersic2_err = t[11]
            self.PA2, self.PA2_err = t[12]
            self.BA2, self.BA2_err = t[13]            
            self.xc3, self.xc3_err = t[14]
            self.yc3, self.yc3_err = t[15]
            self.mag3, self.mag3_err = t[16]
            self.re3, self.re3_err = t[17]
            self.nsersic3, self.nsersic3_err = t[18]
            self.PA3, self.PA3_err = t[19]
            self.BA3, self.BA3_err = t[20]            
            self.sky, self.sky_err = t[21]
            self.error = t[22]
            self.chi2nu = t[23]   
    
    
        if int(self.ncomp) == 4:
            header_keywords= ['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_XC','2_YC','2_MAG','2_RE','2_N','2_AR','2_PA','3_XC','3_YC','3_MAG','3_RE','3_N','3_AR','3_PA','4_XC','4_YC','4_MAG','4_RE','4_N','4_AR','4_PA','5_SKY','CHI2NU','galname']
            self.xc, self.xc_err = t[0]
            self.yc, self.yc_err = t[1]
            self.mag, self.mag_err = t[2]
            self.re, self.re_err = t[3]
            self.nsersic, self.nsersic_err = t[4]
            self.BA, self.BA_err = t[5]
            self.PA, self.PA_err = t[6]
            self.xc2, self.xc2_err = t[7]
            self.yc2, self.yc2_err = t[8]
            self.mag2, self.mag2_err = t[9]
            self.re2, self.re2_err = t[10]
            self.nsersic2, self.nsersic2_err = t[11]
            self.PA2, self.PA2_err = t[12]
            self.BA2, self.BA2_err = t[13]            
            self.xc3, self.xc3_err = t[14]
            self.yc3, self.yc3_err = t[15]
            self.mag3, self.mag3_err = t[16]
            self.re3, self.re3_err = t[17]
            self.nsersic3, self.nsersic3_err = t[18]
            self.PA3, self.PA3_err = t[19]
            self.BA3, self.BA3_err = t[20] 
            self.xc4, self.xc4_err = t[21]
            self.yc4, self.yc4_err = t[22]
            self.mag4, self.mag4_err = t[23]
            self.re4, self.re4_err = t[24]
            self.nsersic4, self.nsersic4_err = t[25]
            self.PA4, self.PA4_err = t[26]
            self.BA4, self.BA4_err = t[27] 
            self.sky, self.sky_err = t[28]
            self.error = t[29]
            self.chi2nu = t[30]    
    
    
    
    
    def write_results(self,printflag=False):
        '''
        GOAL: 
        * Put results from galfit into a logfile by appending values

        PARAMS:
        * none

        OUTPUT: 
        * logfile of parameter types and their associated outputs from galfit 

        '''
        self.get_galfit_results(printflag=printflag)

        output=open(self.logfilename,'a')
        # create string with best-fit parameters
        
        
        if int(self.ncomp) == 1:    
            s = '%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f \n'%(self.xc,self.xc_err,self.yc,self.yc_err,self.mag,self.mag_err, self.re, self.re_err, self.nsersic, self.nsersic_err, self.BA, self.BA_err, self.PA, self.PA_err, self.sky, self.sky_err, self.error,self.chi2nu) + str( self.galname)
        
        #for below, may have to read dummycat.fits for second, third, and/or fourth sersic obj galnames. another option is to simply append a central and off-centered galname column...whichever strikes you as more doable.
        
        
        if int(self.ncomp) == 2:
            s = '%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f'% (self.xc,self.xc_err,self.yc,self.yc_err,self.mag,self.mag_err, self.re, self.re_err, self.nsersic, self.nsersic_err, self.BA, self.BA_err, self.PA, self.PA_err, self.sky, self.sky_err, self.error, self.chi2nu) + ' ' + str(self.galname) + '\n' + '%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f'% (self.xc2,self.xc2_err,self.yc2,self.yc2_err,self.mag2,self.mag2_err, self.re2, self.re2_err, self.nsersic2, self.nsersic2_err, self.BA2, self.BA2_err, self.PA2, self.PA2_err, self.sky, self.sky_err, self.error,self.chi2nu) + ' ' + str(self.galname) + '\n'
        
        
        if int(self.ncomp) == 3:
            s = '%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f'% (self.xc,self.xc_err,self.yc,self.yc_err,self.mag,self.mag_err, self.re, self.re_err, self.nsersic, self.nsersic_err, self.BA, self.BA_err, self.PA, self.PA_err, self.sky, self.sky_err, self.error, self.chi2nu) + ' ' + str(self.galname) + '\n' + '%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f'% (self.xc2,self.xc2_err,self.yc2,self.yc2_err,self.mag2,self.mag2_err, self.re2, self.re2_err, self.nsersic2, self.nsersic2_err, self.BA2, self.BA2_err, self.PA2, self.PA2_err, self.sky, self.sky_err, self.error,self.chi2nu) + ' ' + str(self.galname) + '\n' + '%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f'% (self.xc3,self.xc3_err,self.yc3,self.yc3_err,self.mag3,self.mag3_err, self.re3, self.re3_err, self.nsersic3, self.nsersic3_err, self.BA3, self.BA3_err, self.PA3, self.PA3_err, self.sky, self.sky_err, self.error,self.chi2nu) + ' ' + str(self.galname) + '\n'
        
        
        if int(self.ncomp) == 4:
            s = '%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f'% (self.xc,self.xc_err,self.yc,self.yc_err,self.mag,self.mag_err, self.re, self.re_err, self.nsersic, self.nsersic_err, self.BA, self.BA_err, self.PA, self.PA_err, self.sky, self.sky_err, self.error, self.chi2nu) + ' ' + str(self.galname) + '\n' + '%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f'% (self.xc2,self.xc2_err,self.yc2,self.yc2_err,self.mag2,self.mag2_err, self.re2, self.re2_err, self.nsersic2, self.nsersic2_err, self.BA2, self.BA2_err, self.PA2, self.PA2_err, self.sky, self.sky_err, self.error,self.chi2nu) + ' ' + str(self.galname) + '\n' + '%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f'% (self.xc3,self.xc3_err,self.yc3,self.yc3_err,self.mag3,self.mag3_err, self.re3, self.re3_err, self.nsersic3, self.nsersic3_err, self.BA3, self.BA3_err, self.PA3, self.PA3_err, self.sky, self.sky_err, self.error,self.chi2nu) + ' ' + str(self.galname) + '\n' + '%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f'% (self.xc4,self.xc4_err,self.yc4,self.yc4_err,self.mag4,self.mag4_err, self.re4, self.re4_err, self.nsersic4, self.nsersic4_err, self.BA4, self.BA4_err, self.PA4, self.PA4_err, self.sky, self.sky_err, self.error,self.chi2nu) + ' ' + str(self.galname) + '\n'    
        
        
        output.write(s)
        output.close()
    def display_galfit_model(self,percentile1=.5,percentile2=99.5,p1residual=5,p2residual=99,cmap='viridis',zoom=None):
      '''
      ARGS:
      percentile1 = min percentile for stretch of image and model
      percentile2 = max percentile for stretch of image and model
      p1residual = min percentile for stretch of residual
      p2residual = max percentile for stretch of residual
      cmap = colormap, default is viridis
      '''
      # model name
      self.filename = self.galname+'-unwise-'+'w'+str(self.band)+'-1Comp-galfit-out.fits'
      pngname = self.galname+'-unwise-'+'w'+str(self.band)+'-1Comp-galfit-out.png'
      if self.convflag:
         pngname = self.galname+'-unwise-'+'w'+str(self.band)+'-1Comp-galfit-out-conv.png'
      image,h = fits.getdata(self.filename,1,header=True)
      model = fits.getdata(self.filename,2)
      residual = fits.getdata(self.filename,3)

      if zoom is not None:
         print("who's zoomin' who?")
         # display central region of image
         # figure out how to zoom

         # get image dimensions and center
         xmax,ymax = image.shape
         xcenter = int(xmax/2)
         ycenter = int(ymax/2)

         # calculate new size to display based on zoom factor
         new_xradius = int(xmax/2/(float(zoom)))
         new_yradius = int(ymax/2/(float(zoom)))

         # calculate pixels to keep based on zoom factor
         x1 = xcenter - new_xradius
         x2 = xcenter + new_xradius
         y1 = ycenter - new_yradius
         y2 = ycenter + new_yradius
         
         # check to make sure limits are not outsize image dimensions
         if (x1 < 1):
            x1 = 1
         if (y1 < 1):
            y1 = 1
         if (x2 > xmax):
            x2 = xmax
         if (y2 > ymax):
            y2 = ymax

         # cut images to new size
         image = image[x1:x2,y1:y2]
         model = model[x1:x2,y1:y2]
         residual = residual[x1:x2,y1:y2]         
         pass
      wcs = WCS(h)
      images = [image,model,residual]
      titles = ['image','model','residual']
      v1 = [scoreatpercentile(image,percentile1),
            scoreatpercentile(image,percentile1),
            scoreatpercentile(residual,p1residual)]
      v2 = [scoreatpercentile(image,percentile2),
            scoreatpercentile(image,percentile2),
            scoreatpercentile(residual,p2residual)]
      norms = [simple_norm(image,'asinh',max_percent=percentile2),
               simple_norm(image,'asinh',max_percent=percentile2),
               simple_norm(residual,'linear',max_percent=p2residual)]
               
      plt.figure(figsize=(14,6))
      plt.subplots_adjust(wspace=.0)
      for i,im in enumerate(images): 
         plt.subplot(1,3,i+1,projection=wcs)
         plt.imshow(im,origin='lower',cmap=cmap,vmin=v1[i],vmax=v2[i],norm=norms[i])
         plt.xlabel('RA')
         if i == 0:
            plt.ylabel('DEC')
         else:
            ax = plt.gca()
            ax.set_yticks([])
         plt.title(titles[i],fontsize=16)
      plt.savefig(pngname)
      plt.close()
        
      return xcenter, ycenter


        
    def print_galfit_results(self):
      #self.filename = self.galname+'-unwise-'+'w'+str(self.band)+'-1Comp-galfit-out.fits'
      self.filename = self.gal1.output_image
      rg.print_galfit_results(self.filename)
        
    def run_dmc(self, N=100,convflag=True):
        '''
        GOAL: 
        * Run galfit with monte carlo sampling to find all minima

        OPTIONAL PARAMS:
        * N = number of times to run galfit; default is 100
        * convflag = convolve galfit model with psf; default is True

        OUTPUT: 
        * galfit models for all possible minima

        '''
        #N is Number of random samples
        
        # set up arrays to store galfit output (e.g. xf, yf, rf, etc)
        # 7 output parameters, so Nx7 dimensional arrays
        X = np.empty((0,7))
        
        # uncertainty in fitted parameters
        dX = np.empty((0,7))
        C = np.empty((0)) #create list C (charge) with one vector (not array)

        # download the wise images if the user requests this
        if args.getwise:
            self.get_wise_image()

        # define image names
        self.set_image_names()

        # get the pixel coordinates of the galaxy
        # this uses the image header to translate RA and DEC into pixel coordinates
        self.getpix()

        # set up all of the inputs for galfit
        self.initialize_galfit(convflag=convflag)
        B = 1
        for i in range(N):
            E = 100000
            # this loop selects initial conditions that 
            while(np.random.random()>= np.exp(-B*E)): 
                nX = len(X) # number of unique fits that we already found
                D = np.zeros((1,7))             
                self.set_sersic_params() # select random initial conditions
                X0 = np.array([[self.xc,self.yc,self.mag,self.re,self.nsersic,self.BA,self.PA]])             
                for k in range(nX-1):
                    R = X0-X[k,:]
                    R7 = np.linalg.norm(R)**7
                    deltaD = C[k]*R/R7
                    D = D+deltaD
                E = np.linalg.norm(D)**2
            self.run_galfit_wise(fitBA=1,fitPA=1)
            self.get_galfit_results()
            # append best-fit values and errors to array         
            Qnew=np.array([[self.xc,self.yc,self.mag,self.re,self.nsersic,self.BA,self.PA]])
            dQnew = np.array([[self.xc_err,self.yc_err,self.mag_err,self.re_err,self.nsersic_err,self.BA_err,self.PA_err]])

            # initialize flag
            # the flag is 0 if it does not overlap the list of charges
            # the flag is 1 is it overlaps some other charge
            flag=0    
                    
            for j in range(nX): # loop over fitted variables (xc, yc, mag, re, nsersic, BA, PA)
                #compare the difference in fitted parameter with error in difference
                if (np.linalg.norm(Qnew-X[j,:])<= np.linalg.norm(dQnew+dX[j,:])): 
                    # overlaps
                    flag=1
                    C[j]=C[j]+1 # adding a charge
                    if np.linalg.norm(dQnew)<np.linalg.norm(dX[j,:]):
                        X[j,:]=Qnew
                        dX[j,:]=dQnew
                    break

            if flag==0: # if it does not overlap, then add it to the list
                if self.error==0:
                    X=np.append(X,Qnew,axis=0)
                    dX=np.append(dX,dQnew,axis=0)
                    C = np.append(C,[1],axis=0)
        return X    

    def run_simple(self, convflag=True,zoom=1,sersic_start=None):
        '''
        GOAL: 
        * Run galfit once 

        OPTIONAL PARAMS:
        * convflag = convolve galfit model with psf; default is True

        OUTPUT: 
        * galfit model

        '''
        self.convflag = convflag
        # download the wise images if the user requests this
        self.get_wise_image()

        # define image names
        self.set_image_names()

        # get the pixel coordinates of the galaxy
        # this uses the image header to translate RA and DEC into pixel coordinates
        self.getpix()

        # set up all of the inputs for galfit
        self.initialize_galfit(convflag=convflag)
        if sersic_start is not None:
           xc,yc,mag,re,n,BA,PA =  sersic_start
           self.set_sersic_manual(m=mag,re=re,n=n,BA=BA,PA=PA) # set fixed initial parameters
        self.run_galfit_wise(fitBA=1,fitPA=1)
        self.write_results(printflag=True)
        #self.get_galfit_f
        self.display_galfit_model(zoom=zoom)

       ###########################################################################################################################



def readfile(filename):
    
    fileobj=open(filename,'r')                               #opens the file, of course
    words=fileobj.read().splitlines()                        #splits the header and data into two lists
    header = words[0].split()                                #splits header, data arrays into elements
    data = words[1].split()
    fileobj.close()
    return header,data



def readfile2(filename):                                     
    fileobj=open(filename,'r')
    words=fileobj.read().splitlines()
    header = words[0].split()
    fileobj.close()
    return header    
    


def run_galfit_no_psf(galaxy_sample,WISE_dir,sample_txt_name_nopsf):
    
    homedir = os.getenv('HOME')
    os.chdir(homedir+'/github/'+str(WISE_dir))
    get_ipython().run_line_magic('run', '~/github/virgowise/wisesize.py')

    for n in range(0,len(galaxy_sample)):
        
        try:
            vfid = galaxy_sample['VFID'][n]
            g = galaxy(galaxy_sample['RA'][n], galaxy_sample['DEC'][n],
                      galaxy_sample['radius'][n], name = galaxy_sample['prefix'][n],vfid=galaxy_sample['VFID'][n],band='3')
            print(galaxy_sample['prefix'][n])
            ###
            g.set_sersic_manual()
            ###
            g.run_simple(convflag=False)
            
            #text file with output parameters
            t = homedir+'/github/'+str(WISE_dir)+'/'+galaxy_sample[n]['prefix']+'-unwise-w3-log.txt'
            header,data = readfile(t)
            header.pop(0)                                    #removes the pound_sign from the array
            header.append('prefix')
            header.append('success_flag')
            
            for i in range(0,len(data)):
                data[i] = float(data[i])
            data.append(galaxy_sample[n]['prefix'])
            data.append(1)                                   #success_flag value of one
            
            if n == 0:                                       #if the galaxy is the first entry, then
                file_test = [header,data]                    #append to the list both the header & data lists
                file_plots = [header,data]                   #append to list for corner plots
                
            else:
                file_test2 = [header,data]                   #otherwise, only include the data list
                file_test.append(file_test2[1])
                file_plots.append(file_test2[1])
                    
        except:
            
            t = homedir+'/github/'+WISE_dir+'/'+str(galaxy_sample['prefix'][n])+'-unwise-w3-log.txt'
            header = readfile2(t)
            header.pop(0)                                    #removes the pound_sign from the array
            header.append('prefix')
            header.append('success_flag')
            
            data = []
            for num in range(0,len(header)-2):
                data.append(-999)
            data.append(galaxy_sample[n]['prefix'])
            data.append(0)                                   #success_flag value of zero
               
            file_test2 = [header,data]
            file_test.append(file_test2[1])

            
            print(galaxy_sample['prefix'][n], ' ' , 'was unsuccessful.')
            continue
        
    data_array = np.array(file_test)
    data_array_plots = np.array(file_plots)
    np.savetxt(sample_txt_name_nopsf+'.txt',data_array,fmt="%s")                          #all
    np.savetxt(sample_txt_name_nopsf+'_cornerplots.txt',data_array_plots,fmt="%s")        #for corner plots

    
    
def run_galfit_psf(galaxy_sample,WISE_dir,sample_txt_name_nopsf,sample_txt_name_psf):
    
    os.chdir(homedir+'/github/'+str(WISE_dir))
    get_ipython().run_line_magic('run', '~/github/virgowise/wisesize.py')
    
    for n in range(0,len(galaxy_sample)):
        
        try:
            g = galaxy(galaxy_sample['RA'][n], galaxy_sample['DEC'][n],
                      galaxy_sample['radius'][n], name = galaxy_sample['prefix'][n],band='3')
            print(galaxy_sample['prefix'][n])
            
            ###
            tab = ascii.read(homedir+'/github/'+WISE_dir+'/'+sample_txt_name_nopsf+'.txt')
            sersic_parameters = [tab['xc'][n],tab['yc'][n],tab['mag'][n],tab['re'][n],
                                 tab['nsersic'][n],tab['BA'][n],tab['PA'][n]]
            ###
            
            g.run_simple(convflag=True,sersic_start = sersic_parameters)
            
            t = homedir+'/github/'+WISE_dir+'/'+galaxy_sample[n]['prefix']+'-unwise-w3-log.txt'
            header,data = readfile(t)
            header.pop(0)                                    #removes the pound_sign from the array
            header.append('prefix')
            header.append('success_flag')
            
            for i in range(0,len(data)):
                data[i] = float(data[i])
            data.append(galaxy_sample[n]['prefix'])
            data.append(1)                                   #success_flag value of one
            
            if n == 0:                                       #if the galaxy is the first entry, then
                file_test = [header,data]                    #append to the list both the header & data lists
                file_plots = [header,data]                   #append to list for corner plots
                
            else:
                file_test2 = [header,data]                   #otherwise, only include the data list
                file_test.append(file_test2[1])
                file_plots.append(file_test2[1])
                    
        except:
            
            t = homedir+'/github/'+WISE_dir+'/'+galaxy_sample[n]['prefix']+'-unwise-w3-log.txt'
            header = readfile2(t)
            header.pop(0)                                    #removes the pound_sign from the array
            header.append('prefix')
            header.append('success_flag')
            
            data = []
            for num in range(0,len(header)-2):
                data.append(-999)
            data.append(galaxy_sample[n]['prefix'])
            data.append(0)                                   #success_flag value of zero
               
            file_test2 = [header,data]
            file_test.append(file_test2[1])

            
            print(galaxy_sample['prefix'][n], ' ' , 'was unsuccessful.')
            continue
        
    data_array = np.array(file_test)
    data_array_plots = np.array(file_plots)
    np.savetxt(sample_txt_name_psf+'.txt',data_array,fmt="%s")                          #all
    np.savetxt(sample_txt_name_psf+'_cornerplots.txt',data_array_plots,fmt="%s")        #for corner plots
    
    
 
if __name__ == '__main__':
    homedir = os.getenv("HOME")
    #t_full = Table.read(homedir+'/Desktop/v1-20210308/vf_north_v1_main.fits')
    #obj1 = t_full[481]
    #g1 = galaxy(obj1['RA'],obj1['DEC'],obj1['radius'],name=obj1['prefix'],band=3)

    
    
    
    
    
    
    
    
    
    