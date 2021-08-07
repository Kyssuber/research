# coding: utf-8


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
import warnings
warnings.filterwarnings('ignore',category=DeprecationWarning)
from astropy.wcs import WCS
from astropy.table import Table
from astropy.io import ascii
from matplotlib import pyplot as plt

from scipy.stats import scoreatpercentile
import wget
import tarfile
import glob
import gzip


import astropy.wcs as wcs


homedir = os.getenv("HOME")

os.sys.path.append(homedir+'/github/virgowise/')

dummycat = Table.read(homedir+'/dummycat.fits',format='ascii')

sga_params = Table.read(homedir+'/sga_vf_matched.fits')

sgacut = Table.read(homedir+'/sga_cut.fits',format='ascii')



import rungalfit as rg #This code has all the defined functions that I can use
os.sys.path.append(homedir+'/github/HalphaImaging/python3/')
import plot_cutouts_ha as cutouts #This code has all the defined functions that I can use




vf = Table.read(homedir+'/vfcut.fits',format='ascii')
dummycat = Table.read(homedir+'/dummycat.fits',format='ascii')



       
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
        #variable names for input galaxy
        self.ra = ra#*u.deg
        self.dec = dec#*u.deg
        self.radius = size#*u.arcsec
        self.band = band
        self.galname = name
        self.vfid = vfid
        self.image_rootname = self.galname+'-unwise-w'+str(self.band)
        self.image = self.image_rootname+'-img-m-trim.fits'

        #change directory, grab galaxy's image and mask FITS files
        os.chdir('/mnt/astrophysics/wisesize/'+str(self.vfid))
        im = glob.glob('*w3-img-m-trim.fits')[0]
        im_mask = glob.glob('*mask-trim.fits')[0]
        self.mask_image = im_mask

        self.sigma_image = self.image_rootname+'-std-m-trim.fits'
        self.invvar_image = self.image_rootname+'-invvar-m-trim.fits'

        #change to galfit output directory
        os.chdir(homedir+'/github/gal_output')
        
        # remove the log file from output dir if it exists
        self.logfilename = self.galname+'-unwise-'+'w'+str(self.band)+'-log.txt'
        if os.path.exists(self.logfilename):
           os.remove(self.logfilename)

        # write the header line into the log file
        output=open(self.logfilename,'w')
        output.write('# xc xc_err yc yc_err mag mag_err re re_err nsersic nsersic_err BA BA_err PA PA_err sky sky_err error chi2nu galname \n')
        # close log file
        output.close()    






   def get_wise_image(self,makeplots=False):
      #vfmain = Table.read(homedir+'/vf_north_v1_main.fits')
      #vf defined above
      vfmain = vf
      base_dir = '/mnt/astrophysics/wisesize/'
      #loop through VFID column in vfmain; if the ID matches that of the galaxy defined above, then grab the appropriate files
      for i in vfmain['VFID']:
         if str(self.vfid) == i:
            #os.chdir('/mnt/astrophysics/wisesize/')
            self.image = base_dir + str(self.vfid) + '/unwise-' + str(self.vfid) + '-w3-img-m-trim.fits'
            self.sigma_image = base_dir + str(self.vfid) + '/unwise-' + str(self.vfid) + '-w3-std-m-trim.fits'
            temp = fits.getdata(self.image)
            print(temp.shape)
            self.ximagesize, self.yimagesize = temp.shape
            break
         else:
            continue


   

        
        
   def set_image_names(self):
        '''
        GOAL:
        * Set the psf image name and some additional parameters for galfit 

        PARAMS:
        * none

        OUTPUT: PSF image and optional diffusion kernel, other parameters  

        '''
        dummycat=Table.read(homedir+'/dummycat.fits',format='ascii')
        vfid = self.galname[0:8]

        #just using center til, doesn't matter usually
        os.chdir(homedir+'/github/gal_output/')
        os.system('cp '+homedir+'/github/virgowise/wise_psfs/wise-w3-psf-wpro-09x09-05x05.fits .') 
        self.psf_image = 'wise-w3-psf-wpro-09x09-05x05.fits' 
        self.psf_oversampling = 8

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


        if vfid in dummycat['central galaxy']:
            self.ncomp = len(np.where(dummycat['central galaxy'] == vfid)[0]) + 1
        else:
            self.ncomp=1

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
        self.mag = m
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
        
        if self.vfid in sgacut['VFID']:
            sgaindex = np.where(sga_params['VFID'] == self.vfid)[0]

            BA = sga_params['BA'][sgaindex]
            PA = sga_params['PA'][sgaindex]
            fitBA = 0
            fitPA = 0
        else:
           BA = self.BA
           PA = self.PA
           

        
        #os.system('cp '+self.psf_image+' .')
        self.gal1.set_sersic_params(xobj=self.xc,yobj=self.yc,mag=self.mag,rad=self.re,nsersic=self.nsersic,BA=BA,PA=PA,fitmag=1,fitcenter=1,fitrad=1,fitBA=fitBA,fitPA=fitPA,fitn=1,first_time=0)
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
        t = rg.parse_galfit_1comp(self.filename,ncomp=self.ncomp)
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
        # create string with best-fit parameter

        if int(self.ncomp) == 1:    
            s = '%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f'%(self.xc,self.xc_err,self.yc,self.yc_err,self.mag,self.mag_err, self.re, self.re_err, self.nsersic, self.nsersic_err, self.BA, self.BA_err, self.PA, self.PA_err, self.sky, self.sky_err, self.error,self.chi2nu) + ' ' + str( self.galname)+'\n'
        


        
        if int(self.ncomp) == 2:
           dummycat = Table.read(homedir+'/dummycat.fits',format='ascii')
           vfid = self.galname[0:8]
           index = np.where(dummycat['central galaxy'] == vfid)[0]
           galname2 = dummycat['ID'][index[0]]
           print(galname2)



           s = '%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f'% (self.xc,self.xc_err,self.yc,self.yc_err,self.mag,self.mag_err, self.re, self.re_err, self.nsersic, self.nsersic_err, self.BA, self.BA_err, self.PA, self.PA_err, self.sky, self.sky_err, self.error, self.chi2nu) + ' ' + str(self.galname) + '\n' + '%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f'% (self.xc2,self.xc2_err,self.yc2,self.yc2_err,self.mag2,self.mag2_err, self.re2, self.re2_err, self.nsersic2, self.nsersic2_err, self.BA2, self.BA2_err, self.PA2, self.PA2_err, self.sky, self.sky_err, self.error,self.chi2nu) + ' ' + str(galname2) + '\n'
        
        
        if int(self.ncomp) == 3:
           dummycat = Table.read(homedir+'/dummycat.fits',format='ascii')
           vfid = self.galname[0:8]
           indices = np.where(dummycat['central galaxy'] == vfid)[0]
           ind1 = indices[0]
           ind2 = indices[1]
           galname2 = dummycat['ID'][ind1]
           galname3 = dummycat['ID'][ind2]



           s = '%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f'% (self.xc,self.xc_err,self.yc,self.yc_err,self.mag,self.mag_err, self.re, self.re_err, self.nsersic, self.nsersic_err, self.BA, self.BA_err, self.PA, self.PA_err, self.sky, self.sky_err, self.error, self.chi2nu) + ' ' + str(self.galname) + '\n' + '%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f'% (self.xc2,self.xc2_err,self.yc2,self.yc2_err,self.mag2,self.mag2_err, self.re2, self.re2_err, self.nsersic2, self.nsersic2_err, self.BA2, self.BA2_err, self.PA2, self.PA2_err, self.sky, self.sky_err, self.error,self.chi2nu) + ' ' + str(galname2) + '\n' + '%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f'% (self.xc3,self.xc3_err,self.yc3,self.yc3_err,self.mag3,self.mag3_err, self.re3, self.re3_err, self.nsersic3, self.nsersic3_err, self.BA3, self.BA3_err, self.PA3, self.PA3_err, self.sky, self.sky_err, self.error,self.chi2nu) + ' ' + str(galname3) + '\n'
        
        
        if int(self.ncomp) == 4:
           dummycat = Table.read(homedir+'/dummycat.fits',format='ascii')
           vfid = self.galname[0:8]
           indices = np.where(dummycat['central galaxy'] == vfid)[0]
           ind1 = indices[0]
           ind2 = indices[1]
           ind3 = indices[2]
           galname2 = dummycat['ID'][ind1]
           galname3 = dummycat['ID'][ind2]
           galname4 = dummycat['ID'][ind3]


           s = '%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f'% (self.xc,self.xc_err,self.yc,self.yc_err,self.mag,self.mag_err, self.re, self.re_err, self.nsersic, self.nsersic_err, self.BA, self.BA_err, self.PA, self.PA_err, self.sky, self.sky_err, self.error, self.chi2nu) + ' ' + str(self.galname) + '\n' + '%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f'% (self.xc2,self.xc2_err,self.yc2,self.yc2_err,self.mag2,self.mag2_err, self.re2, self.re2_err, self.nsersic2, self.nsersic2_err, self.BA2, self.BA2_err, self.PA2, self.PA2_err, self.sky, self.sky_err, self.error,self.chi2nu) + ' ' + str(galname2) + '\n' + '%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f'% (self.xc3,self.xc3_err,self.yc3,self.yc3_err,self.mag3,self.mag3_err, self.re3, self.re3_err, self.nsersic3, self.nsersic3_err, self.BA3, self.BA3_err, self.PA3, self.PA3_err, self.sky, self.sky_err, self.error,self.chi2nu) + ' ' + str(galname3) + '\n' + '%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f'% (self.xc4,self.xc4_err,self.yc4,self.yc4_err,self.mag4,self.mag4_err, self.re4, self.re4_err, self.nsersic4, self.nsersic4_err, self.BA4, self.BA4_err, self.PA4, self.PA4_err, self.sky, self.sky_err, self.error,self.chi2nu) + ' ' + str(galname4) + '\n'


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

      self.filename = self.galname+'-unwise-'+'w'+str(self.band)+'-'+str(self.ncomp)+'Comp-galfit-out.fits'
      if self.convflag:
         self.filename = self.galname+'-unwise-'+'w'+str(self.band)+'-'+str(self.ncomp)+'Comp-galfit-out-conv.fits'
      pngname = self.galname+'-unwise-'+'w'+str(self.band)+'-'+str(self.ncomp)+'Comp-galfit-out.png'
      if self.convflag:
         pngname = self.galname+'-unwise-'+'w'+str(self.band)+'-'+str(self.ncomp)+'Comp-galfit-out-conv.png'
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
         ax = plt.subplot(1,3,i+1,projection=wcs)
         plt.imshow(im,origin='lower',cmap=cmap,vmin=v1[i],vmax=v2[i],norm=norms[i])
         ax.set_xlabel('RA')
         if i == 0:
            ax.set_ylabel('DEC')
         else:
            ax = plt.gca()
            ax.set_yticks([])
         plt.title(titles[i],fontsize=16)
      plt.savefig(pngname)
      plt.close()
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
    words=fileobj.read().splitlines()                        #splits the header and data into lists
    header = words[0].split()       #splits header into elements
    header.remove('#')              #removes pound symbol from header
    header.append('success_flag')   #indicates whether galfit runs on galaxy successfully
    header.append('central_flag')   #indicates whether galaxy is central (1) or off-centered (0)
    data1 = words[1].split()        #splits data1 array into elements
    if len(words) > 2:
       if len(words) == 3:
          data2 = words[2].split()
          data=[header,data1,data2]
       if len(words) == 4:
          data2 = words[2].split()
          data3 = words[3].split()
          data = [header,data1,data2,data3]
       if len(words) == 5:
          data2 = words[2].split()
          data3 = words[3].split()
          data4 = words[4].split()
          data = [header,data1,data2,data3,data4]
    else:
       data = [header,data1]
    fileobj.close()
    return data



def readfile2(filename):                                     
    fileobj=open(filename,'r')
    words=fileobj.read().splitlines()
    header = words[0].split()
    fileobj.close()
    return header




def run_galfit_no_psf(galaxy_sample,WISE_dir,sample_txt_name_nopsf):


    #get_ipython().run_line_magic('run', '~/github/virgowise/wisesize.py')

    txt = '%run ~/github/virgowise/wisesize.py'
    os.system(txt)

    #note: galaxy_sample only includes central galaxies by design
    for n in range(0,len(galaxy_sample)):
    
       vfid = galaxy_sample['VFID'][n]
       g = galaxy(galaxy_sample['RA'][n], galaxy_sample['DEC'][n],
                      galaxy_sample['radius'][n], name = galaxy_sample['prefix'][n],vfid=galaxy_sample['VFID'][n],band='3')
       print(galaxy_sample['prefix'][n])
       
       g.set_sersic_manual()
       
       try:
          g.run_simple(convflag=False)
          t = homedir+'/github/'+str(WISE_dir)+'/'+galaxy_sample[n]['prefix']+'-unwise-w3-log.txt'
          data = readfile(t)      #returns list of header,data*



            
          if n == 0:                              #if the galaxy is the first entry, then...
             for i in range(0,len(data)):
                try:
                   m=i+1                          #define index to avoid data[0], which is the header
                   data[m].append(int(1))
                   name = str(data[m][18])        #isolate VFID
                   name = name[0:8]
                   print(name)
                   if name in vf['VFID']:     #if galname is in central galaxy, append 1; else, append 0
                      data[m].append(int(1))
                   else:
                      data[m].append(int(0))
                except:                           #if the index m is beyond the length of the list, then continue (exit loop)
                   continue             

             file_test = data                     #append to the list both the header & data lists
             file_plots = data                    #append to list for corner plots (which will not include defunct galaxies)
   
          else:
             for i in range(0,len(data)-1):
                if (i+1) <  len(data):
                   num=i+1                          #again, skips i=0 (header) index
                   data[num].append(int(1))
                   name = str(data[num][18])        #isolate VFID
                   name = name[0:8]
                   print(name)

                   if name in vf['VFID']:
                      data[num].append(int(1))
                   else:
                      data[num].append(int(0))

                else:
                   print('fin')
                   continue
                file_test.append(data[num])
                #file_plots.append(data[num])
                    
       except:
            
           t = homedir+'/github/'+WISE_dir+'/'+str(galaxy_sample['prefix'][n])+'-unwise-w3-log.txt'
           header = readfile2(t)
           header.remove('#')                                    #removes the pound_sign from the array
           
           header.append('success_flag')
            
           data = []
           for num in range(0,len(header)-2):
              data.append(-999)
           data.append(galaxy_sample[n]['prefix'])
           data.append(-999)                                   #success_flag value of zero, other entries are -999
           data.append(-999)    
           file_test2 = [header,data]
           file_test.append(file_test2[1])

            
           print(galaxy_sample['prefix'][n], ' failed at run_simple.')
           print(galaxy_sample['prefix'][n], ' failed at run_simple.')
           print(galaxy_sample['prefix'][n], ' failed at run_simple.')
           print(galaxy_sample['prefix'][n], ' failed at run_simple.')
           print(galaxy_sample['prefix'][n], ' failed at run_simple.')
           continue
        
    data_array = np.array(file_test)
    #data_array_plots = np.array(file_plots)
    np.savetxt(sample_txt_name_nopsf+'.txt',data_array,fmt="%s")                           #all
    #np.savetxt(sample_txt_name_nopsf+'_cornerplots.txt',data_array_plots,fmt="%s")        #for corner plots
 
    
    
def run_galfit_psf(galaxy_sample,WISE_dir,sample_txt_name_nopsf,sample_txt_name_psf):
   
    #get_ipython().run_line_magic('run', '~/github/virgowise/wisesize.py')

    txt = '%run ~/github/virgowise/wisesize.py'
    os.system(txt)

    try:
       tab = ascii.read(homedir+'/github/'+WISE_dir+'/'+sample_txt_name_nopsf+'.txt')
    except:
       tab = ascii.read('/mnt/astrophysics/kconger_wisesize/gal_output_nopsf/'+str(sample_txt_name_nopsf)+'.txt')
    central_flag = np.asarray(tab['central_flag'])
    central_flag = central_flag.astype(bool)
    success_flag = np.asarray(tab['success_flag'])
    success_flag = success_flag.astype(bool)
    tab_central = tab[central_flag]

    
    for n in range(0,len(galaxy_sample)):
        
       try:
            g = galaxy(galaxy_sample['RA'][n], galaxy_sample['DEC'][n],
                      galaxy_sample['radius'][n], name = galaxy_sample['prefix'][n],vfid=galaxy_sample['VFID'][n],band='3')
            print(galaxy_sample['prefix'][n])
            
            ###
            
            sersic_parameters = [tab_central['xc'][n],tab_central['yc'][n],tab_central['mag'][n],tab_central['re'][n],
                                 tab_central['nsersic'][n],tab_central['BA'][n],tab_central['PA'][n]]
            name = tab_central['galname']
            ###
            
            g.run_simple(convflag=True,sersic_start = sersic_parameters)
            
            t = homedir+'/github/'+WISE_dir+'/'+galaxy_sample[n]['prefix']+'-unwise-w3-log.txt'
            data = readfile(t)      #returns list of header,data1,...header



            if n == 0:                              #if the galaxy is the first entry, then...
               for i in range(0,len(data)):
                  try:
                     m=i+1                          #define index to avoid data[0], which is the header
                     data[m].append(int(1))
                     name = str(data[m][18])        #isolate VFID
                     name = name[0:8]
                     print(name)
                     if name in vf['VFID']:     #if galname is in vf catalog, append 1; else, append 0
                        data[m].append(int(1))
                     else:
                        data[m].append(int(0))
                  except:                           #if the index m is beyond the length of the list, then continue (exit loop)
                     continue             

               file_test = data                     #append to the list both the header & data lists
               file_plots = data                    #append to list for corner plots (which will not include defunct galaxies)
                
            else:
               for i in range(0,len(data)-1):
                  if (i+1) < len(data):
                     num=i+1                          #again, skips i=0 (header) index
                     data[num].append(int(1))
                     name = str(data[num][18])        #isolate VFID
                     name = name[0:8]
                     print(name)
                     if name in vf['VFID']:
                        data[num].append(int(1))
                     else:
                        data[num].append(int(0))
                  else:
                     print('fin')
                     continue
                  file_test.append(data[num])

                  #file_plots.append(dat)


       except:
            
           t = homedir+'/github/'+WISE_dir+'/'+str(galaxy_sample['prefix'][n])+'-unwise-w3-log.txt'
           header = readfile2(t)
           header.remove('#')                                    #removes the pound_sign from the array
           
           header.append('success_flag')
            
           data = []
           for num in range(0,len(header)-2):
              data.append(-999)
           data.append(galaxy_sample[n]['prefix'])
           data.append(-999)                                   #success_flag value of zero, other entries are -999
           data.append(-999)    
           file_test2 = [header,data]
           file_test.append(file_test2[1])

            
           print(galaxy_sample['prefix'][n], ' failed at run_simple.')
           continue

        
    data_array = np.array(file_test)
    #data_array_plots = np.array(file_plots)
    np.savetxt(sample_txt_name_psf+'.txt',data_array,fmt="%s")                           #all
    #np.savetxt(sample_txt_name_psf+'_cornerplots.txt',data_array_plots,fmt="%s")        #for corner plots
    

if __name__ == '__main__':
   print('run_galfit_no_psf(galaxy_sample,WISE_dir,sample_txt_name_nopsf)')
   print('run_galfit_psf(galaxy_sample,WISE_dir,sample_txt_name_nopsf,sample_txt_name_psf)')

   print('catalogs: vfcut (SNR>10 from VF parent sample), sgacut (vfcut matched with SGA)')

   homedir = os.getenv("HOME")
   dummycat = Table.read(homedir+'/dummycat.fits',format='ascii')
   
