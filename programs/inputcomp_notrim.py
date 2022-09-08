


'''
GOAL:
- generate galfit input scripts for sample galaxies
OUTPUT:
- individual scripts for each central galaxy in sample; default output directory is /mnt/astrophysics/kconger_wisesize/github/gal_output
- can directly feed these text files into GALFIT
- Line 31 --> galaxy catalog
- Line 32 --> dummy catalog
- Line 82 --> output directory
- Line 84 --> psf directory
'''

import numpy as np
import warnings
warnings.filterwarnings('ignore')
import os
homedir = os.getenv("HOME")
import glob
from astropy.table import Table
from astropy.wcs import WCS
from astropy.io import ascii
from astropy.io import fits

#os.sys.path.append('/mnt/astrophysics/kconger_wisesize/github/gal_output')
#os.chdir('/mnt/astrophysics/kconger_wisesize/github/gal_output')

#have to call sample
cat = Table.read(homedir+'/sgacut_coadd.fits')
dummycat = Table.read(homedir+'/dummycat.fits')


class galfit:
    def __init__(self, galname=None, vfid=None, vfid_v1=None, r25=None, band=3, image=None, sigma_image=None, psf_image=None, psf_oversampling=None, mask_image=None, xminfit=None, yminfit=None, xmaxfit=None, ymaxfit=None, convolution_size=None, magzp=None, pscale=None, galfile=None, convflag=1, constraintflag=1, fitallflag=0, ncomp=1, xobj=None, yobj=None, mag=None, rad=None, nsersic=None, BA=None, PA=None, fitmag=1, fitcenter=1, fitrad=1, fitBA=1, fitPA=1, fitn=1, first_time=0, asymmetry=0):
        
        self.galname=galname
        self.vfid=vfid
        self.vfid_v1=vfid_v1
        
        #if band is an integer, than prepend with w to indicate that the 'band' is a WISE channel
        try:
            int(band)
            self.band = 'w'+str(band)
        except:
            self.band = band
        
        self.constraintflag=constraintflag
        self.fitallflag=fitallflag
        
        if vfid in dummycat['central galaxy']:
            self.ncomp = len(np.where(dummycat['central galaxy'] == vfid)[0]) + 1
        else:
            self.ncomp=ncomp

        self.asymmetry=0
        self.galfile=str(self.galname)+'-galfit.input.'+str(self.ncomp)+'Comp'

        #self.image_rootname = self.galname+'-unwise-'+str(self.band)

        #change directory, grab galaxy's image and mask FITS filenames
        #NOTE: WISESIZE DIRECTORIES ARE LABELED USING V1 VFIDs!
        os.chdir('/mnt/astrophysics/wisesize/'+str(self.vfid_v1))
        im = glob.glob('*'+str(self.band)+'-img-m.fits')[0]
        self.image = im
        im_mask = glob.glob('*mask.fits')[0]
        self.mask_image = im_mask

        self.sigma_image = glob.glob('*-std-m.fits')[0]
        
        #I think
        try:
            self.invvar_image = glob.glob('*-inv-mask.fits')[0]
        except:
            print('')
        
        temp = fits.getdata(self.image)
        self.ximagesize, self.yimagesize = temp.shape
        
        #default psf:
        #os.system('cp '+homedir+'/github/virgowise/wise_psfs/wise-w3-psf-wpro-09x09-05x05.fits .')
        #self.psf_image = 'wise-w3-psf-wpro-09x09-05x05.fits'

        #'personalized' w3-band psfs according to coadd_id
        #copies w3 psf directory into gal_output...there will be as many as there are galaxies in the vf sample, so be prepared for an influx (pun unintended) of point spread functions.
        os.chdir('/mnt/astrophysics/kconger_wisesize/github/gal_output/')
        os.system('cp '+homedir+'/github/virgowise/sgacut_psfs/'+str(self.vfid)+'* .')
        self.psf_image = glob.glob(str(self.vfid)+'*psf.fits')[0]
        
        #value from original script
        self.psf_oversampling=8
        
        #for WISE
        self.magzp=22.5     #magnitude zeropoint
        self.pscale=0.344
        
        self.convflag=convflag

        self.xobj=self.ximagesize/2
        self.yobj=self.yimagesize/2
        self.mag=7
        self.rad=5
        self.nsersic=2
        self.BA=1
        self.PA=0
        self.fitmag=fitmag
        self.fitn=fitn
        self.fitcenter=fitcenter
        self.fitrad=fitrad
        self.fitBA=fitBA
        self.fitPA=fitPA
        
        #default max and min fits
        #self.xmaxfit=self.ximagesize
        #self.ymaxfit=self.yimagesize
        #self.xminfit=0
        #self.yminfit=0
        #self.convolution_size=self.ximagesize
        
        
        #each cutout from Dustin is 500x500 px; the aim here is instruct GALFIT to only model a certain region about this central galaxy (which I quasi-arbitrarily choose to be a cutout size of 3*d25)
        xc = self.ximagesize/2
        yc = self.yimagesize/2
        self.r25 = r25
        #convert arcseconds to pixels
        r25_px = self.r25/2.75
        size = r25_px * 3
        self.xminfit = xc - size
        self.yminfit = yc - size
        self.xmaxfit = xc + size
        self.ymaxfit = yc + size
        self.convolution_size = self.xmaxfit - self.xminfit


    def create_output_names(self):
        if self.asymmetry == 1:
            self.output_image=str(self.galname)+'-'+ str(self.band)+'-'+ str(self.ncomp)+'Comp-galfit-out-asym.fits'
        elif self.convflag == 1:
            self.output_image=str(self.galname)+'-'+ str(self.band)+'-'+ str(self.ncomp) +'Comp-galfit-out-conv.fits'
        else:
            self.output_image=str(self.galname)+'-'+ str(self.band)+'-'+ str(self.ncomp) +'Comp-galfit-out.fits'
            

    def open_galfit_input(self):
        self.galfit_input=open(self.galfile,'w')


    def write_image_params(self): #,input_image,output_image,sigma_image,psf_image,psf_oversampling,mask_image,xminfit,xmaxfit,yminfit,ymaxfit,convolution_size,magzp,pscale,convflag=1,constraintflag=1,fitallflag=0):
        self.galfit_input.write('# IMAGE PARAMETERS\n')
        self.galfit_input.write('A) '+self.image+'              # Input data image (FITS file)\n')
        self.galfit_input.write('B) '+self.output_image+'       # Name for the output image\n')
        self.galfit_input.write('C) %s                # Sigma image name (made from data if blank or "none") \n'%(self.sigma_image))
        if self.convflag==1:
            self.galfit_input.write('D) '+self.psf_image+'     # Input PSF image and (optional) diffusion kernel\n')
            self.galfit_input.write('E) %i                   # PSF oversampling factor relative to data\n'%(self.psf_oversampling))
        self.galfit_input.write('H) '+str(int(round(self.xminfit)))+' '+str(int(round(self.xmaxfit)))+' '+str(int(round(self.yminfit)))+' '+str(int(round(self.ymaxfit)))+'     # Image region to fit (xmin xmax ymin ymax)\n')
        if self.convflag==1:
            self.galfit_input.write('I) '+str(int(round(self.convolution_size)))+' '+str(int(round(self.convolution_size)))+'             # Size of convolution box (x y)\n')
        self.galfit_input.write('J) %5.2f              # Magnitude photometric zeropoint \n'%(self.magzp))
        self.galfit_input.write('K) %6.5f   %6.5f         # Plate scale (dx dy)  [arcsec/pix]\n'%(self.pscale,self.pscale))
        self.galfit_input.write('O) regular                # Display type (regular, curses, both)\n')
        self.galfit_input.write('P) 0                   # Create output image only? (1=yes; 0=optimize) \n')
        self.galfit_input.write('S) 0                   # Modify/create objects interactively?\n')

        
    def set_sky(self,sky):
        self.sky=sky


    def write_sersic(self,objnumber,profile):

        self.galfit_input.write(' \n')
        self.galfit_input.write('# Object number: %i \n'%(objnumber))
        self.galfit_input.write(' 0) %s             # Object type \n'%(profile))
        self.galfit_input.write(' 1) %8.1f  %8.1f %i %i  # position x, y        [pixel] \n'%(self.xobj,self.yobj,int(self.fitcenter),int(self.fitcenter)))
        self.galfit_input.write(' 3) %5.2f      %i       # total magnitude     \n'%(self.mag,self.fitmag))
        self.galfit_input.write(' 4) %8.2f       %i       #     R_e              [Pixels] \n'%(self.rad,self.fitrad))
        if self.nsersic>5:
            self.nsersic=5
            self.fitn=0
        self.galfit_input.write(' 5) %5.2f       %i       # Sersic exponent (deVauc=4, expdisk=1)   \n'%(self.nsersic,int(self.fitn)))
        self.galfit_input.write(' 9) %5.2f       %i       # axis ratio (b/a)    \n'%(self.BA,int(self.fitBA)))
        self.galfit_input.write('10) %5.2f       %i       # position angle (PA)  [Degrees: Up=0, Left=90] \n'%(self.PA,int(self.fitPA)))
        if self.asymmetry:
            self.galfit_input.write('F1) 0.0001 0.00   1  1     # azim. Fourier mode 1, amplitude & phase angle \n')
        self.galfit_input.write(" Z) 0                  # Output option (0 = residual, 1 = Don't subtract)  \n")
        

    def write_sky(self,objnumber):
        self.galfit_input.write(' \n')
        self.galfit_input.write('# Object number: %i \n'%(objnumber))
        self.galfit_input.write(' 0) sky             # Object type \n')
        self.galfit_input.write(' 1) %8.1f   1  # sky background at center of fitting region [ADUs] \n'%(self.sky))
        self.galfit_input.write(' 2) 0      0       # dsky/dx (sky gradient in x)    \n')
        self.galfit_input.write(' 3) 0      0       # dsky/dy (sky gradient in y) \n')
        self.galfit_input.write(" Z) 0                  # Output option (0 = residual, 1 = Don't subtract)  \n")


    def add_simple_sersic_object(self,objnumber,profile,x,y,mag,rad,nsersic,BA,PA):
        self.galfit_input.write(' \n')
        self.galfit_input.write('# Object number: %i \n'%(objnumber))
        self.galfit_input.write(' 0) %s             # Object type \n'%(profile))
        self.galfit_input.write(' 1) %8.1f  %8.1f 1 1  # position x, y        [pixel] \n'%(x,y))
        self.galfit_input.write(' 3) %5.2f      1       # total magnitude     \n'%(mag))
        self.galfit_input.write(' 4) %8.2f       1       #     R_e              [Pixels] \n'%(rad))
        self.galfit_input.write(' 5) %5.2f       1       # Sersic exponent (deVauc=4, expdisk=1)   \n'%(nsersic))
        self.galfit_input.write(' 9) %5.2f       1       # axis ratio (b/a)    \n'%(BA))
        self.galfit_input.write('10) %5.2f       1       # position angle (PA)  [Degrees: Up=0, Left=90] \n'%(PA))
        self.galfit_input.write(" Z) 0                  # Output option (0 = residual, 1 = Don't subtract)  \n")


    def close_input_file(self):
        self.galfit_input.close()


if __name__ == '__main__':
    
    convflag = input('conv? enter 0 (n) or 1 (y): ')
    #band = input('enter band (1-4 for WISE, r for r-band: ')
    
    if int(convflag) == 1:
        #be sure path to parameters table is correct!
        params = Table.read('/mnt/astrophysics/kconger_wisesize/gal_output/output_params_'+gal.band+'_nopsf.fits')
    
    for i in range(0,len(cat)):
    
        gal = galfit(galname=cat['prefix'][i], vfid=cat['VFID'][i], vfid_v1=cat['VFID_V1'][i], r25 = cat['radius'][i], convflag=convflag, constraintflag=1, fitallflag=0, ncomp=1)
        
        #gal.band=band
         
        #if nopsf already run, then use output params as initial input parameter guesses
        if gal.convflag == 1:
            ind = np.where(cat['galname'][i] == params['galname'])[0]
            gal.xobj=params['xc'][ind]
            gal.yobj=params['yc'][ind]
            gal.mag=params['mag'][ind]
            #help prevent errors associated with unphysical nopsf output params
            if int(params['nsersic'][ind])>5:
                gal.nsersic=5 #restrict nser initial guess to n=5
                gal.rad=5     #revert to default initial guess for Re
            else:
                gal.nsersic=params['nsersic'][ind]
                gal.rad=params['re'][ind]
            gal.BA=params['BA'][ind]
            gal.PA=params['PA'][ind]
        
        print(gal.vfid)
        gal.create_output_names()
        gal.set_sky(0)
        gal.open_galfit_input()
        gal.write_image_params()
        gal.write_sersic(1,'sersic')
        
        if gal.vfid in dummycat['central galaxy']:
            indices = np.where(dummycat['central galaxy'] == gal.vfid)[0]
            for i in range(0,len(indices)):
                index = indices[i]
                n = int(i)+2      #i begins at 0, and 1 is already taken
                x = int(dummycat['x (px)'][index])
                y = int(dummycat['y (px)'][index])
                gal.add_simple_sersic_object(n,'sersic',x=x,y=y,mag=7,rad=5,nsersic=2,BA=1,PA=0)
                
                if n == int(len(indices)+1):
                    gal.write_sky(n+1)
        else:
            gal.write_sky(2)
        if (gal.fitallflag):
            gal.fitall()
        gal.close_input_file()
