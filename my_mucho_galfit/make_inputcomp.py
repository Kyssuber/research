'''
GOAL:
- generate galfit input scripts for sample galaxies

INPUT:
- python make_inputcomp.py -param_file /path/to/paramfile/paramfile (no quotation marks)
- paramfile should be params_make_inputcomp*.txt, where the * is replaced with the band

OUTPUT:
- individual input scripts for each central galaxy in sample
- can directly feed these text files into GALFIT
- directly edit params_inputcomp_notrim.txt to tweak parameters as desired
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
import sys

class galfit:
    def __init__(self, galname=None, vfid=None, vfid_v1=None, r25=None, band=None, directory=None, gal_output_path=None, psf_filepath=None, image=None, sigma_image=None, psf_image=None, psf_oversampling=None, mask_image=None, xminfit=None, yminfit=None, xmaxfit=None, ymaxfit=None, convolution_size=None, magzp=None, pscale=None, galfile=None, convflag=None, constraintflag=1, fitallflag=0, ncomp=1, xobj=None, yobj=None, mag=None, rad=None, nsersic=None, BA=None, PA=None, fitmag=1, fitcenter=1, fitrad=1, fitBA=1, fitPA=1, fitn=1, first_time=0, asymmetry=0):
        
        self.galname=galname
        self.vfid=vfid
        self.vfid_v1=vfid_v1
        
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
        
        self.directory = directory
        
        try:
            os.chdir(self.directory+str(self.vfid_v1))
            path_to_input_im = self.directory+str(self.vfid_v1)+'/'
        except:
            print('self.directory did not work. defaulting to /mnt/astrophysics/wisesize/VFIDxxxx (xxxx is VFID_V1 code)')
            print(self.directory)
            os.chdir('/mnt/astrophysics/wisesize/'+str(self.vfid_v1))
            path_to_input_im = '/mnt/astrophysics/wisesize/'+str(self.vfid_v1)+'/'
        
        im = glob.glob('*'+str(self.band)+'-img-m.fits')[0]
        self.image = path_to_input_im+im
        im_mask = glob.glob('*-m-mask.fits')[0]
        self.mask_image = path_to_input_im+im_mask

        sigma_image = glob.glob('*-std-m.fits')[0]
        self.sigma_image = path_to_input_im+sigma_image
        
        try:
            invvar_image = glob.glob('*-inv-mask.fits')[0]
            self.invvar_image = path_to_input_im+invvar_image
        except:
            print('no invvar mask')
        
        #determining arcsec-to-pixel conversion factor from im header information
        
        hdu1=fits.open(self.image)
        header = hdu1[0].header
        #use header information to generalize conversion from arcseconds to pixels (for 12-micron, 2.75 arcsec per px),
        len_image_arcsec = np.abs(header['NAXIS1']*header['CD1_1'])*3600
        arcsec_per_pixel = len_image_arcsec/header['NAXIS1']       
        
        #determining image shape
        
        temp = fits.getdata(self.image)
        self.ximagesize, self.yimagesize = temp.shape
        
        #default psf:
        #os.system('cp '+homedir+'/github/virgowise/wise_psfs/wise-w3-psf-wpro-09x09-05x05.fits .')
        #self.psf_image = 'wise-w3-psf-wpro-09x09-05x05.fits'

        #'personalized' w3-band psfs according to coadd_id
        #copies w3 psf directory into gal_output...there will be as many as there are galaxies in the vf sample, so be prepared for an influx (pun unintended) of point spread functions.

        self.gal_output_path = gal_output_path
        self.psf_filepath = psf_filepath
        
        try: 
            os.chdir(self.gal_output_path)
        except:
            print('self.gal_output_path did not work. defaulting to /mnt/astrophysics/kconger_wisesize/github/gal_output')
            os.chdir('/mnt/astrophysics/kconger_wisesize/github/gal_output/')
        
        try:
            os.system('cp '+self.psf_filepath+str(self.vfid)+'* .')
        except:
            print('self.psf_filepath did not work')
            os.system('cp '+homedir+'/github/virgowise/sgacut_psfs/'+str(self.vfid)+'* .')
        
        psf_image = glob.glob(str(self.vfid)+'*psf.fits')[0]
        self.psf_image = self.gal_output_path+psf_image
        
        #value from original script
        self.psf_oversampling=8
        
        #for WISE
        self.magzp=magzp     #magnitude zeropoint
        self.pscale=pscale
        
        self.convflag=convflag
        
        self.xobj=self.ximagesize/2
        self.yobj=self.yimagesize/2
        self.mag=mag
        self.rad=rad
        self.nsersic=nsersic
        self.BA=BA
        self.PA=PA
        self.fitmag=fitmag
        self.fitn=fitn
        self.fitcenter=fitcenter
        self.fitrad=fitrad
        self.fitBA=fitBA
        self.fitPA=fitPA
        
        # default max and min fits
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
        r25_px = self.r25/arcsec_per_pixel     #see about 70 lines above for conversion factor
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

    def write_image_params(self): 
        
        self.galfit_input.write('# IMAGE PARAMETERS\n')
        self.galfit_input.write('A) '+self.image+'              # Input data image (FITS file)\n')
        self.galfit_input.write('B) '+self.output_image+'       # Name for the output image\n')
        self.galfit_input.write('C) %s                # Sigma image name (made from data if blank or "none") \n'%(self.sigma_image))
        if self.convflag==1:
            self.galfit_input.write('D) '+self.psf_image+'     # Input PSF image and (optional) diffusion kernel\n')
            self.galfit_input.write('E) %i                   # PSF oversampling factor relative to data\n'%(self.psf_oversampling))
            
        self.galfit_input.write('F) %s                # Bad Pixel Mask (i.e., obj mask; blank or "none") \n'%(self.mask_image))
         
         #if self.constraintflag:
        #    self.galfit_input.write('G) /Users/rfinn/research/LocalClusters/GalfitAnalysis/sersic.constraint        # Parameter constraint file (ASCII)\n')   
            
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
    
    if '-h' in sys.argv or '--help' in sys.argv:
        print("Usage: %s [-param_file (name of parameter file, no single or double quotation marks)]")
        sys.exit(1)
    
    if '-param_file' in sys.argv:
        p = sys.argv.index('-param_file')
        param_file = str(sys.argv[p+1]) 
    
    homedir = os.getenv("HOME")
       
    #create dictionary with keywords and values, from parameter .txt file

    param_dict = {}
    with open(homedir+'/'+param_file) as f:
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
    directory = param_dict['directory']
    psf_filepath = param_dict['psf_filepath']
    psf_filepath = homedir+psf_filepath
    gal_output_path = param_dict['gal_output_path']
    
    #parameters
    band = param_dict['band']
    magzp = float(param_dict['magzp'])
    convflag = int(param_dict['convflag'])
    asymmetry = int(param_dict['asymmetry'])
    pscale = float(param_dict['pscale'])
    psf_oversampling = int(param_dict['psf_oversampling'])
    
    #initial guesses
    mag = int(param_dict['mag'])
    rad = int(param_dict['rad'])
    nsersic = int(param_dict['nsersic'])
    BA = int(param_dict['BA'])
    PA = int(param_dict['PA'])
    
    #fit flags
    fitmag = int(param_dict['fitmag'])
    fitrad = int(param_dict['fitrad'])
    fitn = int(param_dict['fitn'])
    fitBA = int(param_dict['fitBA'])
    fitPA = int(param_dict['fitPA'])
    fitcenter = int(param_dict['fitcenter'])
    
    print('convflag:', str(convflag))
    
    for i in range(0,len(cat)):
    
        gal = galfit(galname=cat['prefix'][i], vfid=cat['VFID'][i], vfid_v1=cat['VFID_V1'][i], 
                     r25 = cat['radius'][i], band=band, directory=directory, gal_output_path=gal_output_path, 
                     psf_filepath=psf_filepath, psf_oversampling=psf_oversampling, magzp=magzp, pscale=pscale, 
                     convflag=convflag, constraintflag=1, fitallflag=0, ncomp=1, mag=mag, rad=rad, nsersic=nsersic, 
                     BA=BA, PA=PA, fitmag=fitmag, fitcenter=fitcenter, fitrad=fitrad, fitBA=fitBA, fitPA=fitPA, 
                     fitn=fitn, first_time=0, asymmetry=asymmetry)

        #if nopsf already run, then use output params as initial input parameter guesses
        #otherwise, use default guesses entered when creating gal class
        
        if gal.convflag == 1:
            print('Checking for noPSF parameter table at',gal_output_path+'output_params_'+gal.band+'_nopsf.fits','...')
            try:
                params = Table.read(gal.gal_output_path+'output_params_'+gal.band+'_nopsf.fits')
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
            except:
                print('ERROR: noPSF parameter file not found. defaulting to initial guesses in param_file.')
        
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
                gal.add_simple_sersic_object(n,'sersic',x=x,y=y,mag=gal.mag,rad=gal.rad,nsersic=gal.nsersic,BA=gal.BA,PA=gal.PA)
                
                if n == int(len(indices)+1):
                    gal.write_sky(n+1)
        else:
            gal.write_sky(2)
        if (gal.fitallflag):
            gal.fitall()
        gal.close_input_file()
        
