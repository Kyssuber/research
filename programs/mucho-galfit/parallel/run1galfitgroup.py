#!/usr/bin/env python

'''
GOAL: 
- this program will run galfit on one of JM's group images
- galfit is run once with generic input parameters, no convolution
- the output is used as input for a second run, w/convolution enabled

INPUT:
- group 

PROCEDURE:
- load virgo catalog (vf_v2_main.fits)
- get virgo galaxies in FOV
- create input file for galfit, adding sersic profile at each location
- run galfit w/out convolution

- parse output from galfit, getting details on all sources
- create new input for galfit, with all sources and convolution
- run galfit with convolution


USAGE:
This should be run from a directory that contains
- galfit input file
- image
- noise image

To run:
python ~/github/virgowise/run1galfitgroup.py galname 

- galname should be the VFID and should correspond to the subdirectories in 
/mnt/astrophysics/rfinn/muchogalfit-output (grawp directory)

/mnt/astrophysics/muchogalfit-output (how directory is mounted on virgo vms)

- each subdirectory has an input file that is created by setup_galfit.py

'''
import os
import sys
#import glob
import numpy as np

from astropy.io import fits
from astropy.wcs import WCS

homedir = os.getenv("HOME")
# add in masking from halphagui
sys.path.append(homedir+'/github/halphagui/')
from maskwrapper import buildmask
import reproject_mask

### DICTIONARIES

pixel_scale = {'FUV':1.5,'NUV':1.5,'g':0.262,'r':0.262,'z':0.262,'W1':2.75,'W2':2.75,'W3':2.75,'W4':2.75}
psf_oversampling = {'FUV':1,'NUV':1,'g':1,'r':1,'g':1,'W1':1,'W2':1,'W3':1,'W4':1}
mag_zeropoint = {'FUV':22.5,'NUV':22.5,'g':22.5,'r':22.5,'g':22.5,'W1':22.5,'W2':22.5,'W3':22.5,'W4':22.5}
image_resolution = {'FUV':6,'NUV':6,'g':1.5,'r':1.5,'z':1.5,'W1':6.1,'W2':6.4,'W3':6.5,'W4':12}
minmag2fit = {'FUV':10,'NUV':10,'g':17,'r':17,'z':17,'W1':10,'W2':10,'W3':10,'W4':10}

# set up a dictionary for the radius to use for the first guess of the sersic profile
# a better way is to use a constant angular size, and then we can translate that into pixels using the pixel scale

class buildgroupmask(buildmask):
    def __init__(self,image):
        self.image_name = image
        self.image, self.imheader = fits.getdata(self.image_name,header = True)
        self.image_wcs = WCS(self.imheader)
        self.ymax,self.xmax = self.image.shape
        
        # SE parameters for masking
        self.threshold=0.005
        self.snr=10
        self.snr_analysis=10        
        self.sepath = homedir+'/github/halphagui/astromatic/'
        self.config='default.sex.HDI.mask'

        # set auto to true to run masking in automatic mode
        self.auto = True
        
        # GAIA catalog, with path set for grawp
        self.gaiapath = '/mnt/astrophysics/rfinn/catalogs/gaia-mask-dr9.virgo.fits'
        self.gaia_mask = None
        self.add_gaia_stars = True        

        t = self.image_name.split('.fit')
        self.mask_image=t[0]+'-mask.fits'
        self.mask_inv_image=t[0]+'-inv-mask.fits'
        print('saving mask image as: ',self.mask_image)

        # SET UP AND RUN SOURCE EXTRACTOR
        self.link_files()
        self.runse()
        self.add_gaia_masks()
        self.write_mask()
        
        self.grow_mask()
        self.grow_mask()
        self.grow_mask()


    def remove_gals(self,xgals,ygals):
        for x,y in zip(xgals,ygals):
            # get mask value at location of galaxy
            maskval = self.maskdat[int(y),int(x)]
            
            
            # get pixels where the mask value = value at location of galaxy

            # and 
            # set those pixels to zero
            # adding commit just because
            self.maskdat[self.maskdat == maskval]=0
        self.write_mask()
        
### FUNCTIONS
def funpack_image(input,output,nhdu=1):
    command = 'funpack -O {} {}'.format(output,input)
    print(command)
    os.system(command)
def funpack_image_astropy(input,output,nhdu=1):
    from astropy.io import fits
    hdu = fits.open(input)
    print('input file = ',input)
    #print(hdu.info())
    fits.writeto(output,data=hdu[nhdu].data, header=hdu[nhdu].header, overwrite=True)
    hdu.close()
    print('finished unpacking image')
    
def parse_galfit_output(galfit_outimage):
    from astropy.io import fits
    #import numpy as np
    numerical_error_flag=0
    header_keywords=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_SKY','CHI2NU']
    fit_parameters=[]
    working_dir=os.getcwd()+'/'
    image_header = fits.getheader(galfit_outimage,2)
    for hkey in header_keywords:
        s=str(image_header[hkey])
        #print hkey
        if s.find('[') > -1:
            s=s.replace('[','')
            s=s.replace(']','')
            t=s.split('+/-')
            values=(float(t[0]),0.)# fit and error
        else:
            t=s.split('+/-')
            try:
                values=(float(t[0]),float(t[1]))# fit and error
            except ValueError:
                # look for * in the string, which indicates numerical problem
                if t[0].find('*') > -1:
                    numerical_error_flag=1
                    t[0]=t[0].replace('*','')
                    t[1]=t[1].replace('*','')
                    values=(float(t[0]),float(t[1]))# fit and error
            except IndexError: # for CHI2NU
                chi2nu=float(t[0])
                continue
        fit_parameters.append(values)
    fit_parameters.append(numerical_error_flag)
    fit_parameters.append(chi2nu)
    return fit_parameters

def convert_invvar_noise(invvar_image, noise_image):
    # convert invvar image to noise
    from astropy.io import fits
    import numpy as np
    # read in invvar image
    print('invvar image = ',invvar_image, os.path.basename(invvar_image))
    hdu = fits.open(invvar_image)
    data = hdu[0].data
    header = hdu[0].header
    hdu.close()
    # operate on pixels to take sqrt(1/x)
    noise_data = np.sqrt(1/data)
    
    # check for bad values, nan, inf
    # set bad pixels to very high value, like 1e6
    noise_data = np.nan_to_num(noise_data,nan=1.e6)

    # write out as noise image
    fits.writeto(noise_image,noise_data,header=header,overwrite=True)

def get_image_size(image):
    from astropy.io import fits
    data = fits.getdata(image)
    data_shape = data.shape
    return data_shape

def get_xy_from_wcs(ra,dec,image):
    from astropy.io import fits
    from astropy.wcs import WCS
    from astropy.coordinates import SkyCoord
    header = fits.getheader(image)
    w = WCS(header)
    c = SkyCoord(ra,dec,unit='deg',frame='icrs')
    xobj,yobj = w.world_to_pixel(c)
    return xobj,yobj

def get_galaxies_in_fov(image,bandpass='W3'):
    """get galaxies in FOV using virgo catalog """
    from astropy.wcs import WCS
    from astropy.coordinates import SkyCoord
    from astropy.table import Table

    # read in virgo catalog
    catalog='/mnt/astrophysics/rfinn/catalogs/Virgo/v2/vf_v2_main.fits'
    try:
        vtab = Table.read(catalog)
    except FileNotFoundError:
        # test to see if running on Virgo VMS
        catalog='/mnt/astrophysics/catalogs/Virgo/v2/vf_v2_main.fits'
        vtab = Table.read(catalog)
    # create a SkyCoord object from RA and DEC of virgo galaxies
    galcoord = SkyCoord(vtab['RA'],vtab['DEC'],frame='icrs',unit='deg')

    # set up image wcs
    image_wcs = WCS(image)

    # get the size of the image
    xmax, ymax = get_image_size(image)
    
    # find galaxies on cutout
    x,y = image_wcs.world_to_pixel(galcoord)

    # create flag to save galaxies on the image
    flag = (x > 0) & (x < xmax) & (y>0) & (y < ymax)        
    vfids = vtab['VFID'][flag]
    x,y = x[flag],y[flag]
    # write out file containing VFID, x, y
    ofilename = f'galsFOV-{bandpass}.txt'
    outfile = open(ofilename,'w')
    for i in range(len(vfids)):
        outfile.write(f'{vfids[i]}, {x[i]:.2f}, {y[i]:.2f} \n')
    outfile.close()
    return x,y

def get_group_mask(image,xgals=None,ygals=None):
    """create the mask for the group image, remove galaxies to be fitted  """
    # create the mask
    # remove objects at the positions of the galaxies to be fitted
    
    m = buildgroupmask(image)
    m.remove_gals(xgals,ygals)
    maskname = get_maskname(image)
    os.rename(m.mask_image,maskname)
    # return mask_image
    return maskname

def get_maskname(image):
    """ return the wise mask for wise images and the r mask for legacy images  """
    # check for wise extenstion
    elist = ['W1','W2','W3','W4']
    
    for e in elist:
        if e in image:
            maskname = image.replace(e+".fits","wise-mask.fits")
            return maskname
    # check for legacy
    llist = ['image-g','image-r','image-z']
    for l in llist:
        if l in image:
            maskname = image.replace(l+".fits","image-r-mask.fits")
            return maskname
        
def write_galfit_input(galdir, output_dir, objname, ra, dec, bandpass, firstpass=True):
    galname = objname
    #print('inside write_galfit_input: ',galdir,galname)
    image = f'{galname}_GROUP-custom-image-{bandpass}.fits.fz'
    invvar_image = f'{galname}_GROUP-custom-invvar-{bandpass}.fits.fz'    
    psf_image = f'{galname}_GROUP-custom-psf-{bandpass}.fits.fz'

    # created images
    sigma_image = f'{galname}_GROUP-custom-std-{bandpass}.fits'    
    if firstpass:
        output_image = f'{galname}_GROUP-{bandpass}-out1.fits'
    else:
        output_image = f'{galname}_GROUP-{bandpass}-out2.fits'    

    
    if firstpass:

        
        # funpack the .fz images so galfit can read them
        # save them in output_dir
        funpack_image(os.path.join(galdir,image),os.path.join(output_dir,image.replace('.fz','')))
        funpack_image(os.path.join(galdir,psf_image),os.path.join(output_dir,psf_image.replace('.fz','')),nhdu=0)
        # unpack invvar image
        funpack_image(os.path.join(galdir,invvar_image),os.path.join(output_dir,invvar_image.replace('.fz','')),nhdu=0)        
    # get information from the image, like the image size and position of gal
    image = image.replace('.fz','')
    psf_image = psf_image.replace('.fz','')
    invvar_image = invvar_image.replace('.fz','')
    # prepend output directory to all images
    sigma_image = output_dir+'/'+sigma_image
    invvar_image = output_dir+'/'+invvar_image
    # check if noise image exists, if not make it from invvar    
    if not os.path.exists(sigma_image):
        convert_invvar_noise(invvar_image,sigma_image)
    sigma_image = os.path.basename(sigma_image)

        
    psf_sampling = psf_oversampling[bandpass]
    pscale = pixel_scale[bandpass]
    magzp = mag_zeropoint[bandpass]
    minmag = minmag2fit[bandpass]    
    # TODO: add mask to galfit input
    # have updated mask wrapper in halpha gui
    maskfound = False

    # skipping masking now for group images
    mask_image = get_maskname(image)

    

    # TODO: need to get xmaxfit,ymaxfit
    xmax,ymax = get_image_size(image) # am I mixing x and y dimensions here?
    print('image dimensions = ',xmax,ymax)
    xminfit = 1
    xmaxfit = xmax 
    yminfit = 1
    ymaxfit = ymax
    
    #if firstpass:
    #    # read in image header and convert RA, DEC to xpixel,ypixel
    
    convolution_size = int(image_resolution[bandpass]/pixel_scale[bandpass]*20)
    if firstpass:
        # remove an straggler galfit files
        # this ensures that the first pass results are in galfit.01
        os.system('rm galfit.??')
        # DONE: need to get (x,y) center of object
        xobj, yobj = get_xy_from_wcs(ra,dec,image)
        xgal, ygal = get_galaxies_in_fov(image,bandpass=bandpass)
        if not os.path.exists(mask_image):
            if bandpass in ['W1','W2','W3','W4']:
                rmask = mask_image.replace(bandpass,'r')
                if os.path.exists(rmask):
                    #reproject mask
                    reproject_mask.reproject_image(rmask,image)
            else:
                mask_image = get_group_mask(image,xgals=xgal,ygals=ygal)
            maskfound = True
        BA=1
        fitBA = 1
        PA=0
        fitPA = 1
        nsersic=2
        fitn = 1
        mag=7
        fitmag = 1
        sky = 0
        # set initial guess to 25 arcsec,
        # and translate into pixels based on pixelscale of bandpass 
        rad = 15/pixel_scale[bandpass]
        fitrad = 1

        
        outfile = open('galfit.input1','w')
        outfile.write('# IMAGE PARAMETERS\n')
        outfile.write('A) '+image+'              # Input data image (FITS file)\n')
        outfile.write('B) '+output_image+'       # Name for the output image\n')
        outfile.write('C) %s                # Sigma image name (made from data if blank or "none") \n'%(sigma_image))
        if maskfound:
            outfile.write(f'F) {mask_image}     # Input PSF image and (optional) diffusion kernel\n')
        outfile.write('H) '+str(int(round(xminfit)))+' '+str(int(round(xmaxfit)))+' '+str(int(round(yminfit)))+' '+str(int(round(ymaxfit)))+'     # Image region to fit (xmin xmax ymin ymax)\n')
        outfile.write('J) %5.2f              # Magnitude photometric zeropoint \n'%(magzp))
        outfile.write('K) %6.5f   %6.5f         # Plate scale (dx dy)  [arcsec/pix]\n'%(pscale,pscale))
        outfile.write('O) regular                # Display type (regular, curses, both)\n')
        outfile.write('P) 0                   # Create output image only? (1=yes; 0=optimize) \n')
        outfile.write('S) 0                   # Modify/create objects interactively?\n')
        outfile.write(' \n')
        # write object
        # github test
        objnum = 1
        for xi, yi in zip(xgal,ygal):
            outfile.write(f'# Object number: {objnum} \n')
            outfile.write(' 0) sersic             # Object type \n')
            outfile.write(' 1) %8.1f  %8.1f 1 1  # position x, y        [pixel] \n'%(xi,yi))
            outfile.write(' 3) %5.2f      %i       # total magnitude     \n'%(mag,fitmag))
            outfile.write(' 4) %8.2f       %i       #     R_e              [Pixels] \n'%(rad,fitrad))
            outfile.write(' 5) %5.2f       %i       # Sersic exponent (deVauc=4, expdisk=1)   \n'%(nsersic,int(fitn)))
            outfile.write(' 9) %5.2f       %i       # axis ratio (b/a)    \n'%(BA,int(fitBA)))
            outfile.write('10) %5.2f       %i       # position angle (PA)  [Degrees: Up=0, Left=90] \n'%(PA,int(fitPA)))
            outfile.write(" Z) 0                  # Output option (0 = residual, 1 = Don't subtract)  \n")
            objnum += 1

        outfile.write(' \n')
        outfile.write(f'# Object number: {objnum} \n')
        outfile.write(' 0) sky             # Object type \n')
        outfile.write(' 1) %8.1f   1  # sky background at center of fitting region [ADUs] \n'%(sky))
        outfile.write(' 2) 0      0       # dsky/dx (sky gradient in x)    \n')
        outfile.write(' 3) 0      0       # dsky/dy (sky gradient in y) \n')
        outfile.write(" Z) 0                  # Output option (0 = residual, 1 = Don't subtract)  \n")
                  
    
        outfile.close()
    else:
        # read in output from first pass

        
        input = open('galfit.01','r')
        all_lines = input.readlines()
        outlines = []
        

        i = 0
        holdfixed = False
        skyobject = True
        for i in range(len(all_lines)):
            
            line = all_lines[i]
            
            if line.startswith('B)'):
                outlines.append(line.replace('out1.fits','out2.fits'))
                holdfixed = False
                
            elif line.startswith('D)'):
                
                outlines.append('D) '+psf_image+'     # Input PSF image and (optional) diffusion kernel\n')
            elif line.startswith('E)'):
                outlines.append('E) %i                   # PSF oversampling factor relative to data\n'%(psf_sampling))
            elif line.startswith('I)'):
                outlines.append('I) '+str(int(round(convolution_size)))+' '+str(int(round(convolution_size)))+'             # Size of convolution box (x y)\n')

            elif line.startswith(' 0)'):
                if line.find('sky') > -1:
                    skyobject = True
                else:
                    skyobject = False
                outlines.append(line)
            
            
                
            elif line.startswith(' 3)'):# check if mag is too faint, then don't fit for new parameters
                t = line.split()
                # check the magnitude returned from no covolution
                #print(f"checking magnitude {float(t[1])} compared to {minmag}",holdfixed)
                if float(t[1]) > minmag and not skyobject:
                    # keep the positions fixed
                    #print(i,len(outlines))
                    outlines[i-2] = outlines[i-2].replace(' 1 1',' 0 0')
                                
                    # set all parameters to fixed
                    # fix the x and y coords
                                
                                
                                
                    holdfixed=True
                    outlines.append(line.replace(' 1 ',' 0 '))
                else:
                    holdfixed=False
                    outlines.append(line)
                #print(f"checking magnitude {float(t[1])} compared to {minmag}",holdfixed)                
            elif (line.startswith(' 4)') or line.startswith(' 5)') or line.startswith(' 9)') or line.startswith('10)')):
                if holdfixed and not skyobject:
                    outlines.append(line.replace(' 1 ',' 0 '))
                else:
                    outlines.append(line)
            else:
                outlines.append(line)
                
        outfile = open('galfit.input2','w')
        for l in outlines:
            outfile.write(l)
        outfile.close()
        input.close()
                                                                                                  
if __name__ == '__main__':
    
    # run this from /mnt/astrophysics

    # DONE: move galfit output to a new destination
    # /mnt/astrophysics/rfinn/muchogalfit-output
    topdir = '/mnt/astrophysics/rfinn/muchogalfit-output/'
    try:
        os.chdir(topdir)
    except FileNotFoundError: # assuming that we are running on virgo vms
        topdir = '/mnt/astrophysics/muchogalfit-output/'
        os.chdir(topdir)
    # take as input the galaxy name
    galname = sys.argv[1]

    # move to muchogalfit-output directory
    output_dir = topdir+galname+'/'
    if not os.path.exists(output_dir):
        print('WARNING: {} does not exist\n Be sure to run setup_galfit.py first')
        os.chdir(topdir)
        sys.exit()
    
    os.chdir(output_dir)

    # define environment variable so funpack can find the correct variables
    os.environ["LD_LIBRARY_PATH"]="/opt/ohpc/pub/compiler/gcc/9.4.0/lib64:/home/siena.edu/rfinn/software/cfitsio-4.2.0/"


    
    # read in input file to get galname, objname, ra, dec, and bandpass
    sourcefile = open(galname+'sourcelist','r')
    galaxies = sourcefile.readlines()
    if len(galaxies) > 1:
        # set the flag to have more than one galaxy in the galfit input file
        multiflag = True
    elif len(galaxies) == 1:
        multiflag = False
    else:
        print('Problem reading sourcelist for {}'.format(galname))
        print('Please check the setup directory')
        os.chdir(topdir)
        sys.exit()

    # parse information from file
    vfid, objname, ra, dec, bandpass = galaxies[0].rstrip().split()
    ra = float(ra)
    dec = float(dec)

    # reset bandpass to the wavelength that is passed in from the command line
    bandpass = sys.argv[2]
    

    # set up path name for image directory
    # directory where galaxy images are
    data_dir = '/mnt/astrophysics/virgofilaments-data/{}/{}_GROUP/'.format(int(ra),objname)



    # TODO: add code to generate galfit input for first run, no convolution, generic starting point
    write_galfit_input(data_dir, output_dir, objname, ra, dec, bandpass)
    
    # code to run galfit
    print('running galfit')
    os.system(f"galfit galfit.input1")


    # TODO: read galfit output, and create new input to run with convolution
    write_galfit_input(data_dir, output_dir, objname, ra, dec, bandpass, firstpass=False)

    # TODO: skipping convolution for now, so that I can test parallel code.  Come back to this.
    # TODO: make sure I am using the correct PSF images
    print('running galfit second time')
    os.system(f"galfit galfit.input2")

    os.chdir(topdir)
