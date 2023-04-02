#!/usr/bin/env python

'''
GOAL: 
- this program will run galfit on one galaxy
- galfit is run once with generic input parameters, no convolution
- the output is used as input for a second run, w/convolution enabled

USAGE:
This should be run from a directory that contains
- galfit input file
- image
- noise image

To run:
python ~/github/virgowise/run1galfit.py galname

- galname should be the VFID and should correspond to the subdirectories in 
/mnt/astrophysics/rfinn/muchogalfit-output (grawp directory)

/mnt/astrophysics/muchogalfit-output (how directory is mounted on virgo vms)

- each subdirectory has an input file that is created by setup_galfit.py

'''
import os
import sys
#import glob
import numpy as np
### DICTIONARIES

pixel_scale = {'FUV':1.5,'NUV':1.5,'g':0.262,'r':0.262,'z':0.262,'W1':2.75,'W2':2.75,'W3':2.75,'W4':2.75}
psf_oversampling = {'FUV':1,'NUV':1,'g':1,'r':1,'g':1,'W1':1,'W2':1,'W3':1,'W4':1}
mag_zeropoint = {'FUV':22.5,'NUV':22.5,'g':22.5,'r':22.5,'g':22.5,'W1':22.5,'W2':22.5,'W3':22.5,'W4':22.5}
image_resolution = {'FUV':6,'NUV':6,'g':1.5,'r':1.5,'z':1.5,'W1':6.1,'W2':6.4,'W3':6.5,'W4':12}

# set up a dictionary for the radius to use for the first guess of the sersic profile
# a better way is to use a constant angular size, and then we can translate that into pixels using the pixel scale

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
            maskname = image.replace(l+".fits","r-mask.fits")
            return maskname
        
def write_galfit_input(galdir, output_dir, objname, ra, dec, bandpass, firstpass=True):
    galname = objname
    #print('inside write_galfit_input: ',galdir,galname)
    image = f'{galname}-custom-image-{bandpass}.fits.fz'
    invvar_image = f'{galname}-custom-invvar-{bandpass}.fits.fz'    
    psf_image = f'{galname}-custom-psf-{bandpass}.fits.fz'

    # created images
    sigma_image = f'{galname}-custom-std-{bandpass}.fits'    
    if firstpass:
        output_image = f'{galname}-{bandpass}-out1.fits'
    else:
        output_image = f'{galname}-{bandpass}-out2.fits'    

    
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
    
    # TODO: add mask to galfit input
    # have updated mask wrapper in halpha gui
    maskfound = False
    mask_image = get_maskname(image)
    if os.path.exists(mask_image):
        maskfound = True
        print("found a mask - will implement masking in galfit")

    # make of values for xminfit, etc for now
    # get image size

    # TODO: need to get xmaxfit,ymaxfit
    xmax,ymax = get_image_size(image) # am I mixing x and y dimensions here?
    print('image dimensions = ',xmax,ymax)
    xminfit = 1
    xmaxfit = xmax 
    yminfit = 1
    ymaxfit = ymax
    
    #if firstpass:
    #    # read in image header and convert RA, DEC to xpixel,ypixel
    

    if firstpass:
        # DONE: need to get (x,y) center of object
        xobj, yobj = get_xy_from_wcs(ra,dec,image)
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

    else:
        # read in output from first pass run of galfit
        t = parse_galfit_output(output_image.replace('out2','out1'))
        print(t)
        # header_keywords=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_SKY','CHI2NU']
        xc, yc, mag, rad, nsersic, BA, PA, sky = t[0][0],t[1][0],t[2][0],t[3][0],t[4][0],t[5][0],t[6][0],t[7][0]
        fitBA = 1
        fitPA = 1
        fitn = 1
        fitmag = 1
        fitrad = 1
        xobj, yobj = xc, yc
        # get convolution size - set to cutout size?
    # TODO: need to decide on the right size for convolution - could be 2x image size?
    # Chien had recommended the full image, but a smaller convolution size should make galfit run faster

    # galfit website FAQ says 40-80 x PSF seeing disk
    # see: https://users.obs.carnegiescience.edu/peng/work/galfit/TFAQ.html#convbox
    convolution_size = int(image_resolution[bandpass]/pixel_scale[bandpass]*50)
    if firstpass:
        outfile = open('galfit.input1','w')
    else:
        outfile = open('galfit.input2','w')        
    outfile.write('# IMAGE PARAMETERS\n')
    outfile.write('A) '+image+'              # Input data image (FITS file)\n')
    outfile.write('B) '+output_image+'       # Name for the output image\n')
    outfile.write('C) %s                # Sigma image name (made from data if blank or "none") \n'%(sigma_image))
    if not firstpass:
        outfile.write('D) '+psf_image+'     # Input PSF image and (optional) diffusion kernel\n')
        outfile.write('E) %i                   # PSF oversampling factor relative to data\n'%(psf_sampling))
    if maskfound:
        outfile.write(f'F) {mask_image}     # Input PSF image and (optional) diffusion kernel\n')
    outfile.write('H) '+str(int(round(xminfit)))+' '+str(int(round(xmaxfit)))+' '+str(int(round(yminfit)))+' '+str(int(round(ymaxfit)))+'     # Image region to fit (xmin xmax ymin ymax)\n')
    if not firstpass:
        outfile.write('I) '+str(int(round(convolution_size)))+' '+str(int(round(convolution_size)))+'             # Size of convolution box (x y)\n')
    outfile.write('J) %5.2f              # Magnitude photometric zeropoint \n'%(magzp))
    outfile.write('K) %6.5f   %6.5f         # Plate scale (dx dy)  [arcsec/pix]\n'%(pscale,pscale))
    outfile.write('O) regular                # Display type (regular, curses, both)\n')
    outfile.write('P) 0                   # Create output image only? (1=yes; 0=optimize) \n')
    outfile.write('S) 0                   # Modify/create objects interactively?\n')
    outfile.write(' \n')
    # write object
    # github test
    outfile.write('# Object number: 1 \n')
    outfile.write(' 0) sersic             # Object type \n')
    outfile.write(' 1) %8.1f  %8.1f 1 1  # position x, y        [pixel] \n'%(xobj,yobj))
    outfile.write(' 3) %5.2f      %i       # total magnitude     \n'%(mag,fitmag))
    outfile.write(' 4) %8.2f       %i       #     R_e              [Pixels] \n'%(rad,fitrad))
    outfile.write(' 5) %5.2f       %i       # Sersic exponent (deVauc=4, expdisk=1)   \n'%(nsersic,int(fitn)))
    outfile.write(' 9) %5.2f       %i       # axis ratio (b/a)    \n'%(BA,int(fitBA)))
    outfile.write('10) %5.2f       %i       # position angle (PA)  [Degrees: Up=0, Left=90] \n'%(PA,int(fitPA)))
    outfile.write(" Z) 0                  # Output option (0 = residual, 1 = Don't subtract)  \n")
    

    outfile.write(' \n')
    outfile.write('# Object number: 2 \n')
    outfile.write(' 0) sky             # Object type \n')
    outfile.write(' 1) %8.1f   1  # sky background at center of fitting region [ADUs] \n'%(sky))
    outfile.write(' 2) 0      0       # dsky/dx (sky gradient in x)    \n')
    outfile.write(' 3) 0      0       # dsky/dy (sky gradient in y) \n')
    outfile.write(" Z) 0                  # Output option (0 = residual, 1 = Don't subtract)  \n")
                  
    
    outfile.close()


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
    data_dir = '/mnt/astrophysics/virgofilaments-data/{}/{}/'.format(int(ra),objname)



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
