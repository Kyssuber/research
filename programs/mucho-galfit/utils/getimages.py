#!/usr/bin/env python

# THIS PROGRAM IS NOT NEEDED BECAUSE JOHN M. MOVED THE IMAGES TO SIENA!!!

def getlegacy(ra,dec,galid='VFID0',pixscale=1,imsize='60',band='g',makeplots=False,subfolder=None):
    """
    Download legacy image for a particular ra, dec
    
    Inputs:
    * ra
    * dec
    * galid = galaxy id (e.g. VFID0001); used for naming the image files
    * imsize = size of cutout in pixe ls
    * band = filter for the fits images that will be returned.  e.g. 'g' or 'r' or 'z'. 
    * pixscale = pixel scale of cutout in arcsec; native is 0.262 for legacy
    * makeplots = boolean, generate plot of image
    * subfolder = default is None; you can specify a name of a subfolder to 
                  save the data in, e.g., subfolder='legacy-images'
    Returns:
    * fits_name = fits image name
    * jpeg_name = jpeg image name
    """
    imsize = int(imsize)

    # make output image names
    if subfolder is not None:
        # check if subfolder exists. if not, make it.
        if not os.path.exists(subfolder):
            os.mkdir(subfolder)
        rootname = subfolder+'/'+str(galid)+'-legacy-'+str(imsize)
    else:
        rootname = str(galid)+'-legacy-'+str(imsize)        
    jpeg_name = rootname+'.jpg'
    fits_name = rootname+'-'+band+'.fits'


    print('legacy imsize = ',imsize)
    
    # check if images already exist
    # if not download images
    if not(os.path.exists(jpeg_name)):
        print('retrieving ',jpeg_name)
        url='http://legacysurvey.org/viewer/jpeg-cutout?ra='+str(ra)+'&dec='+str(dec)+'&layer=dr8&size='+str(imsize)+'&pixscale='+str(pixscale)
        urlretrieve(url, jpeg_name)
    else:
        print('previously downloaded ',jpeg_name)
    if not(os.path.exists(fits_name)):
        print('retrieving ',fits_name)
        url='http://legacysurvey.org/viewer/cutout.fits?ra='+str(ra)+'&dec='+str(dec)+'&layer=dr8&size='+str(imsize)+'&pixscale='+str(pixscale)+'&bands='+band
        urlretrieve(url, fits_name)
    else:
        print('previously downloaded ',fits_name)

    # try to read the data in
    try:
        t,h = fits.getdata(fits_name,header=True)
        
    except IndexError:
        print('problem accessing image')
        print(fits_name)
        url='http://legacysurvey.org/viewer/cutout.fits?ra='+str(ra1)+'&dec='+str(dec1)+'&layer=dr8&size='+str(image_size)+'&pixscale=1.00'
        print(url)
        return None

    # this will trigger if the image is outside the legacy footprint
    if np.mean(t[1]) == 0:
        return None

    return fits_name, jpeg_name


def getunwise(ra,dec,galid='VFID0',pixscale=2.75,imsize='60',bands='1234',makeplots=False,subfolder=None):
    """
    Download unwise image for a particular ra, dec
    
    Inputs:
    * ra
    * dec
    * galid = galaxy id (e.g. VFID0001)
    * imsize = size of cutout in pixels
    * pixscale = pixel scale of cutout in arcsec
      - native is 0.262 for legacy; 
      - 2.75 for wise
    """

    # TODO: where do we get the reprocessed unwise images?
    
    downloadwise = True
    # check if images already exist
    if subfolder is not None:
        image_names = glob.glob(subfolder+'/'+galid+'-unwise*img-m.fits')
    else:
        image_names = glob.glob(galid+'-unwise*img-m.fits')
    if len(image_names) > 3:
        print('unwise images already downloaded')
        print(image_names)
        # should be only one *-img-m.fits image per band
        if len(image_names) > len(bands):
            multiframe=True
        else:
            multiframe = False
        weight_names = glob.glob(galid+'-unwise*std*.fits')
        if not multiframe:
            return image_names,weight_names,multiframe
        else:
            print('going to try new stacking for wise multiframe')
            downloadwise = False
    if downloadwise:
        imsize = int(imsize)
        print('wise image size = ',imsize)
        baseurl = 'http://unwise.me/cutout_fits?version=allwise'
        imurl = baseurl +'&ra=%.5f&dec=%.5f&size=%s&bands=%s'%(ra,dec,imsize,bands)
        print('downloading unwise images')
        print(imurl)
        wisetar = wget.download(imurl)
        tartemp = tarfile.open(wisetar,mode='r:gz') #mode='r:gz'
        wnames = tartemp.getnames()

        print(wnames)
        # check for multiple pointings - means galaxy is split between images
        multiframe = False
        if len(wnames) > 4*len(bands):
            multiframe = True
    
        wmembers = tartemp.getmembers()
        image_names = []
        weight_names = []
        tartemp.extractall()
        for fname in wnames:
            #print(fname)
            t = fname.split('-')
            if subfolder is not None:
                rename = subfolder+'/'+str(galid)+'-'+t[0]+'-'+t[1]+'-'+t[2]+'-'+t[3]+'-'+t[4]
            else:
                rename = str(galid)+'-'+t[0]+'-'+t[1]+'-'+t[2]+'-'+t[3]+'-'+t[4]
            #print('rename = ',rename)
            #print(rename.find('gz'))
            #if os.path.exists(rename): # this should only occur if multiple images are returned from wise
            #    os.remove(rename)
            os.rename(fname, rename)
            if rename.find('.gz') > -1:
                #print('hello????')
                os.system('gunzip '+rename)
                rename = rename.split('.gz')[0]
                #print('after gunzip, rename = ',rename)
            if rename.find('img') > -1:
                image_names.append(rename)
            if rename.find('std') > -1:
                # move ivar images to imagename.weight.fits
                outim = rename.replace('std-m','img-m.std')
                os.rename(rename,outim)
                weight_names.append(outim)
        os.remove(wisetar)
    
    # if multiframe
    # run swarp to create coadded image
    if multiframe:
        image_names=[]
        weight_names=[]
        for b in bands:
            print('running swarp to combine multiple unwise images in band ',b)
            #########################################
            ## COMBINE THE IMAGE FRAMES USING AVERAGE
            #########################################        
            # create default.swarp
            os.system('swarp -d > default.swarp')
            # run swarp
            matchstring = "*w{}-img-m.fits".format(b)            
            if subfolder is not None:
                allfiles = glob.glob(subfolder+'/'+galid+matchstring)


                
            else:
                allfiles = glob.glob(galid+matchstring)

            all_images = " ".join(allfiles)
            output_image = str(galid)+'-'
            s = 'swarp '+all_images+' -COMBINE_TYPE AVERAGE -WEIGHT_SUFFIX .std.fits -SUBTRACT_BACK N'
            print(s)
            os.system(s)

            # rename coadd.fits to the output image name
            outimage = str(galid)+'-unwise-w'+str(b)+'-coadd.fits'
            if subfolder is not None:
                os.rename('coadd.fits',os.path.join(subfolder,outimage))
                image_names.append(os.path.join(subfolder,outimage))                
            else:
                os.rename('coadd.fits',outimage)
                image_names.append(outimage)

            
            #os.rename('coadd.fits','unwise/'+outimage)

            #########################################
            ## COMBINE THE STD FRAMES USING SUM
            #########################################
            # just combine using average and then multiply by sqrt 2
            matchstring = "*w{}-img-m.std.fits".format(b)
            
            if subfolder is not None:
                allfiles = glob.glob(subfolder+'/'+galid+matchstring)
            else:
                allfiles = glob.glob(galid+matchstring)
            #print('allfiles with std images = ',allfiles)
            all_images = " ".join(allfiles)
            
            output_image = str(galid)+'-'
            s = 'swarp '+all_images+' -COMBINE_TYPE MEDIAN -WEIGHT_TYPE NONE -SUBTRACT_BACK N'
            #print(s)
            os.system(s)

            # now take sqrt of image values
            #im,h = fits.getdata('coadd.fits',header=True)
            #im = np.sqrt(2)*im
            # divide by number of images because we took average of data
            #im = im/len(all_images)
            
            weightname = str(galid)+'-unwise-w'+str(b)+'-coadd.std.fits'  
            #fits.writeto(weightname,im,header=h,overwrite=True)


            # just trying average!
            if subfolder is not None:
                os.rename('coadd.fits',os.path.join(subfolder,weightname))
                weight_names.append(os.path.join(subfolder,weightname))
            else:
                os.rename('coadd.fits',weightname)
                weight_names.append(weightname)            
        

    return image_names,weight_names,multiframe

def getgalex(ra,dec,imsize):
    """

    get galex image of a galaxy
    
    Input:
    * ra in deg
    * dec in deg
    * imsize in arcsec
    
    Returns:
    * image
    """
    

    # following procedure outlined here:
    # https://astroquery.readthedocs.io/en/latest/mast/mast.html

    # get data products in region near ra,dec
    obs_table = Observations.query_region("%12.8f %12.8f"%(ra,dec),radius=.1*u.arcmin)
    # create a flag to select galex data
    galexFlag = obs_table['obs_collection'] == 'GALEX'

    # separate out galex data
    data_products = Observations.get_product_list(obs_table[galexFlag])

    # download the observations
    manifest = Observations.download_products(data_products,productType="SCIENCE")

    for m in manifest:
        # choose the first NUV image
        if m['Local Path'].find('nd-int') > -1:
            nuv_path = m['Local Path']
            break
        
    # should be able to construct path from the obs_id in data_products
    # this will let us check if the image is already downloaded
    #
    # DONE: I can also save the cutout in a GALEX folder, and
    # look for the image before calling this function
    #
    # but I should also look for the image, because if I change the image
    # size of the cutout, I don't need to download the big FOV
    
    nuv,nuv_header = fits.getdata(nuv_path,header=True)

    # this is a big image, so we need to get a cutout

    nuv_wcs = WCS(nuv_header)
    position = SkyCoord(ra,dec,unit="deg",frame='icrs')
    try:
        cutout = Cutout2D(nuv,position,(imsize*u.arcsec,imsize*u.arcsec),wcs=nuv_wcs)
    except:
        print('WARNING: problem getting galex cutout')
        cutout = None

    return cutout
    
def display_image(image,percent=99.9,lowrange=False,mask=None,sigclip=True):
    from matplotlib import pyplot as plt
    lowrange=False
    if sigclip:
        clipped_data = sigma_clip(image,sigma_lower=5,sigma_upper=5)#,grow=10)
    else:
        clipped_data = image
    if lowrange:
        norm = simple_norm(clipped_data, stretch='linear',percent=percent)
    else:
        norm = simple_norm(clipped_data, stretch='asinh',percent=percent)

    plt.imshow(image, norm=norm,cmap='gray_r',origin='lower')

