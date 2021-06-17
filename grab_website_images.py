import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
import wget
import os
from matplotlib.image import mpimg
from astropy.wcs import WCS
from scipy.stats import scoreatpercentile
from astropy.visualization import simple_norm



def pull_LS_image(sample=vf,i,initial_pixscale):
    ra = str(sample['RA'][i])
    dec = str(sample['DEC'][i])
    galname = str(sample['VFID'][i])
    initial_pixscale = str(initial_pixscale)
    image_url = 'https://www.legacysurvey.org/viewer/cutout.fits?ra='+ra+'&dec='+dec+'&layer=ls-dr9&pixscale='+initial_pixscale
    filename = '/mnt/astrophysics/kconger_wisesize/LS_cutouts/'+galname+'.fits'
    if os.path.exists(filename):
        os.remove(filename)
        image = wget.download(image_url,out=filename)
    else:
        image = wget.download(image_url,out=filename)



def align_fits(sample=vf,i,initial_pixscale,percentile1=.5,percentile2=99.5,cmap='viridis'):
    image_path = '/mnt/astrophysics/kconger_wisesize/cutouts/unwise-'+galname+'-w3-img-m-trim.fits'
    image_WISE, h_WISE = fits.getdata(image_path,header=True)

    wcs_WISE = WCS(h_WISE)

    v1 = [scoreatpercentile(image_WISE,percentile1),scoreatpercentile(image_WISE,percentile1)]
    v2 = [scoreatpercentile(image_WISE,percentile2),scoreatpercentile(image_WISE,percentile2)]

    pull_LS_image(sample,i,initial_pixscale)
    image_path = '/mnt/astrophysics/kconger_wisesize/LS_cutouts/'+galname+'.fits'
    image_LS, h_LS = fits.getdata(image_path,header=True)


    for scale in np.arange(0.05,10,0.01):
        ratio = (h_WISE['NAXIS1']*h_WISE['CD1_1'])/(h_LS['NAXIS1']*h_LS['CD1_1'])
        if (ratio>0.97) & (ratio<1.03):
            print('complete: pixscale',scale)
            break
        else:
            pull_LS_image(sample,i,scale)
            image_LS, h_LS = fits.getdata(image_path,header=True)

    wcs_LS = WCS(h_LS)

    image_url = 'https://legacysurvey.org/viewer/cutout.jpg?ra='+ra+'&dec='+dec+'&layer=ls-ds9&pixscale='+str(scale)
    ls_jpg = wget.download(image_url)
    image_desired = galname+'.jpg'
    os.rename(ls_jpg,image_desired)

    os.chdir('/mnt/astrophysics/kconger_wisesize/LS_cutouts')
    images = [image_WISE,image_LS,ls_jpg]
    titles = ['WISE Image','LS fits','LS Image']

    plt.figure()
    plt.subplots_adjust(wspace=.0)
    for i,im in enumerate(images):
        plt.xlabe('RA')
        if i == 0:
            plt.subplot(1,3,i+1,projection = image_wcs)
            plt.imshow(im,origin='lower',cmap='viridis',vmin=v1[i],vmax[i])
            plt.ylabel("DEC")
            ax=plt.gca()
            ax.set_yticks([])
            if i == 1:
                plt.subplot(1,3,i+1,projection=wcs_LS[1])
            if i == 2:
                os.chdir('/mnt/astrophysics/kconger_wisesize/LS_cutouts')
                plt.imshow(mpimg.mread(ls_jpg),origin='lower')
                plt.gca().invert_yaxis()
                plt.axis("off")

            plt.title(titles[i],fontsize=16)


    filename = '/mnt/astrophysics/kconger_wisesize/LS_mosaics/'+galname+'-mosaic.png'
    if os.path.exists(filename):
        os.remove(filename)
        plt.savefig(galname+'-mosaic.png')
    else:
        plt.savefig(galname+'-mosaic.png')
    plt.close()




if __name__ == '__main__':
    homedir = os.getenv("HOME")
    vf = Table.read(homedir+'vfcut.fits')
    cutoutdir = homedir+'/github/cutouts/'
    maskdir = homedir+'/github/masks/'
