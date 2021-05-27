import numpy as np
from astropy.nddata import Cutout2D
from matplotlib import pyplot as plt
from astropy import units as u
from astropy.io import fits
import os
from astropy.table import Table
from astropy.wcs import WCS

homedir = os.getenv("HOME")
#read vf table (in home directory)
vfmain = Table.read(homedir+'/vf_north_v1_main.fits')

def trim(gal_sample):
    for i in range(0,len(gal_sample)):
        galname = str(gal_sample['VFID'][i])
        #define filepaths
        filepath_im = '/mnt/astrophysics/wisesize/'+galname+'/unwise-'+galname+'-w3-img-m.fits'
        filepath_n = '/mnt/astrophysics/wisesize/'+galname+'/unwise-'+galname+'-w3-n-m.fits'
        filepath_std = '/mnt/astrophysics/wisesize/'+galname+'/unwise-'+galname+'-w3-std-m.fits'

        #open fits, define WCS coordinates for galaxy
        hdu1 = fits.open(filepath_im)[0]
        wcs1 = WCS(hdu1.header)
        hdu2 = fits.open(filepath_n)[0]
        wcs2 = WCS(hdu2.header)
        hdu3 = fits.open(filepath_std)[0]
        wcs3 = WCS(hdu3.header)

        #extract galaxy's radius (in arcsec) from vf table entry
        vf = gal_sample[i]
        rad = vf['radius']

        #define center of image to be center of galaxy (img is 500x500 px)
        center = (250,250)
        #desired image size, about 2.5x the galaxy size (will be identical for each cutout)
        size = u.Quantity((rad*2.5, rad*2.5),u.arcsec)

        #generating the cutouts
        cutout1 = Cutout2D(hdu1.data, center, size, wcs=wcs1)
        cutout2 = Cutout2D(hdu2.data, center, size, wcs=wcs2)
        cutout3 = Cutout2D(hdu3.data, center, size, wcs=wcs3)

        #updating header information and saving the result
        hdu1.data = cutout1.data
        hdu2.data = cutout2.data
        hdu3.data = cutout3.data
        hdu1.header.update(cutout1.wcs.to_header())
        hdu2.header.update(cutout2.wcs.to_header())
        hdu3.header.update(cutout3.wcs.to_header())
        cutout_filename1 = 'unwise-'+galname+'-w3-img-m-trim.fits'
        cutout_filename2 = 'unwise-'+galname+'-w3-n-m-trim.fits'
        cutout_filename3 = 'unwise-'+galname+'-w3-std-m-trim.fits'

        hdu1.writeto('/mnt/astrophysics/wisesize/'+galname+'/'+cutout_filename1)
        hdu2.writeto('/mnt/astrophysics/wisesize/'+galname+'/'+cutout_filename2)
        hdu3.writeto('/mnt/astrophysics/wisesize/'+galname+'/'+cutout_filename3)
