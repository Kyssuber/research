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
        generic = '/mnt/astrophysics/wisesize/'+galname+'/unwise-'+galname+'-w3-'
        filepath_im = generic + 'img-m.fits'
        filepath_n = generic + 'n-m.fits'
        filepath_std = generic + 'std-m.fits'
        filepath_mask = generic + 'img-m-mask.fits'
        filepath_inv = generic + 'img-m-inv-mask.fits'

        #open fits, define WCS coordinates for galaxy
        hdu1 = fits.open(filepath_im)[0]
        wcs1 = WCS(hdu1.header)
        hdu2 = fits.open(filepath_n)[0]
        wcs2 = WCS(hdu2.header)
        hdu3 = fits.open(filepath_std)[0]
        wcs3 = WCS(hdu3.header)
        hdu4 = fits.open(filepath_mask)[0]
        wcs4 = WCS(hdu4.header)
        hdu5 = fits.open(filepath_inv)[0]
        wcs5 = WCS(hdu5.header)

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
        cutout4 = Cutout2D(hdu4.data, center, size, wcs=wcs4)
        cutout5 = Cutout2D(hdu5.data, center, size, wcs=wcs5)

        #updating header information and saving the result
        hdu1.data = cutout1.data
        hdu2.data = cutout2.data
        hdu3.data = cutout3.data
        hdu4.data = cutout4.data
        hdu5.data = cutout5.data
        hdu1.header.update(cutout1.wcs.to_header())
        hdu2.header.update(cutout2.wcs.to_header())
        hdu3.header.update(cutout3.wcs.to_header())
        hdu4.header.update(cutout4.wcs.to_header())
        hdu5.header.update(cutout5.wcs.to_header())
        cutout_filename1 = 'unwise-'+galname+'-w3-img-m-trim.fits'
        cutout_filename2 = 'unwise-'+galname+'-w3-n-m-trim.fits'
        cutout_filename3 = 'unwise-'+galname+'-w3-std-m-trim.fits'
        cutout_filename4 = 'unwise-'+galname+'-w3-img-m-mask-trim.fits'
        cutout_filename5 = 'unwise-'+galname+'-w3-img-m-inv-mask-trim.fits'

        hdu1.writeto('/mnt/astrophysics/wisesize/'+galname+'/'+cutout_filename1)
        hdu2.writeto('/mnt/astrophysics/wisesize/'+galname+'/'+cutout_filename2)
        hdu3.writeto('/mnt/astrophysics/wisesize/'+galname+'/'+cutout_filename3)
        hdu4.writeto('/mnt/astrophysics/wisesize/'+galname+'/'+cutout_filename4)
        hdu5.writeto('/mnt/astrophysics/wisesize/'+galname+'/'+cutout_filename5)
