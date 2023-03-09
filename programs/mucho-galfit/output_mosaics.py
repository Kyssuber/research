#!/usr/bin/env python
'''
GOAL:
- generate galfit output mosaic .png
'''

from matplotlib import pyplot as plt
import numpy as np
import warnings
warnings.filterwarnings('ignore')
import os
homedir = os.getenv("HOME")
import glob
from astropy.table import Table
from astropy.io import fits
from astropy.visualization import simple_norm
from astropy.wcs import WCS
from astropy.io import ascii
from scipy.stats import scoreatpercentile

#os.sys.path.append('/mnt/astrophysics/kconger_wisesize/github/gal_output')
#os.chdir('/mnt/astrophysics/kconger_wisesize/github/gal_output')

#have to call sample
cat = Table.read(homedir+'/sgacut_coadd.fits')
dummycat = Table.read(homedir+'/dummycat.fits')


class mosaic():
   
    def __init__(galname,band='w3',vfid,ncomp=1):
        self.galname = galname
        self.band = band
        #self.vfid = galname[0:8]
        self.vfid = vfid
        if self.vfid in dummycat['central galaxy']:
            self.ncomp = len(np.where(dummycat['central galaxy'] == vfid)[0]) + 1
        else:
            self.ncomp = ncomp
    
    def createmosaic(self,percentile1=.5,percentile2=99.5,p1residual=5,p2residual=99,cmap='viridis'):
        '''
        ARGS:
        percentile1 = min percentile for stretch of image and model
        percentile2 = max percentile for stretch of image and model
        p1residual = min percentile for stretch of residual
        p2residual = max percentile for stretch of residual
        cmap = colormap, default is viridis
        '''
        
        self.filename = glob.glob(self.galname+'-'+str(self.band)+'*Comp-galfit-out.fits')[0]
        pngname = self.galname+'-'+str(self.band)+'-'+str(self.ncomp)+'Comp-galfit-out.png'
        
        if self.convflag:
            self.filename = glob.glob(self.galname+'-'+str(self.band)+'*Comp-galfit-out-conv.fits')[0]
            pngname = self.galname+'-'+str(self.band)+'-'+str(self.ncomp)+'Comp-galfit-out-conv.png'
        
        image,h = fits.getdata(self.filename,1,header=True)
        model = fits.getdata(self.filename,2)
        residual = fits.getdata(self.filename,3)

        wcs = WCS(h)
        images = [image,model,residual,residual]
        titles = ['image','model','residual (img stretch)','residual (res stretch)']
        #titles = ['image','model','residual','residual (res stretch)']
        v1 = [scoreatpercentile(image,percentile1),
            scoreatpercentile(image,percentile1),
            scoreatpercentile(image,percentile1),
            scoreatpercentile(residual,p1residual)]
        v2 = [scoreatpercentile(image,percentile2),
            scoreatpercentile(image,percentile2),
            scoreatpercentile(image,percentile2),
            scoreatpercentile(residual,p2residual)]
        norms = [simple_norm(image,'asinh',max_percent=percentile2),
               simple_norm(image,'asinh',max_percent=percentile2),
               simple_norm(image,'asinh',max_percent=percentile2),
               simple_norm(residual,'linear',max_percent=p2residual)]
               
        plt.figure(figsize=(14,6))
        plt.subplots_adjust(wspace=.0)
        for i,im in enumerate(images): 
            ax = plt.subplot(1,4,i+1,projection=wcs)
            #ax = plt.subplot(1,4,i+1)

            plt.imshow(im,origin='lower',cmap=cmap,vmin=v1[i],vmax=v2[i],norm=norms[i])
            #if i == 1:
            #   plt.plot(xmax/8,ymax/8,'o',markersize=4,color='red')
            ax.set_xlabel('RA')
            if i == 0:
                ax.set_ylabel('DEC')
            else:
                plt.ylabel(' ')
                ax = plt.gca()
                ax.set_yticks([])
            plt.title(titles[i],fontsize=16)
        plt.savefig(pngname,dpi=300)
        plt.close()
        

if __name__ == '__main__':
    for i in range(0,len(cat)):
        m = mosaic(galname=cat['prefix'][i],vfid=cat['VFID'][i])
        m.createmosaic()
