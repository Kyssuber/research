#!/usr/bin/env python
'''
GOAL:
- generate galfit output mosaic .png
'''

import sys
from matplotlib import pyplot as plt
import numpy as np
import os
import glob
from astropy.table import Table
from astropy.io import fits
from astropy.visualization import simple_norm
from astropy.wcs import WCS
from astropy.io import ascii
from scipy.stats import scoreatpercentile

class mosaic():
   
    def __init__(galname,vfid,band=None,ncomp=1):
        self.galname = galname
        self.band = band
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
        
        filename = glob.glob(self.galname+str(self.band)+'*Comp-galfit-out.fits')[0]
        pngname = self.galname+'-'+str(self.band)+'-'+str(self.ncomp)+'Comp-galfit-out.png'
        
        if self.convflag:
            filename = glob.glob(self.galname+str(self.band)+'*Comp-galfit-out-conv.fits')[0]
            pngname = self.galname+'-'+str(self.band)+'-'+str(self.ncomp)+'Comp-galfit-out-conv.png'
        
        image,h = fits.getdata(filename,1,header=True)
        model = fits.getdata(filename,2)
        residual = fits.getdata(filename,3)

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
    
    if '-h' in sys.argv or '--help' in sys.argv:
        print("Usage: %s [-param_file (name of parameter file, no single or double quotation marks)]" % sys.argv[0])
        sys.exit(1)
    
    if '-param_file' in sys.argv:
        p = sys.argv.index('-param_file')
        param_file = str(sys.argv[p+1])
  
       
    #create dictionary with keywords and values, from parameter .txt file
    
    homedir = os.getenv("HOME")
    
    param_dict = {}
    with open(homedir+'/'+param_file) as f:
        for line in f:
            try:
                key = line.split()[0]
                val = line.split()[1]
                param_dict[key] = val
            except:
                continue
    
    ############
    
    cat_path = param_dict['vf_sample']
    cat = Table.read(cat_path)
    
    dummycat_path = param_dict['dummycat']
    dummycat = Table.read(dummycat_path)
    
    directory = param_dict['directory']
    gal_output_path = param_dict['gal_output_path']
    
    band = param_dict['band']    
    
    #cd to location of all galfit output .fits files
    os.chdir(gal_output_path) 
    
    for i in range(0,len(cat)):
        m = mosaic(galname=cat['prefix'][i],vfid=cat['VFID'][i],band=band)
        try:
            m.createmosaic()
        except:
            print('did not complete for ',m.galname)