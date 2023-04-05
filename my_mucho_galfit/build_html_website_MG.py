'''
GOAL:
- generate .html files that comprise a website interface for viewing vf galaxy images and models
OUTPUT:
- 
General instructions for opening .html locally:
COMMAND+O (that is a capital o, not zero!)
- 
'''
import numpy as np
from matplotlib import pyplot as plt
import os
from astropy.table import Table
from astropy.io import ascii
import glob
from astropy.io import fits
import wget
import matplotlib.image as mpimg
from astropy.wcs import WCS
from scipy.stats import scoreatpercentile
from astropy.visualization import simple_norm
from reproject import reproject_interp
import sys

'''
#NOTES:
#need some way to filter out the group galaxies, or at least disable the hyperlink on homepage that would otherwise lead to galfit params. also dummycat...
#when running .html locally, I'll want to place every folder and .html file into a specific directory. *oRgAnIzAtIoN*
'''

class HomePage():
    def __init__(self,website_title=None,home_color=None,homepage_name=None,catalog=None,dummycat=None,local_path=None,
                 path_to_galhtml=None,path_to_params=None,LS_cutout_folder=None,LS_mosaic_folder=None,mask_folder=None,fits_folder=None,
                 gal_mosaic_folder=None,indices=None):
        
        #index = 0 for w3_nopsf, 1 for w3_psf, 2 for r_nopsf, 3 for r_psf
        if indices != None:
            self.indices = [0,1,2,3]
        else:
            self.indices = indices
        
        self.cat = Table.read(catalog)   #vf subsample catalog
        self.dummycat = Table.read(dummycat)   #catalog flagging galaxies with 2+ Sersic objects in their cutouts
        self.path_to_galhtml = path_to_galhtml   #path to the folder containing each individual galpage html file
        self.local_path = local_path   #path to the folder containing html directories/files when running the website locally
        self.homepage_name = homepage_name
        self.mask_mosaics = mask_folder   #path to the folder holding the .png mask mosaics
        self.LS_mosaics = LS_mosaic_folder   #contains mosaics of w3 cutout, r-band cutout, and LS cutout
        self.htmlpath = self.path_to_galhtml+self.homepage_name   #path to homepage html
        self.fits_folder = fits_folder   #path to all w3, r-band postage stamps
        self.gal_mosaic_folder = gal_mosaic_folder   #.png mosaics of galaxy, model, and residual
        
        self.website_title = website_title
        self.home_color = home_color
        
        self.path_to_params = path_to_params   #path to galfit output parameter files
        self.LS_cutouts = LS_cutout_folder   #will contain the .png files of LS cutouts       

        #if I am using the v2_snrcoadd.fits file, the length is 6780
        if len(self.cat)>702:
            self.cutcat = self.cat[self.cat['sgacut_flag']]
        else:
            self.cutcat = self.cat
         
        #call the remaining parameter files...
        self.w3params_nopsf = Table.read(self.path_to_params+'/output_params_W3_nopsf.fits') 
        self.params_w3_nopsf = self.w3params_nopsf.copy()  #will help to identify groups below...somewhere. This script is becoming unruly.
        self.w3params_psf = Table.read(self.path_to_params+'/output_params_W3_psf.fits')
        self.rparams_nopsf = Table.read(self.path_to_params+'/output_params_r_nopsf.fits')
        self.rparams_psf = Table.read(self.path_to_params+'/output_params_r_psf.fits')
    
    def html_setup(self):
        
        with open(self.htmlpath, 'w') as html:
            
            #initialize title page
            html.write('<html><body>\n')
            html.write('<title>Virgowise Project (masking)</title>\n')
            html.write(f'<body style="background-color:{self.home_color};">\n')   #be gentle with the background color (e.g., no limegreen)
            html.write('<style type="text/css">\n')
            html.write('table, td, th {padding: 5px; text-align: center; border: 2px solid black;}\n')
            html.write('p {display: inline-block;;}\n')
            html.write('</style>\n')
            html.write(f'<font size="40">{self.website_title}</font>\n')
            
            #begin and populate the table. first line is header information; loop creates the data rows.
            html.write('<table><tr><th>Index</th><th>LS Cutout</th><th>Prefix</th><th>RA</th><th>DEC</th><th>Comments</th>\n')
            
            for i in range(len(self.cutcat)):
                html.write('<tr><td>'+str(i)+'</td>\n')   #add index number
                html.write('<td><img src = "' + self.LS_cutouts + self.cutcat['VFID'][i] + '-LS.jpg' + '" height="25%" width = "25%"></img></td>\n')   #cutouts will have the name VFIDxxxx-LS.png, using the v2 IDs
                
                #if galaxy is part of a group, then all parameters will be zeros, and heaps of trouble in terms of generating the galaxy page arises. As such, we effectively disable the hyperlinks in these cases.
                if self.params_w3_nopsf['xc'][i]>0:
                    
                    print('Creating htmlpage for '+self.cutcat['VFID'][i])
                    
                    #CREATE SINGLE GALPAGE using the GalPage class (see below)
                    single_galpage = GalPage(galaxy_index=i, psf_indices=self.indices, 
                                             page_name=self.cutcat['VFID'][i]+'.html', catalog=self.cutcat, 
                                             dummycat=self.dummycat, local_path=self.local_path, 
                                             path_to_galhtml=self.path_to_galhtml, LS_cutout_folder=self.LS_cutouts, 
                                             LS_mosaic_folder=self.LS_mosaics, mask_folder=self.mask_mosaics, 
                                             fits_folder=self.fits_folder, gal_mosaic_folder=self.gal_mosaic_folder, 
                                             w3params_nopsf=self.w3params_nopsf, w3params_psf=self.w3params_psf, 
                                             rparams_nopsf=self.rparams_nopsf, rparams_psf=self.rparams_psf, 
                                             homepage_name=self.homepage_name)
                    
                    pagename = single_galpage.page_name

                    html.write('<td><a href='+self.path_to_galhtml+pagename+'>'+str(self.cutcat['prefix'][i])+'</a></td>\n')   #text hyperlink to galaxy page VFIDxxxx.html (pagename)
                
                else:
                    print('Unable to create htmlpage for '+self.cutcat['VFID'][i])
                    html.write('<td>'+str(self.cutcat['prefix'][i])+'</a></td>\n')
                
                html.write('<td>'+str(self.cutcat['RA_1'][i])+'</td>\n')
                html.write('<td>'+str(self.cutcat['DEC_1'][i])+'</td>\n')

                #if the VFID (v2) is in the dummycat central galaxy column, then write in the Comments column that this particular galaxy is a member of a "Moustakas group"
                if self.cutcat['VFID'][i] in self.dummycat['central galaxy']:
                    html.write('<td>Group Galaxy</td>\n')
                else:
                    html.write('<td>--</td>\n')

            html.write('</tr></table>\n')
            html.write('<br /><br />\n')
            html.write('</html></body>\n')
            html.close()
    
    #creating a few separate functions in order to independently run the PNG creation scripts - if I run them as part of the homepage class loop I use to generate every galpage.html file, then an application memory problem arises. (Why not remove the methods from the class altogether if I won't use them in the initial loop as intended? I think they are organized more nicely as part of the galpage class; and since the variable names are already entangled, I may as well not tinker any further.)           
    
    def create_LS_figures(self):
        
        for i in range(len(self.cutcat)):
            if self.params_w3_nopsf['xc'][i]>0:   #I add this every time, since I use this condition to create the galhtml pages

                #I set test=True to avoid running the automatic execution of the function that creates galhtml pages
                single_galaxy = GalPage(galaxy_index=i, psf_indices=self.indices, 
                                         page_name=self.cutcat['VFID'][i]+'.html', catalog=self.cutcat, 
                                         dummycat=self.dummycat, local_path=self.local_path, 
                                         path_to_galhtml=self.path_to_galhtml, LS_cutout_folder=self.LS_cutouts, 
                                         LS_mosaic_folder=self.LS_mosaics, fits_folder=self.fits_folder, test=True)
                
                print('Creating LS cutout for '+single_galaxy.VFID)
                single_galaxy.compile_LS_cutouts()
                print('Creating LS mosaic for '+single_galaxy.VFID)
                single_galaxy.create_LS_mosaics()
                
                del single_galaxy
                
    def create_galfit_mosaics(self, psf_index):
        
        index_dict = {0:'W3, no conv', 1:'W3, conv', 2:'r-band, no conv', 3:'r-band, conv'}
        
        for i in range(len(self.cutcat)):
            if self.params_w3_nopsf['xc'][i]>0:   #I add this every time, since I use this condition to create the galhtml pages
                
                #I set test=True to avoid running the automatic execution of the function that creates galhtml pages
                single_galaxy = GalPage(galaxy_index=i, psf_indices=self.indices, page_name=self.cutcat['VFID'][i]+'.html', 
                                        catalog=self.cutcat, dummycat=self.dummycat, local_path=self.local_path, 
                                        path_to_galhtml=self.path_to_galhtml, fits_folder=self.fits_folder, 
                                        gal_mosaic_folder=self.gal_mosaic_folder, w3params_nopsf=self.w3params_nopsf, 
                                        w3params_psf=self.w3params_psf, rparams_nopsf=self.rparams_nopsf, 
                                        rparams_psf=self.rparams_psf, test=True)                
        
                single_galaxy.create_model_mosaics_names()
                print('Creating GALFIT mosaics for '+single_galaxy.VFID+f' {index_dict[psf_index]}')
                single_galaxy.create_model_mosaics(psf_index = psf_index)
                
                del single_galaxy
                
    def create_mask_mosaics(self):
        
        for i in range(len(self.cutcat)):
            if self.params_w3_nopsf['xc'][i]>0:   #I add this every time, since I use this condition to create the galhtml pages

                #I set test=True to avoid running the automatic execution of the function that creates galhtml pages
                single_galaxy = GalPage(galaxy_index=i, psf_indices=self.indices, page_name=self.cutcat['VFID'][i]+'.html', 
                                        catalog=self.cutcat, dummycat=self.dummycat, local_path=self.local_path, 
                                        path_to_galhtml=self.path_to_galhtml, fits_folder=self.fits_folder, 
                                        gal_mosaic_folder=self.gal_mosaic_folder, mask_folder=self.mask_mosaics, test=True)                
                
                print('Creating mask mosaic for '+single_galaxy.VFID)
                single_galaxy.create_mask_mosaics()
                
                del single_galaxy
                
class GalPage():
    def __init__(self,galaxy_index=None, psf_indices = [0,1,2,3], page_name=None, catalog=None, dummycat=None, 
                 local_path=None, path_to_galhtml=None, LS_cutout_folder=None, LS_mosaic_folder=None, mask_folder=None, 
                 fits_folder=None, gal_mosaic_folder=None, w3params_nopsf=None, w3params_psf=None, rparams_nopsf=None, 
                 rparams_psf=None, homepage_name=None, test=False):
        
        self.galaxy_index = int(galaxy_index)
        self.psf_indices = psf_indices
        self.catalog = catalog
        self.cutcat = catalog
        self.dummycat = dummycat
        
        self.homepage_name = homepage_name
        self.local_path = local_path
        self.LS_cutouts = LS_cutout_folder
        self.LS_mosaics = LS_mosaic_folder
        self.mask_mosaics = mask_folder
        self.fits_folder = fits_folder
        self.gal_mosaic_folder = gal_mosaic_folder
        self.path_to_galhtml = path_to_galhtml
        
        self.w3params_nopsf = w3params_nopsf
        self.w3params_psf = w3params_psf
        self.rparams_nopsf = rparams_nopsf
        self.rparams_psf = rparams_psf
        
        self.RA = self.cutcat['RA_1'][self.galaxy_index]
        self.RA = str(self.RA)
        self.DEC = self.cutcat['DEC_1'][self.galaxy_index]
        self.DEC = str(self.DEC)
        self.VFID = self.cutcat['VFID'][self.galaxy_index]
        self.VFID = str(self.VFID)
        
        self.page_name = self.VFID+'.html'   #name of galaxy html page
        
        if test==False:
            self.gal_htmlpath = self.path_to_galhtml+self.page_name

        #number of instances of central VFID galaxy appearing in the relevant dummycat column --> ncomp-1
        if self.VFID in self.dummycat['central galaxy']:
            ncomp_flag = (self.dummycat['central galaxy']==self.VFID)
            self.ncomp = len(self.dummycat['central galaxy'][ncomp_flag])+1   #the +1 accounts for the central galaxy itself
            self.ext_list = self.dummycat['ID'][ncomp_flag]
        else:
            self.ncomp = 1  
                
        self.objname = self.cutcat['objname'][self.galaxy_index]
        
        #need w3 header information and w3, r-band image data
        path_to_w3 = glob.glob(self.fits_folder+self.objname+'-custom-image-W3.fits')[0]
        path_to_r = glob.glob(self.fits_folder+self.objname+'-custom-image-r.fits')[0]
        self.wise_im, self.wise_header = fits.getdata(path_to_w3, header=True)
        self.r_im, self.r_header = fits.getdata(path_to_r, header=True)
        
        #WCS header information
        self.wcs_w3 = WCS(self.wise_header)
        self.wcs_r = WCS(self.r_header)
        
        #defining the mask paths
        self.w3mask_path = glob.glob(self.fits_folder+self.objname+'-custom-image-wise-mask.fits')[0]
        self.rmask_path = glob.glob(self.fits_folder+self.objname+'-custom-image-r-mask.fits')[0]
        
        #if not testing the various functions on one galaxy (test==True), then run only the functions that are required for (1) variables and (2) the actual galaxy html pages. Think of these functions as the 'defaults' that enable the user to then generate the PNG files at their liberty.
        if test==False:
            self.create_model_mosaics_names()
            self.tabulate_parameters()
            self.WRITETHEGALPAGE()
        
    #download scaled RGB jpg of VFIDxxxx galaxy; place in self.filename_LS
    def compile_LS_cutouts(self):
        
        pixscale=1   #standard
        pscale=np.abs(float(self.wise_header['CD1_1']))   #grab transformation matrix of WISE image
        xsize=np.abs(int(self.wise_header['NAXIS1']))   #grab length of WISE image
        xsize_arcsec=pscale*3600*xsize   #length convert to arcseconds
        imsize=int(xsize_arcsec/pixscale)   #convert length to an integer
        imsize=str(imsize)   #convert integer length to a...string
        
        self.filename_LS = self.LS_cutouts+self.VFID+'-LS.jpg'
        
        image_url = f'https://www.legacysurvey.org/viewer/cutout.jpg?ra={self.RA}&dec={self.DEC}&layer=ls-dr9&size={imsize}&pixscale={1}'
        if os.path.exists(self.filename_LS):
            os.remove(self.filename_LS)
            image = wget.download(image_url,out=self.filename_LS)
        else:
            image = wget.download(image_url,out=self.filename_LS)
    
    def create_LS_mosaics(self):
        
        #first must be sure that r-band is projected onto w3 coordinates (ensures all images have same size)
        r_scaled, footprint = reproject_interp((self.r_im, self.r_header), self.wise_header)
        titles = ['W3 Image', 'r-band Image', 'LS Image']
        images = [self.wise_im, r_scaled, self.filename_LS]
        
        plt.figure(figsize=(12,6))
        for i,im in enumerate(images):
            plt.xlabel('RA')
            if (i==0)|(i==1):
                if i==0:
                    plt.subplot(1,len(images),i+1,projection = self.wcs_w3)
                if i==1:
                    plt.subplot(1,len(images),i+1,projection = self.wcs_r)
                try:
                    norm = simple_norm(images[i],stretch='asinh',max_percent=99.5)
                    plt.imshow(images[i],origin='lower',cmap='viridis',norm=norm)
                except:
                    norm = simple_norm(self.r_im,stretch='asinh',max_percent=99.5)
                    plt.imshow(self.r_im,origin='lower',cmap='viridis',norm=norm)
                plt.ylabel('DEC')
                ax = plt.gca()
                ax.set_yticks([])
            
            if i==2:
                plt.subplot(1,len(images),i+1)
                plt.imshow(mpimg.imread(images[i]),origin='lower')
                plt.gca().invert_yaxis()
                plt.axis("off")

            plt.subplots_adjust(wspace=0,hspace=0)
            plt.title(titles[i],fontsize=24)
        
        if os.path.exists(self.LS_mosaics+self.VFID+'-LS-mosaic.png'):
            os.remove(self.LS_mosaics+self.VFID+'-LS-mosaic.png')
            plt.savefig(self.LS_mosaics+self.VFID+'-LS-mosaic.png',bbox_inches = 'tight',pad_inches=0.2)
        else:
            plt.savefig(self.LS_mosaics+self.VFID+'-LS-mosaic.png',bbox_inches='tight',pad_inches=0.2)
        plt.close()
    
    def create_model_mosaics_names(self):

        self.file_w3_nopsf = self.fits_folder+self.objname+'-W3-out1.fits'
        self.file_w3_psf = self.fits_folder+self.objname+'-W3-out2.fits'
        self.file_r_nopsf = self.fits_folder+self.objname+'-r-out1.fits'
        self.file_r_psf = self.fits_folder+self.objname+'-r-out2.fits'
        
        self.mosaic_names = [self.file_w3_nopsf,self.file_w3_psf,self.file_r_nopsf,self.file_r_psf]
        
        self.models = [fits.getdata(self.file_w3_nopsf,2),
                 fits.getdata(self.file_w3_psf,2),
                 fits.getdata(self.file_r_nopsf,2),
                 fits.getdata(self.file_r_psf,2)]
        
        self.residuals = [fits.getdata(self.file_w3_nopsf,3),
                    fits.getdata(self.file_w3_psf,3),
                    fits.getdata(self.file_r_nopsf,3),
                    fits.getdata(self.file_r_psf,3)]
        
        self.pngnames = [self.gal_mosaic_folder+self.VFID+'-'+'galfit-model-w3-nopsf.png',
                   self.gal_mosaic_folder+self.VFID+'-'+'galfit-model-w3-psf.png',
                   self.gal_mosaic_folder+self.VFID+'-'+'galfit-model-r-nopsf.png',
                   self.gal_mosaic_folder+self.VFID+'-'+'galfit-model-r-psf.png']
        
        #print('For self.create_model_mosaics(index), index=0 is w3_nopsf, 1 is w3_psf, 2 is r_nopsf, 3 is r_psf.')

    def create_model_mosaics(self, psf_index, percentile1=.5, percentile2=99.5, p1residual=5, p2residual=99, cmap='viridis'):
        
        '''
        ARGS:
        index = 0 for w3_nopsf, 1 for w3_psf, 2 for r_nopsf, 3 for r_psf
        percentile1 = min percentile for stretch of image and model
        percentile2 = max percentile for stretch of image and model
        p1residual = min percentile for stretch of residual
        p2residual = max percentile for stretch of residual
        cmap = colormap, default is viridis
        ''' 
        
        #for index in self.psf_indices:
        if psf_index<2:   #w3 is index=0 or index=1
            images = [self.wise_im,self.models[psf_index],self.residuals[psf_index],self.residuals[psf_index]]
        if psf_index>=2:   #r-band is index=2 or index=3
            images = [self.r_im,self.models[psf_index],self.residuals[psf_index],self.residuals[psf_index]]
        titles = ['Image','Model','Residual (img stretch)','Residual (res stretch)']

        v1 = [scoreatpercentile(images[0],percentile1),
            scoreatpercentile(images[0],percentile1),
            scoreatpercentile(images[0],percentile1),
            scoreatpercentile(images[3],p1residual)]
        v2 = [scoreatpercentile(images[0],percentile2),
            scoreatpercentile(images[0],percentile2),
            scoreatpercentile(images[0],percentile2),
            scoreatpercentile(images[3],p2residual)]

        norms = [simple_norm(images[0],'asinh',max_percent=percentile2,min_cut=v1[0],max_cut=v2[0]),
               simple_norm(images[0],'asinh',max_percent=percentile2,min_cut=v1[1],max_cut=v2[1]),
               simple_norm(images[0],'asinh',max_percent=percentile2,min_cut=v1[2],max_cut=v2[2]),
               simple_norm(images[0],'linear',max_percent=p2residual,min_cut=v1[3],max_cut=v2[3])]

        plt.figure(figsize=(14,6))
        plt.subplots_adjust(wspace=.0)
        for i,im in enumerate(images): 
            ax = plt.subplot(1,4,i+1,projection=self.wcs_w3)
            plt.imshow(im,origin='lower',cmap=cmap,norm=norms[i])  #vmin=v1[i],vmax=v2[i]
            ax.set_xlabel('RA')
            if i == 0:
                ax.set_ylabel('DEC')
            else:
                plt.ylabel(' ')
                ax = plt.gca()
                ax.set_yticks([])
            plt.title(titles[i],fontsize=16)
        plt.savefig(self.pngnames[psf_index],bbox_inches='tight', pad_inches=0.2)   #dpi=200
        plt.close()    

    def create_mask_mosaics(self, percentile1=.5, percentile2=99.5, cmap='viridis'):
        self.w3_mask = fits.getdata(self.w3mask_path)
        self.r_mask = fits.getdata(self.rmask_path)

        titles = ['W3 Image', 'W3 Mask', 'r-band Image', 'r-band Mask']
        images = [self.wise_im, self.w3_mask, self.r_im, self.r_mask]
        
        #norms for images but not for masks
        v1 = [scoreatpercentile(images[0],percentile1),None,scoreatpercentile(images[2],percentile1),None]
        v2 = [scoreatpercentile(images[0],percentile2),None,scoreatpercentile(images[2],percentile2),None]
        norms = [simple_norm(images[0],'asinh',max_percent=percentile2,min_cut=v1[0],max_cut=v2[0]),None, 
                 simple_norm(images[2],'asinh',max_percent=percentile2,min_cut=v1[2],max_cut=v2[2]),None]
              
        plt.figure(figsize=(14,6))
        plt.subplots_adjust(wspace=.0)
        for i,im in enumerate(images): 
            if i<=1:
                ax = plt.subplot(1,4,i+1,projection=self.wcs_w3)
            if i>=2:
                ax = plt.subplot(1,4,i+1,projection=self.wcs_r)
            plt.imshow(im,origin='lower',cmap=cmap,norm=norms[i])
            ax.set_xlabel('RA')
            if i == 0:
                ax.set_ylabel('DEC')
            else:
                plt.ylabel(' ')
                ax = plt.gca()
                ax.set_yticks([])
            plt.title(titles[i],fontsize=16)
        plt.savefig(self.mask_mosaics+self.VFID+'-mask_mosaic.png',bbox_inches='tight', pad_inches=0.2)   #dpi=200
        plt.close() 
    
    def tabulate_parameters(self):
        
        params_list = [self.w3params_nopsf,self.w3params_psf,self.rparams_nopsf,self.rparams_psf]
        
        self.page_params = []
        if self.ncomp>1:
            self.page_extparams = []
        
        #create list of every parameter row for the given indices
        for index in self.psf_indices:
            params=params_list[index]   #selects the correct params table
            param_row=params[params['VFID']==self.VFID]
            self.page_params.append(param_row)
            if self.ncomp>1:
                try:
                    single_ncomp_params=[]  #create empty list for a single external galaxy
                    for num in range(self.ncomp):  #say n=2. then num = 0, 1 (indices of external galaxies in self.ext_list)
                        param_row=params[params['VFID']==self.ext_list[num]]  #find where VFID = external galaxy ID
                        single_ncomp_params.append(param_row)  #append row values to single ncomp list 
                    self.page_extparams.append(single_ncomp_params)  #will comprise list of external galaxy lists. LISTCEPTION.
                except:
                    print('ncomp>1 galaxies not yet implemented.')
    
    #create VFIDxxxx.html for the galaxy!
    def WRITETHEGALPAGE(self):
        
        title_dictionary = {0:'w3 noConv',
                            1:'w3 Conv',
                            2:'r-band noConv',
                            3:'r-band Conv'}
        
        LS_path = self.LS_mosaics+self.VFID+'-LS-mosaic.png'
        mask_path = self.mask_mosaics+self.VFID+'-mask_mosaic.png'

        with open(self.gal_htmlpath,'w') as html:
            
            i = int(self.galaxy_index)
            
            html.write('<html><body> \n')
            html.write(f'<title>{self.VFID}-{self.objname}</title> \n')
            html.write('<style type="text/css"> \n')
            html.write('.img-container{text-align: left;} \n')
            html.write('table, td, th {padding: 5px; text-align: center; border: 1px solid black;} \n')
            html.write('p {display: inline-block;;} \n')
            html.write('</style> \n')
        
            html.write(f'<font size="40">Central Galaxy: {self.VFID}-{self.objname} </font><br /> \n')
            html.write(f'<a href={self.homepage_name}>Return to Homepage</a></br /> \n')

            if i != len(self.cutcat)-1:
                html.write('<a href='+str(self.cutcat['VFID'][i+1])+'.html>Next Galaxy</a></br /> \n') 

            if i != 0:
                html.write('<a href='+str(self.cutcat['VFID'][i-1])+'.html>Previous Galaxy</a></br /> \n')
            
            #add the LS mosaic
            html.write('<div class='+'"'+'img-container'+'"> <!-- Block parent element --> <img src='+'"'+LS_path+'" height="50%" width="55%" /><br /> \n')
            
            for index in self.psf_indices:
                mosaic_path = self.pngnames[index]
                params = self.page_params[index]
                html.write(f'<font size="30"> GALFIT Output Mosaic {title_dictionary[index]}:</font><br /> \n')
                html.write('<div class='+'"'+'img-container'+'"> <!-- Block parent element --> <img src='+'"'+mosaic_path+'" height="70%" width="90%" /><br /> \n')
                
                html.write('<table><tr><th>VFID</th><th>Type</th><th>xc</th><th>xc_err</th><th>yc</th><th>yc_err</th><th>mag</th><th>mag_err</th><th>Re</th><th>Re_err</th><th>nser</th><th>nser_err</th><th>BA</th><th>BA_err</th><th>PA</th><th>PA_err</th><th>err_flag</th></tr> \n')

                html.write('<tr><td>'+self.VFID+'</td> \n')
                html.write('<td>Host</td> \n')
                for p in range(1,15):   #for every parameter value in a row
                    html.write(f'<td>{params[0][p]}</td> \n')   #0th row, which does not change as there is only one row here
                html.write(f'<td>{params[0][17]}</td></tr> \n')   #error flag

                #add the external galaxy parameters, if any.
                if self.ncomp>1:
                    try:
                        param_ext=self.page_extparams[index]  #isolate band&psf parameters of the external galaxies (if there are multiple galaxies, then there will be multiple lists for one index/central galaxy. 
                        for num in range(self.ncomp):  #say n=2. then num = 0, 1 (indices of external galaxies in self.ext_list)
                            single_param_ext = param_ext[num]  #isolate the self.ext_list[num] external galaxy parameters
                            html.write('<tr><td>'+self.ext_list[num]+'</td> \n')  #VFID of external galaxy
                            html.write('<td>External</td> \n')
                            for p in range(1,15):
                                html.write(f'<td>{single_param_ext[p]}</td> \n')
                            if num==self.ncomp-1:
                                html.write(f'<td>{params[17]}</td></tr> \n')
                    except:
                        print('ncomp>1 galaxies not yet implemented.')
                        
                #I want a table for every nopsf/psf pair, so the table will finish if index is psf w3 or psf r-band (OR neither, in the case that I don't include any psf bands...). I would attempt to account for the case where I am only including nopsf tables, but I do not anticipate this condition being met. Ever.
                if (index==1)|(index==3)|(index==np.max(self.psf_indices)):
                    html.write('</tr></table> \n')
                
                #add the mask mosaic
                if index == np.max(self.psf_indices):
                    html.write('<div class='+'"'+'img-container'+'"> <!-- Block parent element --> <img src='+'"'+mask_path+'" height="70%" width="85%" /><br /> \n')
                
            html.write(f'<a href={self.homepage_name}>Return to Homepage</a></br /> \n')
           
            if i != len(self.cutcat)-1:
                html.write('<a href='+str(self.cutcat['VFID'][i+1])+'.html>Next Galaxy</a></br /> \n') 
            if i != 0:
                html.write('<a href='+str(self.cutcat['VFID'][i-1])+'.html>Previous Galaxy</a></br /> \n')

            html.write('<br /><br />\n')    
            html.write('</html></body>\n')     

            html.close()    
            
if __name__ == '__main__':    
    
    print("""*use -h or --help for the list of arguments*
    
    USAGE (iPython/Jupyter):
    ---Running this program automatically initiates the HomePage class (hp)
    ---if the -test arg is False, then replace single_galpage with hp in order to generate all images and such:
            --hp.html_setup() --> create .html homepages for all galaxies in the VF subsample; also initiates GalPage 
            class(single_galpage) for every galaxy in a loop, creating all relevant folders and files.
            --hp.create_LS_figures()
            --hp.create_galfit_mosaics()
            --hp.create_mask_mosaics(psf_index) --> (w3, nopsf), 1 (w3, psf), 2 (r, nopsf), 3 (r, psf)
    ---If the -test arg is True, then the user is wanting to test the GalPage class (this is a stategic idea, 
            since this class contains the bulk of the functions required for the script to run successfully). 
            In this case, the script automatically defines a single galpage class (single_galpage; 
            uses index 0 for a random test galaxy), with which the user can test the following methods:
            --single_galpage.compile_LS_cutouts()
            --single_galpage.create_LS_mosaics()
            --single_galpage.create_model_mosaics_names()
            --single_galpage.create_model_mosaics(psf_index) --> 0 (w3, nopsf), 1 (w3, psf), 2 (r, nopsf), 3 (r, psf)
            --single_galpage.create_mask_mosaics()
            --single_galpage.tabulate_parameters()
            --single_galpage.WRITETHEGALPAGE()            
            
    """)
    print('-----------------------------------------------------')
    print()
    
    if '-h' in sys.argv or '--help' in sys.argv:
        print("Usage: %s [-param_file (name of parameter file, no single or double quotation marks)] \n [-test (True or False; indicates whether user is wanting to either test the GalPage class (True) or just generate galhtml files (False)] \n [-mask (instructs script to create mask mosaics)] \n [-LS (creates LS cutouts, mosaics)] \n [-galfit (creates galfit output mosaics for all psf_indices specified in param_file)] \n [-html (creates .html scripts for galpages and homepage)]")
        sys.exit(1)
    
    if '-param_file' in sys.argv:
        p = sys.argv.index('-param_file')
        param_file = str(sys.argv[p+1])
        
    if '-test' in sys.argv:
        test = True
    else:
        test = False
           
    #create dictionary with keywords and values, from parameter .txt file

    param_dict = {}
    with open(param_file) as f:
        for line in f:
            try:
                key = line.split()[0]
                val = line.split()[1]
                param_dict[key] = val
            except:
                continue    

    #quick function to convert string input of indices (indicating wavelength bands and PSF) to a list of integers
    def str_to_list(string):
        new_list = []
        string=string.strip('[]')
        string=string.replace(',','')
        for num in range(len(string)):
            num=int(num)
            new_list.append(num)
        return new_list
    
    #another quick function to create folders if they do not already exist
    def create_folder(folder_path):
        if os.path.exists(folder_path)==False:
            os.system('mkdir '+folder_path)
            return
    
    indices_str = param_dict['psf_indices']
    psf_indices = str_to_list(indices_str)
    catalog = '/mnt/astrophysics/muchogalfit-input-cats/'+param_dict['catalog']
    dummycat = '/mnt/astrophysics/muchogalfit-input-cats/'+param_dict['dummycat']
    local_path = param_dict['local_path']
    homepage_name = param_dict['homepage_name']
    home_color = param_dict['home_color']
    website_title = param_dict['website_title'].replace('_',' ')
    
    mask_folder = local_path+param_dict['mask_folder']+'/'
    LS_mosaic_folder = local_path+param_dict['LS_mosaic_folder']+'/'
    fits_folder = local_path+param_dict['fits_folder']+'/'
    gal_mosaic_folder = local_path+param_dict['gal_mosaic_folder']+'/'
    LS_cutout_folder = local_path+param_dict['LS_cutout_folder']+'/'
    path_to_params = local_path+param_dict['path_to_params']+'/'
    path_to_galhtml = local_path+param_dict['path_to_galhtml']+'/'
    
    create_folder(LS_mosaic_folder)
    create_folder(gal_mosaic_folder)
    create_folder(LS_cutout_folder)
    create_folder(mask_folder)
    create_folder(path_to_galhtml)
    
    if test=='False':
        hp = HomePage(homepage_name=homepage_name, website_title=website_title, home_color=home_color, catalog=catalog, 
                      dummycat=dummycat, local_path=local_path, path_to_galhtml=path_to_galhtml, path_to_params=path_to_params, 
                      LS_cutout_folder=LS_cutout_folder, LS_mosaic_folder=LS_mosaic_folder, mask_folder=mask_folder, 
                      fits_folder=fits_folder, gal_mosaic_folder=gal_mosaic_folder, indices=psf_indices) 
        
        #additional memory leak precaution measures...
        if '-mask' in sys.argv:
            hp.create_model_mosaics()

        if '-LS' in sys.argv:
            hp.create_LS_figures()

        if '-galfit' in sys.argv:
            for index in psf_indices:   #only create output mosaics for the user-specified indices
                hp.create_galfit_mosaics(psf_index=index)

        if '-html' in sys.argv:
            hp.html_setup()  
        
    if test=='True':
        catalog=Table.read(catalog)
        dummycat=Table.read(dummycat)
        
        #if I am using the v2_snrcoadd.fits file, the length is 6780
        if len(catalog)>702:
            cutcat = catalog[catalog['sgacut_flag']]
        else:
            cutcat = catalog
         
        w3params_nopsf = Table.read(path_to_params+'/output_params_W3_nopsf.fits')        
        w3params_psf = Table.read(path_to_params+'/output_params_W3_psf.fits')
        rparams_nopsf = Table.read(path_to_params+'/output_params_r_nopsf.fits')
        rparams_psf = Table.read(path_to_params+'/output_params_r_psf.fits')
        
        single_galpage = GalPage(galaxy_index=0, psf_indices=psf_indices, page_name='VFID0001.html', catalog=cutcat, 
                                 dummycat=dummycat, local_path=local_path, LS_cutout_folder=LS_cutout_folder, 
                                 LS_mosaic_folder=LS_mosaic_folder, mask_folder=mask_folder, fits_folder=fits_folder, 
                                 gal_mosaic_folder=gal_mosaic_folder, w3params_nopsf=w3params_nopsf, w3params_psf=w3params_psf, 
                                 rparams_nopsf=rparams_nopsf, rparams_psf=rparams_psf, test=True)
        
        print('Single galaxy class single_galpage all set for testing. See printed instructions above.')