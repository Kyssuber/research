'''
GOAL:
- generate .html files that comprise a website interface for viewing vf galaxy images and models
OUTPUT:
- 
General instructions for opening .html locally:
COMMAND+O (that is a capital o, not zero!)
- 
NOTE --> we are currently transitioning to w1 (from r-band) for stellar disk sizes. I keep the variables here for 
completion purposes, but they are indeed commented out. I also am wanting to eliminate the no-psf mosaics on the galpages, but keep the tabulated terms.
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
#need some way to filter out the group galaxies, or at least disable the hyperlink on homepage that would otherwise lead to galfit params.
#when running .html locally, I'll want to place every folder and .html file into a specific directory. *oRgAnIzAtIoN*
'''

class HomePage():
    def __init__(self,website_title=None,home_color=None,homepage_name=None,catalog=None,local_path=None,
                 path_to_galhtml=None,path_to_params=None,LS_cutout_folder=None,LS_mosaic_folder=None,mask_folder=None,fits_folder=None,
                 gal_mosaic_folder=None,indices=None):
        
        #index = 0 for w3_nopsf, 1 for w3_psf, 2 for w1_nopsf, 3 for w1_psf
        if indices is None:
            self.indices = [0,1,2,3]
        else:
            self.indices = indices
        
        self.cat = Table.read(catalog)   #vf subsample catalog
        self.path_to_galhtml = path_to_galhtml   #path to the folder containing each individual galpage html file
        self.local_path = local_path   #path to the folder containing html directories/files when running the website locally
        self.homepage_name = homepage_name
        self.mask_mosaics = mask_folder   #path to the folder holding the .png mask mosaics
        self.LS_mosaics = LS_mosaic_folder   #contains mosaics of w3 cutout, r-band or w1 cutout, and LS cutout
        self.htmlpath = self.path_to_galhtml+self.homepage_name   #path to homepage html
        self.fits_folder = fits_folder   #path to all w3, and r-band or w1, postage stamps
        self.gal_mosaic_folder = gal_mosaic_folder   #.png mosaics of galaxy, model, and residual
        
        self.website_title = website_title
        self.home_color = home_color
        
        self.path_to_params = path_to_params   #path to galfit output parameter files
        self.LS_cutouts = LS_cutout_folder   #will contain the .png files of LS cutouts       

        #if I am using the v2_snrcoadd.fits file, the length is 6780
        if len(self.cat)>496:
            self.cutcat = self.cat[self.cat['subsample_flag']]
        else:
            self.cutcat = self.cat
        
        #define group columns
        self.group_flag = self.cutcat['group_flag']
        self.primaryGroup_flag = self.cutcat['primaryGroup_flag']
        self.group_names = self.cutcat['group_name']

        #call the remaining parameter files...
        #be sure that the size of these tables are compatible with the self.cutcat size. 
        subsample_flag_temp = self.cat['sgacut_flag'] #to trim main catalog to length 702
        flag_for_re = (self.cat.copy()[subsample_flag_temp]['massflag']) & (self.cat.copy()[subsample_flag_temp]['sSFR_flag']) & (self.cat.copy()[subsample_flag_temp]['SFRflag'])
    
        #call the remaining parameter files...
        self.w3params_nopsf = Table.read(self.path_to_params+'output_params_W3_nopsf.fits')[flag_for_re]
        self.params_w3_nopsf = self.w3params_nopsf.copy()  #will help to identify groups below...somewhere. This script is becoming unruly.
        self.w3params_psf = Table.read(self.path_to_params+'output_params_W3_psf.fits')[flag_for_re]
        #self.rparams_nopsf = Table.read(self.path_to_params+'output_params_r_nopsf.fits')[flag_for_re]
        #self.rparams_psf = Table.read(self.path_to_params+'output_params_r_psf.fits')[flag_for_re]
        self.w1params_nopsf = Table.read(self.path_to_params+'output_params_W1_nopsf.fits')[flag_for_re]
        self.w1params_psf = Table.read(self.path_to_params+'output_params_W1_psf.fits')[flag_for_re]

    def html_setup(self):
                
        with open(self.htmlpath, 'w') as html:
            
            #initialize title page
            html.write('<html><body>\n')
            html.write('<title>Virgo WISESize Project</title>\n')
            html.write(f'<body style="background-color:{self.home_color};">\n')   #be gentle with the background color (e.g., no limegreen)
            html.write('<style type="text/css">\n')
            html.write('table, td, th {padding: 5px; text-align: center; border: 2px solid black;}\n')
            html.write('p {display: inline-block;;}\n')
            html.write('</style>\n')
            html.write(f'<font size="40">{self.website_title}</font>\n')
            
            #begin to populate the table. first line is header information; loop creates the data rows.
            html.write('<table><tr><th>VFID</th><th>LS Cutout</th><th>Prefix</th><th>RA</th><th>DEC</th><th>Comments</th>\n')
            
            for i in range(len(self.cutcat)):
                html.write('<tr><td>'+self.cutcat['VFID'][i]+'</td>\n')   #add index number
                html.write('<td><img src = "' + self.LS_cutouts + self.cutcat['VFID'][i] + '-LS.jpg' + '" height="25%" width = "25%"></img></td>\n')   #cutouts will have the name VFIDxxxx-LS.png, using the v2 IDs
                
                #if galfit did not run successfully on a group galaxy then all parameters will be zeros, and heaps of trouble in terms of generating the galaxy page arises. We disable the hyperlinks in these cases. Otherwise, 
                
                if (self.params_w3_psf['xc'][i]>0):
                    if (self.primaryGroup_flag[i])|(~self.group_flag[i]):

                        #CREATE SINGLE GALPAGE using the GalPage class (see below)
                        single_galpage = GalPage(galaxy_index=i, psf_indices=self.indices, 
                                                 page_name=self.cutcat['VFID'][i]+'.html', catalog=self.cutcat, 
                                                 local_path=self.local_path, 
                                                 path_to_galhtml=self.path_to_galhtml, LS_cutout_folder=self.LS_cutouts, 
                                                 LS_mosaic_folder=self.LS_mosaics, mask_folder=self.mask_mosaics, 
                                                 fits_folder=self.fits_folder, gal_mosaic_folder=self.gal_mosaic_folder, 
                                                 w3params_nopsf=self.w3params_nopsf, w3params_psf=self.w3params_psf, 
                                                 #rparams_nopsf=self.rparams_nopsf, rparams_psf=self.rparams_psf,
                                                 w1params_nopsf = self.w1params_nopsf, w1params_psf = self.w1params_psf,
                                                 homepage_name=self.homepage_name)

                        print('Creating htmlpage for '+self.cutcat['VFID'][i])
                        single_galpage.WRITETHEGALPAGE()
                        pagename = self.VFID[i]+'.html'   #name of galaxy html page
                            
                        html.write('<td><a href='+self.path_to_galhtml+pagename+'>'+self.group_names[i]+'</a></td>\n')   #text hyperlink to galaxy page VFIDxxxx.html (pagename)
                    
                    #if neither primary or ncomp=1 galaxy, then determine which galaxy of group *is* the primary galaxy and hyperlink to that htmlpage
                    else:
                        group_rows = self.cutcat[[True if str(x)==str(self.group_names[i]) else False for x in self.group_names]]
                        primary_row = group_rows[group_rows['primaryGroup_flag']]   #all column information for the primary galaxy
                        pagename = primary_row['VFID'][0]+'.html'   #name of galaxy html page
                        print('Linking htmlpage for '+str(self.cutcat['VFID'][i])+' to '+pagename)
                        html.write('<td><a href='+self.path_to_galhtml+pagename+'>'+str(self.group_names[i])+'</a></td>\n')   #hyperlink to galaxy page VFIDxxxx.html (pagename)
                        print('<td><a href='+self.path_to_galhtml+pagename+'>'+str(self.group_names[i])+'</a></td>\n')
                
                #if galfit simply *failed* (or the primary galaxy is not a subsample member), disable hyperlink
                else:
                    print('Unable to create htmlpage for '+self.cutcat['VFID'][i])
                    html.write('<td>'+str(self.group_names[i])+'</a></td>\n')
                
                html.write('<td>'+str(self.cutcat['RA_1'][i])+'</td>\n')
                html.write('<td>'+str(self.cutcat['DEC_1'][i])+'</td>\n')

                #if the VFID (v2) is part of a group galaxy, then write in the Comments column that this particular galaxy is a member of a group
                if (self.group_flag[i]) | ('GROUP' in self.group_names[i]):
                    html.write('<td>Group Galaxy</td>\n')
                else:
                    html.write('<td>--</td>\n')
                
            html.write('</tr></table>\n')
            html.write('<br /><br />\n')
            html.write('</html></body>\n')
            html.close()
    
    #creating a few separate functions in order to independently run the PNG creation scripts - if I run them as part of the homepage class loop I use to generate every galpage.html file, then an application memory problem arises. (Why not remove the functions from the class altogether if I won't use them in the initial loop as intended? I think they are organized more nicely as part of the galpage class; and since the variable names are already entangled, I may as well not tinker any further.)           
    
    def create_LS_figures(self):
        
        for i in range(len(self.cutcat)):
            
            #if galfit ran successfully...
            if (self.params_w3_psf['xc'][i]>0):

                #I set test=True to avoid running the automatic execution of the function that creates galhtml pages
                single_galaxy = GalPage(galaxy_index=i, psf_indices=self.indices, 
                                         page_name=self.cutcat['VFID'][i]+'.html', catalog=self.cutcat, 
                                         local_path=self.local_path, 
                                         path_to_galhtml=self.path_to_galhtml, LS_cutout_folder=self.LS_cutouts, 
                                         LS_mosaic_folder=self.LS_mosaics, fits_folder=self.fits_folder, test=True)
                
                print('Creating LS cutout for '+single_galaxy.VFID)
                single_galaxy.compile_LS_cutouts()
                
                #if this galaxy is either a primary galaxy or not part of a Moustakas group, create the mosaic as well
                if (self.primaryGroup_flag[i])|(~self.group_flag[i]):
                    print('Creating LS mosaic for '+single_galaxy.VFID)
                    single_galaxy.create_LS_mosaics()

                del single_galaxy
                
    def create_galfit_mosaics(self, psf_index):
        
        index_dict = {0:'W3, no conv',
                      1:'W3, conv', 
                      2:'W1, no conv',
                      3:'W1, conv'}
        
        for n in range(len(self.cutcat)):

            #if galfit ran successfully AND this galaxy is either a primary galaxy or not part of a Moustakas group
            if (self.params_w3_psf['xc'][n]>0) & ((self.primaryGroup_flag[n])|(~self.group_flag[n])): 

                #I set test=True to avoid running the automatic execution of the function that creates galhtml pages
                single_galaxy = GalPage(galaxy_index=n, psf_indices=self.indices, page_name=self.cutcat['VFID'][n]+'.html', 
                                        catalog=self.cutcat, local_path=self.local_path, 
                                        path_to_galhtml=self.path_to_galhtml, fits_folder=self.fits_folder, 
                                        gal_mosaic_folder=self.gal_mosaic_folder, w3params_nopsf=self.w3params_nopsf, 
                                        w3params_psf=self.w3params_psf, 
                                        w1params_nopsf = self.w1params_nopsf, w1params_psf = self.w1params_psf,
                                        test=True)                
        
                single_galaxy.create_model_mosaics_names()
                print('Creating GALFIT mosaics for '+single_galaxy.VFID+f' {index_dict[psf_index]}')
                single_galaxy.create_model_mosaics(psf_index = psf_index)
        
                del single_galaxy
    
    def create_mask_mosaics(self):

        for i in range(len(self.cutcat)):
            
            #if galfit ran successfully AND this galaxy is either a primary galaxy or not part of a Moustakas group
            if (self.params_w3_psf['xc'][i]>0) & ((self.primaryGroup_flag[i])|(~self.group_flag[i])): 

                #I set test=True to avoid running the automatic execution of the function that creates galhtml pages
                single_galaxy = GalPage(galaxy_index=i, psf_indices=self.indices, page_name=self.cutcat['VFID'][i]+'.html', 
                                        catalog=self.cutcat, local_path=self.local_path, 
                                        path_to_galhtml=self.path_to_galhtml, fits_folder=self.fits_folder, 
                                        gal_mosaic_folder=self.gal_mosaic_folder, mask_folder=self.mask_mosaics, test=True)   
                
                print('Creating mask mosaic for '+single_galaxy.VFID)
                single_galaxy.create_mask_mosaics()

            del single_galaxy
                
class GalPage():
    def __init__(self,galaxy_index=None, psf_indices = [0,1,2,3], page_name=None, catalog=None, 
                 local_path=None, path_to_galhtml=None, LS_cutout_folder=None, LS_mosaic_folder=None, mask_folder=None, 
                 fits_folder=None, gal_mosaic_folder=None, w3params_nopsf=None, w3params_psf=None, rparams_nopsf=None, 
                 rparams_psf=None, w1params_nopsf=None, w1params_psf=None, homepage_name=None, test=False):

        self.galaxy_index = int(galaxy_index)
        self.psf_indices = psf_indices
        self.catalog = catalog
        self.cutcat = catalog
        
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
        #self.rparams_nopsf = rparams_nopsf
        #self.rparams_psf = rparams_psf
        self.w1params_nopsf = w1params_nopsf
        self.w1params_psf = w1params_psf
        
        self.RA = str(self.cutcat['RA_1'][self.galaxy_index])
        self.DEC = str(self.cutcat['DEC_1'][self.galaxy_index])
        self.VFID = str(self.cutcat['VFID'][self.galaxy_index])
        self.objname = str(self.cutcat['group_name'][self.galaxy_index])   #group name is proxy for objname; if galaxy is not in group, group_name=objname
        
        self.page_name = self.VFID+'.html'   #name of galaxy html page
        
        if not test:
            self.gal_htmlpath = self.path_to_galhtml+self.page_name
        if test:
            self.gal_htmlpath = os.getenv("HOME")+'/'+self.page_name
        
        #find number of galaxies sharing the same group_name
        self.ncomp = len(np.where(self.cutcat['group_name']==self.cutcat['group_name'][self.galaxy_index])[0])
        
        #need w3 header information and w3, r-band image data
        path_to_w3 = glob.glob(self.fits_folder+self.objname+'-custom-image-W3.fits')[0]
        #path_to_r = glob.glob(self.fits_folder+self.objname+'-custom-image-r.fits')[0]
        path_to_w1 = glob.glob(self.fits_folder+self.objname+'-custom-image-W1.fits')[0]

        self.wise_im, self.wise_header = fits.getdata(path_to_w3, header=True)
        #self.r_im, self.r_header = fits.getdata(path_to_r, header=True)
        self.w1_im, self.w1_header = fits.getdata(path_to_w1, header=True)
        
        #WCS header information
        self.wcs_w3 = WCS(self.wise_header)
        #self.wcs_r = WCS(self.r_header)
        self.wcs_w1 = WCS(self.w1_header)

        #defining the mask paths
        try:
            self.w3mask_path = glob.glob(self.fits_folder+self.objname+'-custom-image-wise-mask.fits')[0]
            #self.rmask_path = glob.glob(self.fits_folder+self.objname+'-custom-image-r-mask.fits')[0]
            self.w3_mask = fits.getdata(self.w3mask_path)
            #self.r_mask = fits.getdata(self.rmask_path)
            self.w1_mask = self.w3_mask.copy()   #pixelscales are the same so masks will be also
            
            #define boolean masks for image scaling purposes (so those darn bright stars do not dictate the norm vmin, vmax)
            #when multiplied by image, only the unmasked object (only -- and all of -- the central galaxy, ideally) remains
            self.w3_mask_bool = ~(self.w3_mask>0)
            #self.r_mask_bool = ~(self.r_mask>0)
            self.w1_mask_bool = ~(self.w1_mask>0)
       
        except:
            print(self.objname+' has no mask images.')
        
        #if not testing the various functions on one galaxy (test==True), then run only the functions that are required for (1) variables and (2) the actual galaxy html pages. Think of these functions as the 'defaults' that enable the user to then generate the PNG files at their liberty.
        if test==False:
            self.create_model_mosaics_names()
            self.tabulate_parameters()
        
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
        #r_scaled, footprint = reproject_interp((self.r_im, self.r_header), self.wise_header)
        titles = ['W3 Image', 'W1 Image', 'LS Image']
        images = [self.wise_im, self.w1_im, self.filename_LS]
        
        try:
            bool_masks = [self.w3_mask_bool, self.w1_mask_bool, None]
        except:
            bool_masks = [np.zeros((len(self.wise_im),len(self.wise_im)))+1, np.zeros((len(self.w1_im),len(self.w1_im)))+1, None]
            print(f'{self.VFID} has no mask images.')
        
        plt.figure(figsize=(12,6))
        for i,im in enumerate(images):
            plt.xlabel('RA')
            if (i==0)|(i==1):
                if i==0:
                    plt.subplot(1,len(images),i+1,projection = self.wcs_w3)
                if i==1:
                    plt.subplot(1,len(images),i+1,projection = self.wcs_w1)
                try:
                    norm = simple_norm(images[i]*bool_masks[i],stretch='asinh',max_percent=99.9)
                    plt.imshow(images[i],origin='lower',cmap='viridis',norm=norm)
                except:
                    norm = simple_norm(self.w1_im,stretch='asinh',max_percent=99.9)
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
        self.file_w1_nopsf = self.fits_folder+self.objname+'-W1-out1.fits'
        self.file_w1_psf = self.fits_folder+self.objname+'-W1-out2.fits'
        
        self.mosaic_names = [self.file_w3_nopsf,self.file_w3_psf,self.file_w1_nopsf,self.file_w1_psf]
        
        self.psf_dictionary = {0:self.file_w3_nopsf,
                            1:self.file_w3_psf,
                            2:self.file_w1_nopsf,
                            3:self.file_w1_psf}
        
        #is os.path.exists() --> if not exist, set equal to None
        self.models = []
        self.residuals = []
        self.psf_indices_galaxy = []

        for index in range(4):
            if os.path.exists(self.psf_dictionary[index]):
                self.models.append(fits.getdata(self.psf_dictionary[index],2))
                self.residuals.append(fits.getdata(self.psf_dictionary[index],3))
                self.psf_indices_galaxy.append(index)
            else:
                self.models.append(None)
                self.residuals.append(None)
                self.psf_indices_galaxy.append(None)
        self.psf_indices_galaxy = np.asarray(self.psf_indices_galaxy)
        
        self.pngnames = [self.gal_mosaic_folder+self.VFID+'-'+'galfit-model-w3-nopsf.png',
                   self.gal_mosaic_folder+self.VFID+'-'+'galfit-model-w3-psf.png',
                   self.gal_mosaic_folder+self.VFID+'-'+'galfit-model-w1-nopsf.png',
                   self.gal_mosaic_folder+self.VFID+'-'+'galfit-model-w1-psf.png']
        
        #print('For self.create_model_mosaics(index), index=0 is w3_nopsf, 1 is w3_psf, 2 is w1_nopsf, 3 is w1_psf.')
        
    def create_model_mosaics(self, psf_index, percentile1=.5, percentile2=99.9, p1residual=5, p2residual=99, cmap='viridis'):
        
        '''
        ARGS:
        index = 0 for w3_nopsf, 1 for w3_psf, 2 for w1_nopsf, 3 for w1_psf
        percentile1 = min percentile for stretch of image and model
        percentile2 = max percentile for stretch of image and model
        p1residual = min percentile for stretch of residual
        p2residual = max percentile for stretch of residual
        cmap = colormap, default is viridis
        ''' 
        
        #create boolean masks for scaling purposes (prevents domination of prominent stars or other artifacts)
        try:
            bool_masks = [self.w3_mask_bool, self.w1_mask_bool]
        except:
            #if no im mask, I just create a simple nxn matrix of 1s, so multiplying by bool_mask does not affect the image
            bool_masks = [np.zeros((len(self.wise_im),len(self.wise_im)))+1, np.zeros((len(self.w1_im),len(self.w1_im)))+1]
        
        #for index in self.psf_indices:
        if psf_index<2:   #w3 is index=0 or index=1
            images = [self.wise_im,self.models[psf_index],self.residuals[psf_index],self.residuals[psf_index]]
            bool_mask = bool_masks[0]
        if psf_index>=2:   #w1 is index=2 or index=3
            images = [self.w1_im,self.models[psf_index],self.residuals[psf_index],self.residuals[psf_index]]
            bool_mask = bool_masks[1]
        
        titles = ['Image','Model','Residual (img stretch)','Residual (hard stretch)']
        
        if images[1] is not None:   #if model not None, then galfit ran correctly
            v1 = [scoreatpercentile(images[0]*bool_mask,percentile1),
                scoreatpercentile(images[0]*bool_mask,percentile1),
                scoreatpercentile(images[0]*bool_mask,percentile1),
                scoreatpercentile(images[2]*bool_mask,p1residual)]
            v2 = [scoreatpercentile(images[0]*bool_mask,percentile2),
                scoreatpercentile(images[0]*bool_mask,percentile2),
                scoreatpercentile(images[0]*bool_mask,percentile2),
                scoreatpercentile(images[2]*bool_mask,p2residual)]

            norms = [simple_norm(images[0]*bool_mask,'asinh',max_percent=percentile2,min_cut=v1[0],max_cut=v2[0]),
                   simple_norm(images[0]*bool_mask,'asinh',max_percent=percentile2,min_cut=v1[1],max_cut=v2[1]),
                   simple_norm(images[0]*bool_mask,'asinh',max_percent=percentile2,min_cut=v1[2],max_cut=v2[2]),
                   simple_norm(images[2]*bool_mask,'linear',max_percent=p2residual,min_cut=v1[3],max_cut=v2[3])]
        
        else:
            print(f'GALFIT did not run correctly for {self.VFID} - no model...')
            
            v1 = [scoreatpercentile(images[0]*bool_mask,percentile1),
                None, None, None]
            v2 = [scoreatpercentile(images[0]*bool_mask,percentile2),
                None, None, None]

            norms = [simple_norm(images[0]*bool_mask,'asinh',max_percent=percentile2,min_cut=v1[0],max_cut=v2[0]),
                   None, None, None]

        plt.figure(figsize=(14,6))
        plt.subplots_adjust(wspace=.0)
        for i,im in enumerate(images): 
            ax = plt.subplot(1,4,i+1,projection=self.wcs_w3)
            if v1[i] is None:
                plt.imshow(np.zeros((len(self.wise_im),len(self.wise_im))),origin='lower',cmap='gray')
            else:
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

    def create_mask_mosaics(self, percentile1=.5, percentile2=99.8, cmap='viridis'):

        titles = ['W3 Image', 'W3 Mask', 'W1 Image', 'W1 Mask']
        images = [self.wise_im, self.w3_mask, self.w1_im, self.w1_mask]
        
        #norms for images but not for masks
        v1 = [scoreatpercentile(images[0],percentile1),0,scoreatpercentile(images[2],percentile1),0]
        v2 = [scoreatpercentile(images[0],percentile2),1,scoreatpercentile(images[2],percentile2),1]
        norms = [simple_norm(images[0],'asinh',max_percent=percentile2,min_cut=v1[0],max_cut=v2[0]),None, 
                 simple_norm(images[2],'asinh',max_percent=percentile2,min_cut=v1[2],max_cut=v2[2]),None]
              
        plt.figure(figsize=(14,6))
        plt.subplots_adjust(wspace=.0)
        for i,im in enumerate(images): 
            if i<=1:
                ax = plt.subplot(1,4,i+1,projection=self.wcs_w3)                    
            if i>=2:
                ax = plt.subplot(1,4,i+1,projection=self.wcs_w1)
            if (i==0)|(i==2):
                plt.imshow(im, origin='lower', cmap=cmap, norm=norms[i])
            else:
                plt.imshow(im, origin='lower', cmap=cmap, vmin=v1[i], vmax=v2[i])
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
        
        #list of each table of parameters
        params_list = [self.w3params_nopsf,self.w3params_psf,self.w1params_nopsf,self.w1params_psf]
        
        self.page_params = []
        self.page_extparams = []   #initiate extparams list
        
        if self.cutcat['primaryGroup_flag'][self.galaxy_index]:
            groupname = self.cutcat['group_name'][self.galaxy_index]

        #create list of every parameter row for the given indices
        for index in self.psf_indices:
            params=params_list[index]   #selects the correct parameter table from params_list above
            param_row=params[params['VFID']==self.VFID]   #find the correct parameter row corresponding to the galaxy VFID
            self.page_params.append(param_row)
            single_extparams=[]   #create empty list for external galaxies for single psf_index

            if self.cutcat['primaryGroup_flag'][self.galaxy_index]:   #if primary group galaxy, then find its subsample pals (if any)
                indices = np.where((self.cutcat['group_name']==groupname)&(~self.cutcat['primaryGroup_flag']))[0]
                for num in indices: #for every external galaxy, add row of parameters
                    param_row=params[num]  #find correct row
                    single_extparams.append(param_row)

            self.page_extparams.append(single_extparams)  #will comprise list of external galaxy lists. LISTCEPTION. 
    
    #create VFIDxxxx.html for the galaxy!
    def WRITETHEGALPAGE(self):
        
        title_dictionary = {0:'w3 noConv',
                            1:'w3 Conv',
                            2:'W1 noConv',
                            3:'W1 Conv'}
        
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
                if self.models[index] is not None:
                    mosaic_path = self.pngnames[index]
                    params = self.page_params[index]
                    html.write(f'<font size="30"> GALFIT Output Mosaic {title_dictionary[index]}:</font><br /> \n')
                    html.write('<div class='+'"'+'img-container'+'"> <!-- Block parent element --> <img src='+'"'+mosaic_path+'" height="70%" width="90%" /><br /> \n')

                    html.write('<table><tr><th>VFID</th><th>Type</th><th>xc</th><th>xc_err</th><th>yc</th><th>yc_err</th><th>mag</th><th>mag_err</th><th>Re</th><th>Re_err</th><th>nser</th><th>nser_err</th><th>BA</th><th>BA_err</th><th>PA</th><th>PA_err</th><th>err_flag</th></tr> \n')

                    html.write('<tr><td>'+self.VFID+'</td> \n')
                    html.write('<td>Primary</td> \n')
                    for p in range(1,15):   #for every parameter value in a row
                        html.write(f'<td>{params[0][p]}</td> \n')   #0th row, which does not change as there is only one row here
                    html.write(f'<td>{params[0][17]}</td></tr> \n')   #error flag

                    #add the external galaxy parameters, if any.
                    if len(self.page_extparams)>0:

                        param_ext=self.page_extparams[index]  #isolate band&psf parameters of the external galaxies (if there are multiple galaxies, then there will be multiple lists for one index/central galaxy). 
                        for num in range(len(param_ext)): 
                            single_param_ext = param_ext[num]  #isolate the self.ext_list[num] external galaxy parameters
                            html.write('<tr><td>'+single_param_ext['VFID']+'</td> \n')  #VFID of external galaxy
                            html.write('<td>External</td> \n')
                            for p in range(1,15):
                                html.write(f'<td>{single_param_ext[p]}</td> \n')
                            if num==self.ncomp-1:
                                html.write(f'<td>{params[17]}</td></tr> \n')

                    #I want a table for every nopsf/psf pair, so the table will finish if index is psf w3 or psf w1 (OR neither, in the case that I don't include any psf bands...). I would attempt to account for the case where I am only including nopsf tables, but I do not anticipate this condition being met. Ever.
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
    
    print("""If running in terminal, use -h or --help for the list of arguments. Otherwise...*
    
    USAGE (iPython/Jupyter):
    ---Running this program automatically initiates the HomePage class (hp)
    ---if the -test arg is False, then replace single_galpage with hp in order to generate all images and such:
            --hp.html_setup() --> create .html homepages for all galaxies in the VF subsample; also initiates GalPage 
            class(single_galpage) for every galaxy in a loop, creating all relevant folders and files.
            --hp.create_LS_figures()
            --hp.create_galfit_mosaics()
            --hp.create_mask_mosaics(psf_index) --> (w3, nopsf), 1 (w3, psf), 2 (w1, nopsf), 3 (w1, psf)
    ---If the -test arg is True, then the user is wanting to test the GalPage class (this is a stategic idea, 
            since this class contains the bulk of the functions required for the script to run successfully). 
            In this case, the script automatically defines a single galpage class (single_galpage; 
            uses index 0 for a random test galaxy), with which the user can test the following methods:
            --single_galpage.compile_LS_cutouts()
            --single_galpage.create_LS_mosaics()
            --single_galpage.create_model_mosaics_names()
            --single_galpage.create_model_mosaics(psf_index) --> 0 (w3, nopsf), 1 (w3, psf), 2 (w1, nopsf), 3 (w1, psf)
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
    
    #create the MASTER DIRECTORY. MWUHAHAHAHAHA.
    create_folder(local_path)
    #moves fits folder (created using another script) into local path, if not already done.
    if os.path.exists(fits_folder)==False:
        os.system('mv /mnt/astrophysics/kconger_wisesize/all_input_fits '+local_path)
    #create rest of directories if not already present
    create_folder(LS_mosaic_folder)
    create_folder(gal_mosaic_folder)
    create_folder(LS_cutout_folder)
    create_folder(mask_folder)
    create_folder(path_to_galhtml)
    create_folder(path_to_params)
    
    if '-test' not in sys.argv:
        print('test is false')
        hp = HomePage(homepage_name=homepage_name, website_title=website_title, home_color=home_color, catalog=catalog, 
                      local_path=local_path, path_to_galhtml=path_to_galhtml, path_to_params=path_to_params, 
                      LS_cutout_folder=LS_cutout_folder, LS_mosaic_folder=LS_mosaic_folder, mask_folder=mask_folder, 
                      fits_folder=fits_folder, gal_mosaic_folder=gal_mosaic_folder, indices=psf_indices) 
        
        #additional memory leak precaution measures...
        if '-mask' in sys.argv:
            hp.create_mask_mosaics()

        if '-LS' in sys.argv:
            hp.create_LS_figures()

        if '-galfit' in sys.argv:
            for index in psf_indices:   #only create output mosaics for the user-specified indices
                hp.create_galfit_mosaics(psf_index=index)

        if '-html' in sys.argv:
            hp.html_setup()  
        
    if '-test' in sys.argv:
        print('test is true')
        catalog=Table.read(catalog)
        
        #if I am using the v2_snrcoadd.fits file, the length is 6780
        if len(catalog)>702:
            cutcat = catalog[catalog['sgacut_flag']]
        else:
            cutcat = catalog
         
        w3params_nopsf = Table.read(path_to_params+'/output_params_W3_nopsf.fits')        
        w3params_psf = Table.read(path_to_params+'/output_params_W3_psf.fits')
        #rparams_nopsf = Table.read(path_to_params+'/output_params_r_nopsf.fits')
        #rparams_psf = Table.read(path_to_params+'/output_params_r_psf.fits')
        w1params_nopsf = Table.read(path_to_params+'/output_params_W1_nopsf.fits')
        w1params_psf = Table.read(path_to_params+'/output_params_W1_psf.fits')
        
        single_galpage = GalPage(galaxy_index=59, psf_indices=psf_indices, page_name='VFID0281.html', catalog=cutcat, 
                                 local_path=local_path, LS_cutout_folder=LS_cutout_folder, 
                                 LS_mosaic_folder=LS_mosaic_folder, mask_folder=mask_folder, fits_folder=fits_folder, 
                                 gal_mosaic_folder=gal_mosaic_folder, w3params_nopsf=w3params_nopsf, w3params_psf=w3params_psf, 
                                 w1params_nopsf=w1params_nopsf, w1params_psf=w1params_psf, test=True)
        
        print('Single galaxy class single_galpage all set for testing. See printed instructions above.')
