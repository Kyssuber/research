'''
GOAL:
- generate .html files that comprise a website interface for viewing vf galaxy images and models
OUTPUT:
- 
General instructions for opening .html locally:
COMMAND+O
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
'''
#NOTES:
#need some way to filter out the group galaxies, or at least disable the hyperlink on homepage that would otherwise lead to galfit params. also dummycat...
#when running .html locally, I'll want to place every folder and .html file into a specific directory. *oRgAnIzAtIoN*
'''

class HomePage(homepage_name='home_local.html'catalog=None,dummycat=None,local_path=None,path_to_params=None,
               LS_cutout_folder=None,LS_mosaic_folder=None,mask_folder=None,fits_folder=None,gal_mosaic_folder=None):
    
    def __init__(self):
        
        self.cat = catalog   #vf subsample catalog
        self.dummycat = dummycat   #catalog flagging galaxies with 2+ Sersic objects in their cutouts
        
        self.local_path = local_path   #path to the folder containing html folder when running the website locally
        self.homepage_name = homepage_name
        self.mask_mosaics = self.local_path+mask_folder   #path to the folder holding the mask mosaics
        self.LS_mosaics = self.local_path+LS_mosaic_folder   #contains mosaics of w3 cutout, r-band cutout, and LS cutout
        self.htmlpath = self.local_path+self.homepage_name   #path to homepage html
        self.fits_folder = fits_folder   #path to all w3, r-band postage stamps
        self.gal_mosaic_folder = gal_mosaic_folder   #.png mosaics of galaxy, model, and residual
        
        self.path_to_params = path_to_params   #path to galfit output parameter files
        self.LS_cutouts = self.local_path+LS_cutout_folder   #will contain the .png files of LS cutouts        
        
        self.group_flag = (self.path_to_params['xc'] == 0)   #the first run of galfit did not produce any params if the filename include the word 'group.' this flag will filter out these galaxies (or simply leave them be once problem is resolved).
        
        self.cutcat = self.cat[self.group_flag]     
    
    def html_setup(self):
        
        with open(self.htmlpath, 'w') as html:
            
            #initialize title page
            html.write('<html><body>\n')
            html.write('<title>Virgowise Project</title>\n')
            html.write('<body style="background-color:powderblue;">\n')   #be gentle with the background color (e.g., no limegreen)
            html.write('<style type="text-css">\n>')
            html.write('table, td, th {padding: 5px; text-align: center; border: 2px solid black;}\n')
            html.write('p {display: inline-block;;}\n')
            html.write('</style>\n')
            html.write('<font size="40">GALFIT Data for Virgowise VF Galaxies</font>\n')
            
            #begin and populate the table. first line is header information; loop creates the data rows.
            html.write('<table><tr><th>Index</th><th>LS Cutout</th><th>Prefix</th><th>RA</th><th>DEC</th><th>Comments</th>')
            
            for i in range(0,len(self.cutcat)):
                html.write('<tr><td>'+str(i)+'</td>\n')   #add index number
                
                html.write('<td><img src = "' + self.LS_cutouts + str(self.cutcat['VFID'][i]) + '-LS.jpg' + '" height="50%" width = "50%"></img></td>\n')   #cutouts will have the name VFIDxxxx-LS.png, using the v2 IDs
                
                #CREATE SINGLE GALPAGE using the GalPage class (see below)
                single_galpage = GalPage(galaxy_index=i, page_name=self.cutcat['VFID'][i]+'.html', catalog=self.cutcat, dummycat=self.dummycat, local_path=self.local_path, LS_mosaic_folder=self.LS_mosaics, mask_folder=self.mask_mosaics, fits_folder = self.fits_folder, gal_mosaic_folder=self.gal_mosaic_folder)
                
                pagename = single_galpage.page_name
                
                hmtl.write('<td><a href='+self.local_path+pagename+'>'+str(self.cutcat['prefix'][i])+'</a></td>\n')   #prefix text hyperlink to galaxy page VFIDxxxx.html (pagename)
                
                html.write('<td>'+str(self.cutcat['RA'][i])+'</td>\n')
                html.write('<td>'+str(self.cutcat['DEC'][i])+'</td>\n')

                #if the VFID (v2) is in the dummycat central galaxy column, then write in the Comments column that this particular galaxy is a member of a "Moustakas group"
                if self.cutcat['VFID'] in dummycat['central galaxy']:
                    html.write('<td>Group Galaxy</td>')
                else:
                    html.write('<td>--</td>')

            html.write('</tr></table>\n')
            html.write('<br /><br />\n')
            html.write('</html></body>\n')
            html.close()

class GalPage(galaxy_index=None, page_name=None, catalog=None, dummycat=None, local_path=None, LS_mosaic_folder=None, mask_folder=None, fits_folder=None, gal_mosaic_folder=None):
    
    def __init__:
        self.galaxy_index = int(galaxy_index)
        self.catalog = catalog
        self.dummycat = dummycat
        
        self.local_path = local_path
        self.LS_mosaics = LS_mosaic_folder
        self.mask_mosaics = mask_folder
        self.fits_folder = fits_folder
        self.gal_mosaic_folder = gal_mosaic_folder
        
        self.RA = self.cutcat['RA_1'][self.galaxy_index]
        self.RA = str(self.RA[0])
        self.DEC = self.cutcat['DEC_1'][self.galaxy_index]
        self.DEC = str(self.DEC[0])
        self.VFID = self.cutcat['VFID'][self.galaxy_index]
        self.VFID = str(self.VFID[0])
        
        self.page_name = self.VFID+'.html'   #name of galaxy html page
        
        self.objname = self.cutcat['objname'][self.galaxy_index]
        self.objname = str(self.objname[0])
        
        #need w3 header information and w3, r-band image data
        path_to_w3 = glob.glob(self.fits_folder+self.objname+'*'+'W3.fits')[0]
        path_to_r = glob.glob(self.fits_folder+self.objname+'*'+'r.fits')[0]
        self.wise_im, self.wise_header = fits.getdata(path_to_w3, header=True)
        self.r_im, self.r_header = fits.getdata(path_to_r, header=True)
        
        self.wcs_w3 = WCS(self.wise_header)
        
        self.create_model_mosaics_names()
        
    #download scaled RGB jpg of VFIDxxxx galaxy; place in self.filename_LS
    def compile_LS_cutouts(self):
        
        pixscale=1   #standard
        pscale=np.abs(float(self.wise_header['CD1_1']))   #grab transformation matrix of WISE image
        xsize=np.abs(int(self.wise_header['NAXIS1']))   #grab length of WISE image
        xsize_arcsec=pscale*3600*xsize   #length convert to arcseconds
        imsize=int(xsize_arcsec/pixscale)   #convert length to an integer
        imsize=str(imsize)   #convert integer length to a...string
        
        self.filename_LS = self.LS_cutouts+VFID+'-LS.jpg'
        
        image_url = f'https://www.legacysurvey.org/viewer/cutout.jpg?ra={ra}&dec={dec}&layer=ls-dr9&size={imsize}&pixscale={1}'
        if os.path.exists(self.filename_LS):
            os.remove(self.filename_LS)
            image = wget.download(image_url,out=self.filename_LS)
        else:
            image = wget.download(image_url,out=self.filename_LS)
    
    def create_LS_mosaics(self):
        titles = ['W3 Image', 'r-band Image', 'LS Image']
        images = [self.wise_im, self.r_im, self.filename_LS]
        
        plt.figure(figsize=(12,6))
        for i,im in enumerate(images):
            plt.xlabel('RA')
            if (i==0)|(i==1):
                plt.subplot(1,len(images),i+1,projection = self.wcs_w3)
                norm = simple_norm(images[i],stretch='asinh',max_percent=99.5)
                plt.imshow(images[i],origin='lower',cmap='viridis',norm=norm)
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
        
        if os.path.exists(self.LS_mosaics):
            os.remove(self.LS_mosaics)
            plt.savefig(self.LS_mosaics+self.VFID+'-LS-mosaic.png')
        else:
            plt.savefig(self.LS_mosaics+self.VFID+'-LS-mosaic.png')
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
        
        self.pngnames = [self.galname+'-'+'galfit-model-w3-nopsf.png',
                   self.galname+'-'+'galfit-model-w3-psf.png',
                   self.galname+'-'+'galfit-model-r-nopsf.png',
                   self.galname+'-'+'galfit-model-r-psf.png']
        
        print('For self.create_model_mosaics(index), index=0 is w3_nopsf, 1 is w3_psf, 2 is r_nopsf, 3 is r_psf.')

    def create_model_mosaics(self, index, percentile1=.5, percentile2=99.5, p1residual=5, p2residual=99, cmap='viridis'):
        
        '''
        ARGS:
        index = 0 for w3_nopsf, 1 for w3_psf, 2 for r_nopsf, 3 for r_psf
        percentile1 = min percentile for stretch of image and model
        percentile2 = max percentile for stretch of image and model
        p1residual = min percentile for stretch of residual
        p2residual = max percentile for stretch of residual
        cmap = colormap, default is viridis
        ''' 
        
        if index<2:   #w3 is index=0 or index=1
            images = [self.wise_im,self.models[index],self.residuals[index],self.residuals[index]]
        if index>=2:   #r-band is index=2 or index=3
            images = [self.r_im,self.models[index],self.residuals[index],self.residuals[index]]
        titles = ['Image','Model','Residual (img stretch)','Residual (res stretch)']
        
        v1 = [scoreatpercentile(images[0],percentile1),
            scoreatpercentile(images[0],percentile1),
            scoreatpercentile(images[0],percentile1),
            scoreatpercentile(images[3],p1residual)]
        v2 = [scoreatpercentile(images[0],percentile2),
            scoreatpercentile(images[0],percentile2),
            scoreatpercentile(images[0],percentile2),
            scoreatpercentile(images[3],p2residual)]
        
        norms = [simple_norm(images[0],'asinh',max_percent=percentile2),
               simple_norm(images[0],'asinh',max_percent=percentile2),
               simple_norm(images[0],'asinh',max_percent=percentile2),
               simple_norm(images[0],'linear',max_percent=p2residual)]
               
        plt.figure(figsize=(14,6))
        plt.subplots_adjust(wspace=.0)
        for i,im in enumerate(images): 
            ax = plt.subplot(1,4,i+1,projection=self.wcs_w3)
            plt.imshow(im,origin='lower',cmap=cmap,vmin=v1[i],vmax=v2[i],norm=norms[i])
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
        
        
        #self.gal_mosaic_folder
        
        
        #w3 nopsf
        #w3 psf
        #r-band nopsf
        #r-band psf
        return
    
    
    
    def create_mask_mosaics(self):
        #w3 nopsf
        #w3 psf
        #r-band nopsf
        #r-band psf
        return
    
    def tablulate_parameters(self):
        #w3 nopsf
        #w3 psf
        #r-band nopsf
        #r-band psf
        return
    
    def WRITETHEGALPAGE(self):
        return


def build_gal_html(i, vf_sample, galpath, htmlpath, galfit_params_nopsf, galfit_params_psf, dummycat = None, fixed_galfit_params_nopsf = None, fixed_galfit_params_psf = None):
    
    #VFID (v2) html name of webpage
    htmlname = str(vf_sample[i]['VFID']) + '.html'
    
    #full local path to where this .html is stored
    full_htmlpath = htmlpath + htmlname
    
    #begin to write text in full_htmlpath file
    with open(full_htmlpath, 'w') as html:
        
        #I think much of this is style preparation, similar to the initial lines in a LaTeX document
        html.write('<html><body>\n')
        html.write('<title>' + str(vf_sample['prefix'][i]) + '<\title>\n')  #title that appears on the browser tab
        html.write('<style type="text/css">\n')
        html.write('.img-container{text-align: left;}')
        html.write('table, td, th {padding: 5px; text-align: center; border: 1px solid black;}\n')
        html.write('p {display: inline-block;;}\n')
        html.write('</style>\n')
        
        #name of central galaxy
        html.write('<fontsize="40">Central Galaxy: ' + str(vf_sample['prefix'][i]) + '</font><br />\n')
        
        #add hyperlinks that return user to either homepage, next galaxy (if i != 0), or prev galaxy (if i != len(sample) - 1)
        html.write('<a href=main.html>Return to Homepage</a></br />\n')
        if i != len(vf_sample) - 1:
            html.write('<a href='+str(vf_sample['VFID'][i+1])+'.html>Next Galaxy</a></br />\n')
        if i != 0:
            html.write('<a href='+str(vf_sample['VFID'][i-1])+'.html>Previous Galaxy</a></br />\n')
        
        '''
        unclear how we will treat postage stamps with 2+ Sersic obj, so this and following relevant sections will be commented out until further notice.
        
        #if zero, then ncomp = 1
        ncomp = len(np.where(dummycat['central galaxy'] == vf_sample['VFID'][i])[0]) + 1
        '''
        
        ncomp = 1




if __name__ == '__main__':    
    
    if '-h' in sys.argv or '--help' in sys.argv:
        print("Usage: %s [-param_file (name of parameter file, no single or double quotation marks)]")
        sys.exit(1)
    
    if '-param_file' in sys.argv:
        p = sys.argv.index('-param_file')
        param_file = str(sys.argv[p+1])
    
    
    homedir = os.getenv("HOME")
       
    #create dictionary with keywords and values, from parameter .txt file

    param_dict = {}
    with open(homedir+'/'+param_file) as f:
        for line in f:
            try:
                key = line.split()[0]
                val = line.split()[1]
                param_dict[key] = val
            except:
                continue    
    
    #now...extract parameters and assign to relevantly-named variables
    
    #catalogs and pathnames
    cat_path = param_dict['vf_catalog']
    cat = Table.read(cat_path)
    dummycat_path = param_dict['dummycat']
    dummycat = Table.read(dummycat_path)
    directory = param_dict['directory']
    psf_filepath = param_dict['psf_filepath']
    psf_filepath = homedir+psf_filepath
    gal_output_path = param_dict['gal_output_path']
