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
from reproject import reproject_interp

'''
#NOTES:
#need some way to filter out the group galaxies, or at least disable the hyperlink on homepage that would otherwise lead to galfit params. also dummycat...
#when running .html locally, I'll want to place every folder and .html file into a specific directory. *oRgAnIzAtIoN*
'''

class HomePage(homepage_name='home_local.html'catalog=None,dummycat=None,local_path=None,path_to_params=None,
               LS_cutout_folder=None,LS_mosaic_folder=None,mask_folder=None,fits_folder=None,gal_mosaic_folder=None):
    
    def __init__(self):
        
        #index = 0 for w3_nopsf, 1 for w3_psf, 2 for r_nopsf, 3 for r_psf --> add this key to params.txt
        self.indices = indices
        
        self.cat = catalog   #vf subsample catalog
        self.dummycat = dummycat   #catalog flagging galaxies with 2+ Sersic objects in their cutouts
        
        self.local_path = local_path   #path to the folder containing html folder when running the website locally
        self.homepage_name = homepage_name
        self.mask_mosaics = self.local_path+mask_folder   #path to the folder holding the .png mask mosaics
        self.LS_mosaics = self.local_path+LS_mosaic_folder   #contains mosaics of w3 cutout, r-band cutout, and LS cutout
        self.htmlpath = self.local_path+self.homepage_name   #path to homepage html
        self.fits_folder = fits_folder   #path to all w3, r-band postage stamps
        self.gal_mosaic_folder = gal_mosaic_folder   #.png mosaics of galaxy, model, and residual
        
        self.path_to_params = path_to_params   #path to galfit output parameter files
        self.LS_cutouts = self.local_path+LS_cutout_folder   #will contain the .png files of LS cutouts        
        
        #call one of the parameter textfiles in order to identify galaxies belonging to groups.
        self.w3params_nopsf = Table.read(self.path_to_params+'/output_params_W3_nopsf.fits')        
        self.params_w3_nopsf = Table.read(self.w3params_nopsf)

        #if I am using the v2_snrcoadd.fits file, the length is 6780
        if len(self.cat)>702:
            self.cutcat = self.cat[self.cat['sgacut_flag']]
        else:
            self.cutcat = self.cat
         
        #call the remaining parameter files...
        self.w3params_psf = Table.read(self.path_to_params+'/output_params_W3_psf.fits')
        self.rparams_nopsf = Table.read(self.path_to_params+'/output_params_r_nopsf.fits')
        self.rparams_psf = Table.read(self.path_to_params+'/output_params_r_psf.fits')
        
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
                
                #if galaxy is part of a group, then all parameters will be zeros, and heaps of trouble in terms of generating the galaxy page arises. As such, we effectively disable the hyperlinks in these cases.
                if self.params_w3_nopsf['xc']>0:
                    
                    #CREATE SINGLE GALPAGE using the GalPage class (see below)
                    single_galpage = GalPage(galaxy_index=i, psf_indices=self.indices, 
                                             page_name=self.cutcat['VFID'][i]+'.html', catalog=self.cutcat, 
                                             dummycat=self.dummycat, local_path=self.local_path, 
                                             LS_mosaic_folder=self.LS_mosaics, mask_folder=self.mask_mosaics, 
                                             fits_folder=self.fits_folder, gal_mosaic_folder=self.gal_mosaic_folder, 
                                             w3params_nopsf=self.w3params_nopsf, w3params_psf=self.w3params_psf, 
                                             rparams_nopsf=self.rparams_nopsf, rparams_psf=self.rparams_psf)

                    pagename = single_galpage.page_name

                    hmtl.write('<td><a href='+self.local_path+pagename+'>'+str(self.cutcat['prefix'][i])+'</a></td>\n')   #text hyperlink to galaxy page VFIDxxxx.html (pagename)
                
                else:
                    hmtl.write('<td>'+str(self.cutcat['prefix'][i])+'</a></td>\n')
                
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

class GalPage(galaxy_index=None, page_name=None, catalog=None, dummycat=None, local_path=None, LS_mosaic_folder=None, mask_folder=None, fits_folder=None, gal_mosaic_folder=None, w3params_nopsf=None, w3params_psf=None, rparams_nopsf=None, rparams_psf=None):
    
    def __init__:
        self.galaxy_index = int(galaxy_index)
        self.catalog = catalog
        self.dummycat = dummycat
        
        self.local_path = local_path
        self.LS_mosaics = LS_mosaic_folder
        self.mask_mosaics = mask_folder
        self.fits_folder = fits_folder
        self.gal_mosaic_folder = gal_mosaic_folder
        
        self.page_name = self.VFID+'.html'   #name of galaxy html page
        self.gal_htmlpath = self.local_path+self.page_name
        
        self.w3params_nopsf = w3params_nopsf
        self.w3params_psf = w3params_psf
        self.rparams_nopsf = rparams_nopsf
        self.rparams_psf = rparams_psf
        
        self.RA = self.cutcat['RA_1'][self.galaxy_index]
        self.RA = str(self.RA[0])
        self.DEC = self.cutcat['DEC_1'][self.galaxy_index]
        self.DEC = str(self.DEC[0])
        self.VFID = self.cutcat['VFID'][self.galaxy_index]
        self.VFID = str(self.VFID[0])
        
        #number of instances of central VFID galaxy appearing in the relevant dummycat column --> ncomp-1
        if self.VFID in dummycat['central galaxy']:
            ncomp_flag = (dummycat['central galaxy']==self.VFID)
            self.ncomp = len(dummycat['central galaxy'][ncomp_flag])+1   #the +1 accounts for the central galaxy itself
            self.ext_list = dummycat['ID'][ncomp_flag]
        else:
            self.ncomp = 1  
                
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
        
        #first must be sure that r-band is projected onto w3 coordinates (ensures all images have same size)
        r_scaled, footprint = reproject_interp((self.r_im, self.r_header), wise_hdu.header)
        
        titles = ['W3 Image', 'r-band Image', 'LS Image']
        images = [self.wise_im, r_scaled, self.filename_LS]
        
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
        
        self.pngnames = [self.gal_mosaic_folder+self.galname+'-'+'galfit-model-w3-nopsf.png',
                   self.gal_mosaic_folder+self.galname+'-'+'galfit-model-w3-psf.png',
                   self.gal_mosaic_folder+self.galname+'-'+'galfit-model-r-nopsf.png',
                   self.gal_mosaic_folder+self.galname+'-'+'galfit-model-r-psf.png']
        
        print('For self.create_model_mosaics(index), index=0 is w3_nopsf, 1 is w3_psf, 2 is r_nopsf, 3 is r_psf.')

    def create_model_mosaics(self, percentile1=.5, percentile2=99.5, p1residual=5, p2residual=99, cmap='viridis'):
        
        '''
        ARGS:
        index = 0 for w3_nopsf, 1 for w3_psf, 2 for r_nopsf, 3 for r_psf
        percentile1 = min percentile for stretch of image and model
        percentile2 = max percentile for stretch of image and model
        p1residual = min percentile for stretch of residual
        p2residual = max percentile for stretch of residual
        cmap = colormap, default is viridis
        ''' 
        
        for index in self.index_list:
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
            plt.savefig(self.pngnames[index],dpi=300)
            plt.close()    

    def create_mask_mosaics(self):
        print('Under development, pending completion of masking.')   #may use same scaling technique as with LS_mosaics 
        #w3
        #r-band
        return
    
    def tablulate_parameters(self):
        
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
                single_ncomp_params=[]
                for num in range(ncomp):  #say n=2. then num = 0, 1 (indices of external galaxies in self.ext_list)
                    param_row=params[params['VFID']=self.ext_list[num]]  #find where VFID = external galaxy ID
                    single_ncomp_params.append(param_row)  #append row values to single ncomp list 
                self.page_extparams.append(single_ncomp_params)  #will comprise list of external galaxy lists
                
    #create VFIDxxxx.html for the galaxy!
    def WRITETHEGALPAGE(self):
        
        title_dictionary = {0:'w3 noConv',
                            1:'w3 Conv',
                            2:'r-band noConv',
                            3:'r-band Conv'}
        
        LS_path = self.LS_mosaics+self.VFID+'-LS-mosaic.png'

        with open(self.gal_htmlpath,'w') as html:
            
            html.write('<html><body> \n')
            html.write(f'<title>+{self.VFID}-{self.objname}</title> \n')
            html.write('<style type="text/css"> \n')
            html.write('.img-container{text-align: left;} \n')
            html.write('table, td, th {padding: 5px; text-align: center; border: 1px solid black;} \n')
            html.write('p {display: inline-block;;} \n')
            html.write('</style> \n')
        
            html.write(f'<font size="40">Central Galaxy: {self.VFID}-{self.objname} </font><br /> \n')
            html.write('<a href=main.html>Return to Homepage</a></br /> \n')

            if i != len(self.cutcat)-1:
                html.write('<a href='+str(self.cutcat['VFID'][i+1])+'.html>Next Galaxy</a></br /> \n') 

            if i != 0:
                html.write('<a href='+str(self.cutcat['VFID'][i-1])+'.html>Previous Galaxy</a></br /> \n')

            for index in self.psf_indices:
                mosaic_path = self.pngnames[index]
                params = self.page_params[index]
                html.write(f'<font size="30"> GALFIT Output Mosaic {title_dictionary[index]}:</font><br /> \n')
                html.write('<div class='+'"'+'img-container'+'"> <!-- Block parent element --> <img src='+'"'+mosaic_path+'" height="70%" width="70%" /><br /> \n')
                
                html.write('<table><tr><th>VFID</th><th>Type</th><th>xc</th><th>xc_err</th><th>yc</th><th>yc_err</th><th>mag</th><th>mag_err</th><th>Re</th><th>Re_err</th><th>nser</th><th>nser_err</th><th>BA</th><th>BA_err</th><th>PA</th><th>PA_err</th><th>err_flag</th></tr> \n')

                html.write('<tr><td>'+self.VFID+'</td> \n')
                html.write('<td>Host</td> \n')
                for i in range(1,15):
                    html.write(f'<td>{params[i]}</td> \n')
                html.write(r'<td>{params[17]}</td></tr> \n')

                #add the external galaxy parameters, if any.
                if self.ncomp>1:
                    param_ext=self.page_extparams[index]  #isolate band&psf parameters of the external galaxies
                    for num in range(ncomp):  #say n=2. then num = 0, 1 (indices of external galaxies in self.ext_list)
                        single_param_ext = param_ext[num]
                        html.write('<tr><td>'+self.ext_list[num]+'</td> \n')  #VFID of external galaxy
                        html.write('<td>External</td> \n')
                        for i in range(1,15):
                            html.write(f'<td>{single_param_ext[i]}</td> \n')
                        if num==self.ncomp-1:
                            html.write(r'<td>{params[17]}</td></tr> \n')
                
                #I want a table for every nopsf/psf pair, so the table will finish if index is psf w3 or psf r-band (OR neither, in the case that I don't include any psf bands...). I would attempt to account for the case where I am only including nopsf tables, but I do not anticipate this condition being met. Ever.
                if (index==1)|(index==3)|(index==np.max(self.psf_indices)):
                    html.write('</tr></table> \n')
                
                #add the LS mosaic (and, when completed, the mask mosaic)
                if index == np.max(self.psf_indices):
                    html.write('<div class='+'"'+'img-container'+'"> <!-- Block parent element --> <img src='+'"'LS_path'" height="60%" width="60%" /><br /> \n')
                    #html.write('<div class='+'"'+'img-container'+'"> <!-- Block parent element --> <img src='+'"'+mask_path+'" height="60%" width="50%" /><br /> \n')
                
            html.write('<a href=main.html>Return to Homepage</a></br /> \n')
           
            if i != len(self.cutcat)-1:
                html.write('<a href='+str(self.cutcat['VFID'][i+1])+'.html>Next Galaxy</a></br /> \n') 

            if i != 0:
                html.write('<a href='+str(self.cutcat['VFID'][i-1])+'.html>Previous Galaxy</a></br /> \n')

            html.write('<br /><br />\n')    
            html.write('</html></body>\n')     

            html.close()   


#---tweak param file, be sure the hoopla below is consistent, write a printed comment (similar to re_analysis.py) explaining how to run each component of the classes.
#AND THEN: the testing phase. :-)

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
