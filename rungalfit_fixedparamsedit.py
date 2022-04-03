#!/usr/bin/env python

'''


REQUIRED MODULES:
pyds9


'''
try:
    import pyds9
except ImportError:
    print("WARNING: could not import pyds9")
import os
from astropy.io import fits
from astropy.table import Table
import numpy as np
homedir = os.getenv("HOME")
dummycat = Table.read(homedir+'/dummycat.fits',format='ascii')




def parse_galfit_1comp(galfit_outimage,ncomp,asymflag=0,printflag=False):
    numerical_error_flag=0
    if asymflag:
        header_keywords=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_SKY','1_F1','1_F1PA','CHI2NU']
    else:
        header_keywords=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_SKY','CHI2NU']
    if ncomp == 2:
        header_keywords=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_XC','2_YC','2_MAG','2_RE','2_N','2_AR','2_PA','3_SKY','CHI2NU']
    if ncomp == 3:
        header_keywords=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_XC','2_YC','2_MAG','2_RE','2_N','2_AR','2_PA','3_XC','3_YC','3_MAG','3_RE','3_N','3_AR','3_PA','4_SKY','CHI2NU']

    if ncomp == 4:
        header_keywords=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_XC','2_YC','2_MAG','2_RE','2_N','2_AR','2_PA','3_XC','3_YC','3_MAG','3_RE','3_N','3_AR','3_PA','4_XC','4_YC','4_MAG','4_RE','4_N','4_AR','4_PA','5_SKY','CHI2NU']



        
    fit_parameters=[]
    working_dir=os.getcwd()+'/'
    image_header = fits.getheader(galfit_outimage,2)
    for hkey in header_keywords:
        s=str(image_header[hkey])
        #print hkey
        if s.find('[') > -1:
            s=s.replace('[','')
            s=s.replace(']','')
            t=s.split('+/-')
            values=(float(t[0]),0.)# fit and error
        else:
            t=s.split('+/-')
            try:
                values=(float(t[0]),float(t[1]))# fit and error
            except ValueError:
                # look for * in the string, which indicates numerical problem
                if t[0].find('*') > -1:
                    numerical_error_flag=1
                    t[0]=t[0].replace('*','')
                    t[1]=t[1].replace('*','')
                    values=(float(t[0]),float(t[1]))# fit and error
            except IndexError: # for CHI2NU
                chi2nu=float(t[0])
                continue
        fit_parameters.append(values)
        if printflag:
            print('{:6s}: {:s}'.format(hkey,s))
    fit_parameters.append(numerical_error_flag)
    fit_parameters.append(chi2nu)
    #if printflag:
    #    for i in range(len(header_keywords)):
    #        try:
    #            print('%6s : %5.2f +/- %5.2f'%(header_keywords[i],t[i][0],t[i][1]))
    #        except TypeError:
    #            print(header_keywords[i])

    return fit_parameters










class galfit:
    def __init__(self,galname=None,vfid=None,image=None,sigma_image=None,psf_image=None,psf_oversampling=None,mask_image=None,xminfit=None,yminfit=None,xmaxfit=None,ymaxfit=None,convolution_size=None,magzp=None,pscale=None,convflag=1,constraintflag=1,fitallflag=0,ncomp=1):
        self.galname=galname
        self.image=image
        self.vfid=vfid

        self.sigma_image=sigma_image
        self.psf_image=psf_image
        self.psf_oversampling=psf_oversampling
        self.mask_image=mask_image
        self.xminfit=xminfit
        self.yminfit=yminfit
        self.xmaxfit=xmaxfit
        self.ymaxfit=ymaxfit
        self.convolution_size=convolution_size
        self.magzp=magzp
        self.pscale=pscale
        self.convflag=convflag
        self.constraintflag=constraintflag
        self.fitallflag=fitallflag
        vfid = self.galname[0:8]
        
        if vfid in dummycat['central galaxy']:
            self.ncomp = len(np.where(dummycat['central galaxy'] == vfid)[0]) + 1
        else:
            self.ncomp=ncomp

        self.asymmetry=0

        #print('***%%%%%%%%%%%%%%%%%')
        #print('inside galfit class')
        #print(xminfit,xmaxfit,yminfit,ymaxfit,convolution_size)
        #print(self.xminfit,self.xmaxfit,self.yminfit,self.ymaxfit,self.convolution_size)
        #print('psf_image = ',psf_image)
        #print('self.fitall = ',self.fitallflag)
        #print('***%%%%%%%%%%%%%%%%%')

        


        
    def create_output_names(self):
        if self.asymmetry:
            output_image=str(self.galname)+'-'+ str(self.ncomp) +'Comp-galfit-out-asym.fits'
        elif self.convflag:
            output_image=str(self.galname)+'-'+ str(self.ncomp) +'Comp-galfit-out-conv.fits'
        else:
            output_image=str(self.galname)+'-'+ str(self.ncomp) +'Comp-galfit-out.fits'

        self.output_image=output_image
        # create galfit input file
        self.galfile=str(self.galname)+'-galfit.input.'+str(self.ncomp)+'Comp'


    def open_galfit_input(self):
        self.galfit_input=open(self.galfile,'w')


    def write_image_params(self):#,input_image,output_image,sigma_image,psf_image,psf_oversampling,mask_image,xminfit,xmaxfit,yminfit,ymaxfit,convolution_size,magzp,pscale,convflag=1,constraintflag=1,fitallflag=0):
        self.galfit_input.write('# IMAGE PARAMETERS\n')
        self.galfit_input.write('A) '+self.image+'              # Input data image (FITS file)\n')
        self.galfit_input.write('B) '+self.output_image+'       # Name for the output image\n')
        self.galfit_input.write('C) %s                # Sigma image name (made from data if blank or "none") \n'%(self.sigma_image))
        if self.convflag:
            self.galfit_input.write('D) '+self.psf_image+'     # Input PSF image and (optional) diffusion kernel\n')
            self.galfit_input.write('E) %i                   # PSF oversampling factor relative to data\n'%(self.psf_oversampling))
        #if self.fitallflag:
        #    self.galfit_input.write('F)            # Pixel mask (ASCII file or FITS file with non-0 values)\n')
        #else:
        #    self.galfit_input.write('F) '+self.mask_image+'           # Pixel mask (ASCII file or FITS file with non-0 values)\n')
        #if self.constraintflag:
        #    self.galfit_input.write('G) /Users/rfinn/research/LocalClusters/GalfitAnalysis/sersic.constraint        # Parameter constraint file (ASCII)\n')
        self.galfit_input.write('H) '+str(int(round(self.xminfit)))+' '+str(int(round(self.xmaxfit)))+' '+str(int(round(self.yminfit)))+' '+str(int(round(self.ymaxfit)))+'     # Image region to fit (xmin xmax ymin ymax)\n')
        if self.convflag:
            self.galfit_input.write('I) '+str(int(round(self.convolution_size)))+' '+str(int(round(self.convolution_size)))+'             # Size of convolution box (x y)\n')
        self.galfit_input.write('J) %5.2f              # Magnitude photometric zeropoint \n'%(self.magzp))
        self.galfit_input.write('K) %6.5f   %6.5f         # Plate scale (dx dy)  [arcsec/pix]\n'%(self.pscale,self.pscale))
        self.galfit_input.write('O) regular                # Display type (regular, curses, both)\n')
        self.galfit_input.write('P) 0                   # Create output image only? (1=yes; 0=optimize) \n')
        self.galfit_input.write('S) 0                   # Modify/create objects interactively?\n')


    def set_sersic_params(self,xobj=None,yobj=None,mag=None,rad=None,nsersic=None,BA=None,PA=None,fitmag=1,fitcenter=1,fitrad=1,fitBA=1,fitPA=1,fitn=1,first_time=0):

#run to hold PA, BA fixed

#        sgacut = Table.read(homedir+'/sga_cut.fits',format='ascii')
#        if self.vfid in sgacut['VFID']:
#            sgaindex = np.where(sga_params['VFID'] == self.vfid)[0]

#            BA = sga_params['BA'][sgaindex]
#            PA = sga_params['PA'][sgaindex]
#            fitBA = 0
#            fitPA = 0
            



        self.xobj=xobj
        self.yobj=yobj
        self.mag=mag
        self.rad=rad
#a few lines to ensure that the sersic index does not exceed 6 (which would translate to an unphysical model)
        if nsersic>6:
            self.nsersic = 6
        else:
            self.nsersic=nsersic
        self.BA=BA
        self.PA=PA
        self.fitmag=fitmag
        self.fitn=fitn
        self.fitcenter=fitcenter
        self.fitrad=fitrad
        self.fitBA=fitBA
        self.fitPA=fitPA


        if first_time:
            self.xobj0=xobj
            self.yobj0=yobj
            self.mag0=mag
            self.rad0=rad
            self.nsersic0=nsersic
            self.BA0=BA
            self.PA0=PA
            self.fitmag0=fitmag
            self.fitn0=fitn
            self.fitcenter0=fitcenter
            self.fitrad0=fitrad
            self.fitBA0=fitBA
            self.fitPA0=fitPA
            self.asymmetry0=self.asymmetry

    def reset_sersic_params(self):
        self.xobj=self.xobj0
        self.yobj=self.yobj0
        self.mag=self.mag0
        self.rad=self.rad0
#in case the nsersic can be larger than 6 here, I apply a similar upper limit as above
        if nself.nsersic0>6:
            self.nsersic=6
        else:
            self.nsersic=self.nsersic0
        self.BA=self.BA0
        self.PA=self.PA0
        self.fitmag=self.fitmag0
        self.fitn=self.fitn0
        self.fitcenter=self.fitcenter0
        self.fitrad=self.fitrad0
        self.fitBA=self.fitBA0
        self.fitPA=self.fitPA0
        self.asymmetry=self.asymmetry0
        
    def set_sky(self,sky):
        self.sky=sky

    def write_sersic(self,objnumber,profile):

        
        self.galfit_input.write(' \n')
        self.galfit_input.write('# Object number: %i \n'%(objnumber))
        self.galfit_input.write(' 0) %s             # Object type \n'%(profile))
        self.galfit_input.write(' 1) %8.1f  %8.1f %i %i  # position x, y        [pixel] \n'%(self.xobj,self.yobj,int(self.fitcenter),int(self.fitcenter)))
        self.galfit_input.write(' 3) %5.2f      %i       # total magnitude     \n'%(self.mag,self.fitmag))
        self.galfit_input.write(' 4) %8.2f       %i       #     R_e              [Pixels] \n'%(self.rad,self.fitrad))
        print('sersic n, fitsersicn = ',self.nsersic,self.fitn)
        self.galfit_input.write(' 5) %5.2f       %i       # Sersic exponent (deVauc=4, expdisk=1)   \n'%(self.nsersic,int(self.fitn)))
        print('BA, fitBA = ',self.BA,self.fitBA)
        self.galfit_input.write(' 9) %5.2f       %i       # axis ratio (b/a)    \n'%(self.BA,int(self.fitBA)))
        self.galfit_input.write('10) %5.2f       %i       # position angle (PA)  [Degrees: Up=0, Left=90] \n'%(self.PA,int(self.fitPA)))
        if self.asymmetry:
            self.galfit_input.write('F1) 0.0001 0.00   1  1     # azim. Fourier mode 1, amplitude & phase angle \n')
        self.galfit_input.write(" Z) 0                  # Output option (0 = residual, 1 = Don't subtract)  \n")           

        

    def write_sky(self,objnumber):    
        self.galfit_input.write(' \n')
        self.galfit_input.write('# Object number: %i \n'%(objnumber))
        self.galfit_input.write(' 0) sky             # Object type \n')
        self.galfit_input.write(' 1) %8.1f   1  # sky background at center of fitting region [ADUs] \n'%(self.sky))
        self.galfit_input.write(' 2) 0      0       # dsky/dx (sky gradient in x)    \n')
        self.galfit_input.write(' 3) 0      0       # dsky/dy (sky gradient in y) \n')
        self.galfit_input.write(" Z) 0                  # Output option (0 = residual, 1 = Don't subtract)  \n")




    def add_simple_sersic_object(self,objnumber,profile,x,y,mag,rad,nsersic,BA,PA):
        self.galfit_input.write(' \n')
        self.galfit_input.write('# Object number: %i \n'%(objnumber))
        self.galfit_input.write(' 0) %s             # Object type \n'%(profile))
        self.galfit_input.write(' 1) %8.1f  %8.1f 1 1  # position x, y        [pixel] \n'%(x,y))
        self.galfit_input.write(' 3) %5.2f      1       # total magnitude     \n'%(mag))
        self.galfit_input.write(' 4) %8.2f       1       #     R_e              [Pixels] \n'%(rad))
        self.galfit_input.write(' 5) %5.2f       1       # Sersic exponent (deVauc=4, expdisk=1)   \n'%(nsersic))
        self.galfit_input.write(' 9) %5.2f       1       # axis ratio (b/a)    \n'%(BA))
        self.galfit_input.write('10) %5.2f       1       # position angle (PA)  [Degrees: Up=0, Left=90] \n'%(PA))
        self.galfit_input.write(" Z) 0                  # Output option (0 = residual, 1 = Don't subtract)  \n")

    
    def run_galfit(self,displayflag=False):

        from astropy.table import Table
        import os
        import numpy as np
        
        homedir = os.getenv("HOME")
        dummycat = Table.read(homedir+'/dummycat.fits',format='ascii')
        #print 'self.fitall = ',self.fitall
        self.create_output_names()
        self.open_galfit_input()
        #print('in rungalfit.run_galfit, self.psf_image = ',self.psf_image)
        self.write_image_params()
        #print 'self.fitall = ',self.fitall
        self.write_sersic(1,'sersic')
        print(self.galname[0:8])
        vfid = self.galname[0:8]
        if vfid in dummycat['central galaxy']:
            indices = np.where(dummycat['central galaxy'] == vfid)[0]
            print(indices)
            for i in range(0,len(indices)):
                index = indices[i]
                n = int(i)+2
                #1 already taken, i begins at 0
                x = int(dummycat['x (px)'][index])
                y = int(dummycat['y (px)'][index])
                print(x,y)
                self.add_simple_sersic_object(n,'sersic',x=x,y=y,mag=self.mag,rad=self.rad,nsersic=self.nsersic,BA=self.BA,PA=self.PA)

                
                if n == int(len(indices)+1):
                    self.write_sky(n+1)

        #print 'self.fitall = ',self.fitall
        else:
            self.write_sky(2)
        #print 'self.fitall = ',self.fitall
        if (self.fitallflag):
            print('%%%%%%%%%%%%%% HEY %%%%%%%%%%%%%')
            print('I think fitall is true, just sayin...')
            self.fitall()
        self.close_input_file()
        #print 'self.fitall = ',self.fitall
        s = 'galfit '+self.galfile
        print('run the following: ',s)

        #try:
        #    errno=os.system(s)
        #    self.galfit_flag=1
        #except:
        #    print "PROBLEM RUNNING GALFIT!!!!"
        #    self.galfit_flag=0
        #    return

        errno=os.system(s)
        self.galfit_flag=1

        image_id=str(self.galname)+'-'
        self.galfit_log=image_id+str(self.ncomp)+'Comp-fit.log'
        s='cp fit.log '+self.galfit_log
        os.system(s)
        self.galfit_out=image_id+str(self.ncomp)+'Comp'+'-galfit.01'
        s='mv galfit.01 '+self.galfit_out
        try:
            os.rename('galfit.01',self.galfit_out)
        except:
            print("I hate to be the bearer of bad news, but it appears like galfit did not complete.")
            #galflag[j]=0
            self.galfit_flag=0
            return
        if displayflag:
            self.display_results()

    def display_results(self):
        print('%%%%%%%%%%%%%%%%%%')
        print('inside display_results')
        print('self.galfit_flag = ',self.galfit_flag)
        if (self.galfit_flag < 0.1):
            print('GALFIT did not complete - can not display results')
            return

        subcomp_image=self.galname+'-'+str(self.ncomp)+'Comp'+'-subcomps.fits'
        s='galfit -o3 '+self.galfit_out
        os.system(s)
        os.rename('subcomps.fits',subcomp_image)
        #    - display results (like ds9 multiextension data cube -> use xpa)
        #

        try:
            d.set('frame delete all')
        except NameError:
            d=pyds9.DS9()
            d.set('frame delete all')
        #print 'file to display = ',self.output_image
        s='file new multiframe '+self.output_image
        #print s

        d.set(s)
        d.set('frame delete 1')
        for k in range(2,5):
            s='frame '+str(k)
            d.set(s)
            d.set('zoom to fit')
            #print k
            if k == 2:
                d.set('regions command {text 30 10 #text="Image" font="times 18 bold" color="red"}')
            if k == 3:
                d.set('regions command {text 30 10 #text="Model" font="times 18 bold" color="red"}')
            if k == 4:
                d.set('regions command {text 30 10 #text="Residual" font="times 18 bold" color="red"}')
        d.set('frame match wcs')

        galfit_outimage=self.galname+'-'+'1Comp-galfit-out.fits[2]'
        galfit_outimage_header=self.galname+'-'+'1Comp-galfit-out.fits'
        self.print_galfit_results(galfit_outimage_header)
        print('file to display = ',subcomp_image)
        s='file new multiframe '+subcomp_image
        d.set(s)

        if self.ncomp == 1:
            endframe=8
        if self.ncomp == 2:
            endframe=9
        if self.ncomp == 3:
            endframe=10
        s='frame delete '+str(endframe)
        print(s)
        try:
            d.set(s)
        except ValueError:
            print("couldn't execute the following ds9 command : ",s)
        try:
            d.set('frame delete 5')
        except:
            print("couldn't delete frame 5")
        try:
            d.set('frame delete 6')
        except:
            print("couldn't delete frame 6")
        d.set('frame 7')
        if self.mask_image == None:
            print('no mask available')
        else:
            try:
                d.set('file '+self.mask_image)
            except:
                print("couldn't load mask file ",self.mask_image)

        for k in range(2,endframe):
            if k == 5:
                continue
            if k == 6:
                continue
            s='frame '+str(k)
            d.set(s)
            d.set('scale log')
            d.set('zoom to fit')
        d.set('frame match wcs')
        #d.set('saveimage png image.png')
        #img_name=self.galname+'-'+str(self.ncomp)+'Comp.png'
        #os.rename('image.png',img_name)


    def fitall(self,mindistance=8):
        os.system('cp '+homedir+'research/LocalClusters/sextractor/default.param .')
        os.system('cp '+homedir+'research/LocalClusters/sextractor/default.nnw .')
        s='sex '+self.image+'[1] -c '+homedir+'research/LocalClusters/sextractor/default.sex.24um.galfitsource -WEIGHT_TYPE MAP_RMS -WEIGHT_IMAGE '+self.sigma_image+' -CATALOG_NAME '+self.galname+'test.cat -CATALOG_TYPE ASCII_HEAD'
        os.system(s)
        # read in SE table to get x,y for sources
        #fname=self.galname+'test.fits'
        fname=self.galname+'test.cat'
        print('FITALL CATALOG NAME = ',fname)
        objnumber=2
        profile='sersic'
        try:
            se=atpy.Table(fname,type='ascii')
            print('found ',len(se.X_IMAGE),' sources on the field of ',self.galname)
            nearbyobjflag=sqrt((se.X_IMAGE-self.xobj)**2+(se.Y_IMAGE-self.yobj)**2) > mindistance
            for k in range(len(se.X_IMAGE)):
                if nearbyobjflag[k]:
                    objnumer=objnumber+1
                    self.add_simple_sersic_object(objnumber,profile,se.X_IMAGE[k],se.Y_IMAGE[k],se.MAG_BEST[k],se.FLUX_RADIUS[k,0],2,se.B_IMAGE[k]/se.A_IMAGE[k],se.THETA_IMAGE[k])
        except AttributeError:
            print('WARNING: no sources detected in image!')
        raw_input=('hit any key to continue \n')

    def add_simple_sersic_object(self,objnumber,profile,x,y,mag,rad,nsersic,BA,PA):
        self.galfit_input.write(' \n')
        self.galfit_input.write('# Object number: %i \n'%(objnumber))
        self.galfit_input.write(' 0) %s             # Object type \n'%(profile))
        self.galfit_input.write(' 1) %8.1f  %8.1f 1 1  # position x, y        [pixel] \n'%(x,y))
        self.galfit_input.write(' 3) %5.2f      1       # total magnitude     \n'%(mag))
        self.galfit_input.write(' 4) %8.2f       1       #     R_e              [Pixels] \n'%(rad))
        self.galfit_input.write(' 5) %5.2f       1       # Sersic exponent (deVauc=4, expdisk=1)   \n'%(nsersic))
        self.galfit_input.write(' 9) %5.2f       1       # axis ratio (b/a)    \n'%(BA))
        self.galfit_input.write('10) %5.2f       1       # position angle (PA)  [Degrees: Up=0, Left=90] \n'%(PA))
        self.galfit_input.write(" Z) 0                  # Output option (0 = residual, 1 = Don't subtract)  \n")
                                

    def close_input_file(self):
        self.galfit_input.close()

    def print_params(self):
        print('CURRENT INPUTS: \n mag = %5.2f %i \n Re = %5.2f %i \n n = %5.2f %i\n B/A = %5.2f %i \n PA = %5.2f %i \n fitall = %i \n fitcenter = %i \n'%(self.mag,self.fitmag,self.rad,self.fitrad,self.nsersic,self.fitn,self.BA,self.fitBA,self.PA,self.fitPA,self.fitallflag,self.fitcenter))
                                
    def print_galfit_results(self,image):            
        if self.asymmetry:
            header_keywords=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_SKY','1_F1','1_F1PA','ERROR','CHI2NU']
        else:
            header_keywords=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_SKY','ERROR','CHI2NU']

        if self.ncomp == 2:
            header_keywords=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_XC','2_YC','2_MAG','2_RE','2_N','2_AR','2_PA','3_SKY','ERROR','CHI2NU']

        if self.ncomp == 3:
            header_keywords = ['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_XC','2_YC','2_MAG','2_RE','2_N','2_AR','2_PA','3_XC','3_YC','3_MAG','3_RE','3_N','3_AR','3_PA','4_SKY','ERROR','CHI2NU']
       
        if self.ncomp == 4:
            header_keywords = ['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_XC','2_YC','2_MAG','2_RE','2_N','2_AR','2_PA','3_XC','3_YC','3_MAG','3_RE','3_N','3_AR','3_PA','4_XC','4_YC','4_MAG','4_RE','4_N','4_AR','4_PA','5_SKY','ERROR','CHI2NU']    



            
        t=parse_galfit_1comp(image,ncomp=self.ncomp,printflag=True)
    def edit_params_menu(self):
        flag=str(raw_input('What is wrong?\n n = adjust sersic \n r = reset Re \n o = nearby object (toggle fitall) \n b = B/A \n p = PA \n m = mag \n c = recenter \n f = hold values fixed \n a = toggle asymmetry parameter \n R = reset to original values \n g = go (run galfit) \n x=quit \n '))
        return flag

    def set_n(self):
        n=float(raw_input('sersic exponent = '))
        self.nsersic=n

    def set_r(self):
        r=float(raw_input('radius = '))
        self.rad=r

    def set_BA(self):
        r=float(raw_input('BA = '))
        self.BA=r

    def set_PA(self):
        r=float(raw_input('PA = '))
        self.PA=r

    def set_mag(self):
        r=float(raw_input('mag = '))
        self.mag=r

    def set_center(self):
        r=float(raw_input('xc = '))
        self.xobj=r
        r=float(raw_input('yc = '))
        self.yobj=r

    def toggle_fitall(self):
        self.fitallflag=toggle(self.fitallflag)

    def toggle_asymmetry(self):
        self.asymmetry=toggle(self.asymmetry)

    def print_fix_menu(self):
        self.print_params()
        flag3=str(raw_input('What do you want to hold fixed/toggle?\n n = fix sersic index \n r = fix Re \n b = fix B/A \n p = PA \n c = center \n f = use constraint file \n R = reset to original values \n g = go (run galfit) \n x=quit \n '))
        return flag3

    def fix_n(self):
        n=float(raw_input('sersic exponent = '))
        self.nsersic=n
        self.fitn=toggle(self.fitn)

    def fix_rad(self):
        self.fitrad=toggle(self.fitrad)

    def fix_BA(self):
        self.fitBA=toggle(self.fitBA)
        print(self.fitBA, self.BA)
    def fix_PA(self):
        self.fitPA=toggle(self.fitPA)

    def fix_center(self):
        self.fitcenter=toggle(self.fitcenter)

    def add_constraint_file(self):
        self.constraintflag=toggle(self.constraintflag)
        

if __name__ == '__main__':
    from astropy.table import Table
    from astropy.io import ascii
    homedir = os.getenv("HOME")
    dummycat = Table.read(homedir+'/dummycat.fits',format='ascii')
    vfcut = Table.read(homedir+'/vfcut.fits',format='ascii')
    sgacut = Table.read(homedir+'/sga_cut.fits',format='ascii')
    sga_params = Table.read(homedir+'/sga_vf_matched.fits')
    vf_largeser = Table.read(homedir+'/vf_largeser.fits')
