'''
Sensitivity Study: how many pixels of a galaxy 2D cutout can I block/mask before the data is no longer usable?

AIMS:
--> for a single galaxy postage stamp whose background requires no masking in w3, generate random locations in the image and mask some group of pixels similar to the number needed for a small galaxy or star. 
--> the first mask is an NxN matrix of 0s, and each consecutive iteration adds some circular aperture of 1s until a user-defined percentage of pixels are masked (currently using 70%).
--> I have now added the option of omitting the central galaxy from being masked.
--> additional functions will set up and run galfit on the given galaxy with each mask, recording the output parameters as well as the percentage of masked pixels. 
--> finally, yet another function will output an aesthetically pleasing figure illustrating the dependence of some given parameter on the mask-to-all pixel ratio, having run the above algorithm ten times. 


--currently have two directories located at homedir/ --> NGC3364 and UGC06711
--pathnames are currently tailored to me, myself, and my purposes. because I can. 
--to run...open terminal, 
>ipython
>%run ~/mask_iter.py
then proceed. I can automate the plot-generation as well. 
'''

import os
homedir = os.getenv("HOME")

from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import fits
from astropy.io import ascii
import glob
from astropy.wcs import WCS
import glob
from reproject import reproject_interp

import random


def define_paths(galaxy_name,homedir=homedir):
    path_to_dir = homedir+'/Desktop/maskies_'+galaxy_name+'/'
    cutout_image_path = path_to_dir+galaxy_name+'-custom-image-W3.fits.fz'
    return path_to_dir,cutout_image_path


def set_cutout_var(cutout_image_path):
    cutout_im = fits.getdata(cutout_image_path)
    return cutout_im

    
def show_cutouts(ngc3364_cutout, ugc06711_cutout, homedir=homedir, save=False):
    fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(10,10))
    ax1.imshow(ngc3364_cutout,origin='lower')
    ax2.imshow(ugc06711_cutout,origin='lower')
    ax1.set_title('NGC3364',fontsize=15)
    ax2.set_title('UGC06711',fontsize=15)
    plt.show()
    if save:
        plt.savefig(homedir+'/Desktop/galaxy_cutouts.png',dpi=250)
    plt.close()

    
#convert cutouts from .fz to .fits
#DO FOR BOTH PSF AND IM CUTOUTS
def fz_to_fits(path_to_im, galaxy_name, psf_or_im):
    if psf_or_im=='psf':
        galaxy = path_to_im+galaxy_name+'-custom-psf-W3.fits.fz'
        fz = fits.getdata(galaxy)
        fz_header = (fits.open(galaxy)[0]).header
    if psf_or_im=='im':
        galaxy = path_to_im+galaxy_name+'-custom-image-W3.fits.fz'
        fz = fits.getdata(galaxy)
        fz_header = (fits.open(galaxy)[1]).header
    fits.writeto(galaxy[:-3],fz,header=fz_header,overwrite=True)    #[:-3] removes the .fz from location string
    
    
def set_parameters(path_to_dir,cutout_im,galaxy_name,psf=False):
    
    galname = galaxy_name
    image = path_to_dir+galaxy_name+'-custom-image-W3.fits'
    sigma_image = glob.glob(path_to_dir+'*std*')[0]
    if psf==True:
        psf_image = path_to_dir+galaxy_name+'-custom-psf-W3.fits'
    xminfit = 0
    xmaxfit = np.shape(cutout_im)[0]
    yminfit = 0
    ymaxfit = np.shape(cutout_im)[1]
    xobj = int(xmaxfit/2)
    yobj = int(ymaxfit/2)
    
    #programming syntax gods, please forgive me.
    if psf==True:
        return galname,image,sigma_image,psf_image,xminfit,xmaxfit,yminfit,ymaxfit,xobj,yobj
    if psf==False:
        return galname,image,sigma_image,xminfit,xmaxfit,yminfit,ymaxfit,xobj,yobj
    
        
#adapted from https://stackoverflow.com/questions/44865023/how-can-i-create-a-circular-mask-for-a-numpy-array
def mock_mask(h, w, center=None, radius=None):

    if center is None: # use the middle of the image
        center = (int(w/2), int(h/2))
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]   #create x, y coordinates for every pixel (that otherwise have just an index number)
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = dist_from_center >= radius
    return mask
                        
                        
 #returns number of pixels masked relative to number of pixels in the cutout
def percent_masked(mask,invert=True):
    mask = ~mask if invert else mask
    return(np.sum(mask)/(np.shape(mask)[0]**2))


def isolate_central_galaxy(path_to_dir,galaxy_name,band='w3'):
    
    moustakas_maskname = glob.glob(path_to_dir+galaxy_name+'*maskbits*')[0]
    moustakas_mask = fits.getdata(moustakas_maskname)
    
    if band=='w3':
        
        w3_reference = path_to_dir+galaxy_name+'-custom-image-W3.fits.fz'
        
        rband_hdu = fits.open(moustakas_maskname)[1]
        wise_hdu = fits.open(w3_reference)[1]
        
        wise_mask, footprint = reproject_interp(rband_hdu, wise_hdu.header)
        
        moustakas_mask = wise_mask
    
    '''
    #isolate pixels that DO NOT BELONG to central galaxy, assuming it is in the SGA

    #(First bool --> selects saturated pixels that overlap with the galaxy)
    #(Second bool --> selects pixels that are masked (value > 0) AND have a value smaller than the galaxy flag)

    test = (np.array(rband_hdu.data) > 4096) | ((np.array(rband_hdu.data) > 0) & (np.array(rband_hdu.data) < 4096))

    #result --> test now has a values of one (masked) and zero (unmasked)
    
    #For central_galaxy, simply reverse the > and <.
    '''
    central_galaxy_mask = (np.array(moustakas_mask) < 4096) | ((np.array(moustakas_mask) < 0) & (np.array(moustakas_mask) > 4096)) 
    
    return central_galaxy_mask
  

#if mask_gal=False, then none of the masked pixels will intersect or lie within the bounds of the central galaxy
def create_masks(path_to_dir, cutout_image_path, circ_diameter, central_galaxy, mask_gal=True, ratio_to_terminate=0.70):
    
    if mask_gal==False:
        central_galaxy_mask = isolate_central_galaxy(path_to_dir,central_galaxy)
    
    cutout_image = fits.getdata(cutout_image_path)
    im_header = (fits.open(cutout_image_path)[1]).header
    
    #mock masks will consist of circles of radius circ_diameter/2
    h = np.shape(cutout_image)[0]
    w = np.shape(cutout_image)[1]
    radius = int(circ_diameter/2)
    
    #define initial ratio of masked-to-all pixels
    ratio = 0
    
    #create list into which updated ratios will lie. 0 will be the first ratio.
    ratios = [0]
    
    #create empty list which will hold mask names
    masknames = ['mask-0.fits']
    
    iteration_number = 0
    
    #so long as the ratio is less than the user's 'termination ratio', generate masks
    while ratio < ratio_to_terminate:        
        
        #center of circ_mask will be random each time.
        center = (random.randrange(h-1),random.randrange(w-1))
        
        #create initial mask for ratio=0, append updated ratio, save mask to path_to_dir
        if ratio == 0:
            
            new_mask = mock_mask(h,w,center,radius)
            
            #if we want no galaxy pixels to be masked, then mask_gal==False
            if mask_gal==False:
                mask_to_save = (~new_mask*central_galaxy_mask).astype(int)
                ratio = percent_masked(~new_mask*central_galaxy_mask,invert=False)
            
            else:
                #the ~new_mask ensures that masked px are 1s and unmasked px are 0s
                mask_to_save = (~new_mask).astype(int)
                ratio = percent_masked(new_mask)
            
            #save mask to appropriate directory;
            #and I multiply ratio by 1000 to include the first decimal place in the image name.
            percent = int(ratio*1000)
                        
            fits.writeto(path_to_dir+'mask-{}.fits'.format(str(percent)), mask_to_save, header=im_header, overwrite=True)
            
            ratios.append(ratio)
            masknames.append('mask-{}.fits'.format(str(percent)))
            
            iteration_number+=1
        
        #create mask, append updated ratio, save mask to path_to_dir
        else:
            
            #old mask
            mask = new_mask.copy()
            
            #new_mask is previous mask multiplied by the new circ_mask
            new_mask = mock_mask(h,w,center,radius)*mask
            
            #if we want no galaxy pixels to be masked, then mask_gal==False
            if mask_gal==False:
                mask_to_save = (~new_mask*central_galaxy_mask).astype(int)
                ratio = percent_masked(~new_mask*central_galaxy_mask,invert=False)
    
            else:
                #the ~new_mask ensures that masked px are 1s and unmasked px are 0s
                mask_to_save = (~new_mask).astype(int)
                ratio = percent_masked(new_mask)
                        
            #save mask to appropriate directory
            percent = int(ratio*1000)
            
            #quality control; ensures new ratio is larger than the previous ratio.
            if ratio > ratios[iteration_number]:
                
                ratios.append(ratio)
                                
                fits.writeto(path_to_dir+'mask-{}.fits'.format(str(percent)),
                             mask_to_save,header=im_header,overwrite=True)
                
                masknames.append('mask-{}.fits'.format(str(percent)))

                iteration_number+=1         
    
    print('RATIO OF RATIO LENGTH TO MASKNAME LENGTH: ',len(ratios)/len(masknames))
    return np.asarray(ratios),masknames

    
def write_image_params(galfit_input, i, image, output_image, sigma_image, xminfit, xmaxfit, yminfit, ymaxfit, magzp=22.5, pscale=0.344, convflag=False, psf_image=None, mask_image=None): 

    psf_oversampling=8
    convolution_size = xmaxfit - xminfit
    
    galfit_input.write('# IMAGE PARAMETERS\n')
    galfit_input.write('A) '+image+'              # Input data image (FITS file)\n')
    galfit_input.write('B) '+output_image+'       # Name for the output image\n')
    galfit_input.write('C) %s                  # Sigma image name (made from data if blank or "none") \n'%(sigma_image))
    
    if convflag:
        galfit_input.write('D) '+psf_image+'     # Input PSF image and (optional) diffusion kernel\n')
        galfit_input.write('E) %i                   # PSF oversampling factor relative to data\n'%(psf_oversampling))
    
    if i != 0:
        galfit_input.write('F) %s                # Bad Pixel Mask (i.e., obj mask; blank or "none") \n'%(mask_image))
    else:
        galfit_input.write('F)                   # Bad Pixel Mask (i.e., obj mask; blank or "none") \n')
    
    galfit_input.write('H) '+str(int(round(xminfit)))+' '+str(int(round(xmaxfit)))+' '+str(int(round(yminfit)))+' '+str(int(round(ymaxfit)))+'     # Image region to fit (xmin xmax ymin ymax)\n')
    
    if convflag:
            galfit_input.write('I) '+str(int(round(convolution_size)))+' '+str(int(round(convolution_size)))+'             # Size of convolution box (x y)\n')
    
    galfit_input.write('J) %5.2f              # Magnitude photometric zeropoint \n'%(magzp))
    galfit_input.write('K) %6.5f   %6.5f         # Plate scale (dx dy)  [arcsec/pix]\n'%(pscale,pscale))
    galfit_input.write('O) regular                # Display type (regular, curses, both)\n')
    galfit_input.write('P) 0                   # Create output image only? (1=yes; 0=optimize) \n')
    galfit_input.write('S) 0                   # Modify/create objects interactively?\n')

    
def write_sersic(galfit_input, objnumber, profile, xobj, yobj, fitcenter=1, mag=7, fitmag=1, rad=5, fitrad=1, nsersic=2, fitn=1, BA=1, fitBA=1, PA=0, fitPA=1):

    galfit_input.write(' \n')
    galfit_input.write('# Object number: %i \n'%(objnumber))
    galfit_input.write(' 0) %s             # Object type \n'%(profile))
    galfit_input.write(' 1) %8.1f  %8.1f %i %i  # position x, y        [pixel] \n'%(xobj,yobj,int(fitcenter),int(fitcenter)))
    galfit_input.write(' 3) %5.2f      %i       # total magnitude     \n'%(mag,fitmag))
    galfit_input.write(' 4) %8.2f       %i       #     R_e              [Pixels] \n'%(rad,fitrad))
    galfit_input.write(' 5) %5.2f       %i       # Sersic exponent (deVauc=4, expdisk=1)   \n'%(nsersic,int(fitn)))
    galfit_input.write(' 9) %5.2f       %i       # axis ratio (b/a)    \n'%(BA,int(fitBA)))
    galfit_input.write('10) %5.2f       %i       # position angle (PA)  [Degrees: Up=0, Left=90] \n'%(PA,int(fitPA)))
    galfit_input.write(" Z) 0                  # Output option (0 = residual, 1 = Don't subtract)  \n")

    
def write_sky(galfit_input, objnumber, sky=0):
    galfit_input.write(' \n')
    galfit_input.write('# Object number: %i \n'%(objnumber))
    galfit_input.write(' 0) sky            # Object type \n')
    galfit_input.write(' 1) %8.1f  1       # sky background at center of fitting region [ADUs] \n'%(sky))
    galfit_input.write(' 2) 0      0       # dsky/dx (sky gradient in x)    \n')
    galfit_input.write(' 3) 0      0       # dsky/dy (sky gradient in y) \n')
    galfit_input.write(" Z) 0              # Output option (0 = residual, 1 = Don't subtract)  \n")

    
def close_input_file(galfit_input):
    galfit_input.close()
                        

def parse_galfit_1comp(outimage,printflag=False):
        
        temp=[]
        numerical_error_flag=0

        header_keywords=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_SKY','CHI2NU']

        fit_parameters=[]
        
        working_dir=os.getcwd()+'/'
        image_header = fits.getheader(outimage,2)
        
        for hkey in header_keywords:
            s=str(image_header[hkey])
            
            if s.find('[') > -1:
                s=s.replace('[','')
                s=s.replace(']','')
                t=s.split('+/-')
                values=(float(t[0]),0.)# fit and error
            
            else:
                t=s.split('+/-')
                
                try:
                    values=(float(t[0]),float(t[1]))# fit and error
                    temp.append(values[0])
                    temp.append(values[1])
                
                except ValueError:
                    # look for * in the string, which indicates numerical problem
                    if t[0].find('*') > -1:
                        numerical_error_flag=1
                        t[0]=t[0].replace('*','')
                        t[1]=t[1].replace('*','')
                        values=(float(t[0]),float(t[1]))# fit and error
                        temp.append(values[0])
                        temp.append(values[1])
                
                except IndexError: # for CHI2NU
                    chi2nu=float(t[0])
                    continue
            
            if printflag:
                print('{:6s}: {:s}'.format(hkey,s))
        
        temp.append(numerical_error_flag)
        temp.append(chi2nu)
        
        fit_parameters.append(temp)
        return fit_parameters                   
    
    
def runitallx3(galaxy_name,convflag=False,mask_gal=True,termination_ratio=0.70):
    
    #define paths
    path_to_dir, cutout_image_path = define_paths(galaxy_name)
    
    #set a variable to the read-in image data
    cutout_im = set_cutout_var(cutout_image_path)
    
    #set galfit initial parameters (most of which are simply the default)
    if convflag==False:
        galname, image, sigma_image, xminfit, xmaxfit, yminfit, ymaxfit, xobj, yobj= set_parameters(path_to_dir, cutout_im, galaxy_name, psf=False)   
    
    if convflag==True:
        galname, image, sigma_image, psf_image, xminfit, xmaxfit, yminfit, ymaxfit, xobj, yobj= set_parameters(path_to_dir, cutout_im, galaxy_name, psf=True)   
        fz_to_fits(path_to_dir,galaxy_name,psf_or_im='psf')
        
    #prepare image files in path_to_dir; must run twice for Moustakas' psf and galaxy images
    fz_to_fits(path_to_dir,galaxy_name,psf_or_im='im')
        
    #create a numpy array of ten empty lists, each will hold the ratios for one iteration
    ratio_arrays = np.zeros(10,dtype=list)
    
    #create labels for output table
    labels = ['xc','xc_err','yc','yc_err','mag','mag_err','re','re_err',
                          'nser','nser_err','BA','BA_err','PA','PA_err','sky','sky_err',
                          'error','chi2nu','iteration_number']
    
    
    #initialize the output table...with the above labels.
    output_table = Table(names=labels)

    #loop through iteration ten times
    for iteration in range(10):

        ratio_arrays[iteration],masknames = create_masks(path_to_dir, cutout_image_path, circ_diameter=15, central_galaxy=galaxy_name, mask_gal=mask_gal, ratio_to_terminate=termination_ratio)
        
        #replace empty list with the ratios!
        ratios = ratio_arrays[iteration]
        
        #for every mask, run galfit once
        for i in range(len(masknames)):

            output_image = path_to_dir+galaxy_name+'-{}-out.fits'.format(int(ratios[i]*1000))
            galfile = path_to_dir+galaxy_name+'-mask{}.1Comp'.format(int(ratios[i]*1000))
            mask_image = path_to_dir+masknames[i]
            
            galfit_input=open(galfile,'w')
            
            if convflag==True:
                write_image_params(galfit_input, i, image, output_image, sigma_image, xminfit, xmaxfit, yminfit, ymaxfit, psf_image=psf_image, mask_image=mask_image, convflag=True)
            
            if convflag==False:
                write_image_params(galfit_input, i, image, output_image, sigma_image, xminfit, xmaxfit, yminfit, ymaxfit, mask_image=mask_image, convflag=False)
            
            write_sersic(galfit_input, 1,'sersic',xobj,yobj)
            write_sky(galfit_input, 2)
            
            close_input_file(galfit_input)

        for i in range(len(masknames)):

            os.chdir(path_to_dir)
            os.system('galfit '+galaxy_name+'-mask{}.1Comp'.format(int(ratios[i]*1000)))

        for i in range(len(masknames)):
            out=path_to_dir+galaxy_name+'-{}-out.fits'.format(int(ratios[i]*1000))
            fit_params=parse_galfit_1comp(out,printflag=False)

            fit_params[0].append(iteration)

            output_table.add_row(fit_params[0])

    #combine ratios into one tidy list
    ratio_arrays_collapsed = np.concatenate(ratio_arrays.copy()).ravel()
    output_table.add_column(ratio_arrays_collapsed,name='ratios')

    #write output table
    if (convflag==False)&(mask_gal==True):
        output_table_name = 'output_parameters_nopsf_unprotected.fits'
    
    if (convflag==False)&(mask_gal==False):
        output_table_name = 'output_parameters_nopsf_protected.fits'
    
    if (convflag==True)&(mask_gal==True):
        output_table_name = 'output_parameters_psf_unprotected.fits'
    
    if (convflag==True)&(mask_gal==False):
        output_table_name = 'output_parameters_psf_protected.fits'
    
    output_table.write(path_to_dir+output_table_name,overwrite=True)
    

def read_data_tab(path_to_dir,convflag,mask_gal,ydata_type):
        
    data_type_names = ['re','nser','sky']
    
    if (convflag==False)&(mask_gal==True):
        tab_name = 'output_parameters_nopsf_unprotected.fits'
        png_name = 'paramfig_{}_nopsf_unprotected.png'.format(data_type_names[ydata_type])
    
    if (convflag==False)&(mask_gal==False):
        tab_name = 'output_parameters_nopsf_protected.fits'
        png_name = 'paramfig_{}_nopsf_protected.png'.format(data_type_names[ydata_type])
    
    if (convflag==True)&(mask_gal==True):
        tab_name = 'output_parameters_psf_unprotected.fits'
        png_name = 'paramfig_{}_psf_unprotected.png'.format(data_type_names[ydata_type])
    
    if (convflag==True)&(mask_gal==False):
        tab_name = 'output_parameters_psf_protected.fits'
        png_name = 'paramfig_{}_psf_protected.png'.format(data_type_names[ydata_type])

    tab = Table.read(path_to_dir+tab_name)
    print('Using data table '+tab_name)
    return tab,png_name
    
    
#0 for Re, 1 for nser, 2 for sky
def create_maskiter_plots(path_to_dir, convflag, mask_gal, max_iter=10, ydata_type=0, ymin=None, 
                          ymax=None, ylogscale=False, savefig=False):
    
    data_types = [0,1,2]
    data_type_names = ['re','nser','sky']
    
    data_tab,png_name = read_data_tab(path_to_dir,convflag=convflag, mask_gal=mask_gal,ydata_type=ydata_type)
    
    iteration_col = data_tab['iteration_number']
    
    ydata_all = [data_tab['re'],data_tab['nser'],data_tab['sky']]
    yerr_all = [data_tab['re_err'],data_tab['nser_err'],data_tab['sky_err']]
    ylabels_all = ['Effective Radius (px)','Sersic Index','Sky Value']
    
    ydata = ydata_all[ydata_type]
    yerr = yerr_all[ydata_type]
    ylabel = ylabels_all[ydata_type]

    ratios_all = data_tab['ratios']
    max_ratio = np.max(ratios_all)
    
    poss_iter = np.arange(0,max_iter,1)

    plt.figure(figsize=(10,8))

    for iteration in poss_iter:
        flag = (iteration_col==iteration)
        plt.scatter(ratios_all[flag],ydata[flag],alpha=0.2,s=30)
        plt.plot(ratios_all[flag],ydata[flag],alpha=0.1)
    
    plt.xlabel('Masked-to-All Pixel Ratio',fontsize=20)
    plt.ylabel(ylabel,fontsize=20)
    
    if max_iter==10:
        plt.axhline(np.mean(ydata),color='black',alpha=0.6,
                label='Average = {}'.format(np.round(np.mean(ydata),2)))
    
    else:
        plt.axhline(np.mean(ydata[iteration_col<max_iter]),color='black',alpha=0.6,
                label='Average = {}'.format(np.round(np.mean(ydata[iteration_col<max_iter]),2)))
    
    plt.grid(alpha=0.1)

    #bins --> 0-3, 3-5, 5-7
    #avg err in each bin

    #combine ratios into one tidy list
    #ratio_arrays_collapsed = np.concatenate(ratio_arrays.copy()).ravel()
    
    bin_one_err = yerr[(ratios_all>=0) & (ratios_all<=0.3)]
    bin_one_mean = ydata[(ratios_all>=0) & (ratios_all<=0.3)]
    err_one = np.mean(bin_one_err)
    mean_one = np.mean(bin_one_mean)
    std_one = np.std(bin_one_mean)

    bin_two_err = yerr[(ratios_all>0.3) & (ratios_all<=0.5)]
    bin_two_mean = ydata[(ratios_all>0.3) & (ratios_all<=0.5)]
    err_two = np.mean(bin_two_err)
    mean_two = np.mean(bin_two_mean)
    std_two = np.std(bin_two_mean)
    
    bin_three_err = yerr[(ratios_all>0.5) & (ratios_all<=max_ratio)]
    bin_three_mean = ydata[(ratios_all>0.5) & (ratios_all<=max_ratio)]
    if max_iter!=10:
        bin_three_err = yerr[(ratios_all>0.5) & (ratios_all<=0.7) & (iteration_col<max_iter)]
        bin_three_mean = ydata[(ratios_all>0.5) & (ratios_all<=0.7) & (iteration_col<max_iter)]
        
    err_three = np.mean(bin_three_err)
    mean_three = np.mean(bin_three_mean)
    std_three = np.std(bin_three_mean)
    
    plt.axvline(0,color='r',alpha=0.2)
    
    plt.axvline(0.3,color='r',alpha=0.2)

    plt.axvline(0.5,color='r',alpha=0.2)

    plt.axvline(max_ratio,color='r',alpha=0.2)

    if ylogscale:
        plt.yscale('log')
    
    if (ymin is not None) & (ymax is not None):
        plt.ylim(ymin,ymax)
    
    #boss has also requested a plot of means with standard deviations for each ratio bin.
    plt.scatter(0.15,mean_one,color='black',s=100,zorder=2,marker='*',label='Average bin {}'.format(data_type_names[int(ydata_type)]))
    plt.errorbar(0.15,mean_one,std_one,color='black',capsize=3,linewidth=4,
                 linestyle='None',label='average errorbar size')
    
    plt.scatter(0.40,mean_two,color='black',s=100,zorder=2,marker='*')
    plt.errorbar(0.40,mean_two,std_two,color='black',capsize=3,linewidth=4,
                 linestyle='None')
    
    plt.scatter(0.60,mean_three,color='black',s=100,zorder=2,marker='*')
    plt.errorbar(0.60,mean_three,std_three,color='black',capsize=3,linewidth=4,
                 linestyle='None')
    
    plt.legend(loc='upper left',fontsize=15)
    #plt.show()
    
    if savefig:
        plt.savefig(path_to_dir+png_name,dpi=300)
    
    if not savefig:
        plt.show()
    
    plt.close()
    print('means:',mean_one,mean_two,mean_three)
    print('std:',std_one,std_two,std_three)
    
    
def create_all_plots(path_to_dir, convflag, mask_gal, savefig=False):
    
    data_types = [0,1,2]
    
    for i in data_types:
        create_maskiter_plots(path_to_dir, convflag=convflag, mask_gal=mask_gal, ydata_type=i, savefig=savefig)
        print()
    
if __name__ == '__main__':
    
    #galaxy_name = input('enter galaxy name (with proper capitalization): ')
    galaxy_name='NGC3364'
    print()
    print('runitallx3(galaxy_name,convflag=False,mask_gal=True,termination_ratio=0.70)')
    print()
    print('create_maskiter_plots(path_to_dir, convflag, mask_gal, max_iter=10, ydata_type=0, ymin=None, ymax=None, ylogscale=False, savefig=False)')
    print()
    print('create_all_plots(path_to_dir, convflag, mask_gal, savefig=False)')
    print('--> 0 for Re, 1 for nser, 2 for sky')
    
    #runitallx3(galaxy_name,termination_ratio=0.7,mask_gal=False,convflag=False)