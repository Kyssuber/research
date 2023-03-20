'''
BITMASK DICTIONARY:
https://www.legacysurvey.org/dr9/bitmasks/

INPUT:
--FITS file of optical (r-band) mask

OUTPUT:
--FITS file of 12-micron mask (i.e., the inflated r-band mask)
--mosaic .png displaying the optical mask and the enlarged 12-micron version, for comparison/testing purposes

NOTE: may not apply to cases where central galaxy is not in SGA; may not mask out other galaxies in cutout field (in which instance I think we would run galfit with 2+ sersic obj. anyhow).
'''


#call packages n' libraries n' libraries n' packages
import os
homedir = os.getenv("HOME")

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

#import sys
#np.set_printoptions(threshold=sys.maxsize)


#display cutouts and r-band mask
def diagnostic_one(im_r,im_w3,mask_one,r_vmax=0.8):
    #display r-band image, 12-micron image, r-band mask
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3,figsize=(14,14))
    ax1.imshow(im_r,vmax=r_vmax,origin='lower')
    ax2.imshow(im_w3,origin='lower')
    ax3.imshow(mask_one,origin='lower')
    ax1.set_title('r-band image',fontsize=15)
    ax2.set_title('12-micron image',fontsize=15)
    ax3.set_title('mask image',fontsize=15)
    plt.show()


#create dictionary with exponents (n) and their resulting 2**n value
def create_power_dictionary(min_power=0, max_power=12):

    two_exp = np.arange(min_power,max_power+1,1)
    two_array = two_exp.copy()
    
    for i in two_array:
        two_array[i] = 2**i
        
    powers_of_two = dict(zip(two_exp,two_array))
    
    return(powers_of_two)


#decompose every pixel of the r-band mask into a list of the constituent powers
#this function is run for EVERY pixel. yeesh. :-(
def power_decomposition(old_pixel_val, power_dict):
    
    #create empty decomposition list (which will contain the exponents that comprise the mask)
    decomposition = []
    
    #isolate old pixel value
    n = old_pixel_val
    
    #loop through every exponent in power dictionary
    for i in range(len(power_dict)):
        
        #begin with largest power of 2
        dict_index = len(powers_of_two)-(1+i)
        
        #if n is already zero, finish
        if n == 0:
            break
        
        #if n is *exactly* equal to a power of 2, append to decomposition list and finish
        if n == powers_of_two[dict_index]:
            decomposition.append(dict_index)
            break
        
        #begin with largest power of 2 remaining
        dict_index = len(powers_of_two)-(1+i)
        #determine the number to subtract from n
        byte_to_subtract = powers_of_two[dict_index]

        #in the case where the number to subtract is larger than n, skip to the next power of 2
        if (n-byte_to_subtract)<0:
            pass
        
        #else, subtract number from n and append the corresponding power of 2 to the decomposition list
        if (n-byte_to_subtract)>=0:
            n -= byte_to_subtract
            decomposition.append(dict_index)
    
    #return a list of every powers-of-2 exponent constituting the pixel value
    return(decomposition)


#use decomposition list to determine whether pixel should be masked (neighboring object) or not masked (SGA galaxy or background)
def set_mask_value(power_dict, single_decomposition_list):
    
    if single_decomposition_list:  #checks if list is empty (IF empty, then px must have been zero)
        
        #12 corresponds to SGA galaxy mask, per LegacySurvey website; 11 corresponds to star
        #(in case of tested galaxy IC2389, bright star overlaps galaxy)
        #NOTE: 3 corresponds to saturated r-band pixel...unsure whether to include in this list of conditions
        if (12 in single_decomposition_list) & (11 not in single_decomposition_list): 
            new_pixel_val = 1 
        else:    
            new_pixel_val = 0  #if mask is not SGA galaxy, or IS galaxy as well as star, 
            #ORRR is flagged as some sort of alternative mask, then change value to zero
    else:
        new_pixel_val = 1  #if original pixel value was zero, then no mask --> change to value 1
    
    return(new_pixel_val) #then set this pixel val as the value of whichever px the loop is currently operating on


def diagnostic_two(mask_edit,im_r,vmax=r_vmax):
#display mask_edit and (mask_edit*r-band image), sort of a diagnostic to determine whether the above code functioned correctly
    fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(12,12))
    ax1.imshow(mask_edit,origin='lower')
    ax2.imshow(mask_edit*im_r,vmax=r_vmax,origin='lower')
    ax1.set_title('edited mask',fontsize=15)
    ax2.set_title('edited mask * r-band image',fontsize=15)
    plt.show()


#this was a terrible experience.
def shrink_the_maskies(r_mask, w3_image):
    
    '''
    --create empty 12-micron mask image --> will have same shape as 12-micron image (different number of pixels 
      than r-band image, but covers same area as r-band image)
    --calculate number of r-band pixels that lie inside a w3 pixel
    --convert r_mask row of width N into NxN blocks, each of which eventually translating
      to a single w3 pixel value
    --output --> w3_mask counterpart to r_mask  
    '''
    
    #create empty 12-micron mask --> has same shape as 12-micron image and covers same area as r-band image
    w3_mask=np.zeros(np.shape(w3_image))
    
    #divide width of r-band image into block rows
    #result should be len(w3_mask) rows, all with length len(r_mask)
    r_mask_rows = np.array_split(r_mask,len(w3_mask),axis=0)   
    
    #for each row in empty_w3_mask,
    for row in range(len(w3_mask)):
        
        #isolate corresponding row block in r_mask
        block_row = r_mask_rows[row]
        
        #split block_row into block with dimensions ~nxn
        array_of_blocks = np.array_split(block_row,len(w3_mask),axis=1)
        
        #for every pixel in the w3_mask row
        for i in range(len(w3_mask[row])):
            
            #give the isolated block a snazzy variable name
            single_block = array_of_blocks[i]
            
            #if any pixels within block array_of_blocks[i] are masked, then assign w3_mask pixel to 0
            if 0 in single_block:
                w3_mask[row][i] = 0
            #if none of the pixels in this block are masked, assign w3_mask pixel to 1
            else:
                w3_mask[row][i] = 1
        
        #isolate an n x len(r_mask) block_row, which will thereafter be split into blocks with dimensions nxn
        #if row=0, then block_row width is 0 to n; if row=1, then width is 1*n to 2*n, etc.
        #block_row = r_mask[(row*n):((row+1)*n), :]              
                
    return(w3_mask)


def diagnostic_three(mask_edit,w3_mask,im_w3):
    #oh, yay --> another diagnostic test.
    #displays r-band mask, the new w3_mask, and (w3_mask * w3_image)
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3,figsize=(14,14))
    ax1.imshow(mask_edit,origin='lower')
    ax2.imshow(w3_mask,origin='lower')
    ax3.imshow(im_w3*w3_mask,origin='lower')
    ax1.set_title('r-band mask',fontsize=15)
    ax2.set_title('w3 mask',fontsize=15)
    ax3.set_title('w3-band image',fontsize=15)
    plt.show()


#convert 0 pixels (masked) to 1; and 1 pixels (unmasked) to 0.
def test_to_real(r_mask, w3_mask):
    
    r_mask_real = np.abs(r_mask - 1)
    w3_mask_real = np.abs(w3_mask - 1)
    
    return(r_mask_real, w3_mask_real)





'''
~~~~~EXECUTION EXAMPLE~~~~~
~~~(no diagnostic tests)~~~

#define path to test galaxy
path_to_folder = homedir+'/IC2389/'
im_r_path = path_to_folder+'IC2389-custom-image-r.fits.fz'
im_w3_path = path_to_folder+'IC2389-custom-image-w3.fits.fz'
mask_one_path = path_to_folder+'IC2389-custom-maskbits.fits.fz'

#extract pixel data from .fits files
im_r = fits.getdata(im_r_path)
im_w3 = fits.getdata(im_w3_path)
mask_one = fits.getdata(mask_one_path)  #r-band mask
mask_edit = mask_one.copy()  #r-band mask with 0s for background, 1s for masked pixels

#setup dictionary
powers_of_two = create_power_dictionary()

#create a copy of mask_one, which will be replaced pixel by pixel with updated 0 or 1 values
mask_edit = mask_one.copy()

#the following loop requires ~10 seconds to execute per galaxy. not ideal.
#for every row in the mask image,
for row in range(np.shape(mask_edit)[0]):
    
    #and for every pixel in this row,
    for px in range(len(mask_edit[row])):
        
        #isolate the old pixel value
        old_pixel = mask_edit[row][px]
        
        #extract the 2**n exponent list
        decomposition_list = power_decomposition(old_pixel, powers_of_two)
        
        #calculate the new pixel value, either 0 or 1
        new_px_value = set_mask_value(powers_of_two, decomposition_list)
        
        #assign this new value to the corresponding pixel in mask_edit
        mask_edit[row][px] = new_px_value

#create w3 mask equivalent
w3_mask = shrink_the_maskies(r_mask = im_r, w3_image = im_w3)

#convert to *actual* mask pixel values for which masked pixels have a value of 1 and unmasked pixels have a value of 0
#note that the reverse was true ONLY for diagnostic purposes (so I could easily apply a mask to a cutout with cutout*mask ... #masked pixels would be zero in the cutout and unmasked pixels would retain their original value). 
r_mask, w3_mask = test_to_real(r_mask, w3_mask)

#can be looped through to run for multiple galaxies.

#fin.
'''


