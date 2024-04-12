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
from reproject import reproject_interp


#display cutouts and r-band mask
def diagnostic_one(im_r,im_w3,rband_mask,r_vmax=0.8):
    #display r-band image, 12-micron image, r-band mask
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3,figsize=(14,14))
    ax1.imshow(im_r,vmax=r_vmax,origin='lower')
    ax2.imshow(im_w3,origin='lower')
    ax3.imshow(mask_one,origin='lower')
    ax1.set_title('r-band image',fontsize=15)
    ax2.set_title('12-micron image',fontsize=15)
    ax3.set_title('mask image',fontsize=15)
    plt.show()


def mask_to_bool(path_to_rmask):
    
    '''
    -Isolate pixels that DO NOT BELONG to central galaxy, assuming it is in the SGA; output rmask with 0s corresponding to unmasked px, 1s corresponding to masked px.
    
    -First bool (sat_galpx_flag) --> selects saturated pixels that overlap with the galaxy
    -Second bool (background_flag & non_galaxy_flag) --> selects pixels that are masked (value > 0) AND have a value smaller than the galaxy flag

    '''
    
    r_mask = fits.open(path_to_rmask)[1]
    
    sat_galpx_flag = np.array(r_mask.data) > 4096
    background_flag = np.array(r_mask.data) > 0
    non_galaxy_flag = np.array(r_mask.data) < 4096
    
    new_rmask = (sat_galpx_flag) | (background_flag & non_galaxy_flag)
    
    return(new_rmask)

    
def project_rmask(new_rmask, path_to_imw3):

    w3_hdu = fits.open(path_to_imw3)[1]
    wise_mask, footprint = reproject_interp(new_rmask, wise_hdu.header)
    
    return(w3_hdu, wise_mask)
                  
                  
#oh, yay. another diagnostic test.
#displays r-band mask, the new w3_mask, and (w3_mask * w3_image)
def diagnostic_two(new_rmask, wise_mask, w3_hdu):
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3,figsize=(14,14))
    ax1.imshow(new_rmask,origin='lower')
    ax2.imshow(wise_mask,origin='lower')
    ax3.imshow(w3_hdu*(~wise_mask),origin='lower')
    ax1.set_title('r-band mask',fontsize=15)
    ax2.set_title('w3 mask',fontsize=15)
    ax3.set_title('w3-band image',fontsize=15)
    plt.show()
    
    
    
    