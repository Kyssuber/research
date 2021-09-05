import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.table import Table
from photutils.aperture import CircularAperture
import os

def snr(gal_sample):
    snr15_list = []
    snr30_list = []
    snr15_flag = []
    snr30_flag = []
    for n in gal_sample:
        try:
            vfid=str(n['VFID'])
        except:
            vfid = n['galname']
            vfid = vfid[0:8]
        filepath = '/mnt/astrophysics/wisesize/'+vfid+'/unwise-'+vfid+'-w3-'
        im = filepath+'img-m.fits'
        std = filepath+'std-m.fits'

        hdu_im=fits.open(im)[0]
        hdu_std=fits.open(std)[0]
        #fig,(ax1,ax2) = plt.subplots(1,2)
        #ax1.imshow(hdu_im.data,origin='lower')
        #ax2.imshow(hdu_std.data,origin='lower')

        img_len = len(hdu_im.data)
        position = (int(img_len/2),int(img_len/2))

        len_arcsec = np.abs(hdu_im.header['NAXIS1']*hdu_im.header['CD1_1'])*3600
        arcsec_per_pixel = len_arcsec/hdu_im.header['NAXIS1']
        print(vfid)
        #print('arcsec per pixel: %.2f'%arcsec_per_pixel)

        radius_15 = 15/arcsec_per_pixel
        radius_30 = radius_15*2

        aper15 = CircularAperture(position,radius_15)
        mask15 = aper15.to_mask(method='center')
        aper30 = CircularAperture(position,radius_30)
        mask30 = aper30.to_mask(method='center')

        data15_im = mask15.multiply(hdu_im.data)
        data15_std = mask15.multiply(hdu_std.data)
        data30_im = mask30.multiply(hdu_im.data)
        data30_std = mask30.multiply(hdu_std.data)

        data15_im = data15_im[data15_im!=0]
        data15_std = data15_std[data15_std!=0]
        data30_im = data30_im[data30_im!=0]
        data30_std = data30_std[data30_std!=0]

        noise15 = np.sqrt(np.sum(data15_std**2))
        noise30 = np.sqrt(np.sum(data30_std**2))
        signal15 = np.sqrt(np.sum(data15_im**2))
        signal30 = np.sqrt(np.sum(data30_im**2))

        snr15 = signal15/noise15
        snr15_list.append(snr15)
        snr30 = signal30/noise30
        snr30_list.append(snr30)
        if (snr15 > 10.) & (str(snr15) != 'inf'):
            print(vfid, 'snr15 > 10')
            snr15_flag.append(1)
        if (snr15 < 10.) & (str(snr15) != 'inf'):
            snr15_flag.append(0)
        if str(snr15) == 'inf':
            snr15_flag.append(0)
        if (snr30 > 10.) & (str(snr30) != 'inf'):
            print(vfid, 'snr30 > 10')
            snr30_flag.append(1)
        if (snr30 < 10.) & (str(snr30) != 'inf'):
            snr30_flag.append(0)
        if str(snr30) == 'inf':
            snr30_flag.append(0)

    snr15_list = [round(num,1) for num in snr15_list]
    #print('snr15:')
    #print(snr15_list)
    snr15_array = np.asarray(snr15_list)

    snr30_list = [round(num,1) for num in snr30_list]
    #print('snr30:')
    #print(snr30_list)
    snr30_array = np.asarray(snr30_list)

    snr15_flag = np.asarray(snr15_flag)
    snr30_flag = np.asarray(snr30_flag)

    count15 = 0
    count30 = 0
    count_both = 0
    print(len(snr15_flag),len(snr30_flag))
    for i in range(0,len(gal_sample)):
        if snr15_flag[i] == 1:
            count15 += 1
        if snr30_flag[i] == 1:
            count30 += 1
        if (snr15_flag[i] == 1) & (snr30_flag[i] == 1):
            count_both += 1
    print('SNR>10 for 15arc: ',count15,'of ',len(snr15_flag))
    print('SNR>10 for 30arc: ',count30,'of ',len(snr15_flag))
    print('# galaxies with SNR>10 for both: ',count_both)
            
    os.chdir(homedir)
    np.savetxt('snr15.txt',snr15_array,fmt="%f")
    np.savetxt('snr30.txt',snr30_array,fmt="%f")
    np.savetxt('snr15_flag.txt',snr15_flag,fmt="%s")
    np.savetxt('snr30_flag.txt',snr30_flag,fmt="%s")

    return(snr15_flag)
    return(snr30_flag)



if __name__ == '__main__':
    homedir = os.getenv("HOME")
    vf = Table.read(homedir+'/vf_north_v1_main.fits')
    print('snr(gal_sample)')
    print('sample: vf')
    print("read table: ascii.read('___.txt',format='no_header')")




