''' This script is a rudimentary collection of functions to evaluate new galfit parameters and compare with previous iternations of the output parameter tables. I am adding them here for purposes of maintaining a relatively tidy Jupyter notebook interface. :-)'''

import numpy as np
from matplotlib import pyplot as plt
from astropy.table import Table
import os

homedir = os.getenv("HOME")

#load old, new tables
#create functions for plotting purposes
#etc.

class catalogs:

    def __init__(self):
        
        #magphys table
        self.magphys = Table.read(homedir+'/Desktop/v2-20220820/vf_v2_magphys_legacyExt_final.fits')

        #catalog with WISESize flags
        self.vf = Table.read(homedir+'/Desktop/VF_WISESIZE.fits')

        self.w1_new = Table.read(homedir+'/Desktop/galfit_files/galfit_W1_03112024.fits')
        self.w3_new = Table.read(homedir+'/Desktop/galfit_files/galfit_W3_03112024.fits')

        self.w1_old = Table.read(homedir+'/Desktop/galfit_files/galfit_W1_2023.fits')
        self.w3_old = Table.read(homedir+'/Desktop/galfit_files/galfit_W3_2023.fits')
        
        self.cutcats()
        
    def cutcats(self):
        #wI can apply t-type, SNR
        snr_flag = self.vf['SNRflag']
        t_flag = self.vf['t_flag']

        #certain GALFIT cuts
        n_flag = (self.w3_new['CN']<6) & (self.w1_new['CN']<6)
        fail_flag = ~(self.w3_new['CRE']==0) & ~(self.w1_new['CRE']==0) & ~(self.w3_old['CRE']==0) & ~(self.w1_old['CRE']==0)
        numerical_flag = ~(self.w1_new['CNumerical_Error']) & ~(self.w3_new['CNumerical_Error'])

        #magphys flags (including quality flags from various sources)
        mag_flag = (self.magphys['magphysFlag']) & (self.magphys['logSFR_med']>-1.398) & (self.magphys['logMstar_med']>8.26) & ((self.magphys['logSFR_med']-self.magphys['logMstar_med'])>-11.5)


        all_flags = snr_flag & t_flag & n_flag & fail_flag & numerical_flag & mag_flag

        self.vfcut = self.vf[all_flags]
        self.w1_oldcut = self.w1_old[all_flags]
        self.w3_oldcut = self.w3_old[all_flags]
        self.w1_newcut = self.w1_new[all_flags]
        self.w3_newcut = self.w3_new[all_flags]
        self.magphys_cut = self.magphys[all_flags]

        print('Final subsample size:',len(self.vfcut))
        print('magflag+SNR+t-type only:',len(self.vf[snr_flag&t_flag&mag_flag])) 
        
        print('number of galaxies with w1 galfit errors (old):',len(self.vf[self.w1_old['CERROR']==1.]))
        
        print('number of galaxies with w1 galfit errors (new):',len(self.vf[~self.w1_new['CNumerical_Error']]))
        
        print('number of galaxies with w3 galfit errors (old):',len(self.vf[self.w3_old['CERROR']==1.]))
        
        print('number of galaxies with w3 galfit errors (new):',len(self.vf[~self.w3_new['CNumerical_Error']]))
            
    def comp_oldnew(self,zoom=True):

        plt.figure(figsize=(9,7))

        plt.subplot(2,2,1)
        plt.scatter(self.w1_oldcut['CRE'], self.w1_oldcut['CRE']/self.w1_newcut['CRE'], c=self.magphys_cut['logMstar_med'])
        plt.clim(7.5,10.5)
        plt.axhline(1,color='indigo',ls='dashed')
        plt.xscale('log')
        plt.yscale('log')
        plt.ylabel(r'Re$_{old}$/Re$_{new}$',fontsize=12)
        plt.xlabel(r'Re$_{old}$',fontsize=12)
        plt.title('W1')

        plt.subplot(2,2,2)
        plt.scatter(self.w3_oldcut['CRE'], self.w3_oldcut['CRE']/self.w3_newcut['CRE'], c=self.magphys_cut['logMstar_med'])
        plt.clim(7.5,10.5)
        plt.axhline(1, color='indigo',ls='dashed')
        plt.xscale('log')
        plt.yscale('log')
        plt.title('W3')
        plt.ylabel(r'Re$_{old}$/Re$_{new}$',fontsize=12)
        plt.xlabel(r'Re$_{old}$',fontsize=12)
        plt.colorbar(label='log(Mstar)')

        plt.subplot(2,2,3)
        plt.scatter(self.w1_oldcut['CN'], self.w1_oldcut['CN']/self.w1_newcut['CN'], c=self.magphys_cut['logMstar_med'])
        plt.clim(7.5,10.5)
        plt.axhline(1,color='indigo',ls='dashed')
        plt.xscale('log')
        plt.yscale('log')
        plt.ylabel(r'n$_{old}$/n$_{new}$',fontsize=12)
        plt.xlabel(r'n$_{old}$',fontsize=12)
        plt.title('W1')

        plt.subplot(2,2,4)
        plt.scatter(self.w3_oldcut['CN'], self.w3_oldcut['CN']/self.w3_newcut['CN'], c=self.magphys_cut['logMstar_med'])
        plt.clim(7.5,11)
        plt.axhline(1,color='indigo',ls='dashed')
        plt.xscale('log')
        plt.yscale('log')
        plt.title('W3')
        plt.ylabel(r'n$_{old}$/n$_{new}$',fontsize=12)
        plt.xlabel(r'n$_{old}$',fontsize=12)
        plt.colorbar(label='log(Mstar)')

        plt.subplots_adjust(wspace=0.3, hspace=0.4)
        plt.show()

        if zoom:
            plt.figure(figsize=(9,7))

            plt.subplot(2,2,1)
            plt.scatter(self.w1_oldcut['CRE'], self.w1_oldcut['CRE']/self.w1_newcut['CRE'], c=self.magphys_cut['logMstar_med'])
            plt.clim(7.5,10.5)
            plt.axhline(1,color='indigo',ls='dashed')
            plt.ylabel(r'Re$_{old}$/Re$_{new}$',fontsize=12)
            plt.xlabel(r'Re$_{old}$',fontsize=12)
            plt.title('W1')
            plt.ylim(0,2)

            plt.subplot(2,2,2)
            plt.scatter(self.w3_oldcut['CRE'], self.w3_oldcut['CRE']/self.w3_newcut['CRE'], c=self.magphys_cut['logMstar_med'])
            plt.clim(7.5,10.5)
            plt.axhline(1, color='indigo',ls='dashed')
            plt.title('W3')
            plt.ylabel(r'Re$_{old}$/Re$_{new}$',fontsize=12)
            plt.xlabel(r'Re$_{old}$',fontsize=12)
            plt.colorbar(label='log(Mstar)')
            plt.ylim(0,2)
            plt.xlim(-0.5,50)

            plt.subplot(2,2,3)
            plt.scatter(self.w1_oldcut['CN'], self.w1_oldcut['CN']/self.w1_newcut['CN'], c=self.magphys_cut['logMstar_med'])
            plt.clim(7.5,10.5)
            plt.axhline(1,color='indigo',ls='dashed')
            plt.ylabel(r'n$_{old}$/n$_{new}$',fontsize=12)
            plt.xlabel(r'n$_{old}$',fontsize=12)
            plt.title('W1')
            plt.ylim(0,2)
            plt.xlim(0,6)

            plt.subplot(2,2,4)
            plt.scatter(self.w3_oldcut['CN'], self.w3_oldcut['CN']/self.w3_newcut['CN'], c=self.magphys_cut['logMstar_med'])
            plt.clim(7.5,10.5)
            plt.axhline(1,color='indigo',ls='dashed')
            plt.title('W3')
            plt.ylabel(r'n$_{old}$/n$_{new}$',fontsize=12)
            plt.xlabel(r'n$_{old}$',fontsize=12)
            plt.colorbar(label='log(Mstar)')
            plt.ylim(0,2)
            plt.xlim(0,6)

            plt.subplots_adjust(wspace=0.3, hspace=0.4)
            plt.show()

    def newbands(self):
        return
    
    
if __name__=='__main__':
    print('''USAGE: ''')

#plt.axline((0, 0), slope=1, color='indigo',ls='dashed')
