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
        
        self.r_new = Table.read(homedir+'/Desktop/galfit_files/galfit_r_03112024.fits')
        
        #self.w1_new = Table.read(homedir+'/Desktop/galfit_files/vf_v2_galfit_W1-fixBA.fits')
        #self.w3_new = Table.read(homedir+'/Desktop/galfit_files/vf_v2_galfit_W3-fixBA.fits')
        
        self.cutcats()
        
    def cutcats(self):
        #wI can apply t-type, SNR
        snr_flag = self.vf['SNRflag']
        t_flag = self.vf['t_flag']

        #certain GALFIT cuts (must NOT have n>6 and NOT have zeros in rows and NOT have a numerical error)
        n_flag = (self.w3_new['CN']<6) & (self.w1_new['CN']<6)
        fail_flag = ~(self.w3_new['CRE']==0) | ~(self.w1_new['CRE']==0) | ~(self.w3_old['CRE']==0) | ~(self.w1_old['CRE']==0)
        numerical_flag = ~(self.w1_new['CNumerical_Error']) & ~(self.w3_new['CNumerical_Error'])

        #magphys flags (including quality flags from various sources)
        mag_flag = (self.magphys['magphysFlag']) & (self.magphys['logSFR_med']>-1.398) & (self.magphys['logMstar_med']>8.26) & ((self.magphys['logSFR_med']-self.magphys['logMstar_med'])>-11.5)


        all_flags = snr_flag & t_flag & n_flag & fail_flag & numerical_flag & mag_flag

        self.vfcut = self.vf[all_flags]
        self.w1_oldcut = self.w1_old[all_flags]
        self.w3_oldcut = self.w3_old[all_flags]
        self.w1_newcut = self.w1_new[all_flags]
        self.w3_newcut = self.w3_new[all_flags]
        
        self.rcut = self.r_new[all_flags]
        self.magphys_cut = self.magphys[all_flags]

        print('Final subsample size (inc. SNR, t-type, magphys, and W1+W3 galfit flags):',len(self.vfcut))
        print('magflag+SNR+t-type only:',len(self.vf[snr_flag&t_flag&mag_flag])) 
                
        print('')
        print('for full catalog (galfit errors only):')
        print(f"fraction galaxies with w1 galfit errors (old): {np.sum(self.w1_old['CERROR']==1.)/len(self.vf):.03f}")
                
        print(f"fraction galaxies with w1 galfit errors (new): {np.sum(self.w1_new['CNumerical_Error'])/len(self.vf):.03f}")
        
        print(f"fraction galaxies with w3 galfit errors (old): {len(self.vf[self.w3_old['CERROR']==1.])/len(self.vf):.03f}")
        
        print(f"fraction galaxies with w3 galfit errors (new): {np.sum(self.w3_new['CNumerical_Error'])/len(self.vf):.03f}")
        
        print('')
        print('for subsample galaxies...(applying SNR flag and galfit error flag):')
        
        num_snr_w1old_flag = (self.w1_old['CERROR']==1.) & (self.vf['SNRflag'])
        num_snr_w3old_flag = (self.w3_old['CERROR']==1.) & (self.vf['SNRflag'])
        num_snr_w1new_flag = (self.w1_new['CNumerical_Error']) & (self.vf['SNRflag'])
        num_snr_w3new_flag = (self.w3_new['CNumerical_Error']) & (self.vf['SNRflag'])
        
        print(f'fraction galaxies with w1 galfit errors (old): {np.sum(num_snr_w1old_flag)/len(self.vf):.03f}')
        print(f'fraction galaxies with w1 galfit errors (new): {np.sum(num_snr_w1new_flag)/len(self.vf):.03f}')
        print(f'fraction galaxies with w3 galfit errors (old): {np.sum(num_snr_w3old_flag)/len(self.vf):.03f}')
        print(f'fraction galaxies with w3 galfit errors (new): {np.sum(num_snr_w3new_flag)/len(self.vf):.03f}')
        
        '''
        galfit_W1 = self.w3_new
        print("W1 results:")
        print("------------")
        W1galfitflag = galfit_W1['MAG'] > 0
        nflag = galfit_W1['Numerical_Error']
        cnflag = galfit_W1['CNumerical_Error']
        print(f"percent success w/out convolution: {np.sum(W1galfitflag & ~nflag)/len(galfit_W1):.03f}")
        print(f"percent success with convolution: {np.sum(W1galfitflag & ~cnflag)/len(galfit_W1):.03f}")
        '''
        
    def comp_oldnew(self,zoom=False):

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
        plt.clim(7.5,10.5)
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

    def rbandcomp(self):
            
        rflag = ~(self.rcut['CNumerical_Error']) & (self.rcut['CN']<6) & ~(self.rcut['CRE']==0)
        print(f'frac points without r-band errors (/518): {len(self.rcut[rflag])/len(self.rcut):.03f}')
        
        rvals = [self.rcut['CRE'][rflag]*0.262,self.rcut['CN'][rflag],
                 self.rcut['CPA'][rflag],self.rcut['CAR'][rflag]]
        w3vals = [self.w3_newcut['CRE'][rflag]*2.75,self.w3_newcut['CN'][rflag],
                  self.w3_newcut['CPA'][rflag],self.w3_newcut['CAR'][rflag]]
        w1vals = [self.w1_newcut['CRE'][rflag]*2.75,self.w1_newcut['CN'][rflag],
                  self.w1_newcut['CPA'][rflag],self.w1_newcut['CAR'][rflag]]
        labels = ['Re (arcsec)','nser','PA','B/A']
        
        plt.figure(figsize=(14,2.5))
        for n in range(len(rvals)):
            plt.subplot(1,4,n+1)
            plt.scatter(rvals[n],w3vals[n],color='b',s=40,alpha=0.7,label='W3')
            plt.scatter(rvals[n],w1vals[n],color='r',s=40,alpha=0.7,label='W1')
            plt.axline((0, 0), slope=1, color='indigo',ls='dashed')
            plt.xlabel(f'r-band {labels[n]}',fontsize=12)
            plt.ylabel(f'WISE {labels[n]}',fontsize=12)
            if (n==0):
                plt.xscale('log')
                plt.yscale('log')
            if (n==0)|(n==3):
                plt.legend(framealpha=0.2)
            
        plt.subplots_adjust(wspace=0.45)
        plt.show()

        return
    
    
if __name__=='__main__':
    print('''USAGE: ''')
