'''
GOAL: analyze GALFIT output parameters for Virgowise subsample galaxies.

pixel_scale = {'FUV':1.5,'NUV':1.5,'g':0.262,'r':0.262,'z':0.262,'W1':2.75,'W2':2.75,'W3':2.75,'W4':2.75}
psf_oversampling = {'FUV':1,'NUV':1,'g':1,'r':1,'g':1,'W1':1,'W2':1,'W3':1,'W4':1}
mag_zeropoint = {'FUV':22.5,'NUV':22.5,'g':22.5,'r':22.5,'g':22.5,'W1':22.5,'W2':22.5,'W3':22.5,'W4':22.5}
image_resolution = {'FUV':6,'NUV':6,'g':1.5,'r':1.5,'z':1.5,'W1':6.1,'W2':6.4,'W3':6.5,'W4':12}
'''

import numpy as np
from matplotlib import pyplot as plt
from astropy.table import Table

import os
homedir = os.getenv("HOME")

class catalogs:
    
    def __init__(self,conv=False):
        path_to_dir = homedir+'/Desktop/v2-20220820/'
        self.v2_env = Table.read(path_to_dir+'vf_v2_environment.fits')
        self.v2_main = Table.read(homedir+'/v2_snrcoadd.fits')
        self.magphys = Table.read(path_to_dir+'vf_v2_magphys_07-Jul-2022.fits')
        self.z0mgs = Table.read(path_to_dir+'vf_v2_z0mgs.fits')
        
        self.conv = conv
        
        if self.conv==False:
            self.rdat = Table.read(homedir+'/output_params_r_nopsf.fits')
            self.w3dat = Table.read(homedir+'/output_params_W3_nopsf.fits')
        if self.conv==True:
            self.rdat = Table.read(homedir+'/output_params_r_psf.fits')
            self.w3dat = Table.read(homedir+'/output_params_W3_psf.fits')
            
        self.cut_cats()
            
    def cut_cats(self):
        subsample_flag = self.v2_main['sgacut_flag']
        
        self.v2_envcut = self.v2_env[subsample_flag]
        self.v2_maincut = self.v2_main[subsample_flag]
        self.magphyscut = self.magphys[subsample_flag]
        self.z0mgscut = self.z0mgs[subsample_flag]
        
        self.re_rband = self.rdat['re']
        self.re_w3band = self.w3dat['re']

        self.PA_rband = self.rdat['PA']
        self.PA_w3band = self.w3dat['PA']

        self.BA_rband = self.rdat['BA']
        self.BA_w3band = self.w3dat['BA']

        group_flag = (self.re_rband == 0.0) | (self.re_w3band == 0.0)
        r_flag = (self.re_rband == 0.0)
        w3_flag = (self.re_w3band == 0.0)
        err_flag = (self.rdat['err_flag']==1) | (self.w3dat['err_flag']==1)
        
        n_fails_w3 = len(self.re_rband[r_flag])
        n_fails_r = len(self.re_w3band[w3_flag])

        self.re_rband_cut = self.re_rband[~group_flag & ~err_flag]
        self.re_w3band_cut = self.re_w3band[~group_flag & ~err_flag]

        self.PA_rband_cut = self.PA_rband[~group_flag & ~err_flag]
        self.PA_w3band_cut = self.PA_w3band[~group_flag & ~err_flag]

        self.BA_rband_cut = self.BA_rband[~group_flag & ~err_flag]
        self.BA_w3band_cut = self.BA_w3band[~group_flag & ~err_flag]
        
        #apply final-ish cut to envcut and maincut catalogs
        self.v2_envcut = self.v2_envcut[~group_flag & ~err_flag]
        self.v2_maincut = self.v2_maincut[~group_flag & ~err_flag]
        self.magphyscut = self.magphyscut[~group_flag & ~err_flag]
        self.z0mgscut = self.z0mgscut[~group_flag & ~err_flag]
        
        #define env flags
        self.clusflag = self.v2_envcut['cluster_member']
        self.rgflag = self.v2_envcut['rich_group_memb']
        self.pgflag = self.v2_envcut['poor_group_memb']
        self.filflag = self.v2_envcut['filament_member']
        self.fieldflag = self.v2_envcut['pure_field']

        n_tot = len(self.re_rband)     #should be the same for both r-band and w3

        print(f'No GALFIT data for {n_fails_w3} w3 galaxies and {n_fails_r} r galaxies.')
        print(f'Total number of galaxies with GALFIT error flags: {int(np.sum(np.ones(len(err_flag))*err_flag))}')
        print(f'Total number of galaxies: {n_tot}')

        self.sizerats = (self.re_w3band_cut*2.75)/(self.re_rband_cut*0.262)
        self.PArats = self.PA_w3band_cut/self.PA_rband_cut

    def envbins(self, savefig=False):

        #plot num of galaxies in each bin
        fig = plt.figure(figsize=(10,6))
        ax = fig.add_subplot(111)
        ind=np.arange(0,5,1)
        width=0.15
        
        x = np.arange(0,5,1)
        env_bin_dat = np.array([len(self.v2_envcut[self.clusflag]), len(self.v2_envcut[self.rgflag]), 
                                len(self.v2_envcut[self.pgflag]), len(self.v2_envcut[self.filflag]), 
                                len(self.v2_envcut[self.fieldflag])])

        plt.scatter(x,env_bin_dat,color='purple')
        plt.title('Number of galaxies in each bin (tot='+str(len(self.v2_envcut))+')',fontsize=18)
        plt.ylabel('# Galaxies',fontsize=18)

        ax.set_xlim(-width,len(ind)-1.5*width)
        xTickMarks= ['Cluster','Rich \n Group','Poor \n Group','Filament','Field']
        plt.xticks(ind, xTickMarks, rotation=10, fontsize=18)
        plt.xlim(-0.5,4.5)
        plt.ylim(0,200)
        plt.grid(color='purple',alpha=0.2)
        plt.show()

        if savefig==True:
            plt.savefig(homedir+'/Desktop/overleaf_figures/envbins.png', dpi=300)

    def env_means(self, combine_mid=False, convflag=False, savefig=False):    
        
        ratios = self.sizerats
        
        re_data = [ratios[self.clusflag],ratios[self.rgflag],ratios[self.pgflag],
                   ratios[self.filflag],ratios[self.fieldflag]]
        mean = []
        err = []
        index = np.arange(1,6,1)
        env_names = ['Cluster','Rich \n Group','Poor \n Group','Filament','Field']

        for j,i in enumerate(re_data):      #j==index, i==value
            mean.append(np.mean(i))
            err.append(np.std(i)/np.sqrt(len(i)))
        
        err_color = 'orangered'
        plt.figure(figsize=(10,6))
        plt.scatter(index,mean,color='blue',s=40,zorder=2)
        plt.errorbar(index,mean,yerr=err,fmt='None',color=err_color,zorder=1)

        xmin,xmax = plt.xlim()
        xfield = np.linspace(xmin,xmax,50)
        ymax = np.ones(50)*(mean[-1] + err[-1])
        ymin = np.ones(50)*(mean[-1] - err[-1])
        plt.fill_between(xfield,ymax,ymin,color=err_color,alpha=.1)

        plt.xticks(index, env_names, rotation=10, fontsize=20)
        plt.tick_params(axis='both', which='major', labelsize=15)
        plt.grid(alpha=0.2)
        plt.ylabel(r'Size Ratio (12$\mu m/$optical)',fontsize=20)
        plt.show()
        
        if savefig==True:
            plt.savefig('env_ratio.png', dpi=300)
    
    def mass_hist(self, z0mgs_comp=False, savefig=False):
        
        z0mgs_mass = self.z0mgscut['logmass']
        magphys_mass = self.magphyscut['logMstar']
        
        #remove entries where there is no magphys data available for that galaxy
        err_flag = (self.magphyscut['magphysFlag'])
    
        z0mgs_env_mass = [z0mgs_mass[self.clusflag&err_flag], z0mgs_mass[self.rgflag&err_flag], 
                          z0mgs_mass[self.pgflag&err_flag], 
                          z0mgs_mass[self.filflag&err_flag], z0mgs_mass[self.fieldflag&err_flag]]
        magphys_env_mass = [magphys_mass[self.clusflag&err_flag], magphys_mass[self.rgflag&err_flag], 
                            magphys_mass[self.pgflag&err_flag], 
                            magphys_mass[self.filflag&err_flag], magphys_mass[self.fieldflag&err_flag]]

        env_names = ['cluster','rich group','poor group','filament','field']
        
        #n_bins_mag = [12,12,12,12,12]
        #n_bins_z0 = [12,12,12,12,12]
        mybins=np.linspace(7.5,11.5,12)
        
        fig = plt.figure(figsize=(15,10))
        plt.subplots_adjust(hspace=.4,wspace=.2)
        
        for i in range(1,6):
            ax=fig.add_subplot(2,3,i)
            if i == 4:
                ax.set_position([0.24,0.125,0.228,0.3])
            if i == 5:
                ax.set_position([0.55,0.125,0.228,0.3])
            
            plt.hist(magphys_env_mass[i-1],bins=mybins,alpha=0.2,cumulative=True,density=True,label='Magphys')
            
            if z0mgs_comp==True:
                plt.hist(z0mgs_env_mass[i-1],bins=mybins,alpha=0.2,cumulative=True,density=True,label='z0mgs')
            
            if i == 1:
                plt.legend(fontsize=18)
            plt.xlabel(r'log(M*/$M_\odot$)',fontsize=22)
            plt.title(env_names[i-1],fontsize=22)
            
        plt.show()
        
        if savefig==True:
            plt.savefig('mass_hist.png',dpi=300)
 
    
if __name__ == '__main__':
    cat = catalogs()
    cat.envbins()
    cat.env_means()
    cat.mass_hist(z0mgs_comp=True)