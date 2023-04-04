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
from scipy import stats
from scipy.stats import median_abs_deviation as MAD
from scipy.stats import kstest
from astropy.stats import bootstrap

import os
homedir = os.getenv("HOME")

#bootstrap function courtesy of Rose Finn; 
def get_bootstrap_confint(d,bootfunc=np.median,nboot=100):
    # use astropy.stats bootstrap function
    # it will create nboot resamplings of the data
    # and calculate the bootfunc of each resample.
    # it will return the e.g. median for each of the nboot resamples
    bootsamp = bootstrap(d,bootfunc=bootfunc,bootnum=nboot)

    # sort the bootstrap sampled medians
    bootsamp.sort()

    # get indices corresponding to 68% confidence interval
    ilower = int(((nboot - .68*nboot)/2))
    iupper = nboot-ilower

    # return the e.g. median at the 68% confidence interval
    # need to subtract from median to get the actual errorbars
    # like err_lower = actual_median - bootsamp[ilower]
    # and err_upper = bootsamp[iupper] - actual_median
    return bootsamp[ilower],bootsamp[iupper]

class catalogs:
    
    def __init__(self,conv=False,MeanMedian='mean',MADmultiplier=5):
        self.MADmultiplier = MADmultiplier
        path_to_dir = homedir+'/Desktop/v2-20220820/'
        self.v2_env = Table.read(path_to_dir+'vf_v2_environment.fits')
        self.v2_main = Table.read(homedir+'/v2_snrcoadd.fits')
        self.magphys = Table.read(path_to_dir+'vf_v2_magphys_07-Jul-2022.fits')
        self.z0mgs = Table.read(path_to_dir+'vf_v2_z0mgs.fits')
        
        self.MeanMedian = MeanMedian  #whether I plot median or mean size ratios for the self.env_means() figure
        
        self.conv = conv
        print('Convolution: ',str(self.conv))
        
        if self.conv==False:
            self.rdat = Table.read(homedir+'/output_params_r_nopsf.fits')
            self.w3dat = Table.read(homedir+'/output_params_W3_nopsf.fits')
            self.kimparams = Table.read(homedir+'/kimparams_nopsf.fits')   #my nopsf galfit parameters
        if self.conv==True:
            self.rdat = Table.read(homedir+'/output_params_r_psf.fits')
            self.w3dat = Table.read(homedir+'/output_params_W3_psf.fits')
            self.kimparams = Table.read(homedir+'/kimparams_psf.fits')   #my psf galfit parameters
    
        self.roseparams = self.w3dat.copy()        
        self.cut_cats()
            
    def cut_cats(self):
        subsample_flag = self.v2_main['sgacut_flag']
        
        self.v2_env = self.v2_env[subsample_flag]
        self.v2_maincut = self.v2_main[subsample_flag]
        self.magphys = self.magphys[subsample_flag]
        self.z0mgs = self.z0mgs[subsample_flag]
        
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
        
        self.cut_flags = (~group_flag)&(~err_flag)
        #self.cut_flags = (~group_flag)
        
        self.re_rband_cut = self.re_rband[self.cut_flags]
        self.re_w3band_cut = self.re_w3band[self.cut_flags]

        self.PA_rband_cut = self.PA_rband[self.cut_flags]
        self.PA_w3band_cut = self.PA_w3band[self.cut_flags]

        self.BA_rband_cut = self.BA_rband[self.cut_flags]
        self.BA_w3band_cut = self.BA_w3band[self.cut_flags]
        
        #apply final cut to envcut and maincut catalogs
        self.v2_envcut = self.v2_env[self.cut_flags]
        self.v2_maincut = self.v2_maincut[self.cut_flags]
        self.magphyscut = self.magphys[self.cut_flags]
        self.z0mgscut = self.z0mgs[self.cut_flags]
        
        self.kimparams_cut = self.kimparams[self.cut_flags]
        self.roseparams_cut = self.roseparams[self.cut_flags]
        
        #ratio of Rose's w3 effective radius to my own...both cases are without convolution.
        self.comp_ratios = self.roseparams_cut['re']/self.kimparams_cut['re']
        
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

        if savefig==True:
            plt.savefig(homedir+'/Desktop/envbins.png', dpi=300)
        
        plt.show()

    def env_means(self, trimOutliers=False, combine_mid=False, savefig=False):    
        
        #will generate the self.outlier_flag variable needed to, well, trim the outliers.
        if trimOutliers==True:
            self.compareKim(savefig=False)
            ratios = self.sizerats[self.outlier_flag]
            clusflag = self.clusflag.copy()[self.outlier_flag]
            rgflag = self.rgflag.copy()[self.outlier_flag]
            pgflag = self.pgflag.copy()[self.outlier_flag]
            filflag = self.filflag.copy()[self.outlier_flag]
            fieldflag = self.fieldflag.copy()[self.outlier_flag]
        else:
            ratios = self.sizerats
            clusflag = self.clusflag.copy()
            rgflag = self.rgflag.copy()
            pgflag = self.pgflag.copy()
            filflag = self.filflag.copy()
            fieldflag = self.fieldflag.copy()
        
        re_data = [ratios[clusflag],ratios[rgflag],ratios[pgflag],
                   ratios[filflag],ratios[fieldflag]]
        central_pts = []
        err = []
        index = np.arange(1,6,1)
        env_names = ['Cluster','Rich \n Group','Poor \n Group','Filament','Field']

        for j,i in enumerate(re_data):      #j==index, i==value
            if self.MeanMedian=='mean':
                central_pts.append(np.mean(i))
            if self.MeanMedian=='median':
                central_pts.append(np.median(i))
            err.append(np.std(i)/np.sqrt(len(i)))
        
        err_color = 'orangered'
        plt.figure(figsize=(10,6))
        plt.scatter(index,central_pts,color='blue',s=40,zorder=2,label=self.MeanMedian)
        plt.errorbar(index,central_pts,yerr=err,fmt='None',color=err_color,zorder=1)

        xmin,xmax = plt.xlim()
        xfield = np.linspace(xmin,xmax,50)
        ymax = np.ones(50)*(central_pts[-1] + err[-1])
        ymin = np.ones(50)*(central_pts[-1] - err[-1])
        plt.fill_between(xfield,ymax,ymin,color=err_color,alpha=.1)
        
        plt.xticks(index, env_names, rotation=10, fontsize=20)
        plt.tick_params(axis='both', which='major', labelsize=15)
        plt.grid(alpha=0.2)
        plt.ylabel(r'Size Ratio (12$\mu m/$optical)',fontsize=20)
        
        plt.ylim(0.6,0.96)
        
        plt.legend(fontsize=15)
        
        if savefig==True:
            plt.savefig(homedir+'/Desktop/env_ratio.png', dpi=300)
            
        plt.show()
    
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
        
        fig = plt.figure(figsize=(10,7))
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
            
            if (i == 1)&(z0mgs_comp==True):
                plt.legend(fontsize=18)
            plt.xlabel(r'log(M*/$M_\odot$)',fontsize=22)
            plt.title(env_names[i-1],fontsize=22)
                    
        if savefig==True:
            plt.savefig(homedir+'/Desktop/mass_hist.png',dpi=300)
        
        plt.show()
 
    def compareSGA(self,savefig=False):
        
        r50_sga_r = self.v2_maincut['SGA_r50_arcsec']  #arcsec
        r50_gal_r = self.re_rband_cut.copy()*0.262 #arcsec
        
        plt.figure(figsize=(8,6))
        plt.scatter(r50_sga_r,r50_gal_r,color='orange')
        plt.axline([0, 0], [1, 1], color='black',label='1-to-1')
        
        plt.xlabel('SGA r50 (r-band)',fontsize=18)
        plt.ylabel('GALFIT r50 (r-band)',fontsize=18)
        if self.conv==True:
            plt.title('Re Comparison (PSF)',fontsize=20)
        if self.conv==False:
            plt.title('Re Comparison (noPSF)',fontsize=20)
        
        slope, intercept, r_value, p_value, std_err = stats.linregress(r50_sga_r,r50_gal_r)
        # Create empty plot with blank marker containing the extra label
        plt.plot([], [], ' ', label=f'Slope (linregress) {round(slope,3)}')
                
        plt.legend(fontsize=16)
                
        plt.xscale('log')
        plt.yscale('log')
                
        if savefig==True:
            plt.savefig(homedir+'/Desktop/SGA_r50_comparison.png',dpi=300)
        
        plt.show()
    
    def comparePSF(self,savefig=False):
        
        if self.conv==False:
            self.rdat_comp = Table.read(homedir+'/output_params_r_psf.fits')
            self.w3dat_comp = Table.read(homedir+'/output_params_W3_psf.fits')
            xlabels = ['psf Re (px)','psf Re (px)']
            ylabels = ['nopsf Re (px)','nopsf Re (px)']
        if self.conv==True:
            self.rdat_comp = Table.read(homedir+'/output_params_r_nopsf.fits')
            self.w3dat_comp = Table.read(homedir+'/output_params_W3_nopsf.fits')
            xlabels = ['nopsf Re (px)','nopsf Re (px)']
            ylabels = ['psf Re (px)','psf Re (px)']
        
        titles = ['w3 Re Comparison','r-band Re Comparison']
                        
        rdat_comp = self.rdat_comp[self.cut_flags]
        w3dat_comp = self.w3dat_comp[self.cut_flags]
                
        err_flag_comp_r = (rdat_comp['err_flag']==1)
        err_flag_comp_w3 = (w3dat_comp['err_flag']==1)
        err_flag_comp = err_flag_comp_r | err_flag_comp_w3
        
        re_rband_comp = rdat_comp['re'][~err_flag_comp]
        re_w3band_comp = w3dat_comp['re'][~err_flag_comp]
        re_rband = self.re_rband_cut.copy()[~err_flag_comp]
        re_w3band = self.re_w3band_cut.copy()[~err_flag_comp]
        
        re_data_init = [re_w3band*2.75,re_rband*0.262]
        re_data_comp = [re_w3band_comp*2.75,re_rband_comp*0.262]
                
        plt.figure(figsize=(14,6))
        for panel in range(2):
            ax = plt.subplot(1,2,panel+1)
            plt.scatter(re_data_comp[panel], re_data_init[panel], color='crimson', alpha=0.5)
            plt.title(titles[panel],fontsize=16)
            ax.set_xlabel(xlabels[panel],fontsize=16)
            ax.set_ylabel(ylabels[panel],fontsize=16)
            ax.axline((0, 0), slope=1, label='1-to-1',color='black')
            
            slope, intercept, r_value, p_value, std_err = stats.linregress(re_data_comp[panel], re_data_init[panel])
            # Create empty plot with blank marker containing the extra label
            plt.plot([], [], ' ', label=f'Slope (linregress) {round(slope,3)}')
            
            plt.xscale('log')
            plt.yscale('log')
            plt.legend(fontsize=14)
 
        if savefig==True:
            plt.savefig(homedir+'/Desktop/PSF_r50_comparison.png',dpi=300)
        
        plt.show()
    
    def compareKim(self, savefig=False):
        
        lower_bound = np.median(self.comp_ratios) - MAD(self.comp_ratios)*int(self.MADmultiplier)
        upper_bound = np.median(self.comp_ratios) + MAD(self.comp_ratios)*int(self.MADmultiplier)
        self.outlier_flag = (self.comp_ratios > lower_bound)&(self.comp_ratios < upper_bound)
        
        plt.figure(figsize=(8,6))
        plt.scatter(self.kimparams_cut['re'][self.outlier_flag],
                    self.roseparams_cut['re'][self.outlier_flag],alpha=0.2,color='crimson',
                    label='Points within +/- MAD*{}'.format(self.MADmultiplier))
        plt.scatter(self.kimparams_cut['re'][~self.outlier_flag],
                    self.roseparams_cut['re'][~self.outlier_flag],alpha=0.2,color='black',s=20,label='All Points')

        xline=np.linspace(np.min(self.kimparams_cut['re']),np.max(self.kimparams_cut['re']),100)
        plt.plot(xline,xline,color='black')
        
        plt.xscale('log')
        plt.yscale('log')
        plt.legend(fontsize=14,loc='upper left')
        
        plt.xlabel('Kim w3 (px)',fontsize=18)
        plt.ylabel('Rose w3 (px) ',fontsize=18)
        if self.conv==True:
            plt.title('Re Comparison (PSF)',fontsize=20)
        if self.conv==False:
            plt.title('Re Comparison (noPSF)',fontsize=20)
        
        print('List of Galaxy Outliers')
        print(self.kimparams_cut['VFID'][~self.outlier_flag])
        
        if savefig==True:
            plt.savefig(homedir+'/Desktop/Re_comparison_Kim.png',dpi=300)
        
        plt.show()

    def recreate_LCS_hist(self, keep_errs=False, savefig=False):
        re_r = self.re_rband_cut.copy()
        re_w3 = self.re_w3band_cut.copy()
        clusflag = self.clusflag.copy()
        if keep_errs:
            re_r = self.re_rband[self.re_w3band>0]
            re_w3 = self.re_w3band[self.re_w3band>0]
            clusflag = self.v2_env['cluster_member'][self.re_w3band>0]
        
        re_r = re_r*0.262
        re_w3 = re_w3*2.75
        
        data_clus = re_w3[clusflag]/re_r[clusflag]
        data_ext = re_w3[~clusflag]/re_r[~clusflag]
        data = [data_clus, data_ext]
        labels = ['Cluster Galaxies', 'External Galaxies']
        titles = ['Size Ratio Distribution', '']
        colors = ['crimson', 'blue']
        xlabels = ['',r'R$_{12}$/R$_r$']
        
        mybins=np.linspace(0,3,40)
        
        fig = plt.figure(figsize=(8,8))
        plt.subplots_adjust(hspace=.1)
        
        for panel in range(2):
            ax = fig.add_subplot(2,1,panel+1)
            plt.hist(data[panel],bins=mybins,color=colors[panel], alpha=0.7, label=labels[panel],density=False,cumulative=False)
            plt.hist([], color='white', label='0 < R12/Rr < 2')
            plt.title(titles[panel],fontsize=16)
            ax.set_xlabel(xlabels[panel],fontsize=20)
            ax.set_ylabel(r'N$_{gal}$',fontsize=20)
            plt.xlim(-0.1,2)
            ax.legend(fontsize=15)
        
        if savefig==True:
            plt.savefig(homedir+'/Desktop/LCS_hists.png',dpi=300)

        plt.show()
        
        #ks-test statistics
        print('p-value (> 0.05, "same distribution"):')
        print('%.3f'%(kstest(data_clus[data_clus<=2],data_ext[data_ext<=2])[1]))
        print()
        print('Mean cluster size ratio: %.3f'%(np.mean(data_clus[data_clus<=2])))
        print('Median cluster size ratio: %.3f'%(np.median(data_clus[data_clus<=2])))
        print('Mean external size ratio: %.3f'%(np.mean(data_ext[data_ext<=2])))
        print('Median external size ratio: %.3f'%(np.median(data_ext[data_ext<=2])))
        print()
        print('Cluster mean/median uncertainty: %.3f'%(np.std(data_clus[data_clus<=2])/np.sqrt(len(data_clus[data_clus<=2]))))
        print('External mean/median uncertainty: %.3f'%(np.std(data_ext[data_ext<=2])/np.sqrt(len(data_ext[data_ext<=2]))))
     
    def recreate_LCS_mass(self, keep_errs=False, savefig=False):
        re_r = self.re_rband_cut.copy()
        re_w3 = self.re_w3band_cut.copy()
        clusflag = self.clusflag.copy()
        logmass = self.z0mgscut['logmass'].copy()
        if keep_errs:
            re_r = self.re_rband[self.re_w3band>0]
            re_w3 = self.re_w3band[self.re_w3band>0]
            clusflag = self.v2_env['cluster_member'].copy()[self.re_w3band>0]
            logmass = self.z0mgs['logmass'].copy()[self.re_w3band>0]
            
        re_r = re_r*0.262
        re_w3 = re_w3*2.75
        re_ratio = re_w3/re_r

        #define mass bins --> min<log(Mstar)<max, 6 bins
        mass_bounds = np.linspace(np.min(logmass),np.max(logmass),7)
        bin_bounds = []
        for bound in range(6):
            bin_bounds.append([mass_bounds[bound],mass_bounds[bound+1]])
        
        #now for the tedious task of creating flags which will separate the ratios into these mass bins :}
        re_mass_clus=[]
        re_mass_ext=[]
        mass_coord=[]
        err_clus=[]
        err_ext=[]
        
        #note --> bootstrap errors correspond to 68% confidence interval
        for bound in bin_bounds: 
            bound_flag = (logmass>bound[0])&(logmass<=bound[1])
            if self.MeanMedian=='mean':
                avg_re_clus = np.mean(re_ratio[bound_flag & clusflag])
                avg_re_ext = np.mean(re_ratio[bound_flag & ~clusflag])
                lower_clus, upper_clus = get_bootstrap_confint(re_ratio[bound_flag & clusflag],bootfunc=np.mean,nboot=100)
                lower_ext, upper_ext = get_bootstrap_confint(re_ratio[bound_flag & ~clusflag],bootfunc=np.mean,nboot=100)
            if self.MeanMedian=='median':
                avg_re_clus = np.median(re_ratio[bound_flag & clusflag])
                avg_re_ext = np.median(re_ratio[bound_flag & ~clusflag])
                lower_clus, upper_clus = get_bootstrap_confint(re_ratio[bound_flag & clusflag],bootfunc=np.median,nboot=100)
                lower_ext, upper_ext = get_bootstrap_confint(re_ratio[bound_flag & ~clusflag],bootfunc=np.median,nboot=100)
            
            mass_coord.append(np.mean(logmass[bound_flag]))
            re_mass_clus.append(avg_re_clus)
            re_mass_ext.append(avg_re_ext)
            
            #err_clus.append([avg_re_clus-lower_clus, upper_clus-avg_re_clus])
            #err_ext.append([avg_re_ext-lower_ext, upper_ext-avg_re_ext])
            err_clus.append([lower_clus, avg_re_clus])
            err_ext.append([lower_ext, upper_ext])
        
        plt.figure(figsize=(8,6))
        
        plt.scatter(logmass[clusflag],re_ratio[clusflag],color='crimson',s=15,alpha=0.1,label='Cluster',zorder=1)
        plt.scatter(logmass[~clusflag],re_ratio[~clusflag],color='blue',s=15,alpha=0.1,label='External',zorder=1)
        plt.scatter(mass_coord,re_mass_clus,color='crimson',s=250,edgecolors='black',label='<Cluster>',zorder=3)
        plt.scatter(mass_coord,re_mass_ext,color='blue',s=250,edgecolors='black',label='<External>',zorder=3)
        
        for n in range(6):
            plt.plot([mass_coord[n],mass_coord[n]], [err_clus[n][0],err_clus[n][1]],color='crimson',zorder=2)
            #create lower, upper caps on errorbars
            plt.plot([mass_coord[n]-0.08,mass_coord[n]+0.08],[err_clus[n][0],err_clus[n][0]],color='crimson',zorder=2)  
            plt.plot([mass_coord[n]-0.08,mass_coord[n]+0.08],[err_clus[n][1],err_clus[n][1]],color='crimson',zorder=2)
            
            plt.plot([mass_coord[n],mass_coord[n]], [err_ext[n][0],err_ext[n][1]],color='blue',zorder=2)
            #create lower, upper caps on errorbars
            plt.plot([mass_coord[n]-0.08,mass_coord[n]+0.08],[err_ext[n][0],err_ext[n][0]],color='blue',zorder=2)
            plt.plot([mass_coord[n]-0.08,mass_coord[n]+0.08],[err_ext[n][1],err_ext[n][1]],color='blue',zorder=2)
            
        plt.ylabel(r'R$_{12}$/R$_r$',fontsize=18)
        plt.xlabel(r'log$_{10}$(M$_*$/M$_\odot$)',fontsize=18)
        plt.title('Median Disk Size Ratios vs. Stellar Masses',fontsize=20)
        plt.ylim(-0.1,2.5)
        plt.legend()
        
        if savefig==True:
            plt.savefig(homedir+'/Desktop/LCS_mass.png',dpi=300)

        plt.show()
    
    
if __name__ == '__main__':
    print("""USAGE:
    cat = catalogs(conv=False,MeanMedian='mean',MADmultiplier=5) --> initiate catalog class. 
        MeanMedian will propagate to all plots.
    cat.envbins(savefig=False) --> plots number of subsample galaxies in each environment bin
    cat.env_means(trimOutliers=False, combine_mid=False, savefig=False) --> plots either mean 
        or median size ratio (w3/r) in each environment bin; trimOutliers will output an additional plot which compares
        my no PSF parameters to Rose's parameters, allowing the user to visualize which points are omitted in the 
        trimmed env_means plot; combine_mid is currently in development, but the aim is to merge the three middle
        environment values in order to, in part, help reduce error bars.
    cat.mass_hist(z0mgs_comp=True,savefig=False) --> generate mass histogram subplots per environment bin; 
        will compare MAGPHYS stellar masses with z0mgs values if True
    cat.compareSGA(savefig=False) --> compares Rose's GALFIT r-band Re values with SGA's non-parametric r50
    cat.comparePSF(savefig=False) --> plots noPSF vs. PSF Re values for w3, r-band (one subplot per band)
    cat.compareKim(savefig=False) --> compares my noPSF w3-band Re values with Rose's noPSF values
    cat.recreate_LCS_hist(keep_errs=False,savefig=False) --> generates vertically-oriented histogram subplots of R12/Rr 
        distribution, separated into cluster vs. all else (external). keep_errs governs whether the final counts include 
        galaxies with GALFIT error flags (True if yes, False if no).
    cat.recreate_LCS_mass(keep_errs=False,savefig=False) --> generates scatterplot of size ratio vs. mass bin, 
        the format being similar to a skeleton version of Figure 13 from Finn+18. I use z0mgs masses here 
        (odd results with MAGPHYS for min/max, also many galaxies with nan or masked values).
    """)
    print('-----------------------------------------------------')
    print()