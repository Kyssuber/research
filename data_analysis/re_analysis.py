'''
GOAL: analyze GALFIT output parameters for Virgowise subsample galaxies.

pixel_scale = {'FUV':1.5,'NUV':1.5,'g':0.262,'r':0.262,'z':0.262,'W1':2.75,'W2':2.75,'W3':2.75,'W4':2.75}   arcsec/px
psf_oversampling = {'FUV':1,'NUV':1,'g':1,'r':1,'g':1,'W1':1,'W2':1,'W3':1,'W4':1}
mag_zeropoint = {'FUV':22.5,'NUV':22.5,'g':22.5,'r':22.5,'g':22.5,'W1':22.5,'W2':22.5,'W3':22.5,'W4':22.5}
image_resolution = {'FUV':6,'NUV':6,'g':1.5,'r':1.5,'z':1.5,'W1':6.1,'W2':6.4,'W3':6.5,'W4':12}   arcsec
'''

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams
from astropy.table import Table
from scipy.stats import median_abs_deviation as MAD
from astropy.stats import bootstrap
from itertools import combinations
from scipy.stats import kstest
from scipy.stats import binned_statistic
from scipy.stats import spearmanr
from scipy.stats import linregress
from scipy.stats import ttest_1samp
from scipy.stats import mannwhitneyu
from scipy.stats import anderson_ksamp
from scipy.stats import PermutationMethod
from scipy.special import gamma, gammainc, gammaincinv
from matplotlib import ticker

import os
homedir = os.getenv("HOME")
path_to_dir = homedir+'/Desktop/v2-20220820/'

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

#need to calculate L50
def L(r,re,n):
    I_e = 1     #just a brightness normalization term; set equal to 1.
    b_n = gammaincinv(2*n, 0.5)
    r_ratio = r/re
    constant = 2*np.pi*n*(re**2)/b_n**(2*n)
    gamma_term = gammainc(2*n,b_n*(r_ratio**(1/n))) #note that scipy.special.gammainc is NORMALIZED. 
                                                    #need to multiple gamma_term * gamma(2*n) in cases
                                                    #where I'm not calculating an L ratio or I care
                                                    #about units, etc. 
    return constant*I_e*np.exp(b_n)*gamma_term

#need to calculate r90 using L50, so isolate r from L equation above...as Dr. Tosun would exclaim, "this is a nice problem." It was not a nice problem.
def r(L,re,n):
    I_e = 1
    b_n = gammaincinv(2*n,0.5)          #also approximated by 1.999 * n-0.327
    const_numer = 2*np.pi*n*(re**2)*I_e*np.exp(b_n)
    const_denom = b_n**(2*n)
    constant = const_numer/const_denom
    frac = L/constant     #oddly enough, this fraction is 0.90 for L90 and 0.50 for L50...I guess constant = Lmax
    gamma_term = gammaincinv(2*n, L/constant)   #gamma(2*n) corrects the normalization
    
    exp_term = (gamma_term/b_n)**n
    
    return re*exp_term

def calculate_r90(re, n):
    return r(1.8*L(re,re,n),re,n)

class catalogs:
    
    def __init__(self,conv=False,MeanMedian='mean',MADmultiplier=5,cutAGN=False,W1=False):
        self.MADmultiplier = MADmultiplier
        self.v2_env = Table.read(path_to_dir+'vf_v2_environment.fits')
        self.v2_main = Table.read(homedir+'/Desktop/galfit_files/VF_WISESIZE_v2.fits')  #has my flagshh
        self.magphys = Table.read(path_to_dir+'/vf_v2_magphys_legacyExt_mergedNS.fits')
        self.z0mgs = Table.read(path_to_dir+'vf_v2_z0mgs.fits')
        self.HI_tab = Table.read(path_to_dir+'vf_v2_CO_HI.fits')
        self.hyp_tab = Table.read(path_to_dir+'vf_v2_hyperleda.fits')
        self.altmagphys = Table.read(homedir+'/Desktop/galfit_files/vf-altphot.fits')
        #self.sgaparams = Table.read(homedir+'/sgacut_SGAparams.fits')
        #self.sgaparams.sort('VFID_1')   #sort in ascending order all rows according to VFID
        self.MeanMedian = MeanMedian  #whether I plot median or mean size ratios for the self.env_means() figure
        
        self.conv = conv
        self.cutAGN = cutAGN
        self.W1 = W1
        print('Convolution: ',str(self.conv))
        if self.cutAGN:
            print('AGN cut from sample.')
        
        
        self.rdat = Table.read(homedir+'/Desktop/galfit_files/galfit_r_03112024.fits')
        #self.w3dat = Table.read(homedir+'/Desktop/galfit_files/galfit_W3_03112024.fits')
        #self.w1dat = Table.read(homedir+'/Desktop/galfit_files/galfit_W1_03112024.fits')
        
        self.w3dat = Table.read(homedir+'/Desktop/galfit_files/vf_v2_galfit_W3-fixBA.fits')
        self.w1dat = Table.read(homedir+'/Desktop/galfit_files/vf_v2_galfit_W1-fixBA.fits')
        
        
        self.cut_cats()
    
    def cut_cats(self):
        
        subsample_flag=self.v2_main['subsample_flag']
        
        self.v2_env = self.v2_env[subsample_flag]
        self.v2_maincut = self.v2_main[subsample_flag]
        self.magphyscut = self.magphys[subsample_flag]
        self.z0mgs = self.z0mgs[subsample_flag]
        self.HI_tab = self.HI_tab[subsample_flag]
        self.hyp_tab = self.hyp_tab[subsample_flag]
        self.altmagphys_cut = self.altmagphys[subsample_flag]
                
        if self.W1:
            self.re_w1band = self.w1dat['CRE'][subsample_flag]
            self.n_w1band = self.w1dat['CN'][subsample_flag]
        
        self.n_rband = self.rdat['CN'][subsample_flag]
        self.re_rband = self.rdat['CRE'][subsample_flag]
        
        self.n_w3band = self.w3dat['CN'][subsample_flag]
        self.re_w3band = self.w3dat['CRE'][subsample_flag]
        
        self.PA_rband = self.rdat['CPA'][subsample_flag]
        self.PA_w3band = self.w3dat['CPA'][subsample_flag]

        self.BA_rband = self.rdat['CAR'][subsample_flag]
        self.BA_w3band = self.w3dat['CAR'][subsample_flag]

        fail_flag = (self.rdat['CXC'][subsample_flag] == 0.0) | (self.w3dat['CXC'][subsample_flag] == 0.0)
        r_flag = (self.rdat['CXC'][subsample_flag] == 0.0)
        w3_flag = (self.w3dat['CXC'][subsample_flag] == 0.0)
        err_flag = (self.rdat['CNumerical_Error'][subsample_flag]) | (self.w3dat['CNumerical_Error'][subsample_flag])
        sersic_flag = (self.rdat['CN'][subsample_flag]>6) | (self.w3dat['CN'][subsample_flag]>6) #unphysical
        
        nsersic_fails = len(self.re_w3band[sersic_flag])

        n_fails_r = len(self.re_rband[r_flag])
        n_fails_w3 = len(self.re_w3band[w3_flag])
        
        if self.W1: 
            fail_flag = (self.w1dat['CXC'][subsample_flag] == 0.0) | (self.w3dat['CXC'][subsample_flag] == 0.0)  #no entries; galfit failed.
            w1_flag = (self.w1dat['CXC'][subsample_flag] == 0)  #no w1 entries; galfit failed.
            err_flag = (self.w1dat['CNumerical_Error'][subsample_flag]) | (self.w3dat['CNumerical_Error'][subsample_flag]) | (self.w1dat['CXC'][subsample_flag]==0.0) | (self.w3dat['CXC'][subsample_flag]==0.0) #ALL ERRORS
            
            sersic_flag = (self.w1dat['CN'][subsample_flag]>6) | (self.w3dat['CN'][subsample_flag]>6) #unphysical 
            
            n_fails_w1 = len(self.re_w1band[w1_flag])
            nsersic_fails = len(self.re_w3band[sersic_flag])
            
        self.cut_flags = (~fail_flag)&(~err_flag)&(~sersic_flag)   #only keep galaxies that are NOT errors or failures. only winners in this house.
        
        #apply optional AGN cut
        if self.cutAGN:

            AGN_flags = (self.v2_maincut['WISE_AGN_flag'])|(self.v2_maincut['kauffman_AGN_flag'])
            
            print(f'fraction AGN in VF subsample (before filtering out GALFIT errors): {np.round(len(self.v2_maincut[AGN_flags])/len(self.v2_maincut),3)}')
            print(f'This fracion is calculated from {len(self.v2_maincut[AGN_flags])}/{len(self.v2_maincut)}')

            self.cut_flags = (~AGN_flags) & (~fail_flag) & (~err_flag) & (~sersic_flag)
            print(f'Number of galaxies flagged with AGN *AND* a GALFIT error: {len(self.v2_env[(AGN_flags) & (err_flag)])}')
            print()
        
        self.re_rband_cut = self.re_rband[self.cut_flags]
        self.re_w3band_cut = self.re_w3band[self.cut_flags]
        if self.W1:
            self.re_w1band_cut = self.re_w1band[self.cut_flags]
            self.n_w1band_cut = self.n_w1band[self.cut_flags]

        self.n_rband_cut = self.n_rband[self.cut_flags]
        self.n_w3band_cut = self.n_w3band[self.cut_flags]
            
        self.PA_rband_cut = self.PA_rband[self.cut_flags]
        self.PA_w3band_cut = self.PA_w3band[self.cut_flags]

        self.BA_rband_cut = self.BA_rband[self.cut_flags]
        self.BA_w3band_cut = self.BA_w3band[self.cut_flags]
        
        #apply final cut to envcut and maincut catalogs
        self.v2_envcut = self.v2_env[self.cut_flags]
        self.v2_maincut = self.v2_maincut[self.cut_flags]
        self.magphyscut = self.magphyscut[self.cut_flags]
        self.z0mgscut = self.z0mgs[self.cut_flags]
        self.HI_tab_cut = self.HI_tab[self.cut_flags]
        self.hyp_tab_cut = self.hyp_tab[self.cut_flags]
        self.altmagphys_cut = self.altmagphys_cut[self.cut_flags]
        
        #define env flags
        self.clusflag = self.v2_envcut['cluster_member']
        self.rgflag = self.v2_envcut['rich_group_memb']
        self.pgflag = self.v2_envcut['poor_group_memb']
        self.filflag = self.v2_envcut['filament_member']
        self.fieldflag = self.v2_envcut['pure_field']
        
        if self.W1:
            print(f'No GALFIT data for {n_fails_w3} w3 galaxies, {n_fails_w1} w1 galaxies, and {n_fails_r} r galaxies.')
        if not self.W1:
            print(f'No GALFIT data for {n_fails_w3} w3 galaxies and {n_fails_r} r galaxies.')
        print(f'Total number of galaxies with GALFIT errors or error flags: {int(np.sum(np.ones(len(err_flag))*err_flag))}')
        print(f'Total number of galaxies with nser>6: {nsersic_fails}')
        print()
        print(f'Total number of subsample galaxies remaining: {len(self.re_w3band_cut)}')

        if self.W1:
            self.sizerats = (self.re_w3band_cut*2.75)/(self.re_w1band_cut*2.75)
        else:
            self.sizerats = (self.re_w3band_cut*2.75)/(self.re_rband_cut*0.262)
            self.PArats = self.PA_w3band_cut/self.PA_rband_cut 
        
    def mass_matching(self):
        
        #select one cluster galaxy
        #create pool of field galaxies with logmass within 0.3 dex relative to the one cluster galaxy
        #from this pool, randomly draw int(len(field)/len(cluster)) galaxies
            #try without replacement first.
        #repeat (2) and (3) for rich group, filament, and poor group galaxies
        #repeat all for every cluster galaxy
        
        n_dex = 0.3
        
        n_to_draw = [int(round(len(self.v2_envcut[self.fieldflag])/len(self.v2_envcut[self.clusflag]),0)),
                     int(round(len(self.v2_envcut[self.filflag])/len(self.v2_envcut[self.clusflag]),0)),
                     int(round(len(self.v2_envcut[self.pgflag])/len(self.v2_envcut[self.clusflag]),0)),
                     int(round(len(self.v2_envcut[self.rgflag])/len(self.v2_envcut[self.clusflag]),0))]
        
        #remove entries where there is no magphys data available for that galaxy    
        err_flag_cut = self.magphyscut['magphysFlag']
        
        #isolate the rows of galaxies in each environment type
        m_clus = self.magphyscut[(err_flag_cut) & (self.clusflag)]
        
        m_rg = self.magphyscut[(err_flag_cut) & (self.rgflag)]
        m_pg = self.magphyscut[(err_flag_cut) & (self.pgflag)]
        m_fil = self.magphyscut[(err_flag_cut) & (self.filflag)]
        m_field = self.magphyscut[(err_flag_cut) & (self.fieldflag)]
        
        #create empty lists in which the mass-matched environment galaxies will go
        comp_galaxies = [[],[],[],[]]
        
        for index in range(len(m_clus)):
            
            #for every cluster galaxy, randomly choose two galaxies from each environment with masses within 0.3 dex of the cluster galaxy mass
            field_pool = m_field[(m_field['logMstar_med'] < (m_clus['logMstar_med'][index]+n_dex)) & (m_field['logMstar_med'] > (m_clus['logMstar_med'][index]-n_dex))]
            rg_pool = m_rg[(m_rg['logMstar_med'] < (m_clus['logMstar_med'][index]+n_dex)) & (m_rg['logMstar_med'] > (m_clus['logMstar_med'][index]-n_dex))]
            pg_pool = m_pg[(m_pg['logMstar_med'] < (m_clus['logMstar_med'][index]+n_dex)) & (m_pg['logMstar_med'] > (m_clus['logMstar_med'][index]-n_dex))]
            fil_pool = m_fil[(m_fil['logMstar_med'] < (m_clus['logMstar_med'][index]+n_dex)) & (m_fil['logMstar_med'] > (m_clus['logMstar_med'][index]-n_dex))]
            
            #the 'pools' of galaxies from which the mass-matched samples will be drawn
            pools = [field_pool, fil_pool, pg_pool, rg_pool]
            
            for env in range(4):
                if n_to_draw[env] <= len(pools[env]):
                    choices = np.random.choice(pools[env]['VFID'],n_to_draw[env])
                else:    
                    try:
                        choices = np.random.choice(pools[env]['VFID'],n_to_draw[env]-1)
                    except ValueError:
                        choices = np.random.choice(pools[env]['VFID'],0)
                for gal in choices:
                    comp_galaxies[env].append(gal.astype(str))
            
        #print('n_field, n_rg, n_pg, n_filament:', len(comp_galaxies[0]), len(comp_galaxies[1]), 
        #      len(comp_galaxies[2]), len(comp_galaxies[3]))
        #print('n_clus:',len(m_clus))

        self.comp_galaxies = comp_galaxies
        
        #for every VFID in each comp_galaxies array, find the index at which this VFID appears in v2main_cut and
        #extract the corresponding size ratio in self.sizerats. do likewise for masses and SFRs. 
        self.matched_sizes = [[],[],[],[]]
        self.matched_Mstar = [[],[],[],[]]
        self.matched_SFR = [[],[],[],[]]
        for env in range(4):
            for VFID in self.comp_galaxies[env]:
                index_sizes = np.where(VFID==self.v2_maincut['VFID'])[0]
                index_sizes = index_sizes[0]
                self.matched_sizes[env].append(self.sizerats[index_sizes])
                self.matched_SFR[env].append(self.magphyscut['logSFR_med'][index_sizes])
                self.matched_Mstar[env].append(self.magphyscut['logMstar_med'][index_sizes])
                
    def inclination_check(self, savefig=False):
        
        #aim: check dependence of size ratio, SFR, Mstar on optical axis ratio
        #also looks for whether trends exist between size ratio and Sérsic index or T-type
        print('Spearman rank test: -1<p<1, with 0 indicating no correlation.')
        
        #remove entries where there is no magphys data available for that galaxy        
        err_flag_cut = (self.magphyscut['magphysFlag'])
        
        #cut logsfr, logmass
        logsfr_cut = self.magphyscut['logSFR_med'][err_flag_cut]
        logmass_cut = self.magphyscut['logMstar_med'][err_flag_cut]
        
        #cut size ratio, optical axis ratio arrays
        sizerats = self.sizerats[err_flag_cut]
        BA = self.BA_rband_cut[err_flag_cut]
        
        x = [BA, BA[(logsfr_cut>-6)], BA[(logmass_cut>6)]]
        y = [sizerats, logsfr_cut[(logsfr_cut>-6)], logmass_cut[(logmass_cut>6)]]
        colors = ['crimson','green','blue']
        y_labels = [r'R$_{12}$/R$_r$', 'logSFR_med', r'log$M_*$']
        png_titles = ['Sizerat vs BA', 'SFR vs BA', 'Mstar vs BA']
        
        fig=plt.figure(figsize=(20,4))
        plt.subplots_adjust(hspace=.1,wspace=.2)
        
        for index in range(len(y)):
            ax = fig.add_subplot(1, len(y), index+1)
            plt.scatter(x[index], y[index], color=colors[index], s=20)
            plt.xlabel('Axis Ratio',fontsize=18)
            plt.ylabel(y_labels[index],fontsize=18)
            
            #spearman-rank test
            sr = spearmanr(x[index],y[index])
            rho = sr[0]
            p = sr[1]
            plt.plot([],[],color="None",label=f'rho={np.round(rho,4)}')
            plt.plot([],[],color="None",label=f'p-value={np.round(p,4)}')
            
            plt.legend(handlelength=0,fontsize=17)
            plt.title(png_titles[index],fontsize=15)
            
            if index==0:
                plt.yscale('log')
            
        if savefig:
            plt.savefig(homedir+'/Desktop/inclination_check.png', dpi=100, bbox_inches='tight', pad_inches=0.2)

        plt.show()
        
        vfid_list = self.v2_maincut['VFID'][err_flag_cut]
        print(f'Size ratio outliers: {vfid_list.data[(sizerats>2)]}, \n logSFR outliers: {vfid_list.data[(logsfr_cut<-6)]}, \n logMstar outliers: {vfid_list.data[(logmass_cut<6)]}')
        
        #also looks for whether trends exist between size ratio and Sérsic index or T-type
        #ttype_dict = {'Sa':1,'Sab':2,'Sb':3,'Sbc':4,'Sc':5,'Scd':6,'Sd':7,'Sdm':8,'Sm':9,'Ir':10}
        
        nser = self.rdat['CN'][self.cut_flags]
        
        y2 = [self.sizerats, self.sizerats[(nser<4)]]
        x2 = [self.hyp_tab_cut['t'], nser[(nser<4)]]
        x_labels2 = ['Morphological T-Type', 'Sérsic Index (r-band)']
        colors2 = ['orangered','blueviolet']
        
        fig=plt.figure(figsize=(20,4))
        plt.subplots_adjust(hspace=.1,wspace=.2)
        
        for index in range(2):
            ax = fig.add_subplot(1, 2, index+1)
            plt.scatter(x2[index][y2[index]<10], y2[index][y2[index]<10], color=colors2[index], s=20)
            plt.xlabel(x_labels2[index],fontsize=18)
            plt.ylabel(r'R$_{12}$/R$_r$',fontsize=18)
            
            #spearman-rank test
            sr = spearmanr(x2[index],y2[index])
            rho = sr[0]
            p = sr[1]
            plt.plot([],[],color="None",label=f'rho={np.round(rho,4)}')
            plt.plot([],[],color="None",label=f'p-value={np.round(p,4)}')
            
            plt.yscale('log')
            plt.legend(handlelength=0,fontsize=17)
        
        if savefig:
            plt.savefig(homedir+'/Desktop/morphology_check.png', dpi=100, bbox_inches='tight', pad_inches=0.2)

        plt.show()

    def mass_completeness_check(self, mass_match=False, savefig=False):
        
        #aim: Make histogram of full VFS sample with quality flags. Also make histogram with T-type and 12um SNR cut. Take the ratio of these histograms (WISESize sample/VFS) and plot them vs. mass. Below some mass this ratio will decline. We can use that point as our mass completeness limit.
        
        mass_bins = np.arange(6,11.5,0.15)
        
        logmass = self.magphys['logMstar_med']
        
        #remove entries where there is no magphys data available for that galaxy
        err_flag = (self.magphys['magphysFlag'])
        
        logmass_cut = logmass[self.v2_main['SNRflag'] & self.v2_main['t_flag'] & (err_flag)]        
        logmass = logmass[err_flag]
        if mass_match:
            self.mass_matching()
            logmass_cut = np.concatenate([self.matched_Mstar[0], self.matched_Mstar[1], self.matched_Mstar[2], 
                                     self.matched_Mstar[3]])
        
        counts_full = plt.hist(logmass,bins=mass_bins,range=[6, 11.5],alpha=0.5,label='VFS')
        counts_cut = plt.hist(logmass_cut,bins=mass_bins,range=[6, 11.5],alpha=0.6,label='WISESize')
        plt.title('Mass Histograms',fontsize=14)
        plt.legend(handlelength=0.3,fontsize=14)
        counts_ratio = counts_cut[0]/counts_full[0]
        
        xedges = counts_full[1]
        xcenters = (xedges[:-1] + xedges[1:]) / 2
                
        plt.figure(figsize=(10,6))
        plt.scatter(xcenters,counts_ratio)    
        plt.xlabel(r'logM$_*$',fontsize=18)
        plt.ylabel('WISE / VFS Ratio',fontsize=18)
        plt.grid(alpha=0.1)
        
        if savefig:
            plt.savefig(homedir+'/Desktop/mass_completeness_hist.png', dpi=100, bbox_inches='tight', 
                        pad_inches=0.2)
    
    def sfrmstar_magphys(self, show_HI=False, show_D25=False, show_sizerat=True, savefig=False):

        sizerats_all = (self.w3dat['CRE']*2.75)/(self.w1dat['CRE']*2.75)
        
        subsample_cut = (self.v2_main['SNRflag']) & (self.v2_main['t_flag'])
        
        logsfr = self.altmagphys['combined_logSFR_med']
        logmass = self.altmagphys['combined_logMstar_med']

        #remove entries where there is no magphys data available for that galaxy. this includes values == 0.
        ######boss also requested removal of SFR==0 galaxies, which will eliminate that horizontal stripe.
        err_flag = (self.altmagphys['magphys_flag']) &  (self.altmagphys['combined_logSFR_med']!=0.)
        err_flag_cut = (self.altmagphys_cut['magphys_flag']) & (self.altmagphys_cut['combined_logSFR_med']!=0.)
     
        #if self.cutAGN:
        #    AGN_flags = (self.v2_main['WISE_AGN_flag']) | (self.v2_main['kauffman_AGN_flag'])
        #    err_flag = (self.altmagphys['magphys_flag']) & (self.altmagphys['combined_logSFR_med']!=0.) & (~AGN_flags)
    
        logmass_cut = logmass[subsample_cut & err_flag]  
        logsfr_cut = logsfr[subsample_cut & err_flag]

        logsfr = logsfr[err_flag]
        logmass = logmass[err_flag]

        plt.figure(figsize=(10,6))
        plt.scatter(logmass,logsfr,color='gray',s=3,alpha=0.05)
        
#############################################################
        if show_HI:
            HI_tab = Table.read(path_to_dir+'vf_v2_CO_HI.fits')
        
            #prepare MHI_to_Mstar data for colorbar
            MHI_to_Mstar = HI_tab['MHI']/(10**self.magphys['logMstar_med'])

            for n in range(len(MHI_to_Mstar)):
                if HI_tab['MHI'].mask[n]:     #if value is masked, set to be -999
                    MHI_to_Mstar[n] = -999
                else:
                    MHI_to_Mstar[n] = MHI_to_Mstar[n] if ((np.log(HI_tab['MHI'][n])>0) & (MHI_to_Mstar[n]<10) & (MHI_to_Mstar[n]!=1)) else -999
        
            MHI_to_Mstar_cut = MHI_to_Mstar[subsample_cut & (err_flag)]
            plt.scatter(logmass_cut,logsfr_cut,marker='^',color='red',s=30,alpha=0.3,label='VF subsample')
            plt.scatter(logmass_cut[(MHI_to_Mstar_cut>-999)],logsfr_cut[(MHI_to_Mstar_cut>-999)],
                    c=MHI_to_Mstar_cut[(MHI_to_Mstar_cut>-999)], cmap='viridis', s=60, alpha=0.9,
                    label='Subsample with HI measurements',edgecolor='black')
            cb = plt.colorbar()
            cb.set_label(label=r'M$_{HI}$/M$_*$',size=25)
            plt.clim(0,1)
           
        if show_D25:
            hyp_tab = Table.read(path_to_dir+'vf_v2_hyperleda.fits')
            #prepare D25 data for colorbar
            d25 = 10**(hyp_tab['logd25'])       
            d25_cut = d25[subsample_cut & (err_flag)]
            plt.scatter(logmass_cut,logsfr_cut, c=d25_cut, cmap='viridis', s=60, 
                        alpha=0.6, label='VF Subsample')
            cb = plt.colorbar()
            cb.set_label(label=r'Optical D25',size=23)
            plt.clim(15,40)
#############################################################
            
        if show_sizerat:
            plt.scatter(self.altmagphys_cut['combined_logMstar_med'][err_flag_cut], self.altmagphys_cut['combined_logSFR_med'][err_flag_cut], c=self.sizerats[err_flag_cut], cmap='viridis', s=60, alpha=0.6, label='VF Subsample',zorder=2)
            cb = plt.colorbar()
            if self.W1:
                cb.set_label(label=r'R$_{12}$/R$_{3.4}$',size=23)
                plt.clim(0.65,1.2)
            if not self.W1:
                cb.set_label(label=r'R$_{12}$/R$_r$',size=23)
                plt.clim(0.1,1)

        #add z0mgs main sequence line
        #y = logmass*0.47 - 4.88
        #plt.plot(logmass,y,label=('Main Sequence'),color='black',alpha=0.6)
        
        #fit main-sequence line to data
    
        #fit to sSFR constraints
        p1 = np.poly1d(np.polyfit(logmass[((logsfr-logmass)>-11.5)],logsfr[((logsfr-logmass)>-11.5)],1)) 
        #p2 = np.poly1d(np.polyfit(logmass[(logmass-logsfr>-11.5)],logsfr[(logmass-logsfr>-11.5)],2))

        x_p = np.linspace(np.min(logmass),np.max(logmass),1000)

        plt.plot(x_p,p1(x_p),color='black',label='Main Sequence fit',alpha=0.7,linewidth=3,zorder=2)
        #plt.plot(x_p,p2(x_p),color='blue',label='2d fit')
        print(f'p1 parameters: m={np.round(p1[1],2)}, b={np.round(p1[0],2)}')
        
        xplot = np.sort(logmass,axis=None)
        
        plt.plot([xplot[0],np.max(xplot)], [-11.5+xplot[0],-11.5+np.max(xplot)], color='red', linestyle=':', alpha=0.5, label='log(sSFR)>-11.5 limit',zorder=3)
        
        plt.axhline(-1.318,color='blue',linestyle='-.',alpha=0.5,label='log(SFR)>-1.32 limit',zorder=3)
        plt.axvline(8.26,color='green',linestyle='--',alpha=0.5,label=r'log(M$_*$)>8.26 limit',zorder=3)
        
        #plt.scatter(logmass_cut[((logsfr_cut-logmass_cut)<-11.5)], logsfr_cut[((logsfr_cut-logmass_cut)<-11.5)] color='crimson', facecolor='None',s=100,zorder=3)
        #plt.scatter(logmass_cut[((logsfr_cut)<-1.318)], logsfr_cut[((logsfr_cut)<-1.318)], color='crimson', facecolor='None',s=100,zorder=3)
        #plt.scatter(logmass_cut[((logmass_cut)<8.26)], logsfr_cut[((logmass_cut)<8.26)], color='crimson', facecolor='None',s=100,zorder=3)
    
        print()
        print(f'ngal below Mstar limit: {len(logsfr_cut[((logmass_cut)<8.26)])}')
        print(f'ngal below SFR limit: {len(logsfr_cut[((logsfr_cut)<-1.318)])}')
        print(f'ngal below sSFR limit: {len(logsfr_cut[((logsfr_cut-logmass_cut)<-11.5)])}')
        print(f'total ngal circled (note that this is NOT the sum of the above three values): {len(logmass_cut[(logmass_cut<8.26) | (logsfr_cut<-1.318) | ((logsfr_cut-logmass_cut)<-11.5)])}')
        
        cb.ax.tick_params(labelsize=15)
        
        plt.xlabel(r'log($M_*$/$M_\odot$)',fontsize=22)
        plt.ylabel(r'log(SFR/($M_\odot$/yr))',fontsize=22)
        
        plt.xlim(6,11.6)
        plt.ylim(-5,2.3)
        
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        
        plt.legend(fontsize=12)

        if savefig==True:
            plt.savefig(homedir+'/Desktop/sfrmstar_magphys.png', dpi=100, bbox_inches='tight', pad_inches=0.2)
        
        plt.show()
        
    def sfrmstar_z0mgs(self, show_HI=False, show_D25=False, show_sizerat=True, savefig=False):
        
        #prepare D25 data for colorbar
        d25 = 10**(self.hyp_tab_cut['logd25'])
        
        #prepare MHI_to_Mstar data for colorbar
        MHI_to_Mstar = self.HI_tab['MHI']/(10**self.z0mgs['logmass'])
        
        for n in range(len(MHI_to_Mstar)):
            if self.HI_tab['MHI'].mask[n]:     #if value is masked, set to be -999
                MHI_to_Mstar[n] = -999
            else:
                MHI_to_Mstar[n] = MHI_to_Mstar[n] if ((np.log(self.HI_tab['MHI'][n])>0) & (MHI_to_Mstar[n]<10) & (MHI_to_Mstar[n]!=1)) else -999
        
        path_to_dir = homedir+'/Desktop/v2-20220820/'
        logsfr = Table.read(path_to_dir+'vf_v2_z0mgs.fits')['logsfr']
        logmass = Table.read(path_to_dir+'vf_v2_z0mgs.fits')['logmass']
        
        logsfr_cut = self.z0mgscut['logsfr']
        logmass_cut = self.z0mgscut['logmass']
        
        MHI_to_Mstar_cut = MHI_to_Mstar
        d25_cut = d25

        plt.figure(figsize=(10,6))
        plt.scatter(logmass,logsfr,color='gray',s=3,alpha=0.05,label='VF sample')
        
        if show_HI:
            plt.scatter(logmass_cut,logsfr_cut,marker='^',color='red',s=30,alpha=0.3,label='VF subsample')
            plt.scatter(logmass_cut[(MHI_to_Mstar_cut>-999)],logsfr_cut[(MHI_to_Mstar_cut>-999)],
                    c=MHI_to_Mstar_cut[(MHI_to_Mstar_cut>-999)], cmap='viridis', s=60, alpha=0.9,
                    label='Subsample with HI measurements',edgecolor='black')
            cb = plt.colorbar()
            cb.set_label(label=r'M$_{HI}$/M$_*$',size=25)
            plt.clim(0,1)
        
        if show_sizerat:
            plt.scatter(self.magphyscut['logMstar_med'][err_flag_cut],self.magphyscut['logSFR_med'][err_flag_cut], 
                        c=self.sizerats[err_flag_cut], cmap='viridis', s=60, alpha=0.6, label='VF Subsample')
            cb = plt.colorbar()
            cb.set_label(label=r'R$_{12}$/R$_r$',size=23)
            plt.clim(0.1,1)
            
        if show_D25:
            plt.scatter(logmass_cut, logsfr_cut, c=d25_cut, cmap='viridis', s=60, 
                        alpha=0.6, label='VF Subsample')
            cb = plt.colorbar()
            cb.set_label(label=r'Optical D25',size=23)
            plt.clim(15,40)
            
        #add main sequence line
        y = logmass*0.47 - 4.88
        plt.plot(logmass,y,label=('Main Sequence'),color='black',alpha=0.6)
        
        cb.ax.tick_params(labelsize=15)
        
        plt.xlabel(r'log(M*/$M_\odot$)',fontsize=22)
        plt.ylabel(r'log(SFR/($M_\odot$/yr))',fontsize=22)
        
        plt.xlim(6,11.6)
        plt.ylim(-5,2.3)
        
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        
        plt.legend(fontsize=12)

        if savefig==True:
            plt.savefig(homedir+'/Desktop/sfrmstar_z0mgs.png', dpi=100, bbox_inches='tight', pad_inches=0.2)
        
        plt.show()
     
    def r12_vs_rstar(self, savefig=False):
        
        logsfr = self.altmagphys_cut['combined_logSFR_med']
        logmass = self.altmagphys_cut['combined_logMstar_med']
        
        #remove entries where there is no magphys data available for that galaxy
        err_flag = (self.altmagphys_cut['magphys_flag']) & (logmass!=0.)
        
        logsfr = logsfr[err_flag]
        logmass = logmass[err_flag]

        if self.W1:
            r_arcsec = self.re_w1band_cut[err_flag]*2.75
        if not self.W1:
            r_arcsec = self.re_rband_cut[err_flag]*0.262
        w3_arcsec = self.re_w3band_cut[err_flag]*2.75
        
        rr_dat = [r_arcsec,
                  r_arcsec[self.clusflag[err_flag]],
                  r_arcsec[self.rgflag[err_flag]],
                   r_arcsec[self.pgflag[err_flag]],
                  r_arcsec[self.filflag[err_flag]],
                    r_arcsec[self.fieldflag[err_flag]]]
        w3_dat = [w3_arcsec,
                  w3_arcsec[self.clusflag[err_flag]],
                  w3_arcsec[self.rgflag[err_flag]],
                   w3_arcsec[self.pgflag[err_flag]],
                  w3_arcsec[self.filflag[err_flag]],
                  w3_arcsec[self.fieldflag[err_flag]]]
        
        logmass_dat = [logmass,
                       logmass[self.clusflag[err_flag]],
                       logmass[self.rgflag[err_flag]],
                       logmass[self.pgflag[err_flag]],
                       logmass[self.filflag[err_flag]],
                       logmass[self.fieldflag[err_flag]]]
        
        logsfr_dat = [logsfr,
                       logsfr[self.clusflag[(err_flag)]],
                       logsfr[self.rgflag[(err_flag)]],
                       logsfr[self.pgflag[(err_flag)]],
                       logsfr[self.filflag[(err_flag)]],
                       logsfr[self.fieldflag[(err_flag)]]]
        
        labels = ['Full Subsample','Cluster','Rich Group','Poor Group','Filament','Field']
        
        colors=['purple','crimson','black','magenta','green','blue']
        
        #take upper left panel of mstar figure and create one row of three panels, separated according to different mass bins; color-code by sSFR. this plot will be generated regardless of whether the user selects mstar or sfr.
        
        logsfr = logsfr[logmass>=8.26]
        r_arcsec = r_arcsec[logmass>=8.26]
        w3_arcsec = w3_arcsec[logmass>=8.26]
        envflags = self.v2_envcut[err_flag]
        envflags = envflags[logmass>=8.26]
        
        logmass = logmass[logmass>=8.26]

        max_mass = np.max(logmass)
        min_mass = np.min(logmass)
        divisor = (max_mass - min_mass)/3
                        
        mass_bin1 = (logmass>=min_mass) & (logmass<(min_mass+divisor))
        mass_bin2 = (logmass>=(min_mass+divisor)) & (logmass<(min_mass+2*divisor))
        mass_bin3 = (logmass>=(min_mass+2*divisor)) & (logmass<(min_mass+3*divisor))
        
        #mass bin flags
        mass_bins = [mass_bin1,mass_bin2,mass_bin3,mass_bin1,mass_bin2,mass_bin3,mass_bin1,mass_bin2,mass_bin3]
                
        titles = [f'{round(min_mass,2)} <'+r' $M_*$ '+f'< {round(min_mass+divisor,2)}',
                  f'{round(min_mass+divisor,2)} <'+r' $M_*$ '+f'< {round(min_mass+2*divisor,2)}',
                 f'{round(min_mass+2*divisor,2)} <'+r' $M_*$ '+f'< {round(min_mass+3*divisor,2)}']
        
        labels = ['Cluster','Cluster','Cluster','Filament/Groups','Filament/Groups','Filament/Groups',
                  'Field','Field','Field']
        
        filgroupflag = ((envflags['rich_group_memb'])|(envflags['poor_group_memb'])|(envflags['filament_member']))
        
        fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(34,24))
        plt.subplots_adjust(hspace=.05,wspace=.05)
        
        for n, ax in enumerate(axes.flat):      #j==index, i==value          
            n = n+1

            #color-code by log(sSFR)
            if n<4:
                im = ax.scatter(r_arcsec[mass_bins[n-1]&envflags['cluster_member']],
                            w3_arcsec[mass_bins[n-1]&envflags['cluster_member']],edgecolors='black',
                        c=(logsfr[mass_bins[n-1]&envflags['cluster_member']]-logmass[mass_bins[n-1]&envflags['cluster_member']]),s=100,cmap='plasma')
                
                print('std for cluster:',np.std(w3_arcsec[mass_bins[n-1]&envflags['cluster_member']]/r_arcsec[mass_bins[n-1]&envflags['cluster_member']]))
            
            elif (n>=4) & (n<=6):
                im = ax.scatter(r_arcsec[mass_bins[n-1]&filgroupflag],w3_arcsec[mass_bins[n-1]&filgroupflag],
                            edgecolors='black',s=100,
                            c=(logsfr[mass_bins[n-1]&filgroupflag]-logmass[mass_bins[n-1]&filgroupflag]),
                            cmap='plasma')        

                print('std for filament/groups:',np.std(w3_arcsec[mass_bins[n-1]&filgroupflag]/r_arcsec[mass_bins[n-1]&filgroupflag]))
                
            else:
                im = ax.scatter(r_arcsec[mass_bins[n-1]&envflags['pure_field']],
                            w3_arcsec[mass_bins[n-1]&envflags['pure_field']],edgecolors='black',
                        c=(logsfr[mass_bins[n-1]&envflags['pure_field']]-logmass[mass_bins[n-1]&envflags['pure_field']]),s=100,cmap='plasma')

                print('std for field:', np.std(w3_arcsec[mass_bins[n-1]&envflags['pure_field']]/r_arcsec[mass_bins[n-1]&envflags['pure_field']]))
            
            im.set_clim(-11,-9)
            
            ax.axline((0, 0), slope=1, color='indigo')
            
            ax.set_xlim(2, 2e2)
            ax.set_ylim(1, 130)
            
            ax.set_xscale('log')
            ax.set_yscale('log')
            if n in [1,2,3]:
                ax.set_title(titles[n-1],fontsize=35,pad=20)
            if n in [7,8,9]:
                if not self.W1:
                    ax.set_xlabel(r'R$_{r}$ (arcsec)',fontsize=30)
                else:
                    ax.set_xlabel(r'R$_{3.4}$ (arcsec)',fontsize=30)
            if n in [1,4,7]:
                ax.set_ylabel(r'R$_{12}$ (arcsec)',fontsize=30)
            
            ax.set_xticks([])
            ax.set_yticks([])
            plt.xticks(fontsize=25)
            plt.yticks(fontsize=25)
            
            plt.text(.05, .95, labels[n-1], ha='left', va='top',transform=ax.transAxes,fontsize=35)
        
        cb = fig.colorbar(im, ax=axes.ravel().tolist(), pad=0.02, aspect=30)
        cb.set_label(label=r'log(sSFR)',size=30)
        cb.ax.tick_params(labelsize=25)
        #tick information from 
        #https://stackoverflow.com/questions/22012096/how-to-set-number-of-ticks-in-plt-colorbar
        tick_locator = ticker.MaxNLocator(nbins=5)
        cb.locator = tick_locator
        cb.update_ticks()
        
        if savefig==True:
            plt.savefig(homedir+'/Desktop/r12_rstar_all.png', dpi=100, bbox_inches='tight', pad_inches=0.2)
        
        plt.show()
        
    def r12_vs_rstar_v2(self, savefig=False):
        
        logsfr = self.altmagphys_cut['combined_logSFR_med']
        logmass = self.altmagphys_cut['combined_logMstar_med']
        
        #remove entries where there is no magphys data available for that galaxy
        err_flag = (self.altmagphys_cut['magphys_flag']) & (logmass!=0.)
        
        logsfr = logsfr[err_flag]
        logmass = logmass[err_flag]
        envflags = self.v2_envcut[err_flag]

        w3_arcsec = self.re_w3band_cut[err_flag]*2.75
        r_arcsec = self.re_w1band_cut[err_flag]*2.75
                        
        #take upper left panel of mstar figure and create one row of three panels, separated according to different mass bins; color-code by sSFR. this plot will be generated regardless of whether the user selects mstar or sfr.
        
        logsfr = logsfr[logmass>=8.26]
        r_arcsec = r_arcsec[logmass>=8.26]
        w3_arcsec = w3_arcsec[logmass>=8.26]
        envflags = envflags[logmass>=8.26]
        logmass = logmass[logmass>=8.26]
                        
        mass_bin1 = (logmass < 9.5)
        mass_bin2 = (logmass > 9.5)
        
        x = [r_arcsec[mass_bin1&envflags['cluster_member']], r_arcsec[mass_bin1&envflags['rich_group_memb']],
            r_arcsec[mass_bin1&envflags['poor_group_memb']], r_arcsec[mass_bin1&envflags['filament_member']],
            r_arcsec[mass_bin1&envflags['pure_field']], r_arcsec[mass_bin2&envflags['cluster_member']],
            r_arcsec[mass_bin2&envflags['rich_group_memb']], r_arcsec[mass_bin2&envflags['poor_group_memb']],
            r_arcsec[mass_bin2&envflags['filament_member']], r_arcsec[mass_bin2&envflags['pure_field']]]
        
        y = [w3_arcsec[mass_bin1&envflags['cluster_member']], w3_arcsec[mass_bin1&envflags['rich_group_memb']],
            w3_arcsec[mass_bin1&envflags['poor_group_memb']], w3_arcsec[mass_bin1&envflags['filament_member']], 
            w3_arcsec[mass_bin1&envflags['pure_field']], w3_arcsec[mass_bin2&envflags['cluster_member']],
            w3_arcsec[mass_bin2&envflags['rich_group_memb']], w3_arcsec[mass_bin2&envflags['poor_group_memb']], 
            w3_arcsec[mass_bin2&envflags['filament_member']], w3_arcsec[mass_bin2&envflags['pure_field']]]
        
        cs = [logsfr[mass_bin1&envflags['cluster_member']]-logmass[mass_bin1&envflags['cluster_member']], 
              logsfr[mass_bin1&envflags['rich_group_memb']]-logmass[mass_bin1&envflags['rich_group_memb']], 
              logsfr[mass_bin1&envflags['poor_group_memb']]-logmass[mass_bin1&envflags['poor_group_memb']], 
              logsfr[mass_bin1&envflags['filament_member']]-logmass[mass_bin1&envflags['filament_member']], 
              logsfr[mass_bin1&envflags['pure_field']]-logmass[mass_bin1&envflags['pure_field']], 
              logsfr[mass_bin2&envflags['cluster_member']]-logmass[mass_bin2&envflags['cluster_member']], 
              logsfr[mass_bin2&envflags['rich_group_memb']]-logmass[mass_bin2&envflags['rich_group_memb']], 
              logsfr[mass_bin2&envflags['poor_group_memb']]-logmass[mass_bin2&envflags['poor_group_memb']], 
              logsfr[mass_bin2&envflags['filament_member']]-logmass[mass_bin2&envflags['filament_member']], 
              logsfr[mass_bin2&envflags['pure_field']]-logmass[mass_bin2&envflags['pure_field']]]

        titles = ['Cluster', 'Rich Group', 'Poor Group', 'Filament', 'Field', 'Cluster', 'Rich Group', 'Poor Group', 'Filament', 'Field']
                
        fig, axes = plt.subplots(nrows=2, ncols=5, figsize=(30,10))
        plt.subplots_adjust(hspace=0.,wspace=0.)
        
        for n, ax in enumerate(axes.flat):      #j==index, i==value          
            n = n+1   
            #color-code by log(sSFR)
            im = ax.scatter(x[n-1], y[n-1], edgecolors='black', c=cs[n-1], s=100, cmap='plasma')
            
            if n in [1,2,3,4,5]:
                print(f'std for {titles[n-1]} (logmass < 9.5):',np.std(y[n-1]/x[n-1]))
            else:
                print(f'std for {titles[n-1]} (logmass > 9.5):',np.std(y[n-1]/x[n-1]))
            
            im.set_clim(-11,-9)
            
            ax.axline((0, 0), slope=1, color='indigo')
            
            ax.set_xlim(2, 2e2)
            ax.set_ylim(1, 130)
            
            ax.set_xscale('log')
            ax.set_yscale('log')
            
            if n in [1,2,3,4,5]:
                ax.set_title(titles[n-1],fontsize=25,pad=10)
                plt.text(.05, .95, str(np.round(np.min(logmass),1))+r'<log(M$_*$)<9.5', ha='left', va='top', 
                         transform=ax.transAxes, fontsize=16)
            else:
                ax.set_xlabel(r'R$_{3.4}$ (arcsec)',fontsize=20)
                plt.text(.05, .95, r'9.5<log(M$_*$)<'+str(np.round(np.max(logmass),1)), ha='left', 
                         va='top', transform=ax.transAxes, fontsize=16)
                
            if n in [1,6]:   #leftmost column
                ax.set_ylabel(r'R$_{12}$ (arcsec)',fontsize=20)
                ax.tick_params(axis='both', labelsize=15)
                if n==1:
                    ax.tick_params(labelbottom=False,labelright=False,labeltop=False)
                    
            elif n in [7,8,9,10]:   #all of bottom row
                ax.tick_params(axis='x',labelsize=15)
                ax.tick_params(labelleft=False)
  
            else:   #rest of top row
                ax.tick_params(labelbottom=False,labelleft=False)
            
            #add tick marks to top, right of each panel
            ax2 = ax.secondary_xaxis('top')
            ax2.tick_params(labeltop=False)
            ax3 = ax.secondary_yaxis('right')
            ax3.tick_params(labelright=False)
        
        #add colorbar spanning entire height of subplot figure
        cb = fig.colorbar(im, ax=axes.ravel().tolist(), pad=0.02, aspect=30)
        cb.set_label(label=r'log(sSFR/yr$^{-1}$)',size=30)
        cb.ax.tick_params(labelsize=25)
        #tick info from https://stackoverflow.com/questions/22012096/how-to-set-number-of-ticks-in-plt-colorbar
        tick_locator = ticker.MaxNLocator(nbins=5)
        cb.locator = tick_locator
        cb.update_ticks()
        
        if savefig==True:
            plt.savefig(homedir+'/Desktop/r12_rstar_all.png', dpi=100, bbox_inches='tight', pad_inches=0.2)
        
        plt.show()
        
    def ratio_MS(self, savefig=False):
        
        m = 0.77
        b = -8.03
                
        logsfr = self.altmagphys_cut['combined_logSFR_med']
        logmass = self.altmagphys_cut['combined_logMstar_med']
        
        #fit main-sequence line to data (1D --> p1; parameters taken from sfrmstar_magphys)
        Mstar_full = logmass
        
        #remove entries where there is no magphys data available for that galaxy
        err_flag_cut = (self.altmagphys_cut['magphys_flag']) & (logmass>0.)

        #Now calculate distance of each point to the best-fit line
        #Distance = (| a*x1 + b*y1 + c |) / (sqrt( a*a + b*b))
        dist = ((m*logmass[err_flag_cut]) + (-1*logsfr[err_flag_cut]) + b) / (np.sqrt((m)**2 + (1)**2))
        dist = -1*dist

        plt.figure(figsize=(10,6))
        plt.axhline(1,linestyle='--',color='r',alpha=0.4)
        
        t=self.re_w3band_cut[err_flag_cut]
        t1=self.re_w1band_cut[err_flag_cut]
        t2=self.v2_envcut[err_flag_cut]
        
        plt.scatter(dist,self.sizerats[err_flag_cut],color='blue',s=40,alpha=0.7)   #all points
        plt.xlabel(r'Distance from Main Sequence Line',fontsize=17)
        
        if self.W1:
            plt.ylabel(r'Size Ratio ($R_{12}/R_{3.4}$)',fontsize=17)
            #plt.yscale('log')
            plt.ylim(0,2)    
            
            print(f'# galaxies within 0.5 dex of MS: {len(self.sizerats[err_flag_cut & (np.abs(dist<0.5))])}')              
        if not self.W1:
            plt.ylabel(r'Size Ratio ($R_{12}/R_{r}$)',fontsize=17)
            plt.ylim(0, 1.75)   #artificially trim outliers with size ratios>5 

        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)

        res = spearmanr(dist,self.sizerats[err_flag_cut])

        print('np.polyfit m,b:',np.polyfit(dist,self.sizerats[err_flag_cut],1))
        print('spearmanr:',res)
        
        self.clusflag = self.v2_envcut['cluster_member']
        self.rgflag = self.v2_envcut['rich_group_memb']
        self.pgflag = self.v2_envcut['poor_group_memb']
        self.filflag = self.v2_envcut['filament_member']
        self.fieldflag = self.v2_envcut['pure_field']
        
        #######Testing whether an environment correlation exists#######
        
        distclus = ((m*logmass[err_flag_cut&self.clusflag]) + (-1*logsfr[err_flag_cut&self.clusflag]) + b) / (np.sqrt((m)**2 + (1)**2))
        distrg = ((m*logmass[err_flag_cut&self.rgflag]) + (-1*logsfr[err_flag_cut&self.rgflag]) + b) / (np.sqrt((m)**2 + (1)**2))
        distpg = ((m*logmass[err_flag_cut&self.pgflag]) + (-1*logsfr[err_flag_cut&self.pgflag]) + b) / (np.sqrt((m)**2 + (1)**2))
        distfil = ((m*logmass[err_flag_cut&self.filflag]) + (-1*logsfr[err_flag_cut&self.filflag]) + b) / (np.sqrt((m)**2 + (1)**2))
        distfield = ((m*logmass[err_flag_cut&self.fieldflag]) + (-1*logsfr[err_flag_cut&self.fieldflag]) + b) / (np.sqrt((m)**2 + (1)**2))
        
        
        resclus = spearmanr(distclus,self.sizerats[err_flag_cut&self.clusflag])
        resrg = spearmanr(distrg,self.sizerats[err_flag_cut&self.rgflag])
        respg = spearmanr(distpg,self.sizerats[err_flag_cut&self.pgflag])
        resfil = spearmanr(distfil,self.sizerats[err_flag_cut&self.filflag])
        resfield = spearmanr(distfield,self.sizerats[err_flag_cut&self.fieldflag])
        
        print()
        print('cluster spearmanr:',resclus)
        print('rich group spearmanr:',resrg)
        print('poor group spearmanr:',respg)
        print('filament spearmanr:',resfil)
        print('field spearmanr:',resfield)
        ###############################################################
        
        if savefig==True:
            plt.savefig(homedir+'/Desktop/MS_dist.png', dpi=100, bbox_inches='tight', pad_inches=0.2)
        
        plt.show()

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
        
        print('Fractional Breakdown:')
        print(f'Cluster: {round(len(self.v2_envcut[self.clusflag])/len(self.v2_envcut)*100,1)}%')
        print(f'Rich Group: {round(len(self.v2_envcut[self.rgflag])/len(self.v2_envcut)*100,1)}%')
        print(f'Poor Group: {round(len(self.v2_envcut[self.pgflag])/len(self.v2_envcut)*100,1)}%')
        print(f'Filament: {round(len(self.v2_envcut[self.filflag])/len(self.v2_envcut)*100,1)}%')
        print(f'Field: {round(len(self.v2_envcut[self.fieldflag])/len(self.v2_envcut)*100,1)}%')
        print(f'RG+Filament: {round(len(self.v2_envcut[self.rgflag&self.filflag])/len(self.v2_envcut)*100,1)}%')
        print(f'PG+Filament: {round(len(self.v2_envcut[self.pgflag&self.filflag])/len(self.v2_envcut)*100,1)}%')

        if savefig==True:
            plt.savefig(homedir+'/Desktop/envbins.png', dpi=100)
        
        plt.show()

    def env_means_comp(self, mass_match=False, errtype='bootstrap', savefig=False):    
        index = np.arange(1,6,1)
        env_names = ['Cluster','Rich \n Group','Poor \n Group','Filament','Field']
        
        #prepare WISESize projected environment sizes for standard errors
        #   [ratio_of_VFS_galaxies_in_environment] * N_WISESIZE
        ntot_wise = 13756
        n_wise = np.asarray([len(self.v2_envcut[self.clusflag])*ntot_wise/len(self.v2_envcut), 
                  len(self.v2_envcut[self.rgflag])*ntot_wise/len(self.v2_envcut), 
                  len(self.v2_envcut[self.pgflag])*ntot_wise/len(self.v2_envcut),
                  len(self.v2_envcut[self.filflag])*ntot_wise/len(self.v2_envcut), 
                  len(self.v2_envcut[self.fieldflag])*ntot_wise/len(self.v2_envcut)])
        
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
        err_wisesize = []
        err_upper_bootstrap = []
        err_lower_bootstrap = []
        
        for j,i in enumerate(re_data):      #j==index, i==value
            if self.MeanMedian=='mean':
                central_pts.append(np.mean(i))
                lower_err, upper_err = get_bootstrap_confint(i,bootfunc=np.mean,nboot=200)

            if self.MeanMedian=='median':
                central_pts.append(np.median(i))
                lower_err, upper_err = get_bootstrap_confint(i,bootfunc=np.median,nboot=200)
            
            err.append(np.std(i)/np.sqrt(len(i)))
            err_wisesize.append(np.std(i)/np.sqrt(n_wise[j]))
            err_upper_bootstrap.append(upper_err)
            err_lower_bootstrap.append(lower_err)
        
        err_color = 'orangered'
        plt.figure(figsize=(7,7))
        #plt.scatter(index,central_pts,color='blue',s=50,zorder=2,edgecolors='black',label=self.MeanMedian)
       
        if errtype!='bootstrap':
            plt.errorbar(index,central_pts,yerr=err,fmt='None',color='orange',linewidth=8,label='Current Precision - VFS')
        
        plt.errorbar(index,central_pts,yerr=err_wisesize,fmt='None',color='black',linewidth=8,label='Projected Precision - WISESize',zorder=4)
        
        xmin,xmax = plt.xlim()
        xfield = np.linspace(xmin,xmax,50)
        ymax = np.ones(50)*(central_pts[-1] + err[-1])
        ymin = np.ones(50)*(central_pts[-1] - err[-1])
        ymax_wise = np.ones(50)*(central_pts[-1] + err_wisesize[-1])
        ymin_wise = np.ones(50)*(central_pts[-1] - err_wisesize[-1])
        
        if errtype=='bootstrap':
            labels = [None,None,None,None,'Current Precision - VFS']
            for n in range(5):
                plt.plot([index[n],index[n]], [err_lower_bootstrap[n],err_upper_bootstrap[n]],
                         color='orange',linewidth=8,label=labels[n],zorder=3)
            ymax = np.ones(50)*(err_upper_bootstrap[-1])
            ymin = np.ones(50)*(err_lower_bootstrap[-1])
        
        plt.fill_between(xfield,ymax,ymin,color=err_color,alpha=.1,zorder=1)
        plt.fill_between(xfield,ymax_wise,ymin_wise,color='black',alpha=.1,zorder=2)
        
        plt.xticks(index, env_names, rotation=10, fontsize=20)
        plt.tick_params(axis='both', which='major', labelsize=15)
        plt.grid(alpha=0.2)
        
        if self.W1:
            plt.ylabel(r'R$_{12}$/R$_{3.4}$',fontsize=20)
        else:
            plt.ylabel(r'R$_{12}$/R$_r$',fontsize=20)
            plt.ylim(0.6,.85)
        
        leg = plt.legend(fontsize=18)
        
        if savefig==True:
            plt.savefig(homedir+'/Desktop/env_ratio_comp.png', bbox_inches='tight', pad_inches=0.2, dpi=100)
            
        plt.show()
        
    def env_means(self, mass_match=False, trimOutliers=False, errtype='bootstrap', r90=False, savefig=False):    
        index = np.arange(1,6,1)
        env_names = ['Cluster','Rich \n Group','Poor \n Group','Filament','Field']
        
        #will generate the self.outlier_flag variable needed to, well, trim the outliers.
        ratios = self.sizerats
        clusflag = self.clusflag.copy()
        rgflag = self.rgflag.copy()
        pgflag = self.pgflag.copy()
        filflag = self.filflag.copy()
        fieldflag = self.fieldflag.copy()
        
        #if r90, then change ratio type accordingly. all else should be the same.
        if r90:
            ratios = np.zeros(len(self.sizerats))
            for i in range(len(self.sizerats)):
                ratios[i] = calculate_r90(self.re_w3band_cut[i],self.n_w3band_cut[i])/calculate_r90(self.re_w1band_cut[i],self.n_w1band_cut[i])
        
        re_data = [ratios[clusflag],ratios[rgflag],ratios[pgflag],ratios[filflag],
                   ratios[fieldflag]]
        
        for i in range(len(ratios)):
            if (ratios[i]==0) | (ratios[i]==np.nan):
                print(ratios[i],self.re_w3band_cut[i],self.re_w1band_cut[i],calculate_r90(self.re_w3band_cut[i],self.n_w3band_cut[i]),
                     calculate_r90(self.re_w1band_cut[i],self.n_w1band_cut[i]))
        
        #create 100 iterations of size ratio v. environment bin; plot median and STD, 
        #compare with non-mass matching result.
        if mass_match:
            print('Completing 500 mass-match iterations...')
            
            meds = [[],[],[],[],[]]
            err = [[],[],[],[],[]]
            meds_all = [[],[],[],[],[]]  #cluster, rg, pg, filament, field
            err_all = [[],[],[],[],[]]
            for n in range(500):
                self.mass_matching()    #field rich group poor group filament
                re_data = [ratios[clusflag], np.asarray(self.matched_sizes[3]), np.asarray(self.matched_sizes[2]), 
                           np.asarray(self.matched_sizes[1]), np.asarray(self.matched_sizes[0])]
                for i in range(len(re_data)):
                    meds[i].append(np.median(re_data[i]))

            for i in range(5):
                meds_all[i].append(np.median(meds[i]))
                #err_all[i].append(np.std(re_data[i])/np.sqrt(len(re_data[i])))
                err_all[i].append(np.std(meds[i]))
                
            #remove list of lists; need one list of five floats.
            err_all_arr = []
            for item in err_all:
                err_all_arr.append(item[0])
            meds_all_arr = []
            for item in meds_all:
                meds_all_arr.append(item[0])
                        
            plt.figure(figsize=(10,6))
            plt.scatter(index,meds_all_arr,color='blue',s=50,zorder=2,edgecolors='black',label='Median')
            plt.errorbar(index,meds_all_arr,yerr=err_all_arr,fmt='None',color='orangered',zorder=1)
            
            xmin,xmax = plt.xlim()
            xfield = np.linspace(xmin,xmax,50)
            ymax = np.ones(50)*(meds_all_arr[-1] + err_all_arr[-1])
            ymin = np.ones(50)*(meds_all_arr[-1] - err_all_arr[-1])
            
            plt.fill_between(xfield,ymax,ymin,color='orangered',alpha=.1)
        
            plt.xticks(index, env_names, rotation=10, fontsize=20)
            plt.tick_params(axis='both', which='major', labelsize=15)
            plt.grid(alpha=0.2)

            plt.show()

            #now re-define ra_data in order to plot 'just one' of the mass-matching iterations
            re_data = [ratios[clusflag], np.asarray(self.matched_sizes[3]), np.asarray(self.matched_sizes[2]), 
                           np.asarray(self.matched_sizes[1]), np.asarray(self.matched_sizes[0])]

        central_pts = []
        low_25 = []   #25% lowest value in bin
        high_75 = []  #75% highest value in bin
        
        err = []
        err_upper_bootstrap = []
        err_lower_bootstrap = []
        
        for j,i in enumerate(re_data):      #j==index, i==value
            if self.MeanMedian=='mean':
                central_pts.append(np.mean(i))
                lower_err, upper_err = get_bootstrap_confint(i,bootfunc=np.mean,nboot=200)
                reordered_i = np.sort(i)
                low_25.append(reordered_i[int((len(i)-1)*0.25)])
                high_75.append(reordered_i[int((len(i)-1)*0.75)])

            if self.MeanMedian=='median':
                central_pts.append(np.median(i))
                lower_err, upper_err = get_bootstrap_confint(i,bootfunc=np.median,nboot=200)
                reordered_i = np.sort(i)
                low_25.append(reordered_i[int((len(i)-1)*0.25)])
                high_75.append(reordered_i[int((len(i)-1)*0.75)])
            
            err.append(np.std(i)/np.sqrt(len(i)))
            err_upper_bootstrap.append(upper_err)
            err_lower_bootstrap.append(lower_err)
        
        err_color = 'orangered'
        plt.figure(figsize=(10,6))
        plt.scatter(index,central_pts,color='blue',s=50,zorder=2,edgecolors='black',label=self.MeanMedian)
        
        #plt.scatter(index,low_25,color='green',s=30,edgecolors='black',label='25% value')
        #plt.plot(index,low_25,color='green',alpha=0.3,ls='--')
        #plt.scatter(index,high_75,color='green',s=30,edgecolors='black',label='75% value')
        #plt.plot(index,high_75,color='green',alpha=0.3,ls='--')
       
        if errtype!='bootstrap':
            plt.errorbar(index,central_pts,yerr=err,fmt='None',color=err_color,zorder=1)

        xmin,xmax = plt.xlim()
        xfield = np.linspace(xmin,xmax,50)
        ymax = np.ones(50)*(central_pts[-1] + err[-1])
        ymin = np.ones(50)*(central_pts[-1] - err[-1])
        
        if errtype=='bootstrap':
            for n in range(5):
                plt.plot([index[n],index[n]], [err_lower_bootstrap[n],err_upper_bootstrap[n]],color='blue',zorder=1)
            ymax = np.ones(50)*(err_upper_bootstrap[-1])
            ymin = np.ones(50)*(err_lower_bootstrap[-1])
        
        plt.fill_between(xfield,ymax,ymin,color=err_color,alpha=.1)
        
        plt.xticks(index, env_names, rotation=10, fontsize=20)
        plt.tick_params(axis='both', which='major', labelsize=15)
        plt.grid(alpha=0.2)
        plt.ylabel(r'R$_{12}$/R$_r$',fontsize=20)
    
        if self.W1:
            plt.ylabel(r'R$_{12}$ / R$_{3.4}$',fontsize=20)
            if r90:
                plt.ylabel(r'R90$_{12}$ / R90$_{3.4}$',fontsize=20)
        else:
            plt.ylabel(r'R50$_{12}$ / R50$_r$',fontsize=20)
            if r90:
                plt.ylabel(r'R90$_{12}$ / R90$_r$',fontsize=20)
        
        #plt.legend(fontsize=15)
        
        if savefig==True:
            plt.savefig(homedir+'/Desktop/env_ratio.png', bbox_inches='tight', pad_inches=0.2, dpi=100)
            
        plt.show()    
        
        t_stat, p_val = ttest_1samp(a=re_data[0], popmean = central_pts[-1])
        
        print(f'cluster median: {np.round(central_pts[0],3)} +{np.round(err_upper_bootstrap[0]-central_pts[0],3)} -{np.round(central_pts[0]-err_lower_bootstrap[0],3)}')
        print(f'field median: {np.round(central_pts[-1],3)} +{np.round(err_upper_bootstrap[-1]-central_pts[-1],3)} -{np.round(central_pts[-1]-err_lower_bootstrap[-1],3)}')
        
        print(t_stat,p_val)
        
    
    def sfr_hist(self, z0mgs_comp=False, savefig=False):
        
        z0mgs_sfr = self.z0mgscut['logsfr']
        magphys_sfr = self.altmagphys_cut['combined_logSFR_med']
        
        #remove entries where there is no magphys data available for that galaxy
        err_flag = (self.altmagphys_cut['magphys_flag']) & (magphys_sfr!=0)  #ALSO INCLUDE != 0 values. 
    
        z0mgs_env_sfr = [z0mgs_sfr[self.clusflag&err_flag], z0mgs_sfr[self.rgflag&err_flag], 
                          z0mgs_sfr[self.pgflag&err_flag], 
                          z0mgs_sfr[self.filflag&err_flag], z0mgs_sfr[self.fieldflag&err_flag]]
        magphys_env_sfr = [magphys_sfr[self.clusflag&err_flag], magphys_sfr[self.rgflag&err_flag], 
                            magphys_sfr[self.pgflag&err_flag], 
                            magphys_sfr[self.filflag&err_flag], magphys_sfr[self.fieldflag&err_flag]]

        env_names = ['Cluster','Rich Group','Poor Group','Filament','Field']

        mybins=np.linspace(np.min(magphys_sfr),np.max(magphys_sfr),1000)
        #mybins=1000
        
        fig = plt.figure(figsize=(11,8))
        plt.subplots_adjust(hspace=.4,wspace=.3)
        
        for i in range(1,6):
            ax=fig.add_subplot(2,3,i)
            if i == 4:
                ax.set_position([0.24,0.125,0.228,0.3])
            if i == 5:
                ax.set_position([0.55,0.125,0.228,0.3])
            
            plt.hist(magphys_env_sfr[i-1],bins=mybins,alpha=0.8,cumulative=True,density=True,
                     histtype='step',label='Magphys')
            
            if z0mgs_comp==True:
                plt.hist(z0mgs_env_sfr[i-1],bins=mybins,alpha=0.8,cumulative=True,density=True,
                         histtype='step',label='z0mgs')
            
            if (i == 1)&(z0mgs_comp==True):
                plt.legend(fontsize=13,loc='upper left')
            plt.xlabel(r'log(SFR/[$M_\odot$ yr$^{-1}$])',fontsize=16)
            plt.xlim(-1.4,0.85)
            plt.title(env_names[i-1],fontsize=16)
            
            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)
        
        print('K-S p-value (> 0.003 (3sigma), "same distribution"):')
        for n in range(len(env_names)):
            print(env_names[n])
            print('%.5f'%(kstest(z0mgs_env_sfr[n],magphys_env_sfr[n])[1]))
        
        if savefig==True:
            plt.savefig(homedir+'/Desktop/sfr_hist.png',bbox_inches='tight', pad_inches=0.2, dpi=100)
        
        plt.show()
        
        fig2 = plt.figure(figsize=(6,4))
        plt.scatter(magphys_sfr,z0mgs_sfr,s=10)
        plt.axline((0, 0), slope=1, color='indigo')
        plt.xlabel('[MAGPHYS logSFR]',fontsize=12)
        plt.ylabel('[MAGPHYS logSFR] / [z0MGS logSFR]',fontsize=12) 
        #plt.ylim(-5,5)
        plt.show()
        
    def mass_hist_oneplot(self, savefig=False):
        
        #magphys_mass = self.magphyscut['logMstar_med']
        magphys_mass = self.altmagphys_cut['combined_logMstar_med']
        
        #remove entries where there is no magphys data available for that galaxy
        err_flag = (self.altmagphys_cut['magphys_flag']) & (magphys_mass!=0)  #ALSO INCLUDE != 0 values. 
                        
                        
        env_mass = [magphys_mass[self.clusflag&err_flag], 
                            magphys_mass[self.rgflag&err_flag], 
                            magphys_mass[self.pgflag&err_flag], 
                            magphys_mass[self.filflag&err_flag], 
                            magphys_mass[self.fieldflag&err_flag]]
                                
        env_names = ['cluster','rich group','poor group','filament','field']
        mybins=np.linspace(7.5,11,12000)
        
        fig = plt.figure(figsize=(8,6))
        
        for i in range(5):
            plt.hist(env_mass[i],bins=mybins,cumulative=True,density=True,linewidth=2,
                     alpha=0.6,histtype='step',label=env_names[i]+f' ({len(env_mass[i])})')

            plt.xlabel(r'log(M*/$M_\odot$)',fontsize=22)
            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)
        plt.xlim(7.5,11.)
        plt.legend(loc='upper left',fontsize=14)
        
        print('K-S p-value (> 0.003 (3sigma), "same distribution"):')
        pairs = list(combinations(env_mass,2))
        names = list(combinations(env_names,2))
        for n in range(len(pairs)):
            pair = pairs[n]
            print(names[n])
            print('%.5f'%(kstest(np.ndarray.tolist(pair[0]),np.ndarray.tolist(pair[1]))[1]))
         
        plt.ylabel('Total Fraction of Galaxies',fontsize=22)
        if savefig==True:
            plt.savefig(homedir+'/Desktop/mass_hist.png',bbox_inches='tight', pad_inches=0.2, dpi=100)
        
        plt.show()
    
    def ndensity_hist_oneplot(self, fullsample=False, savefig=False):
        
        ndensity = self.v2_envcut['n5th']
            
        ndensity_env = [ndensity[self.clusflag],ndensity[self.rgflag], ndensity[self.pgflag],
                        ndensity[self.filflag],ndensity[self.fieldflag]]
        
        if fullsample:
            path_to_dir = homedir+'/Desktop/v2-20220820/'
            v2_env = Table.read(path_to_dir+'vf_v2_environment.fits')
            clusflag = v2_env['cluster_member']
            rgflag = v2_env['rich_group_memb']
            pgflag = v2_env['poor_group_memb']
            filflag = v2_env['filament_member']
            fieldflag = v2_env['pure_field']
            magphys_mass = self.magphys['logMstar_med']
            ndensity = v2_env['n5th']
            ndensity_env = [ndensity[clusflag&(magphys_mass>8.25)], 
                                ndensity[rgflag&(magphys_mass>8.25)],
                                ndensity[pgflag&(magphys_mass>8.25)], 
                                ndensity[filflag&(magphys_mass>8.25)],
                                ndensity[fieldflag&(magphys_mass>8.25)]]                                
        
        env_names = ['cluster','rich group','poor group','filament','field']
        mybins=np.linspace(np.min(ndensity),21,12000)
        
        fig = plt.figure(figsize=(8,6))
        
        for i in range(5):
            plt.hist(ndensity_env[i],bins=mybins,cumulative=True,density=True,linewidth=2,
                     alpha=0.6,histtype='step',label=env_names[i]+f' ({len(ndensity_env[i])})')

            plt.xlabel(r'Volume Density',fontsize=22)
            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)
        plt.xlim(0,20.)
        plt.legend(loc='lower right',fontsize=14)
        
        print('K-S p-value (> 0.003 (3sigma), "same distribution"):')
        pairs = list(combinations(ndensity_env,2))
        names = list(combinations(env_names,2))
        for n in range(len(pairs)):
            pair = pairs[n]
            print(names[n])
            print('%.5f'%(kstest(np.ndarray.tolist(pair[0]),np.ndarray.tolist(pair[1]))[1]))
                
        if savefig==True:
            plt.savefig(homedir+'/Desktop/ndensity_hist.png',bbox_inches='tight', pad_inches=0.2, dpi=100)
        
        plt.show()
    
    def hist_dist_rats(self, savefig=False):
        
        mybins=np.linspace(0,3,40)
        if self.W1:
            xlabels=['', '', '', '', r'R$_{12}$/R$_{3.4}$']
        if not self.W1:
            xlabels=['', '', '', '', r'$R_{12}$/$R_r$']
        colors=['crimson','orange','green','blue','violet']
        labels=['Cluster','Rich Group','Poor Group','Filament','Field']
        
        ratios = self.sizerats.copy()
        clusflag = self.clusflag.copy()
        rgflag = self.rgflag.copy()
        pgflag = self.pgflag.copy()
        filflag = self.filflag.copy()
        fieldflag = self.fieldflag.copy()
        
        re_data = [ratios[clusflag],ratios[rgflag],ratios[pgflag],
                   ratios[filflag],ratios[fieldflag]]
        
        err_upper_bootstrap = []
        err_lower_bootstrap = []
        
        for i in re_data:
            lower_err, upper_err = get_bootstrap_confint(i,bootfunc=np.median,nboot=1000)
            
            err_upper_bootstrap.append(upper_err)
            err_lower_bootstrap.append(lower_err)
        
        fig = plt.figure(figsize=(8,18))
        plt.subplots_adjust(hspace=.2)        
        
        for panel in range(5):
            ax = fig.add_subplot(5,1,panel+1)
            plt.hist(re_data[panel], bins=mybins, color=colors[panel], alpha=0.7, label=labels[panel], density=False, cumulative=False)
            
            med_val = np.median(re_data[panel])
            plt.axvline(med_val,color='black',linestyle='-')
            plt.axvline(err_upper_bootstrap[panel],color='black',linestyle='--')
            plt.axvline(err_lower_bootstrap[panel],color='black',linestyle='--')

            ax.set_xlabel(xlabels[panel],fontsize=20)
            ax.set_ylabel(r'N$_{gal}$',fontsize=20)
            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)
            plt.xlim(-0.05,2.2)
            ax.legend(fontsize=14)
        
        if savefig==True:
            plt.savefig(homedir+'/Desktop/hist_dist.png', bbox_inches='tight', pad_inches=0.2, dpi=100)

        plt.show()
        
        #test statistics
        
        pairs = list(combinations(re_data,2))
        names = list(combinations(labels,2))
        
        print('K-S p-value (> 0.003 (3sigma), "same distribution"):')
        for n in range(len(pairs)):
            pair = pairs[n]
            print(names[n])
            print('%.5f'%(kstest(np.ndarray.tolist(pair[0]),np.ndarray.tolist(pair[1]))[1]))
        
        print()
        print('Mann-Whitney U (null hypothesis: the distribution underlying x is same as the distribution underlying y): ')
        for n in range(len(pairs)):
            pair = pairs[n]
            print(names[n])
            print('%.5f'%(mannwhitneyu(np.ndarray.tolist(pair[0]),np.ndarray.tolist(pair[1]))[1]))
        
        print()
        print('Anderson-Darling (tests the null hypothesis that a sample is drawn from a population that follows a particular distribution): ')
        method = PermutationMethod(n_resamples=9999)
        for n in range(len(pairs)):
            pair = pairs[n]
            print(names[n])
            print(anderson_ksamp([np.ndarray.tolist(pair[0]),np.ndarray.tolist(pair[1])],method=method))
            #print('%.5f'%(anderson_ksamp([np.ndarray.tolist(pair[0]),np.ndarray.tolist(pair[1])])[1]))
    
    def wisesize_mass(self, nbins=5, savefig=False):

        self.data = [self.sizerats[(self.clusflag) | (self.rgflag)], 
                     self.sizerats[(self.pgflag) | (self.filflag) | (self.fieldflag)]]
        self.mass_data = [self.altmagphys_cut['combined_logMstar_med'][(self.clusflag) | (self.rgflag)], 
                          self.altmagphys_cut['combined_logMstar_med'][(self.pgflag) | (self.filflag) | (self.fieldflag)]]
        
        #lastly...some code does not function (pun intended) if I do not explicltly remove masked rows
        self.data[0] = self.data[0][self.mass_data[0]>0]
        self.data[1] = self.data[1][self.mass_data[1]>0]
        self.mass_data[0] = self.mass_data[0][self.mass_data[0]>0]
        self.mass_data[1] = self.mass_data[1][self.mass_data[1]>0]
        
        mass_data = self.mass_data.copy()   #logmass data
        data = self.data.copy()   #Re data
        logmass = self.altmagphys_cut['combined_logMstar_med'].copy()
        
        nbins=nbins
        
        bin_med_clus, bin_edges_clus, binnumber_clus = binned_statistic(mass_data[0],data[0],statistic='median', bins=nbins,range=[(8.6, 11)])
        bin_centers_clus = .5*(bin_edges_clus[:-1]+bin_edges_clus[1:])
        
        bin_med_fall, bin_edges_fall, binnumber_fall = binned_statistic(mass_data[1],data[1],statistic='median', bins=nbins,range=[(8.6, 11)])
        bin_centers_fall = .5*(bin_edges_fall[:-1]+bin_edges_fall[1:])
        
        bin_bounds_clus = []
        bin_bounds_fall = []
        for bound in range(nbins):
            bin_bounds_clus.append([bin_edges_clus[bound],bin_edges_clus[bound+1]])
            bin_bounds_fall.append([bin_edges_fall[bound],bin_edges_fall[bound+1]])
        
        err_clus = []
        err_fall = []
        for bound in bin_bounds_clus:
            if self.MeanMedian=='mean':
                bound_flag_clus = (mass_data[0]>bound[0])&(mass_data[0]<=bound[1])
                lower_clus, upper_clus = get_bootstrap_confint(data[0][bound_flag_clus],bootfunc=np.mean,nboot=1000)
                err_clus.append([lower_clus,upper_clus])
            if self.MeanMedian=='median':
                bound_flag_clus = (mass_data[0]>bound[0])&(mass_data[0]<=bound[1])
                lower_clus, upper_clus = get_bootstrap_confint(data[0][bound_flag_clus],bootfunc=np.median,nboot=1000)
                err_clus.append([lower_clus,upper_clus])
        for bound in bin_bounds_fall:
            if self.MeanMedian=='mean':
                bound_flag_fall = (mass_data[1]>bound[0])&(mass_data[1]<=bound[1])
                lower_fall, upper_fall = get_bootstrap_confint(data[1][bound_flag_fall],bootfunc=np.mean,nboot=1000)
                err_fall.append([lower_fall,upper_fall])
            if self.MeanMedian=='median':
                bound_flag_fall = (mass_data[1]>bound[0])&(mass_data[1]<=bound[1])
                lower_fall, upper_fall = get_bootstrap_confint(data[1][bound_flag_fall],bootfunc=np.median,nboot=1000)
                err_fall.append([lower_fall,upper_fall])

        plt.figure(figsize=(8,6))
        
        plt.scatter(mass_data[0],data[0],color='crimson',s=15,alpha=0.2,label='Cluster & RG',zorder=1)
        plt.scatter(mass_data[1],data[1],color='blue',s=15,alpha=0.2,label='PG, Filament, & Field',zorder=1)
        plt.scatter(bin_centers_clus, bin_med_clus, color='crimson', s=250, edgecolors='black', label='<Cluster & RG>',zorder=3)
        plt.scatter(bin_centers_fall, bin_med_fall, color='blue', s=250, edgecolors='black', label='<PG, Filament, & Field>',zorder=3)

        for n in range(nbins):
            plt.plot([bin_centers_clus[n],bin_centers_clus[n]], [err_clus[n][0],err_clus[n][1]],color='crimson',zorder=2)
            #create lower, upper caps on errorbars
            plt.plot([bin_centers_clus[n]-0.08,bin_centers_clus[n]+0.08],[err_clus[n][0],err_clus[n][0]],color='crimson',zorder=2)  
            plt.plot([bin_centers_clus[n]-0.08,bin_centers_clus[n]+0.08],[err_clus[n][1],err_clus[n][1]],color='crimson',zorder=2)

            plt.plot([bin_centers_fall[n],bin_centers_fall[n]], [err_fall[n][0],err_fall[n][1]],color='blue',zorder=2)
            #create lower, upper caps on errorbars
            plt.plot([bin_centers_fall[n]-0.08,bin_centers_fall[n]+0.08],[err_fall[n][0],err_fall[n][0]],color='blue',zorder=2)
            plt.plot([bin_centers_fall[n]-0.08,bin_centers_fall[n]+0.08],[err_fall[n][1],err_fall[n][1]],color='blue',zorder=2)
            
        plt.ylabel(r'R$_{12}$/R$_{3.4}$',fontsize=18)
        plt.xlabel(r'log$_{10}$(M$_*$/M$_\odot$)',fontsize=18)
        
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        
        plt.title('Median Re Size Ratios vs. Stellar Mass',fontsize=20)
        plt.ylim(0.5,1)
        plt.xlim(8.5,11)
        plt.legend(fontsize=14)
        
        if savefig==True:
            plt.savefig(homedir+'/Desktop/wisesize_comp_mass_env.png', bbox_inches='tight', pad_inches=0.2, dpi=100)

        plt.show()    
    
if __name__ == '__main__':
    '''
    print("""USAGE:
    ***NOTE: The mass-matching flag, which appears every so often, dictates whether the subsample is 
        'tailored' to the masses of the cluster galaxies.***"""
    print('-----------------------------------------------------')
    
    print("""
    cat = catalogs(conv=False,MeanMedian='mean',MADmultiplier=5,cutAGN=False,W1=False) --> initiate 
        catalog class. 
        MeanMedian will propagate to all plots. MADmultip[lier is a proxy for the extent to which a 
        galaxy must be an 'outlier' in order to be removed from calculations.
        cutAGN is precisely that and will remove flagged AGN from the plots below (both WISE and BPT AGN)
        W1 switches r-band with WISE 3.4-micron GALFIT results.
    cat.inclination_check(savefig=False) --> checks for size ratio, SFR, and Mstar dependencies on
        optical axis ratio (found using GALFIT 1-Sersic profile fits on r-band images of the galaxies).
    cat.mass_completeness_check(mass_match=False,savefig=False) --> compares mass distribution of galaxies 
        for full VFS and the cut sample.
    cat.sfrmstar_magphys(show_HI=False,show_sizerat=True,show_D25=False,savefig=False) --> generates
        SFR v. Mstar plot, color-coded according to either available HI gas mass, D25 (Hyperleda), 
        or GALFIT size ratios. SFR and Mstar values from MAGPHYS.
    cat.sfrmstar_z0mgs(show_HI=False,show_sizerat=True,savefig=False) --> same as above, but with z0mgs.
    cat.r12_vs_rstar(savefig=False) --> plots R_12 vs. R_r effective radii for the full
        subsample as well as according to each environment bin (six total panels).
    cat.r12_vs_rstar_v2(savefig=False) --> plots R_12 vs. R_r effective radii for the full
        subsample as well as according to each environment bin (5 env columns, 2 mass ranges). 
    cat.envbins(savefig=False) --> plots number of subsample galaxies in each environment bin
    cat.ratio_MS(showHI=False,savefig=False) --> plots size ratios as a function of MS distance, 
        colored according to normalized HI mass in the case where showHI=True.
    cat.hist_dist_rats(savefig=False) --> size ratio distributions for each of the five defined 
        environment bins (plt.subplots with one column, five rows). Rainbow colors, very aesthestic.
        Outputs a lot of statistical data, very fun.
    cat.env_means(mass_match=False, trimOutliers=False, errtype='bootstrap', W1=False, 
                  r90=False, savefig=False) --> plots either mean or median size ratio (w3/r or W3/W1) 
                  in each environment bin; trimOutliers will output an additional plot which compares my 
                  no PSF parameters to Rose's parameters, allowing the user to visualize which points are 
                  omitted in the trimmed env_means plot; errtype either 'bootstrap' or 'std_err;' r90
                  simply switches r50 data for, well, r90.
    cat.env_means_comp(savefig=True) --> same as above, but with projected WISESize uncertainties
    cat.sfr_hist(z0mgs_comp=True,savefig=False) --> generate sfr histogram subplots per environment bin; 
        will compare MAGPHYS SFRs with z0mgs values if True
    cat.mass_hist_oneplot(savefig=False) --> generate same MAGPHYS mass histograms 
        but as a single plot; fullsample=True will use all VFS galaxies.
    cat.ndensity_hist_oneplot(fullsample=False, savefig=False) --> same as above, but with volume
        density (unsure of the units presently).
    cat.wisesize_mass(nbins=3,savefig=False) --> creates median/mean size ratios vs. mass bin plot for two
        environment classifications: low-density (PG, filament, field) and high-density (cluster, RG).
    """)
    print('-----------------------------------------------------')
    print()
    '''