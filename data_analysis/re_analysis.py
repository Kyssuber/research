'''
GOAL: analyze GALFIT output parameters for Virgowise subsample galaxies.

pixel_scale = {'FUV':1.5,'NUV':1.5,'g':0.262,'r':0.262,'z':0.262,'W1':2.75,'W2':2.75,'W3':2.75,'W4':2.75}
psf_oversampling = {'FUV':1,'NUV':1,'g':1,'r':1,'g':1,'W1':1,'W2':1,'W3':1,'W4':1}
mag_zeropoint = {'FUV':22.5,'NUV':22.5,'g':22.5,'r':22.5,'g':22.5,'W1':22.5,'W2':22.5,'W3':22.5,'W4':22.5}
image_resolution = {'FUV':6,'NUV':6,'g':1.5,'r':1.5,'z':1.5,'W1':6.1,'W2':6.4,'W3':6.5,'W4':12}
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
    
    def __init__(self,conv=False,MeanMedian='mean',MADmultiplier=5,cutAGN=False):
        self.MADmultiplier = MADmultiplier
        path_to_dir = homedir+'/Desktop/v2-20220820/'
        self.v2_env = Table.read(path_to_dir+'vf_v2_environment.fits')
        self.v2_main = Table.read(homedir+'/v2_snrcoadd.fits')
        self.magphys = Table.read(path_to_dir+'vf_v2_magphys_07-Jul-2022.fits')
        self.z0mgs = Table.read(path_to_dir+'vf_v2_z0mgs.fits')
        self.HI_tab = Table.read(path_to_dir+'vf_v2_CO_HI.fits')
        self.hyp_tab = Table.read(path_to_dir+'vf_v2_hyperleda.fits')
        self.sgaparams = Table.read(homedir+'/sgacut_SGAparams.fits')
        self.sgaparams.sort('VFID_1')   #sort in ascending order all rows according to VFID
        
        self.MeanMedian = MeanMedian  #whether I plot median or mean size ratios for the self.env_means() figure
        
        self.conv = conv
        self.cutAGN = cutAGN
        print('Convolution: ',str(self.conv))
        if self.cutAGN:
            print('AGN cut from sample.')
        
        if not self.conv:
            self.rdat = Table.read(homedir+'/output_params_r_nopsf.fits')
            self.w3dat = Table.read(homedir+'/output_params_W3_nopsf.fits')
            self.kimparams = Table.read(homedir+'/kimparams_nopsf.fits')   #my nopsf galfit parameters
        if self.conv:
            self.rdat = Table.read(homedir+'/output_params_r_psf.fits')
            self.w3dat = Table.read(homedir+'/output_params_W3_psf.fits')
            self.kimparams = Table.read(homedir+'/kimparams_psf.fits')   #my psf galfit parameters
    
        self.roseparams = self.w3dat.copy()        
        self.cut_cats()
    
    def cut_cats(self):
        subsample_flag = self.v2_main['sgacut_flag']
        
        self.v2_env = self.v2_env[subsample_flag]
        self.v2_maincut = self.v2_main[subsample_flag]
        self.magphyscut = self.magphys[subsample_flag]
        self.z0mgs = self.z0mgs[subsample_flag]
        self.HI_tab = self.HI_tab[subsample_flag]
        self.hyp_tab = self.hyp_tab[subsample_flag]
        
        self.re_rband = self.rdat['re']
        self.re_w3band = self.w3dat['re']

        self.PA_rband = self.rdat['PA']
        self.PA_w3band = self.w3dat['PA']

        self.BA_rband = self.rdat['BA']
        self.BA_w3band = self.w3dat['BA']

        fail_flag = (self.re_rband == 0.0) | (self.re_w3band == 0.0)
        r_flag = (self.re_rband == 0.0)
        w3_flag = (self.re_w3band == 0.0)
        err_flag = (self.rdat['err_flag']==1) | (self.w3dat['err_flag']==1)
        
        n_fails_w3 = len(self.re_rband[r_flag])
        n_fails_r = len(self.re_w3band[w3_flag])
        
        self.cut_flags = (~fail_flag)&(~err_flag)
        
        #apply optional AGN cut
        if self.cutAGN:
            #WISE color magnitudes
            wise_mag_cut=Table.read(homedir+'/Desktop/v2-20220820/vf_v2_unwise.fits')[subsample_flag]
            #spectral line strengths
            bpt_lines_cut=Table.read(homedir+'/Desktop/v2-20220820/vf_v2_nsa_v0.fits')[subsample_flag]
            wise_agn = (wise_mag_cut['w1_mag'] - wise_mag_cut['w2_mag']>0.65) & (wise_mag_cut['w2_mag']-wise_mag_cut['w3_mag']<4.0)
            agn_kauffman = (np.log10(bpt_lines_cut['O3FLUX']/bpt_lines_cut['HBFLUX']) > (.61/(np.log10(bpt_lines_cut['N2FLUX']/bpt_lines_cut['HAFLUX']-.05))+1.3)) | (np.log10(bpt_lines_cut['N2FLUX']/bpt_lines_cut['HAFLUX']) > 0.)
        
            self.wise_agn_flag = np.asarray(wise_agn)
            self.agn_kauffman_flag = np.asarray(agn_kauffman)
            AGN_flags = (self.wise_agn_flag)|(self.agn_kauffman_flag)
            
            print(f'# WISE AGN in VF subsample: {len(self.v2_maincut[self.wise_agn_flag])}')
            print(f'# BPT AGN in VF subsample: {len(self.v2_maincut[self.agn_kauffman_flag])}')

            self.cut_flags = (~AGN_flags) & (~fail_flag) & (~err_flag)
            print(f'Number of galaxies flagged with AGN and a GALFIT error: {len(self.v2_env[(AGN_flags) & (err_flag)])}')
        
        self.re_rband_cut = self.re_rband[self.cut_flags]
        self.re_w3band_cut = self.re_w3band[self.cut_flags]

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
        self.sgaparams_cut = self.sgaparams[self.cut_flags]
        
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
    
    def mass_matching(self):
        
        #select one cluster galaxy
        #create pool of field galaxies with logmass within 0.3 dex relative to the one cluster galaxy
        #from this pool, randomly draw int(len(field)/len(cluster)) galaxies
            #try without replacement first.
        #repeat (2) and (3) for rich group, filament, and poor group galaxies
        #repeat all for every cluster galaxy
        
        n_dex = 0.3
        
        n_to_draw = [int(round(len(self.v2_envcut[self.fieldflag])/len(self.v2_envcut[self.clusflag]),0)),
                     int(round(len(self.v2_envcut[self.rgflag])/len(self.v2_envcut[self.clusflag]),0)),
                     int(round(len(self.v2_envcut[self.pgflag])/len(self.v2_envcut[self.clusflag]),0)),
                     int(round(len(self.v2_envcut[self.filflag])/len(self.v2_envcut[self.clusflag]),0))]
        
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
            field_pool = m_field[(m_field['logMstar'] < (m_clus['logMstar'][index]+n_dex)) & (m_field['logMstar'] > (m_clus['logMstar'][index]-n_dex))]
            rg_pool = m_rg[(m_rg['logMstar'] < (m_clus['logMstar'][index]+n_dex)) & (m_rg['logMstar'] > (m_clus['logMstar'][index]-n_dex))]
            pg_pool = m_pg[(m_pg['logMstar'] < (m_clus['logMstar'][index]+n_dex)) & (m_pg['logMstar'] > (m_clus['logMstar'][index]-n_dex))]
            fil_pool = m_fil[(m_fil['logMstar'] < (m_clus['logMstar'][index]+n_dex)) & (m_fil['logMstar'] > (m_clus['logMstar'][index]-n_dex))]
            
            #the 'pools' of galaxies from which the mass-matched samples will be drawn
            pools = [field_pool, rg_pool, pg_pool, fil_pool]
            
            for env in range(4):
                if n_to_draw[env] <= len(pools[env]):
                    choices = np.random.choice(pools[env]['VFID'],n_to_draw[env])
                else:    
                    try:
                        choices = np.random.choice(pools[env]['VFID'],n_to_draw[env]-1)
                    except ValueError:
                        choices = np.random.choice(pools[env]['VFID'],0)
                for gal in choices:
                    comp_galaxies[env].append(gal)
        
        print('n_field, n_rg, n_pg, n_filament:',len(comp_galaxies[0]),len(comp_galaxies[1]),len(comp_galaxies[2]),len(comp_galaxies[3]))
        print('n_clus:',len(m_clus))
        
        #once process is confirmed to be correct, I could replace comp_galaxies[env] arrays with T, F as opposed to 
        #just inserting the galaxy VFID.
        #define the result as self.comp_galaxies, add flag kwarg to certain functions below, then I can re-generate
        #certain plots (e.g., cat.env_means()). 
    
    def inclination_check(self, savefig=False):
        
        #aim: check dependence of size ratio, SFR, Mstar on optical axis ratio
        #also looks for whether trends exist between size ratio and Sérsic index or T-type
        print('Spearman rank test: -1<p<1, with 0 indicating no correlation.')
        
        #remove entries where there is no magphys data available for that galaxy        
        err_flag_cut = (self.magphyscut['magphysFlag'])
        
        #cut logsfr, logmass
        logsfr_cut = self.magphyscut['logSFR'][err_flag_cut]
        logmass_cut = self.magphyscut['logMstar'][err_flag_cut]
        
        #cut size ratio, optical axis ratio arrays
        sizerats = self.sizerats[err_flag_cut]
        BA = self.BA_rband_cut[err_flag_cut]
        
        x = [BA, BA[(logsfr_cut>-6)], BA[(logmass_cut>6)]]
        y = [sizerats, logsfr_cut[(logsfr_cut>-6)], logmass_cut[(logmass_cut>6)]]
        colors = ['crimson','green','blue']
        y_labels = [r'R$_{12}$/R$_r$', 'logSFR', r'log$M_*$']
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
            plt.savefig(homedir+'/Desktop/inclination_check.png', dpi=150, bbox_inches='tight', pad_inches=0.2)

        plt.show()
        
        vfid_list = self.v2_maincut['VFID'][err_flag_cut]
        print(f'Size ratio outliers: {vfid_list.data[(sizerats>2)]}, \n logSFR outliers: {vfid_list.data[(logsfr_cut<-6)]}, \n logMstar outliers: {vfid_list.data[(logmass_cut<6)]}')
        
        #also looks for whether trends exist between size ratio and Sérsic index or T-type
        #ttype_dict = {'Sa':1,'Sab':2,'Sb':3,'Sbc':4,'Sc':5,'Scd':6,'Sd':7,'Sdm':8,'Sm':9,'Ir':10}
        
        y2 = self.sizerats
        x2 = [self.hyp_tab_cut['t'], self.rdat['nsersic'][self.cut_flags]]
        x_labels2 = ['Morphological T-Type', 'Sérsic Index (r-band)']
        colors2 = ['orangered','blueviolet']
        
        fig=plt.figure(figsize=(20,4))
        plt.subplots_adjust(hspace=.1,wspace=.2)
        
        for index in range(2):
            ax = fig.add_subplot(1, 2, index+1)
            plt.scatter(x2[index][y2<20], y2[y2<20], color=colors2[index], s=20)
            plt.xlabel(x_labels2[index],fontsize=18)
            plt.ylabel(r'R$_{12}$/R$_r$',fontsize=18)
            
            #spearman-rank test
            sr = spearmanr(x2[index],y2)
            rho = sr[0]
            p = sr[1]
            plt.plot([],[],color="None",label=f'rho={np.round(rho,4)}')
            plt.plot([],[],color="None",label=f'p-value={np.round(p,4)}')
            
            plt.yscale('log')
            plt.legend(handlelength=0,fontsize=17)
        
        if savefig:
            plt.savefig(homedir+'/Desktop/morphology_check.png', dpi=150, bbox_inches='tight', pad_inches=0.2)

        plt.show()
        
        
        
        
        
        
        
        
        
    def mass_completeness_check(self, savefig=False):
        
        #aim: Make histogram of full VFS sample with quality flags. Also make histogram with T-type and 12um SNR cut. Take the ratio of these histograms (WISESize sample/VFS) and plot them vs. mass. Below some mass this ratio will decline. We can use that point as our mass completeness limit.
        
        mass_bins = np.arange(6,11.5,0.15)
        
        logmass = self.magphys['logMstar']
        
        #remove entries where there is no magphys data available for that galaxy
        err_flag = (self.magphys['magphysFlag'])
        
        logmass_cut = logmass[self.v2_main['sgacut_flag'] & (err_flag)]        
        logmass = logmass[err_flag]
        
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
            plt.savefig(homedir+'/Desktop/mass_completeness_hist.png', dpi=150, bbox_inches='tight', 
                        pad_inches=0.2)
    
    
    
    
    
    def sfrmstar_magphys(self, show_HI=False, show_D25=False, show_sizerat=True, savefig=False):
        
        #prepare D25 data for colorbar
        d25 = 10**(self.hyp_tab['logd25'])
        
        #prepare MHI_to_Mstar data for colorbar
        MHI_to_Mstar = self.HI_tab['MHI']/(10**self.z0mgs['logmass'])
        
        for n in range(len(MHI_to_Mstar)):
            if self.HI_tab['MHI'].mask[n]:     #if value is masked, set to be -999
                MHI_to_Mstar[n] = -999
            else:
                MHI_to_Mstar[n] = MHI_to_Mstar[n] if ((np.log(self.HI_tab['MHI'][n])>0) & (MHI_to_Mstar[n]<10) & (MHI_to_Mstar[n]!=1)) else -999
        
        logsfr = self.magphys['logSFR']
        logmass = self.magphys['logMstar']
        
        #remove entries where there is no magphys data available for that galaxy
        err_flag = (self.magphys['magphysFlag'])
        err_flag_cut = (self.magphyscut['magphysFlag'])
        
        logsfr_cut = logsfr[(self.v2_main['sgacut_flag']) & (err_flag)]
        logmass_cut = logmass[self.v2_main['sgacut_flag'] & (err_flag)]
        MHI_to_Mstar_cut = MHI_to_Mstar[(err_flag[self.v2_main['sgacut_flag']])]
        d25_cut = d25[(err_flag[self.v2_main['sgacut_flag']])]
        
        logsfr = logsfr[err_flag]
        logmass = logmass[err_flag]

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
            plt.scatter(self.magphyscut['logMstar'][err_flag_cut],self.magphyscut['logSFR'][err_flag_cut], 
                        c=self.sizerats[err_flag_cut], cmap='viridis', s=60, alpha=0.6, label='VF Subsample')
            cb = plt.colorbar()
            cb.set_label(label=r'R$_{12}$/R$_r$',size=23)
            plt.clim(0.1,1)
            
        if show_D25:
            plt.scatter(logmass_cut,logsfr_cut, c=d25_cut, cmap='viridis', s=60, 
                        alpha=0.6, label='VF Subsample')
            cb = plt.colorbar()
            cb.set_label(label=r'Optical D25',size=23)
            plt.clim(15,40)
            
        #add z0mgs main sequence line
        #y = logmass*0.47 - 4.88
        #plt.plot(logmass,y,label=('Main Sequence'),color='black',alpha=0.6)
        
        #fit main-sequence line to data
        p1 = np.poly1d(np.polyfit(logmass[(logsfr-logmass)>-11.5],logsfr[(logsfr-logmass)>-11.5],1))
        #p2 = np.poly1d(np.polyfit(logmass[(logmass-logsfr>-12)],logsfr[(logmass-logsfr>-12)],2))

        x_p = np.linspace(np.min(logmass),np.max(logmass),1000)

        plt.plot(x_p,p1(x_p),color='black',label='MS fit',alpha=0.7,linewidth=3)
        #plt.plot(x_p,p2(x_p),color='blue',label='2d fit')
        print(f'p1 parameters: m={np.round(p1[1],2)}, b={np.round(p1[0],2)}')
        
        #plt.plot(logmass, -11.5+logmass, color='red', linestyle='-', alpha=0.5, label='log(sSFR)>-11 limit')
        
        cb.ax.tick_params(labelsize=15)
        
        plt.xlabel(r'log(M*/$M_\odot$)',fontsize=22)
        plt.ylabel(r'log(SFR/($M_\odot$/yr))',fontsize=22)
        
        plt.xlim(6,11.6)
        plt.ylim(-5,2.3)
        
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        
        plt.legend(fontsize=12)

        if savefig==True:
            plt.savefig(homedir+'/Desktop/sfrmstar_magphys.png', dpi=150, bbox_inches='tight', pad_inches=0.2)
        
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
            plt.scatter(self.magphyscut['logMstar'][err_flag_cut],self.magphyscut['logSFR'][err_flag_cut], 
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
            plt.savefig(homedir+'/Desktop/sfrmstar_z0mgs.png', dpi=150, bbox_inches='tight', pad_inches=0.2)
        
        plt.show()
    
    
    def r12_vs_rstar(self, sfr_mstar='mstar', savefig=False):
        
        logsfr = self.magphyscut['logSFR']
        logmass = self.magphyscut['logMstar']
        
        #remove entries where there is no magphys data available for that galaxy
        err_flag = (self.magphyscut['magphysFlag'])
        
        logsfr = logsfr[err_flag]
        logmass = logmass[err_flag]
        
        r_arcsec = self.re_rband_cut[(self.v2_maincut['VFID']!='VFID1984')&(err_flag)]*0.262
        w3_arcsec = self.re_w3band_cut[(self.v2_maincut['VFID']!='VFID1984')&(err_flag)]*2.75
        
        rr_dat = [r_arcsec,
                  r_arcsec[self.clusflag[(self.v2_maincut['VFID']!='VFID1984')&(err_flag)]],
                  r_arcsec[self.rgflag[(self.v2_maincut['VFID']!='VFID1984')&(err_flag)]],
                   r_arcsec[self.pgflag[(self.v2_maincut['VFID']!='VFID1984')&(err_flag)]],
                  r_arcsec[self.filflag[(self.v2_maincut['VFID']!='VFID1984')&(err_flag)]],
                    r_arcsec[self.fieldflag[(self.v2_maincut['VFID']!='VFID1984')&(err_flag)]]]
        w3_dat = [w3_arcsec,
                  w3_arcsec[self.clusflag[(self.v2_maincut['VFID']!='VFID1984')&(err_flag)]],
                  w3_arcsec[self.rgflag[(self.v2_maincut['VFID']!='VFID1984')&(err_flag)]],
                   w3_arcsec[self.pgflag[(self.v2_maincut['VFID']!='VFID1984')&(err_flag)]],
                  w3_arcsec[self.filflag[(self.v2_maincut['VFID']!='VFID1984')&(err_flag)]],
                  w3_arcsec[self.fieldflag[(self.v2_maincut['VFID']!='VFID1984')&(err_flag)]]]
        
        logmass_dat = [logmass,
                       logmass[self.clusflag[(self.v2_maincut['VFID']!='VFID1984')&(err_flag)]],
                       logmass[self.rgflag[(self.v2_maincut['VFID']!='VFID1984')&(err_flag)]],
                       logmass[self.pgflag[(self.v2_maincut['VFID']!='VFID1984')&(err_flag)]],
                       logmass[self.filflag[(self.v2_maincut['VFID']!='VFID1984')&(err_flag)]],
                       logmass[self.fieldflag[(self.v2_maincut['VFID']!='VFID1984')&(err_flag)]]]
        
        logsfr_dat = [logsfr,
                       logsfr[self.clusflag[(self.v2_maincut['VFID']!='VFID1984')&(err_flag)]],
                       logsfr[self.rgflag[(self.v2_maincut['VFID']!='VFID1984')&(err_flag)]],
                       logsfr[self.pgflag[(self.v2_maincut['VFID']!='VFID1984')&(err_flag)]],
                       logsfr[self.filflag[(self.v2_maincut['VFID']!='VFID1984')&(err_flag)]],
                       logsfr[self.fieldflag[(self.v2_maincut['VFID']!='VFID1984')&(err_flag)]]]
        
        labels = ['Full Subsample','Cluster','Rich Group','Poor Group','Filament','Field']
        
        colors=['purple','crimson','black','magenta','green','blue']
        
        fig = plt.figure(figsize=(30,16))
        plt.subplots_adjust(hspace=.2,wspace=.3)
        
        for n in range(1,7):
            ax = fig.add_subplot(2,3,n)
            if sfr_mstar=='mstar':
                plt.scatter(rr_dat[n-1],w3_dat[n-1],alpha=1,c=logmass_dat[n-1],edgecolors='black',
                            s=150,cmap='viridis')
            if sfr_mstar=='sfr':
                plt.scatter(rr_dat[n-1],w3_dat[n-1],alpha=0.9,c=logsfr_dat[n-1],edgecolors='black',
                            s=150,cmap='plasma')
            plt.axline((0, 0), slope=1, color='indigo')
            plt.text(.05, .95, labels[n-1], ha='left', va='top',transform=ax.transAxes,fontsize=25)
            
            plt.xlim(0.7,800)
            plt.ylim(0.7,340)
            
            plt.xscale('log')
            plt.yscale('log')

            plt.xlabel(r'R$_{r}$ (arcsec)',fontsize=25)
            plt.ylabel(r'R$_{12}$ (arcsec)',fontsize=25)
        
            plt.xticks(fontsize=20)
            plt.yticks(fontsize=20)
            
            cb = plt.colorbar()
            if sfr_mstar=='mstar':
                cb.set_label(label=r'log(M$_{star}$/M$_{\odot}$)',size=25)
                plt.clim(8.5,10)
            if sfr_mstar=='sfr':
                cb.set_label(label=r'log(SFR/($M_\odot$/yr))',size=25)
                plt.clim(-1,0.6)
            cb.ax.tick_params(labelsize=20)
        
        if savefig==True:
            plt.savefig(homedir+'/Desktop/r12_rstar.png', dpi=200, bbox_inches='tight', pad_inches=0.2)
        
        plt.show()
        
    def ratio_MS(self, showHI=False, savefig=False):
        
        m = 0.86
        b = -8.79
        
        #prepare MHI_to_Mstar data for colorbar
        MHI_to_Mstar = self.HI_tab_cut['MHI']/(10**self.magphyscut['logMstar'])
        
        logsfr = self.magphyscut['logSFR']
        logmass = self.magphyscut['logMstar']
        
        #fit main-sequence line to data (1D --> p1; parameters taken from sfrmstar_magphys)
        Mstar_full = logmass
        
        #main sequence best-fit line
        #slope: 0.47 +/- 0.01
        #y-intercept: -4.88 +/- 0.10
        #Mstar_full = logmass
        #y = Mstar_full*0.47 - 4.88
        
        #remove entries where there is no magphys data available for that galaxy
        err_flag_cut = (self.magphyscut['magphysFlag'])

        #Now calculate distance of each point to the best-fit line
        #Distance = (| a*x1 + b*y1 + c |) / (sqrt( a*a + b*b))
        dist = (m*logmass[err_flag_cut] + logsfr[err_flag_cut] + b) / (np.sqrt((m)**2 + (1)**2))
        
        MHI_to_Mstar_cut = MHI_to_Mstar[err_flag_cut]
        
        plt.figure(figsize=(10,6))
        plt.axhline(1,linestyle='--',color='r',alpha=0.4)
        
        if showHI:
            plt.scatter(dist,self.sizerats[err_flag_cut],color='gray',s=10,alpha=0.3)   #all points
            plt.scatter(dist,self.sizerats[err_flag_cut],c=MHI_to_Mstar_cut)   #will only plot points with MHI data
            plt.colorbar().set_label(label='MHI_to_Mstar',size=15)
            plt.clim(0,1)
        if not showHI:
            plt.scatter(dist,self.sizerats[err_flag_cut],color='blue',s=40,alpha=0.7)   #all points
        
        #plt.yscale('log')
        plt.ylim(0, 1.75)   #artificially trim outliers with size ratios>5 (there is indeed one with ratio~100)
        plt.xlim(-3.8,1.5)
        
        plt.xlabel(r'Distance from Main Sequence Line',fontsize=17)
        plt.ylabel(r'Size Ratio ($R_{12}/R_{r}$)',fontsize=17)
        
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        
        print(np.polyfit(dist,self.sizerats[err_flag_cut],1))
        
        if savefig==True:
            plt.savefig(homedir+'/Desktop/MS_dist.png', dpi=200, bbox_inches='tight', pad_inches=0.2)
        
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

        if savefig==True:
            plt.savefig(homedir+'/Desktop/envbins.png', dpi=200)
        
        plt.show()

    def env_means(self, trimOutliers=False, combine_mid=False, errtype='bootstrap', savefig=False):    
        
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
        err_upper_bootstrap = []
        err_lower_bootstrap = []
        
        index = np.arange(1,6,1)
        env_names = ['Cluster','Rich \n Group','Poor \n Group','Filament','Field']
        
        for j,i in enumerate(re_data):      #j==index, i==value
            if self.MeanMedian=='mean':
                central_pts.append(np.mean(i))
                lower_err, upper_err = get_bootstrap_confint(i,bootfunc=np.mean,nboot=100)

            if self.MeanMedian=='median':
                central_pts.append(np.median(i))
                lower_err, upper_err = get_bootstrap_confint(i,bootfunc=np.median,nboot=100)
            
            err.append(np.std(i)/np.sqrt(len(i)))
            err_upper_bootstrap.append(upper_err)
            err_lower_bootstrap.append(lower_err)
            
        err_color = 'orangered'
        plt.figure(figsize=(10,6))
        plt.scatter(index,central_pts,color='blue',s=50,zorder=2,edgecolors='black',label=self.MeanMedian)
       
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
    
        plt.ylim(0.6,.85)
        
        #plt.legend(fontsize=15)
        
        if savefig==True:
            plt.savefig(homedir+'/Desktop/env_ratio.png', bbox_inches='tight', pad_inches=0.2, dpi=200)
            
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

        env_names = ['Cluster','Rich Group','Poor Group','Filament','Field']

        mybins=np.linspace(7.5,11.5,1000)
        
        fig = plt.figure(figsize=(11,8))
        plt.subplots_adjust(hspace=.4,wspace=.3)
        
        for i in range(1,6):
            ax=fig.add_subplot(2,3,i)
            if i == 4:
                ax.set_position([0.24,0.125,0.228,0.3])
            if i == 5:
                ax.set_position([0.55,0.125,0.228,0.3])
            
            plt.hist(magphys_env_mass[i-1],bins=mybins,alpha=0.8,cumulative=True,density=True,
                     histtype='step',label='Magphys')
            
            if z0mgs_comp==True:
                plt.hist(z0mgs_env_mass[i-1],bins=mybins,alpha=0.8,cumulative=True,density=True,
                         histtype='step',label='z0mgs')
            
            if (i == 1)&(z0mgs_comp==True):
                plt.legend(fontsize=13,loc='upper left')
            plt.xlabel(r'log(M*/$M_\odot$)',fontsize=22)
            plt.xlim(7.5,11.5)
            plt.title(env_names[i-1],fontsize=22)
            
            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)
        
        print('K-S p-value (> 0.003 (3sigma), "same distribution"):')
        for n in range(len(env_names)):
            print(env_names[n])
            print('%.5f'%(kstest(z0mgs_env_mass[n],magphys_env_mass[n])[1]))
        
        if savefig==True:
            plt.savefig(homedir+'/Desktop/mass_hist.png',bbox_inches='tight', pad_inches=0.2, dpi=200)
        
        plt.show()
 
    def mass_hist_oneplot(self, savefig=False):
        
        magphys_mass = self.magphyscut['logMstar']
        
        #remove entries where there is no magphys data available for that galaxy
        err_flag = (self.magphyscut['magphysFlag'])
    
        magphys_env_mass = [magphys_mass[self.clusflag&err_flag], magphys_mass[self.rgflag&err_flag], 
                            magphys_mass[self.pgflag&err_flag], 
                            magphys_mass[self.filflag&err_flag], magphys_mass[self.fieldflag&err_flag]]

        env_names = ['cluster','rich group','poor group','filament','field']
        mybins=np.linspace(7.5,11,12000)
        
        fig = plt.figure(figsize=(8,6))
        
        for i in range(5):
            plt.hist(magphys_env_mass[i],bins=mybins,cumulative=True,density=True,
                     histtype='step',label=env_names[i])

            plt.xlabel(r'log(M*/$M_\odot$)',fontsize=22)
            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)
        plt.xlim(7.5,11.)
        plt.legend(loc='upper left',fontsize=16)
        
        print('K-S p-value (> 0.003 (3sigma), "same distribution"):')
        pairs = list(combinations(magphys_env_mass,2))
        names = list(combinations(env_names,2))
        for n in range(len(pairs)):
            pair = pairs[n]
            print(names[n])
            print('%.5f'%(kstest(pair[0],pair[1])[1]))
                
        if savefig==True:
            plt.savefig(homedir+'/Desktop/mass_hist.png',bbox_inches='tight', pad_inches=0.2, dpi=200)
        
        plt.show()
    
    def hist_dist_rats(self, savefig=False):
        
        mybins=np.linspace(0,3,30)
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
            plt.xlim(-0.1,2)
            ax.legend(fontsize=15)
        
        if savefig==True:
            plt.savefig(homedir+'/Desktop/hist_dist.png', bbox_inches='tight', pad_inches=0.2, dpi=200)

        plt.show()
        
        #ks-test statistics
        
        print('K-S p-value (> 0.003 (3sigma), "same distribution"):')
        pairs = list(combinations(re_data,2))
        names = list(combinations(labels,2))
        for n in range(len(pairs)):
            pair = pairs[n]
            print(names[n])
            print('%.5f'%(kstest(pair[0],pair[1])[1]))
    
    def compareSGA(self,savefig=False):
        
        r50_sga_r = self.sgaparams_cut['r50 (rband)']  #arcsec
        r50_gal_r = self.re_rband_cut.copy()*0.262 #arcsec
        
        plt.figure(figsize=(8,6))
        plt.scatter(r50_sga_r,r50_sga_r/r50_gal_r,color='crimson',edgecolors='black',alpha=0.4)
        
        plt.xlabel(r'r50$_{SGA}$',fontsize=18)
        plt.ylabel(r'r50$_{SGA}$/r50$_{GALFIT}$',fontsize=18)
        if self.conv==True:
            plt.title('r-band r50 Comparison (PSF)',fontsize=20)
        if self.conv==False:
            plt.title('r-band r50 comparison (noPSF)',fontsize=20)
        
        m, b = np.polyfit(r50_sga_r,r50_sga_r/r50_gal_r,1)
#        slope, intercept, r_value, p_value, std_err = stats.linregress(r50_sga_r,r50_gal_r)
        # Create empty plot with blank marker containing the extra label
        
        plt.axhline(1,color='black',label='1-to-1')
        
        plt.plot([np.min(r50_sga_r),np.max(r50_sga_r)], [np.min(m*r50_sga_r+b),np.max(m*r50_sga_r+b)], color='crimson', alpha=0.5,
                 label=f'Best-fit slope (np.polyfit) {round(m,3)}')
                
        plt.legend(fontsize=16)
                
        plt.xscale('log')
        #plt.yscale('log')
        #plt.xlim(0,100)
        #plt.ylim(0,2)
        
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        
        if savefig==True:
            plt.savefig(homedir+'/Desktop/SGA_r50_comparison.png', bbox_inches='tight', pad_inches=0.2, dpi=200)
        
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
            plt.savefig(homedir+'/Desktop/PSF_r50_comparison.png', bbox_inches='tight', pad_inches=0.2, dpi=200)
        
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
            plt.savefig(homedir+'/Desktop/Re_comparison_Kim.png', bbox_inches='tight', pad_inches=0.2, dpi=200)
        
        plt.show()    

    def wisesize_mass(self, nbins=5, savefig=False):

        self.data = [self.sizerats[(self.clusflag) | (self.rgflag)], 
                     self.sizerats[(self.pgflag) | (self.filflag) | (self.fieldflag)]]
        self.mass_data = [self.z0mgscut['logmass'][(self.clusflag) | (self.rgflag)], 
                          self.z0mgscut['logmass'][(self.pgflag) | (self.filflag) | (self.fieldflag)]]
        
        #lastly...some code does not function (pun intended) if I do not explicltly remove masked rows
        self.data[0] = self.data[0][self.mass_data[0]>0]
        self.data[1] = self.data[1][self.mass_data[1]>0]
        self.mass_data[0] = self.mass_data[0][self.mass_data[0]>0]
        self.mass_data[1] = self.mass_data[1][self.mass_data[1]>0]
        
        mass_data = self.mass_data.copy()   #logmass data
        data = self.data.copy()   #Re data
        logmass = self.z0mgscut['logmass'].copy()
        
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
            
        plt.ylabel(r'R$_{12}$/R$_r$',fontsize=18)
        plt.xlabel(r'log$_{10}$(M$_*$/M$_\odot$)',fontsize=18)
        
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        
        plt.title('Median Re Size Ratios vs. Stellar Mass',fontsize=20)
        plt.ylim(0.5,1)
        plt.xlim(8.5,11)
        plt.legend(fontsize=14)
        
        if savefig==True:
            plt.savefig(homedir+'/Desktop/wisesize_comp_mass_env.png', bbox_inches='tight', pad_inches=0.2, dpi=200)

        plt.show()    
    
if __name__ == '__main__':
    '''
    print("""USAGE:
    cat = catalogs(conv=False,MeanMedian='mean',MADmultiplier=5,cutAGN=False) --> initiate catalog class. 
        MeanMedian will propagate to all plots. MADmultip[lier is a proxy for the extent to which a 
        galaxy must be an 'outlier' in order to be removed from calculations.
        cutAGN is precisely that and will remove flagged AGN from the plots below (both WISE and BPT AGN)
    cat.inclination_check(savefig=False) --> checks for size ratio, SFR, and Mstar dependencies on
        optical axis ratio (found using GALFIT 1-Sersic profile fits on r-band images of the galaxies).
    cat.sfrmstar_magphys(show_HI=False,show_sizerat=True,show_D25=False,savefig=False) --> generates
        SFR v. Mstar plot, color-coded according to either available HI gas mass, D25 (Hyperleda), 
        or GALFIT size ratios. SFR and Mstar values from MAGPHYS.
    cat.sfrmstar_z0mgs(show_HI=False,show_sizerat=True,savefig=False) --> same as above, but with z0mgs.
    cat.r12_vs_rstar(sfr_mstar='sfr', savefig=False) --> plots R_12 vs. R_r effective radii for the full
        subsample as well as according to each environment bin (six total panels). User can dictate
        the colorbar, namely whether data are colored according to SFR or to Mstar (both from MAGPHYS)
    cat.envbins(savefig=False) --> plots number of subsample galaxies in each environment bin
    cat.ratio_MS(showHI=False,savefig=False) --> plots size ratios as a function of MS distance, 
        colored according to normalized HI mass in the case where showHI=True.
    cat.hist_dist_rats(savefig=False) --> size ratio distributions for each of the five defined 
        environment bins (plt.subplots with one column, five rows). Rainbow colors, very aesthestic.
    cat.env_means(trimOutliers=False, errtype='bootstrap', savefig=False) --> plots either mean 
        or median size ratio (w3/r) in each environment bin; trimOutliers will output an additional plot which 
        compares my no PSF parameters to Rose's parameters, allowing the user to visualize which points are omitted 
        in the trimmed env_means plot; errtype either 'bootstrap' or 'std_err'
    cat.mass_hist(z0mgs_comp=True,savefig=False) --> generate mass histogram subplots per environment bin; 
        will compare MAGPHYS stellar masses with z0mgs values if True
    cat.mass_hist_oneplot(savefig=False) --> generate same MAGPHYS mass histograms but as a single plot
    cat.compareSGA(savefig=False) --> compares Rose's GALFIT r-band Re values with SGA's non-parametric r50
    cat.comparePSF(savefig=False) --> plots noPSF vs. PSF Re values for w3, r-band (one subplot per band)
    cat.compareKim(savefig=False) --> compares my noPSF w3-band Re values with Rose's noPSF values
    cat.wisesize_mass(nbins=3,savefig=False) --> creates median/mean size ratios vs. mass bin plot for two
        environment classifications: low-density (PG, filament, field) and high-density (cluster, RG).
    """)
    print('-----------------------------------------------------')
    print()
    '''