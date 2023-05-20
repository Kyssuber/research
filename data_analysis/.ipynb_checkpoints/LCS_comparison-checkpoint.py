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
from scipy.stats import median_abs_deviation as MAD
from scipy.stats import kstest
from astropy.stats import bootstrap
from scipy.stats import binned_statistic

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
    
    def __init__(self,conv=False,MeanMedian='mean',cutAGN=False):
        path_to_dir = homedir+'/Desktop/v2-20220820/'
        self.v2_env = Table.read(path_to_dir+'vf_v2_environment.fits')
        self.v2_main = Table.read(homedir+'/v2_snrcoadd.fits')
        self.magphys = Table.read(path_to_dir+'vf_v2_magphys_07-Jul-2022.fits')
        self.z0mgs = Table.read(path_to_dir+'vf_v2_z0mgs.fits')
        
        self.MeanMedian = MeanMedian  #whether I plot median or mean size ratios for the self.env_means() figure
        
        self.conv = conv
        print('Convolution: ',str(self.conv))
        self.cutAGN=cutAGN
        
        if self.conv==False:
            self.rdat = Table.read(homedir+'/output_params_r_nopsf.fits')
            self.w3dat = Table.read(homedir+'/output_params_W3_nopsf.fits')
        if self.conv==True:
            self.rdat = Table.read(homedir+'/output_params_r_psf.fits')
            self.w3dat = Table.read(homedir+'/output_params_W3_psf.fits')
    
        self.cut_cats()
    
    def cut_cats(self):
        subsample_flag = self.v2_main['sgacut_flag']
        
        self.v2_env = self.v2_env[subsample_flag]
        self.v2_maincut = self.v2_main[subsample_flag]
        self.magphyscut = self.magphys[subsample_flag]
        self.z0mgscut = self.z0mgs[subsample_flag]
        
        self.re_rband = self.rdat['re']
        self.re_w3band = self.w3dat['re']
        self.nser = self.rdat['nsersic']

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
            wise_agn = (wise_mag_cut['w1_mag'] - wise_mag_cut['w2_mag']>0.75) & (wise_mag_cut['w2_mag']-wise_mag_cut['w3_mag']<4.0)
            agn_kauffman = (np.log10(bpt_lines_cut['O3FLUX']/bpt_lines_cut['HBFLUX']) > (.61/(np.log10(bpt_lines_cut['N2FLUX']/bpt_lines_cut['HAFLUX']-.05))+1.3)) | (np.log10(bpt_lines_cut['N2FLUX']/bpt_lines_cut['HAFLUX']) > 0.)
        
            self.wise_agn_flag = np.asarray(wise_agn)
            self.agn_kauffman_flag = np.asarray(agn_kauffman)
            AGN_flags = (self.wise_agn_flag)|(self.agn_kauffman_flag)

            print(f'# WISE AGN in VF subsample: {len(self.v2_maincut[self.wise_agn_flag])}')
            print(f'# BPT AGN in VF subsample: {len(self.v2_maincut[self.agn_kauffman_flag])}')
            
            self.cut_flags = (~AGN_flags) & (~fail_flag) & (~err_flag)
        
        self.re_rband_cut = self.re_rband[self.cut_flags]
        self.re_w3band_cut = self.re_w3band[self.cut_flags]
        self.nser = self.nser[self.cut_flags]
        
        #apply final cut to envcut and maincut catalogs
        self.v2_envcut = self.v2_env[self.cut_flags]
        self.v2_maincut = self.v2_maincut[self.cut_flags]
        self.magphyscut = self.magphyscut[self.cut_flags]
        self.z0mgscut = self.z0mgscut[self.cut_flags]

        n_tot = len(self.re_rband)     #should be the same for both r-band and w3

        print(f'No GALFIT data for {n_fails_w3} w3 galaxies and {n_fails_r} r galaxies.')
        print(f'Total number of galaxies with GALFIT error flags: {int(np.sum(np.ones(len(err_flag))*err_flag))}')
        print(f'Total number of galaxies: {n_tot}')
        print()
        self.sizerats = (self.re_w3band_cut*2.75)/(self.re_rband_cut*0.262)
        
        #introduce some important quantities for infall/core distinction
        #define core as being within 1 R200 of the cluster center, infalling as between 1 and 5 R200 from the center

        self.RA = self.v2_maincut['RA_1']
        self.DEC = self.v2_maincut['DEC_1']
        self.RA_center = 187.70
        self.DEC_center = 12.34
        virial_radius = 5.8224 #(degrees)
        five_vir = 5*virial_radius
        
        dist = np.sqrt((self.RA-self.RA_center)**2 + (self.DEC-self.DEC_center)**2)
        virial_1flag = (dist<=virial_radius)
        virial_5flag = (dist>virial_radius)&(dist<=five_vir)
        
        self.data_core = self.sizerats[self.v2_envcut['cluster_member'] & (self.nser<2)]
        self.data_fall = self.sizerats[virial_5flag & (self.nser<2) & (~self.v2_envcut['cluster_member'])]
        
        self.data = [self.data_core, self.data_fall]
        self.mass_data = [self.z0mgscut['logmass'][self.v2_envcut['cluster_member'] & (self.nser<2)], 
                          self.z0mgscut['logmass'][virial_5flag & (self.nser<2) & (~self.v2_envcut['cluster_member'])]]
        
        #lastly...some code does not function (pun intended) if I do not explicltly remove masked rows
        self.data[0] = self.data[0][self.mass_data[0]>0]
        self.data[1] = self.data[1][self.mass_data[1]>0]
        self.mass_data[0] = self.mass_data[0][self.mass_data[0]>0]
        self.mass_data[1] = self.mass_data[1][self.mass_data[1]>0]
        
    def LCS_hist(self, ReDisk='Re', savefig=False):
        
        lcs_tab = Table.read(homedir+'/Desktop/misc_cats/LCS_final_sample.fits')
        disk_sizes = lcs_tab['sizeratio_disk']
        re_sizes = lcs_tab['sizeratio_re']
        
        #if disk sizes, restrict B/T ratio; if Re, restrict sersic index
        if ReDisk=='Re':
            infall_flag = (lcs_tab['infall']) & (lcs_tab['sersic_n']<2)
            core_flag = (lcs_tab['core']) & (lcs_tab['sersic_n']<2)
        if ReDisk=='Disk':
            infall_flag = (lcs_tab['infall']) & (lcs_tab['BT']<0.3)
            core_flag = (lcs_tab['core']) & (lcs_tab['BT']<0.3)
                    
        data = [disk_sizes[core_flag],disk_sizes[infall_flag]]
        titles = ['LCS Disk Size Ratio Distribution (B/T<0.3)', '']
        if ReDisk=='Re':
            data = [re_sizes[core_flag],re_sizes[infall_flag]]
            titles = ['LCS Re Size Ratio Distribution (nser<2)', '']
        
        labels = ['Core Galaxies', 'Infalling Galaxies']
        
        colors = ['crimson', 'blue']
        xlabels = ['',r'R$_{24}$/R$_r$']
        
        mybins=np.linspace(0,3,30)
        
        fig = plt.figure(figsize=(8,8))
        plt.subplots_adjust(hspace=.2)
        
        for panel in range(2):
            ax = fig.add_subplot(2,1,panel+1)
            plt.hist(data[panel],bins=mybins,color=colors[panel], alpha=0.7, label=labels[panel],density=False,cumulative=False)
            #plt.hist([], color='white', label=r'0 < R$_{24}$/R$_{r}$ < 2')
            plt.title(titles[panel],fontsize=16)
            ax.set_xlabel(xlabels[panel],fontsize=20)
            ax.set_ylabel(r'N$_{gal}$',fontsize=20)
            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)
            plt.xlim(-0.1,2)
            ax.legend(fontsize=15)
        
        if savefig==True:
            plt.savefig(homedir+'/Desktop/lcs_comp_hists.png', bbox_inches='tight', pad_inches=0.2, dpi=200)

        plt.show()
        
        #ks-test statistics
        uncertainty_core = np.std(data[0][data[0]<=3])/np.sqrt(len(data[0][data[0]<=3]))
        uncertainty_infall = np.std(data[0][data[0]<=3])/np.sqrt(len(data[0][data[0]<=3]))
        
        print('K-S p-value (> 0.003 (3sigma), "same distribution"):')
        print('%.5f'%(kstest(data[0],data[1])[1]))
        print()
        print('Mean core size ratio: %.3f'%(np.mean(data[0])))
        print('Median core size ratio: %.3f'%(np.median(data[0])))
        print('Mean infalling size ratio: %.3f'%(np.mean(data[1])))
        print('Median infalling size ratio: %.3f'%(np.median(data[1])))
        print()
        print('Core mean/median uncertainty: %.3f'%uncertainty_core)
        print('Infalling mean/median uncertainty: %.3f'%uncertainty_infall)
    
    def LCS_mass(self, ReDisk='Re', savefig=False):
        
        lcs_tab = Table.read(homedir+'/Desktop/misc_cats/LCS_final_sample.fits')
        disk_sizes = lcs_tab['sizeratio_disk']
        re_sizes = lcs_tab['sizeratio_re']
        logmass = lcs_tab['logMstar']

        #if disk sizes, restrict B/T ratio; if Re, restrict sersic index
        if ReDisk=='Re':
            infall_flag = (lcs_tab['infall']) & (lcs_tab['sersic_n']<2)
            core_flag = (lcs_tab['core']) & (lcs_tab['sersic_n']<2)
            data = [re_sizes[core_flag],re_sizes[infall_flag]]
        if ReDisk=='Disk':
            infall_flag = (lcs_tab['infall']) & (lcs_tab['BT']<0.3)
            core_flag = (lcs_tab['core']) & (lcs_tab['BT']<0.3)
            data = [disk_sizes[core_flag],disk_sizes[infall_flag]]
        mass_data = [lcs_tab['logMstar'][core_flag],lcs_tab['logMstar'][infall_flag]]
        
        nbins=6
        
        #define mass bins --> min<log(Mstar)<max, nbins
        mass_bounds = np.linspace(9.1,10.6,nbins+1)
        bin_bounds = []
        for bound in range(nbins):
            bin_bounds.append([mass_bounds[bound],mass_bounds[bound+1]])
        
        #now for the tedious task of creating flags which will separate the ratios into these mass bins :}
        re_mass_clus=[]
        re_mass_ext=[]
        mass_coord=[]
        err_clus=[]
        err_ext=[]
        
        #note --> bootstrap errors correspond to 68% confidence interval
        for bound in bin_bounds: 
            bound_flag_clus = (mass_data[0]>bound[0])&(mass_data[0]<=bound[1])
            bound_flag_fall = (mass_data[1]>bound[0])&(mass_data[1]<=bound[1])
            bound_flag = (logmass>bound[0])&(logmass<=bound[1])
            if self.MeanMedian=='mean':
                avg_re_clus = np.mean(data[0][bound_flag_clus])
                avg_re_ext = np.mean(data[1][bound_flag_fall])
                lower_clus, upper_clus = get_bootstrap_confint(data[0][bound_flag_clus],bootfunc=np.mean,nboot=1000)
                lower_ext, upper_ext = get_bootstrap_confint(data[1][bound_flag_fall],bootfunc=np.mean,nboot=1000)
            if self.MeanMedian=='median':
                avg_re_clus = np.median(data[0][bound_flag_clus])
                avg_re_ext = np.median(data[1][bound_flag_fall])
                lower_clus, upper_clus = get_bootstrap_confint(data[0][bound_flag_clus],bootfunc=np.median,nboot=1000)
                lower_ext, upper_ext = get_bootstrap_confint(data[1][bound_flag_fall],bootfunc=np.median,nboot=1000) 
            
            mass_coord.append(np.mean(logmass[bound_flag]))
            re_mass_clus.append(avg_re_clus)
            re_mass_ext.append(avg_re_ext)
            
            #err_clus.append([avg_re_clus-lower_clus, upper_clus-avg_re_clus])
            #err_ext.append([avg_re_ext-lower_ext, upper_ext-avg_re_ext])
            err_clus.append([lower_clus, upper_clus])
            err_ext.append([lower_ext, upper_ext])
        
        plt.figure(figsize=(8,6))
        
        plt.scatter(mass_data[0],data[0],color='crimson',s=15,alpha=0.1,label='Core',zorder=1)
        plt.scatter(mass_data[1],data[1],color='blue',s=15,alpha=0.1,label='Infalling',zorder=1)
        plt.scatter(mass_coord,re_mass_clus,color='crimson',s=250,edgecolors='black',label='<Core>',zorder=3)
        plt.scatter(mass_coord,re_mass_ext,color='blue',s=250,edgecolors='black',label='<Infalling>',zorder=3)
        
        for n in range(nbins):
            plt.plot([mass_coord[n],mass_coord[n]], [err_clus[n][0],err_clus[n][1]],color='crimson',zorder=2)
            #create lower, upper caps on errorbars
            plt.plot([mass_coord[n]-0.04,mass_coord[n]+0.04],[err_clus[n][0],err_clus[n][0]],color='crimson',zorder=2)  
            plt.plot([mass_coord[n]-0.04,mass_coord[n]+0.04],[err_clus[n][1],err_clus[n][1]],color='crimson',zorder=2)
            
            plt.plot([mass_coord[n],mass_coord[n]], [err_ext[n][0],err_ext[n][1]],color='blue',zorder=2)
            #create lower, upper caps on errorbars
            plt.plot([mass_coord[n]-0.04,mass_coord[n]+0.04],[err_ext[n][0],err_ext[n][0]],color='blue',zorder=2)
            plt.plot([mass_coord[n]-0.04,mass_coord[n]+0.04],[err_ext[n][1],err_ext[n][1]],color='blue',zorder=2)
            
        plt.ylabel(r'R$_{24}$/R$_r$',fontsize=18)
        plt.xlabel(r'log$_{10}$(M$_*$/M$_\odot$)',fontsize=18)
        
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        
        plt.title('LCS Median Re Size Ratios (nser<2)',fontsize=20)
        if ReDisk=='Disk':
            plt.title('LCS Median Disk Size Ratios (B/T<0.3)',fontsize=20)
        plt.ylim(0,1.6)
        plt.xlim(9,10.6)
        plt.legend(fontsize=14)
        
        if savefig==True:
            plt.savefig(homedir+'/Desktop/LCS_comp_mass.png', bbox_inches='tight', pad_inches=0.2, dpi=200)

        plt.show()    

    def wisesize_hist(self, savefig=False):
        
        labels = ['Cluster Galaxies', 'Infalling Galaxies']
        titles = ['WISESize Size Ratio Distribution (r-band nser<2)', '']
        colors = ['crimson', 'blue']
        xlabels = ['',r'R$_{12}$/R$_r$']

        mybins=np.linspace(0,3,40)

        fig = plt.figure(figsize=(8,8))
        plt.subplots_adjust(hspace=.2)

        for panel in range(2):
            ax = fig.add_subplot(2,1,panel+1)
            plt.hist(self.data[panel],bins=mybins,color=colors[panel], alpha=0.7,label=labels[panel],
                     density=False,cumulative=False)
            plt.hist([], color='white', label=r'0 < R$_{12}$/R$_{r}$ < 2')
            plt.title(titles[panel],fontsize=16)
            ax.set_xlabel(xlabels[panel],fontsize=20)
            ax.set_ylabel(r'N$_{gal}$',fontsize=20)
            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)
            #ax.set_yticks([0,])
            plt.xlim(-0.1,2)
            ax.legend(fontsize=15)

        if savefig==True:
            plt.savefig(homedir+'/Desktop/wisesize_comp_hists.png', bbox_inches='tight', pad_inches=0.2, dpi=200)

        plt.show()
        
        data_core = self.data[0]
        data_fall = self.data[1]
        
        #ks-test statistics
        print('Note: the following results are for galaxies with size ratios <= 3')
        print()
        print('K-S p-value (> 0.003, "same distribution"):')
        print('%.5f'%(kstest(data_core[data_core<=3],data_fall[data_fall<=3])[1]))
        print()
        print('Mean core size ratio: %.3f'%(np.mean(data_core[data_core<=3])))
        print('Median core size ratio: %.3f'%(np.median(data_core[data_core<=3])))
        print('Mean infalling size ratio: %.3f'%(np.mean(data_fall[data_fall<=3])))
        print('Median infalling size ratio: %.3f'%(np.median(data_fall[data_fall<=3])))
        print()
        print('Core mean/median uncertainty: %.3f'%(np.std(data_core[data_core<=3])/np.sqrt(len(data_core[data_core<=3]))))
        print('Infalling mean/median uncertainty: %.3f'%(np.std(data_fall[data_fall<=3])/np.sqrt(len(data_fall[data_fall<=3])))) 
    
    def wisesize_mass(self, nbins=3, savefig=False):
        
        mass_data = self.mass_data.copy()   #logmass data
        data = self.data.copy()   #Re data
        logmass = self.z0mgscut['logmass'].copy()
        
        nbins=nbins
        
        bin_med_clus, bin_edges_clus, binnumber_clus = binned_statistic(mass_data[0],data[0],statistic='median', bins=nbins)
        bin_centers_clus = .5*(bin_edges_clus[:-1]+bin_edges_clus[1:])
        
        bin_med_fall, bin_edges_fall, binnumber_fall = binned_statistic(mass_data[1],data[1],statistic='median', bins=nbins)
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
        
        plt.scatter(mass_data[0],data[0],color='crimson',s=15,alpha=0.2,label='Cluster',zorder=1)
        plt.scatter(mass_data[1],data[1],color='blue',s=15,alpha=0.2,label='Infall',zorder=1)
        plt.scatter(bin_centers_clus, bin_med_clus, color='crimson', s=250, edgecolors='black', label='<Cluster>',zorder=3)
        plt.scatter(bin_centers_fall, bin_med_fall, color='blue', s=250, edgecolors='black', label='<Infall>',zorder=3)

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
        
        plt.title('Median Re Size Ratios vs. Stellar Mass (VFS)',fontsize=20)
        #plt.ylim(0.5,0.92)
        plt.xlim(8,11.5)
        plt.legend(fontsize=14)
        
        if savefig==True:
            plt.savefig(homedir+'/Desktop/wisesize_comp_mass.png', bbox_inches='tight', pad_inches=0.2, dpi=200)

        plt.show()    
    


        
if __name__ == '__main__':
    print("""USAGE:
    lcs = catalogs(conv=False, MeanMedian='mean',cutAGN=False) --> initiate catalog class. 
        MeanMedian will propagate to all plots; can either be 'mean' or 'median.' cutAGN
        will do precisely that, flagging both BPT AGN and AGN based on WISE magnitude colors.
    lcs.wisesize_hist(savefig=False) --> generates vertically-oriented histogram subplots of 
        R12/Rr distribution, separated into cluster vs. all else (external).
    lcs.LCS_hist(ReDisk='Re',savefig=False) --> generates same type of histogram plot as above but with 
        the actual LCS values. 'Re' or 'Disk' to instruct which column to use for size ratios.
    lcs.wisesize_mass(savefig=False) --> generates scatterplot of size ratio vs. mass bin using 
        LCS data, the format being similar to a skeleton version of Figure 13 from Finn+18. 
        Masses are from John Moustakas' SED fitting. 
    lcs.wisesize_mass(nbins=3, savefig=False) --> generates scatterplot of size ratio vs. mass bin 
        using VFS data, the format being similar to a skeleton version of Figure 13 from Finn+18. 
        I use z0mgs masses here (odd results with MAGPHYS for min/max, also many galaxies with nan 
        or masked values). 
    """)
    print('-----------------------------------------------------')
    print()