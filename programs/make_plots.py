#!/usr/bin/env python

import os
from astropy.io import fits
import numpy as np
from astropy.table import Table
from astropy.io import ascii
from matplotlib import pyplot as plt
from scipy import stats
from numpy import percentile
from scipy.stats import spearmanr


def hist_dist(xlimit_lower=False,xlimit_upper=False,savefig=False):
    ratio_bins = [re_ratio_1,re_ratio_2]
    colors = ['red','blue']
    xlabels = ['Re Ratio (nopsf)','Re Ratio (psf)']
    xlimits = [xlimit_lower,xlimit_upper]
    fig = plt.figure(figsize=(15,10))
    plt.subplots_adjust(hspace=.3,wspace=.2)
    for i in range(1,3):
        ax=fig.add_subplot(2,2,i)
        plt.hist(ratio_bins[i-1],bins=100,color=colors[i-1])
        plt.xlabel(xlabels[i-1],fontsize=18)
        #plt.xlim(xlims[i-1])
        if i == 1:
            plt.text(1.8,32,'mean: '+'%.3f'%(np.mean(re_ratio_1)),fontsize=16)
            plt.text(1.8,28,'median: '+'%.3f'%(np.median(re_ratio_1)),fontsize=16)
            plt.text(1.8,24,'# galaxies: '+str(len(re_ratio_1)),fontsize=16)        
        if i == 2:
            plt.text(1.8,32,'mean: '+'%.3f'%(np.mean(re_ratio_2)),fontsize=16)
            plt.text(1.8,28,'median: '+'%.3f'%(np.median(re_ratio_2)),fontsize=16)
            plt.text(1.8,24,'# galaxies: '+str(len(re_ratio_2)),fontsize=16) 
    if savefig==True:
        plt.savefig('hist_dist.png', dpi=300)
               
            
def sfr_mstar(trimgray=False,savefig=True):
    #size ratio v. offset from main sequence
    mgstab_full = Table.read(homedir+'/Desktop/v1-20210308/vf_north_v1_z0mgs.fits')
    mgstab = mgs

    Mstar = mgstab['logmass'] #x-axis
    SFR = mgstab['logsfr'] #y-axis; frequently found with W4 data
    D25 = 0.1 * 10**(hyp['logd25']) * 60 #arcsec
    D25 = D25[t2_error]
    D25 = D25[sersic_flag_t2]

    #full catalog...
    Mstar_full = mgstab_full['logmass']
    SFR_full = mgstab_full['logsfr']

    #cuts
    Mstar = Mstar[t2_error]
    Mstar = Mstar[sersic_flag_t2]
    SFR = SFR[t2_error]
    SFR = SFR[sersic_flag_t2]

    plt.figure(figsize=(11,6))
    plt.scatter(Mstar_full,SFR_full,color='gray',s=2,alpha=0.3,label='Full SFR Sample')
    #plt.scatter(Mstar,SFR,c=re_ratio_2,cmap='cool',linewidths=1,vmin=0,vmax=2,s=25,alpha=0.8,label='Size Subsample')
    plt.scatter(Mstar,SFR,c=D25,cmap='cool',linewidths=1,vmin=0,vmax=600,s=25,alpha=0.8,label='Size Subsample')

    cb = plt.colorbar()

    cb.set_label(r'Optical D25 (arcsec)',fontsize=20)
    plt.xlabel(r'log(M*/$M_\odot$)',fontsize=22)
    plt.ylabel(r'log(SFR/($M_\odot$/yr))',fontsize=22)

    #best-fit line test...
    #slope: 0.47 +/- 0.01
    #y-intercept: -4.88 +/- 0.10
    y = Mstar_full*0.47 - 4.88

    plt.plot(Mstar_full,y,label=('Main Sequence'),color='black',alpha=0.7)

    plt.legend(loc='upper left',fontsize=18)

    #trim the "gray" perimeter to magnify on the central, colored points
    if trimgray==True:
        plt.ylim(-3,2)
        plt.xlim(7.45,11.5)
        
    if savefig == True:
        plt.savefig('sfr_mstar.png', dpi=300)
        

def ratio_ms(savefig=False):

    #Distance = (| a*x1 + b*y1 + c |) / (sqrt( a*a + b*b))
    d2 = (0.47*Mstar - SFR - 4.88) / (np.sqrt(0.47**2 + 1**2))
    
    plt.figure(figsize=(9,6))
    plt.scatter(d2,re_ratio_2,s=6,alpha=0.5,label='sample')
    plt.yscale('log')
    plt.ylabel(r'Size ratio (12$\mu$m/optical)',fontsize=15)
    plt.xlabel('Distance from MS (log(SFR)/log(M*))',fontsize=15)

    plt.scatter(d2[kauffman_agn_flag2], re_ratio_2[kauffman_agn_flag2], s=85, facecolors='none', 
                edgecolors='g',label='BPT AGN')
    plt.scatter(d2[~wise_agn_flag2], re_ratio_2[~wise_agn_flag2], s=85, facecolors='none', 
                edgecolors='r',label='WISE AGN')
    #plt.plot([np.min(d2),np.max(d2)],[mean_re-std_re,mean_re-std_re],label=r'-1$\sigma$')
    plt.legend(loc='upper left',fontsize=12)
    
    if savefig == True:
        plt.savefig('snr_mstar.png', dpi=300)
        
    #Spearman rank correlation test
    #~kauffman_agn_flag excludes BPT AGN
    #wise_agn_flag exclude WISE AGN

    print('Free Correlation (rho value, no outliers) '+'%0.4f'%spearmanr(d2,re_ratio_2)[0])
    print('Free Correlation (p value, no outliers) '+'%0.4f'%spearmanr(d2,re_ratio_2)[1])
    print(' ')


def envbins(savefig=False):
    

    #first plot num of galaxies in each bin
    fig = plt.figure(figsize=(8,4))
    ax = fig.add_subplot(111)
    ind=np.arange(0,8,1)
    width=0.15

    x = np.array([1,2,3,4,5,6,7,8])-1
    data1 = np.array([len(cluster_re_free),len(fil_rg_re_free),len(rg_all_re_free),len(rg_re_free),
                     len(filament_re_free),len(fil_iso_re_free),len(pg_re_free),len(field_re_free)])

    plt.scatter(x,data1)
    plt.title('Number of galaxies in each bin (tot='+str(len(re_ratio_2))+')',fontsize=15)
    plt.ylabel('# Galaxies',fontsize=15)
    plt.ylim(0,250)

    ax.set_xlim(-width,len(ind)-1.5*width)
    xTickMarks= ['Cluster','Fil and RG','RG','RG (no Fil)','Fil','Fil no RG','PG','Field']
    ax.set_xticks(ind+width)
    xtickNames=ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames,rotation=45,fontsize=13)
    plt.xlim(-0.5,7.5)
    plt.grid(alpha=0.5)
        
    if savefig==True:
        plt.savefig('snr_mstar.png', dpi=300)
        
        
def violin_onion(savefig=False):     
    env_names = ['Cluster','Fil and RG','RG','RG (no Fil)','Fil','Fil no RG','PG','Field']
    #x axis --> env names; y axis --> Re ratios

    #set up the canvas
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15, 4))

    #data corresponding, I assume, to each env bin. I will populate the list as such
    data =[cluster_re_free, fil_rg_re_free, rg_all_re_free, rg_re_free, filament_re_free, fil_iso_re_free,
           pg_re_free, field_re_free]

    # mod = mode?
    mod = []

    # median
    med = []

    # inter-quartile ranges; IQR = stats.iqr(data, interpolation = 'midpoint')
    iqr = []

    # upper and lower adjacent values
    # UAV = Q3 + 1.5*IQR 
    # LAV = Q1 - 1.5*IQR
    avs = []

    #taken from https://matplotlib.org/2.0.2/examples/statistics/customized_violin_demo.html
    def adjacent_values(vals):  
        q1 = np.percentile(vals,25)
        q3 = np.percentile(vals,75)
        upper_adjacent_value = q3 + ((q3 - q1) * 1.5)
        upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])
        lower_adjacent_value = q1 - ((q3 - q1) * 1.5)
        lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
        return lower_adjacent_value, upper_adjacent_value   

    for arr in data:
        sarr = sorted(arr)
        med.append(percentile(sarr, 50))   #used 50 instead of 0.5
        mod.append(stats.mode(np.asarray(sarr)[np.isfinite(sarr)].astype(int)))
        iqr.append([percentile(sarr, 32), percentile(sarr, 68)])   #used 32, 68 instead of 0.32, 0.68
        avs.append(adjacent_values(sarr))

    plt.rcParams['axes.axisbelow'] = True
    plt.rc('grid', c='0.5', ls='-', lw=0.2)
    plt.grid(True)

    ticks=np.arange(0,8,1)
    parts = ax.violinplot(data, positions=ticks,widths =0.5, showmeans=False, showmedians=False,
                      showextrema=True)

    med_col = ['tomato','coral','orange','yellow','chartreuse','cyan','fuchsia','orchid']
    avs_col = ['tomato','coral','orange','yellow','chartreuse','cyan','fuchsia','orchid']

    for i in range(len(med)):
        ax.plot([i, i], avs[i], '-', lw=1, zorder = 7,color=avs_col[i])
        ax.plot([i, i], iqr[i], '-', lw=5, zorder = 7,color='black')
        ax.plot(i, med[i], 'o', mec='none', ms=8, zorder = 7,color=med_col[i],label=str(env_names[i])+' Median')

    med_ld = med
    err_ld = iqr

    plt.xticks(ticks,env_names, rotation=10, fontsize=20)

    plt.ylabel(r"R$_e$ Ratios (12$\mu$m/optical)",fontsize=20)      
    
    if savefig==True:
        plt.savefig('snr_mstar.png', dpi=300)
    
    
def env_means(savefig=False):    
    data2 = [cluster_re_free, rg_re_free, pg_re_free, filament_re_free, field_re_free]
    mean2 = []
    err2 = []
    index = [1,2,3,4,5]
    env_names2 = ['Cluster','RG','PG','Fil','Field']

    for j,i in enumerate(data2):
        mean2.append(np.mean(i))
        err2.append(np.std(i)/np.sqrt(len(i)))

    plt.figure(figsize=(10,6))
    plt.scatter(index,mean2)
    plt.errorbar(index,mean2,yerr=err2,fmt='o')
    
    xmin,xmax = plt.xlim()
    xfield = np.linspace(xmin,xmax,50)
    ymax = np.ones(50)*(mean2[-1] + err2[-1])
    ymin = np.ones(50)*(mean2[-1] - err2[-1])
    plt.fill_between(xfield,ymax,ymin,color='0.5',alpha=.2)

    plt.xticks(index,env_names2, rotation=10, fontsize=20)
    plt.grid()
    plt.ylabel(r'Mean (12$\mu$m/optical)',fontsize=20)
    print(' ')    
    
    if savefig==True:
        plt.savefig('snr_mstar.png', dpi=300)
    
        
if __name__ == '__main__':
    
    homedir = os.getenv("HOME")
    
    #note: sgacut is already selected according to SNR>10 (photutils elliptical apertures) and t-type>=1 (Hyperleda classifications); see notebook for guidance regarding the latter.
    #sga matched to vf catalog
    sgacut = Table.read(homedir+'/sgacut_coadd.fits')
    #sga r-band r50 parameters
    sga_params = Table.read(homedir+'/sga_params')
    #hyperleda catalog for t-type
    hyp = Table.read(homedir+'/Desktop/v1-20210308/vf_north_v1_hyperleda.fits')
    #for AGN cuts
    wise = Table.read(homedir+'/Desktop/v1-20210308/vf_north_v1_unwise.fits')
    wise.add_column(hyp['VFID'],name='VFID')
    #for BPT AGN cuts
    bpt = Table.read(homedir+'/Desktop/v1-20210308/vf_north_v1_nsa_v0.fits')
    #environment classifications
    env = Table.read(homedir+'/Desktop/v1-20210308/vf_north_v1_main_envsummary_Kourkchi.fits')
    #z0mgs from Leroy+19
    mgstab = Table.read(homedir+'/Desktop/v1-20210308/vf_north_v1_z0mgs.fits')
    #general Virgo v1 table
    vf = Table.read(homedir+'/Desktop/v1-20210308/vf_north_v1_main.fits')
    
    #trimming the files to all possess compatible lengths and galaxies
    
    #sga_params --> sga_params table only containing sgacut entries
    #sgacut --> resized sgacut to be compatible with sga_p
    flag_p = []
    for i in sga_params['name']:
        if i not in sgacut['VFID']:
            flag_p.append(False)
        else:
            flag_p.append(True)
    sga_params = sga_params[flag_p]
    flag_p2 = []
    for i in sgacut['VFID']:
        if i not in sga_params['name']:
            flag_p2.append(False)
        else:
            flag_p2.append(True)
    sgacut = sgacut[flag_p2]
    
    #use hyp['VFID'] as the reference column, since all catalogs contain these entries
    flag_all = []
    for i in hyp['VFID']:
        if i not in sga_params['name']:
            flag_all.append(False)
        else:
            flag_all.append(True)
    
    hyp = hyp[flag_all]
    wise = wise[flag_all]
    bpt = bpt[flag_all]
    env = env[flag_all]
    mgs = mgstab[flag_all]
    vf = vf[flag_all]
    
    #**large** size ratio outliers, update as applicable
    outlier_vfid = ['VFID0432','VFID1533','VFID1932','VFID2627','VFID2762','VFID2977','VFID3289',
               'VFID4122','VFID4272','VFID4877','VFID5918','VFID6357','VFID6454','VFID0721','VFID4454',
               'VFID5360']
    
    #trim all catalogs accordingly
    flag_outliers = []
    for i in hyp['VFID']:
        if i in outlier_vfid:
            flag_outliers.append(False)
        else:
            flag_outliers.append(True)

    ### FINAL CUTS ###
    hyp = hyp[flag_outliers]
    wise = wise[flag_outliers]
    bpt = bpt[flag_outliers]
    env = env[flag_outliers]
    mgs = mgs[flag_outliers]
    vf = vf[flag_outliers]
    
    flag_sga = []
    flag_sga2 = []
    
    for i in sga_params['name']:
        if i not in hyp['VFID']:
            flag_sga.append(False)
        else:
            flag_sga.append(True)
    for i in sgacut['name']:
        if i not in hyp['VFID']:
            flag_sga2.append(False)
        else:
            flag_sga2.append(True)
    ### FINAL CUTS ###
    sga_params = sga_params[flag_sga]
    sgacut = sgacut[flag_sga2]
    
    #cut GALFIT output parameter tables
    t1 = Table.read(homedir+'/nopsf_params.txt',format='ascii')
    t2 = Table.read(homedir+'/psf_params.txt',format='ascii')
    
    #generate flag that only permit galaxies appearing in SGA to be filtered through
    flag_t1 = []
    flag_t2 = []
    for i in t1['galname']:
        if str(i) not in sgacut['prefix']:
            flag_t1.append(False)
        else:
            flag_t1.append(True)
        
    for i in t2['galname']:
        if str(i) not in sgacut['prefix']:
            flag_t2.append(False)
        else:
            flag_t2.append(True)
    
    t1 = t1[flag_t1]
    t2 = t2[flag_t2]

    t1_error = t1['error']
    t2_error = t2['error']
    
    #isolate and convert Re to arcseconds, for compatibility purposes with the optical Re
    #for W3 camera, one px = 2.75 arcsec

    re_free_nopsf = t1['re']*2.75
    re_free_psf = t2['re']*2.75
    
    #convert sga params to optical half-light radius for galaxy subsample
    #(values are for r-band)
    #note: parameters are dimensionless, so multiplying by 10 arcsec yields an r50 in arcsec

    r0 = 10 #arcsec
    alpha_1 = sga_params['R_COG_PARAMS_ALPHA1']
    alpha_2 = sga_params['R_COG_PARAMS_ALPHA2']
    m0 = sga_params['R_COG_PARAMS_M0']

    comp_one = np.log10(0.5)/(0.4*m0)
    comp_two = np.exp(-comp_one) - 1
    comp_three = (1/alpha_1) * comp_two
    comp_four = comp_three**(-1/alpha_2)
    #label as tentative, since we must rearrange the order 
    r50_tentative = r0 * comp_four
    
    r50 = []
    for i in range(0,len(sgacut)):
        galname = sgacut['VFID'][i]
        #find index at which sga_params (corresponding to r50_tentative index) fits the sgacut galname 
        index = np.where(sga_params['name'] == galname)[0][0]
        #use that index to extract r50_tentative element and place in r50
        r50.append(r50_tentative[index])

        #calculating ratios --> (12-micron Re) / (r-band Re)

    re_ratio_free_nopsf = re_free_nopsf/r50
    re_ratio_free_psf = re_free_psf/r50

    #and filter errors --> galaxies with numerical errors are excluded from analysis

    #first convert 0,1 to True, False

    t1_error_list = []
    t2_error_list = []

    for i in range(0,len(sgacut)):
        if t1_error[i] == 0:
            t1_error_list.append(True)
        if (t1_error[i] == 1) | (t1_error[i] == -999):
            t1_error_list.append(False)
        if t2_error[i] == 0:
            t2_error_list.append(True)
        if (t2_error[i] == 1) | (t2_error[i] == -999):
            t2_error_list.append(False)

    #redefining errors such that int values become bools.

    t1_error = np.asarray(t1_error_list)
    t2_error = np.asarray(t2_error_list)
    
    #mark AGN in sample, according to either WISE colors or BPT diagram cuts
    #  can use WISE colors to isolate galaxies with AGN from sample
    #      W1-W2 > 0.75, W2-W3 < 4.0
    #Asmus+20, Jarret+17 in agreement with these cut selections
    wise_agn = wise['VFID'][np.where((wise['w1_mag'] - wise['w2_mag']>0.75) & (wise['w2_mag']-wise['w3_mag']<4.0))[0]]
    #generate WISE flag (removes WISE AGN from sample)
    wise_agn_flag = []
    for i in wise['VFID']:
        wise_agn_flag.append(i not in wise_agn)
    
    #do likewise for BPT spectral data!
    #code courtesy of Rose Finn
    agn_kauffman = (np.log10(bpt['O3FLUX']/bpt['HBFLUX']) > (.61/(np.log10(bpt['N2FLUX']/bpt['HAFLUX']-.05))+1.3)) | (np.log10(bpt['N2FLUX']/bpt['HAFLUX']) > 0.)
    agn_kauffman_flag = np.asarray(agn_kauffman)

    
    #then treat the previously defined GALFIT errors as flags...no error behinds left.
    #(i.e., remove galaxies with error flags)

    #strictly 12-micron Re
    re_1 = re_free_nopsf[t1_error]
    re_2 = re_free_psf[t2_error]

    #strictly r-band Re
    r50_1 = np.asarray(r50)[t1_error]
    r50_2 = np.asarray(r50)[t2_error]

    #Re ratios (recall: 12-micron/r-band)
    re_ratio_1 = re_ratio_free_nopsf[t1_error]
    re_ratio_2 = re_ratio_free_psf[t2_error]

    #also trim WISE AGN flags...
    wise_agn_flag1 = np.asarray(wise_agn_flag)[t1_error]
    wise_agn_flag2 = np.asarray(wise_agn_flag)[t2_error]

    #annnnd, BPT AGN flags...
    kauffman_agn_flag1 = np.asarray(agn_kauffman_flag)[t1_error]
    kauffman_agn_flag2 = np.asarray(agn_kauffman_flag)[t2_error]
    
    
    #AND FILTER n==5 or n>6 galaxies:
    nser1 = t1['nsersic'][t1_error]
    nser2 = t2['nsersic'][t2_error]

    sersic_flag1_t1 = (nser1!=5) 
    sersic_flag2_t1 = (nser1 < 6)
    sersic_flag1_t2 = (nser2!=5)
    sersic_flag2_t2 = (nser2 < 6)

    sersic_flag_t1 = (nser1!=5) & (nser1<6)
    sersic_flag_t2 = (nser2!=5) & (nser2<6)

    #strictly 12-micron Re
    re_1 = re_1[sersic_flag_t1]
    re_2 = re_2[sersic_flag_t2]

    #strictly r-band Re
    r50_1 = r50_1[sersic_flag_t1]
    r50_2 = r50_2[sersic_flag_t2]

    #Re ratios (recall: 12-micron/r-band)
    re_ratio_1 = re_ratio_1[sersic_flag_t1]
    re_ratio_2 = re_ratio_2[sersic_flag_t2]

    #also trim WISE AGN flags...
    wise_agn_flag1 = wise_agn_flag1[sersic_flag_t1]
    wise_agn_flag2 = wise_agn_flag2[sersic_flag_t2]

    #annnnd, BPT AGN flags...
    kauffman_agn_flag1 = kauffman_agn_flag1[sersic_flag_t1]
    kauffman_agn_flag2 = kauffman_agn_flag2[sersic_flag_t2]
    
    envcut = env
    #isolate columns; convert 0,1 to False,True

    cluster = envcut['flag_clus']
    cluster = np.asarray(list(map(bool,cluster)))

    filament = envcut['flag_fil']
    filament = np.asarray(list(map(bool,filament)))

    field = envcut['flag_isol']
    field = np.asarray(list(map(bool,field)))

    #Environment cuts
    
    group = envcut['flag_gro']
    #poor group
    pg_flag = []
    #rich group
    rg_flag = []

    #define rich, poor groups 

    for n in range(0,len(envcut)):
        if int(group[n]) == 1:
            rg_flag.append(False)
            pg_flag.append(True)
        if int(group[n]) == 2:
            rg_flag.append(True)
            pg_flag.append(False)
        if int(group[n]) == 0:
            rg_flag.append(False)
            pg_flag.append(False)

    rg_flag = np.asarray(rg_flag)
    pg_flag = np.asarray(pg_flag)
        
    #CL
    cluster = cluster

    #Filament and rich group
    fil_rg = (rg_flag & filament)

    #rich group (all, meaning isolated rich groups or those in filaments)
    rg_all = (rg_flag)

    #rich group only (include if the row only has a "2")
    rg = (rg_flag) & (~fil_rg)

    #filament (all galaxies with filament marker)
    filament = filament

    #filament galaxies that do not belong to a rich group
    fil_only = (filament) & (~rg_flag)

    #poor group
    pg_flag = pg_flag

    #pure field
    field = field
        
    #trimming bool parameters according to galaxies with GALFIT numerical errors, for free case

    cluster_free = cluster[t2_error]
    cluster_free = cluster_free[sersic_flag_t2]
    fil_rg_free = fil_rg[t2_error]
    fil_rg_free = fil_rg_free[sersic_flag_t2]
    rg_all_free = rg_all[t2_error]
    rg_all_free = rg_all_free[sersic_flag_t2]
    rg_free = rg[t2_error]
    rg_free = rg_free[sersic_flag_t2]
    filament_free = filament[t2_error]
    filament_free = filament_free[sersic_flag_t2]
    filament_iso_free = fil_only[t2_error]
    filament_iso_free = filament_iso_free[sersic_flag_t2]
    pg_free = pg_flag[t2_error]
    pg_free = pg_free[sersic_flag_t2]
    field_free = field[t2_error]
    field_free = field_free[sersic_flag_t2]
        
    #isolate re_ratios for each category

    cluster_re_free = re_ratio_2[cluster_free]
    fil_rg_re_free = re_ratio_2[fil_rg_free]
    rg_all_re_free = re_ratio_2[rg_all_free]
    rg_re_free = re_ratio_2[rg_free]
    filament_re_free = re_ratio_2[filament_free]
    fil_iso_re_free = re_ratio_2[filament_iso_free]
    pg_re_free = re_ratio_2[pg_free]
    field_re_free = re_ratio_2[field_free]
    
    
    
    
    
    
   