import numpy as np
from numpy import argsort
from matplotlib import pyplot as plt
from astropy.table import Table
import os
homedir=os.getenv("HOME")

class altvals():
    def __init__(self):
        #first load the tables.
        path_to_dir = homedir+'/Desktop/v2-20220820'
        self.phot_all = Table.read(f'{path_to_dir}/virgofilaments-v3b-legacyphot.fits',hdu=2)
        self.phot_names_all = Table.read(f'{path_to_dir}/virgofilaments-v3b-legacyphot.fits',hdu=1)['VFID']
        self.phot_maskbit_all = Table.read(f'{path_to_dir}/virgofilaments-v3b-legacyphot.fits',hdu=3)['MASKBITS']
        self.vf_all = Table.read(f'{path_to_dir}/vf_v2_environment.fits')
        self.magphys_all = Table.read(f'{path_to_dir}/vf_v2_magphys_legacyExt_mergedNS.fits')
        
        self.trim_tables()
        
#now that we are all initiated...

    def trim_tables(self):
        #John's photometry catalog does not include the most massive galaxies due to runtime restrictions and other problematic characteristics (difficulties with establishing where the galaxy ends, greater likelihood of intersections with foreground stars, usw). the nameflag ensures vf and magphys begin at the same size as phot_*
        nameflag = np.zeros(len(self.vf_all),dtype=bool)
        for i in range(len(self.vf_all)):
            if self.vf_all[i]['VFID'] in self.phot_names_all:
                nameflag[i]=True
            else:
                nameflag[i]=False

        #trim the two tables accordingly
        self.vf = self.vf_all[nameflag]
        self.magphys = self.magphys_all[nameflag]

        #for whatever reason, the vf output table is not entirely row-matched to the phot table. as such, the magphys output table will also not be entirely row-matched. I use np.argsort to correct this problem.
        
        ##I sort photometry tables in the same fashion...just to be absolutely certain of consistency. magphys and vf SHOULD be row-matched by virtue of their being part of the VF catalogs.
        #organizing in this way row-matches to vf and magphys, since their VFIDs ARE in ascending order.
        ind = np.argsort(self.phot_names_all)
        self.phot_names_all = self.phot_names_all[ind]
        self.phot_all = self.phot_all[ind]
        self.phot_maskbit_all = self.phot_maskbit_all[ind]

        #I will also be removing galaxies with negative AP06 values in either W1, W4, or NUV.
        #in W4, there are 1716 such galaxies; in NUV, 128 (though these galaxies ALSO have errant W4 values).
        #in future iterations of this code we will likely have a solution that does not involve flagging.
        flux_flag = (self.phot_all['FLUX_AP06_NUV']>0) & (self.phot_all['FLUX_AP06_W4']>0) & (self.phot_all['FLUX_AP06_W1']>0) & (self.phot_all['FLUX_AP06_W2']>0)

        #also...some galaxies with vcosmic<0. I'll just pluck these for now. we are also unsure about whether to proceed with vcosmic or use recession velocity. who knows what the future will hold? I sure as heck do not.
        vel_flag = (self.vf['Vcosmic']>0)

        self.vf = self.vf[flux_flag&vel_flag]
        self.phot = self.phot_all[flux_flag&vel_flag]
        self.phot_names = self.phot_names_all[flux_flag&vel_flag]
        self.magphys = self.magphys[flux_flag&vel_flag]
        self.phot_maskbit = self.phot_maskbit_all[flux_flag&vel_flag] #need to track [g,r,z] saturation, "medium" bits, prominent W1 foreground stars

    def logsfr_NUV(self):

            #COG_MTOT_* or FLUX_AP06_*, WHERE * = [FUV, NUV, g, r, z, W1, W2, W3, W4]
            # corresponding wavelengths --> [.15 um, .22 um, .472 um, .642 um, .926 um, 3.4 um, 4.6 um, 12 um, 22 um]

            self.vcosmic = self.vf['Vcosmic']  #km/s          #vcosmic? vr? will decide which to use later.
            dist_km = (self.vcosmic/68.)*3.08e19         #assume H0 = 68 km/s/Mpc; 1 Mpc = 3.08e19 km 
            dist_meters = dist_km*1.e3

            #and L_nu = f_nu * 4piD, where D is the flow-corrected distance
            wave_W4 = 22.e-6
            freq_W4 = 3.e8/wave_W4   #units are Hz

            w4_AP06_flux = self.phot['FLUX_AP06_W4']   #nanomaggies
            fW4_Jy = (w4_AP06_flux * 3.631e-6)  #Jy
            fW4_mks = fW4_Jy * 1.e-26   #W/m^2/Hz

            nuL_22 = freq_W4 * fW4_mks * (4*np.pi*dist_meters**2) #units are Watts
            nuL_22 = nuL_22 * 1.e7   #erg/s

            logSFR_W4 = np.log10(nuL_22) - 42.69

            NUV_AP06_flux = self.phot['FLUX_AP06_NUV']
            fNUV_Jy = NUV_AP06_flux * 3.631e-6   #Jy
            fNUV_mks = fNUV_Jy * 1.e-26   #W/m^2/Hz

            wave_NUV = 0.22e-6
            freq_NUV = 3.e8/wave_NUV   #units are Hz

            nuL_NUV = freq_NUV * fNUV_mks * (4*np.pi*dist_meters**2)*1.e7 #units are egs/s/Hz
            nuL_NUV_corrected = (nuL_NUV) + (2.26*nuL_22)

            logSFR_NUV_corrected = np.log10(nuL_NUV_corrected) - 43.17
            
            self.logSFR_color = logSFR_NUV_corrected

    def logmstar_C12(self):

        w1_AP06_flux = self.phot['FLUX_AP06_W1']    #nanomaggies
        w2_AP06_flux = self.phot['FLUX_AP06_W2']    #nanomaggies

        w1_magAB = 22.5 - 2.5 * np.log10(w1_AP06_flux)
        w2_magAB = 22.5 - 2.5 * np.log10(w2_AP06_flux)

        #m_Vega = m_AB - mag_offset, where
        #mag_offset_W1 = 2.699
        #mag_offset_W2 = 3.339
        #(source: https://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html#conv2flux)

        w1_magVega = w1_magAB - 2.699
        w2_magVega = w2_magAB - 3.339

        dist_mpc = (self.vcosmic/68.)
        dist_pc = dist_mpc*1.e6
        w1_absmag_AB = w1_magAB - 5*np.log10(dist_pc/10.)

        sun_absmag_AB = 5.92

        L_W1 = 10**(-0.4 * (w1_absmag_AB - sun_absmag_AB))   #units of solar luminosities
        logL_W1 = np.log10(L_W1)

        #use the equation which incorporates the W1-W2 WISE colors
        A0 = -0.376
        A1 = -1.053
        C_12 = w1_magVega - w2_magVega

        logMstar = A0 + A1*C_12 + logL_W1

        self.logMstar_C12 = logMstar

    def find_fitted_coefficients(self):

        gsaturation = (self.phot_maskbit & 2**2) == 2**2
        rsaturation = (self.phot_maskbit & 2**3) == 2**3
        zsaturation = (self.phot_maskbit & 2**4) == 2**4

        self.sat_flags = (gsaturation | rsaturation | zsaturation)
        self.medium_flag = (self.phot_maskbit & 2**11) == 2**11
        self.wisestar_flag = (self.phot_maskbit & 2**8) == 2**8

        fitflags = (self.sat_flags | self.medium_flag | self.wisestar_flag)  #flag galaxies which are either saturated, treated as a point source, or have a prominent W1 foreground star.

        #finding quantitative relationships and parameter uncertainties, but only for unflagged galaxies!
        p_sfr, cov_sfr = np.polyfit(self.magphys['logSFR_med'][~fitflags & (self.magphys['logMstar_med']>6.)],
                                    (self.logSFR_color[~fitflags & (self.magphys['logMstar_med']>6.)]),1,cov=True)

        p_mstar, cov_mstar = np.polyfit(self.magphys['logMstar_med'][~fitflags & (self.magphys['logMstar_med']>6.)],
                                        self.logMstar_C12[~fitflags & (self.magphys['logMstar_med']>6.)],1,cov=True)
        
        self.m_mass = p_mstar[0]
        self.m_mass_err = np.sqrt(cov_mstar[0][0])
        self.b_mass = p_mstar[1]
        self.b_mass_err = np.sqrt(cov_mstar[1][1])
        self.m_sfr = p_sfr[0]
        self.m_sfr_err = np.sqrt(cov_sfr[0][0])
        self.b_sfr = p_sfr[1]
        self.b_sfr_err = np.sqrt(cov_sfr[1][1])
        
        print('CurveFit Parameters log(Mstar_color) = mlog(Mstar_MAGPHYS) + b')
        print(f'm={np.round(self.m_mass,3)} +/- {np.round(self.m_mass_err,3)}, b={np.round(self.b_mass,3)} +/- {np.round(self.b_mass_err,3)}')
        print()
        print('CurveFit Parameters log(SFR_color) = mlog(SFR_MAGPHYS) + b')
        print(f'm={np.round(self.m_sfr,3)} +/- {np.round(self.m_sfr_err,3)}, b={np.round(self.b_sfr,3)} +/- {np.round(self.b_sfr_err,3)}')

    def output_altmagphys_table(self):

        self.magphys_mstar_alt = (self.logMstar_C12 - self.b_mass)/self.m_mass  #alternative magphys values
        self.magphys_sfr_alt = (self.logSFR_color - self.b_sfr)/self.m_sfr      #alternative magphys SFR values
        
        self.magphys_mstar_alt_err = np.sqrt((self.m_mass_err*self.logMstar_C12/self.m_mass)**2 + (np.sqrt(self.m_mass_err**2 + self.b_mass_err**2)*self.b_sfr/self.m_sfr)**2)
        self.magphys_sfr_alt_err = np.sqrt((self.m_sfr_err*self.logSFR_color/self.m_sfr)**2 + (np.sqrt(self.m_sfr_err**2 + self.b_sfr_err**2)*self.b_sfr/self.m_sfr)**2)
        
        good_phot_flag = np.zeros(len(self.vf_all),dtype=bool)
        
        #only true if vcosmic>0 and FLUX_AP06>0 in the relevant bands
        for i in range(len(self.vf_all)):
            good_phot_flag[i] = True if self.vf_all['VFID'][i] in self.vf['VFID'] else False   
            #(self.vf is so trimmed that only the good phot VFID matches are present)
        
        columns = names=['VFID', 'logMstar_C12', 'logSFR_color', 'MAGPHYS_logMstar_med', 'alt_logMstar_med', 
                         'alt_logMstar_med_err', 'MAGPHYS_logSFR_med', 'alt_logSFR_med', 'alt_logSFR_med_err',
                          'magphys_flag', 'unsaturated_flag', 'nowisestar_flag', 'nomedstar_flag', 'good_phot_flag']
        
        dtype=['S8',float,float,float,float,float,float,float,float,bool,bool,bool,bool,bool]
        
        #create empty table
        self.out_tab = Table(np.zeros((len(self.vf_all),len(columns))),dtype=dtype,names=columns)
        
        #add VFID column
        self.out_tab['VFID'] = self.vf_all['VFID']
        
        #add magphys columns
        self.out_tab['MAGPHYS_logMstar_med'] = self.magphys_all['logMstar_med']
        self.out_tab['MAGPHYS_logSFR_med'] = self.magphys_all['logSFR_med']
        self.out_tab['magphys_flag'] = self.magphys_all['magphysFlag']
        
        #add good_phot_flag column...
        self.out_tab['good_phot_flag'] = good_phot_flag
        
        #for every entry in trimmed phot catalog, find index of out_tab at which VFID = phot_names
        #then, assign appropriate values...
        for i in range(len(self.vf)):
            try:
                index = np.where(self.out_tab['VFID'] == self.vf['VFID'][i])[0][0]
                self.out_tab['logMstar_C12'][index] = self.logMstar_C12[i]
                self.out_tab['logSFR_color'][index] = self.logSFR_color[i]
                self.out_tab['alt_logMstar_med'][index] = self.magphys_mstar_alt[i]
                self.out_tab['alt_logMstar_med_err'][index] = self.magphys_mstar_alt_err[i]
                self.out_tab['alt_logSFR_med'][index] = self.magphys_sfr_alt[i]
                self.out_tab['alt_logSFR_med_err'][index] = self.magphys_sfr_alt_err[i]
                self.out_tab['unsaturated_flag'][index] = ~self.sat_flags[i]
                self.out_tab['nowisestar_flag'][index] = ~self.medium_flag[i]
                self.out_tab['nomedstar_flag'][index] = ~self.wisestar_flag[i]
                
            except:
                pass
                
                
        #then create table mask? use bool flag (self.out_tab==0.)

if __name__ == '__main__': 
    
    phot_table = altvals()
    phot_table.logsfr_NUV()
    phot_table.logmstar_C12()
    phot_table.find_fitted_coefficients()
    phot_table.output_altmagphys_table()
    
    print()
    phot_table.out_tab.write(homedir+'/Desktop/vf-altphot.fits',overwrite=True)
    print(f'Table written to {homedir}/Desktop/vf-altphot.fits')