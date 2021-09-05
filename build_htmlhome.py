import numpy as np
from matplotlib import pyplot as plt
import os
from astropy.table import Table
from astropy.io import ascii
import glob


def build_htmlhome(sample,ku_or_siena,htmlpage=False):

    if str(ku_or_siena) == 'ku':
        path = '/Users/k215c316/'
    if str(ku_or_siena) == 'siena':
        path = '/mnt/astrophysics/kconger_wisesize/'
        
    htmlpath = '/mnt/astrophysics/kconger_wisesize/main_local.html'
    maskpath = path+'mask_mosaics/'
    cutoutpath = path+'LS_mosaics/'

    with open(htmlpath,'w') as html:
        html.write('<html><body>\n')
        html.write('<title>WISESize Project</title>\n')
        html.write('<body style="background-color:powderblue;">\n')
        html.write('<style type="text/css">\n')
        html.write('table, td, th {padding: 5px; text-align: center; border: 2px solid black;}\n')
        html.write('p {display: inline-block;;}\n')
        html.write('</style>\n')
        html.write('<font size="40"> WISESize GALFIT Data for VF Galaxies (SNR>10) </font>\n')

        html.write('<table><tr><th>Index</th><th>Prefix</th><th>RA</th><th>DEC</th><th>Image</th><th>Mask</th><th>2+ Sersic</th><th>Re_ratio (fixed/free PA,BA)</th><th>Comments</th></tr>\n')

        for i in range(0,len(sample)):
            html.write('<tr><td>'+str(i)+'</td>\n')
            html.write('<td>'+str(sample['prefix'][i])+'</td>\n')
            html.write('<td>'+str(sample['RA'][i])+'</td>\n')
            html.write('<td>'+str(sample['DEC'][i])+'</td>\n')
            html.write('<td><img src = "' + cutoutpath + str(sample['VFID'][i]) + '-mosaic.png' + '" height="auto" width = "100%"></img></td>\n')
            html.write('<td><img src = "' + maskpath + str(sample['VFID'][i]) + '-mask-mosaic.png' + '" height="auto" width = "100%"></img></td>\n')
            if vf['VFID'] in dummycat['central galaxy']:
                html.write('<td>Yes</td>')
            else:
                html.write('<td>No</td>')

            if vf['prefix'][i] in re_ratio['prefix']:
                index = np.where(re_ratio['prefix'] == vf['prefix'][i])[0]
                html.write('<td>'+ str('%.5f'%re_ratio['re_ratio'][index]) +' </td>\n')
            else:
                html.write('<td>N/A</td>\n')
            
            html.write('<td> </td>\n')

        html.write('</tr></table>\n')
        html.write('<br /><br />\n')
        html.write('</html></body>\n')
        html.close()







#once psf is properly established, may no longer be necessary to include no_psf row. Until then, do include, and be sure both are the correct length (which will not be so if GALFIT fails with ncomp>1 for psf and not non_psf, or vice versa)


def build_html_one(sample,i,ku_or_siena,paba_comparison=False):
    if str(ku_or_siena) == 'ku':
        path = '/Users/k215c316/'
    if str(ku_or_siena) == 'siena':
        path = '/Users/kconger/'

    galpath = '/mnt/astrophysics/kconger_wisesize/'


    htmlpath = '/mnt/astrophysics/kconger_wisesize/gal_html/'+str(sample[i]['VFID'])+'.html'

    mosaicpath = galpath+'/output_mosaics/'

    tab_nopsf = ascii.read(galpath+'gal_output_nopsf/set_one.txt')
    #print('tab_nopsf:',len(tab_nopsf))
    tab_psf = ascii.read(galpath+'/gal_output_psf/set_one.txt')
    #print('tab_psf:',len(tab_psf))


    #index329 not in tab_psf...
    index=np.where(tab_nopsf['galname'] == 'index329')[0]
    index=int(index)
    tab_nopsf.remove_row(index)
    
        
    
    if paba_comparison == True:
        tab_sga_psf = ascii.read(galpath+'gal_output_fixed_psf/set_one.txt')
        tab_sga_nopsf = ascii.read(galpath+'gal_output_fixed_nopsf/set_one.txt')

    with open(htmlpath,'w') as html:
        html.write('<html><body>\n')
        html.write('<title>'+str(sample['prefix'][i])+'</title>')
        html.write('<style type="text/css">\n')
        html.write('.img-container{text-align: left;}')
        html.write('table, td, th {padding: 5px; text-align: center; border: 1px solid black;}\n')
        html.write('p {display: inline-block;;}\n')
        html.write('</style>\n')
        
        html.write('<font size="40">Central Galaxy: '+str(sample['prefix'][i]) +'</font><br /> \n')

        html.write('<a href=main.html>Return to Homepage</a></br /> \n')

        if i != len(sample)-1:
            html.write('<a href='+str(sample['VFID'][i+1])+'.html>Next Galaxy</a></br /> \n') 

        if i != 0:
            html.write('<a href='+str(sample['VFID'][i-1])+'.html>Previous Galaxy</a></br /> \n')

        #define ncomp, since filename is predicated on number of Sersic obj in postage stamp
        if sample['VFID'][i] in dummycat['central galaxy']:
            ncomp = len(np.where(dummycat['central galaxy'] == sample['VFID'][i])[0]) + 1
        else:
            ncomp = 1


        mosaicfile = path + 'output_mosaics/'+str(sample['prefix'][i])+'-unwise-w3-{}Comp-galfit-out-conv.png'.format(ncomp)
        mosaicfile_sga = path + 'output_mosaics_fixed/'+str(sample['prefix'][i])+'-unwise-w3-{}Comp-galfit-out-conv.png'.format(ncomp)
        print(mosaicfile)


        html.write('<font size="30">GALFIT Output (PSF) with no parameters held fixed:</font><br /> \n')

        html.write('<div class='+'"'+'img-container'+'"> <!-- Block parent element --> <img src='+'"'+mosaicfile+'" /><br /> \n')

        html.write('<font size="30">GALFIT Output (PSF) with PA, BA parameters held fixed:</font><br /> \n')

        html.write('<div class='+'"'+'img-container'+'"> <!-- Block parent element --> <img src='+'"'+mosaicfile_sga+'" /><br /> \n') 
        
        html.write('<div class='+'"'+'img-container'+'"> <!-- Block parent element --> <img src='+'"'+path+'LS_mosaics/'+str(sample['VFID'][i])+'-mosaic.png" height="50%" width="50%" /><br /> \n')

        html.write('<div class='+'"'+'img-container'+'"> <!-- Block parent element --> <img src='+'"'+path+'mask_mosaics/'+str(sample['VFID'][i])+'-mask-mosaic'+'.png'+'" /><br /> \n')

        index = np.where(sample['prefix'][i] == tab_nopsf['galname'])[0]
        index = int(index)

        
        if tab_psf['success_flag'][index] == 0:
            html.write('<font size="40">GALFIT Failed!</font>\n')

        else:

            #if the galaxy is a central galaxy, then determine the number of "perimeter" galaxies
            if str(tab_nopsf['galname'][index])[0:8] in dummycat['central galaxy']:
                
                external_number = len(np.where(str(tab_nopsf['galname'][index])[0:8] == dummycat['central galaxy'])[0])

            #write parameters of central galaxy as first row
                html.write('<font size="30">GALFIT Parameters</font><br />')
                html.write('<table><tr><th>Galname</th><th>Type</th><th>Fixed/Free PA,BA</th><th>PSF?</th><th>xc</th><th>xc_err</th><th>yc</th><th>yc_err</th><th>mag</th><th>mag_err</th><th>Re</th><th>Re_err</th><th>nser</th><th>nser_err</th><th>BA</th><th>BA_err</th><th>PA</th><th>PA+err</th><th>Error</th></tr> \n')

                html.write('<tr><td>'+str(sample['VFID'][i])+'</td> \n')
                html.write('<td>Host</td> \n')
                html.write('<td>Free</td> \n')
                html.write('<td>No PSF</td> \n')
                html.write('<td>'+str(tab_nopsf[index][0])+'</td> \n')
                html.write('<td>'+str(tab_nopsf[index][1])+'</td> \n')
                html.write('<td>'+str(tab_nopsf[index][2])+'</td> \n')
                html.write('<td>'+str(tab_nopsf[index][3])+'</td> \n')
                html.write('<td>'+str(tab_nopsf[index][4])+'</td> \n')
                html.write('<td>'+str(tab_nopsf[index][5])+'</td> \n')
                html.write('<td>'+str(tab_nopsf[index][6])+'</td> \n')
                html.write('<td>'+str(tab_nopsf[index][7])+'</td> \n')
                html.write('<td>'+str(tab_nopsf[index][8])+'</td> \n')
                html.write('<td>'+str(tab_nopsf[index][9])+'</td> \n')
                html.write('<td>'+str(tab_nopsf[index][10])+'</td> \n')
                html.write('<td>'+str(tab_nopsf[index][11])+'</td> \n')
                html.write('<td>'+str(tab_nopsf[index][12])+'</td> \n')
                html.write('<td>'+str(tab_nopsf[index][13])+'</td> \n')
                html.write('<td>'+str(tab_nopsf[index][16])+'</td></tr> \n')

                html.write('<tr><td>--</td> \n')
                html.write('<td>--</td> \n')
                html.write('<td>-- </td> \n')
                html.write('<td>PSF</td> \n')
                html.write('<td>'+str(tab_psf[index][0])+'</td> \n')
                html.write('<td>'+str(tab_psf[index][1])+'</td> \n')
                html.write('<td>'+str(tab_psf[index][2])+'</td> \n')
                html.write('<td>'+str(tab_psf[index][3])+'</td> \n')
                html.write('<td>'+str(tab_psf[index][4])+'</td> \n')
                html.write('<td>'+str(tab_psf[index][5])+'</td> \n')
                html.write('<td>'+str(tab_psf[index][6])+'</td> \n')
                html.write('<td>'+str(tab_psf[index][7])+'</td> \n')
                html.write('<td>'+str(tab_psf[index][8])+'</td> \n')
                html.write('<td>'+str(tab_psf[index][9])+'</td> \n')
                html.write('<td>'+str(tab_psf[index][10])+'</td> \n')
                html.write('<td>'+str(tab_psf[index][11])+'</td> \n')
                html.write('<td>'+str(tab_psf[index][12])+'</td> \n')
                html.write('<td>'+str(tab_psf[index][13])+'</td> \n')
                html.write('<td>'+str(tab_psf[index][16])+'</td></tr> \n')

                if paba_comparison == True:

                    try:
                        index2 = np.where(str(tab_sga_psf['galname']) == str(sample['prefix'][0]))
                        print(index2)
                        index2 = int(index2[0])

                        index3 = np.where(tab_sga_nopsf['galname'] == sample['prefix'][i])[0]
                        index3 = int(index3)

                        html.write('<tr><td>'+str(sample['VFID'][i])+'</td> \n')
                        html.write('<td>Host</td> \n')
                        html.write('<td>Fixed</td> \n')
                        html.write('<td>NoPSF</td> \n')
                        html.write('<td>'+str(tab_sga_nopsf[index3][0])+'</td> \n')
                        html.write('<td>'+str(tab_sga_nopsf[index3][1])+'</td> \n')
                        html.write('<td>'+str(tab_sga_nopsf[index3][2])+'</td> \n')
                        html.write('<td>'+str(tab_sga_nopsf[index3][3])+'</td> \n')
                        html.write('<td>'+str(tab_sga_nopsf[index3][4])+'</td> \n')
                        html.write('<td>'+str(tab_sga_nopsf[index3][5])+'</td> \n')
                        html.write('<td>'+str(tab_sga_nopsf[index3][6])+'</td> \n')
                        html.write('<td>'+str(tab_sga_nopsf[index3][7])+'</td> \n')
                        html.write('<td>'+str(tab_sga_nopsf[index3][8])+'</td> \n')
                        html.write('<td>'+str(tab_sga_nopsf[index3][9])+'</td> \n')
                        html.write('<td>'+str(tab_sga_nopsf[index3][10])+'</td> \n')
                        html.write('<td>'+str(tab_sga_nopsf[index3][11])+'</td> \n')
                        html.write('<td>'+str(tab_sga_nopsf[index3][12])+'</td> \n')
                        html.write('<td>'+str(tab_sga_nopsf[index3][13])+'</td> \n')
                        html.write('<td>'+str(tab_sga_nopsf[index3][16])+'</td></tr> \n')

                        html.write('<tr><td>--</td> \n')
                        html.write('<td>--</td> \n')
                        html.write('<td>--</td> \n')
                        html.write('<td>PSF</td> \n')
                        html.write('<td>'+str(tab_sga_psf[index2][0])+'</td> \n')
                        html.write('<td>'+str(tab_sga_psf[index2][1])+'</td> \n')
                        html.write('<td>'+str(tab_sga_psf[index2][2])+'</td> \n')
                        html.write('<td>'+str(tab_sga_psf[index2][3])+'</td> \n')
                        html.write('<td>'+str(tab_sga_psf[index2][4])+'</td> \n')
                        html.write('<td>'+str(tab_sga_psf[index2][5])+'</td> \n')
                        html.write('<td>'+str(tab_sga_psf[index2][6])+'</td> \n')
                        html.write('<td>'+str(tab_sga_psf[index2][7])+'</td> \n')
                        html.write('<td>'+str(tab_sga_psf[index2][8])+'</td> \n')
                        html.write('<td>'+str(tab_sga_psf[index2][9])+'</td> \n')
                        html.write('<td>'+str(tab_sga_psf[index2][10])+'</td> \n')
                        html.write('<td>'+str(tab_sga_psf[index2][11])+'</td> \n')
                        html.write('<td>'+str(tab_sga_psf[index2][12])+'</td> \n')
                        html.write('<td>'+str(tab_sga_psf[index2][13])+'</td> \n')
                        html.write('<td>'+str(tab_sga_psf[index2][16])+'</td></tr> \n')

                    except:

                        print('no sga equivalent')
                    
                for f in range(0,external_number):
                    #num will ensure central galaxy is skipped
                    num=f+1
                    #only VFID galaxies (i.e., galaxies included in VF catalog) 
                    if tab_psf[index+num]['galname'] in galvf['VFID']:
                        html.write('<tr><td>'+str(tab_psf[index+num]['galname'])+'</td>')
                        html.write('<td>External</td> \n')
                        html.write('<td>Free</td> \n')
                        html.write('<td>NoPSF</td> \n')
                        html.write('<td>'+str(tab_nopsf[index+num][0])+'</td> \n')
                        html.write('<td>'+str(tab_nopsf[index+num][1])+'</td> \n')
                        html.write('<td>'+str(tab_nopsf[index+num][2])+'</td> \n')
                        html.write('<td>'+str(tab_nopsf[index+num][3])+'</td> \n')
                        html.write('<td>'+str(tab_nopsf[index+num][4])+'</td> \n')
                        html.write('<td>'+str(tab_nopsf[index+num][5])+'</td> \n')
                        html.write('<td>'+str(tab_nopsf[index+num][6])+'</td> \n')
                        html.write('<td>'+str(tab_nopsf[index+num][7])+'</td> \n')
                        html.write('<td>'+str(tab_nopsf[index+num][8])+'</td> \n')
                        html.write('<td>'+str(tab_nopsf[index+num][9])+'</td> \n')
                        html.write('<td>'+str(tab_nopsf[index+num][10])+'</td> \n')
                        html.write('<td>'+str(tab_nopsf[index+num][11])+'</td> \n')
                        html.write('<td>'+str(tab_nopsf[index+num][12])+'</td> \n')
                        html.write('<td>'+str(tab_nopsf[index+num][13])+'</td> \n')
                        html.write('<td>'+str(tab_nopsf[index+num][16])+'</td></tr> \n')

                        html.write('<tr><td>'+str(tab_nopsf[index+num]['galname'])+'</td>')
                        html.write('<td>--</td> \n')
                        html.write('<td>--</td> \n')
                        html.write('<td>PSF</td> \n')
                        html.write('<td>'+str(tab_psf[index+num][0])+'</td> \n')
                        html.write('<td>'+str(tab_psf[index+num][1])+'</td> \n')
                        html.write('<td>'+str(tab_psf[index+num][2])+'</td> \n')
                        html.write('<td>'+str(tab_psf[index+num][3])+'</td> \n')
                        html.write('<td>'+str(tab_psf[index+num][4])+'</td> \n')
                        html.write('<td>'+str(tab_psf[index+num][5])+'</td> \n')
                        html.write('<td>'+str(tab_psf[index+num][6])+'</td> \n')
                        html.write('<td>'+str(tab_psf[index+num][7])+'</td> \n')
                        html.write('<td>'+str(tab_psf[index+num][8])+'</td> \n')
                        html.write('<td>'+str(tab_psf[index+num][9])+'</td> \n')
                        html.write('<td>'+str(tab_psf[index+num][10])+'</td> \n')
                        html.write('<td>'+str(tab_psf[index+num][11])+'</td> \n')
                        html.write('<td>'+str(tab_psf[index+num][12])+'</td> \n')
                        html.write('<td>'+str(tab_psf[index+num][13])+'</td> \n')
                        html.write('<td>'+str(tab_psf[index+num][16])+'</td> \n')
                    
            
            else:
                #if the galaxy is a central galaxy with no "perimeter" galaxies, then go ahead and add one row to the page.
                html.write('<font size="30">GALFIT Parameters</font><br />')
                html.write('<table><th>Galname</th><th>Type</th><th>Fixed/Free PA,BA</th><th>PSF?</th><th>xc</th><th>xc_err</th><th>yc</th><th>yc_err</th><th>mag</th><th>mag_err</th><th>Re</th><th>Re_err</th><th>nser</th><th>nser_err</th><th>BA</th><th>BA_err</th><th>PA</th><th>PA_err</th><th></th><th>Error</th></tr> \n')

                html.write('<tr><td>'+str(sample['VFID'][i])+'</td> \n')
                html.write('<td>Host</td> \n')
                html.write('<td>Free</td> \n')
                html.write('<td><NoPSF></td> \n')
                html.write('<td>'+str(tab_nopsf[index][0])+'</td> \n')
                html.write('<td>'+str(tab_nopsf[index][1])+'</td> \n')
                html.write('<td>'+str(tab_nopsf[index][2])+'</td> \n')
                html.write('<td>'+str(tab_nopsf[index][3])+'</td> \n')
                html.write('<td>'+str(tab_nopsf[index][4])+'</td> \n')
                html.write('<td>'+str(tab_nopsf[index][5])+'</td> \n')
                html.write('<td>'+str(tab_nopsf[index][6])+'</td> \n')
                html.write('<td>'+str(tab_nopsf[index][7])+'</td> \n')
                html.write('<td>'+str(tab_nopsf[index][8])+'</td> \n')
                html.write('<td>'+str(tab_nopsf[index][9])+'</td> \n')
                html.write('<td>'+str(tab_nopsf[index][10])+'</td> \n')
                html.write('<td>'+str(tab_nopsf[index][11])+'</td> \n')
                html.write('<td>'+str(tab_nopsf[index][12])+'</td> \n')
                html.write('<td>'+str(tab_nopsf[index][13])+'</td> \n')
                html.write('<td>'+str(tab_nopsf[index][16])+'</td> \n')

                html.write('<tr><td>--</td> \n')
                html.write('<td>--</td> \n')
                html.write('<td>--</td> \n')
                html.write('<td><PSF></td> \n')
                html.write('<td>'+str(tab_psf[index][0])+'</td> \n')
                html.write('<td>'+str(tab_psf[index][1])+'</td> \n')
                html.write('<td>'+str(tab_psf[index][2])+'</td> \n')
                html.write('<td>'+str(tab_psf[index][3])+'</td> \n')
                html.write('<td>'+str(tab_psf[index][4])+'</td> \n')
                html.write('<td>'+str(tab_psf[index][5])+'</td> \n')
                html.write('<td>'+str(tab_psf[index][6])+'</td> \n')
                html.write('<td>'+str(tab_psf[index][7])+'</td> \n')
                html.write('<td>'+str(tab_psf[index][8])+'</td> \n')
                html.write('<td>'+str(tab_psf[index][9])+'</td> \n')
                html.write('<td>'+str(tab_psf[index][10])+'</td> \n')
                html.write('<td>'+str(tab_psf[index][11])+'</td> \n')
                html.write('<td>'+str(tab_psf[index][12])+'</td> \n')
                html.write('<td>'+str(tab_psf[index][13])+'</td> \n')
                html.write('<td>'+str(tab_psf[index][16])+'</td> \n')
                
                if paba_comparison == True:

                    try:
                        index2 = np.where(tab_sga_psf['galname'] == sample['prefix'][i])[0]
                        index2 = int(index2)
                        index3 = np.where(tab_sga_nopsf['galname'] == sample['prefix'][i])[0]
                        index3 = int(index3)

                        html.write('<tr><td>'+str(sample['VFID'][i])+'</td>')
                        html.write('<td>Host</td> \n')
                        html.write('<td>Fixed</td> \n')
                        html.write('<td>NoPSF</td> \n')
                        html.write('<td>'+str(tab_sga_nopsf[index3][0])+'</td>')
                        html.write('<td>'+str(tab_sga_nopsf[index3][1])+'</td>')
                        html.write('<td>'+str(tab_sga_nopsf[index3][2])+'</td>')
                        html.write('<td>'+str(tab_sga_nopsf[index3][3])+'</td>')
                        html.write('<td>'+str(tab_sga_nopsf[index3][4])+'</td>')
                        html.write('<td>'+str(tab_sga_nopsf[index3][5])+'</td>')
                        html.write('<td>'+str(tab_sga_nopsf[index3][6])+'</td>')
                        html.write('<td>'+str(tab_sga_nopsf[index3][7])+'</td>')
                        html.write('<td>'+str(tab_sga_nopsf[index3][8])+'</td>')
                        html.write('<td>'+str(tab_sga_nopsf[index3][9])+'</td>')
                        html.write('<td>'+str(tab_sga_nopsf[index3][10])+'</td>')
                        html.write('<td>'+str(tab_sga_nopsf[index3][11])+'</td>')
                        html.write('<td>'+str(tab_sga_nopsf[index3][12])+'</td>')
                        html.write('<td>'+str(tab_sga_nopsf[index3][13])+'</td> \n')
                        html.write('<td>'+str(tab_sga_nopsf[index3][16])+'</td></tr> \n')

                        html.write('<tr><td>'+str(sample['VFID'][i])+'</td>')
                        html.write('<td>Host</td> \n')
                        html.write('<td>Fixed</td> \n')
                        html.write('<td>NoPSF</td> \n')
                        html.write('<td>'+str(tab_sga_psf[index2][0])+'</td>')
                        html.write('<td>'+str(tab_sga_psf[index2][1])+'</td>')
                        html.write('<td>'+str(tab_sga_psf[index2][2])+'</td>')
                        html.write('<td>'+str(tab_sga_psf[index2][3])+'</td>')
                        html.write('<td>'+str(tab_sga_psf[index2][4])+'</td>')
                        html.write('<td>'+str(tab_sga_psf[index2][5])+'</td>')
                        html.write('<td>'+str(tab_sga_psf[index2][6])+'</td>')
                        html.write('<td>'+str(tab_sga_psf[index2][7])+'</td>')
                        html.write('<td>'+str(tab_sga_psf[index2][8])+'</td>')
                        html.write('<td>'+str(tab_sga_psf[index2][9])+'</td>')
                        html.write('<td>'+str(tab_sga_psf[index2][10])+'</td>')
                        html.write('<td>'+str(tab_sga_psf[index2][11])+'</td>')
                        html.write('<td>'+str(tab_sga_psf[index2][12])+'</td>')
                        html.write('<td>'+str(tab_sga_psf[index2][13])+'</td> \n')
                        html.write('<td>'+str(tab_sga_psf[index2][16])+'</td></tr> \n')

                    except:
                        print('no sga equivalent')
                
                
            html.write('</tr></table> \n')

            html.write('<a href=main.html>Return to Homepage</a></br /> \n')

        
            
        if i != len(sample)-1:
            html.write('<a href='+str(sample['VFID'][i+1])+'.html>Next Galaxy</a></br /> \n')         

        if i != 0:
            html.write('<a href='+str(sample['VFID'][i-1])+'.html>Previous Galaxy</a></br /> \n')
    
        html.write('<br /><br />\n')    
        html.write('</html></body>\n')     
            
        html.close()
        



def build_html_all(sample,ku_or_siena,paba_comparison):
    for i in range(0,len(sample)):
        build_html_one(sample,i,ku_or_siena,paba_comparison)
        print(sample['VFID'][i])
    build_htmlhome_galfit(sample,ku_or_siena)







def build_htmlhome_galfit(sample,ku_or_siena):
    
    if str(ku_or_siena) == 'ku':
        path = '/Users/k215c316/'
        htmlpath = '/mnt/astrophysics/kconger_wisesize/gal_html/main.html'
    if str(ku_or_siena) == 'siena':
        path = '/Users/kconger/'
        htmlpath = '/mnt/astrophysics/kconger_wisesize/gal_html/main.html'
    
        

    stamppath = path+'LS_cutouts/'
    galhtmlpath = path+'gal_html/'

    with open(htmlpath,'w') as html:
        html.write('<html><body>\n')
        html.write('<title>WISESize Project</title>\n')
        html.write('<body style="background-color:powderblue;">\n')
        html.write('<style type="text/css">\n')
        html.write('table, td, th {padding: 5px; text-align: center; border: 2px solid black;}\n')
        html.write('p {display: inline-block;;}\n')
        html.write('</style>\n')
        html.write('<font size="40"> WISESize GALFIT Data for VF Galaxies (SNR>10) </font>\n')

        html.write('<table><tr><th>Index</th><th>LS Cutout</th><th>Prefix</th><th>Galaxy</th><th>RA</th><th>DEC</th><th>2+ Sersic</th><th>Re_ratio (fixed/free PA,BA)</th><th>Comments</th></tr>\n')

        for i in range(0,len(sample)):
            html.write('<tr><td>'+str(i)+'</td>\n')
            html.write('<td><img src = "' + stamppath + str(sample['VFID'][i]) + '.jpg' + '" height="auto" width = "100%"></img></td>\n')
            html.write('<td>'+str(sample['prefix'][i])+'</td>\n')
############for the below entry, do not run until build_htmlpages runs
            html.write('<td><a href='+galhtmlpath+str(sample['VFID'][i])+'.html>'+str(sample['objname'][i])+'</a></td>\n')
            html.write('<td>'+str(sample['RA'][i])+'</td>\n')
            html.write('<td>'+str(sample['DEC'][i])+'</td>\n')

            if vf['VFID'][i] in dummycat['central galaxy']:
                html.write('<td>Yes</td>')
            else:
                html.write('<td>--</td>')

            if vf['prefix'][i] in re_ratio['prefix']:
                index = np.where(re_ratio['prefix'] == vf['prefix'][i])[0]
                html.write('<td>'+ str('%.5f'%re_ratio['re_ratio'][index]) +' </td>\n')
            else:
                html.write('<td>N/A</td>\n')

            html.write('<td> </td>\n')

        html.write('</tr></table>\n')
        html.write('<br /><br />\n')
        html.write('</html></body>\n')
        html.close()






if __name__ == '__main__':
    homedir = os.getenv("HOME")
    vf = Table.read(homedir+'/vfcut.fits',format='ascii')
    dummycat = Table.read(homedir+'/dummycat.fits',format='ascii')
    galvf = Table.read('/mnt/astrophysics/kconger_wisesize/github/research/galfitcut.fits',format='ascii')
    re_ratio = Table.read(homedir+'/re_ratio.fits',format='ascii')
    print('build_htmlhome(sample,ku_or_siena,htmlpage=False)')
    print('build_html_one(sample,i,ku_or_siena,paba_comparison=False)')
    print('build_htmlhome_galfit(sample,ku_or_siena)')
    print('build_html_all(sample,ku_or_siena,paba_comparison=False)')
    print(' ')
    print('paba_comparison: True, False')
    print(' ')
    print('Be sure gal_html directory exists before running.')
