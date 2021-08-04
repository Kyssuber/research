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
        html.write('<font size="40"> WISESize GALFIT Data for VF Galaxies (SNR>15) </font>\n')

        html.write('<table><tr><th>Index</th><th>Prefix</th><th>RA</th><th>DEC</th><th>Image</th><th>Mask</th><th>2+ Sersic</th><th>Comments</th></tr>\n')

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
            html.write('<td> </td>\n')

        html.write('</tr></table>\n')
        html.write('<br /><br />\n')
        html.write('</html></body>\n')
        html.close()










def build_html_one(sample,i,ku_or_siena,psf_or_nopsf):
    if str(ku_or_siena) == 'ku':
        path = '/Users/k215c316/'
        mainpath = path+'main_ku.html'
    if str(ku_or_siena) == 'siena':
        path = '/Users/kconger/'
        mainpath = path+'main_siena.html'
    if str(psf_or_nopsf) == 'nopsf':
        galpath = '/mnt/astrophysics/kconger_wisesize/gal_output_nopsf/'
    if str(psf_or_nopsf) == 'psf':
        galpath = '/mnt/astrophysics/kconger_wisesize/gal_output/'
    htmlpath = '/mnt/astrophysics/kconger_wisesize/gal_html/'+str(sample[i]['VFID'])+'.html'

    mosaicpath = galpath+'/output_mosaics/'
    tab = ascii.read(galpath+'/set_one.txt')




    with open(htmlpath,'w') as html:
        html.write('<html><body>\n')
        html.write('<title>'+str(sample['prefix'][i])+'</title>')
        html.write('<style type="text/css">\n')
        html.write('.img-container{text-align: left;}')
        html.write('table, td, th {padding: 5px; text-align: center; border: 1px solid black;}\n')
        html.write('p {display: inline-block;;}\n')
        html.write('</style>\n')
        
        html.write('<font size="40">Central Galaxy: '+str(sample['prefix'][i]) +'</font><br /> \n')
        html.write('<a href='+mainpath+'>Return to Homepage</a></br /> \n')

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
        print(mosaicfile)


        
        if str(psf_or_nopsf) == 'nopsf':
            html.write('<div class='+'"'+'img-container'+'"> <!-- Block parent element --> <img src='+'"'+path+'output_mosaics/'+str(sample['prefix'][i])+'-unwise-w3-1Comp-galfit-out.png'+'" /><br /> \n')

        if str(psf_or_nopsf) == 'psf':
            html.write('<div class='+'"'+'img-container'+'"> <!-- Block parent element --> <img src='+'"'+mosaicfile+'" /><br /> \n')

        
        html.write('<div class='+'"'+'img-container'+'"> <!-- Block parent element --> <img src='+'"'+path+'LS_mosaics/'+str(sample['VFID'][i])+'-mosaic'+'.png'+'" /><br /> \n')

        html.write('<div class='+'"'+'img-container'+'"> <!-- Block parent element --> <img src='+'"'+path+'mask_mosaics/'+str(sample['VFID'][i])+'-mask-mosaic'+'.png'+'" /><br /> \n')

        index = np.where(sample['prefix'][i] == tab['galname'])[0]
        index = int(index)

        
        if tab['success_flag'][index] == 0:
            html.write('<font size="40">GALFIT Failed!</font>\n')

        else:

            #if the galaxy is a central galaxy, then determine the number of "perimeter" galaxies
            if str(tab['galname'][index])[0:8] in dummycat['central galaxy']:
                
                external_number = len(np.where(str(tab['galname'][index])[0:8] == dummycat['central galaxy'])[0])

            #write parameters of central galaxy as first row
                html.write('<font size="30">GALFIT Parameters</font><br />')
                html.write('<table><tr><th>Galname</th><th>Type</th><th>xc</th><th>xc_err</th><th>yc</th><th>yc_err</th><th>mag</th><th>mag_err</th><th>Re</th><th>Re_err</th><th>nser</th><th>nser_err</th><th>BA</th><th>BA_err</th><th>PA</th><th>PA+err</th></tr> \n')
                html.write('<tr><td>'+str(sample['VFID'][i])+'</td>')
                html.write('<td>Host</td> \n')
                html.write('<td>'+str(tab[index][0])+'</td>')
                html.write('<td>'+str(tab[index][1])+'</td>')
                html.write('<td>'+str(tab[index][2])+'</td>')
                html.write('<td>'+str(tab[index][3])+'</td>')
                html.write('<td>'+str(tab[index][4])+'</td>')
                html.write('<td>'+str(tab[index][5])+'</td>')
                html.write('<td>'+str(tab[index][6])+'</td>')
                html.write('<td>'+str(tab[index][7])+'</td>')
                html.write('<td>'+str(tab[index][8])+'</td>')
                html.write('<td>'+str(tab[index][9])+'</td>')
                html.write('<td>'+str(tab[index][10])+'</td>')
                html.write('<td>'+str(tab[index][11])+'</td>')
                html.write('<td>'+str(tab[index][12])+'</td>')
                html.write('<td>'+str(tab[index][13])+'</td></tr> \n')

                for f in range(0,external_number):
                    #num will ensure central galaxy is skipped
                    num=f+1
                    #only VFID galaxies (i.e., galaxies included in VF catalog) 
                    if tab[index+num]['galname'] in galvf['VFID']:
                        html.write('<tr><td>'+str(tab[index+num]['galname'])+'</td>')
                        html.write('<td>External</td> \n')
                        html.write('<td>'+str(tab[index+num][0])+'</td>')
                        html.write('<td>'+str(tab[index+num][1])+'</td>')
                        html.write('<td>'+str(tab[index+num][2])+'</td>')
                        html.write('<td>'+str(tab[index+num][3])+'</td>')
                        html.write('<td>'+str(tab[index+num][4])+'</td>')
                        html.write('<td>'+str(tab[index+num][5])+'</td>')
                        html.write('<td>'+str(tab[index+num][6])+'</td>')
                        html.write('<td>'+str(tab[index+num][7])+'</td>')
                        html.write('<td>'+str(tab[index+num][8])+'</td>')
                        html.write('<td>'+str(tab[index+num][9])+'</td>')
                        html.write('<td>'+str(tab[index+num][10])+'</td>')
                        html.write('<td>'+str(tab[index+num][11])+'</td>')
                        html.write('<td>'+str(tab[index+num][12])+'</td>')
                        html.write('<td>'+str(tab[index+num][13])+'</td> \n')
                    
            
            else:
                #if the galaxy is a central galaxy no "perimeter" galaxies, then go ahead and add one row to the page.
                html.write('<font size="30">GALFIT Parameters</font><br />')
                html.write('<table><th>Galname</th><th>Type</th><th>xc</th><th>xc_err</th><th>yc</th><th>yc_err</th><th>mag</th><th>mag_err</th><th>Re</th><th>Re_err</th><th>nser</th><th>nser_err</th><th>BA</th><th>BA_err</th><th>PA</th><th>PA_err</th></tr>')
                html.write('<tr><td>'+str(sample['VFID'][i])+'</td>')
                html.write('<td>Host</td> \n')
                html.write('<td>'+str(tab[index][0])+'</td>')
                html.write('<td>'+str(tab[index][1])+'</td>')
                html.write('<td>'+str(tab[index][2])+'</td>')
                html.write('<td>'+str(tab[index][3])+'</td>')
                html.write('<td>'+str(tab[index][4])+'</td>')
                html.write('<td>'+str(tab[index][5])+'</td>')
                html.write('<td>'+str(tab[index][6])+'</td>')
                html.write('<td>'+str(tab[index][7])+'</td>')
                html.write('<td>'+str(tab[index][8])+'</td>')
                html.write('<td>'+str(tab[index][9])+'</td>')
                html.write('<td>'+str(tab[index][10])+'</td>')
                html.write('<td>'+str(tab[index][11])+'</td>')
                html.write('<td>'+str(tab[index][12])+'</td>')
                html.write('<td>'+str(tab[index][13])+'</td>')
            
            html.write('</tr></table> \n')

            html.write('<a href=main.html>Return to Homepage</a></br /> \n')

        
            
        if i != len(sample)-1:
            html.write('<a href='+str(sample['VFID'][i+1])+'.html>Next Galaxy</a></br /> \n')         

        if i != 0:
            html.write('<a href='+str(sample['VFID'][i-1])+'.html>Previous Galaxy</a></br /> \n')
    
        html.write('<br /><br />\n')    
        html.write('</html></body>\n')     
            
        html.close()
        



def build_html_all(sample,ku_or_siena,psf_or_nopsf):
    for i in range(0,len(sample)):
        build_html_one(sample,i,ku_or_siena,psf_or_nopsf)
        print(sample['VFID'][i])
    build_htmlhome_galfit(sample,ku_or_siena,psf_or_nopsf)







def build_htmlhome_galfit(sample,ku_or_siena,psf_or_nopsf):
    
    if str(ku_or_siena) == 'ku':
        path = '/Users/k215c316/'
        htmlpath = '/mnt/astrophysics/kconger_wisesize/gal_html/main.html'
    if str(ku_or_siena) == 'siena':
        path = '/Users/kconger/'
        htmlpath = '/mnt/astrophysics/kconger_wisesize/gal_html/main.html'
    if str(psf_or_nopsf) == 'nopsf':
        galpath = '/mnt/astrophysics/kconger_wisesize/gal_output_nopsf/'
    if str(psf_or_nopsf) == 'psf':
        galpath = '/mnt/astrophysics/kconger_wisesize/gal_output/'
        

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
        html.write('<font size="40"> WISESize GALFIT Data for VF Galaxies (SNR>15) </font>\n')

        html.write('<table><tr><th>Index</th><th>LS Cutout</th><th>Prefix</th><th>Galaxy</th><th>RA</th><th>DEC</th><th>2+ Sersic</th><th>Comments</th></tr>\n')

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
    print('build_htmlhome(sample,ku_or_siena,htmlpage=False)')
    print('build_html_one(sample,i,ku_or_siena,psf_or_nopsf)')
    print('build_htmlhome_galfit(sample,ku_or_siena,psf_or_nopsf)')
    print('build_html_all(sample,ku_or_siena,psf_or_nopsf)')
    
