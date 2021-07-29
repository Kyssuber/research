import numpy as np
from matplotlib import pyplot as plt
import os
from astropy.table import Table


def build_htmlhome(sample,siena_or_ku,htmlpage=False):

    if str(siena_or_ku) == 'ku':
        path = '/Users/k215c316/'
    if str(siena_or_ku) == 'siena':
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
            html.write('<td>+str(sample['RA'][i])+'</td>\n')
            html.write('<td>+str(sample['DEC'][i])+'</td>\n')
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










def build_htmlpage_one(sample,i,ku_or_siena,psf_or_nopsf):
    if str(siena_or_ku) == 'ku':
        path = '/Users/k215c316/'
    if str(siena_or_ku) == 'siena':
        path = '/Users/kconger/'
    if str(psf_or_nopsf) == 'nopsf':
        galpath = '/mnt/astrophysics/kconger_wisesize/gal_output_nopsf/'
    if str(psf_or_nopsf) == 'psf':
        galpath = '/mnt/astrophysics/kconger_wisesize/gal_output/'

    mosaicpath = galpath+'/output_mosaics/'
    tab = ascii.read(galpath+'/set_one.txt')

    with open(path,'w') as html:
        html.write('<html><body>\n')
        html.write('<title>'+str(sample['prefix'][i])+'</title>')
        html.write('<style type="text/css">\n')
        html.write('.img-container{text-align: left;}')
        html.write('table, td, th {padding: 5px; text-align: center; border: 1px solid black;}\n')
        html.write('p {display: inline-block;;}\n')
        html.write('</style>\n')
        
        html.write('<font size="40">Central Galaxy: '+str(vfmain['prefix'][i]) +'</font><br /> \n')
        html.write('<a href=main.html>Return to Homepage</a></br /> \n')

        if i != len(sample)-1:
            html.write('<a href='+str(sample['VFID'][i+1])+'.html>Next Galaxy</a></br /> \n') 

        if i != 0:
            html.write('<a href='+str(sample['VFID'][i-1])+'.html>Previous Galaxy</a></br /> \n')



        html.write('<div class='+'"'+'img-container'+'"> <!-- Block parent element --> <img src='+'"'+path+'output_mosaics/'+str(sample['prefix'][i])+'-unwise-w3-1Comp-galfit-out.png'+'" /><br /> \n')
        html.write('<div class='+'"'+'img-container'+'"> <!-- Block parent element --> <img src='+'"'+path+'LS_mosaics/'+str(sample['VFID'][i])+'-mosaic'+'.png'+'" /><br /> \n')

        index = np.where(t['galname'] == sample['VFID'])[0]


        
        if tab['success_flag'][index] == 0:
            html.write('<font size="40">GALFIT Failed!</font>\n')



        else:

            if sample['VFID'][index] in dummycat['central galaxy']:
                
                external_number = len(np.where(sample['VFID'][index] == dummycat['central galaxy'])[0])
                
                html.write('<font size="30">GALFIT Parameters</font><br />')
                html.write('<table><tr><th>Galname</th><th>Type</th><th>xc</th><th>xc_err</th><th>yc</th><th>yc_err</th><th>mag</th><th>mag_err</th><th>Re</th><th>Re_err</th><th>nser</th><th>nser_err</th><th>BA</th><th>BA_err</th><th>PA</th><th>PA+err</th></th> \n')
                html.write('<td>'+str(sample['VFID'][i])+'</td>')
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
                html.write('<td>'+str(tab[index][13])+'</td> \n')

                for i in range(0,external_number):
                    num=i+1
                    if tab[index+1]['galname'] in galvf['VFID']:
                        html.write('<td>'+str(tab[index+i]['galname'])+'</td>')
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
                html.write('<font size="30">GALFIT Parameters</font><br />')
                html.write('<table><tr><th>xc</th><th>xc_err</th><th>yc</th><th>yc_err</th><th>mag</th><th>mag_err</th><th>Re</th><th>Re_err</th><th>nser</th><th>nser_err</th><th>BA</th><th>BA_err</th><th>PA</th><th>PA_err</th></tr>')
                html.write('<td>'+str(sample['VFID'][i])+'</td>')
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
        







if __name__ == '__main__':
    homedir = os.getenv("HOME")
    vf = Table.read(homedir+'/vfcut.fits',format='ascii')
    dummycat = Table.read(homedir+'/dummycat.fits',format='ascii')
    galvf = Table.read('/mnt/astrophysics/kconger_wisesize/github/research/galfitcut.fits',format='ascii')
    print('build_htmlhome(sample,ku_or_siena,htmlpage=False)')
    print('build_html_one(sample,ku_or_siena,psf_or_nopsf)')
    
    
