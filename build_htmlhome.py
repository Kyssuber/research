import numpy as np
from matplotlib import pyplot as plt
import os
from astropy.table import Table


def build_htmlhome(sample):
    path = '/Users/kconger/'
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
        html.write('<font size="40">Galaxy (SNR>10) Img/Mask Mosaics </font>\n')

        html.write('<table><tr><th>Index</th><th>VFID</th><th>Image</th><th>Mask</th></tr>\n')

        for i in range(0,len(sample)):
            html.write('<tr><td>'+str(i)+'</td>\n')
            html.write('<td>'+str(sample['VFID'][i])+'</td>\n')
            html.write('<td><img src = "' + cutoutpath + str(sample['VFID'][i]) + '-mosaic.png' + '" height="auto" width = "100%"></img></td>\n')
            html.write('<td><img src = "' + maskpath + str(sample['VFID'][i]) + '-mask-mosaic.png' + '" height="auto" width = "100%"></img></td>\n')

        html.write('</tr></table>\n')
        html.write('<br /><br />\n')
        html.write('</html></body>\n')
        html.close()

if __name__ == '__main__':
    homedir = os.getenv("HOME")
    vf = Table.read(homedir+'/vfcut.fits')
    print('build_htmlhome(sample)')
    
    
