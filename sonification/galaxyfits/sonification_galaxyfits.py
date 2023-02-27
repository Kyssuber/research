'''
GOAL: 
-- sonify a single 2D galaxy cutout using an input parameter file.
-- create vertical bar animation that sweeps across and WILL BE in sync (in principle) with overlayed sonified data (separate script)
Output: 
-- in desired directory, will receive 
-----.mid and data.csv, both comprising versions of the data "trimmed" in the define_strips() function below. 
-----.mp4 of vertical bar animation
-----The .mid converted to a .wav file via FluidSynth (midi2audio)
'''

import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
from astropy.table import Table
from astropy.nddata import Cutout2D

import csv
import pandas as pd

from midiutil import MIDIFile
from audiolazy import str2midi

from matplotlib import animation, rc
from IPython.display import HTML
from matplotlib.lines import Line2D
from mido import MidiFile
from midi2audio import FluidSynth

from scipy.stats import scoreatpercentile
from astropy.visualization import simple_norm


def define_strips(fits_type, filepath, gal_name, filename, vf, r25_multiple_x, r25_multiple_y):

    hdu1=fits.open(filepath+filename)
    im = fits.getdata(filepath+filename)
    
    #if name of galaxy is VFID (in 2D cutout filename), taken from the original filenames /mnt/astrophysics/wisesize on virtual machine, then run the following:
    if fits_type == 'VFID':
        header = hdu1[0].header
        vf_index = np.where(gal_name == vf['VFID_V1'])[0]  #remember that these img name VFIDs are from the v1 catalog!
        r25 = vf['radius'][vf_index]  #optical, in arcseconds
        print('all savefiles will include the galaxy VFID',gal_name)
    
    #if instead the name of the galaxy is the objname (in 2D cutout filename), taken from /mnt/virgofilaments-data on virtual machine, then run the following:
    elif fits_type == 'objname':
        header = hdu1[1].header
        vf_index = np.where(gal_name == vf['objname'])[0]
        r25 = vf['radius'][vf_index]  #optical, in arcseconds
        print('all savefiles will include the galaxy objname',gal_name)
    
    #if neither...then the default bar length will be 3/4 of the image height.
    else:
        return
    
    if (fits_type=='VFID') | (fits_type=='objname'):
        #use header information to generalize conversion from arcseconds to pixels (for 12-micron, 2.75 arcsec per px),
        len_image_arcsec = np.abs(header['NAXIS1']*header['CD1_1'])*3600
        arcsec_per_pixel = len_image_arcsec/header['NAXIS1']
        r25 = r25/arcsec_per_pixel

        lower_bound_y = int(header['CRPIX2']-(r25_multiple_y*r25))
        upper_bound_y = int(header['CRPIX2']+(r25_multiple_y*r25))

        lower_bound_x = int(header['CRPIX1']-(r25_multiple_x*r25))
        upper_bound_x = int(header['CRPIX1']+(r25_multiple_x*r25))
    
    else:
        lower_bound_y = int(im.shape[1]-(im.shape[1]*3/4))
        upper_bound_y = int(im.shape[1]+(im.shape[1]*3/4))
        
        lower_bound_x = int(im.shape[0]-(im.shape[0]*3/4))
        upper_bound_x = int(im.shape[0]+(im.shape[0]*3/4))

    #first isolate pixels in this APPROXIMATE band.
    band = im[:,lower_bound_y:upper_bound_y]
    #then each ith vertical strip will be band[i,:] (all y pixels in ith strip)

    #create empty array for 1px strips
    strips=[]
    #for all vertical strips spanning the desired x-axis of the image
    for i in range(lower_bound_x,upper_bound_x):
        strips.append(band[i,:])
    
    return(im,lower_bound_y,upper_bound_y,lower_bound_x,upper_bound_x,strips,gal_name)

def create_csv_file(strips, pathname, csv_savename):
    # open the file in the write mode
    with open(pathname+csv_savename, 'w') as f:
        # create the csv writer
        writer = csv.writer(f)
        # write a row to the csv file
        for i in strips:
            writer.writerow(i)
    f.close()
    
def load_csv_file(csv_savename,pathname):
    #load data as a pandas dataframe
    df = pd.read_csv(pathname+csv_savename,header=None)
    #add column to df representing mean of each pixel strip
    df = df.assign(mean=df.mean(axis=1))
    return(df)

'''a typical sonification mapping function
mapping value(s) from one range to another range
returns float(s)
e.g., want to map 5, which came from data ranging from 1 to 10, to a range from 100 to 200
map_value(5,1,10,100,200) '''

def map_value(value, min_value, max_value, min_result, max_result):
    result = min_result + (value - min_value)/(max_value - min_value)*(max_result - min_result)
    return result

'''now. compress time. how many strips per beat of music?
recall that, for the example .fits, there will be 500 strips per band

ONE OPTION: SET CONVERSION FACTOR TO SUPPRESS "TIME"
(we have # strips per beat, but what will be the beat duration?)'''

#y_data will be df['mean']
def compress_time(y_data,strips_per_beat):
    #rescale strip number to beats
    t_data = np.arange(0,len(y_data),1) / strips_per_beat
    duration_beats=np.max(t_data)
    print('Duration:',duration_beats, 'beats')
    #one beat = one quarter note
    return(t_data)

#many, MANY points are around or below 0.1 once normalized, and they may likely be mapped to the same music note
#one option: scale data points by **0.5
def normalize_scale(y_data,y_scale):
    y_data = map_value(y_data,min(y_data),max(y_data),0,1)
    y_data_scaled = y_data**y_scale
    return(y_data,y_data_scaled)

def visualize_data(t_data,y_data):
    plt.figure(figsize=(6,4))
    plt.scatter(t_data,y_data,s=y_data*100,alpha=0.5,edgecolors='black')
    plt.scatter(x=[],y=[],color='white',label='data scaled by **0.5')
    plt.xlabel('x-coordinate pixel',fontsize=12)
    plt.ylabel('normalized pixel value (12'+r'$\mu$m)',fontsize=12)
    plt.grid(alpha=0.5)
    plt.legend()
    plt.show()

def note_to_midi(y_data_scaled,note_names):
    note_midis = [str2midi(n) for n in note_names]  #list of midi note numbers
    n_notes = len(note_midis)
    print('Resolution:',n_notes,'notes')

    #MAPPING DATA TO MIDIS!

    midi_data = []
    #for every data point (which I represent with len(t_data))
    for i in range(len(t_data)):
        #map y_data_scaled values such that smallest px is lowest note, largest px is highest note
        #assigns midi note number to whichever y_data_scaled[i], once mapped, is nearest to
        note_index = round(map_value(y_data_scaled[i],0,1,0,len(note_midis)-1))
        midi_data.append(note_midis[note_index])

    return(midi_data)

def map_note_velocities(y_data_scaled,vel_min,vel_max):
    #map data to note velocities
    #like midi numbers, ranges from 0 to 127
    #don't want vel_min=0, as this value will generate a *very* faint sound
    vel_data = []
    for i in range(len(y_data_scaled)):
        note_velocity = round(map_value(y_data_scaled[i],0,1,vel_min,vel_max)) #larger values, heavier sound
        vel_data.append(note_velocity)
    
    return(vel_data)

#save data as midi file
def save_midi_file(t_data, midi_data, volume, bpm, midi_savename, VFID=None):

    #create midi file object, add tempo
    my_midi_file = MIDIFile(1) #one track
    my_midi_file.addTempo(track=0,time=0,tempo=bpm) #only one track, so track=0th track; begin at time=0, tempo is bpm

    #add midi notes to file
    for i in range(len(t_data)):
        my_midi_file.addNote(track=0,channel=0,pitch=midi_data[i],time=t_data[i],duration=2,volume=vel_data[i])

    #write file
    with open(midi_savename,"wb") as f:
        my_midi_file.writeFile(f)

        
########################################################################


def get_midi_length(path_to_midi,midi_savename):    
    mid = MidiFile(midi_savename)
    time = mid.length
    #there is a slight time delay between end of animation and end of audio
    #'duration' of each midi note is 2 sec, but I think time-3 seconds is preferable (at least for VFID0047)
    #but still not flawless - a bit of asynchrony remains.
    return(time-3)

def create_animation_midi(t_data, midi_data, time, vid_name_midi):
        
    fig=plt.figure(figsize=(10,6))

    plt.scatter(t_data,midi_data,s=y_data*100,alpha=0.5,edgecolors='black')
    plt.xlabel('RA scaled (deg)',fontsize=12)
    plt.ylabel('Midi note numbers',fontsize=12)

    X_MIN = int(np.min(t_data))
    X_MAX = int(np.max(t_data))
    Y_MIN = int(np.min(midi_data))
    Y_MAX = int(np.max(midi_data))
    X_VALS = np.arange(X_MIN, X_MAX+1,0.05) # possible pixel values for the line, increments 0.05
    
    def update_line(num, line):
        i = X_VALS[num]
        line.set_data([i, i], [Y_MIN-5, Y_MAX+5])
        return line, 
    
    line, = plt.plot([], [], lw=2)

    l,v = plt.plot(X_MIN, Y_MIN, X_MAX, Y_MAX, linewidth=2, color= 'red')

    line_anim = animation.FuncAnimation(fig, update_line, frames=len(X_VALS), fargs=(l, ),
                                        blit=True) 
    FFWriter = animation.FFMpegWriter()
    line_anim.save(vid_name_midi,fps=len(X_VALS)/time)
    plt.close()

def create_animation_cutout(fits_type, filepath, filename, im, lower_bound_x, upper_bound_x, lower_bound_y, upper_bound_y, time, vid_name_cutout):
    
    '''
    percentile1 = min percentile for stretch of image and model
    percentile2 = max percentile for stretch of image and model
    cmap = colormap, default is viridis
    '''
    
    hdu1=fits.open(filepath+filename)
    im = fits.getdata(filepath+filename)
    
    if fits_type == 'VFID':
        header = hdu1[0].header

    if fits_type == 'objname':
        header = hdu1[1].header
    
    X_MIN = int(lower_bound_x)
    X_MAX = int(upper_bound_x)
    Y_MIN = int(lower_bound_y)
    Y_MAX = int(upper_bound_y)
    X_VALS = np.arange(X_MIN, X_MAX+1,0.05) # possible pixel values for the line, increments 0.05   
    #X_VALS = np.flip(X_VALS)
    
    #create cutout; image center is galaxy center for these 500x500 px iamges; 
    #length and width are xmax-xmin, ymax-ymin respectively
    cutout = Cutout2D(im,(header['CRPIX1'],header['CRPIX2']), (X_MAX-X_MIN,Y_MAX-Y_MIN))
   
    percentile1=.5
    percentile2=99.5
    cmap='viridis'
    
    v1 = scoreatpercentile(cutout.data,percentile1)
    v2 = scoreatpercentile(cutout.data,percentile2)
    norm =  simple_norm(cutout.data,'asinh',min_cut=v1,max_cut=v2)
    
    plt.figure()
    plt.imshow(im,origin='lower',cmap=cmap,norm=norm)
    plt.xlim(X_MIN-10, X_MAX+10)
    plt.ylim(Y_MIN-10, Y_MAX+10)

    circle1 = plt.Circle((250, 250), upper_bound_x-250, color='r',alpha=0.2)
    fig = plt.gcf()
    ax = fig.gca()
    ax.add_patch(circle1)

    def update_line(num, line):
        i = X_VALS[num]
        line.set_data([i, i], [Y_MIN, Y_MAX])
        return line, 

    line, = plt.plot([], [], lw=2)

    l,v = plt.plot(X_MIN, Y_MIN, X_MAX, Y_MAX, linewidth=2, color= 'red')

    line_anim = animation.FuncAnimation(fig, update_line, frames=len(X_VALS), fargs=(l, ),
                                        blit=True) 
    FFWriter = animation.FFMpegWriter()
    line_anim.save(vid_name_cutout,fps=len(X_VALS)/time)
    plt.close()

def convert_to_wav(midi_savename,wav_savename,soundfont,gain=2):
    #the .sf2 file was downloaded from internet; can select alternative .sf2 files which act as converter maps between midi values and notes.
    #gain = 2; helps control volume. if run without, then resulting file is too soft. :-( 
    
    fs = FluidSynth(sound_font=soundfont, gain=gain)
    fs.midi_to_audio(midi_savename, wav_savename)  
    
if __name__ == '__main__':
    
    
    if '-h' in sys.argv or '--help' in sys.argv:
        print("Usage: %s [-param_path (path to location of the parameter file, no single or double quotation marks)] [-param_file (name of parameter file, no single or double quotation marks)]")
        sys.exit(1)
    
    if '-param_path' in sys.argv:
        p = sys.argv.index('-param_path')
        path_to_paramfile = str(sys.argv[p+1])
    
    if '-param_file' in sys.argv:
        p = sys.argv.index('-param_file')
        param_file = str(sys.argv[p+1])
    
    homedir = os.getenv("HOME")
    
    #create dictionary with keywords and values, from parameter .txt file

    param_dict = {}
    with open(homedir+path_to_paramfile+param_file) as f:
        for line in f:
            try:
                key = line.split()[0]
                val = line.split()[1]
                param_dict[key] = val
            except:
                continue
    
    #extracting parameters and assigning to variables
    
    VFID_or_objname = param_dict['VFID_or_objname']
    pathname = homedir+param_dict['pathname']
    filename = param_dict['filename']
    
    if VFID_or_objname == 'VFID':
        gal_name = filename.split('-')[1]  #isolates VFID from filename string
    elif VFID_or_objname == 'objname':
        gal_name = filename.split('-')[0]  #isolates objname from filename string
    else:
        gal_name = 'central_obj'
    
    gal_VFID = param_dict['galaxy_VFID']
    
    vf_path = homedir+param_dict['vf_catalog']
    vf = Table.read(vf_path)
    csv_savename = gal_name+param_dict['csv_savename']
    
    vid_name_cutout = pathname + gal_VFID+ param_dict['vid_name_cutout']
    vid_name_midi = pathname + param_dict['vid_name_midi']
    
    r25_multiple_x = float(param_dict['r25_multiple_x'])
    r25_multiple_y = float(param_dict['r25_multiple_y'])
    
    strips_per_beat = int(param_dict['strips_per_beat'])
    y_scale = float(param_dict['y_scale'])
    note_names = param_dict['note_names']   #calls which note_names I added to params, either note_names1, note_names2, etc.
    note_names = param_dict[note_names].split("-")   #converts note_names1, note_names2, etc. to list of notes
    vel_min = int(param_dict['vel_min'])
    vel_max = int(param_dict['vel_max'])
    bpm = int(param_dict['bpm'])
    midi_or_2D = param_dict['midi_or_2D']
    
    wav_savename = pathname+gal_VFID+param_dict['wav_savename']
    soundfont = homedir+param_dict['soundfont']

    #create strips
    #will need first four terms für später
    im, lower_bound_y, upper_bound_y, lower_bound_x, upper_bound_x, strips, VFID = define_strips(VFID_or_objname, pathname, gal_name, filename, vf, r25_multiple_x, r25_multiple_y)
    
    #write and save csv file
    create_csv_file(strips, pathname, csv_savename)
    
    #load csv file as pd table
    df = load_csv_file(csv_savename,pathname)
    
    #define y_data as the mean pixel value of each strip
    y_data = df['mean']
    
    #'compress the x-axis data to define the duration of each beat
    t_data = compress_time(y_data,strips_per_beat)
    
    #normalize y_data, then scale according to **0.5
    y_data, y_data_scaled = normalize_scale(y_data,y_scale)
    
    #generate a plot of the data, if desired
    #visualize_data(t_data,y_data)
    
    #create midi data list (notes to midi values)
    midi_data = note_to_midi(y_data_scaled,note_names)
    
    #determine volume of each datum
    vel_data = map_note_velocities(y_data_scaled,vel_min,vel_max)
    
    #save the midi file!
    midi_savename = pathname+VFID+'.mid'
    save_midi_file(t_data, midi_data, vel_data, bpm, midi_savename, VFID = VFID)
    
    #extract length of midi file...
    time = get_midi_length(pathname,midi_savename)
        
    #generate and save animation
    if midi_or_2D == 'midi':
        create_animation_midi(t_data,midi_data,time,vid_name_midi)
        print('.mp4 saved to',vid_name_midi)
    
    if midi_or_2D == '2D':
        create_animation_cutout(VFID_or_objname, pathname, filename, im, lower_bound_x, upper_bound_x, lower_bound_y, upper_bound_y, time, vid_name_cutout)
        print('.mp4 saved to',vid_name_cutout)
    
    if midi_or_2D == 'both':
        create_animation_midi(t_data,midi_data,time,vid_name_midi)
        print('.mp4 saved to',vid_name_midi)
        create_animation_cutout(VFID_or_objname, pathname, filename, im, lower_bound_x, upper_bound_x, lower_bound_y, upper_bound_y, time, vid_name_cutout)
        print('.mp4 saved to',vid_name_cutout)
    
    #lastly, generate and save the wavfile
    convert_to_wav(midi_savename,wav_savename,soundfont,gain=2)
    print('.wav saved to',wav_savename)
    
    print('fin.')


