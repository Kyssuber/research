'''
Goal: 
-- create vertical bar animation that sweeps across and is in sync with overlayed sonified data
Requirements: 
-- .wav of audio file. software to .mid to .wav is available here at https://www.zamzar.com/convert/midi-to-wav/
-- parameter textfile, currently params_animation_galaxyfits.txt; update accordingly
Outputs:
-- .mp4 animation saved to Desktop, combined with sound
'''

import os
import ffmpeg
import sys
import numpy as np

def combine_save(video_pathname,audio_pathname,output_pathname):
    input_video = ffmpeg.input(video_pathname)
    input_audio = ffmpeg.input(audio_pathname)
    ffmpeg.concat(input_video, input_audio, v=1, a=1).output(output_pathname).run()

    
if __name__ == '__main__':
    
    if '-h' in sys.argv or '--help' in sys.argv:
        print("Usage: %s [-param_file (name of parameter file, no single or double quotation marks)]")
        sys.exit(1)
    
    if '-param_file' in sys.argv:
        p = sys.argv.index('-param_file')
        param_file = str(sys.argv[p+1])
     
    homedir = os.getenv("HOME")
       
    #create dictionary with keywords and values, from parameter .txt file

    param_dict = {}
    with open(homedir+'/sonification/galaxyfits/'+param_file) as f:
        for line in f:
            try:
                key = line.split()[0]
                val = line.split()[1]
                param_dict[key] = val
            except:
                continue    
    
    #now...extract parameters and assign to relevantly-named variables
    
    gal_VFID = param_dict['galaxy_VFID']
    
    path_to_wav = param_dict['path_to_wav']
    wavname = gal_VFID+param_dict['wavname']
    
    pathname = param_dict['pathname']
    video_name_midi = param_dict['vid_name_midi']
    video_name_2D = gal_VFID+param_dict['vid_name_cutout']
    
    outputdir = param_dict['savedir']
    outputvid_midi = gal_VFID+param_dict['outputvid_midi']
    outputvid_cutout = gal_VFID+param_dict['outputvid_cutout']
    
    audio_pathname = homedir + path_to_wav + wavname
    video_pathname_midi = homedir + pathname + video_name_midi
    video_pathname_cutout = homedir + pathname + video_name_2D
    output_pathname_midi = homedir + outputdir + outputvid_midi
    output_pathname_cutout = homedir + outputdir + outputvid_cutout
    
    midi_or_2D = param_dict['midi_or_2D']
    
    if midi_or_2D == 'midi':
        combine_save(video_pathname_midi,audio_pathname,output_pathname_midi)
        print('.mp4 saved at',output_pathname_midi)
    
    if midi_or_2D == '2D':
        combine_save(video_pathname_cutout,audio_pathname,output_pathname_cutout)
        print('.mp4 saved at',output_pathname_cutout)
        
    if midi_or_2D == 'both':
        combine_save(video_pathname_midi,audio_pathname,output_pathname_midi)
        print('.mp4 saved at',output_pathname_midi)
        combine_save(video_pathname_cutout,audio_pathname,output_pathname_cutout)
        print('.mp4 saved at',output_pathname_cutout)
    
    
    print('fin.')
    
    
    
    
    
    

