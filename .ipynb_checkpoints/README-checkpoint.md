# research
- Local Macbook notebooks and .py programs.

## data_analysis
- a dedicated repository to organize code which outputs various relevant analyses of the model fit parameters. VERY unique to my directory paths, filenames, etc.

## notebooks
- somewhat unwieldy directory featuring various coding tests and tutorials. 'unwieldy' is not an embellishment - you may think of this section as a modest programmer's junk drawer.

## programs

### virtual-machine
- currently relatively vacant. all 'virtual-machine' code is found in the Master branch.

##### output_mosaics.py -->
GOAL: 
- generate galfit output mosaic .png files for each virgo-filament galaxy 


### generating_html_locally
- note that all can be run outside of ipython environment. 
- e.g., in terminal, python make_inputcomp.py -h 
- the aim is to write a series of routines that enable me to run galfit on vigo-filament galaxies using Rose Finn's virtual machine

#### build_html_website.py
- script that writes an html file in order to generate a locally-hosted website to 'easily' navigate galfit output mosaics and parameters for the Virgo subsample. Instructions are physically taped to a desk monitor that is out of commission. 

#### params_build_html_website.txt
- list of script parameters to tailor the build_html_website.py outputs to the user's dreams and desires.


### mucho-galfit
- comprises .py files to facilitate the running of galfit in parallel on the Siena HPC machine

##### make_input_file.py -->
- generates input .txt file for running galfit in parallel 
- each directory is labeled as the galaxy's rounded RA coordinate, so this input file contains one row per galaxy with the galaxy objname and corresponding round(RA). 

##### gather_output_params.py -->
GOAL:
- gather all GALFIT output parameters, as well as the galnames and 'central' designations, and compile into an Astropy table. 
- essentially a repurposing of /virtual-machine/output_params.py, but with dir names unique to our purposes. 

## sonification
- a novice side project which aims to convert visual data into scientifically meaningful audio. I ultimately hope to contribute even modestly to the mainstreaming of sonfication into the larger scientific community as a tool to enhance sensitivity to pattern detection. it is additionally a nifty outreach pursuit, if I ever become proficient enough to manipulate code to generate aesthetically appealing media for the public.

### galaxyfits
- the first (and presently only) addition to this larger blanket project. galaxyfits consists of code that takes the fits image of a 2D galaxy cutout and traces the length of the image with a 1-px wide red bar whose height is based on the central galaxy's (OPTICAL) diameter. For each step, this bar occupies a column of pixels, with a mean value that the code maps onto a set of user-defined midi notes. The string of these notes is then converted into a .wav file, plastered onto an animation of this bar sweeping the image and synchronized as well as my wittle program can manage, and saved as a .mp4. As this text is a miserable attempt at an explanation, I include a YouTube link below to an example galaxy (NGC2273A - see video description).
- the textfiles allow for convenient tweaking of certain parameters, such as the desired soundfont as well as the chords/midi notes. 
- NOTE: the 'optical diameter' aspect of the code currently relies on our collaboration's VF catalog as well as the galaxies' unique VFIDs in order to extract the radii from said catalog. I include the catalog along with an example fits file in order to test the code (I cannot supply the soundfont because it is a rather hefty bugger), but users will have to manipulate the source code a bit to align with their own purposes. (I find is comedic how I am assuming that any human aside from myself will desire to use this.)
- example: https://www.youtube.com/watch?v=ff8R6JEHp7c
> python sonification_galaxyfits.py -h

> python animiation_galaxyfits.py -h

- For an in-progress GUI I am in the midst of developing in my heaps of recreation time, please have a glimpse at my 'GUI Odyssey' GitHub repository. Thank you, and good night.




































