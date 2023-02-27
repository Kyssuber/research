# research
- Local Macbook notebooks and .py programs.

## data_analysis
- a dedicated repository to organize code which outputs various relevant analyses of the model fit parameters. VERY unique to my directory paths, filenames, etc.

## notebooks
- somewhat unwieldy directory featuring various coding tests and tutorials. 'unwieldy' is not an embellishment - you may think of this section as a modest programmer's junk drawer.

## programs

### virtual-machine
- note that all can be run outside of ipython environment. 
- e.g., in terminal, python make_inputcomp.py -h 
- the aim is to write a series of routines that enable me to run galfit on vigo-filament galaxies using Rose Finn's virtual machine

##### params_inputcomp*txt -->
- NOT CURRENTLY IN USE.
- the * is a placeholder for the wavelength band and whether or not convolution is involved 
- w3 is currently the only iteration. cp and change name and parameters accordingly, as needed. 
- contains all necessary parameters to correctly run make_inputcomp.py 
- can edit as needed, especially the pathnames to certain files and directories. 
- note that "output" filenames are embedded in make_inputcomp.py itself. I do not recommend altering, since these names are propagated through the entire virtual-machine routine. 

##### make_inputcomp.py --> 
GOAL: 
- generate galfit input scripts for sample galaxies 

INPUT: 
- python make_inputcomp.py -param_file /path/to/paramfile/paramfile (no quotation marks) 
- paramfile should be params_inputcomp*.txt, where the * is replaced with the band 

OUTPUT: 
- individual input scripts for each central galaxy in sample 
- can directly feed these text files into GALFIT 
- directly edit params_inputcomp*.txt to tweak parameters as desired 

##### params_galfit_mosaics_output*.txt -->
- the * is a placeholder for the wavelength band and whether or not convolution is involved 
- w3 is currently the only iteration. cp and change name and parameters accordingly, as needed. 
- parameters file for run_galfit.py, output_params.py, and output_mosaics.py 
- (I merged the separate text files due to immense overlap and the relatively few number of params for each 

##### run_galfit.py -->
GOAL: 
- run galfit on all scripts generated using make_inputcomp.py 
OUTPUT: 
- .FITS output files containing mosaic and output parameters (via header information); output_params.py and output_mosaic.py will then attend to these clunky buggers. 
- the 'organize' parameter will, more often than not, be True. In this scenario, all output .fits will be funneled into the gal_output_path of your choice. 

##### output_params.py -->
GOAL:
- gather all GALFIT output parameters, as well as the galnames and 'central' designations, and compile into an Astropy table. 
- only includes galaxies with VFIDs. entries in 'dummy catalog' which do not have a VFID are not added to the table. 
- be sure directories contain all GALFIT output .fits files before running. 
- output will be multiple .fits tables, one for each central galaxy and one that combines parameters of all galaxies. 

##### output_mosaics.py -->
GOAL: 
- generate galfit output mosaic .png files for each virgo-filament galaxy 

### mucho-galfit
- comprises .py files to facilitate the running of galfit in parallel on the Siena HPC machine

##### make_input_file.py -->
- generates input .txt file for running galfit in parallel 
- each directory is labeled as the galaxy's rounded RA coordinate, so this input file contains one row per galaxy with the galaxy objname and corresponding round(RA). 

##### gather_output_params.py -->
GOAL:
- gather all GALFIT output parameters, as well as the galnames and 'central' designations, and compile into an Astropy table. 
- essentially a repurposing of /virtual-machine/output_params.py, but with dir names unique to our purposes. 

#### build_html_website*.py
- script in development that will write an html file in order to generate a locally-hosted website to 'easily' navigate galfit output mosaics and parameters for the Virgo subsample.

### make_plots.py
- in development, potentially indefinitely.

## sonification
- a novice side project which aims to convert visual data into scientifically meaningful audio. I ultimately hope to contribute even modestly to the mainstreaming of sonfication into the larger scientific community as a tool to enhance sensitivity to pattern detection. it is additionally a nifty outreach pursuit, if I ever become proficient enough to manipulate code to generate aesthetically appealing media for the public.

### galaxyfits
- the first (and presently only) addition to this larger blanket project. galaxyfits consists of code that takes the fits image of a 2D galaxy cutout and traces the length of the image with a 1-px wide red bar whose height is based on the central galaxy's (OPTICAL) diameter. For each step, this bar occupies a column of pixels, with a mean value that the code maps onto a set of user-defined midi notes. The string of these notes is then converted into a .wav file, plastered onto an animation of this bar sweeping the image and synchronized as well as my wittle program can manage, and saved as a .mp4. As this text is a miserable attempt at an explanation, I include a YouTube link below to an example galaxy (NGC2273A - see video description).
- the textfiles allow for convenient tweaking of certain parameters, such as the desired soundfont as well as the chords/midi notes. 
- NOTE: the 'optical diameter' aspect of the code currently relies on our collaboration's VF catalog as well as the galaxies' unique VFIDs in order to extract the radii from said catalog. I include the catalog along with an example fits file in order to test the code (I cannot supply the soundfont because it is a rather hefty bugger), but users will have to manipulate the source code a bit to align with their own purposes. (I find is comedic how I am assuming that any human aside from myself will desire to use this.)
- example: https://www.youtube.com/watch?v=ff8R6JEHp7c
> python sonification_galaxyfits.py -param_file params_sonification_galaxyfits.txt

> python animiation_galaxyfits.py -param_file params_animation_galaxyfits.txt




































