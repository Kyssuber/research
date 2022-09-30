# research
- Local Macbook notebooks and .py programs.

## notebooks
- somewhat unwieldy directory featuring various coding tests and tutorials; largely for my own benefit.

## programs

### virtual-machine
- note that all can be run outside of ipython environment. 
- e.g., in terminal, python make_inputcomp.py -h 
- the aim is to write a series of routines that enable me to run galfit on vigo-filament galaxies using Rose Finn's virtual machine

##### params_inputcomp*txt -->
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

### make_plots.py
- in development, potentially indefinitely.






































