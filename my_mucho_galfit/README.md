#INSTRUCTIONS FOR build_html_website_MG.py

##USAGE (iPython/Jupyter):

- Running this program automatically initiates the HomePage class (hp)
    
- if the -test arg is False, then replace single_galpage with hp in order to generate all images and such:
           
Markup : * hp.html_setup() --> create .html homepages for all galaxies in the VF subsample; also initiates GalPage class(single_galpage) for every galaxy in a loop, creating all relevant folders and files.
           
-- hp.create_LS_figures()
           
-- hp.create_galfit_mosaics(psf_index)--> (w3, nopsf), 1 (w3, psf), 2 (r, nopsf), 3 (r, psf)
           
-- hp.create_mask_mosaics()
   
- If the -test arg is True, then the user is wanting to test the GalPage class (this is a stategic idea, since this class contains the bulk of the functions required for the script to run successfully). In this case, the script automatically defines a single galpage class (single_galpage; uses index 0 for a random test galaxy), with which the user can test the following methods:
            
-- single_galpage.compile_LS_cutouts()
            
-- single_galpage.create_LS_mosaics()
            
-- single_galpage.create_model_mosaics_names()
            
-- single_galpage.create_model_mosaics(psf_index) --> 0 (w3, nopsf), 1 (w3, psf), 2 (r, nopsf), 3 (r, psf)
            
-- single_galpage.create_mask_mosaics()
            
-- single_galpage.tabulate_parameters()
            
-- single_galpage.WRITETHEGALPAGE()
            
 - **The current functionality is such that testing will have to occur in an iPython or Jupyter/notebook environment. ApoloGEEZ.**

 - **BEWARE OF MEMORY LEAKS** when running hp.create_galfit_mosaics(). I recommend taking care before executing that the user uses one psf_index at a time, then quit() iPython, then repeat. The additional (psf_index) arg is the best failsafe I could manage. I am a humble README.md. 

##USAGE (Terminal python):

*use -h or --help for the list of arguments.*

   - python build_html_website_MG.py -h
   
- -param_file path/to/paramfile (name of parameter file, no single or double quotation marks)
   
- -test (indicates that user is wanting to test the GalPage class functions)
   
- -mask (instructs script to create mask mosaics)
   
- -LS (creates LS cutouts, mosaics)
   
- -galfit (creates galfit output mosaics for all psf_indices specified in param_file)
   
- -html (creates .html scripts for galpages and homepage)
