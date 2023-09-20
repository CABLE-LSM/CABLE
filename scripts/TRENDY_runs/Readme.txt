This folder contains the following subfolders:

v11: scripts associated with TRENDY v11 runs (2022). This setup uses the mpi CABLE version as in previous TRENDY runs.

v12: scripts associated with TRENDY v12 runs (2023). This setup differs from v11 and older versions as it runs CABLE in serial mode running multiple instances of the model in parallel.
     This folder contains the following scripts:

     - run_TRENDY.sh: This is the main run script in which basic settings (e.g. domain, number of landmasks, model paths, etc.) are set.
                      This script parses arguments through to other scripts located in this folder and runs them.
     - split_landmask.R: R script that creates landmasks for the individual runs as specified in the run_TRENDY.sh script. It works for 
                         both global and regional domains.
     - run_cable.sh: Script that runs all steps required for a complete CABLE run (including spinup). This script needs to be checked 
                     carefully before starting CABLE runs. It contains the run sequence and steps to run (e.g. number of spinup iterations),
                     as well as the PBS settings. This script also requires 'aux/run_cable-pop_lib.sh'.
     - merge_outputs.sh: in the serial setup (TRENDY v12 onwards), CABLE is run multiple times separately and independently. The outputs will
                         be located in folders named run1,run2,...,runx. The purpose of this script is to merge the output files spatially.
                         This is done with python script 'merge_to_output2d.py' written by Matthias Cuntz. For larger files this can be time consuming.
                         Therefore, one can specify the steps one would like to merge spatially in the main script run_TRENDY.sh ('mergesteps' variable).
                         For example, it makes sense to not merge the spinup files for most applications. Note that this script waits until 
                         all runs are finished and only starts if all runs are finished successfully.
     - cleanup.sh: cleans up the folder structure. Copies the exectuable and all restart, log, landmask files into one folder and deletes the rest.
     
     For debugging purposes, the individual steps described above can be switched on or off in the main script. 

     - postprocessing/postprocess_TRENDY.sh: main script to postprocess raw CABLE outputs to formats required by the TRENDY-MIP. This includes regridding,
                                             time stamps, unit conversion, attributes, etc. This script calls create_zonal_nbp_tables.r, which creates 
                                             overview tables of NBP values in different latitudinal zones. It also has the option to run ILAMB on some
                                             key variables, the functionality of which is to be improved in future versions.

     one more note: the scripts under aux also exist under the main scripts directory but may have been modified to work explicitly for gadi/TRENDY. 
                    Might need some harmonising in the future.
     

For questions please contact Juergen Knauer (J.Knauer@westernsydney.edu.au)
                         