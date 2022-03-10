# cheby_checker
Python code related to orbital ephemeris calculations (for solar system objects) using fast Chebyshev polynomial representations. 

## Installation / prerequisites
 - Git LFS should be instantiated in your clone of this repo:
   ```bash
   # sudo apt install git-lfs # if needed 
   git lfs install
   git lfs pull
   # sanity check
   head -n 1 ./dev_data/de440.bsp
   ```
   As long as this output doesn't start with "version", you have the file and not the pointer from Git LFS.
 - A containerized approach has been set up to take care of
   - Installing all of the REBOUND/REBOUNDx NBody stuff
   - Installing all other python dependencies 
   - See `./cheby_container/README.md` for more details.
 - For a local installation, follow the steps in the Dockerfile from ln. 25 down. Also run `pip install -e .` to install this package itself.

## Functionalities exist/will-exist to:
 - Perform n-body integrations (using a custom REBOUND wrapper) based on inital orbital states provided from an ORBFIT orbital fit to observations;   
 - Convert the output from n-body integrations (in the form of periodic output of posn, vel, covar, etc) into chebyshev polynomial representations;
 - Conduct ephemeris calculations; 
 - Check the residuals between observations and expected positions from previous orbit extrapolations (pChecker-v2);
 - Calculate the objects which are likely to be within the Field-of-View of a telescope (MPChecker-v2);
 - Attribute observations (of unknown designation) to known objects (CheckID-v2); 
 - Identify archival observations (ITF or designated) of new objects (s9m-v2);

## Architecture of Data Flow
 - Checking_MPChecker_PCheck_CheckIDX_s9m_sifter.png
 ![Github](https://github.com/matthewjohnpayne/cheby_checker/blob/master/Checking_MPChecker_PCheck_CheckIDX_s9m_sifter.png) 


## to-do : 2021-11-08
 - See README in cheby_checker/cheby_checker code repo
 - See Google Doc in https://docs.google.com/document/d/1SsKRwc8nC_oxPtODlQOwGAzCq42Y4AvWFe-Fu6Fk9bI/edit#
