# cheby_checker
Python code related to orbital ephemeris calculations (for solar system objects) using fast chebysheb polynomial representations. 

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

## Installation / prerequisites
In order to get the NBODY part of the code working, one needs to have REBOUND/REBOUNDx installed. 
The proceedure pasted below was required to get the code working on machines operated by MA & MJP.
Obviously this is non-optimal, but one would hope it could be automated (via setup.py) in the future. 
 - cd /home/mikea/Github
 - git clone https://github.com/hannorein/rebound.git
 - cd rebound
 - make clean; make
 - export REB_DIR=/home/mikea/Github/rebound
 - cd ../

 - git clone -b holman https://github.com/matthewholman/reboundx.git

 - cd reboundx/examples/central_force/
 - perl -pi -e 's/OPENGL=1/OPENGL=0/g' Makefile
 - make clean; make
 - export REBX_DIR=/Users/matthewjohnpayne/Envs/reboundx



## to-do : 2021-06-08
 - See README in cheby_checker/cheby_checker code repo
 - See Google Doc in ... 
