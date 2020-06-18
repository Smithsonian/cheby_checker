# cheby_checker
Python code related to orbital ephemeris calculations (for solar system objects) using fast chebysheb polynomial representations. 

Functionalities exist/will-exist to:
 - Perform n-body integrations (using a custom REBOUND wrapper) based on inital orbital states provided from an ORBFIT orbital fit to observations;   
 - Convert the output from n-body integrations (in the form of periodic output of posn, vel, covar, etc) into chebyshev polynomial representations;
 - Conduct ephemeris calculations; 
 - Check the residuals between observations and expected positions from previous orbit extrapolations (pChecker-v2);
 - Calculate the objects which are likely to be within the Field-of-View of a telescope (MPChecker-v2);
 - Attribute observations (of unknown designation) to known objects (CheckID-v2); 

# installation / prerequisites
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

 - The following step is necessary until holman merges pull request from MA
 - vi /Users/matthewjohnpayne/Envs/reboundx/examples/ephem_forces/ephem_forces.py
   # This file is in : ../../../reboundx/examples/ephem_forces/ephem_forces.py
   # Assume that the library file will be in the same directory : ../../../reboundx/examples/ephem_forces/libreboundx.so
   rebx_location = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'libreboundx.so' )

 - cd reboundx/examples/central_force/
 - perl -pi -e 's/OPENGL=1/OPENGL=0/g' Makefile
 - make clean; make
 - export REBX_DIR=/Users/matthewjohnpayne/Envs/reboundx



# to-do : 2020-06-16
 - Get reboundx working (can't import ephem_forces from elsewhere) 
 - Draw together disparate repos (nbody, mpchecker, orbit_cheby) into this repo
 - Compile a list of extent and missing functionality
 - Compile a list of extent and missing tests
 - Compile a list of extent and missing demos (notebooks)

