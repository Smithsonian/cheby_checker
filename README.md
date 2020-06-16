# cheby_checker
Python code related to orbital ephemeris calculations (for solar system objects) using fast chebysheb polynomial representations. 

Functionalities exist/will-exist to:
 - Perform n-body integrations (using a custom REBOUND wrapper) based on inital orbital states provided from an ORBFIT orbital fit to observations;   
 - Convert the output from n-body integrations (in the form of periodic output of posn, vel, covar, etc) into chebyshev polynomial representations;
 - Conduct ephemeris calculations; 
 - Check the residuals between observations and expected positions from previous orbit extrapolations (pChecker-v2);
 - Calculate the objects which are likely to be within the Field-of-View of a telescope (MPChecker-v2);
 - Attribute observations (of unknown designation) to known objects (CheckID-v2); 

# to-do : 2020-06-16
 - Draw together disparate repos (nbody, mpchecker, orbit_cheby) into this repo
 - Compile a list of extent and missing functionality
 - Compile a list of extent and missing tests
 - Compile a list of extent and missing demos (notebooks)

