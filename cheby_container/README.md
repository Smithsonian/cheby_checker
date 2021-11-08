# Installing prereqs for cheby_checker

MJP: 2021-11-08

## Docker approach

The scripts in this directory provide a means to set up a
reasonably minimal container that can get you up-and-running 
with a development version of the cheby_checker code. 

The important / difficult part was/is ensuring that the rebound & 
reboundx (holman-branch) codes are installed and compiled, along
with the required ephemeris files. 

Note that I am deliberately *NOT* cloning in the cheby_checker 
repo from github, but rather am binding in a local volume (see 
the build_container.py file) to allow interactive development
of an assumed local version. 

### To build & execute container ...
>>> python3 build_container.py

This builds and runs the image and deposit you on the command line of the running container. 

### To test / run some of the cheby_checker code ...
>>> cd /cheby_checker/tests
>>> pytest test_sql.py

This should successfully run a bunch of tess of the sql module. 

MJP 
 
