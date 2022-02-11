# cheby_checker tests go in this directory

## 2022 : MJP: The following are WIP, attempting to get tests to pass... 



### test_orbit_cheby_horizons.py :
- orbit_cheby.py has an INTERNAL dependance on the nbody.py module
- sql.py has an INTERNAL dependance on the cheby_checker.py module 
- Here I am starting to develo`p tests that explicitly do integrations, from the EXACT same starting coords as Horizons
- This allows me to double-check the accuracy of the rebound-to-cheby conversions, the XYZ-to-RADEC conversions, ..
- I have *NOT* yet been able to run these in a container, due to compilation problems with rebound/reboundx ... 

### test_orbit_cheby.py :
- orbit_cheby.py has an INTERNAL dependance on the nbody.py module
- sql.py has an INTERNAL dependance on the cheby_checker.py module 
- a number of tests have been developed that establish working & accurate code for multiple functions around input STATES
- few/no tests have yet been done to incorporate COVARIANCE MATRIX input / fitting / etc




## 2022 : MJP: The following have been checked to pass tests within a containerized environment initialized/built using the code in "cheby_container"

### test_convenience_Horizons.py
 - convenience_Horizons.py only depends on EXTERNAL packages
 - pytest test_convenience_Horizons.py
 - passes within containerized environment 

### test_obs_pos.py
 - obs_pos.py depends on INTERNAL MPC_library.py: eventually want to change dependancy to wis.py
 - pytest test_obs_pos.py
 - passes within containerized environment

### test_cheby_checker.py
 - cheby_checker.py only depends on EXTERNAL packages
 - pytest test_cheby_checker.py
 - passes within containerized environment

### test_sql.py :
- sql.py has an INTERNAL dependance on the cheby_checker.py module (above)
- pytest test_sql.py
- passes within containerized environment

### test_coco.py :
- coco.py has an INTERNAL dependance on the cheby_checker.py module (above)
- coco.py also depends on INTERNAL MPC_library.py: should shift the rotn & ecliptic angle variables to coco & Base respectively.
- pytest test_coco.py
- passes within containerized environment

### test_nbody_parse.py :
- nbody.py has an INTERNAL dependance on the coco.py module (above)
- nbody.py also depends on INTERNAL MPC_library.py: should shift away from this if/when possible
- nbody.py also depends on flaky REBOUNDX:EPHEM library
- pytest test_nbody_parse.py
- passes within containerized environment

### test_nbody_run.py :
- nbody.py has an INTERNAL dependance on the coco.py module (above)
- nbody.py also depends on INTERNAL MPC_library.py: should shift away from this if/when possible
- nbody.py also depends on the s REBOUNDX:EPHEM library
- There are tests in place which get the same cartesian coords to ~10km when comparing an orbfit result to a JPL result. This seems good enough.

### test_orbit_cheby_locations.py :
- orbit_cheby.py has an INTERNAL dependance on the nbody.py module
- sql.py has an INTERNAL dependance on the cheby_checker.py module 
- The orbit_cheby module has a lot going on, so I am sub-dividing the tests: Here I am splitting out tests of the "_define_locations", "_take_triangular", & "_make_square" functions within orbit_cheby. 
- pytest test_orbit_cheby_locations.py
- passes within containerized environment

### test_orbit_cheby_creation.py :
- orbit_cheby.py has an INTERNAL dependance on the nbody.py module
- sql.py has an INTERNAL dependance on the cheby_checker.py module 
- The orbit_cheby module has a lot going on, so I am sub-dividing the tests: Here I am splitting out tests of the instantiation of MSC & MCS_Loader objects within orbit_cheby, as well as the *from_coord_arrays* function(s) used to populate the objects with data at start up.
 - Detailed tests of the *accuracy* of the chebyshevs that get created as a by-product of calling *from_coord_arrays* will promarily be performed elsewhere
- pytest test_orbit_cheby_creation.py
- passes within containerized environment


