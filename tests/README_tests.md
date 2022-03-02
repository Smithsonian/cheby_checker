# cheby_checker tests go in this directory

[//]: # (TODO: Use Pytest skip/skipf decorators to skip WIP Tests. https://docs.pytest.org/en/latest/how-to/skipping.html#skipping-test-functions)
[//]: # (TODO: Get tests running in Docker/Local: get dev data, dependent packages mounted/installed)

## 2022 : MJP: The following are WIP, attempting to get tests to pass...

### test_orbit_cheby.py :
- orbit_cheby.py has an INTERNAL dependence on the nbody.py module
- orbit_cheby.py has an INTERNAL dependence on the cheby_checker.py module 
- Not sure whether any of the tests in here are still required. 
- Many tests have been established in test_orbit_cheby_locations / test_orbit_cheby_creation / test_orbit_cheby_horizons / test_orbit_cheby_covariance
- Definitely need to write tests of the propagation of RA, Dec covariances, but they might best belong in test_orbit_cheby_covariance

### test_precalc.py :
- orbit_precalc.py has an INTERNAL dependence on the nbody.py module
- orbit_precalc.py has an INTERNAL dependence on the sql.py module 
- orbit_precalc.py has an INTERNAL dependence on the orbit_cheby.py module 
- orbit_precalc.py has an INTERNAL dependence on the obs_pos.py module 
- Have started developing the required tests ... 

## 2022 : MJP: The following have been checked to pass tests within a containerized environment initialized/built using the code in `../cheby_container`

### `pytest test_convenience_Horizons.py`
 - `convenience_Horizons.py` only depends on EXTERNAL packages
 - passes within and without containerized environment 

### test_obs_pos.py
 - obs_pos.py depends on INTERNAL MPC_library.py: eventually want to change dependency to wis.py
 - pytest test_obs_pos.py
 - passes within containerized environment

### test_cheby_checker.py
 - cheby_checker.py only depends on EXTERNAL packages
 - pytest test_cheby_checker.py
 - passes within containerized environment

### `pytest test_sql.py`
- `sql.py` has an INTERNAL dependence on the `cheby_checker.py` module (above)
- passes within and without containerized environment

### test_coco.py :
- coco.py has an INTERNAL dependence on the cheby_checker.py module (above)
- coco.py also depends on INTERNAL MPC_library.py: should shift the rotn & ecliptic angle variables to coco & Base respectively.
- pytest test_coco.py
- passes within containerized environment

### test_nbody_parse.py :
- nbody.py has an INTERNAL dependence on the coco.py module (above)
- nbody.py also depends on INTERNAL MPC_library.py: should shift away from this if/when possible
- nbody.py also depends on flaky REBOUNDX:EPHEM library
- pytest test_nbody_parse.py
- passes within containerized environment

### test_nbody_run.py :
- nbody.py has an INTERNAL dependence on the coco.py module (above)
- nbody.py also depends on INTERNAL MPC_library.py: should shift away from this if/when possible
- nbody.py also depends on the s REBOUNDX:EPHEM library
- There are tests in place which get the same cartesian coords to ~10km when comparing an orbfit result to a JPL result. This seems good enough.

### test_orbit_cheby_locations.py :
- orbit_cheby.py has an INTERNAL dependence on the nbody.py module
- sql.py has an INTERNAL dependence on the cheby_checker.py module 
- The orbit_cheby module has a lot going on, so I am sub-dividing the tests: Here I am splitting out tests of the "_define_locations", "_take_triangular", & "_make_square" functions within orbit_cheby. 
- pytest test_orbit_cheby_locations.py
- passes within containerized environment

### test_orbit_cheby_creation.py :
- orbit_cheby.py has an INTERNAL dependence on the nbody.py module
- sql.py has an INTERNAL dependence on the cheby_checker.py module 
- The orbit_cheby module has a lot going on, so I am sub-dividing the tests: Here I am splitting out tests of the instantiation of MSC & MCS_Loader objects within orbit_cheby, as well as the *from_coord_arrays* function(s) used to populate the objects with data at start up.
 - Detailed tests of the *accuracy* of the chebyshevs that get created as a by-product of calling *from_coord_arrays* will primarily be performed elsewhere
- pytest test_orbit_cheby_creation.py
- passes within containerized environment

### test_orbit_cheby_horizons.py :
- orbit_cheby.py has an INTERNAL dependence on the nbody.py module
- sql.py has an INTERNAL dependence on the cheby_checker.py module 
- Here I am starting to develop tests that explicitly do integrations, from the EXACT same starting coords as Horizons
- This allows me to double-check the accuracy of the rebound-to-cheby conversions & the XYZ-to-RADEC conversions, ..
- pytest test_orbit_cheby_horizons.py
- passes within containerized environment


### test_nbody_covariance.py :
- nbody.py has an INTERNAL dependence on the coco.py module (above)
- nbody.py also depends on INTERNAL MPC_library.py: should shift away from this if/when possible
- nbody.py also depends on the s REBOUNDX:EPHEM library
- Here I am focusing on tests of the propagation of the cartesian covariance matrix by the run_mpcorb routine (which uses *_get_covariance_from_tangent_vectors* under-the-hood)
- pytest test_nbody_covariance.py
- passes within containerized environment
