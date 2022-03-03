# cheby_checker tests go in this directory

[//]: # (TODO: 2 tests failing. From `mpcpp` being unavailable. cf. Dockerfile:59)

[//]: # (    MAY NOT BE ACTIVE BUGS ANY LONGER: The first three are from reboundx and may be fixable via https://github.com/dtamayo/reboundx/issues/26, or more probably https://github.com/dtamayo/reboundx/issues/39)

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

## The following test suites pass within a containerized environment built using the code in `../cheby_container`. They also pass in a local environment (cf. `../README.md`).

### `pytest test_convenience_Horizons.py`
 - `convenience_Horizons.py` only depends on EXTERNAL packages
 - passes within and without containerized environment 

### `pytest test_obs_pos.py`
 - `obs_pos.py` depends on INTERNAL `MPC_library.py`: eventually want to change dependency to `wis.py`
 - passes within and without containerized environment
   
[//]: # (TODO: no wis.py in repo)

### `pytest test_cheby_checker.py`
 - cheby_checker.py only depends on EXTERNAL packages

[//]: # (TODO: Run this suite first.)

### `pytest test_sql.py`
- `sql.py` has an INTERNAL dependence on the `cheby_checker.py` module (above)

### `pytest test_coco.py`
- `coco.py` has an INTERNAL dependence on the `cheby_checker.py` module (above)
- `coco.py` also depends on INTERNAL `MPC_library.py`: should shift the rotn & ecliptic angle variables to coco & Base respectively.

### `pytest test_nbody_parse.py`
- `nbody.py` has an INTERNAL dependence on the coco.py module (above)
- `nbody.py` also depends on INTERNAL `MPC_library.py`: should shift away from this if/when possible
- `nbody.py` also depends on flaky `REBOUNDX:EPHEM` library

### `pytest test_nbody_run.py`
- nbody.py has an INTERNAL dependence on the coco.py module (above)
- nbody.py also depends on INTERNAL MPC_library.py: should shift away from this if/when possible
- nbody.py also depends on the s REBOUNDX:EPHEM library
- There are tests in place which get the same cartesian coords to ~10km when comparing an orbfit result to a JPL result. This seems good enough.

[//]: # (TODO: test_run_mpcorb_A is not passing on my local setup, but is on the container.)

### `pytest test_orbit_cheby_locations.py`
- `orbit_cheby.py` has an INTERNAL dependence on the `nbody.py` module
- `sql.py` has an INTERNAL dependence on the `cheby_checker.py` module 
- The `orbit_cheby` module has a lot going on, so I am sub-dividing the tests: Here I am splitting out tests of the "_define_locations", "_take_triangular", & "_make_square" functions within `orbit_cheby`. 

### `pytest test_orbit_cheby_creation.py`
- orbit_cheby.py has an INTERNAL dependence on the nbody.py module
- sql.py has an INTERNAL dependence on the cheby_checker.py module 
- The orbit_cheby module has a lot going on, so I am sub-dividing the tests: Here I am splitting out tests of the instantiation of MSC & MCS_Loader objects within orbit_cheby, as well as the *from_coord_arrays* function(s) used to populate the objects with data at start up.
- Detailed tests of the *accuracy* of the chebyshevs that get created as a by-product of calling *from_coord_arrays* will primarily be performed elsewhere
- passes within containerized environment (not yet local)

[//]: # (TODO: 4 tests failing locally.. cf. TODOs in orbit_cheby.py. Likewise with test_orbit_cheby_horizons.py. These are passing on the container now.)

### `pytest test_orbit_cheby_horizons.py`
- orbit_cheby.py has an INTERNAL dependence on the nbody.py module
- sql.py has an INTERNAL dependence on the cheby_checker.py module 
- Here I am starting to develop tests that explicitly do integrations, from the EXACT same starting coords as Horizons
- This allows me to double-check the accuracy of the rebound-to-cheby conversions & the XYZ-to-RADEC conversions, ..
- passes within containerized environment (not yet local)

### `pytest test_nbody_covariance.py`
- nbody.py has an INTERNAL dependence on the coco.py module (above)
- nbody.py also depends on INTERNAL MPC_library.py: should shift away from this if/when possible
- nbody.py also depends on the s REBOUNDX:EPHEM library
- Here I am focusing on tests of the propagation of the cartesian covariance matrix by the run_mpcorb routine (which uses *_get_covariance_from_tangent_vectors* under-the-hood)
- passes within containerized environment

## The following suites were identified during a March 2022 review

The following suites are passing as-is: 
- `test_data_classes.py` 

The following tests pass with new changes:
- `test_malloc_reboundx.py` (after using a similar `integration_function(.)` call as in `test_nbody_run.py`)
- `test_nbody_NbodySim.py` (with a provision to skip `test_initialize_integration_function_A()` if there's not enough memory to do the simulation (Need > 58.1 GB).)
- `test_sifter_sqlite.py` passing after updating the tests to use the proper `SQLSifter` object. Also made several updates to `sql.py` to ensure the `SQLSifter` class routines make the appropriate references to `DB`.

## The following suites have issues I cannot easily resolve
### `test_convenience_functions.py`
- Not passing. The `nbody` library has changed from what it used to be in `archaic/nbody20220131`. On `nbody.py:385` there's a reference to the function that can ostensibly load an orbit file from text, not JSON, which is the only current function that runs. Perhaps this test isn't necessary anymore?

### `test_ephem.py`
- Skipping the one test for now since I'm not sure how to instantiate the Ephemeris object. 

### `test_nbody_ParseElements.py`
- Skipping 5 tests as the fixtures appear to be unavailable.

### `test_orbit_cheby.py`
- see todo in file:427.

### `test_sifter_base.py` & `test_sifter_tracklets.py`
- Need `mpcpp`?

### `test_sockets.py`
- As the tests mention, the latter two tests pass after running `sockets_server_starter_DO_NOT_DELETE.py` in the background.
- Not sure what's wrong with the first test yet.

### `tmp.py`
- Has a call to `integration_function(.)` out of a function, and the test appears to be looking for fixtures perhaps related to this call. Skipping for now.