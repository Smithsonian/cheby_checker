# -*- coding: utf-8 -*-
# /tests/test_nbody.py

"""
----------------------------------------------------------------------------
tests for cheby_checker/nbody
 - Here I focus on the functions that RUN the nbody integrations.
 - In the accompanying test_nbody_parse script I test the functions for PARSING the input data

Jan 2022
Matthew Payne

Prev Work:
Mike Alexandersen, Matthew Payne & Matthew Holman

This code simplified as of Jan 2022
Removing many tests of non-json input
 - The non-json input methods *may* still work, but for now I just want to ensure that the json inputs work
 - Old tests of the non-json input can be found in the tests/archiaic/ttestt_nbody* file(s)
 tests for mpc_nbody

----------------------------------------------------------------------------
"""

# import third-party packages
# -----------------------------------------------------------------------------
import sys
import os
import numpy as np
import pytest
from astropy.time import Time
import pytest
from filecmp import cmp
import getpass
import json
import time
import glob


# Import neighbouring packages
# -----------------------------------------------------------------------------
sys.path.append(os.environ['REBX_DIR'])
from examples.ephem_forces import ephem_forces

# The main nbody code we are trying to test
from cheby_checker import nbody

# old conversion library that may be useful for cross-comparison of various tests ...
from cheby_checker import MPC_library as mpc

import convenience_Horizons as Horizons

this_dir = os.path.abspath(os.path.dirname( __file__ ))
repo_dir = os.path.abspath(os.path.dirname( this_dir ))
data_dir = os.path.join(repo_dir, 'dev_data')
json_dir = os.path.join(data_dir, 'mpc_orb_jsons')
std_json_dir = os.path.join(json_dir, 'standard_mp') # Standard grav-only fits


# Utility functions to help with testing
# -----------------------------------------------------------------------------

def similar_xyzuvw(xyzv0, xyzv1, threshold_xyz=1e-13, threshold_v=1e-14): # 15 mm, 1.5 mm/day
    """
    Calculate the difference between two sets of cartesian coordinates.
    """
    if isinstance(xyzv0, list):
        xyzv0 = np.array(xyzv0)
    if isinstance(xyzv1, list):
        xyzv1 = np.array(xyzv1)
    error = xyzv0 - xyzv1
    good_tf = np.abs(error) < np.array([threshold_xyz] * 3 + [threshold_v] * 3)
    return np.all(good_tf), error


def is_nbody_output_good_enough(times, data, target='30102'):
    """
    Mike Alexandersen & Matt Payne
    Helper function for determining whether the saved output from an nbody
    integration is good enough.
    """
    # Check 20 timesteps
    # NB If len(times)<20, indicees will initially contain repeats, but we subsequently force uniqueness using "set"
    indicees = np.linspace(0, len(times) - 1, 20).astype(int)
    for i in set(indicees):
    
        # Get Horizons "state"" for that time
        horizons_xyzv = nice_Horizons(target, '500@0', times[i], 'smallbody')
        
        # The "state" from our integration
        mpc_xyzv = data[i, 0, :]
        
        # Check whether our integration agrees with Horizons ( within threshold )
        error, good_tf = compare_xyzv(horizons_xyzv,
                                        mpc_xyzv,
                                        1e-7, 1e-8) # MJP 2020-09-03 : Artificially increased thresholds to allow me to make progress while waiting for Holman to debug
                                        #5e-11, 2e-13)  # 7.5m, 30 mm/day
        if np.all(good_tf):
            pass # print('Awesome!')
        else:
            print(f'\n Discrepancy in *is_nbody_output_good_enough()* ...')
            print(f'Time, timestep: {times[i]:}, {i:}')
            print(f'Horizons : {["%18.15e" % _ for _ in horizons_xyzv]}')
            print(f'N-body   : {["%18.15e" % _ for _ in mpc_xyzv]}')
            print(f'Position off by [au]: {error[:3]:}')
            print(f'Velocity off by [au/day]: {error[3:6]:}')
        assert np.all(good_tf)


# Tests of NbodySim
# -----------------------------------------------------------------------------

def test_nbody_A():
    """
    Test instantiation of NbodySim object
    """
    # Instantiate
    N = nbody.NbodySim()

    # Check that the expected attributes exist
    assert \
        N.tstart             == None  and \
        N.tstop              == None  and \
        N.geocentric         == False  and \
        N.integration_epoch  == None  and \
        N.mpcorb_list        == []  and \
        N.input_n_particles  == None  and \
        N.save_output        == False  and \
        N.save_output_file   == None  and \
        N.verbose            == False  and \
        N.CHECK_EPOCHS       == True  and \
        N.unpacked_primary_provisional_designation == None and \
        N.helio_ecl_vec_EXISTS   == False  and \
        N.helio_ecl_vec          == None  and \
        N.helio_ecl_cov_EXISTS   == False  and \
        N.helio_ecl_cov          == None  and \
        N.bary_eq_vec_EXISTS     == False  and \
        N.bary_eq_vec            == None  and \
        N.bary_eq_cov_EXISTS     == False  and \
        N.bary_eq_cov            == None  and \
        N.non_grav_EXISTS        == False  and \
        N.non_grav_dict_list     == []  and \
        N.output_times       == None  and \
        N.output_states      == None  and \
        N.output_covar       == None
    

def test_integration_function_A():
    """
    If we put ANYTHING into the ephem_forces.integration_function,
    will it work or crash and burn?

    Most likely if there is a problem, it'll cause pytest to crash entirely,
    so might as well start with this.

    Note that *ephem_forces.integration_function* is NOT the main function that will
    be called by cheby_checker.nbody, but it IS called by
    *ephem_forces.production_integration_function_wrapper*, which IS called by
    cheby_checker.nbody

    Note that this does NOT perform any tests of integration accuracy:
    I am purely trying to get a simple call to execute

    """

    # Define the initial positions. We use the state for Asteroid (3666) Holman from DE441
    # 2458849.500000000 = A.D. 2020-Jan-01 00:00:00.0000 TDB [del_T=     69.183900 s]
    # X = 3.338875350265349E+00 Y =-9.176518267602161E-01 Z =-5.038590677470149E-01
    # VX= 2.805663315227095E-03 VY= 7.550408688437705E-03 VZ= 2.980028207454247E-03
    tstart = 2458849.5
    instates = np.array([[3.338875349745594E+00, -9.176518281675284E-01, -5.038590682977396E-01, 2.805663319000732E-03, 7.550408687780768E-03, 2.980028206579994E-03]])

    # Define the length of the integration (days)
    tstep = 20.0
    tend = tstart + 1.0


    # Define the quantities associated with the variational equaitons
    n_var = 6
    n_particles = 7
    geocentric = 0
    invar_part = np.zeros(6, dtype=int)
    invar = np.identity(6)

    # Run the integration
    times, states, var, var_ng, status = ephem_forces.integration_function(tstart, tend, tstep, geocentric, n_particles, instates, n_var, invar_part, invar)


    assert isinstance(states, np.ndarray)
    assert isinstance(times, np.ndarray)
    assert isinstance(var, np.ndarray)


def test_integration_function_B():
    """
    ... ephem_forces.integration_function,
    ...

    """

    # Define the variables that will be used in the query
    target  = '12345'
    centre  = '500@0'
    epochs  = ['2458850.0','2458880.0']
    id_type = 'smallbody'
    refplane='earth'

    # Call the *nice_Horizons* function to get the cartesian states at the first time
    # This is returning EQUATORIAL BERYCENTRIC coordinates
    horizons_zero = Horizons.nice_Horizons(target, centre, epochs[0], id_type, refplane=refplane)

    # Run the integration
    tstart = float(epochs[0])
    tend   = float(epochs[-1])
    tstep  = 0.01
    geocentric = 0
    n_particles = 1
    instates = np.array([horizons_zero])
    n_var = 6
    invar_part = np.zeros(6, dtype=int)
    invar = np.identity(6)
    times, states, var, var_ng, status = ephem_forces.integration_function(tstart,
                                                                            tend,
                                                                            tstep,
                                                                            geocentric,
                                                                            n_particles,
                                                                            instates,
                                                                            n_var,
                                                                            invar_part,
                                                                            invar)

    # Now call horizons again at some of the output times at which the integration_function produced output
    for n, t in enumerate(times):

        # Check horizons at the simulation-time output
        h = Horizons.nice_Horizons(target, centre, t, id_type, refplane=refplane)

        # Check similarity: threshold_xyz=1e-11, threshold_v=1e-12 : # 150 cm, 15 cm/day
        similar_bool , error = similar_xyzuvw(h, states[n][0], threshold_xyz=1e-11, threshold_v=1e-12)
        assert similar_bool
        

def test_production_integration_function_wrapper_A():
    """
    First test of *ephem_forces.production_integration_function_wrapper*
    Doing a 2-particle integration
    Testing the STRUCTURE of the returned arrays
    *NOT* testing the numerical accuracy

    """
    tstart  = 2456184.7
    tstop   = tstart + 600
    epoch   = tstart
    
    geocentric = 0
    n_particles = 2
    instates = np.array([
                        [-3.1, 2.7, 3.6, -0.006, -0.004, -0.002] ,
                        [-4.1, 3.7, 5.6, -0.004, -0.003, -0.002]
                        ])
        
    # Call the function
    outtimes, states, partial_derivatives_wrt_state, partial_derivatives_wrt_NG, return_value = ephem_forces.production_integration_function_wrapper(  tstart,
                                                                tstop,
                                                                epoch,
                                                                instates,
                                                                non_grav_dict_list = None,
                                                                geocentric = 0 ,
                                                                epsilon = 1e-8,
                                                                tstep_min = 0.02,
                                                                tstep_max = 32.)
    # check that we get the expected variable-types back
    assert isinstance(states, np.ndarray)
    assert isinstance(outtimes, np.ndarray)
    assert isinstance(partial_derivatives_wrt_state, np.ndarray)
    assert isinstance(return_value, tuple)

    # check that the SHAPES of the returned arrays are as expected ...
    # outtimes.shape                        ~ (161,)
    # states.shape                          ~ (161, 2, 6)
    # partial_derivatives_wrt_state.shape   ~ (161, 12, 6)
    assert  outtimes.shape[0] == states.shape[0] == partial_derivatives_wrt_state.shape[0]
    assert  states.ndim ==3 and \
            states.shape[1] == n_particles and \
            states.shape[2] == 6
    assert  partial_derivatives_wrt_state.ndim ==3 and \
            partial_derivatives_wrt_state.shape[1] == 6*n_particles and \
            partial_derivatives_wrt_state.shape[2] == 6
    
    
def test_production_integration_function_wrapper_B():
    """
    Doing a timing test of the production_integration_function_wrapper

    When running directly from a file in ephem_forcs, this test takes
    ~1 secs to run (20,000 day integration)

    During development there were issues in which the same query run from a
    different directory would take ~3mins (i.e. 100 times longer)

    I want to ensure that we don't have this problem
    """
    start_time = time.time()
    #print(__file__, ':test_production_integration_function_wrapper_B')
    
    # Define the inputs
    instates = np.array([[3.338875349745594E+00, -9.176518281675284E-01, -5.038590682977396E-01, 2.805663319000732E-03, 7.550408687780768E-03, 2.980028206579994E-03]])
    tstart, trange = 2458849.5, 20000
    epoch = tstart
    tend = tstart + trange
    
    # Run the integration
    times, states, var, var_ng, status = ephem_forces.production_integration_function_wrapper(tstart, tend, epoch, instates)

    # Check that the execution time is less than 5 seconds (when it gets bad it ballons to 3 minutes)
    end_time = time.time()
    allowed_time = 5 #seconds
    print('test_production_integration_function_wrapper_B: end_time-start_time',end_time-start_time)
    assert end_time-start_time < allowed_time


def test_production_integration_function_wrapper_C():
    """
    Another test of timing...
    In practice 10 particles seemed to take ~2.5 secs (end-start below)
    """
    start_time = time.time()
    #print(__file__, ':test_production_integration_function_wrapper_B')
    
    # Define the inputs
    instates = np.array([
        [3.338875349745594E+00, -9.176518281675284E-01, -5.038590682977396E-01, 2.805663319000732E-03, 7.550408687780768E-03, 2.980028206579994E-03],
        [3.388875349745594E+00, -9.176518281675284E-01, -5.038590682977396E-01, 2.805663319000732E-03, 7.540408687780768E-03, 2.980028206579994E-03],
        [3.438875349745594E+00, -9.176518281675284E-01, -5.038590682977396E-01, 2.805663319000732E-03, 7.530408687780768E-03, 2.980028206579994E-03],
        [3.488875349745594E+00, -9.176518281675284E-01, -5.038590682977396E-01, 2.805663319000732E-03, 7.520408687780768E-03, 2.980028206579994E-03],
        [3.538875349745594E+00, -9.176518281675284E-01, -5.038590682977396E-01, 2.805663319000732E-03, 7.510408687780768E-03, 2.980028206579994E-03],
        [3.588875349745594E+00, -9.176518281675284E-01, -5.038590682977396E-01, 2.805663319000732E-03, 7.550408687780768E-03, 2.980028206579994E-03],
        [3.638875349745594E+00, -9.176518281675284E-01, -5.038590682977396E-01, 2.805663319000732E-03, 7.540408687780768E-03, 2.980028206579994E-03],
        [3.688875349745594E+00, -9.176518281675284E-01, -5.038590682977396E-01, 2.805663319000732E-03, 7.530408687780768E-03, 2.980028206579994E-03],
        [3.738875349745594E+00, -9.176518281675284E-01, -5.038590682977396E-01, 2.805663319000732E-03, 7.520408687780768E-03, 2.980028206579994E-03],
        [3.788875349745594E+00, -9.176518281675284E-01, -5.038590682977396E-01, 2.805663319000732E-03, 7.510408687780768E-03, 2.980028206579994E-03],
    ])
    tstart, trange = 2458849.5, 20000
    epoch = tstart
    tend = tstart + trange
    
    # Run the integration
    times, states, var, var_ng, status = ephem_forces.production_integration_function_wrapper(tstart, tend, epoch, instates)

    # Check that the execution time is less than 15 seconds (when it gets bad it ballons to 3 minutes)
    end_time = time.time()
    allowed_time = 15 #seconds
    print('test_production_integration_function_wrapper_C: end_time-start_time',end_time-start_time)
    assert end_time-start_time < allowed_time


def test_production_integration_function_wrapper_D():
    """

    Testing the numerical ACCURACY of the results returned from a
    single run of the *ephem_forces.production_integration_function_wrapper*
    function.

    Note that many more tests of the accuracy of the reboundx integrator
    need to be performed, but they need to be done in the REBOUNDX
    package itself
    """
    
    # Define the variables that will be used in the query
    target  = '719' # Asteroid #123456
    centre  = '500@0'  # <<-- Barycentric
    epochs  = ['2458000.0','2458200.0']
    id_type = 'smallbody'
    refplane='earth' # <<--Equatorial

    # Call the *nice_Horizons* function to get the cartesian states at the first time
    # This is returning EQUATORIAL BERYCENTRIC coordinates
    horizons_zero = Horizons.nice_Horizons(target, centre, epochs[0], id_type, refplane=refplane)
    print('horizons_zero=',horizons_zero)
    
    # Call the production_integration_function_wrapper that we want to test
    # I believe that the integration is performed in EQUATORIAL BERYCENTRIC coordinates
    tstart = epoch = epochs[0]
    tstop  = epochs[-1]
    instates = np.array([ horizons_zero ])
    outtimes, states, partial_derivatives_wrt_state, partial_derivatives_wrt_NG, return_value = ephem_forces.production_integration_function_wrapper(   float(tstart),
                                                                float(tstop),
                                                                float(epoch),
                                                                instates,
                                                                non_grav_dict_list = None,
                                                                geocentric = 0 ,
                                                                epsilon = 1e-10,
                                                                tstep_min = 0.02,
                                                                tstep_max = 32.)
                                                                
    # Now call horizons again at some of the output times at which the *production_integration_function_wrapper()* produced output
    for n, t in enumerate(outtimes):
    
        # Check horizons at the simulation-time output
        h = Horizons.nice_Horizons(target, centre, t, id_type, refplane=refplane)

        # Check similarity: threshold_xyz=1e-11, threshold_v=1e-13 : # 15 m, 1.5 m/day
        similar_bool , error = similar_xyzuvw(h, states[n][0], threshold_xyz=1e-10, threshold_v=1e-11)
        assert similar_bool, f'test_production_integration_function_wrapper_D:n={n}, error={error}'
        

def test_run_integration_A():
    """
    Test the NbodySim.__run_integration() function
    This
        calls *production_integration_function_wrapper*
        then performs various functions to reshape the partial derivs,
        then calculates the cov-matrix at each timestep

    We want to check that the output has the expected type & shape
    """
    
    # Instantiate
    N = nbody.NbodySim()
    
    # Declare inputs
    tstart = 2459200
    tstop  = 2459300
    data_file = '2000SR210.json'
    mpcorb_list = [ os.path.join(std_json_dir , data_file) ]

    # Parse the inputs
    N._parse_inputs_run_mpcorb(tstart = tstart , tstop = tstop , mpcorb_list = mpcorb_list)
    
    # Convert to barycentric equatorial
    N.make_bary_equatorial()

    # Check some quantities BEFORE running the test function ...
    assert N.input_n_particles == 1
    assert N.output_times   == None
    assert N.output_states  == None
    assert N.output_covar   == None

    # ## ### #### ##### ######
    # Run the *__run_integration* function we want to test
    # ## ### #### ##### ######
    N.verbose = True
    N._run_integration()
    
    # Check the output has the expected type & shape
    assert isinstance(N.output_times, np.ndarray)
    assert isinstance(N.output_states, np.ndarray)
    assert isinstance(N.output_covar, np.ndarray)

    # output_times.shape   ~ (161,)
    # output_states.shape  ~ (161, 1, 6)
    # output_covar.shape   ~ (161, 1, 6, 6)
    n_particles = 1
    assert  N.output_times.shape[0] == N.output_states.shape[0] == N.output_covar.shape[0]
    assert  N.output_states.ndim ==3 and \
            N.output_states.shape[1] == n_particles and \
            N.output_states.shape[2] == 6
    assert  N.output_covar.ndim ==4 and \
            N.output_covar.shape[1] == n_particles and \
            N.output_covar.shape[2] == 6 and \
            N.output_covar.shape[3] == 6


@pytest.mark.skip(reason="Not Passing: similar_bool == False")
def test_run_mpcorb_A():
    """
    Test the overall *run_mpcorb* function

    Note that here we are testing that
    (a) the cartesian positions from the orbfit fit are "close enough" to the JPL fit (10's of km)
    (b) the cartesian positions remain close to the JPL fit throughout the integration by reboundx
    """

    # Instantiate
    N = nbody.NbodySim()
    
    # Declare inputs
    tstart = 2459200
    tstop  = 2459300
    data_file = '2000SR210.json'
    mpcorb_list = [ os.path.join(std_json_dir , data_file) ]

    # Run the *run_mpcorb* that we want to test
    N.run_mpcorb( tstart = tstart , tstop = tstop , mpcorb_list = mpcorb_list )


    # Define the variables that will be used in the query
    target  = '56986' # 2000SR210 == Asteroid #56986
    centre  = '500@0'  # <<-- Barycentric
    id_type = 'smallbody'
    refplane='earth' # <<--Equatorial


    # Now call horizons again at the output times at which the *production_integration_function_wrapper()* produced output
    for n, t in enumerate(N.output_times):
    
        # Check horizons at the simulation-time output
        h = Horizons.nice_Horizons(target, centre, t, id_type, refplane=refplane)

        # Check similarity: threshold_xyz=1e-6, threshold_v=1e-8
        # The expected differences are at the few-times 10^-7 level => few-times 10km
        #  -- This seems entirely reasonable for the differences in the results of orbit fits
        similar_bool , error = similar_xyzuvw(h, N.output_states[n][0], threshold_xyz=1e-6, threshold_v=1e-8)
        assert similar_bool


def test_run_mpcorb_B():
    """
    Test the overall *run_mpcorb* function on MULTIPLE INPUT FILES

    """

    # Instantiate
    N = nbody.NbodySim()
    
    # Declare inputs
    tstart = 2459200
    tstop  = 2459300
    mpc_orb_json_files = glob.glob(std_json_dir + '/*.json' )[:3]

    # Run the *run_mpcorb* that we want to test
    N.run_mpcorb( tstart = tstart , tstop = tstop , mpcorb_list = mpc_orb_json_files )

    # Add some tests !!!
        
    

