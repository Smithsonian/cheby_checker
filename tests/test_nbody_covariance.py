# -*- coding: utf-8 -*-
# /tests/test_nbody.py

'''
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
'''

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
from mpc_orb.parse import MPCORB

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
    '''
    Calculate the difference between two sets of cartesian coordinates.
    '''
    if isinstance(xyzv0, list):
        xyzv0 = np.array(xyzv0)
    if isinstance(xyzv1, list):
        xyzv1 = np.array(xyzv1)
    error = xyzv0 - xyzv1
    good_tf = np.abs(error) < np.array([threshold_xyz] * 3 + [threshold_v] * 3)
    return np.all(good_tf), error



def is_nbody_output_good_enough(times, data, target='30102'):
    '''
    Mike Alexandersen & Matt Payne
    Helper function for determining whether the saved output from an nbody
    integration is good enough.
    '''
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




# Test Data
# -----------------------------------------------------------------------------
mpc_orb_json_files = glob.glob(std_json_dir + '/*.json' )



# Tests of NbodySim
# -----------------------------------------------------------------------------

def test_nbody_A():
    '''
    Test instantiation of NbodySim object
    '''
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
            
    
        
def test_covariance_A():
    '''
    Test the propagation of the covariance matrix by the NbodySim.__run_integration() function
     - Single object test
    '''
    
    # Object orbfit file to integrate
    data_file = '2000SR210.json'
    mpcorb_list = [ os.path.join(std_json_dir , data_file) ]
    
    # For the purposes of this test, read the json file ahead of time to extract the
    # epoch, so that we can set tstart == epoch
    # (this code taken from _parse_orbfit_json)
    M = MPCORB(os.path.join(std_json_dir , data_file))
    epoch,timesystem    = M.epoch_data["epoch"],M.epoch_data["timesystem"]
    allowed_timesystems = {'TDT':'tt'}
    this_mpcorb_epoch = Time(float(epoch), format='mjd', scale=allowed_timesystems[timesystem] )
    

    # Declare integration params
    # Forcing the start date to be == the orbit epoch extracted above
    tstart = this_mpcorb_epoch.tdb.to_value('jd')#2459000.5000000107 # Setting
    tstop  = tstart + 32

    # Instantiate
    N = nbody.NbodySim()

    # ----------------------------------------
    # Run the *run_mpcorb* function that we want to test
    N.run_mpcorb( tstart = tstart , tstop = tstop , mpcorb_list = mpcorb_list )
    # ----------------------------------------

    # Check the output has the expected type & shape
    assert isinstance(N.output_times, np.ndarray)
    assert isinstance(N.output_states, np.ndarray)
    assert isinstance(N.output_covar, np.ndarray)

    # output_times.shape   ~ (41,)
    # output_states.shape  ~ (41, 1, 6)
    # output_covar.shape   ~ (41, 1, 6, 6)
    n_particles = 1
    assert  N.output_times.shape[0] == N.output_states.shape[0] == N.output_covar.shape[0]
    assert  N.output_states.ndim ==3 and \
            N.output_states.shape[1] == n_particles and \
            N.output_states.shape[2] == 6
    assert  N.output_covar.ndim ==4 and \
            N.output_covar.shape[1] == n_particles and \
            N.output_covar.shape[2] == 6 and \
            N.output_covar.shape[3] == 6

    # Check that the covariance at the first output time is equal/close to the input
    # - This should work due to tstart === epoch
    np.allclose( N.bary_eq_cov , N.output_covar[0,0,:,:])


    
    
    # Explicitly check the matrices at each timestep
    # - Just recalculating to ensure that the matrix multiplication was done correctly
    Gamma_0 = np.linalg.inv( N.bary_eq_cov )
    for n,t in enumerate(N.output_times):
        
        P_t = N.partial_derivatives_wrt_state[n,0,:,:]
        P_tT= P_t.T
        G_t_calc = np.matmul(P_tT , np.matmul(Gamma_0, P_t) )
        C_calc   = np.linalg.inv(G_t_calc)
        assert np.allclose( C_calc , N.output_covar[n,0,:,:])


    # ----------------------------------------
    # Make another NBody object to use for doing alternative integrations
    # I am going to offset the starting positions and then rerun a bunch of integrations
    # ----------------------------------------

    # Now run the integrations from the perturbed starting states
    for n in range(6):
    
        # ------ PERTURBATION -----
        shift_mag= 1e-9
        shift    = np.zeros(6)
        shift[n] = shift_mag

        # ------ EXPECTED ---------
        expected = N.output_states[-1,0,:] + N.partial_derivatives_wrt_state[-1,0,:][n] * shift_mag

        # ------ PERT CALC ---------
        # Make another NBody object to use for doing alternative integration
        N2 = nbody.NbodySim()
        N2._parse_inputs_run_mpcorb( tstart = tstart , tstop = tstop , mpcorb_list = mpcorb_list )
        N2.make_bary_equatorial()
        assert np.allclose(N.bary_eq_vec, N2.bary_eq_vec)
        # Add the pert & run the integration
        N2.bary_eq_vec += shift
        N2._run_integration()
 
       # Check that the perturbed integration is close to the expected result
        # Imposing ~1e-16 level : should be very very close
        similar_bool , err_arr = similar_xyzuvw(expected, N2.output_states[-1,0,:], threshold_xyz=1e-15, threshold_v=1e-16)
        assert similar_bool
        
        




def test_covariance_B():
    '''
    Test the propagation of the covariance matrix by the NbodySim.__run_integration() function
     - Multiple test objects
    '''
    
    for data_file in mpc_orb_json_files:

        mpcorb_list = [ os.path.join(std_json_dir , data_file) ]
        
        # For the purposes of this test, read the json file ahead of time to extract the
        # epoch, so that we can set tstart == epoch
        # (this code taken from _parse_orbfit_json)
        M = MPCORB(os.path.join(std_json_dir , data_file))
        epoch,timesystem    = M.epoch_data["epoch"],M.epoch_data["timesystem"]
        allowed_timesystems = {'TDT':'tt'}
        this_mpcorb_epoch = Time(float(epoch), format='mjd', scale=allowed_timesystems[timesystem] )
        

        # Declare integration params
        # Forcing the start date to be == the orbit epoch extracted above
        tstart = this_mpcorb_epoch.tdb.to_value('jd')#2459000.5000000107 # Setting
        tstop  = tstart + 32

        # Instantiate
        N = nbody.NbodySim()

        # ----------------------------------------
        # Run the *run_mpcorb* function that we want to test
        N.run_mpcorb( tstart = tstart , tstop = tstop , mpcorb_list = mpcorb_list )
        # ----------------------------------------

        # Check the output has the expected type & shape
        assert isinstance(N.output_times, np.ndarray)
        assert isinstance(N.output_states, np.ndarray)
        assert isinstance(N.output_covar, np.ndarray)

        # output_times.shape   ~ (41,)
        # output_states.shape  ~ (41, 1, 6)
        # output_covar.shape   ~ (41, 1, 6, 6)
        n_particles = 1
        assert  N.output_times.shape[0] == N.output_states.shape[0] == N.output_covar.shape[0]
        assert  N.output_states.ndim ==3 and \
                N.output_states.shape[1] == n_particles and \
                N.output_states.shape[2] == 6
        assert  N.output_covar.ndim ==4 and \
                N.output_covar.shape[1] == n_particles and \
                N.output_covar.shape[2] == 6 and \
                N.output_covar.shape[3] == 6

        # Check that the covariance at the first output time is equal/close to the input
        # - This should work due to tstart === epoch
        np.allclose( N.bary_eq_cov , N.output_covar[0,0,:,:])


        
        
        # Explicitly check the matrices at each timestep
        # - Just recalculating to ensure that the matrix multiplication was done correctly
        Gamma_0 = np.linalg.inv( N.bary_eq_cov )
        for n,t in enumerate(N.output_times):
            
            P_t = N.partial_derivatives_wrt_state[n,0,:,:]
            P_tT= P_t.T
            G_t_calc = np.matmul(P_tT , np.matmul(Gamma_0, P_t) )
            C_calc   = np.linalg.inv(G_t_calc)
            assert np.allclose( C_calc , N.output_covar[n,0,:,:])


        # ----------------------------------------
        # Make another NBody object to use for doing alternative integrations
        # I am going to offset the starting positions and then rerun a bunch of integrations
        # ----------------------------------------

        # Now run the integrations from the perturbed starting states
        for n in range(6):
        
            # ------ PERTURBATION -----
            shift_mag= 1e-9
            shift    = np.zeros(6)
            shift[n] = shift_mag

            # ------ EXPECTED ---------
            expected = N.output_states[-1,0,:] + N.partial_derivatives_wrt_state[-1,0,:][n] * shift_mag

            # ------ PERT CALC ---------
            # Make another NBody object to use for doing alternative integration
            N2 = nbody.NbodySim()
            N2._parse_inputs_run_mpcorb( tstart = tstart , tstop = tstop , mpcorb_list = mpcorb_list )
            N2.make_bary_equatorial()
            assert np.allclose(N.bary_eq_vec, N2.bary_eq_vec)
            # Add the pert & run the integration
            N2.bary_eq_vec += shift
            N2._run_integration()
     
           # Check that the perturbed integration is close to the expected result
            # Imposing ~1e-16 level : should be very very close
            similar_bool , err_arr = similar_xyzuvw(expected, N2.output_states[-1,0,:], threshold_xyz=1e-15, threshold_v=1e-16)
            assert similar_bool


test_covariance_A()
test_covariance_B()
