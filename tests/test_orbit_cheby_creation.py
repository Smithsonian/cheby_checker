# -*- coding: utf-8 -*-

'''
    --------------------------------------------------------------
    tests of orbit_cheby's base class
    ../tests/test_orbit_cheby_creation
    
    Just focusing on tests of the MSC/MSC_Loader instantiation and
    population-with-data
    
    Tests of other detailed functionality & accuracy will be done
    elsewhere

    Feb 2022
    Matt Payne
    
    --------------------------------------------------------------
'''


# Import third-party packages
# --------------------------------------------------------------
import sys
import os
import numpy as np
import pytest
import json
import glob
import math
from datetime import datetime
from astropy.time import Time

# Import MPC packages
# -----------------------------------------------------------------------------
from mpc_orb.parse import MPCORB

# Import neighboring packages
# --------------------------------------------------------------
from cheby_checker import nbody
from cheby_checker import orbit_cheby
from cheby_checker import cheby_checker
from cheby_checker import obs_pos
from cheby_checker import coco

this_dir = os.path.abspath(os.path.dirname( __file__ ))
repo_dir = os.path.abspath(os.path.dirname( this_dir ))
data_dir = os.path.join(repo_dir, 'dev_data')
json_dir = os.path.join(data_dir, 'mpc_orb_jsons')
std_json_dir = os.path.join(json_dir, 'standard_mp')


# Constants & Test Data
# -----------------------------------------------------------------------------
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'dev_data')

FLAT_FILES = [  os.path.join(DATA_DIR, '2022AA_demo.txt') ,
                os.path.join(DATA_DIR, 'simulation_states.dat')]
orbfit_filenames = [os.path.join(DATA_DIR, file) for file in ['30101.eq0_horizons', '30102.eq0_horizons']]
mpc_orb_json_files = glob.glob(std_json_dir + '/*.json' )


# Convenience data / functions to aid testing
# --------------------------------------------------------------
def similar_xyzuvw(xyzv0, xyzv1, threshold_xyz=1e-13, threshold_v=1e-14): # 15 mm, 1.5 mm/day
    '''
    Calculate the difference between two sets of cartesian coordinates.
    '''
    if isinstance(xyzv0, list):
        xyzv0 = np.array(xyzv0)
    if isinstance(xyzv1, list):
        xyzv1 = np.array(xyzv1)
    error = xyzv0 - xyzv1
    if len(error) == 3:
        good_tf = np.abs(error) < np.array([threshold_xyz] * 3 )
    if len(error) == 6:
        good_tf = np.abs(error) < np.array([threshold_xyz] * 3 + [threshold_v] * 3)

    return np.all(good_tf), error
    


def convenience_call_to_nbody_run_mpcorb( json_filepath_or_list, tstart=2459200 , tstop=2459295 ):
    ''' 
        A convenience function to return a simulation object
        Proper testing of mpc_nbody is done elsewhere (test_nbody_run.py)
    '''

    # Instantiate
    N = nbody.NbodySim()
    #N.verbose = True
    
    # Now run the integrator
    # NB: tstart = 2459200 = First date in Sector # 600 (check using B.map_JD_to_sector_number_and_sector_start_JD(2459200, B.standard_MJDmin )
    #     tstop = 2459295 = :Last day in Sector # 602   (check using B.map_JD_to_sector_number_and_sector_start_JD(2459295, B.standard_MJDmin )
    # This is useful to know for later testing ...
    mpcorb_list = [ json_filepath_or_list ] if isinstance(json_filepath_or_list,str) else json_filepath_or_list
    N.run_mpcorb( tstart = tstart , tstop = tstop , mpcorb_list = mpcorb_list )

    # Quick tests ...
    for attrib in ["output_times", "output_states"]:
        assert attrib in N.__dict__
        assert isinstance( N.__dict__[attrib] , np.ndarray )

    return N
    
def convenience_call_to_get_MSCs_from_file_via_nbody_run( json_filepath_or_list, tstart=2459200 , tstop=2459295):
    '''
        A convenience function to return a list of MSCs
        Starts from mpc_orb_json file(s)
        Runs NBody integration via nbody.NbodySim.run_mpcorb
        Uses MSC_Loader to do all of the work to declare and populate a list of MSC objects
        Tests of MSC_Loader are performed below (e.g. *test_create_loader* & *test_loader_from_nbodysim*)
        
    '''

    # Instantiate Nbody Object, Run NBody Simulation, and return populated Object
    N = convenience_call_to_nbody_run_mpcorb(json_filepath_or_list , tstart=tstart , tstop=tstop)
    
    # Use the MSC_Loader to do all of the work to declare and populate a list of MSC objects
    MSCs = orbit_cheby.MSC_Loader(NbodySim = N).MSCs
    
    return MSCs
    





# Actual tests : MSC Object [ tests of MSC_Loader object are below]
# -----------------------------------------------------------------

def test_convenience_call_to_nbody_run_mpcorb():
    ''' just calling the functions to ensure they do something ...'''
    # Call with one file
    convenience_call_to_nbody_run_mpcorb( mpc_orb_json_files[0])
    # Call with two files
    convenience_call_to_nbody_run_mpcorb( mpc_orb_json_files[:2])

    assert True


def test_create_empty_MSC():
    ''' Create an empty MSC class object'''
    
    # Initialize the multi_sector_cheby
    result = orbit_cheby.MSC()
    
    # Check the result is as expected
    assert isinstance( result , orbit_cheby.MSC )
    assert 'sector_coeffs' in result.__dict__
    assert isinstance( result.sector_coeffs , dict )




def test_from_coord_arrays_A(  ):
    '''
        Use MSC functionality to load from numpy arrays
        
        Significant reliance on *generate_cheb_for_sector()* under-the-hood
         - We test the ACCURACY of the returned coefficients, so we are testing
           both *from_coord_arrays*  &  *from_coord_arrays*
         
        Here we test whether various quanities are populated as expected
        when we pass in data from a SINGLE PARTICLE integration
        
        
        
    '''
    # The file that we will use as the source of the data
    mpc_orb_json_filepath = mpc_orb_json_files[0]
    
    # Do a local read of the json (using MPCORB) as it's useful to be able to read the name of the object
    M = MPCORB(mpc_orb_json_filepath)
    primary_unpacked_provisional_designation = M.designation_data["unpacked_primary_provisional_designation"]
    
    # Instantiate Nbody Object, Run NBody Simulation, and return populated Object
    N = convenience_call_to_nbody_run_mpcorb(mpc_orb_json_filepath)
    
    # Instantiate MSC object
    M=orbit_cheby.MSC()

    # ## ### #### #####
    # Run the *from_coord_arrays* function that we want to test ...
    # NB1: This attempts to load from arrays ...
    # NB2: Need to do array slicing that will be taken care of within MSC_Loader in practical operation
    # ## ### #### #####
    M.from_coord_arrays(primary_unpacked_provisional_designation, N.output_times , N.output_states[:,0,:] )

    # check that the expected attributes have been set
    for attr in ["TDB_init", "TDB_final", "sector_init", "sector_final", "sector_coeffs"]:
        assert attr in M.__dict__

    # check that the "supported_sector_numbers" are as expected
    # NB: convenience_call_to_nbody_run_mpcorb uses
    #     tstart = 2459200 = First date in Sector # 600 (check using B.map_JD_to_sector_number_and_sector_start_JD(2459200, B.standard_MJDmin )
    #     tstop = 2459295 = :Last day in Sector # 602   (check using B.map_JD_to_sector_number_and_sector_start_JD(2459295, B.standard_MJDmin )
    # Hence expected_supported_sector_numbers = [600,601,602]
    B = cheby_checker.Base()
    assert M.sector_init  == B.map_JD_to_sector_number_and_sector_start_JD(N.tstart, B.standard_MJDmin )[0]
    assert M.sector_final == B.map_JD_to_sector_number_and_sector_start_JD(N.tstop, B.standard_MJDmin )[0]
    assert M.TDB_init     == B.map_JD_to_sector_number_and_sector_start_JD(N.tstart, B.standard_MJDmin )[1]
    assert M.TDB_final    == B.map_JD_to_sector_number_and_sector_start_JD(N.tstop, B.standard_MJDmin )[1] + B.sector_length_days - B.epsilon

    # Check the settings used for the cheby-fitting
    assert M.minorder == 5 , 'The default is expected to be *5*, albeit there is little reasoning for why ...'

    # Check that the sector_coeff dict is of the correct shape
    for sector_number , cheb_coeffs in M.sector_coeffs.items():

        # The # of coefficients for each sector should be > minorder
        #      ( if order = 7 , N_coefficients = 8)
        # NB> v.shape = ( N_coefficients , N_coords  ) , and N_coords = 6
        assert cheb_coeffs.shape[0] > M.minorder
        assert cheb_coeffs.shape[1] , N.output_states[:,0,:].shape[1]
        
        
    # Check the accuracy of the sector_coeff dict coefficients ...
    for n, t in enumerate(N.output_times): # Loop over the times from the NBody integration output ...
        
        # Input state ...
        inputState = N.output_states[:,0,:][n]
        
        # Evaluate cheby at given time ...
        # NB: returns array of shape (N_times, N_components)
        chebyEval = M.evaluate_components([t])[0]

        # Check that the cheby-evaluation is "close-enough" to the input state
        # NB: threshold_xyz=1e-11 => ~150cm == 1.5m
        similarBool, err_arr = similar_xyzuvw(inputState,chebyEval , threshold_xyz=1e-11, threshold_v=1e-11)
        assert similarBool , f'inputState={inputState} , chebyEval={chebyEval} , err_arr={err_arr}'



def test_from_coord_arrays_B(  ):
    '''
        Use MSC functionality to load from numpy arrays
        
        Significant reliance on *generate_cheb_for_sector()* under-the-hood
         - We test the ACCURACY of the returned coefficients, so we are testing
           both *from_coord_arrays*  &  *from_coord_arrays*
         
        Here we test whether various quanities are populated as expected
        when we pass in data from MANY SINGLE-PARTICLE integrations
        
        Follows tthe same pattern as test_from_coord_arrays_A, above
        
    '''

    # Loop over the files that we will use as the source of the data
    for mpc_orb_json_filepath in mpc_orb_json_files:
        print('\n','--'*22)
        print(mpc_orb_json_filepath)
        # Instantiate Nbody Object, Run NBody Simulation, and return populated Object
        N = convenience_call_to_nbody_run_mpcorb(mpc_orb_json_filepath)
        #print('N.output_times=',N.output_times)
        
        # Do a local read of the json (using MPCORB) as it's useful to be able to read the name of the object
        M = MPCORB(mpc_orb_json_filepath)
        primary_unpacked_provisional_designation = M.designation_data["unpacked_primary_provisional_designation"]

        # Instantiate MSC object
        M=orbit_cheby.MSC()

        # ## ### #### #####
        # Run the *from_coord_arrays* function that we want to test ...
        # NB1: This attempts to load from arrays ...
        # NB2: Need to do array slicing that will be taken care of within MSC_Loader in practical operation
        # ## ### #### #####
        M.from_coord_arrays(primary_unpacked_provisional_designation, N.output_times , N.output_states[:,0,:] )

        # check that the expected attributes have been set
        for attr in ["TDB_init", "TDB_final", "sector_init", "sector_final", "sector_coeffs"]:
            assert attr in M.__dict__

        # check that the "supported_sector_numbers" are as expected
        # NB: convenience_call_to_nbody_run_mpcorb uses
        #     tstart = 2459200 = First date in Sector # 600 (check using B.map_JD_to_sector_number_and_sector_start_JD(2459200, B.standard_MJDmin )
        #     tstop = 2459295 = :Last day in Sector # 602   (check using B.map_JD_to_sector_number_and_sector_start_JD(2459295, B.standard_MJDmin )
        # Hence expected_supported_sector_numbers = [600,601,602]
        B = cheby_checker.Base()
        assert M.sector_init  == B.map_JD_to_sector_number_and_sector_start_JD(N.tstart, B.standard_MJDmin )[0]
        assert M.sector_final == B.map_JD_to_sector_number_and_sector_start_JD(N.tstop, B.standard_MJDmin )[0]
        assert M.TDB_init     == B.map_JD_to_sector_number_and_sector_start_JD(N.tstart, B.standard_MJDmin )[1]
        assert M.TDB_final    == B.map_JD_to_sector_number_and_sector_start_JD(N.tstop, B.standard_MJDmin )[1] + B.sector_length_days - B.epsilon

        # Check the settings used for the cheby-fitting
        #print('B.standard_MJDmin =', B.standard_MJDmin )
        assert M.minorder == 5 , 'The default is expected to be *5*, albeit there is little reasoning for why ...'

        # Check that the sector_coeff dict is of the correct shape
        for sector_number , cheb_coeffs in M.sector_coeffs.items():

            # The # of coefficients for each sector should be > minorder
            #      ( if order = 7 , N_coefficients = 8)
            # NB> v.shape = ( N_coefficients , N_coords  ) , and N_coords = 6
            assert cheb_coeffs.shape[0] > M.minorder
            assert cheb_coeffs.shape[1] , N.output_states[:,0,:].shape[1]
            
            
        # Check the accuracy of the sector_coeff dict coefficients ...
        for n, t in enumerate(N.output_times): # Loop over the times from the NBody integration output ...

            # Input state ...
            inputState = N.output_states[:,0,:][n]
            
            # Evaluate cheby at given time ...
            # NB: returns array of shape (N_times, N_components)
            chebyEval = M.evaluate_components([t])[0]
            
            # Check that the cheby-evaluation is "close-enough" to the input state
            # NB1: threshold_xyz=1e-11 => ~150cm == 1.5m
            # NB2: **** 2022-02-07 *** : For reasons related to lack of control over the number of output points, I have changed
            #      min number of points-per-secttor to 5 (instead of 7) and worsened the allowed accuracy on this tast
            #      from 1e-11 to 1e-10 (1,500cm). I would like to go back to the more stringent consttraint. But in order
            #      to go back, I/we need better control over the output points returned from the production_integration_function_wrapper
            similarBool, err_arr = similar_xyzuvw(inputState,chebyEval , threshold_xyz=1e-10, threshold_v=1e-11)
            assert similarBool , f'inputState={inputState} , chebyEval={chebyEval} , err_arr={err_arr}'






# Actual tests : MSC_Loader
# --------------------------------------------------------------
def test_create_loader():
    ''' Create an empty MSC_Loader class object'''
    
    # Initialize the multi_sector_cheby
    result = orbit_cheby.MSC_Loader()
    
    # Check the result is as expected
    assert isinstance( result , orbit_cheby.MSC_Loader )
    assert 'MSCs' in result.__dict__
    assert isinstance( result.MSCs , list )


def test_loader_from_arrays():
    ''' 
        Use loader to create an array of MSCs, starting from numpy arrays
        Pretty High-Level function, lots of dependencies ... 
    '''

    # The file that we will use as the source of the data
    mpc_orb_json_filepath = mpc_orb_json_files[0]
    
    # Do a local read of the json (using MPCORB) as it's useful to be able to read the name of the object
    M = MPCORB(mpc_orb_json_filepath)
    primary_unpacked_provisional_designation = M.designation_data["unpacked_primary_provisional_designation"]
    
    # Instantiate Nbody Object, Run NBody Simulation, and return populated Object
    N = convenience_call_to_nbody_run_mpcorb(mpc_orb_json_filepath)

    
    # Attempt to instantiate MSC via MSC_Loader
    MSCs = orbit_cheby.MSC_Loader(
                                  primary_unpacked_provisional_designations = primary_unpacked_provisional_designation,
                                  times_TDB = N.output_times,
                                  statearray = N.output_states).MSCs
    assert isinstance(MSCs , list)
    assert len(MSCs) == 1
    
    M = MSCs[0]
    assert isinstance(M, orbit_cheby.MSC)

    return M


def test_loader_from_nbodysim():
    ''' 
        Use loader to create an array of MSCs, starting from mpc_nbody.nbodysim object
        Pretty High-Level function, lots of dependencies ... 
    '''

    # The file that we will use as the source of the data
    mpc_orb_json_filepath = mpc_orb_json_files[0]

    # Instantiate Nbody Object, Run NBody Simulation, and return populated Object
    N = convenience_call_to_nbody_run_mpcorb(mpc_orb_json_filepath)
    
    # Use the MSC_Loader to do all of the work to declare and populate a list of MSC objects
    MSCs = orbit_cheby.MSC_Loader(NbodySim = N).MSCs
    assert isinstance(MSCs , list)
    assert len(MSCs) == 1
    
    M = MSCs[0]
    assert isinstance(M, orbit_cheby.MSC)
    
    return MSCs

