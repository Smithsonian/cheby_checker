# -*- coding: utf-8 -*-
# sifter/tests/test_orbit_cheby

'''
    --------------------------------------------------------------
    tests of orbit_cheby's base class
    
    Jan 2020
    Matt Payne
    
    --------------------------------------------------------------
'''


# Import third-party packages
# --------------------------------------------------------------
import sys
import os
import numpy as np
import pytest

# Import neighboring packages
# --------------------------------------------------------------
sys.path.append(os.path.dirname(os.path.dirname(
                                                os.path.realpath(__file__))))
from cheby_checker import orbit_cheby
from cheby_checker import nbody_reader
from cheby_checker import mpc_nbody



# Constants & Test Data
# -----------------------------------------------------------------------------
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'dev_data')
FLAT_FILES = [os.path.join(DATA_DIR, '2022AA_demo.txt') , os.path.join(DATA_DIR, 'simulation_states.dat')]
orbfit_filenames = [os.path.join(DATA_DIR, file) for file in ['30101.eq0_horizons', '30102.eq0_horizons']]


# Convenience data / functions to aid testing
# --------------------------------------------------------------


@pytest.mark.parametrize( ('orbfit_file' ) , orbfit_filenames  )
def test_convenience_call_to_mpc_nbody_nbodysim(orbfit_file):
    ''' 
        A convenience function to return a simulation object
        Proper testing of mpc_nbody is done elsewhere (test_run_nbody.py)
    '''

    # First, let's initiate the class with an input file:
    Sim = mpc_nbody.NbodySim(orbfit_file, 'eq')

    # Now run the integrator, by calling the object.
    Sim(tstep=20, trange=1000)
    
    # Quick tests ...
    for attrib in ["output_times", "output_vectors"]:
        assert attrib in Sim.__dict__

    return Sim




# Actual tests : MSC Object
# --------------------------------------------------------------

def test_create_empty_MSC():
    ''' Create an empty MSC class object'''
    
    # Initialize the multi_sector_cheby
    result = orbit_cheby.MSC()
    
    # Check the result is as expected
    assert isinstance( result , orbit_cheby.MSC )
    assert 'sector_coeffs' in result.__dict__
    assert isinstance( result.sector_coeffs , dict )


@pytest.mark.parametrize( ('flat_file' , 'expected_sectors') ,
                         [
                            (FLAT_FILES[0] , (0,624) ),
                            (FLAT_FILES[1] , (503,522) )
                          ]
                         )
def test_from_coord_arrays( flat_file , expected_sectors ):
    '''
        Use MSC functionality to load from numpy arrays
        Significant reliance on *generate_cheb_for_sector()*
    '''
    
    # Get the data from supplied test files
    name, times, states  = nbody_reader.parse_nbody_txt( flat_file )

    # Instantiate
    M=orbit_cheby.MSC()

    # attempt to load from arrays ...
    primary_unpacked_provisional_designation = name
    M.from_coord_arrays(primary_unpacked_provisional_designation,times , states)

    # check that the expected attributes have been set
    for attr in ["TDB_init", "TDB_final", "sector_init", "sector_final", "sector_coeffs"]:
        assert attr in M.__dict__

    # check that the expected sector numbers are in the sector_coeffs dict
    assert isinstance( M.sector_coeffs , dict )
    assert M.sector_init == expected_sectors[0]
    assert M.sector_final == expected_sectors[1]
    assert expected_sectors[0] in M.sector_coeffs
    assert expected_sectors[1] in M.sector_coeffs
    assert expected_sectors[1] + 1 not in M.sector_coeffs
    assert expected_sectors[0] - 1 not in M.sector_coeffs


@pytest.mark.parametrize( ('flat_file' ) , FLAT_FILES  )
def test_generate_cheb_for_sector(flat_file):
    ''' 
        Use MSC function to generate the chebyshev coefficients for a given (32 day) sector
        ***Need to test the hell out of this, as it underlies everything***
        *** *** Current tests insufficient !!!!!!!!!!!!!!!!!!!!!!!! *** ***
    '''

    # Get the data from file
    name, times, states  = nbody_reader.parse_nbody_txt( flat_file )

    # Instantiate
    M=orbit_cheby.MSC()

    # Loop over some different starting times:
    # - want to test different parts of the data
    delta = 10
    for n in range(delta):
    
        # indicees to use to select from data (ensure that its in allowed range)
        ind = np.min( [n * int(len(times)/(delta+1)), len(times)-32] ) + np.arange(32)
        
        # make times relative
        # (we require 0<t<M.sector_length_days in M.generate_cheb_for_sector)
        rel_times = times[ind] - times[ind][0]
        assert np.all(rel_times < M.sector_length_days )
    
        # Call *generate_cheb_for_sector* using sub-arrays
        cheb_coeffs, maxErr = M.generate_cheb_for_sector(rel_times , states[ind])
        
        # Test the shape of the returned coefficients
        # We expect shape = (No, Nc), where Nc is the number of components being fitted,
        # and No is 1+order (and note that "order" can iterate upwards)
        assert M.minorder == 7 , 'The default is expected to be *7*, albeit there is little reasoning for why ...'
        assert cheb_coeffs.shape[0] >= M.minorder
        assert cheb_coeffs.shape[1] == states.shape[1]
    
        # Evaluate returned coeffs ...
        quickEval = np.polynomial.chebyshev.chebval(rel_times , cheb_coeffs).T
        print(flat_file)
        print(cheb_coeffs.shape, maxErr)
        # Test the accuracy of the returned coefficients
        assert M.maxerr <= 1e-8
        assert np.max( np.abs(quickEval - states[ind]) ) < M.maxerr

@pytest.mark.parametrize( ('flat_file' ) ,FLAT_FILES  )
def test_evaluate_components(flat_file):
    """
        Use MSC function to evaluate cheby components at set of supplied times
        ***Need to test the hell out of this, as it underlies everything***
    """

    # Get the data from file
    name, times, states  = nbody_reader.parse_nbody_txt( flat_file )
    
    # Instantiate MSC via MSC_Loader
    MSCs = orbit_cheby.MSC_Loader(
                                  primary_unpacked_provisional_designations = name,
                                  times_TDB = times,
                                  statearray = states).MSCs

    # Test the MSCs & M returned (repeating *test_loader_from_arrays* ) ...
    assert isinstance(MSCs , list)
    assert len(MSCs) == 1
    M = MSCs[0]
    assert isinstance(M, orbit_cheby.MSC)
    
    # Indicees for data supported by this MSC:
    init, final = M.get_valid_range_of_dates()
    ind = np.where( (times >= init) & (times <= final ) )

    # Make the call to the low-level function, *evaluate_components*
    evaluatedComponents      = M.evaluate_components( times[ind] )
    evaluatedComponentsSlice = M.evaluate_components( times[ind] , component_slice_spec=slice(0,3))

    # Make a call to the higher level function, *generate_XYZ*
    XYZs = M.generate_XYZ( times[ind] )

    # Check that the shape of the returned components is as expected ...
    assert states.T.shape == evaluatedComponents.shape
    assert evaluatedComponentsSlice.shape == XYZs.shape

    # Check that the values in evaluatedComponentsSlice == those in XYZs
    assert np.all( evaluatedComponentsSlice == XYZs )

    # Check that the values in evaluatedComponentsSlice are all close to the input state
    np.all(np.abs(evaluatedComponents[0] - states[:, 0][ind]) < M.maxerr)




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


@pytest.mark.parametrize( ('flat_file' ) ,FLAT_FILES  )
def test_loader_from_arrays(flat_file):
    ''' 
        Use loader to create an array of MSCs, starting from numpy arrays
        Pretty High-Level function, lots of dependencies ... 
    '''

    # Get the data from file
    name, times, states  = nbody_reader.parse_nbody_txt( flat_file )
    
    # Attempt to instantiate MSC via MSC_Loader
    MSCs = orbit_cheby.MSC_Loader(
                                  primary_unpacked_provisional_designations = name,
                                  times_TDB = times,
                                  statearray = states).MSCs
    assert isinstance(MSCs , list)
    assert len(MSCs) == 1
    
    M = MSCs[0]
    assert isinstance(M, orbit_cheby.MSC)

    return M


@pytest.mark.parametrize( ('orbfit_file' ) , orbfit_filenames  )
def test_loader_from_nbodysim(orbfit_file):
    ''' 
        Use loader to create an array of MSCs, starting from mpc_nbody.nbodysim object
        Pretty High-Level function, lots of dependencies ... 
    '''

    # Get the data from Sim ...
    Sim = test_convenience_call_to_mpc_nbody_nbodysim(orbfit_file)
    
    # Use the MSC_Loader to do all of the work to declare and populate a list of MSC objects
    MSCs = orbit_cheby.MSC_Loader(NbodySim = Sim).MSCs
    assert isinstance(MSCs , list)
    assert len(MSCs) == 1
    
    M = MSCs[0]
    assert isinstance(M, orbit_cheby.MSC)
    
    return MSCs

