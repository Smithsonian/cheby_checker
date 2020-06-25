# -*- coding: utf-8 -*-
# sifter/tests/test_base

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
from tests.test_parse_input import is_parsed_good_enough, compare_xyzv
from cheby_checker import orbit_cheby
from cheby_checker import nbody_reader
from cheby_checker import mpc_nbody



# Constants & Test Data
# -----------------------------------------------------------------------------
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'dev_data')
FLAT_FILES = [os.path.join(DATA_DIR, '2022AA_demo.txt') , os.path.join(DATA_DIR, 'simulation_states.dat')]


# Convenience data / functions to aid testing
# --------------------------------------------------------------

def test_convenience_call_to_mpc_nbody_nbodysim():
    ''' 
        A convenience function to return a simulation object
        Proper testing of mpc_nbody is done elsewhere (test_run_nbody.py)
    '''
    # Define some files that have data in them
    filenames = [os.path.join(DATA_DIR, file)
                 for file in ['30101.eq0_horizons', '30102.eq0_horizons']]

    # First, let's initiate the class with an input file:
    Sim = mpc_nbody.NbodySim(filenames[0], 'eq')

    # Now run the integrator, by calling the object.
    Sim(tstep=20, trange=1000)
    
    # Quick tests ...
    for attrib in ["output_times", "output_vectors"]:
        assert attrib in Sim.__dict__

    return Sim


def create_single_MSC_from_array() :
    return orbit_cheby.MSC_Loader(FROM_ARRAY = True ,
                                 unpacked_provisional_designations = name,
                                 times_TDB = times,
                                 statearray = states).MSCs[0]

def create_list_of_MSCs_from_arrays(n=5) :
    arrays, names = [],[]
    for i in range(1,n):
        name, times, states  = nbody_reader.parse_nbody_txt( text_filepath )
        names.append(name+"_"+str(i))
        arrays.append(states)
    states_3D = np.stack(arrays, axis=2)

    return orbit_cheby.MSC_Loader(FROM_ARRAY = True ,
                                     unpacked_provisional_designations = names,
                                     times_TDB = times,
                                     statearray = states_3D).MSCs


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
                            (FLAT_FILES[1] , (0,624) )
                          ]
                         )
def test_from_coord_arrays( flat_file , expected_sectors ):
    '''
        Use MSC functionality to load from numpy arrays
        Significant reliance on *generate_cheb_for_sector()*
        Need to test separately
    '''
    print(f'flat_file:{flat_file}')
    assert False
    
    
    
    '''
    # Get the data from supplied test files
    name, times, states  = nbody_reader.parse_nbody_txt( flat_file )

    # Instantiate
    M=orbit_cheby.MSC()

    # attempt to load from arrays ...
    primary_unpacked_provisional_designation, TDB_init , TDB_final = name, 2440000, 2459999
    M.from_coord_arrays(primary_unpacked_provisional_designation, TDB_init , TDB_final , times , states)

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

    print(flat_file,  len( M.sector_coeffs) , M.TDB_init, M.TDB_final, M.sector_init, M.sector_final , expected_sectors)

    assert False
    '''

def test_generate_cheb_for_sector():
    ''' 
        Use MSC function to generate the chebyshev coefficients for a given (32 day) sector
        Need to test the hell out of this, as it underlies everything
        Currently only testing it on 2 input files
    '''

    # Get the data from some different test files
    for f in FLAT_FILES :

        # - Note the use of RELATIVE TIMES
        name, times, states  = nbody_reader.parse_nbody_txt( f )
        rel_times = times - times[0]
    
        # Instantiate
        M=orbit_cheby.MSC()

        # Loop over some different starting times:
        # - want to test different parts of the data
        delta = 10
        for n in range(delta):
        
            # indicees to use to select from data (ensure that its in allowed range)
            ind = np.min( [n * int(len(times)/(delta+1)), len(times)-32] ) + np.arange(32)
        
            # Call *generate_cheb_for_sector* using sub-arrays
            cheb_coeffs = M.generate_cheb_for_sector(rel_times[ind], states[ind])
        
            # Evaluate returned coeffs ...
            quickEval = np.polynomial.chebyshev.chebval(rel_times[ind]  , cheb_coeffs).T
            assert np.max( np.abs(quickEval - states[ind]) ) < M.maxerr

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
        ... tests need breaking down ...
    '''

    # Get the data from file
    name, times, states  = nbody_reader.parse_nbody_txt( FLAT_FILES[0] )
    
    # Attempt to instantiate MSC  via MSC_Loder
    MSCs = orbit_cheby.MSC_Loader(
                                  primary_unpacked_provisional_designations = name,
                                  times_TDB = times,
                                  statearray = states).MSCs
    assert isinstance(MSCs , list)
    assert len(MSCs) == 1
    
    M = MSCs[0]
    assert isinstance(M, orbit_cheby.MSC)


