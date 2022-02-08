# -*- coding: utf-8 -*-

'''
    --------------------------------------------------------------
    tests of orbit_cheby's base class
    ../tests/test_orbit_cheby_locations.py
    
    Testing only the _define_locations()
    funuction / quantities within orbit_cheby.MSC

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
# cheby_checker/                 # <<-- repo
# cheby_checker/cheby_checker    # <<-- python
# cheby_checker/tests            # <<-- tests
this_dir = os.path.abspath(os.path.dirname( __file__ ))
repo_dir = os.path.abspath(os.path.dirname( this_dir ))
data_dir = os.path.join(repo_dir, 'dev_data')
json_dir = os.path.join(data_dir, 'mpc_orb_jsons')
std_json_dir = os.path.join(json_dir, 'standard_mp')
test_dir = os.path.join(repo_dir, 'tests')
code_dir = os.path.join(repo_dir, 'cheby_checker')
for d in [test_dir, code_dir]:
    sys.path.append( d )

import nbody
import orbit_cheby
import cheby_checker
import obs_pos
import coco

# Constants & Test Data
# -----------------------------------------------------------------------------
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'dev_data')

FLAT_FILES = [  os.path.join(DATA_DIR, '2022AA_demo.txt') ,
                os.path.join(DATA_DIR, 'simulation_states.dat')]
orbfit_filenames = [os.path.join(DATA_DIR, file) for file in ['30101.eq0_horizons', '30102.eq0_horizons']]
mpc_orb_json_files = glob.glob(std_json_dir + '/*.json' )



# Actual tests : MSC Object [ tests of MSC_Loader object are below]
# -----------------------------------------------------------------

def test_create_empty_MSC():
    ''' Create an empty MSC class object'''
    
    # Initialize the multi_sector_cheby
    result = orbit_cheby.MSC()
    
    # Check the result is as expected
    assert isinstance( result , orbit_cheby.MSC )
    assert 'sector_coeffs' in result.__dict__
    assert isinstance( result.sector_coeffs , dict )




def test_define_locations_A(  ):
    '''
        
    '''
    
    # Instantiate MSC object
    M=orbit_cheby.MSC()

    # Define number of coordinates (reqd in order to allow *_define_locations_A* to work)
    M.n_coordinates = 6
    
    # Call the *_define_locations* function that we want to test
    M._define_locations()
    
    # Check various attributes have been populated
    assert hasattr(M,'triangular_mapping') and isinstance( M.triangular_mapping , dict)
    assert hasattr(M,'coord_map')          and isinstance( M.coord_map , dict)
    assert hasattr(M,'inv_map')            and isinstance( M.inv_map , dict)
    assert hasattr(M,'tri_map')            and isinstance( M.tri_map , dict)
    assert hasattr(M,'combi_map')          and isinstance( M.combi_map , dict)
    assert hasattr(M,'XYZ_slice_spec')
    assert hasattr(M,'covXYZ_slice_spec')

    
def test_define_locations_B(  ):
    '''
        ...
    '''
    
    for n in range(6,10):
    
        # Instantiate MSC object
        M=orbit_cheby.MSC()

        # Define number of coordinates (reqd in order to allow *_define_locations_A* to work)
        M.n_coordinates = n
        
        # Call the *_define_locations* function that we want to test
        M._define_locations()
        
        # Check various attributes have been populated
        assert hasattr(M,'triangular_mapping') and isinstance( M.triangular_mapping , dict)
        assert hasattr(M,'coord_map')          and isinstance( M.coord_map , dict)
        assert hasattr(M,'inv_map')            and isinstance( M.inv_map , dict)
        assert hasattr(M,'tri_map')            and isinstance( M.tri_map , dict)
        assert hasattr(M,'combi_map')          and isinstance( M.combi_map , dict)
        assert hasattr(M,'XYZ_slice_spec')
        assert hasattr(M,'covXYZ_slice_spec')

    
def test_define_locations_C(  ):
    '''
        ...
    '''
    
    # Instantiate MSC object
    M=orbit_cheby.MSC()

    # Define number of coordinates (reqd in order to allow *_define_locations_A* to work)
    M.n_coordinates = 6
    
    # Call the *_define_locations* function that we want to test
    M._define_locations()
    
    # Check various attributes are as expected
    assert M.coord_map == { 'x':0,'y':1,'z':2,'vx':3,'vy':4,'vz':5}
    assert M.inv_map   == { 0:'x', 1:'y',2:'z',3:'vx',4:'vy',5:'vz'}
    assert np.array_equal(M.covXYZ_slice_spec , np.array([6,7,8,12,13,17]) )



def test_take_triangular_A():
    '''
    '''
    
    for n in range(6,10):
    
        # Instantiate MSC object
        M=orbit_cheby.MSC()

        # Define number of coordinates (reqd in order to allow *_define_locations_A* to work)
        M.n_coordinates = n

        # Call the *_define_locations* function that creates some useful quantities that we need
        M._define_locations()

        # Create a stack of symmetric covariance matrices of the correct shape to use as test input
        Nt = 41
        c = np.arange( Nt*n*n ).reshape(Nt,n,n)
        i, j = np.triu_indices( n )
        for k,_ in enumerate(c):
            c[k] = _ + _.T - np.diag(_.diagonal())   ### This makes each time-slice symmetric

        ####################
        # Call the conversion function tthat we want to test
        t = M._take_triangular( c )
        ####################
        
        # Check the output ...
        No = [k for k,v in M.triangular_mapping.items() if v == n][0]
        assert t.shape == (Nt,No)


def test__make_square_A():
    '''
    '''
    
    for n in range(6,10):
    
        # Instantiate MSC object
        M=orbit_cheby.MSC()

        # Define number of coordinates (reqd in order to allow *_define_locations_A* to work)
        M.n_coordinates = n

        # Call the *_define_locations* function that creates some useful quantities that we need
        M._define_locations()

        # Create a stack of symmetric covariance matrices of the correct shape to use as test input
        Nt = 41
        c = np.arange( Nt*n*n ).reshape(Nt,n,n)
        i, j = np.triu_indices( n )
        for k,_ in enumerate(c):
            c[k] = _ + _.T - np.diag(_.diagonal())   ### This makes each time-slice symmetric

        ####################
        # Call the two conversion functions that we want to test
        t = M._take_triangular( c )
        s = M._make_square( t )
        t2= M._take_triangular( s )
        ####################
        
        # Check the output ...
        # (i) As above, check that the triangulr part from _take_triangular is of the expected shape
        No = [k for k,v in M.triangular_mapping.items() if v == n][0]
        assert t.shape == (Nt,No)
        # (ii) Check that the "round trip" of _take_triangular + _make_square leaves the array unchanged
        assert s.shape == c.shape
        assert np.allclose(c,s)
        # (iii) And check again that a further iteration gets the same triangular result ...
        assert t.shape == t2.shape
        assert np.allclose(t,t2)

