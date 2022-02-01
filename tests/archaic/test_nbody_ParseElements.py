    # -*- coding: utf-8 -*-
# /tests/test_nbody.py

'''
----------------------------------------------------------------------------
tests for mpc_nbody

Dec 2021
Matthew Payne

Prev Work:
Mike Alexandersen, Matthew Payne & Matthew Holman

This code simplified as of Dec 2021
Removing many tests of non-json input
 - The non-json input methods *may* still work, but for now I just want to ensure that the json inputs work
 
 
----------------------------------------------------------------------------
'''

# import third-party packages
# -----------------------------------------------------------------------------
import sys
import os
import numpy as np
import pytest
from astroquery.jplhorizons import Horizons
from astropy.time import Time
import pytest
from filecmp import cmp
import getpass
import json



# Import neighbouring packages
# -----------------------------------------------------------------------------
try:  # Import ephem_forces from whereever REBX_DIR is set to live
    sys.path.append(os.environ['REBX_DIR'])
    from examples.ephem_forces import ephem_forces
except (KeyError, ModuleNotFoundError):
    from reboundx.examples.ephem_forces import ephem_forces

# cheby_checker/                 # <<-- repo
# cheby_checker/cheby_checker    # <<-- python
# cheby_checker/tests            # <<-- tests
this_dir = os.path.abspath(os.path.dirname( __file__ ))
repo_dir = os.path.abspath(os.path.dirname( this_dir ))
test_dir = os.path.join(repo_dir, 'tests')
code_dir = os.path.join(repo_dir, 'cheby_checker')
for d in [test_dir, code_dir]:
    sys.path.append( d )

# import the main mnbody code that we want to test ...
#from cheby_checker
import nbody

# old conversion library that may be useful for cross-comparison of various tests ...
from code_dir import MPC_library as mpc




# Constants & Test Data
# -----------------------------------------------------------------------------
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__))), 'dev_data')



# Utility functions to help with testing
# -----------------------------------------------------------------------------
def _get_and_set_junk_data(P, BaryEqDirect=False ):
    """
    For an input P = ParseElements-class object,
    Just populate some junk data within the obbject
    Of use for testing some save function(s).
    """
    P.time = Time(2458849.5, format='jd', scale='tdb')
    v      = np.array( [[3., 2., 1., 0.3, 0.2, 0.1]] )
    CoV    = 0.01 * np.ones((1,6,6))
    
    # Default is to make helio-ecl, then calc bary-eq from that
    if not BaryEqDirect:
        P.helio_ecl_vec              = v
        P.helio_ecl_vec_EXISTS       = True
        
        P.helio_ecl_cov              = CoV
        P.helio_ecl_cov_EXISTS       = True
    
        P.make_bary_equatorial()
        
    # Alternative is to directly set bary-eq
    else:
        P.bary_eq_vec                = v
        P.bary_eq_vec_EXISTS         = True
        
        P.bary_eq_cov                = CoV
        P.bary_eq_cov_EXISTS         = True


def is_parsed_good_enough(new_results_file, expected_results_file):
    '''
    Helper function to help test whether a just-created "new_results_file" file matches
    the "expected_results_file" in the "dev_data" directory
    '''
    
    if cmp(new_results_file, expected_results_file):
        assert True  # If files are identical, no further testing needed.
        
    else:  # If files not identical, investigate further:
        with open(new_results_file, 'r') as fileA, open(expected_results_file, 'r') as fileB :
            five_tf = []
            for _ in range(0, 5):  # First five lines should be identical
                lineA = fileA.readline()
                lineB = fileB.readline()
                five_tf.append(lineA == lineB)
            xyzA = np.array(fileA.readline().split(), dtype=float)
            xyzB = np.array(fileB.readline().split(), dtype=float)
            vA = np.array(fileA.readline().split(), dtype=float)
            vB = np.array(fileB.readline().split(), dtype=float)
            error, good_tf = compare_xyzv(np.concatenate([xyzA, vA]),
                                          np.concatenate([xyzB, vB]),
                                          1e-13, 1e-14)  # 15 mm, 1.5 mm/day
            if np.all(good_tf) & np.all(five_tf):
                pass # print('Awesome!')
            else:
                print(f'\n Problem detected in *is_parsed_good_enough* ... ')
                print(f'new_results_file={new_results_file}, expected_results_file={expected_results_file}')
                print(f'First five lines identical: {five_tf:}')
                print(f'Position off by: {error[:3]:}')
                print(f'Velocity off by: {error[3:6]:}')
            assert np.all(good_tf) & np.all(five_tf)


def compare_xyzv(xyzv0, xyzv1, threshold_xyz, threshold_v):
    '''
    Calculate the difference between two sets of cartesian coordinates.
    '''
    if isinstance(xyzv0, list):
        xyzv0 = np.array(xyzv0)
    if isinstance(xyzv1, list):
        xyzv1 = np.array(xyzv1)
    error = xyzv0 - xyzv1
    good_tf = np.abs(error) < np.array([threshold_xyz] * 3 + [threshold_v] * 3)
    return error, good_tf




# Tests of ParseElements
# -----------------------------------------------------------------------------

"""
@pytest.mark.parametrize(   ('data_file'),
                         [  '30101.eq0_postfit',
                            '30102.eq0_postfit',
                            '30101.eq0_horizons',
                            '30102.eq0_horizons'][:1])
def test_parse_orbfit_felfile_txt(data_file):

    '''
    Test that OrbFit files get parsed correctly.
    NB: The ...eq0... files passed in (above) are
        the OLD text-file output from ORBFIT orbit-fitting
        They have filetypes like .eq0/.eq1
    '''
    P = nbody.ParseElements()
    
    # Check that the expected attributes exist
    # and that they are initiated == None
    assert P.helio_ecl_vec_EXISTS   is False
    assert P.helio_ecl_vec          is None
    assert P.helio_ecl_cov_EXISTS   is False
    assert P.helio_ecl_cov          is None
    
    # Read the contents of the test file
    # We are doing this here because we are explicitly testing ONLY the
    #    *parse_orbfit_felfile_txt* function below
    with open(os.path.join(DATA_DIR, data_file),'r') as fh:
        file_contents=fh.readlines()

    # call parse_orbfit_felfile_txt
    P.parse_orbfit_felfile_txt(file_contents, CHECK_EPOCHS=False)
    
    # Check that the expected attributes exist
    # and that they are populated
    
    assert P.helio_ecl_vec_EXISTS   is True
    
    assert isinstance(P.helio_ecl_vec, np.ndarray)
    assert P.helio_ecl_vec.ndim == 2
    assert P.helio_ecl_vec.shape == (1,6)
    
    assert P.helio_ecl_cov_EXISTS   is True
    
    assert isinstance(P.helio_ecl_cov, np.ndarray)

    assert P.helio_ecl_cov.ndim == 3
    assert P.helio_ecl_cov.shape == (1,6,6)
"""


@pytest.mark.parametrize(   ('data_file'),
                         [  '10199fel_num.json',
                            '1566fel_num.json',
                            '2003AF23fel_num.json',
                            '2017AP4fel_num.json',
                            '545808fel_num.json'])
def test_parse_orbfit_json_A(data_file):

    '''
    Test that OrbFit files get parsed correctly.
    NB(1): The ...json... files passed in (above) are
        the mpcorb format jsons derived from ORBFIT orbit-fitting
    NB(2): This test deliberately only works for 6-dimension stuff, i.e. gravity-only
    '''
    P = nbody.ParseElements()
    
    # Check that the expected attributes exist
    assert P.helio_ecl_vec_EXISTS   is False
    assert P.helio_ecl_vec          is None
    assert P.helio_ecl_cov_EXISTS   is False
    assert P.helio_ecl_cov          is None
    
    # Read the contents of the test file
    # We are doing this here because we are explicitly testing ONLY the
    #    *parse_orbfit_json* function below
    with open(os.path.join(DATA_DIR, data_file),'r') as json_file:
        file_contents = json.load(json_file)

    
    # call parse_orbfit_json
    P.parse_orbfit_json(file_contents, CHECK_EPOCHS=False)
    
    # Check that the expected attributes exist
    # and that they are populated
    
    assert P.helio_ecl_vec_EXISTS   is True
    
    assert isinstance(P.helio_ecl_vec, np.ndarray)
    assert P.helio_ecl_vec.ndim == 2
    assert P.helio_ecl_vec.shape == (1,6)
    
    assert P.helio_ecl_cov_EXISTS   is True
    
    assert isinstance(P.helio_ecl_cov, np.ndarray)

    assert P.helio_ecl_cov.ndim == 3
    assert P.helio_ecl_cov.shape in [(1,6,6),(1,7,7),(1,8,8),(1,9,9)]
    


"""
def test_save_elements():
    '''Test that saving input-elements to an outpuut-file works correctly.'''
    # Get rid of an save_file.tmp file in the test directory
    if os.path.isfile('save_file.tmp'):
        os.remove('save_file.tmp')
        
    # Instantiate ...
    P = nbody.ParseElements()
    
    # Populate variables (junk data)
    _get_and_set_junk_data(P , BaryEqDirect=True)
    
    # Save to file
    P.save_elements()
    
    # Check contents of file are as expected
    assert cmp('./save_file.tmp', os.path.join(DATA_DIR, 'expected_junk_save.dat'))
    
    # Get rid of an save_file.tmp file in the test directory
    if os.path.isfile('save_file.tmp'):
        os.remove('save_file.tmp')





names_of_variables     = ('data_file', 'file_type', 'test_result_file')
values_for_each_test   = [
    pytest.param('30101.ele220', 'ele220', 'holman_ic_30101',
                 marks=pytest.mark.xfail(reason='Not implemented yet.')),
    pytest.param('30102.ele220', 'ele220', 'holman_ic_30102',
                 marks=pytest.mark.xfail(reason='Not implemented yet.')),
    ('30101.eq0_postfit', 'eq', 'holman_ic_30101'),
    ('30102.eq0_postfit', 'eq', 'holman_ic_30102'),
    ('30101.eq0_horizons', 'eq', 'holman_ic_30101_horizons'),
    ('30102.eq0_horizons', 'eq', 'holman_ic_30102_horizons'),
    ('10199fel_num.json',  'json', 'expected_10199fel_num.txt'),
 ]
@pytest.mark.parametrize( names_of_variables, values_for_each_test )
def test_instantiation_from_data_files(data_file, file_type, test_result_file):
    '''
    Test that instantiation with data works
    (by doing so we essentially test most/all functionalities).
    '''
    # Where we will save a local test file
    save_file='save_file.tmp'
    if os.path.isfile(save_file) : os.remove(save_file)
    
    # Instantiate from file (which calls *make_bary_equatorial*)
    # and then save to save_file='save_file.tmp'
    nbody.ParseElements( input = os.path.join(DATA_DIR, data_file),
                         filetype = file_type,
                         save_parsed=True,
                         save_file=save_file,
                         CHECK_EPOCHS=False )

    # Check the output
    is_parsed_good_enough( save_file , os.path.join(DATA_DIR, test_result_file) )
    
    # Tidy-up by removing any local test file
    if os.path.isfile(save_file) : os.remove(save_file)

"""

# End 
