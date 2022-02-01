    # -*- coding: utf-8 -*-

'''
----------------------------------------------------------------------------
tests for cheby_checker/nbody
 - Here I focus on the functions that PARSE the input data.
 - In the accompanying test_nbody_run script I test the functions for RUNNING the nbody integrations.

Jan 2022
Matthew Payne

Prev Work:
Mike Alexandersen, Matthew Payne & Matthew Holman

This code simplified as of Jan 2022
Removing many tests of non-json input
 - The non-json input methods *may* still work, but for now I just want to ensure that the json inputs work
 - Old tests of the non-json input can be found in the tests/archiaic/ttestt_nbody* file(s)
 
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
data_dir = os.path.join(repo_dir, 'dev_data')
json_dir = os.path.join(data_dir, 'mpc_orb_jsons')
test_dir = os.path.join(repo_dir, 'tests')
code_dir = os.path.join(repo_dir, 'cheby_checker')
for d in [test_dir, code_dir]:
    sys.path.append( d )

# import the main mnbody code that we want to test ...
#from cheby_checker
import nbody

# old conversion library that may be useful for cross-comparison of various tests ...
import MPC_library as mpc




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




# Tests of ParseElements
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
        N.non_grav_array         == []  and \
        N.output_times       == None  and \
        N.output_vectors     == None  and \
        N.output_n_times     == None  and \
        N.output_n_particles == None
            
    


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
    # Instantiate
    N = nbody.NbodySim()
        
    # Check that a few key attributes have the expected values before running the test
    # (these are the variables we expect to get populated by the function call below)
    assert N.helio_ecl_vec_EXISTS   is False
    assert N.helio_ecl_vec is None
    assert N.helio_ecl_cov_EXISTS   is False
    assert N.helio_ecl_cov is None
    assert N.integration_epoch  == None
    
    # call _parse_orbfit_json [this is the function we are testing]
    N._parse_orbfit_json( os.path.join(json_dir , data_file) )
    
    # Check that the expected attributes have now been populated
    assert N.helio_ecl_vec_EXISTS   is True
    assert isinstance(N.helio_ecl_vec, np.ndarray)
    assert N.helio_ecl_vec.ndim == 2
    assert N.helio_ecl_vec.shape == (1,6)
    
    assert N.helio_ecl_cov_EXISTS   is True
    assert isinstance(N.helio_ecl_cov, np.ndarray)
    assert N.helio_ecl_cov.ndim == 3
    assert N.helio_ecl_cov.shape in [(1,6,6),(1,7,7),(1,8,8),(1,9,9)]
    
    # Check that some of the numerical values are as expeected...
    with open(os.path.join(json_dir , data_file)) as f:
        data_dict = json.load(f)

        # - We expect the helio_ecl_vec to be the same as the data_file -> data_dict["CAR"]["elements"].values()
        SIMILAR, error = similar_xyzuvw( N.helio_ecl_vec , list(data_dict["CAR"]["elements"].values()) )
        assert SIMILAR, f'helio_ecl_vec != data_dict["CAR"]["elements"].values(): Diff=={error}'

        # - We expect the output time-object to be an astropy-Time object instantiated using data_file -> data_dict["epoch_data"]["epoch"]
        T = Time(float(data_dict["epoch_data"]["epoch"]) , format='mjd', scale='tt')
        assert N.integration_epoch.tdb.to_value("jd") == T.tdb.to_value("jd")


def test_parse_inputs_run_mpcorb_A():
    '''
    Test the *parse_inputs_run_mpcorb* routine
    This should
    (1) test that OrbFit files get parsed correctly (as in test_parse_orbfit_json_A, above)
    (2) test that the other variables gett populated properly
        (tstart, tstop, input_n_particles, ...)
    '''
    # Instantiate
    N = nbody.NbodySim()

    # Declare inputs
    tstart = 50000
    tstop  = 55000
    data_file = '545808fel_num.json'
    mpcorb_list = [ os.path.join(json_dir , data_file) ]

    # Check that a few key attributes have the expected values before running the test
    # (these are the variables we expect to get populated by the function call below)
    assert N.helio_ecl_vec_EXISTS   is False
    assert N.helio_ecl_vec is None
    assert N.helio_ecl_cov_EXISTS   is False
    assert N.helio_ecl_cov is None
    assert N.tstart             == None
    assert N.tstop              == None
    assert N.input_n_particles  == None
    assert N.integration_epoch  == None
    assert N.mpcorb_list        == []

    # Call the function
    N._parse_inputs_run_mpcorb(tstart = tstart , tstop = tstop , mpcorb_list = mpcorb_list)


    # Check that the expected attributes have now been populated
    assert N.tstart == tstart
    assert N.tstop  == tstop
    assert N.mpcorb_list == mpcorb_list
    assert N.input_n_particles == len(mpcorb_list)
    
    assert N.helio_ecl_vec_EXISTS   is True
    assert isinstance(N.helio_ecl_vec, np.ndarray)
    assert N.helio_ecl_vec.ndim == 2
    assert N.helio_ecl_vec.shape == (1,6)
    
    assert N.helio_ecl_cov_EXISTS   is True
    assert isinstance(N.helio_ecl_cov, np.ndarray)
    assert N.helio_ecl_cov.ndim == 3
    assert N.helio_ecl_cov.shape in [(1,6,6),(1,7,7),(1,8,8),(1,9,9)]

    # Check that some of the numerical values are as expected...
    # - We expect the helio_ecl_vec to be the same as the data_file -> data_dict["CAR"]["elements"].values()
    with open(os.path.join(json_dir , data_file)) as f:
        data_dict = json.load(f)
        SIMILAR, error = similar_xyzuvw( N.helio_ecl_vec , list(data_dict["CAR"]["elements"].values()) )
        assert SIMILAR, f'helio_ecl_vec != data_dict["CAR"]["elements"].values(): Diff=={error}'


def test_make_bary_equatorial_A():
    '''
    Test the creation of the barycentric elements / covar-components
    This requries that we initially populate heliocentric components, ...
    ... and then do the transformation using *make_bary_equatorial*
    '''

    # Instantiate
    N = nbody.NbodySim()

    # Declare inputs
    tstart = 50000
    tstop  = 55000
    data_file = '545808fel_num.json'
    mpcorb_list = [ os.path.join(json_dir , data_file) ]

    # Check that a few key attributes have the expected values before running the test
    # (these are the variables we expect to get populated by the function call below)
    assert N.bary_eq_vec_EXISTS     is False
    assert N.bary_eq_vec == None
    assert N.bary_eq_cov_EXISTS     is False
    assert N.bary_eq_cov == None

    # Call the parse function (already tested above) to get the data set-up
    N._parse_inputs_run_mpcorb(tstart = tstart , tstop = tstop , mpcorb_list = mpcorb_list)

    # NOW CALL THE ROUTINE WE WANT TO TEST: make_bary_equatorial
    N.make_bary_equatorial()
    
    # Check that the expected attributes have now been populated
    assert N.bary_eq_vec_EXISTS     is True
    assert isinstance(N.bary_eq_vec, np.ndarray)
    assert N.bary_eq_cov_EXISTS     is True
    assert isinstance(N.bary_eq_cov, np.ndarray)


def test_make_bary_equatorial_B():
    '''
    Trying to test the accuracy of the numerical components of the vectors & covariance-matrix
    produced by the *make_bary_equatorial* function
    '''
    pass
    
    

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

