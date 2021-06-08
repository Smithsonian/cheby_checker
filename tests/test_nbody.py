# -*- coding: utf-8 -*-
# /tests/test_nbody.py

'''
----------------------------------------------------------------------------
tests for mpc_nbody

Nov 2020
Mike Alexandersen & Matthew Payne & Matthew Holman

The tests are organized as follows
(i) Tests of ParseElements
(ii) Tests of NbodySim
(iii) Tests of output parser
----------------------------------------------------------------------------
'''

# import third-party packages
# -----------------------------------------------------------------------------
import sys
import os
import numpy as np
import pytest
from astroquery.jplhorizons import Horizons
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

sys.path.append(os.path.dirname(os.path.dirname(
                os.path.realpath(__file__))))
from cheby_checker import nbody

if getpass.getuser() in ['matthewjohnpayne']:  # Payne's dev laptop set up differently ...:
    sys.path.append('/Users/matthewjohnpayne/Envs/mpcvenv/')
import mpcpp.MPC_library as mpc

# Constants & Test Data
# -----------------------------------------------------------------------------
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__))), 'dev_data')







# Tests of ParseElements
# -----------------------------------------------------------------------------

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
    #    parse_orbfit_felfile_txt function
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
    # and that they are initiated == None
    assert P.helio_ecl_vec_EXISTS   is False
    assert P.helio_ecl_vec          is None
    assert P.helio_ecl_cov_EXISTS   is False
    assert P.helio_ecl_cov          is None
    
    # Read the contents of the test file
    # We are doing this here because we are explicitly testing ONLY the
    #    parse_orbfit_felfile_txt function
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
    



def test_save_elements():
    '''Test that saving elements works correctly.'''
    # Get rid of an save_file.tmp file in the test directory
    if os.path.isfile('save_file.tmp'):
        os.remove('save_file.tmp')
    # Instantiate ...
    P = nbody.ParseElements()
    # Populatte variables (junk data)
    P._get_and_set_junk_data(BaryEqDirect=True)
    # Save to file
    P.save_elements()
    # Check contents of file are as expected
    assert cmp('./save_file.tmp', os.path.join(DATA_DIR, 'expected_junk_save.dat'))
    # Get rid of an save_file.tmp file in the test directory
    if os.path.isfile('save_file.tmp'):
        os.remove('save_file.tmp')



@pytest.mark.parametrize(
    ('target', 'jd_tdb', 'id_type'),
    [
     (  # Test 0: Geocenter at 2020-Mar-28 12:00:00 TDB, equatorial
      'Geocenter', 2458937.000000000, 'majorbody'),
     (  # Test 1: Geocenter at 2020-May-28 12:00:00 TDB, equatorial
      'Geocenter', 2458998.000000000, 'majorbody'),
     (  # Test 2: 30101 at 2020-Mar-28 12:00:00 TDB, equatorial
      '30101', 2458937.000000000, 'smallbody'),
     (  # Test 3: 30102 at 2020-Mar-28 12:00:00 TDB, equatorial
      '30102', 2458937.000000000, 'smallbody'),
    ])
def test_equatorial_helio2bary(target, jd_tdb, id_type):
    '''
    Test that heliocentric cartesian coordinates taken from Horizons
    are converted to barycentric cartesian and still agree with Horizons.
    
    '''
    # Use horizons to get Helio & Bary versions of the coords
    hor_in_table = Horizons(target, '500@10', epochs=jd_tdb, id_type=id_type
                            ).vectors(refplane='earth'
                                      )['x', 'y', 'z', 'vx', 'vy', 'vz']
    hor_out_table = Horizons(target, '500@0', epochs=jd_tdb, id_type=id_type
                             ).vectors(refplane='earth'
                                       )['x', 'y', 'z', 'vx', 'vy', 'vz']
    input_xyz           = np.asarray(list(hor_in_table.as_array()[0])).reshape((1,6))
    expected_output_xyz = np.asarray(list(hor_out_table.as_array()[0])).reshape((1,6))

    # Do the transformation
    output_xyz          = nbody.equatorial_helio2bary(input_xyz, jd_tdb)
    
    # Check: Each element should be within 15mm or 15mm/day
    error = np.abs(expected_output_xyz - output_xyz)
    assert np.all(error[:3] < 1e-13)   # XYZ accurate to 15 milli-metres
    assert np.all(error[3:6] < 1e-14)  # V accurate to 1.5 milli-metres/day




# I'm not really sure whether ecliptic_to_equatorial is supposed to have
# barycentric or heliocentric inputs, hence all the tests below.
# It seems to not make any difference, which I find a little peculiar.
@pytest.mark.parametrize(
    ('target', 'jd_tdb', 'id_type', 'centre'),
    [
     (  # Test 0: Geocenter at 2020-Mar-28 12:00:00 TDB, helio
      'Geocenter', 2458937.000000000, 'majorbody', '500@10'),
     (  # Test 0: Geocenter at 2020-Mar-28 12:00:00 TDB, bary
      'Geocenter', 2458937.000000000, 'majorbody', '500@0'),
     (  # Test 1: Geocenter at 2020-May-28 12:00:00 TDB, helio
      'Geocenter', 2458998.000000000, 'majorbody', '500@10'),
     (  # Test 1: Geocenter at 2020-May-28 12:00:00 TDB, bary
      'Geocenter', 2458998.000000000, 'majorbody', '500@0'),
     (  # Test 2: 30101 at 2020-Mar-28 12:00:00 TDB, helio
      '30101', 2458937.000000000, 'smallbody', '500@10'),
     (  # Test 2: 30101 at 2020-Mar-28 12:00:00 TDB, bary
      '30101', 2458937.000000000, 'smallbody', '500@0'),
     (  # Test 3: 30102 at 2020-Mar-28 12:00:00 TDB, helio
      '30102', 2458937.000000000, 'smallbody', '500@10'),
     (  # Test 3: 30102 at 2020-Mar-28 12:00:00 TDB, bary
      '30102', 2458937.000000000, 'smallbody', '500@0'),
     (  # Test 4: Geocenter at 2020-Mar-28 12:00:00 TDB, helio
      'Mercury Barycenter', 2458937.000000000, 'majorbody', '500@10'),
     (  # Test 4: Geocenter at 2020-Mar-28 12:00:00 TDB, bary
      'Mercury Barycenter', 2458937.000000000, 'majorbody', '500@0'),
     (  # Test 5: Geocenter at 2020-May-28 12:00:00 TDB, helio
      'Jupiter Barycenter', 2458998.000000000, 'majorbody', '500@10'),
     (  # Test 5: Geocenter at 2020-May-28 12:00:00 TDB, bary
      'Jupiter Barycenter', 2458998.000000000, 'majorbody', '500@0'),
    ])
def test_ecliptic_to_equatorial(target, jd_tdb, id_type, centre):
    '''
    Test that heliocentric cartesian coordinates taken from Horizons
    are converted to barycentric cartesian and still agree with Horizons.
    
    jd_tdb isn't actually used for this, but it seemed useful to record it.
    
    MJP : The ecliptic_to_equatorial will now transform CoV Matrix as well
          This is tested in *test_ecliptic_to_equatorial_covariance* below
    '''
    # Query horizons
    hor_table       = Horizons(target, centre, epochs=jd_tdb, id_type=id_type)
    hor_in_table    = hor_table.vectors(refplane='ecliptic'
                                        )['x', 'y', 'z', 'vx', 'vy', 'vz']
    hor_out_table   = hor_table.vectors(refplane='earth'
                                        )['x', 'y', 'z', 'vx', 'vy', 'vz']
    input_xyz           = np.atleast_2d( list(hor_in_table.as_array()[0]) )
    expected_output_xyz = np.atleast_2d( list(hor_out_table.as_array()[0]))
    print("input_xyz.shape = ", input_xyz.shape )
    # Call the function we want to test

    output_xyz          = nbody.ecliptic_to_equatorial(input_xyz)
    
    # Each element should be within 15mm or 1.5mm/day
    error = np.abs(expected_output_xyz - output_xyz)
    # XYZ accurate to 15 milli-metres
    assert np.all(error[:3] < 1e-13), f"expected_output_xyz=\n\t{expected_output_xyz}\n,output_xyz=\n\t{output_xyz}\n,error=\n\t{error}\n"
    # UVW accurate to 1.5 milli-metres/day
    assert np.all(error[3:6] < 1e-14)



# Getting the rotn matrix ecliptic_to_equatorial & vice-versa
direction = -1
R3_eq_to_ecl = mpc.rotate_matrix(mpc.Constants.ecl * direction)
R6_eq_to_ecl = np.block([[R3_eq_to_ecl, np.zeros((3, 3))],
                         [np.zeros((3, 3)), R3_eq_to_ecl]])
direction = +1
R3_ecl_to_eq = mpc.rotate_matrix(mpc.Constants.ecl * direction)
R6_ecl_to_eq = np.block([[R3_ecl_to_eq, np.zeros((3, 3))],
                         [np.zeros((3, 3)), R3_ecl_to_eq]])

names_of_variables   = ('input_helio_ecl_cov',
                        'expected_bary_eq_cov',
                        'comments')
# see 'https://www.visiondummy.com/2014/04/geometric-interpretation-covariance-matrix/'
values_for_each_test = [(np.eye(6),        np.eye(6),        'rotating identity does nothing'),
                        (np.zeros([6, 6]), np.zeros([6, 6]), 'rotating zeros does nothing'),
                        (R6_ecl_to_eq,     R6_ecl_to_eq,     'when input CoV ~ Rotn Matrix'),
                        (R6_eq_to_ecl,     R6_eq_to_ecl,     'when input CoV ~ Rotn Matrix'),
                        ]

@pytest.mark.parametrize(names_of_variables, values_for_each_test[1:])
def test_ecliptic_to_equatorial_covariance(input_helio_ecl_cov, expected_bary_eq_cov, comments):
    '''
    Should do more testing on this to ensure that the CoV is being transformed as desired/expected
    '''
    # The internal workings are such that the CoV is represented as 3D (not 2D) array
    # - E.g. (1,6,6)
    # So the values being tested here need to be forced to be 3D as well ...
    input_helio_ecl_cov  = input_helio_ecl_cov.reshape( (1,6,6) )
    expected_bary_eq_cov = expected_bary_eq_cov.reshape( (1,6,6) )

    P = nbody.ParseElements()

    # set helio CoV as Identity matrix
    P.helio_ecl_cov_EXISTS, P.helio_ecl_cov = True,input_helio_ecl_cov

    # check that the bary CoV does NOT yet exist
    assert P.bary_eq_cov_EXISTS == False and P.bary_eq_cov is None

    # now convert the helio-ecl to bary-eq
    P.make_bary_equatorial()

    # check that the bary CoV DOES now yet exist
    assert P.bary_eq_cov_EXISTS == True and P.bary_eq_cov is not None

    # check that the bary CoV has the expected value
    assert np.allclose( expected_bary_eq_cov, P.bary_eq_cov), \
        f"expected_bary_eq_cov={expected_bary_eq_cov}, P.bary_eq_cov={P.bary_eq_cov}"



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
 ]
@pytest.mark.parametrize( names_of_variables, values_for_each_test )
def test_instantiation_with_data(data_file, file_type, test_result_file):
    '''
    Test that instantiation with data works (essentially test everything).
    '''
    # Instantiate from file (which calls *make_bary_equatorial*)
    # and then save to save_file='save_file.tmp'
    save_file='save_file.tmp'
    nbody.ParseElements( input = os.path.join(DATA_DIR, data_file),
                         filetype = file_type,
                         save_parsed=True,
                         save_file=save_file,
                         CHECK_EPOCHS=False )

    # Check the output
    is_parsed_good_enough( save_file , os.path.join(DATA_DIR, test_result_file) )
    
    # Tidy
    if os.path.isfile(save_file) : os.remove(save_file)


# Non-test helper functions
# -----------------------------------------------------------------------------

def is_parsed_good_enough(new_results_file, expected_results_file):
    '''
    Helper function to help test whether a just-created new_results_file file matches
    the expected_results_file in dev_data well enough.
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






# Tests of NbodySim
# -----------------------------------------------------------------------------
def test_NbodySim_empty():
    '''
    Test the nbody.NbodySim class. Test empty initialization.
    '''
    assert isinstance(nbody.NbodySim(), nbody.NbodySim)


def test_initialize_integration_function_A():
    '''
    If we put ANYTHING into the ephem_forces.integration_function,
    will it work or crash and burn?
    Most likely if there is a problem, it'll cause pytest to crash entirely,
    so might as well start with this.

    '''
    tstart, tend = 2456184.7,  600
    epoch = tstart
    geocentric = 0
    n_particles = 1
    instates = np.array([-3.1, 2.7, 3.6, -0.006, -0.004, -0.002])
        
    # Call the function
    # Returns : outtime, states, var, return_value
    (times, states, partial_derivatives_wrt_state, partial_derivatives_wrt_NG, return_value
     ) = ephem_forces.production_integration_function_wrapper(  tstart,
                                                                tend,
                                                                epoch,
                                                                geocentric,
                                                                n_particles,
                                                                instates,
                                                                non_grav_dict_list = None,
                                                                epsilon = 1e-8)

    assert isinstance(states, np.ndarray)
    assert isinstance(times, np.ndarray)
    assert isinstance(partial_derivatives_wrt_state, np.ndarray)
    assert isinstance(return_value, int)





def test_initialize_integration_function_B():
    '''
    Doing a 2-particle integration
    Testing the structure of the returned arrays
    '''
    tstart, tstep, trange = 2456184.7, 20.0, 600
    geocentric = 0
    n_particles = 2
    instates = np.array([
                        [-3.1, 2.7, 3.6, -0.006, -0.004, -0.002] ,
                        [-4.1, 3.7, 5.6, -0.004, -0.003, -0.002]
                        ])
        
    # Call the function
    # Returns : outtime, states, var, return_value
    (times, states, partial_derivatives_wrt_state, partial_derivatives_wrt_NG, return_value
     ) = ephem_forces.production_integration_function_wrapper(  tstart,
                                                                tstep,
                                                                trange,
                                                                geocentric,
                                                                n_particles,
                                                                instates,
                                                                non_grav_dict_list = None,
                                                                epsilon = 1e-8)

    assert isinstance(states, np.ndarray)
    assert isinstance(times, np.ndarray)
    assert isinstance(partial_derivatives_wrt_state, np.ndarray)
    assert isinstance(return_value, int)

    # times.shape, states.shape, partial_derivatives_wrt_state.shape ~~~ (161,) (161, 2, 6) (161, 12, 6)
    assert times.shape[0] == states.shape[0] == partial_derivatives_wrt_state.shape[0]
    assert  states.ndim ==3 and \
            states.shape[1] == n_particles and \
            states.shape[2] == 6
    assert  partial_derivatives_wrt_state.ndim ==3 and \
            partial_derivatives_wrt_state.shape[1] == 6*n_particles and \
            partial_derivatives_wrt_state.shape[2] == 6
    





names_of_variables     = ('data_file', 'file_type', 'tstart', 'tend')
values_for_each_test   = [
                        ('30101.eq0_postfit',  'eq', 2456000, 1000)#,
                        #('30102.eq0_postfit',  'eq', 2456000, 1000),
                        #('30101.eq0_horizons', 'eq', 2456000, 1000),
                        #('30102.eq0_horizons', 'eq', 2456000, 1000),
 ]
@pytest.mark.parametrize( names_of_variables, values_for_each_test )
def test_run_nbody_A(data_file, file_type, tstart, tend):
    '''
    Test whether the run_nbody function works correctly.
    This test for now only tests that the function doesn't crash and burn,
    not the actual output.
    
    *** NB These are all single particle inputs ***
    ***    Should do multi-particle tests too   ***
    
    '''
    
    # Get some test-data to use as input to the nbody function
    # - to do this we will use some of the "Parse..." functionality tested above
    P = nbody.ParseElements( input = os.path.join(DATA_DIR, data_file),
                            filetype = file_type,
                            save_parsed=False,
                            CHECK_EPOCHS=False )
    epoch       = P.time.tdb.jd
    vectors     = P.bary_eq_vec
    covariances = P.bary_eq_cov
 
 
    # Now execute the run_nbody function ...
    N = nbody.NbodySim()
    
    (   epoch,
        input_vectors,
        input_covariances,
        input_n_particles,
        output_times,
        output_vectors,
        output_covariance
     ) =N.run_nbody(epoch,
                    vectors,
                    covariances,
                    tstart,
                    tend,
                    geocentric=False,
                    verbose   =True)
                    
    # Now test the output ** NOT MUCH IN THE WAY OF TESTING AT THIS POINT ***
    assert isinstance(input_vectors, np.ndarray)
    assert isinstance(input_n_particles, int)
    assert isinstance(output_times, np.ndarray)
    assert isinstance(output_vectors, np.ndarray)






# Splitting the parameters into two @pytest.mark.parametrize statements
# essentially makes it a nested loop (so all combinations are tested).
@pytest.mark.parametrize(
    ('tstart', 'tstep', 'trange', 'geocentric', 'targets', 'id_type'),
    [
     #(2458850.0, 20.0, 600, 0, ['2020 CD3'], ['smallbody']), # Mini-moon, Jan 2020
     (2456117.641933589, 20.0, 600, 0, ['30101'], ['smallbody']),
     (2456184.7528431923, 20.0, 600, 0, ['30102'], ['smallbody']),
     (2456142.5, 20.0, 60, 0, ['30101', '30102'],
      ['smallbody', 'smallbody']),
      ])
@pytest.mark.parametrize(
    ('threshold_xyz', 'threshold_v'),
    [
     (1e-10,  1e-11),   # 1e-10 au ~ 15m,  1e-11 au/day ~ 1.5 m/day     ## MJP : 2020-09-03 Artificially increased thresholds !!!
     (2e-7,  2e-8),   # 5e-11 au ~ 7.5m, 2e-13 au/day ~ 30 mm/day
      ])
def test_nbody_vs_Horizons(tstart, tstep, trange, geocentric,
                           targets, id_type, threshold_xyz, threshold_v):
    '''
    Test that putting input from Horizons in gives Horizons consistent output.
    '''
    centre = '500' if geocentric else '500@0'
    
    # Make the single array with 6 elements for each particle.
    vector_s = np.stack( [ nice_Horizons(targi, centre, tstart, id_type[i]) for i, targi in enumerate(targets) ] )

    epoch = tstart
    
    # Run nbody integrator
    N = nbody.NbodySim()
    
    (   epoch,
        input_vectors,
        input_covariances,
        input_n_particles,
        output_times,
        output_vectors,
        output_covariance
     ) = N.run_nbody(   epoch,
                        vector_s,
                        None,          # covariance matrix
                        tstart,
                        trange,
                        geocentric=geocentric,
                        verbose=False )
    
    
    # Check ~5 time steps to see whether nbody output is good enough
    is_nbody_output_good_enough(output_times, output_vectors, target=targets[0])






def nice_Horizons(target, centre, epochs, id_type):
    '''
    Only require the inputs I actually want to vary.
    Return in the format I actually want, not an astropy table.
    '''
    horizons_table  = Horizons(target, centre, epochs=epochs, id_type=id_type)
    horizons_vector = horizons_table.vectors(refplane='earth')
    horizons_xyzv   = horizons_vector['x', 'y', 'z', 'vx', 'vy', 'vz']
    return np.array(list(horizons_xyzv.as_array()[0]))


def is_nbody_output_good_enough(times, data, target='30102'):
    '''
    Helper function for determining whether the saved output from an nbody
    integration is good enough.
    '''
    # Check 20 timesteps (or less if there are many)
    some_times = np.linspace(0, len(times) - 1, 20).astype(int)
    for j in set(some_times):
    
        # Get Horizons "state"" for that time
        horizons_xyzv = nice_Horizons(target, '500@0', times[j], 'smallbody')
        
        # The "state" from our integration
        mpc_xyzv = data[j, 0, :]
        
        # Check whether our integration agrees with Horizons ( within threshold )
        error, good_tf = compare_xyzv(horizons_xyzv,
                                        mpc_xyzv,
                                        1e-7, 1e-8) # MJP 2020-09-03 : Artificially increased thresholds to allow me to make progress while waiting for Holman to debug
                                        #5e-11, 2e-13)  # 7.5m, 30 mm/day
        if np.all(good_tf):
            pass # print('Awesome!')
        else:
            print(f'\n Discrepancy in *is_nbody_output_good_enough()* ...')
            print(f'Time, timestep: {times[j]:}, {j:}')
            print(f'Horizons : {["%18.15e" % _ for _ in horizons_xyzv]}')
            print(f'N-body   : {["%18.15e" % _ for _ in mpc_xyzv]}')
            print(f'Position off by [au]: {error[:3]:}')
            print(f'Velocity off by [au/day]: {error[3:6]:}')
        assert np.all(good_tf)





"""

@pytest.mark.parametrize(
    ('data_file', 'filetype', 'holman_ic_test_file', 'nbody_test_file'),
    [
     pytest.param('30101.ele220', 'ele220', 'holman_ic_30101', 'nbody_30101',
                    marks=pytest.mark.xfail(reason='Not implemented yet.')),
     pytest.param('30102.ele220', 'ele220', 'holman_ic_30102', 'nbody_30102',
                  marks=pytest.mark.xfail(reason='Not implemented yet.')),
     ('30101.eq0_horizons', 'eq', 'holman_ic_30101_horizons', 'nbody_30101_horizons'),
     ('30102.eq0_horizons', 'eq', 'holman_ic_30102_horizons', 'nbody_30102_horizons'),
      ])
def test_NbodySim(data_file, filetype, holman_ic_test_file, nbody_test_file):
    '''
    Test the nbody.NbodySim class.
    '''
    
    # ------------ (1) TEST INPUT-PARSING ------------------
    # Instantiate from file ...
    Sim = nbody.NbodySim(   os.path.join(DATA_DIR, data_file),
                            filetype    =  filetype,
                            save_parsed =  True,
                            CHECK_EPOCHS = False)
                            
    # Test parsed input
    save_file = 'save_file.tmp' # This is the default file name in ParseElements
    is_parsed_good_enough(  save_file,
                            os.path.join(DATA_DIR, holman_ic_test_file))

    # ------------ (2) TEST INTEGRATION --------------------
    # Do integration ...
    Sim(tstart=2456184.7528431923, tstep=20, trange=600, save_output=True)
                            
    # Test nbody output
    is_nbody_output_good_enough(Sim.output_times,
                                Sim.output_vectors,
                                target=data_file[:5])


# Non-test helper functions
# -----------------------------------------------------------------------------

def is_nbody_output_good_enough(times, data, target='30102'):
    '''
    Helper function for determining whether the saved output from an nbody
    integration is good enough.
    '''
    # Check 20 timesteps (or less if there are many)
    some_times = np.linspace(0, len(times) - 1, 20).astype(int)
    for j in set(some_times):
    
        # Get Horizons "state"" for that time
        horizons_xyzv = nice_Horizons(target, '500@0', times[j], 'smallbody')
        
        # The "state" from our integration
        mpc_xyzv = data[j, 0, :]
        
        # Check whether our integration agrees with Horizons ( within threshold )
        error, good_tf = compare_xyzv(horizons_xyzv,
                                        mpc_xyzv,
                                        1e-7, 1e-8) # MJP 2020-09-03 : Artificially increased thresholds to allow me to make progress while waiting for Holman to debug
                                        #5e-11, 2e-13)  # 7.5m, 30 mm/day
        if np.all(good_tf):
            pass # print('Awesome!')
        else:
            print(f'\n Discrepancy in *is_nbody_output_good_enough()* ...')
            print(f'Time, timestep: {times[j]:}, {j:}')
            print(f'Horizons : {["%18.15e" % _ for _ in horizons_xyzv]}')
            print(f'N-body   : {["%18.15e" % _ for _ in mpc_xyzv]}')
            print(f'Position off by [au]: {error[:3]:}')
            print(f'Velocity off by [au/day]: {error[3:6]:}')
        assert np.all(good_tf)



"""
"""


# Tests of NBody Reader
# -----------------------------------------------------------------------------

# Set up a filepath (file will be created during testing)
test_filepath = os.path.join(os.path.dirname(os.getcwd() ), 'dev_data', '2022AA_demo.txt')


# Actual tests ...
# --------------------------------------------------------------

@pytest.mark.parametrize(('test_filepath'), [test_filepath])
def test_text_file_creation(test_filepath):
    
    # Remove test file if it exists
    if os.path.isfile(test_filepath):
        os.remove(test_filepath)

    # Use convenience func in nbody_reader to create a text file
    nbody_reader.create_nbody_txt(test_filepath)

    # Check that the test file has been created
    assert os.path.isfile(text_filepath)



@pytest.mark.parametrize(('test_filepath'), [test_filepath])
def test_text_file_creation(test_filepath):
    
    # Use convenience func in nbody_reader to create a text file
    nbody_reader.create_nbody_txt(test_filepath)
    
    # Parse the text file
    result = nbody_reader.parse_nbody_txt(test_filepath)

    # Check that the parsed result is as expected
    assert len(result) == 2
    name, a = result
    assert isinstance(name, str), isinstance(a, np.ndarray)
    assert a.shape == (20000,28)



"""




# End

