"""
    Tests of the coco (coordinate conversion) module

    MJP: 20220128
"""

# Import third-party packages
# --------------------------------------------------------------
import numpy as np
import sys, os
import pytest

# Import neighboring packages
# --------------------------------------------------------------
# cheby_checker/                 # <<-- repo
# cheby_checker/cheby_checker    # <<-- python
# cheby_checker/tests            # <<-- tests
this_dir = os.path.abspath(os.path.dirname( __file__ ))
repo_dir = os.path.abspath(os.path.dirname( this_dir ))
test_dir = os.path.join(repo_dir, 'tests')
code_dir = os.path.join(repo_dir, 'cheby_checker')
for d in [test_dir, code_dir]:
    sys.path.append( d )

import coco
import convenience_Horizons as Horizons
import MPC_library as mpc

'''
     (  # Test 1: Geocenter at 2020-May-28 12:00:00 TDB, equatorial
      '399', 2458998.000000000, 'majorbody'),
     (  # Test 2: 30101 at 2020-Mar-28 12:00:00 TDB, equatorial
      '30101', 2458937.000000000, 'smallbody'),
     (  # Test 3: 30102 at 2020-Mar-28 12:00:00 TDB, equatorial
      '30102', 2458937.000000000, 'smallbody'),
'''
# -- (1) -------------------------------------------------------------------
@pytest.mark.parametrize(
    ('target', 'jd_tdb', 'id_type'),
    [
     ( # Test 0: Geocenter at 2020-Mar-28 12:00:00 TDB, equatorial
      '399', 2458937.000000000, 'majorbody'),
    ])

def test_equatorial_helio2bary(target, jd_tdb, id_type):
    '''
    Test that heliocentric cartesian coordinates taken from Horizons
    are converted to barycentric cartesian and still agree with Horizons.
    
    '''
    # Use horizons to get Helio & Bary versions of the coords
    hor_helio = Horizons.nice_Horizons(target, '500@10', epochs=jd_tdb, id_type=id_type, refplane='earth')#.vectors(refplane='earth')['x', 'y', 'z', 'vx', 'vy', 'vz']
    hor_bary  = Horizons.nice_Horizons(target, '500@0', epochs=jd_tdb, id_type=id_type, refplane='earth')#.vectors(refplane='earth')['x', 'y', 'z', 'vx', 'vy', 'vz']

    input_xyz           = hor_helio.reshape((1,6))
    
    expected_output_xyz = hor_bary.reshape((1,6))

    # Do the transformation
    output_xyz          = coco.equatorial_helio2bary(input_xyz, jd_tdb)


    # Check: Each element should be within 15mm or 15mm/day
    error = np.abs(expected_output_xyz - output_xyz)
    assert np.all(error[:3] < 1e-13), f'XYZ:{error[:3]}'   # XYZ accurate to 15 milli-metres
    assert np.all(error[3:6] < 1e-14),f'UVW:{error[3:6]}'   # V accurate to 1.5 milli-metres/day




# -- (2) -------------------------------------------------------------------
# MA: I'm not really sure whether ecliptic_to_equatorial is supposed to have
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
    hor_table       = Horizons.Horizons(target, centre, epochs=jd_tdb, id_type=id_type)
    hor_in_table    = hor_table.vectors(refplane='ecliptic'
                                        )['x', 'y', 'z', 'vx', 'vy', 'vz']
    hor_out_table   = hor_table.vectors(refplane='earth'
                                        )['x', 'y', 'z', 'vx', 'vy', 'vz']
    input_xyz           = np.atleast_2d( list(hor_in_table.as_array()[0]) )
    expected_output_xyz = np.atleast_2d( list(hor_out_table.as_array()[0]))
 
    # Call the function we want to test
    output_xyz          = coco.ecliptic_to_equatorial(input_xyz)
    
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

names_of_variables   = ('input_ecl_cov',
                        'expected_eq_cov',
                        'comments')
# see 'https://www.visiondummy.com/2014/04/geometric-interpretation-covariance-matrix/'
values_for_each_test = [(np.eye(6),        np.eye(6),        'rotating identity does nothing'),
                        (np.zeros([6, 6]), np.zeros([6, 6]), 'rotating zeros does nothing'),
                        (R6_ecl_to_eq,     R6_ecl_to_eq,     'when input CoV ~ Rotn Matrix'),
                        (R6_eq_to_ecl,     R6_eq_to_ecl,     'when input CoV ~ Rotn Matrix'),
                        ]

@pytest.mark.parametrize(names_of_variables, values_for_each_test[1:])
def test_ecliptic_to_equatorial_covariance(input_ecl_cov, expected_eq_cov, comments):
    '''
    Should do more testing on this to ensure that the CoV is being transformed as desired/expected
    '''
    # The internal workings are such that the CoV is represented as 3D (not 2D) array
    # - E.g. (1,6,6)
    # So the values being tested here need to be forced to be 3D as well ...
    input_ecl_cov   = input_ecl_cov.reshape( (1,6,6) )
    expected_eq_cov = expected_eq_cov.reshape( (1,6,6) )

    # now convert the helio-ecl to bary-eq
    output_eq_cov = coco.ecliptic_to_equatorial(input_ecl_cov)

    # check that the bary CoV has the expected value
    assert np.allclose( expected_eq_cov, output_eq_cov), \
        f"expected_eq_cov={expected_eq_cov}, output_eq_cov={output_eq_cov}"


