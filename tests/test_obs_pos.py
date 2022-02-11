"""
	Tests of the obs_pos module

    MJP: 20220128
"""

# Import third-party packages
# --------------------------------------------------------------
import numpy as np
import sys, os
from astropy.time import Time
import pytest

# Import neighboring packages
# --------------------------------------------------------------
from cheby_checker import obs_pos
from cheby_checker.cheby_checker import Base
import convenience_Horizons as Horizons
from cheby_checker import coco


# Files / Directories
# --------------------------------------------------------------
HEAD_DIR = os.path.dirname(os.path.realpath(os.getcwd()))
DATA_DIR = os.path.join(HEAD_DIR, 'dev_data')
#filenames = [os.path.join(DATA_DIR, file)
#              for file in ['30101.eq0_horizons', '30102.eq0_horizons']]


# Tests ...
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
    elif len(error) == 6:
        good_tf = np.abs(error) < np.array([threshold_xyz] * 3 + [threshold_v] * 3)
    else:
        sys.exit( f'len(error)={len(error)} : error={error} : ***EXITING!!!*** ')
    return np.all(good_tf), error
    
    

def test_ObsPos():
    '''
    Test the basic instantiation of the "ObsPos" class object.

    '''
    
    # Instantiate ObsPos object
    op = obs_pos.ObsPos()
    
    # Check has the expected function-attributes
    assert hasattr(op,'get_heliocentric_equatorial_xyz')
    assert hasattr(op,'get_heliocentric_ecliptic_xyz')
    assert hasattr(op,'equatorial_to_ecliptic')
    assert hasattr(op,'check_obsCode')



def test_geocenter():
    '''
    '''
    
    # Instantiate
    op = obs_pos.ObsPos()

    # Query Horizons for heliocentric_equatorial posn of the geocenter
    target  = '399'       # Earth Geocenter
    centre  = '500@10'    # Heliocentre
    epochs  = '2458850.0'
    timescale = 'tdb'
    id_type = 'majorbody'

    refplane='earth'
    horizons_helio_eq = Horizons.nice_Horizons(target, centre, epochs, id_type, refplane=refplane )
    refplane='ecliptic'
    horizons_helio_ec = Horizons.nice_Horizons(target, centre, epochs, id_type, refplane=refplane )
    

    # Astropy-Time object (for easy conversion)
    t = Time( epochs , format='jd', scale=timescale )
    t = t.utc  # <<-- This converts to utc (from tdb, above)

    # Call the function(s) to be tested
    helio_eq = op.get_heliocentric_equatorial_xyz( t.jd , obsCode='500' )
    helio_ec = op.get_heliocentric_ecliptic_xyz( t.jd , obsCode='500' )


    # Perform the checks ...
    # These both pass at the 1e-11 ( ~1.5m ) level, but not at 1e-12
    similar_bool, err = similar_xyzuvw(horizons_helio_eq[:3], helio_eq, threshold_xyz=1e-11, threshold_v=1e-11)
    assert similar_bool, f'NOT SIMILAR!:\n\t err={err}, horizons_helio_eq[:3]={horizons_helio_eq[:3]},helio_eq={helio_eq}'
    similar_bool, err = similar_xyzuvw(horizons_helio_ec[:3], helio_ec, threshold_xyz=1e-11, threshold_v=1e-11)
    assert similar_bool, f'NOT SIMILAR!:\n\t err={err}, horizons_helio_ec[:3]={horizons_helio_ec[:3]},helio_ec={helio_ec}'
    print(err)


def test_get_heliocentric_equatorial_xyz_A():
    '''
    Test the ObsPos.get_heliocentric_equatorial_xyz() function
    
    Here we examine the HELIOCENTRIC EQUATORIAL position of the F51 observatory

    '''
    
    # Instantiate ObsPos object
    op = obs_pos.ObsPos()
    
    # Query Horizons for heliocentric_equatorial posn of F51
    # NB we use a hack, setting the target as the Sun, and the center as the observatory)
    # - as such, we need to multiply the values by -1
    target  = '10'       # Heliocentre
    centre  = 'F51'      # PanSTARRS
    epochs  = '2458850.0'
    timescale = 'tdb'
    id_type = 'majorbody'
    refplane='earth'
    
    horizons_state = Horizons.nice_Horizons(target, centre, epochs, id_type, refplane=refplane )
    horizons_state = -1*horizons_state
    
    # Astropy-Time object (for easy conversion)
    t = Time( epochs , format='jd', scale=timescale )
    t = t.utc  # <<-- This converts to utc (from tdb, above)

    # Call the function to be tested
    result = op.get_heliocentric_equatorial_xyz( t.jd , obsCode=centre )
    
    # Compare the results from get_heliocentric_equatorial_xyz with Horizons ...
    # *** THESE ONLY PASS AT THE 1e-9 LEVEL, NOT AT 1e-10 or 1e-11 ***
    assert np.allclose(horizons_state[:3], result, rtol=1e-9, atol=1e-9)

def test_get_heliocentric_equatorial_xyz_B():
    '''
    Multiple tests of the ObsPos.get_heliocentric_equatorial_xyz() function
    '''
    
    # Instantiate DB object
    op = obs_pos.ObsPos()
    
    # Query Horizons for heliocentric_equatorial posn of F51
    # NB we use a hack, setting the target as the Sun, and the center as the observatory)
    # - as such, we need to multiply the values by -1
    target  = '10'       # Heliocentre
    centres  = ['F51','I41','500', '000', '703']      # Various obsCodes
    epochs  = '2458850.0'
    timescale = 'tdb'
    id_type = 'majorbody'
    refplane='earth'

    # Astropy-Time object (for easy conversion)
    t = Time( epochs , format='jd', scale=timescale )
    t = t.utc  # <<-- This converts to utc (from tdb, above)

    for centre in centres:
        horizons_state = Horizons.nice_Horizons(target, centre, epochs, id_type, refplane=refplane )
        horizons_state = -1*horizons_state
    
        # Call the function to be tested
        result = op.get_heliocentric_equatorial_xyz( t.jd , obsCode=centre )
    
        # Compare the results from get_heliocentric_equatorial_xyz with Horizons ...
        # *** THESE ONLY PASS AT THE 1e-9 LEVEL, NOT AT 1e-10 or 1e-11 ***
        assert np.allclose(horizons_state[:3], result, rtol=1e-9, atol=1e-9)


def test_get_heliocentric_ecliptic_xyz_A():
    '''
    Test the ObsPos.get_heliocentric_ecliptic_xyz() function
    
    Here we examine the HELIOCENTRIC ECLIPTIC position of the F51 observatory

    '''
    
    # Instantiate DB object
    op = obs_pos.ObsPos()
    
    # Query Horizons for heliocentric_ECLIPTTIC posn of F51
    # NB we use a hack, setting the target as the Sun, and the center as the observatory)
    # - as such, we need to multiply the values by -1
    target  = '10'       # <<-- Asteroid number 54321 == 2000 JA81
    centre  = 'F51'      # Heliocentric
    epochs  = '2458850.0'
    timescale = 'tdb'
    id_type = 'majorbody'
    refplane='ecliptic'
    
    horizons_state = Horizons.nice_Horizons(target, centre, epochs, id_type, refplane=refplane )
    horizons_state = -1*horizons_state
    
    # Astropy-Time object (for easy conversion)
    t = Time( epochs , format='jd', scale=timescale )
    t = t.utc  # <<-- This converts to utc (from tdb, above)

    # Call the function to be tested
    result = op.get_heliocentric_ecliptic_xyz( t.jd , obsCode=centre )
    
    # Compare the results from get_heliocentric_equatorial_xyz with Horizons ...
    # *** THESE ONLY PASS AT THE 1e-9 LEVEL, NOT AT 1e-10 or 1e-11 ***
    assert np.allclose(horizons_state[:3], result, rtol=1e-9, atol=1e-9)


def test_get_heliocentric_ecliptic_xyz_B():
    '''
    Multiple tests of the ObsPos.get_heliocentric_ecliptic_xyz() function
    '''
    
    # Instantiate DB object
    op = obs_pos.ObsPos()
    
    # Query Horizons for heliocentric_equatorial posn of F51
    # NB we use a hack, setting the target as the Sun, and the center as the observatory)
    # - as such, we need to multiply the values by -1
    target  = '10'       # Heliocentre
    centres  = ['F51','I41','500', '000', '703']      # Various obsCodes
    epochs  = '2458850.0'
    timescale = 'tdb'
    id_type = 'majorbody'
    refplane='ecliptic'

    # Astropy-Time object (for easy conversion)
    t = Time( epochs , format='jd', scale=timescale )
    t = t.utc  # <<-- This converts to utc (from tdb, above)

    for centre in centres:
        horizons_state = Horizons.nice_Horizons(target, centre, epochs, id_type, refplane=refplane )
        horizons_state = -1*horizons_state
    
        # Call the function to be tested
        result = op.get_heliocentric_ecliptic_xyz( t.jd , obsCode=centre )
    
        # Compare the results from get_heliocentric_equatorial_xyz with Horizons ...
        # *** THESE ONLY PASS AT THE 1e-9 LEVEL, NOT AT 1e-10 or 1e-11 ***
        assert np.allclose(horizons_state[:3], result, rtol=1e-9, atol=1e-9)

    

def test_get_barycentric_A():
    '''
    Test the ObsPos.get_heliocentric_ecliptic_xyz()
    and the  coco.ecliptic_to_equatorial(input_xyz) functions
    
    Here we examine the BARYCENTRIC positions (ECLIPTIC & EQUATORIAL) of the F51 observatory
    
    I am ensuring that the MPC's positions are consistent with those from JPL

    '''
    
    # Instantiate DB object
    op = obs_pos.ObsPos()
    
    # Query Horizons for heliocentric posn of F51
    # NB we use a hack, setting the target as the Sun, and the center as the observatory)
    # - as such, we need to multiply the values by -1
    target  = '10'       #
    centre  = 'F51'      # Heliocentric
    epochs  = '2458853.0'#2458850.0'
    timescale = 'tdb'
    id_type = 'majorbody'
    
    refplane='ecliptic'
    horizons_helio_eclip = Horizons.nice_Horizons(target, centre, epochs, id_type, refplane=refplane )
    horizons_helio_eclip = -1*horizons_helio_eclip
    #print('horizons_helio_eclip=', horizons_helio_eclip)
    
    refplane='earth'
    horizons_helio_eq = Horizons.nice_Horizons(target, centre, epochs, id_type, refplane=refplane )
    horizons_helio_eq = -1*horizons_helio_eq
    #print('horizons_helio_eq=', horizons_helio_eq)

    # Query Horizons for barycentric posn of F51
    # Won't work directly, so have to put this hack in place ...
    # ... get posn of heliocentre w.r.t. barycenter, and then add to posn of F51 w.r.t. heliocenter
    refplane='ecliptic'
    b2h = Horizons.nice_Horizons(target, '@0', epochs, id_type, refplane=refplane )
    horizons_bary_eclip = b2h + horizons_helio_eclip
    #print('horizons_bary_eclip=', horizons_bary_eclip)

    refplane='earth'
    b2h = Horizons.nice_Horizons(target, '@0', epochs, id_type, refplane=refplane )
    horizons_bary_eq = b2h + horizons_helio_eq
    #print('horizons_bary_eq=', horizons_bary_eq)

    # Astropy-Time object (for easy conversion)
    t = Time( epochs , format='jd', scale=timescale )
    t = t.utc  # <<-- This converts to utc (from tdb, above)

    # Call the functions to be tested
    helio_eclip = op.get_heliocentric_ecliptic_xyz( t.jd , obsCode=centre )
    helio_eq    = coco.ecliptic_to_equatorial( np.atleast_2d(helio_eclip) )
    helio_eq2   = op.get_heliocentric_equatorial_xyz( t.jd , obsCode=centre )
    bary_eq     = coco.equatorial_helio2bary( helio_eq , float(epochs) )
    bary_eclip  = coco.ecliptic_to_equatorial( np.atleast_2d(bary_eq) , backwards=True )

    # Test the results ...
    #
    # (1) Compare the helio_eclip with Horizons ...
    # These pass at the 1e-10 ( ~15m ) level, but not at 1e-11
    similar_bool, err = similar_xyzuvw(horizons_helio_eclip[:3], helio_eclip, threshold_xyz=1e-10, threshold_v=1e-10)
    assert similar_bool, f'NOT SIMILAR!:\n\t err={err}'
    
    # (2) helio_eq
    # These pass at the 1e-10 ( ~15m ) level, but not at 1e-11
    similar_bool, err = similar_xyzuvw(horizons_helio_eq[:3], helio_eq[0], threshold_xyz=1e-10, threshold_v=1e-10)
    assert similar_bool, f'NOT SIMILAR!:\n\t err={err}, horizons_helio_eq[:3]={horizons_helio_eq[:3]},helio_eq={helio_eq}'
    similar_bool, err = similar_xyzuvw(horizons_helio_eq[:3], helio_eq2, threshold_xyz=1e-10, threshold_v=1e-10)
    assert similar_bool, f'NOT SIMILAR!:\n\t err={err}, horizons_helio_eq[:3]={horizons_helio_eq[:3]},helio_eq2={helio_eq2}'

    # (3) bary_eq
    # These pass at the 1e-10 ( ~15m ) level, but not at 1e-11
    similar_bool, err = similar_xyzuvw(horizons_bary_eq[:3], bary_eq[0], threshold_xyz=1e-10, threshold_v=1e-10)
    assert similar_bool, f'NOT SIMILAR!:\n\t err={err}, horizons_bary_eq[:3]={horizons_bary_eq[:3]},bary_eq={bary_eq}'


    # (4) bary_eclip
    # These pass at the 1e-10 ( ~15m ) level, but not at 1e-11
    similar_bool, err = similar_xyzuvw(horizons_bary_eclip[:3], bary_eclip[0], threshold_xyz=1e-10, threshold_v=1e-10)
    assert similar_bool, f'NOT SIMILAR!:\n\t err={err}, horizons_bary_eclip[:3]={horizons_bary_eclip[:3]},bary_eclip={bary_eclip}'

    print("Passed test_get_barycentric_A")
    
test_geocenter()
