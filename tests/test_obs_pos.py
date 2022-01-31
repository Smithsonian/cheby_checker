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
# cheby_checker/                 # <<-- repo
# cheby_checker/cheby_checker    # <<-- python
# cheby_checker/tests            # <<-- tests
this_dir = os.path.abspath(os.path.dirname( __file__ ))
repo_dir = os.path.abspath(os.path.dirname( this_dir ))
test_dir = os.path.join(repo_dir, 'tests')
code_dir = os.path.join(repo_dir, 'cheby_checker')
for d in [test_dir, code_dir]:
    sys.path.append( d )
import obs_pos
from cheby_checker import Base

import convenience_Horizons as Horizons


# Files / Directories
# --------------------------------------------------------------
HEAD_DIR = os.path.dirname(os.path.realpath(os.getcwd()))
DATA_DIR = os.path.join(HEAD_DIR, 'dev_data')
#filenames = [os.path.join(DATA_DIR, file)
#              for file in ['30101.eq0_horizons', '30102.eq0_horizons']]


# Tests ...
# --------------------------------------------------------------




# --------------------------------------------------------------
# Tests of class instantiation ...
# --------------------------------------------------------------

def test_ObsPos():
    '''
    Test the basic instantiation of the "ObsPos" class object.

    '''
    
    # Instantiate DB object
    op = obs_pos.ObsPos()
    
    # Check has the expected function-attributes
    assert hasattr(op,'get_heliocentric_equatorial_xyz')
    assert hasattr(op,'get_heliocentric_ecliptic_xyz')
    assert hasattr(op,'equatorial_to_ecliptic')
    assert hasattr(op,'check_obsCode')





def test_get_heliocentric_equatorial_xyz_A():
    '''
    Test the ObsPos.get_heliocentric_equatorial_xyz() function
    
    Here we examine the HELIOCENTRIC EQUATORIAL position of the F51 observatory

    '''
    
    # Instantiate DB object
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
    Multiple tests of the ObsPos.get_heliocentric_eclipttic_xyz() function
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

    
#test_get_heliocentric_equatorial_xyz()
#test_get_heliocentric_ecliptic_xyz()
