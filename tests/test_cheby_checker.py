"""
	Tests of the cheby_checker module
"""

# Import third-party packages
# --------------------------------------------------------------
import numpy as np
import sys, os
from astropy_healpix import HEALPix as hp
import pytest

# Import neighboring packages
# --------------------------------------------------------------
sys.path.append(os.path.dirname(os.path.dirname(
    os.path.realpath(__file__))))
from cheby_checker import cheby_checker


# Tests ...
# --------------------------------------------------------------

# --------------------------------------------------------------
# Tests of class instantiation ...
# --------------------------------------------------------------

def test_Base():
    '''
    Test the basic operation of the "Base" class object.
    '''
    
    # Instantiate Base object
    b = cheby_checker.Base()
    
    # Check has the expected variable-attributes
    variables = ['secsPerDay' ,
        'epsilon',
        'HP_nside',
        'HP_order',
        'HPix',
        'HP_npix',
        'sector_length_days',
        'sector_gap',
        'standard_MJDmin',
        'standard_MJDmax',
        'JDlist',
        'db_dir',
        'db_filename',
        ]
    for variable in variables:
        assert hasattr(b,variable)

    # Check has the expected function-attributes
    functions = [   '_fetch_data_directory' ,
                    'map_JD_to_sector_number',
                    'map_sector_number_to_sector_start_JD',
                    'map_JDtimes_to_relative_times',
                    'get_required_sector_dict'
            ]
    for f in functions:
        assert hasattr(b,f)

    # Check variables have reasonable values
    assert b.sector_length_days < 365. # Actually expect == 32, but I guess this may change...
    assert b.secsPerDay == 86400
    assert np.all( b.JDlist == np.arange(b.standard_MJDmin, b.standard_MJDmax+1, 1) )
    assert b.standard_MJDmin < b.standard_MJDmax, f"b.standard_MJDmin={b.standard_MJDmin}, b.standard_MJDmax={b.standard_MJDmax}"
    assert isinstance(b.db_dir , str)
    assert isinstance(b.db_filename , str)

    
    
# --------------------------------------------------------------
# Tests of class functions  ...
# --------------------------------------------------------------

def test_fetch_data_directory():
    '''
    Test the function to fetch the name of the data directory
    NB. This function also *MAKES* the directory if it doesn't exist
    '''
    # Instantiate Base object
    b = cheby_checker.Base()
    
    # Call function
    data_dir = b._fetch_data_directory()
    
    # Check dir exists
    assert os.path.isdir(data_dir)


def test_map_JD_to_sector_number():
    '''
    Test the function to map a date to a sector-number
    '''
    # Instantiate Base object
    b = cheby_checker.Base()

    # Parameters to use in test
    JD0                 = b.standard_MJDmin
    JD_TDB              = np.array( [JD0 , JD0 + b.sector_length_days/2 , JD0 + b.sector_length_days] )
    expected_sectors    = np.array( [0,    0,                             1] )
    
    # Call function
    sectors = b.map_JD_to_sector_number( JD_TDB , JD0)
    
    # Check results
    assert np.all(sectors==expected_sectors )
    

def test_map_sector_number_to_sector_start_JD():
    '''
    Test the function to map a sector number to the starting julian date of that sector
    '''
    # Instantiate Base object
    b = cheby_checker.Base()

    # Parameters to use in test
    JD0          = b.standard_MJDmin
    sectors      = np.array( [0,    1,                             2] )
    expected_JDs = np.array( [JD0 , JD0+b.sector_length_days , JD0+2*b.sector_length_days] )
    
    # Call function
    start_JDs = b.map_sector_number_to_sector_start_JD( sectors , JD0)
    
    # Check results
    assert np.all(start_JDs==expected_JDs ), f"start_JDs={start_JDs}, expected_JDs={expected_JDs}"
    
def test_map_JD_to_sector_number_and_sector_start_JD():
    '''
    Test the function to map a date to a sector-number & a sector-start-JD
    '''
    # Instantiate Base object
    b = cheby_checker.Base()

    # Parameters to use in test
    JD0                 = b.standard_MJDmin
    JD_TDB              = np.array( [JD0 , JD0 + 0.5*b.sector_length_days , JD0 + 1.5*b.sector_length_days] )
    expected_sectors    = np.array( [0,    0,                              1] )
    expected_JDs        = np.array( [JD0 , JD0,                            JD0+b.sector_length_days] )

    # Call function
    sectors, start_JDs = b.map_JD_to_sector_number_and_sector_start_JD( JD_TDB , JD0)
    
    # Check results
    assert np.all(sectors==expected_sectors )
    assert np.all(start_JDs==expected_JDs ), f"start_JDs={start_JDs}, expected_JDs={expected_JDs}"

def test_map_JDtimes_to_relative_times():
    '''
    Test the function that maps a given JD to the relative time within its sector
    I.e. how far the time is from time-zero for that sector
    '''
    # Instantiate Base object
    b = cheby_checker.Base()

    # Parameters to use in test
    JD0                 = b.standard_MJDmin
    JD_TDB              = np.array( [JD0 , JD0 + 0.5*b.sector_length_days , JD0 + 1.25*b.sector_length_days] )
    expected_sectors    = np.array( [0,    0,                              1] )
    expected_times      = np.array( [0,          0.5*b.sector_length_days,        0.25*b.sector_length_days] )
    
    # Call function
    sectors, relative_times = b.map_JDtimes_to_relative_times( JD_TDB)
    
    # Check results
    assert np.all(sectors==expected_sectors )
    assert np.all(relative_times==expected_times )
    
    
def test_get_required_sector_dict():
    '''
    Test the function used to define a complete list of sectors spanning the JDs in JDlist
    Will look like {0: 2440000, 1: 2440032, ..., 750: 2464000}
    '''

    # Instantiate Base object
    b = cheby_checker.Base()

    # Call function
    sector_dict = b.get_required_sector_dict()

    # Check the result looks like {0: 2440000, 1: 2440032, ...}
    assert 0 in sector_dict and sector_dict[0] == b.standard_MJDmin
    assert 1 in sector_dict and sector_dict[1] == b.standard_MJDmin + b.sector_length_days
    
    JD0                 = b.standard_MJDmin
    sector, start_JD = b.map_JD_to_sector_number_and_sector_start_JD( b.standard_MJDmax , JD0)
    assert sector in sector_dict and sector_dict[sector] == start_JD
