# -*- coding: utf-8 -*-
# sifter/tests/test_tracklet

'''
    --------------------------------------------------------------
    tests of sifter's Tracklet class
    
    Jan 2020
    Matt Payne & Mike Alexandersen
    
    --------------------------------------------------------------
    '''



# Import third-party packages
# --------------------------------------------------------------
import sys, os
from astropy import units as u
from astropy_healpix import healpy
from obs80.obs80 import parse80
import numpy as np

# Import neighboring packages
# --------------------------------------------------------------
sys.path.append( os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'sifter') )
import precalc
import sql




def test_instantiation():
    assert isinstance( precalc.Tracklet() , precalc.Tracklet)


def test_parse_observations():
    '''
        Test that precalc.Tracklet.parse_observations is working correctly.
        parse_observations should return:
        integer healpix,
        integer Julian date,
        rate of motion in angle per time units,
        angle of motion relative to RA+ axis (between -180 and 180 deg).
        Rate of motion and angle should be astropy Quantities with units. 
    '''
    T = precalc.Tracklet()  # Create tracklet object

    # define observations  (test pair, snagged from 2011 QF99)
    observation_pair = ['     K11Q99F*~C2011 08 29.52378 01 57 34.729+14 35 44.64         22.8 rc~0qBd568', '     K11Q99F ~C2011 08 29.61470 01 57 34.343+14 35 42.59         22.9 rc~0qBd568']

    # call parse_observations
    # *** AT PRESENT THIS JUST RETURNS RANDOM VALUES
    result = T.parse_observations( observation_pair )
    
    # check that the returned results are formatted as expected
    assert len(result) == 4
    JD, HP, tracklet_name, tracklet_dictionary = result
    assert isinstance(JD, int)
    assert isinstance(HP, np.int64)
    assert isinstance(tracklet_name, str)
    assert len(tracklet_name) == 26
    assert isinstance(tracklet_dictionary, dict)
    assert 'RoM' in tracklet_dictionary.keys()
    assert isinstance(tracklet_dictionary['RoM'], u.quantity.Quantity)
    assert tracklet_dictionary['RoM'].unit.is_equivalent(u.Unit("arcsec / h"))
    assert tracklet_dictionary['RoM'] >= 0
    assert 'AoM' in tracklet_dictionary.keys()
    assert isinstance(tracklet_dictionary['AoM'], u.quantity.Quantity)
    assert tracklet_dictionary['AoM'].unit.is_equivalent(u.deg)
    print(tracklet_dictionary['AoM'])
    assert 0 * u.deg <= tracklet_dictionary['AoM'] < 360 * u.deg
    # The above assertions should be true for any observation pair. 
    # The below assertions are only for the specific example pair. 
    assert JD == 2455803   # Calculated independently
    test_HP = healpy.ang2pix(16, 1.31605272338514, 0.51303437035803, nest=True)
    assert HP == test_HP   # Calculated using ang2pix
    assert HP == 34        # Calculated manually
    assert tracklet_name == 'K11Q99F_2455803.02378__568'


def test_save_tracklet():

    # Set up a Tracklet and use the parse_observations routine to get JD, HP, ...
    T = precalc.Tracklet()
    observation_pair = ['     K11Q99F*~C2011 08 29.52378 01 57 34.729+14 35 44.64         22.8 rc~0qBd568', '     K11Q99F ~C2011 08 29.61470 01 57 34.343+14 35 42.59         22.9 rc~0qBd568']
    JD, HP, tracklet_name, tracklet_dictionary = T.parse_observations( observation_pair )

    # In order to save data, we require sql-db to exist, so let's set that up...
    conn = sql.create_connection( sql.fetch_db_filepath() )
    cur = conn.cursor()
    sql.create_specific_table(conn)

    # Now save the data
    T.save_tracklet(JD, HP, tracklet_name, tracklet_dictionary)

    # Test the data was uploaded and can be downloaded
    cur.execute('SELECT * from tracklets')
    f = cur.fetchone()
    assert( len(f)>3 and f[3] == tracklet_name), 'data not uploaded'

    # Drop the table to facilitate future testing
    cur.execute("DROP TABLE tracklets;")
    cur.execute('SELECT name from sqlite_master WHERE type = "table" AND name = "tracklets"')
    assert cur.fetchone() == None , 'table still exists'

# Call the tests while developing (should really learn how to use pytest... )
test_instantiation()
#test_instantiation_with_observations()
test_parse_observations()
test_save_tracklet()
print('All tests of Tracklet class passed')

