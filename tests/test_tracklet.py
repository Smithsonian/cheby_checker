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

# Import neighboring packages
# --------------------------------------------------------------
sys.path.append( os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'sifter') )
import precalc
import sql




def test_instantiation():
    assert isinstance( precalc.Tracklet() , precalc.Tracklet)


def test_parse_observations():

    T = precalc.Tracklet()
    
    # define observations
    # *** AT PRESENT THESE ARE JUST DUMMY OBS ***
    observations = []

    # call parse_observations
    # *** AT PRESENT THIS JUST RETURNS RANDOM VALUES
    result = T.parse_observations( observations )
    
    # check that the returned results are as expected
    assert len(result) == 4
    JD, HP, tracklet_name, tracklet_dictionary = result
    assert isinstance(JD, int)
    assert isinstance(HP, int)
    assert isinstance(tracklet_name, str)
    assert isinstance(tracklet_dictionary, dict)

def test_save_tracklet():

    # Set up a Tracklet and use the parse_observations routine to get JD, HP, ...
    T = precalc.Tracklet()
    observations = []
    JD, HP, tracklet_name, tracklet_dictionary = T.parse_observations( observations )

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


