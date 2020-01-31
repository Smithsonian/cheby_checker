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
import numpy as np 

# Import neighboring packages
# --------------------------------------------------------------
sys.path.append( os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'sifter') )
import precalc
import sql



def convenience_func_create_db_and_tables():
    
    # In order to save data, we require sql-db to exist, so let's set that up...
    # Force deletion then creation of db...
    if os.path.isfile( sql.fetch_db_filepath() ):
        os.remove(sql.fetch_db_filepath())
    conn = sql.create_connection( sql.fetch_db_filepath() )
    cur = conn.cursor()

    # Create required table(s)
    sql.create_specific_table(conn)
    
    # Double-check that this worked by getting the count of tables with the name
    # - if the count is 1, then table exists
    cur.execute('SELECT name from sqlite_master WHERE type = "table" AND name = "tracklets"')
    res =cur.fetchone()
    assert len(res) == 1 , 'table does not exist'
    conn.close()



def test_instantiation():
    assert isinstance( precalc.Tracklets() , precalc.Tracklets )


def test_parse_observations():

    T = precalc.Tracklets()
    
    # define observations
    # *** AT PRESENT THESE ARE JUST DUMMY/BLANK OBS ***
    observation_pair = [ [],[] ]

    # call parse_observations
    # *** AT PRESENT THIS JUST RETURNS RANDOM VALUES
    result = T.parse_observation_lists( observation_pair )
    
    # check that the returned results are as expected
    assert len(result) == 4
    JD_list, HP_list, tracklet_name_list, tracklet_dictionary_list = result
    assert isinstance(JD_list, list)
    assert isinstance(HP_list, list)
    assert isinstance(tracklet_name_list, list)
    assert isinstance(tracklet_dictionary_list, list)
    assert len(JD_list) == len(HP_list) == len(tracklet_name_list) == len(tracklet_dictionary_list)
    for JD, HP, tracklet_name, tracklet_dictionary in zip(JD_list, HP_list, tracklet_name_list, tracklet_dictionary_list):
        assert isinstance(JD, int)
        assert isinstance(HP, int)
        assert isinstance(tracklet_name, str)
        assert isinstance(tracklet_dictionary, dict)


def test_save_tracklets():

    # Create db from scratch
    convenience_func_create_db_and_tables()
    

    # Set up a Tracklet and use the parse_observations routine to get JD, HP, ...
    T = precalc.Tracklets()
    observation_pair = [ [],[] ]
    JD_list, HP_list, tracklet_name_list, tracklet_dictionary_list = T.parse_observation_lists( observation_pair )

    # Now save the data in the db
    T.save_tracklets(JD_list, HP_list, tracklet_name_list, tracklet_dictionary_list)

    # Test the data was uploaded and can be downloaded
    cur = T.conn.cursor()
    cur.execute('SELECT * from tracklets')
    f = cur.fetchall()
    assert( len(f)==2 and np.all([ len(_)>3 for _ in f]) ), 'data not uploaded'

    # Completely delete db to facilitate future testing
    os.remove(sql.fetch_db_filepath())

def test_instantiation_with_observations():
    
    # Create db from scratch
    convenience_func_create_db_and_tables()
    
    # define observations
    # *** AT PRESENT THESE ARE JUST DUMMY/BLANK OBS ***
    observation_pair = [ [],[] ]
    
    # instantiate with observation_pair
    T = precalc.Tracklets( observation_pair )
    
    # test that the above caused the tracklet to be uploaded to db
    cur = T.conn.cursor()
    cur.execute('SELECT * from tracklets')
    f = cur.fetchone()
    assert( len(f)>3 ), 'data not uploaded'

    # Completely delete db to facilitate future testing
    os.remove(sql.fetch_db_filepath())



# Call the tests while developing (should really learn how to use pytest... )
test_instantiation()
test_parse_observations()
test_save_tracklets()
#test_instantiation_with_observations()
print('All tests of Tracklets* class passed')

