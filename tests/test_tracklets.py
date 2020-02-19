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
    observation_pairs= [[ '     K11Q99F*~C2011 08 29.52378 01 57 34.729+14 35 44.64         22.8 rc~0qBd568',
                          '     K11Q99F ~C2011 08 29.61470 01 57 34.343+14 35 42.59         22.9 rc~0qBd568'],
                        [ '     K11Q99F*~C2012 08 29.52378 01 57 34.729+14 35 44.64         22.8 rc~0qBd568',
                          '     K11Q99F ~C2012 08 29.61470 01 57 34.343+14 35 42.59         22.9 rc~0qBd568']]

    # call parse_observations
    # *** AT PRESENT THIS JUST RETURNS RANDOM VALUES
    tracklet_dictionary_list = T.parse_observation_lists( observation_pairs )

    # check that the returned results are as expected
    assert isinstance(tracklet_dictionary_list, list)
    assert len(observation_pairs) == len(tracklet_dictionary_list)
    for tracklet_dictionary in tracklet_dictionary_list:
        assert 'JD' in tracklet_dictionary
        assert 'HP' in tracklet_dictionary
        assert 'tracklet_name' in tracklet_dictionary


def test_save_tracklets():

    # Create db from scratch
    convenience_func_create_db_and_tables()
    

    # Set up a Tracklet and use the parse_observations routine to get JD, HP, ...
    T = precalc.Tracklets()
    observation_pairs= [[ '     K11Q99F*~C2011 08 29.52378 01 57 34.729+14 35 44.64         22.8 rc~0qBd568',
                         '     K11Q99F ~C2011 08 29.61470 01 57 34.343+14 35 42.59         22.9 rc~0qBd568'],
                        [ '     K11Q99F*~C2012 08 29.52378 01 57 34.729+14 35 44.64         22.8 rc~0qBd568',
                         '     K11Q99F ~C2012 08 29.61470 01 57 34.343+14 35 42.59         22.9 rc~0qBd568']]
    tracklet_dictionary_list = T.parse_observation_lists( observation_pairs )

    # Now save the data in the db
    T.save_tracklets(tracklet_dictionary_list)

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
    observation_pairs= [[ '     K11Q99F*~C2011 08 29.52378 01 57 34.729+14 35 44.64         22.8 rc~0qBd568',
                         '     K11Q99F ~C2011 08 29.61470 01 57 34.343+14 35 42.59         22.9 rc~0qBd568'],
                        [ '     K11Q99F*~C2012 08 29.52378 01 57 34.729+14 35 44.64         22.8 rc~0qBd568',
                         '     K11Q99F ~C2012 08 29.61470 01 57 34.343+14 35 42.59         22.9 rc~0qBd568']]
    
    # instantiate with observation_pair
    T = precalc.Tracklets( observation_pairs )
    
    # test that the above caused the tracklet to be uploaded to db
    cur = T.conn.cursor()
    cur.execute('SELECT * from tracklets')
    f = cur.fetchone()
    assert( len(f)>3 ), 'data not uploaded'

    # Completely delete db to facilitate future testing
    os.remove(sql.fetch_db_filepath())



