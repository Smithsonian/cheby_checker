# -*- coding: utf-8 -*-
# sifter/tests/test_sqlite

'''
    --------------------------------------------------------------
    tests of sifter's sqlite functions
    
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
import sql




def test_db_creation():

    # Where do we want the db to live
    assert 'sifter' in sql.fetch_db_filepath()

    # Does a db get created
    conn = sql.create_connection( sql.fetch_db_filepath() )
    assert os.path.isfile( os.path.join( sql.fetch_db_filepath() ) )


def test_table_creation():
    expected_table_name = 'tracklets'
    
    conn = sql.create_connection( sql.fetch_db_filepath() )
    cur = conn.cursor()

    # Create the table
    sql.create_specific_table(conn)

    # Should actually implement tests rather than just running ...
    # - get the count of tables with the name
    cur.execute('SELECT name from sqlite_master WHERE type = "table" AND name = "tracklets"')

    # - if the count is 1, then table exists
    assert len(cur.fetchone()) == 1 , 'table does not exist'
    
    # Drop the table to facilitate future testing
    cur.execute("DROP TABLE tracklets;")
    cur.execute('SELECT name from sqlite_master WHERE type = "table" AND name = "tracklets"')
    assert cur.fetchone() == None , 'table still exists'

def test_tracklet_upsert():

    # set up db & table
    conn = sql.create_connection( sql.fetch_db_filepath() )
    cur = conn.cursor()
    sql.create_specific_table(conn)
    
    # create some data and then upload it ...
    jd =123
    hp = 456
    tracklet_name = 'kjhdfasdf'
    tracklet_dict = {'asd':'fgh' , 'ghfgh':987}
    sql.upsert_tracklet(conn, jd, hp, tracklet_name, tracklet_dict)

    # test that the data was actually uploaded
    cur.execute('SELECT * from tracklets')
    f = cur.fetchone()
    assert( len(f)>3 and f[3] == tracklet_name), 'data not uploaded'

    # Drop the table to facilitate future testing
    cur.execute("DROP TABLE tracklets;")
    cur.execute('SELECT name from sqlite_master WHERE type = "table" AND name = "tracklets"')
    assert cur.fetchone() == None , 'table still exists'

def test_tracklet_query():
    
    # set up db & table
    conn = sql.create_connection( sql.fetch_db_filepath() )
    cur = conn.cursor()
    sql.create_specific_table(conn)
    
    # create some data and then upload it ...
    jd =123
    hp = 456
    
    tracklet_name1 = 'kjhdfasdf'
    tracklet_dict1 = {'asd':'fgh' , 'ghfgh':888}
    sql.upsert_tracklet(conn, jd, hp, tracklet_name1, tracklet_dict1)

    tracklet_name2 = 'uituyiu'
    tracklet_dict2 = {'asd':'gdhjdhjdhj' , 'ghfgh':985421}
    sql.upsert_tracklet(conn, jd, hp, tracklet_name2, tracklet_dict2)

    # query the data & check that two dictionaries are returned
    result = sql.query_tracklets_jdhp(conn, jd, hp)
    assert isinstance(result, dict) and len(result) == 2
    assert tracklet_name1 in result and tracklet_name2 in result
    assert isinstance(result[tracklet_name1] , dict)
    assert isinstance(result[tracklet_name2] , dict)

    # Drop the table to facilitate future testing
    cur.execute("DROP TABLE tracklets;")
    cur.execute('SELECT name from sqlite_master WHERE type = "table" AND name = "tracklets"')
    assert cur.fetchone() == None , 'table still exists'

def test_tracklet_delete():
    
    # set up db & table
    conn = sql.create_connection( sql.fetch_db_filepath() )
    cur = conn.cursor()
    sql.create_specific_table(conn)
    
    # create some data and then upload it ...
    jd =123
    hp = 456
    
    tracklet_name1 = 'kjhdfasdf'
    tracklet_dict1 = {'asd':'fgh' , 'ghfgh':888}
    sql.upsert_tracklet(conn, jd, hp, tracklet_name1, tracklet_dict1)
    
    tracklet_name2 = 'uituyiu'
    tracklet_dict2 = {'asd':'gdhjdhjdhj' , 'ghfgh':985421}
    sql.upsert_tracklet(conn, jd, hp, tracklet_name2, tracklet_dict2)
    
    # query the data & check that two dictionaries are returned
    result = sql.query_tracklets_jdhp(conn, jd, hp)
    assert isinstance(result, dict) and len(result) == 2
    assert tracklet_name1 in result and tracklet_name2 in result
    assert isinstance(result[tracklet_name1] , dict)
    assert isinstance(result[tracklet_name2] , dict)
    
    # now delete a tracklet & check that only one dictionary is subsequently returned
    sql.delete_tracklet(conn, tracklet_name1)
    result = sql.query_tracklets_jdhp(conn, jd, hp)
    assert isinstance(result, dict) and len(result) == 1
    assert tracklet_name1 not in result and tracklet_name2 in result
    assert isinstance(result[tracklet_name2] , dict)
    
    # Drop the table to facilitate future testing
    cur.execute("DROP TABLE tracklets;")
    cur.execute('SELECT name from sqlite_master WHERE type = "table" AND name = "tracklets"')
    assert cur.fetchone() == None , 'table still exists'



# Call the tests while developing (should really learn how to use pytest... ) 
test_db_creation()
test_table_creation()
test_tracklet_upsert()
test_tracklet_query()
test_tracklet_delete()
print('All tests of SQL-LITE class passed')
