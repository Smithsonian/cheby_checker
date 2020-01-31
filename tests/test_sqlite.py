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

    if os.path.isfile( sql.fetch_db_filepath() ):
        os.remove(sql.fetch_db_filepath())

    # Does a db get created
    conn = sql.create_connection( sql.fetch_db_filepath() )
    assert os.path.isfile( os.path.join( sql.fetch_db_filepath() ) )


def test_table_creation():
    expected_table_name = 'tracklets'
    
    # set up db & table
    if os.path.isfile( sql.fetch_db_filepath() ):
        os.remove(sql.fetch_db_filepath())
    conn = sql.create_connection( sql.fetch_db_filepath() )
    cur = conn.cursor()

    # Create the table
    sql.create_specific_table(conn)

    # - get the count of tables with the name
    cur.execute('SELECT name from sqlite_master WHERE type = "table" AND name = "tracklets"')

    # - if the count is 1, then table exists
    assert len(cur.fetchone()) == 1 , 'table does not exist'
    
    # Delete the db to facilitate future testing
    os.remove(sql.fetch_db_filepath())

def test_tracklet_upsert():

    # set up db & table
    if os.path.isfile( sql.fetch_db_filepath() ):
        os.remove(sql.fetch_db_filepath())
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

    # Delete the db to facilitate future testing
    os.remove(sql.fetch_db_filepath())


def test_tracklets_upsert():
    """ Here we are updating/inserting **lists** of tracklet data """
    
    # set up db & table
    if os.path.isfile( sql.fetch_db_filepath() ):
        os.remove(sql.fetch_db_filepath())
    conn = sql.create_connection( sql.fetch_db_filepath() )
    cur = conn.cursor()
    sql.create_specific_table(conn)
    
    # create some data and then upload it ...
    jd_list = [123, 234]
    hp_list = [456, 567]
    tracklet_name_list = ['kjhdfasdf', 'iuyeruy']
    tracklet_dict_list = [{'asd':'fgh' , 'ghfgh':987} , {'asd':'klhj' , 'ghfgh':4563} ]
    sql.upsert_tracklets(conn, jd_list, hp_list, tracklet_name_list, tracklet_dict_list)
    
    # test that the data was actually uploaded
    cur.execute('SELECT * from tracklets')
    f = cur.fetchall()
    assert( len(f) == 2 ) , 'data not uploaded'
    for i in range(len(f)):
        assert f[i][3] == tracklet_name_list[i], 'data not uploaded'


    # Delete the db to facilitate future testing
    os.remove(sql.fetch_db_filepath())




def test_tracklet_query():
    
    # set up db & table
    if os.path.isfile( sql.fetch_db_filepath() ):
        os.remove(sql.fetch_db_filepath())
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
    list_of_tuples = sql.query_tracklets_jdhp(conn, jd, hp)
    assert isinstance(list_of_tuples, list) and len(list_of_tuples) == 2
    for tup in list_of_tuples:
        assert isinstance( tup, (tuple, list,))
        assert tup[0] in [tracklet_name1, tracklet_name2]
        assert tup[1] in [tracklet_dict1, tracklet_dict2]
 
    # Delete the db to facilitate future testing
    os.remove(sql.fetch_db_filepath())

def test_delete_tracklet():
    
    # set up db & table
    if os.path.isfile( sql.fetch_db_filepath() ):
        os.remove(sql.fetch_db_filepath())
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
    list_of_tuples = sql.query_tracklets_jdhp(conn, jd, hp)
    assert isinstance(list_of_tuples, list) and len(list_of_tuples) == 2
    for tup in list_of_tuples:
        assert isinstance( tup, (tuple, list,))
        assert tup[0] in [tracklet_name1, tracklet_name2]
        assert tup[1] in [tracklet_dict1, tracklet_dict2]
    
    # now delete a tracklet & check that only one dictionary is subsequently returned
    sql.delete_tracklet(conn, tracklet_name1)
    list_of_tuples = sql.query_tracklets_jdhp(conn, jd, hp)
    assert isinstance(list_of_tuples, list) and len(list_of_tuples) == 1
    for tup in list_of_tuples:
        assert isinstance( tup, (tuple, list,))
        assert tup[0] in [ tracklet_name2]
        assert tup[1] in [ tracklet_dict2]
    
    # Delete the db to facilitate future testing
    os.remove(sql.fetch_db_filepath())

def test_delete_tracklets():
    
    # set up db & table
    if os.path.isfile( sql.fetch_db_filepath() ):
        os.remove(sql.fetch_db_filepath())
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

    tracklet_name3 = 'sdfhsdfgh'
    tracklet_dict3 = {'asd':'zbczbzfb' , 'ghfgh':34563}
    sql.upsert_tracklet(conn, jd, hp, tracklet_name3, tracklet_dict3)

    tracklet_name4 = 'oiippipi'
    tracklet_dict4 = {'asd':'aeraer' , 'ghfgh':56856}
    sql.upsert_tracklet(conn, jd, hp, tracklet_name4, tracklet_dict4)


    # query the data & check that 4 dictionaries are returned
    list_of_tuples = sql.query_tracklets_jdhp(conn, jd, hp)
    assert isinstance(list_of_tuples, list) and len(list_of_tuples) == 4
    for tup in list_of_tuples:
        assert isinstance( tup, (tuple, list,))
        assert tup[0] in [tracklet_name1, tracklet_name2, tracklet_name3, tracklet_name4]
        assert tup[1] in [tracklet_dict1, tracklet_dict2, tracklet_dict3, tracklet_dict4]
    
    # now delete two tracklets & check that two dictionaries remain
    sql.delete_tracklets(conn, [tracklet_name1 , tracklet_name3 ] )
    list_of_tuples = sql.query_tracklets_jdhp(conn, jd, hp)
    assert isinstance(list_of_tuples, list) and len(list_of_tuples) == 2
    for tup in list_of_tuples:
        assert isinstance( tup, (tuple, list,))
        assert tup[0] in [ tracklet_name2,  tracklet_name4]
        assert tup[1] in [ tracklet_dict2,  tracklet_dict4]

    # Delete the db to facilitate future testing
    os.remove(sql.fetch_db_filepath())




# Call the tests while developing (should really learn how to use pytest... ) 
test_db_creation()
test_table_creation()
test_tracklet_upsert()
test_tracklets_upsert()
test_tracklet_query()
test_delete_tracklet()
test_delete_tracklets()
print('All tests of SQL-LITE class passed')
