# -*- coding: utf-8 -*-
# cheby_checker/tests/test_sqlite

'''
    --------------------------------------------------------------
    tests of sifter's sqlite functions

    Jan 2020
    Matt Payne & Mike Alexandersen

    --------------------------------------------------------------
    '''


# Import third-party packages
# --------------------------------------------------------------
import sys
import os
import pytest


# Import neighboring packages
# --------------------------------------------------------------
from cheby_checker import sifter_sql as sql


# Convenience data / functions to aid testing
# --------------------------------------------------------------
test_tracklet_dict_list = []
for i in range(4):
    test_tracklet_dict_list.append({'JD': 123 + i,
                                    'HP': 456 + i,
                                    'tracklet_name': 'kjhdfasdf' + str(i),
                                    'asd': 'fgh',
                                    'ghfgh': 987
                                    }
                                   )


# Actual tests ...
# --------------------------------------------------------------
def test_db_creation():
    '''Test that an empty database can be created.'''

    # Where do we want the db to live
    assert 'sifter' in sql.fetch_db_filepath()

    if os.path.isfile(sql.fetch_db_filepath()):
        os.remove(sql.fetch_db_filepath())

    # Does a db get created
    conn = sql.create_connection(sql.fetch_db_filepath())
    assert os.path.isfile(os.path.join(sql.fetch_db_filepath()))


def test_table_creation():
    '''Test table creation.'''
    expected_table_name = 'tracklets'

    # set up db & table
    if os.path.isfile(sql.fetch_db_filepath()):
        os.remove(sql.fetch_db_filepath())
    conn = sql.create_connection(sql.fetch_db_filepath())
    cur = conn.cursor()

    # Create the table
    sql.create_specific_table(conn)

    # - get the count of tables with the name
    cur.execute('SELECT name from sqlite_master WHERE type = "table" AND name = "tracklets"')

    # - if the count is 1, then table exists
    assert len(cur.fetchone()) == 1, 'table does not exist'

    # Delete the db to facilitate future testing
    os.remove(sql.fetch_db_filepath())


@pytest.mark.parametrize('tracklet_dict_list', [test_tracklet_dict_list])
def test_tracklet_upsert(tracklet_dict_list):
    '''Test tracklet upsertion into the database.'''
    # set up db & table
    if os.path.isfile(sql.fetch_db_filepath()):
        os.remove(sql.fetch_db_filepath())
    conn = sql.create_connection(sql.fetch_db_filepath())
    cur = conn.cursor()
    sql.create_specific_table(conn)

    # create some data and then upload it ...
    tracklet_dict = tracklet_dict_list[0]
    sql.upsert_tracklet(conn, tracklet_dict['JD'], tracklet_dict['HP'],
                        tracklet_dict['tracklet_name'], tracklet_dict)

    # test that the data was actually uploaded
    cur.execute('SELECT * from tracklets')
    f = cur.fetchone()
    assert (len(f) > 3 and f[3] == tracklet_dict['tracklet_name']),\
        'data not uploaded'

    # Delete the db to facilitate future testing
    os.remove(sql.fetch_db_filepath())


@pytest.mark.parametrize(('tracklet_dict_list'), [test_tracklet_dict_list])
def test_tracklets_upsert(tracklet_dict_list):
    """ Here we are updating/inserting **lists** of tracklet data """

    # set up db & table
    if os.path.isfile(sql.fetch_db_filepath()):
        os.remove(sql.fetch_db_filepath())
    conn = sql.create_connection(sql.fetch_db_filepath())
    cur = conn.cursor()
    sql.create_specific_table(conn)

    # upload data ...
    JD = [tracklet_dic['JD'] for tracklet_dic in tracklet_dict_list]
    HP = [tracklet_dic['HP'] for tracklet_dic in tracklet_dict_list]
    tracklet_name = [tracklet_dic['tracklet_name']
                     for tracklet_dic in tracklet_dict_list]

    sql.upsert_tracklets(conn, JD, HP, tracklet_name, tracklet_dict_list)

    # test that the data was actually uploaded
    cur.execute('SELECT * from tracklets')
    f = cur.fetchall()
    assert(len(f) == len(tracklet_dict_list)), 'data not uploaded'
    for ii, fi in enumerate(f):
        assert fi[3] == tracklet_name[ii], 'data not uploaded'

    # Delete the db to facilitate future testing
    os.remove(sql.fetch_db_filepath())


@pytest.mark.parametrize(('tracklet_dict_list'), [test_tracklet_dict_list])
def test_tracklet_query(tracklet_dict_list):
    '''Test querying a tracklet.'''
    # set up db & table
    if os.path.isfile(sql.fetch_db_filepath()):
        os.remove(sql.fetch_db_filepath())
    conn = sql.create_connection(sql.fetch_db_filepath())
    cur = conn.cursor()
    sql.create_specific_table(conn)

    # upload data ...
    JD = [tracklet_dic['JD'] for tracklet_dic in tracklet_dict_list]
    HP = [tracklet_dic['HP'] for tracklet_dic in tracklet_dict_list]
    tracklet_name = [tracklet_dic['tracklet_name']
                     for tracklet_dic in tracklet_dict_list]
    sql.upsert_tracklets(conn, JD, HP, tracklet_name, tracklet_dict_list)

    # query the data & check that requisite dictionaries are returned
    list_of_tuples = sql.query_tracklets_jdhp(conn, JD[0], HP[0])
    assert isinstance(list_of_tuples, list) and len(list_of_tuples) == 1

    # Delete the db to facilitate future testing
    os.remove(sql.fetch_db_filepath())


@pytest.mark.parametrize(('tracklet_dict_list'), [test_tracklet_dict_list])
def test_tracklet_query_mutiple_HP(tracklet_dict_list):
    '''Test querying multiple Heal Pix.'''
    # set up db & table
    if os.path.isfile(sql.fetch_db_filepath()):
        os.remove(sql.fetch_db_filepath())
    conn = sql.create_connection(sql.fetch_db_filepath())
    cur = conn.cursor()
    sql.create_specific_table(conn)

    # upload data ...
    JD = [tracklet_dic['JD'] for tracklet_dic in tracklet_dict_list]
    HP = [tracklet_dic['HP'] for tracklet_dic in tracklet_dict_list]
    tracklet_name = [tracklet_dic['tracklet_name']
                     for tracklet_dic in tracklet_dict_list]
    sql.upsert_tracklets(conn, JD, HP, tracklet_name, tracklet_dict_list)

    # query the data & check that requisite dictionaries are returned
    list_of_tuples = sql.query_tracklets_jd_hplist(conn, JD[0], HP)
    assert isinstance(list_of_tuples, list) and len(list_of_tuples) == 1

    # Delete the db to facilitate future testing
    os.remove(sql.fetch_db_filepath())


@pytest.mark.parametrize(('tracklet_dict_list'), [test_tracklet_dict_list])
def test_delete_tracklet(tracklet_dict_list):
    '''Test deletion of a tracklet.'''
    # set up db & table
    if os.path.isfile(sql.fetch_db_filepath()):
        os.remove(sql.fetch_db_filepath())
    conn = sql.create_connection(sql.fetch_db_filepath())
    cur = conn.cursor()
    sql.create_specific_table(conn)

    # upload data ...
    JD = [tracklet_dic['JD'] for tracklet_dic in tracklet_dict_list]
    HP = [tracklet_dic['HP'] for tracklet_dic in tracklet_dict_list]
    tracklet_name = [tracklet_dic['tracklet_name']
                     for tracklet_dic in tracklet_dict_list]
    sql.upsert_tracklets(conn, JD, HP, tracklet_name, tracklet_dict_list)

    # query the data & check that required # of dictionaries are returned
    list_of_tuples = sql.query_tracklets_jdhp(conn, JD[0], HP[0])
    assert isinstance(list_of_tuples, list) and len(list_of_tuples) == 1

    # now delete a tracklet & check that one less dictionary is subsequently returned
    sql.delete_tracklet(conn, tracklet_dict_list[0]['tracklet_name'])
    list_of_tuples = sql.query_tracklets_jdhp(conn, JD[0], HP[0])
    assert isinstance(list_of_tuples, list) and len(list_of_tuples) == 0

    # Delete the db to facilitate future testing
    os.remove(sql.fetch_db_filepath())


@pytest.mark.parametrize(('tracklet_dict_list'), [test_tracklet_dict_list])
def test_delete_tracklets(tracklet_dict_list):
    '''Test deleting multiple tracklets.'''
    # set up db & table
    if os.path.isfile(sql.fetch_db_filepath()):
        os.remove(sql.fetch_db_filepath())
    conn = sql.create_connection(sql.fetch_db_filepath())
    cur = conn.cursor()
    sql.create_specific_table(conn)

    # upload data ...
    JD = [tracklet_dic['JD'] for tracklet_dic in tracklet_dict_list]
    HP = [tracklet_dic['HP'] for tracklet_dic in tracklet_dict_list]
    tracklet_name = [tracklet_dic['tracklet_name'] for tracklet_dic in tracklet_dict_list]
    sql.upsert_tracklets(conn, JD, HP, tracklet_name, tracklet_dict_list)

    # query the data & check that required # of dictionaries are returned
    list_of_tuples = sql.query_tracklets_jdhp(conn, JD[0], HP[0])
    assert isinstance(list_of_tuples, list) and len(list_of_tuples) == 1
    list_of_tuples = sql.query_tracklets_jdhp(conn, JD[1], HP[1])
    assert isinstance(list_of_tuples, list) and len(list_of_tuples) == 1

    # now delete two tracklets & check that two fewer dictionaries remain
    sql.delete_tracklets(conn, [tracklet_dict_list[0]['tracklet_name'],
                                tracklet_dict_list[1]['tracklet_name']])

    list_of_tuples = sql.query_tracklets_jdhp(conn, JD[0], HP[0])
    assert isinstance(list_of_tuples, list) and len(list_of_tuples) == 0
    list_of_tuples = sql.query_tracklets_jdhp(conn, JD[1], HP[1])
    assert isinstance(list_of_tuples, list) and len(list_of_tuples) == 0

    # Delete the db to facilitate future testing
    os.remove(sql.fetch_db_filepath())


# End of file.
