# -*- coding: utf-8 -*-
# sifter/tests/test_tracklet

"""
--------------------------------------------------------------
tests of sifter's Tracklet class

Jan 2020
Matt Payne & Mike Alexandersen

--------------------------------------------------------------
"""


# Import third-party packages
# --------------------------------------------------------------
import os
import numpy as np
import pytest 

# Import neighboring packages
# --------------------------------------------------------------

#@pytest.mark.skip(reason="Sifter Code/Tests have NOT been reviewed in 2022 by MJP")
#from cheby_checker import sifter_precalc as precalc
from cheby_checker import sql
DB = sql.SQLSifter()
DATA_DIR = os.path.join(os.path.dirname(__file__), 'dev_data')


# Convenience data / functions to aid testing
# --------------------------------------------------------------
test_tracklet = ['     K11Q99F*~C2011 08 29.52378 01 57 34.729+14 35 44.64         22.8 rc~0qBd568',
                 '     K11Q99F ~C2011 08 29.61470 01 57 34.343+14 35 42.59         22.9 rc~0qBd568']


def convenience_func_create_db_and_tables():
    """
    In order to save data, we require sql-db to exist,
    so let's set that up...
    Force deletion then creation of db...
    """
    if os.path.isfile(DB.fetch_db_filepath()):
        os.remove(DB.fetch_db_filepath())
    conn = DB.create_connection()
    cur = conn.cursor()

    # Create required table(s)
    DB.create_sifter_tables(conn)

    # Double-check that this worked by getting the count of tables with the name
    # - if the count is 1, then table exists
    cur.execute('SELECT name from sqlite_master WHERE '
                'type = "table" AND name = "tracklets"')
    res = cur.fetchone()
    assert len(res) == 1, 'table does not exist'
    conn.close()


# Actual tests ...
# --------------------------------------------------------------
@pytest.mark.skip(reason="Sifter Code/Tests have NOT been reviewed in 2022 by MJP")
def test_instantiation():
    """Test instantiation of the Tracklets class with no observations."""
    assert isinstance(precalc.Tracklets(), precalc.Tracklets)


@pytest.mark.skip(reason="Sifter Code/Tests have NOT been reviewed in 2022 by MJP")
#@pytest.mark.parametrize(('tracklet_obs'), [test_tracklet])
def test_parse_tracklet_observations(tracklet_obs):
    """Test that observations get parsed correctly."""
    T = precalc.Tracklets()

    # call parse_tracklet_observations
    tracklet_dictionary = T.parse_tracklet_observations(tracklet_obs)

    # check that the returned results are as expected
    assert isinstance(tracklet_dictionary, dict)
    assert tracklet_obs == tracklet_dictionary['observations']
    assert 'JD' in tracklet_dictionary
    assert 'HP' in tracklet_dictionary
    assert 'tracklet_name' in tracklet_dictionary


@pytest.mark.skip(reason="Sifter Code/Tests have NOT been reviewed in 2022 by MJP")
#@pytest.mark.parametrize(('observation_pair_list'), [[test_tracklet, test_tracklet]])
def test_save_tracklets(observation_pair_list):
    """Test that creating a db works and saving stuff to it works."""
    # Create db from scratch
    convenience_func_create_db_and_tables()

    # Set up a Tracklet & use parse_tracklet_observations to get JD, HP, ...
    T = precalc.Tracklets()
    tracklet_dictionary_list = [T.parse_tracklet_observations(obs_pair)
                                for obs_pair in observation_pair_list]

    # Now save the data in the db
    T.save_tracklets(tracklet_dictionary_list)

    # Test the data was uploaded and can be downloaded
    cur = T.conn.cursor()
    cur.execute('SELECT * from tracklets')
    f = cur.fetchall()
    assert(len(f) == 2 and np.all([len(_) > 3 for _ in f])), 'data not uploaded'

    # Completely delete db to facilitate future testing
    os.remove(DB.fetch_db_filepath())


@pytest.mark.skip(reason="Sifter Code/Tests have NOT been reviewed in 2022 by MJP")
#@pytest.mark.parametrize(('observation_pairs'), [test_tracklet])
def test_instantiation_with_observations(observation_pairs):
    """Test instantiation of the Tracklets class with some observations."""
    # Create db from scratch
    convenience_func_create_db_and_tables()

    # instantiate with observation_pair
    T = precalc.Tracklets(observation_pairs)

    # test that the above caused the tracklet to be uploaded to db
    cur = T.conn.cursor()
    cur.execute('SELECT * from tracklets')
    f = cur.fetchone()
    assert(len(f) > 3), 'data not uploaded'

    # Completely delete db to facilitate future testing
    os.remove(DB.fetch_db_filepath())


# End of file.
