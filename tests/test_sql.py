"""
	Tests of the sql module
	Currently (20200805) incomplete: MJP transfered most of the sql demo/tests from Demonstrate_SQLandPreCalc.ipynb
"""

# Import third-party packages
# --------------------------------------------------------------
import numpy as np
import sys, os
from astropy_healpix import healpy
import pytest

# Import neighboring packages
# --------------------------------------------------------------
sys.path.append(os.path.dirname(os.path.dirname(
    os.path.realpath(__file__))))
from cheby_checker import sql
from tests import test_orbit_cheby
from cheby_checker import orbit_cheby
from cheby_checker import obs_pos

# Files / Directories
# --------------------------------------------------------------
HEAD_DIR = os.path.dirname(os.path.realpath(os.getcwd()))
DATA_DIR = os.path.join(HEAD_DIR, 'dev_data')
filenames = [os.path.join(DATA_DIR, file)
              for file in ['30101.eq0_horizons', '30102.eq0_horizons']]


# Tests ...
# --------------------------------------------------------------
def test_func_create_db_and_tables():
    
    # In order to save data, we require sql-db to exist, so let's set that up...
    # Force deletion then creation of db...
    if os.path.isfile( sql.fetch_db_filepath() ):
        os.remove(sql.fetch_db_filepath())
    conn = sql.create_connection( sql.fetch_db_filepath() )
    cur  = conn.cursor()
    
    # Test creation of db
    assert os.path.isfile( os.path.join( sql.fetch_db_filepath() ) ), 'no db'

    # Create required table(s) for cheby-coeff storage
    sql.create_object_coefficients_table(conn)
    sql.create_objects_by_jdhp_table(conn)
    sql.create_object_desig_table(conn)

    # Double-check that this worked by getting the count of tables with the name
    # - if the count is 1, then table exists
    cur.execute('SELECT name from sqlite_master WHERE type = "table" AND name = "objects_by_jdhp"')
    assert len(cur.fetchone()) == 1 , 'jdhp table does not exist'
    cur.execute('SELECT name from sqlite_master WHERE type = "table" AND name = "object_coefficients"')
    assert len(cur.fetchone()) == 1 , 'coeff table does not exist'
    cur.execute('SELECT name from sqlite_master WHERE type = "table" AND name = "object_desig"')
    assert len(cur.fetchone()) == 1 , 'coeff table does not exist'

    # Test that the expected column names are in the *object_coefficients* table
    cur.execute("SELECT * FROM object_coefficients ")
    names = [description[0] for description in cur.description]
    expected_names = ['coeff_id', 'object_id'] + sql.generate_sector_field_names()
    assert np.all( [n in expected_names for n in names ])
    
    
    # Return connection
    return conn

def test_inserts():
    
    # (0) Create demo MSC(s) & create empty db using functions from test_orbit_cheby
    MSCs = test_orbit_cheby.test_loader_from_nbodysim(test_orbit_cheby.orbfit_filenames[0])
    
    # (1) Ensure db exists
    conn = test_func_create_db_and_tables()

    # (2) Test Low level direct designation insert: result = object_id = 1
    result = sql.insert_desig(conn ,MSCs[0].primary_unpacked_provisional_designation )
    assert result == 1

    # (3) Query functionality: result = object_id = 1
    result = sql.query_number_by_desig(conn, MSCs[0].primary_unpacked_provisional_designation)
    assert result == 1
    
    # (4) Test MSC upsert
    object_id = result
    sql.upsert_MSC(conn ,MSCs[0] , object_id )

    # (5) Test that there is something in the database by querying it ...
    cur = conn.cursor()
    cur.execute("SELECT * FROM object_coefficients WHERE object_id=?", ( object_id, ))
    result = cur.fetchall()
    assert len(result) ==1

    # (6) Test that the expected sectors are populated
    M = MSCs[0]
    # names of the sectors we expect to have (binary blob) numerical data in them ...
    expected_sector_field_names = sql.generate_sector_field_names( sector_dict = \
        { sector_num: orbit_cheby.Base().map_sector_number_to_sector_start_JD(sector_num , orbit_cheby.Base().standard_MJDmin) \
            for sector_num in M.sector_coeffs.keys()})
    # names of the fields in the object_coefficients table (getting from query)
    names = [description[0] for description in cur.description]
    for n,name in enumerate(names):
        # We expect these to be populated
        if name in expected_sector_field_names or name in ['coeff_id', 'object_id']:
            assert result[0][n] is not None
        # We expect these fields to be empty because the input MSC has data for only a subset of sectors
        else:
            assert result[0][n] is None



def test_queries():

    # (0) Create demo MSC(s) & create empty db using functions from test_orbit_cheby
    #  NB: test_loader_from_nbodysim() calls test_convenience_call_to_mpc_nbody_nbodysim() 
    MSCs = test_orbit_cheby.test_loader_from_nbodysim(test_orbit_cheby.orbfit_filenames[0])
    M = MSCs[0]
    
    # (1) Ensure db exists
    conn = test_func_create_db_and_tables()

    # (2) Low level direct inserts to put MSC into db so that we have something to query against!
    object_id = sql.insert_desig(conn ,MSCs[0].primary_unpacked_provisional_designation )
    sql.upsert_MSC(conn ,M , object_id )

    # (3) Execute the query
    result = sql.query_object_coefficients(conn , object_id)

    # (4) Examine the characteristics of the output
    assert isinstance(result, dict)
    
    # (5) Check the sectors we expect to have (binary blob) numerical data in them ...
    expected_sector_field_names = sql.generate_sector_field_names( sector_dict = \
        { sector_num: orbit_cheby.Base().map_sector_number_to_sector_start_JD(sector_num , orbit_cheby.Base().standard_MJDmin) \
            for sector_num in M.sector_coeffs.keys()})
    assert len(result) == len(expected_sector_field_names)
    for k,v in result.items():
        assert k in expected_sector_field_names
        assert v.ndim == 2
        assert np.all(M.sector_coeffs[ int(k.split("_")[1]) ] == v )


    # (6) Query checker db for object-coefficients over a limited subset of sectors
    #- This would likely *NOT* be performed by the user directly (but would instead be called by some other function)
    n=4
    first_few_populated_sector_numbers = list(M.sector_coeffs.keys())[:n]
    result = sql.query_object_coefficients(conn ,object_id , sector_numbers = first_few_populated_sector_numbers )
    
    # (7) Examine the characteristics of the output
    assert isinstance(result, dict)
    assert len(result) == n
    for k,v in result.items():
        assert k in expected_sector_field_names
        assert v.ndim == 2
        assert np.all(M.sector_coeffs[ int(k.split("_")[1]) ] == v )


def test_healpix():

    # (0) Create demo MSC(s) & create empty db using functions from test_orbit_cheby
    MSCs = test_orbit_cheby.test_loader_from_nbodysim(test_orbit_cheby.orbfit_filenames[0])
    M = MSCs[0]
    
    # (1) Ensure db exists
    conn = test_func_create_db_and_tables()

    # (2) Low level direct inserts to put MSC into db so that we have something to query against!
    object_id = sql.insert_desig(conn ,MSCs[0].primary_unpacked_provisional_designation )
    sql.upsert_MSC(conn ,M , object_id )

    # (3) Explicitly evaluate the HP functionality to get the appropriate HP for each date
    # - See "Demonstrate_Orbital_Chebyshev_Functionality" demo for more details
    # - NB We need to provide dates that are supported by the MSC
    n=10
    JDlist = np.linspace(M.TDB_init, 0.5*(M.TDB_final+M.TDB_init), num=n )
    observatoryXYZ = np.array( [obs_pos.ObsPos().get_heliocentric_equatorial_xyz(jd , obsCode='500') for jd in JDlist] ).T
    HPlist= M.generate_HP( JDlist ,
                    observatoryXYZ,
                    APPROX = True )
                    
    assert len(HPlist) == len(JDlist) == n
    
    # (4) Now the code to do the insert ...
    sql.insert_HP(conn, JDlist, HPlist, object_id)

    # (5) sql query to get the *object_ids* that are in the listed HP on the given day
    #May be of little use in practice, but helpful for development
    result = sql.query_jd_hp(conn, JDlist[0], HPlist[:5])
    
    # (6) Examine the characteristics of the output
    print('result=',result)
    assert isinstance(result, list)
    assert result == [object_id]
