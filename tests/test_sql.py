"""
	Tests of the sql module
    
	Currently (20200805) incomplete:
    GOOD POINTS:
        MJP transfered most of the sql demo/tests from
        Demonstrate_SQLandPreCalc.ipynb
    BAD POINTS
        (i) Some later tests are dependent on NBODY Integrations
        (ii) Sifter tests not in here
"""

# Import third-party packages
# --------------------------------------------------------------
import numpy as np
import sys, os
from astropy_healpix import healpy
import pytest
import sqlite3

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

# Using decorators to do the db cleaning ...
# - This deletes any extant db before the test
# - And deletes any extand db after the test
def db_handler(func):
    def inner_function(*args, **kwargs):
    
        # Delete db if any exists for any reason
        filepath = sql.DB.fetch_db_filepath()
        if os.path.isfile(filepath):
            os.remove(filepath)
            
        # Run test function
        result = func(*args, **kwargs)

        # Delete the db
        if os.path.isfile(filepath):
            os.remove(filepath)
        
        # Return any result
        return result

    return inner_function

@db_handler
def test_DB():
    
    # Instantiate DB object & check it has the expected attributes
    db = sql.DB()
    assert hasattr(db,'create_table')
    assert hasattr(db,'fetch_db_filepath')
    assert hasattr(db,'create_connection')
    assert hasattr(db,'db_file')
    assert hasattr(db,'conn')
    cur  = db.conn.cursor()
    assert isinstance( cur, sqlite3.Cursor )
        
    # Check db exists
    # ( if no db exists, the above instantiation
    #   should have created one)
    assert os.path.isfile( db.db_file )
    


@db_handler
def test_SQLChecker_TableCreation():

    # Instantiate
    C = sql.SQLChecker()
    assert hasattr(C,'create_table')
    assert hasattr(C,'fetch_db_filepath')
    assert hasattr(C,'create_connection')
    assert hasattr(C,'db_file')
    assert hasattr(C,'conn')
    cur  = C.conn.cursor()
    assert isinstance( cur, sqlite3.Cursor )
    
    # Check db exists
    assert os.path.isfile( C.db_file )
    
    # Create all checker tables
    C.create_all_checker_tables()

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
    expected_names = ['coeff_id', 'object_id'] + C.generate_sector_field_names()
    assert np.all( [n in expected_names for n in names ])
  
  
"""

@db_handler
def test_inserts():
    
    # (0) Create demo MSC(s) & create empty db using functions from test_orbit_cheby
    MSCs = test_orbit_cheby.test_loader_from_nbodysim(test_orbit_cheby.orbfit_filenames[0])
    
    # (1) Instantiate
    C = sql.SQLChecker()

    # (2) Test Low level direct designation insert: result = object_id = 1
    result = C.insert_desig(conn ,MSCs[0].primary_unpacked_provisional_designation )
    assert result == 1

    # (3a) Test that there is something in the db by executing a query
    cur = conn.cursor()
    cur.execute("SELECT * FROM object_coefficients WHERE object_id=?", ( 1, ))
    result = cur.fetchall()
    assert len(result) ==1
    # (3b) There is also an appropriate query functionality: result = object_id = 1
    result = C.query_number_by_desig( MSCs[0].primary_unpacked_provisional_designation )
    assert result == 1
    
    # (4) Test MSC upsert
    object_id = result
    C.upsert_MSC( MSCs[0] , object_id )

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

"""
