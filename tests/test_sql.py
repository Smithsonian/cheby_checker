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
import pickle

# Import neighboring packages
from cheby_checker import sql
from cheby_checker.cheby_checker import Base


# Files / Directories
# --------------------------------------------------------------
HEAD_DIR = os.path.dirname(os.path.realpath(os.getcwd()))
DATA_DIR = os.path.join(HEAD_DIR, 'dev_data')
filenames = [os.path.join(DATA_DIR, file)
              for file in ['30101.eq0_horizons', '30102.eq0_horizons']]


# Tests ...
# --------------------------------------------------------------


@pytest.fixture()
def db():
    db = sql.DB()
    db.clear_database()
    return db


# --------------------------------------------------------------
# Tests of class instantiation ...
# --------------------------------------------------------------
def test_DB(db):
    '''
    Test the basic operation of the "DB" class object.
    
    The sql.DB class handles basic database connections & locations
    
    N.B. (1)
    Instantiation invoves the following two functions,
    so we test that a database-filepath is created and
    that a connection is available
        .fetch_db_filepath()
        .create_connection()
        
    N.B. (2)
    The only other functionality available in the "DB"
    class object is the convenience function,
    *create_table*, so we test that here too.

    '''
    
    # Check has the expected function-attributes
    assert hasattr(db,'create_connection')
    assert hasattr(db,'create_table')
    
    # Check has the expected variable-attributes
    assert hasattr(db,'conn')
    assert hasattr(db, 'cur')
        
    # Check db is connected.
    assert not db.conn.closed
    
    # Attempt to create a table in the database using the *create_table* function
    sql_statement = """
        CREATE TABLE IF NOT EXISTS test_table_name (
        object_id integer PRIMARY KEY,
        unpacked_primary_provisional_designation TEXT UNIQUE);
    """
    db.create_table( sql_statement)
    
    # Check that table creation worked by getting the count of tables with the name
    # - if the count is 1, then table exists
    db.cur.execute("SELECT table_name FROM information_schema.tables WHERE table_schema='public' AND table_type='BASE TABLE' AND table_name='test_table_name';")
    assert len(db.cur.fetchone()) == 1 , 'test_table_name table does not exist'
    

def test_SQLChecker(db):
    """
    Test the basic operation of the "SQLChecker" class object:
     - SQLChecker handles all database interactions required by cheby_checker
    
    In this first test of SQLChecker, we just test instantiation.
    
    N.B. (1)
    Instantiating "SQLChecker" is the same as instantiating "DB"
    So this test is the same as above
    
    N.B. (2)
    More detailed tests of the additional functionality in
    "SQLChecker" are performed in other test-fuunctions below

    """
    
    # Instantiate DB object
    db = sql.SQLChecker()
    
    # Check has the expected function-attributes
    assert hasattr(db,'create_connection')
    assert hasattr(db,'create_table')
    
    # Check has the expected variable-attributes
    assert hasattr(db,'conn')
    
    # Attempt to create a table in the database using the ** function
    sql_statement = """
            CREATE TABLE IF NOT EXISTS test_table_name (
            object_id integer PRIMARY KEY,
            unpacked_primary_provisional_designation TEXT UNIQUE);
    """
    db.create_table(sql_statement)
    
    # Check that table creation worked by getting the count of tables with the name
    # - if the count is 1, then table exists
    db.cur.execute("SELECT table_name from information_schema.tables WHERE table_name = 'test_table_name'")
    assert len(db.cur.fetchone()) == 1, 'test_table_name table does not exist'
    


    
# --------------------------------------------------------------
# Tests of sector field names ...
# --------------------------------------------------------------
def test_SQLChecker_generate_sector_field_names(db):
    """
    Test the *generate_sector_field_names* function(s) in "SQLChecker"
    
    The sectors are given "names" / "labels" / "field-headings"
     - By default the names are derived from the keys in Base().get_required_sector_dict()
     - Base().get_required_sector_dict() has been tested in "test_cheby_checker"
     
    Here I explicitly test / demonstrate that I know what the sector-names look like ...
    
    """
    # Instantiate
    C = sql.SQLChecker()
    cur = C.conn.cursor()
    
    # Create all checker tables (creates two tables)
    C.create_all_checker_tables()
    
    # Call the function
    sector_name_list = C.generate_sector_field_names( sector_dict = Base().get_required_sector_dict() )
    
    # Check results
    required_sector_dict = Base().get_required_sector_dict()
    assert sector_name_list[0] == f'sector_{0}_{required_sector_dict[0]}'
    assert sector_name_list[1] == f'sector_{1}_{required_sector_dict[1]}'


# --------------------------------------------------------------
# Tests of table creation ...
# --------------------------------------------------------------

def test_SQLChecker_TableCreation(db):
    """
    Test the table creation function(s) in "SQLChecker"
     - These are convenience functions that create the three tables I/we expect to be required to
       support all "cheby checker" operations
    """

    # Instantiate
    C = sql.SQLChecker()
    cur = C.conn.cursor()
    
    # Create all checker tables (creates *TWO* tables)
    C.create_all_checker_tables()

    # Double-check that this worked by getting the count of tables with the name
    # - if the count is 1, then table exists

    cur.execute("SELECT table_name from information_schema.tables WHERE table_name = 'objects_by_jdhp'")
    assert len(cur.fetchone()) == 1, 'jdhp table does not exist'
    cur.execute("SELECT table_name from information_schema.tables WHERE table_name = 'object_coefficients'")
    assert len(cur.fetchone()) == 1, 'coeff table does not exist'

    # Test that the expected column names are in the *object_coefficients* table
    cur.execute("SELECT * FROM object_coefficients ")
    names = [description[0] for description in cur.description]
    expected_names = ['object_coeff_id', 'unpacked_primary_provisional_designation'] + C.generate_sector_field_names()
    assert np.all( [n in expected_names for n in names ])
  
    # Test that the expected column names are in the *objects_by_jdhp* table
    cur.execute("SELECT * FROM objects_by_jdhp ")
    names = [description[0] for description in cur.description]
    expected_names = ['jdhp_id', 'jd', 'hp', 'object_coeff_id']
    assert np.all( [n in expected_names for n in names ])


# --------------------------------------------------------------
# Tests of basic data-insert routines ...
# --------------------------------------------------------------

def test_upsert_coefficients(db):
    """
    Test the insertion of coefficients into the coeff table
    NB:
     - We want this to be an upsert/replace with complete removal
    """
    
    # Instantiate & create all checker tables (creates two tables)
    C = sql.SQLChecker()
    C.create_all_checker_tables()
    
    # ------------- 1 ----------------
    # Simple test of insert
    # --------------------------------
    # Define some data to be inserted
    # NB: These will be used to store pickled dictionaries, so let's test that ...
    unpacked_primary_provisional_designation = '2020 AB'
    sector_field_names = C.generate_sector_field_names()[:2]
    sector_values_raw  = [ np.array([[n,n],[n,n]]) for n, _ in enumerate(sector_field_names) ]
    sector_values      = [ pickle.dumps( _, pickle.HIGHEST_PROTOCOL) for _ in sector_values_raw]
    
    # Call upsert function
    object_coeff_id = C.upsert_coefficients(unpacked_primary_provisional_designation , sector_field_names, sector_values)

    # Check that the returned id is an integer:
    # - Because we have a new db & table for this test, the indexing should start from one (apparently), ...
    assert object_coeff_id == 1
    
    # Query the database and see what the returned data looks like
    sqlstr = "SELECT * FROM object_coefficients WHERE unpacked_primary_provisional_designation=?"
    cur = C.conn.cursor()
    cur.execute(sqlstr , ( unpacked_primary_provisional_designation, ))
    result_raw     = cur.fetchall()[0]
    
    # Check that the returned results are the same as the inserted values
    result_object_coeff_id  = result_raw[0]
    result_desig            = result_raw[1]
    result_sectors          = [ pickle.loads( _ ) for _ in result_raw if _ != None and not isinstance(_,(str,int))]
    assert result_object_coeff_id == object_coeff_id
    assert result_desig == unpacked_primary_provisional_designation
    assert np.all( [ a==b for a,b in zip(sector_values_raw , result_sectors) ] )




    # ------------- 2 ----------------
    # Simple test of additional data insert
    # --------------------------------
    unpacked_primary_provisional_designation = '2021 XY'
    sector_field_names = C.generate_sector_field_names()[:3]
    sector_values_raw  = [ np.array([[n,n],[n,n]]) for n, _ in enumerate(sector_field_names) ]
    sector_values      = [ pickle.dumps( _, pickle.HIGHEST_PROTOCOL) for _ in sector_values_raw]
    
    # Call upsert function
    object_coeff_id = C.upsert_coefficients(unpacked_primary_provisional_designation , sector_field_names, sector_values)

    # Check that the returned id is an integer:
    # - This should be the second line ...
    assert object_coeff_id == 2
    
    # Query the database and see what the returned data looks like
    sqlstr = "SELECT * FROM object_coefficients WHERE unpacked_primary_provisional_designation=?"
    cur = C.conn.cursor()
    cur.execute(sqlstr , ( unpacked_primary_provisional_designation, ))
    result_raw     = cur.fetchall()[0]
    
    # Check that the returned results are the same as the inserted values
    result_object_coeff_id  = result_raw[0]
    result_desig            = result_raw[1]
    result_sectors          = [ pickle.loads( _ ) for _ in result_raw if _ != None and not isinstance(_,(str,int))]
    assert result_object_coeff_id == object_coeff_id
    assert result_desig == unpacked_primary_provisional_designation
    assert np.all( [ a==b for a,b in zip(sector_values_raw , result_sectors) ] )


    

    # ------------- 3 ----------------
    # Now we try to re-insert / update the same '2021 XY' object as in ---2--- above
    # But note that there are a different (smaller) number of elements being inserted now
    # So there should still only be 2 rows in the table (but the index will have updated because that's how it works ...)
    # --------------------------------

    unpacked_primary_provisional_designation = '2021 XY'
    sector_field_names = C.generate_sector_field_names()[:2]
    sector_values_raw  = [ np.array([[17*n,17*n],[17*n,17*n]]) for n, _ in enumerate(sector_field_names) ]
    sector_values      = [ pickle.dumps( _, pickle.HIGHEST_PROTOCOL) for _ in sector_values_raw]
    
    # Call upsert function
    object_coeff_id = C.upsert_coefficients(unpacked_primary_provisional_designation , sector_field_names, sector_values)
    
    # Check that the returned id is as expected (the index updates)
    assert object_coeff_id == 3
    
    # Check that there's 2 (not 3) rows in the database
    sqlstr = "SELECT object_coeff_id,unpacked_primary_provisional_designation FROM object_coefficients"
    cur    = C.conn.cursor()
    cur.execute(sqlstr)
    assert len(cur.fetchall()) == 2
    
    # Query the database and see what the returned data looks like
    sqlstr = "SELECT * FROM object_coefficients WHERE unpacked_primary_provisional_designation=?"
    cur = C.conn.cursor()
    cur.execute(sqlstr , ( unpacked_primary_provisional_designation, ))
    result_raw     = cur.fetchall()[0]
    
    # Check that the retuurned coefficients are the same as this update and NOT like those inserted in ---2---
    result_object_coeff_id  = result_raw[0]
    result_desig            = result_raw[1]
    result_sectors          = [ pickle.loads( _ ) for _ in result_raw if _ != None and not isinstance(_,(str,int))]
    assert len(result_sectors) == 2
    assert np.all( [ a==b for a,b in zip(sector_values_raw , result_sectors) ] )




def test_insert_HP(db):
    """
    TEST THE FUNCTION THAT INSERTS LISTS OF JD & HP
    
    Note that this implicitly tests *delete_JDHP_by_object_coeff_id*, as a delete
    step is used to ensure old data for objects with the same object_coeff_id are removed
    
    """

    # Instantiate & create all checker tables (creates three tables)
    C = sql.SQLChecker()
    C.create_all_checker_tables()

    # -------------------1---------------------------------
    # Simple test of  data insert
    # -----------------------------------------------------
    # Create sample inputs
    JDlist = [2440000, 2440001, 2440002]
    HPlist = [17,      17,      18]
    object_coeff_id = 123
    
    # Call the insert function
    C.insert_HP( JDlist, HPlist, object_coeff_id )

    # Query the db and examine the results
    sqlstr = "SELECT * FROM objects_by_jdhp WHERE object_coeff_id=?"
    cur = C.conn.cursor()
    cur.execute(sqlstr , ( object_coeff_id, ))
    result_raw     = cur.fetchall()

    # Check the inputs are the same as the outputs
    for n,r in enumerate(result_raw):
        id,jd,hp,obj_co_id = r
        assert jd==JDlist[n]
        assert hp==HPlist[n]
        assert obj_co_id==object_coeff_id

    # -------------------2---------------------------------
    # Simple test of additional data insert
    # -----------------------------------------------------
    JDlist = [2450000, 2450001, 2450002]
    HPlist = [57,      57,      58]
    object_coeff_id = 456
    
    # Call the insert function
    C.insert_HP( JDlist, HPlist, object_coeff_id )

    # Query the db and examine the results
    sqlstr = "SELECT * FROM objects_by_jdhp WHERE object_coeff_id=?"
    cur = C.conn.cursor()
    cur.execute(sqlstr , ( object_coeff_id, ))
    result_raw     = cur.fetchall()

    # Check the inputs are the same as the outputs
    for n,r in enumerate(result_raw):
        id,jd,hp,obj_co_id = r
        assert jd==JDlist[n]
        assert hp==HPlist[n]
        assert obj_co_id==object_coeff_id


    # -------------------3---------------------------------
    # *** THESE ARE FOR THE SAME object_coeff_id AS IN ---1--- ABOVE
    # *** BUT WE HAVE FEWER ENTRIES, AND THEY HAVE DIFFERENT VALUES
    # *** WE WANT TO TEST THAT ALL THE OLD STUFF IS DELETED
    # -----------------------------------------------------
    JDlist = [2440000, 2440001]
    HPlist = [19,      19]
    object_coeff_id = 123

    # Call the insert function
    C.insert_HP( JDlist, HPlist, object_coeff_id )

    # Query the db and examine the results
    sqlstr = "SELECT * FROM objects_by_jdhp WHERE object_coeff_id=?"
    cur = C.conn.cursor()
    cur.execute(sqlstr , ( object_coeff_id, ))
    result_raw     = cur.fetchall()

    # Check the inputs are the same as the outputs
    for n,r in enumerate(result_raw):
        id,jd,hp,obj_co_id = r
        assert jd==JDlist[n]
        assert hp==HPlist[n]
        assert obj_co_id==object_coeff_id





"""

@db_handler
def test_inserts():
    
    # (0) Create demo MSC(s) & create empty db using functions from test_orbit_cheby
    MSCs = test_orbit_cheby.test_loader_from_nbodysim(test_orbit_cheby.orbfit_filenames[0])
    
    # (1) Instantiate
    C = sql.SQLChecker()

    # (2) Test Low level direct designation insert: result = object_id = 1
    result = C.insert_desig(conn ,MSCs[0].unpacked_primary_provisional_designation )
    assert result == 1

    # (3a) Test that there is something in the db by executing a query
    cur = conn.cursor()
    cur.execute("SELECT * FROM object_coefficients WHERE object_id=?", ( 1, ))
    result = cur.fetchall()
    assert len(result) ==1
    # (3b) There is also an appropriate query functionality: result = object_id = 1
    result = C.query_number_by_desig( MSCs[0].unpacked_primary_provisional_designation )
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

"""

# --------------------------------------------------------------
# Tests of basic query routines ...
# --------------------------------------------------------------

def test_query_object_coefficients(db):
    """
    Test the convenience query that searches the object_coefficients and returns the (unpickled) data
    """

    # Instantiate & create all checker tables (creates two tables)
    C = sql.SQLChecker()
    C.create_all_checker_tables()
    
    # ------------- 1 ----------------
    # Simple test of query
    # --------------------------------
    # Define some data to be inserted
    # NB: The *query_object_coefficient* routine is designed to un-pickle the stored data,
    #     so we need to supply pickled data ...
    unpacked_primary_provisional_designation = '2020 AB'
    sector_field_names = C.generate_sector_field_names()[:2]
    sector_values_raw  = [ np.array([[n,n],[n,n]]) for n, _ in enumerate(sector_field_names) ]
    sector_values      = [ pickle.dumps( _, pickle.HIGHEST_PROTOCOL) for _ in sector_values_raw]
    
    # Call the convenient upsert function (tested above)
    object_coeff_id = C.upsert_coefficients(unpacked_primary_provisional_designation , sector_field_names, sector_values)

    # Now call the query function that we want to test
    result_dict = C.query_object_coefficients(unpacked_primary_provisional_designation)
    
    # Check the results are as expected
    assert isinstance(result_dict, dict)
    for n,v in zip(sector_field_names , sector_values_raw):
        assert n in result_dict
        assert np.all( v == result_dict[n] )
    for k,v in result_dict.items():
        if k not in sector_field_names:
            assert v is None

    # ------------- 2 ----------------
    # Simple test of query
    # NB: THIS IS A DIFFERENT DESIGNATION
    # --------------------------------
    # Define some data to be inserted
    unpacked_primary_provisional_designation = '2021 XY'
    sector_field_names = C.generate_sector_field_names()[:2]
    sector_values_raw  = [ np.array([[2*n,2*n],[2*n,2*n]]) for n, _ in enumerate(sector_field_names) ]
    sector_values      = [ pickle.dumps( _, pickle.HIGHEST_PROTOCOL) for _ in sector_values_raw]
    
    # Call the convenient upsert function (tested above)
    object_coeff_id = C.upsert_coefficients(unpacked_primary_provisional_designation , sector_field_names, sector_values)

    # Now call the query function that we want to test
    result_dict = C.query_object_coefficients(unpacked_primary_provisional_designation)
    
    # Check the results are as expected
    assert isinstance(result_dict, dict)
    for n,v in zip(sector_field_names , sector_values_raw):
        assert n in result_dict
        assert np.all( v == result_dict[n] )
    for k,v in result_dict.items():
        if k not in sector_field_names:
            assert v is None


    # ------------- 3 ----------------
    # NB: THIS IS THE SAME DESIGNATIONS AS IN ----1----
    # So we expect the data to have been updated/replaced
    # --------------------------------
    # Define some data to be inserted
    unpacked_primary_provisional_designation = '2020 AB'
    sector_field_names = C.generate_sector_field_names()[:2]
    sector_values_raw  = [ np.array([[3*n,3*n],[3*n,3*n]]) for n, _ in enumerate(sector_field_names) ]
    sector_values      = [ pickle.dumps( _, pickle.HIGHEST_PROTOCOL) for _ in sector_values_raw]
    
    # Call the convenient upsert function (tested above)
    object_coeff_id = C.upsert_coefficients(unpacked_primary_provisional_designation , sector_field_names, sector_values)

    # Now call the query function that we want to test
    result_dict = C.query_object_coefficients(unpacked_primary_provisional_designation)
    
    # Check the results are as expected
    assert isinstance(result_dict, dict)
    for n,v in zip(sector_field_names , sector_values_raw):
        assert n in result_dict
        assert np.all( v == result_dict[n] )
    for k,v in result_dict.items():
        if k not in sector_field_names:
            assert v is None



def test_query_desig_by_object_coeff_id(db):
    """
    Test the query_desig_by_object_coeff_id
    NB(1): When *query_desig_by_object_coeff_id* is called with a list of object_ids,
           it returns the associated unpacked_primary_provisional_designations
        
    NB(2): Not sure *query_desig_by_object_coeff_id* is ever used any more.
           May be pointtless to keep it.
    
    """

    # Instantiate & create all checker tables (creates two tables)
    C = sql.SQLChecker()
    C.create_all_checker_tables()
    
    # Define some data to be inserted
    expected_desig_by_id = {}
    desigs = ['2020 AB', '2020 XY', '2021 PQ' ]
    for i, desig in enumerate( desigs ):
    
        sector_field_names = C.generate_sector_field_names()[:2]
        sector_values_raw  = [ np.array([[i + 2*n,i + 2*n],[i + 2*n,i + 2*n]]) for n, _ in enumerate(sector_field_names) ]
        sector_values      = [ pickle.dumps( _, pickle.HIGHEST_PROTOCOL) for _ in sector_values_raw]
    
        # Call convenience upsert function to get data into the system
        # NB: *upsert_coefficients* returns the object_coeff_id of the inserted coefficients
        object_coeff_id = C.upsert_coefficients(desig , sector_field_names, sector_values)
        
        # Save object_coeff_id into dict, so that we use the lengthening dict as the basis for a query below
        expected_desig_by_id[object_coeff_id]=desig
        
        # Call the *query_desig_by_object_coeff_id* function that we are trying to test
        # NB Returned dict has : key is object_coeff_id, value is desig
        result_dict = C.query_desig_by_object_coeff_id( list(expected_desig_by_id.keys())  )
        
        # Check the results
        assert isinstance(result_dict ,dict)
        assert len(result_dict) == len(expected_desig_by_id)
        for id in result_dict:
            assert result_dict[id] == expected_desig_by_id[id],\
                f"{i},{desig}\nresult_dict={result_dict}"


def test_query_coefficients_by_jd_hp(db):
    """
    The coefficients_by_jd_hp function is designed to report back the coefficients for
    any objects in the given JD & HPlist
    
        returns
        -------
        dictionary-of-dictionaries
        - keys   = unpacked_primary_provisional_designation
        - values = list of coeff-dictionaries for each unpacked_primary_provisional_designation

    So we will need to pre-populate the tables with values that can give differing results
    depending on the search performed
    
    """

    # Instantiate & create all checker tables (creates two tables)
    C = sql.SQLChecker()
    C.create_all_checker_tables()
    
    # Define some data to be inserted ...
    desigs = [  '2020 AB',
                '2020 XY',
                '2021 PQ' ]
    JDlist = [2440000, 2440001, 2440002, 2440003]
    HPlist = [  [17,      17,      18,      18],
                [17,      18,      18,      19],
                [17,      18,      19,      19] ]

    # Loop to do the data insert (this is NOT the function being tested)
    save_dict = {}
    for i, desig in enumerate(desigs):
    
        # Continue to define the data to be inserted ...
        sector_field_names = C.generate_sector_field_names()[:2]
        sector_values_raw  = [ np.array([[i + 2*n,i + 2*n],[i + 2*n,i + 2*n]]) for n, _ in enumerate(sector_field_names) ]
        sector_values      = [ pickle.dumps( _, pickle.HIGHEST_PROTOCOL) for _ in sector_values_raw]

        # Call the convenient upsert functions (tested above) to insert the data into the object_coefficients & objects_by_jdhp tables
        object_coeff_id = C.upsert_coefficients(desig , sector_field_names, sector_values)
        C.insert_HP( JDlist, HPlist[i], object_coeff_id )  # <<-- E.g. HPlist[1] == [17,      18,      18,      19]

        # Save the inputs in a convenient dictionary to help with the query-verification below
        save_dict[desig] = { sfn:svr for sfn,svr in zip(sector_field_names, sector_values_raw) }



    # Query the data & verify the returned data is as expected
    # Loop over the possible days
    # For each day, I set out (by hand) the expected designations (and hence data) for different HP queries
    for n, JD in enumerate(JDlist):
    
        # -------------------1--------------------
        # Searching for objects in a single HP (17)
        # - In *expected_desigs* I define (by-hand) the expected desigs by day for the defined HPlist, and then use [n] for the JD-loop
        HPlist = [17]
        expected_desigs = [     [  '2020 AB','2020 XY','2021 PQ' ], # <<-- Expected for n = 0 [JD=2440000] & HP = 17
                                ['2020 AB'],                        # <<-- Expected for n = 1 [JD=2440001] & HP = 17
                                [],                                 # <<-- Expected for n = 2 [JD=2440002] & HP = 17
                                []                                  # <<-- Expected for n = 3 [JD=2440003] & HP = 17
                            ][n]
        
        #query
        # NB:
        #returns
        #-------
        #dictionary-of-dictionaries
        #- keys   = unpacked_primary_provisional_designation
        #- values = list of coeff-dictionaries for each unpacked_primary_provisional_designation
        result_dict = C.query_coefficients_by_jd_hp(JD, HPlist , sector_numbers = [0,1])
        
        #check
        # (a) same number of returns
        assert len(expected_desigs) == len(result_dict), f"expected_desigs={expected_desigs} , result_dict={result_dict}"
        # (b) same ids
        assert np.all( [ ed in result_dict for ed in expected_desigs ])
        # (c) same coefficients [complicated due to the structure of the dictionaries passed around ...]
        for ed in expected_desigs:
            expected_svr_dict = save_dict[ed]
            returned_svr_dict = result_dict[ed]
            for esfn,esvr in expected_svr_dict.items():
                rsvr = returned_svr_dict[esfn]
                assert np.all( esvr == rsvr )



        # -------------------2--------------------
        # Searching for objects in a list of HPs (17 & 18)
        # - In *expected_desigs* I define (by-hand) the expected desigs by day for the defined HPlist, and then use [n] for the JD-loop
        HPlist = [17, 18]
        expected_desigs = [
                    [ '2020 AB','2020 XY','2021 PQ' ],      # <<-- Expected for n = 0 [JD=2440000] & HP = [17, 18]
                    ['2020 AB','2020 XY','2021 PQ'],        # <<-- Expected for n = 1 [JD=2440001] & HP = [17, 18]
                    ['2020 AB','2020 XY'],                  # <<-- Expected for n = 2 [JD=2440002] & HP = [17, 18]
                    ['2020 AB']                             # <<-- Expected for n = 3 [JD=2440003] & HP = [17, 18]
                        ][n]
        
        #query
        # NB:
        #returns
        #-------
        #dictionary-of-dictionaries
        #- keys   = unpacked_primary_provisional_designation
        #- values = list of coeff-dictionaries for each unpacked_primary_provisional_designation
        result_dict = C.query_coefficients_by_jd_hp(JD, HPlist , sector_numbers = [0,1])


        #check
        # (a) same number of returns
        assert len(expected_desigs) == len(result_dict), f"expected_desigs={expected_desigs} , result_dict={result_dict}"
        # (b) same ids
        assert np.all( [ ed in result_dict for ed in expected_desigs ])
        # (c) same coefficients [complicated due to the structure of the dictionaries passed around ...]
        for ed in expected_desigs:
            expected_svr_dict = save_dict[ed]
            returned_svr_dict = result_dict[ed]
            for esfn,esvr in expected_svr_dict.items():
                rsvr = returned_svr_dict[esfn]
                assert np.all( esvr == rsvr )


        # -------------------3--------------------
        # Searching for objects in a list of HPs (17 & 19)
        # - In *expected_desigs* I define (by-hand) the expected desigs by day for the defined HPlist, and then use [n] for the JD-loop
        HPlist = [17, 19]
        expected_desigs = [ [ '2020 AB','2020 XY','2021 PQ' ], ['2020 AB'], ['2021 PQ'], ['2020 XY','2021 PQ' ] ][n]
        
        #query
        # NB:
        #returns
        #-------
        #dictionary-of-dictionaries
        #- keys   = unpacked_primary_provisional_designation
        #- values = list of coeff-dictionaries for each unpacked_primary_provisional_designation
        result_dict = C.query_coefficients_by_jd_hp(JD, HPlist , sector_numbers = [0,1])


        #check
        # (a) same number of returns
        assert len(expected_desigs) == len(result_dict), f"n={n}, JD={JD}, HPlist={HPlist}, expected_desigs={expected_desigs} , result_dict={result_dict}"
        # (b) same ids
        assert np.all( [ ed in result_dict for ed in expected_desigs ])
        # (c) same coefficients [complicated due to the structure of the dictionaries passed around ...]
        for ed in expected_desigs:
            expected_svr_dict = save_dict[ed]
            returned_svr_dict = result_dict[ed]
            for esfn,esvr in expected_svr_dict.items():
                rsvr = returned_svr_dict[esfn]
                assert np.all( esvr == rsvr )




    

"""
def test_queries():

    # (0) Create demo MSC(s) & create empty db using functions from test_orbit_cheby
    #  NB: test_loader_from_nbodysim() calls test_convenience_call_to_mpc_nbody_nbodysim() 
    MSCs = test_orbit_cheby.test_loader_from_nbodysim(test_orbit_cheby.orbfit_filenames[0])
    M = MSCs[0]
    
    # (1) Ensure db exists
    conn = test_func_create_db_and_tables()

    # (2) Low level direct inserts to put MSC into db so that we have something to query against!
    object_id = sql.insert_desig(conn ,MSCs[0].unpacked_primary_provisional_designation )
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
    object_id = sql.insert_desig(conn ,MSCs[0].unpacked_primary_provisional_designation )
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
