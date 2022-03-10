# -*- coding: utf-8 -*-
# /tests/test_nbody.py

'''
----------------------------------------------------------------------------
tests for cheby_checker/precalc
 - Here I focus on the convenience funtions that allow us to perform end-to-end-precalculations after an orit-fit has been run by orbfit
 
Feb 2022
Matthew Payne

----------------------------------------------------------------------------
'''

# import third-party packages
# -----------------------------------------------------------------------------
import sys
import os
import numpy as np
import pytest
from astropy.time import Time
import pytest
from filecmp import cmp
import getpass
import json
import time
import glob
import sqlite3
import pickle

# Import MPC packages
# -----------------------------------------------------------------------------
from mpc_orb.parse import MPCORB

# Import neighbouring packages
# -----------------------------------------------------------------------------
from cheby_checker     import precalc
from cheby_checker     import sql
from cheby_checker     import cheby_checker

this_dir = os.path.abspath(os.path.dirname( __file__ ))
repo_dir = os.path.abspath(os.path.dirname( this_dir ))
data_dir = os.path.join(repo_dir, 'dev_data')
json_dir = os.path.join(data_dir, 'mpc_orb_jsons')
std_json_dir = os.path.join(json_dir, 'standard_mp') # Standard grav-only fits
mpc_orb_json_files = glob.glob(std_json_dir + '/*.json' )


# Utility functions to help with testing
# -----------------------------------------------------------------------------

def similar_xyzuvw(xyzv0, xyzv1, threshold_xyz=1e-13, threshold_v=1e-14): # 15 mm, 1.5 mm/day
    '''
    Calculate the difference between two sets of cartesian coordinates.
    '''
    if isinstance(xyzv0, list):
        xyzv0 = np.array(xyzv0)
    if isinstance(xyzv1, list):
        xyzv1 = np.array(xyzv1)
    error = xyzv0 - xyzv1
    good_tf = np.abs(error) < np.array([threshold_xyz] * 3 + [threshold_v] * 3)
    return np.all(good_tf), error

def db_handler(func):
    '''
    Using decorators to do the db set-up ...
    
    NB: This is a similar approach to test_sql.py

     - This deletes any extant db before the test
     - Populates the database with empty tables
     - And deletes any extant db after the test
     
    '''
    
    def inner_function(*args, **kwargs):
    
        # Delete db if any exists for any reason
        filepath = sql.DB.fetch_db_filepath()
        if os.path.isfile(filepath):
            os.remove(filepath)
            
        # Create the required empty tables
        # - See test_sql for tests of this functionality
        sql.SQLChecker().create_all_checker_tables()
            
        # Run test function
        result = func(*args, **kwargs)

        # Delete the db
        if os.path.isfile(filepath):
            os.remove(filepath)
        
        # Return any result
        return result

    return inner_function




# Tests of PreCalc
# -----------------------------------------------------------------------------

@db_handler
def test_PreCalc_A():
    '''
    Test instantiation of PreCalc object
    '''
    # Instantiate
    P = precalc.PreCalc()

    # Check that the expected attributes from Base exist
    assert \
        hasattr(P, 'standard_MJDmin') and\
        hasattr(P, 'standard_MJDmax')

    # Check that the expected attributes from ObsPos exist
    assert \
        hasattr(P, 'obsCodes')

    # Check that the expected attributes from SQLChecker exist
    assert \
        hasattr(P, 'db_file') and\
        hasattr(P, 'conn')

    # Check that we can get a sqlite cursor
    cur  = P.conn.cursor()
    assert isinstance( cur, sqlite3.Cursor )

    # Satisfy ourselves that the @db_handler is working and has created a database with
    # the required empty tables ()
    # - if the count is 1, then table exists
    cur = P.conn.cursor()
    cur.execute('SELECT name from sqlite_master WHERE type = "table" AND name = "objects_by_jdhp"')
    assert len(cur.fetchone()) == 1 , 'jdhp table does not exist'
    cur.execute('SELECT name from sqlite_master WHERE type = "table" AND name = "object_coefficients"')
    assert len(cur.fetchone()) == 1 , 'coeff table does not exist'




@db_handler
def test_end_to_end_precalc_A():
    '''
    Test end_to_end_precalc
    
    A high level function to handle ...
    (i) calling nbody on 1-or-many ORBFIT files
    (ii) calling MSCLoader on the results of (i)
    (iii) calling PreCalc.upsert() on the results of (ii)
    
    
    '''

    # Instantiate PreCalc
    P = precalc.PreCalc()

    # Specify input orbfit file
    mpc_orb_json_filepath = mpc_orb_json_files[0]

    # Call the *end_to_end_precalc()* func that we want to test
    # ---------------------------------------------------
    P.end_to_end_precalc( [mpc_orb_json_filepath] )
    # ---------------------------------------------------

    # Check that data has been inserted into the database
    # - NB at this point we are NOT testing for accuracy
    # ---------------------------------------------------
    # Use MPCORB to get the expected designation
    M = MPCORB(mpc_orb_json_filepath)
    desig = M.designation_data["unpacked_primary_provisional_designation"]

    # Query the object_coefficients table in the db and see what the returned data looks like
    # ---------------------------------------------------
    C   = sql.SQLChecker()
    cur = C.conn.cursor()
    sqlstr = "SELECT * FROM object_coefficients WHERE unpacked_primary_provisional_designation=?"
    cur.execute(sqlstr , ( desig, ))
    result_raw              = cur.fetchall()[0]
    result_object_coeff_id  = result_raw[0]
    result_desig            = result_raw[1]
    result_sectors          = [ pickle.loads( _ ) for _ in result_raw if _ != None and not isinstance(_,(str,int))]
    
    # Check the designation is as expected
    assert desig == result_desig
    
    # Check that the expected number of sectors is returned
    array_of_days = cheby_checker.Base().JDlist
    count_per_sector = {}
    for s in cheby_checker.Base().map_JD_to_sector_number( array_of_days , cheby_checker.Base().standard_MJDmin):
        n = 1 if s not in count_per_sector else count_per_sector[s] + 1
        count_per_sector[s] = n
    expected_supported_sectors = { k:v for k,v in count_per_sector.items() if v > 4 }
    assert len(result_sectors) == len(expected_supported_sectors)
    
    # Now query the HEALPIX table
    # ---------------------------------------------------
    # Query the db and examine the results
    sqlstr = "SELECT * FROM objects_by_jdhp WHERE object_coeff_id=?"
    cur = C.conn.cursor()
    cur.execute(sqlstr , ( result_object_coeff_id, ))  #<<-- result_object_coeff_id from coefficient query above
    result_raw     = cur.fetchall()

    # Check the inputs are the same as the outputs
    # NB: Not testing the healpix here...
    for n,r in enumerate(result_raw):
        id,jd,hp,obj_co_id = r
        assert jd==array_of_days[n]
        assert obj_co_id==result_object_coeff_id
        
    print("Passed: test_end_to_end_precalc_A" )

@db_handler
def test_end_to_end_precalc_B():
    '''
    Test end_to_end_precalc
    
    As per test_end_to_end_precalc_B (above), but
    now we test multiple single integrations & reads

    
    '''
    
    # Specify input orbfit file
    for mpc_orb_json_filepath in mpc_orb_json_files:
        print( 'test_end_to_end_precalc_B\n\t', mpc_orb_json_filepath)

        # Instantiate PreCalc
        P = precalc.PreCalc()

        # Call the *end_to_end_precalc()* func that we want to test
        # ---------------------------------------------------
        P.end_to_end_precalc( [mpc_orb_json_filepath] )
        # ---------------------------------------------------

        # Check that data has been inserted into the database
        # - NB at this point we are NOT testing for accuracy
        # ---------------------------------------------------
        # Use MPCORB to get the expected designation
        M = MPCORB(mpc_orb_json_filepath)
        desig = M.designation_data["unpacked_primary_provisional_designation"]

        # Query the object_coefficients table in the db and see what the returned data looks like
        # ---------------------------------------------------
        C   = sql.SQLChecker()
        cur = C.conn.cursor()
        sqlstr = "SELECT * FROM object_coefficients WHERE unpacked_primary_provisional_designation=?"
        cur.execute(sqlstr , ( desig, ))
        result_raw              = cur.fetchall()[0]
        result_object_coeff_id  = result_raw[0]
        result_desig            = result_raw[1]
        result_sectors          = [ pickle.loads( _ ) for _ in result_raw if _ != None and not isinstance(_,(str,int))]
        
        # Check the designation is as expected
        assert desig == result_desig
        
        # Check that the expected number of sectors is returned
        array_of_days = cheby_checker.Base().JDlist
        count_per_sector = {}
        for s in cheby_checker.Base().map_JD_to_sector_number( array_of_days , cheby_checker.Base().standard_MJDmin):
            n = 1 if s not in count_per_sector else count_per_sector[s] + 1
            count_per_sector[s] = n
        expected_supported_sectors = { k:v for k,v in count_per_sector.items() if v > 4 }
        assert len(result_sectors) == len(expected_supported_sectors)
        
        # Now query the HEALPIX table
        # ---------------------------------------------------
        # Query the db and examine the results
        sqlstr = "SELECT * FROM objects_by_jdhp WHERE object_coeff_id=?"
        cur = C.conn.cursor()
        cur.execute(sqlstr , ( result_object_coeff_id, ))  #<<-- result_object_coeff_id from coefficient query above
        result_raw     = cur.fetchall()

        # Check the inputs are the same as the outputs
        # NB: Not testing the healpix here...
        for n,r in enumerate(result_raw):
            id,jd,hp,obj_co_id = r
            assert jd==array_of_days[n]
            assert obj_co_id==result_object_coeff_id
            
    print("Passed: test_end_to_end_precalc_B" )

    
"""    *** DELIBERATELY COMMENTING OUT THIS TEST [MJP: 2022-03-10] DUE TO PROBLEMS WITH MULTI-PARTICLE INTEGRATIONS ***
@db_handler
def test_end_to_end_precalc_C():
    '''
    Test end_to_end_precalc
    
    As per test_end_to_end_precalc_B (above), but
    now we test multiple integration & reads
    
    As part of this, we move try and pass in the
    observatoryXYZ positions for the healpix calculation
 
        *** NOTES ON PROBLEMS WITH MULTI-PARTICLE INPUT ***
        When testing on 4-particle input, I found that the integration ran wild,
        with the solutions zooming off to 10,000au (individually they stayed @ 3au)
        The problem was found to be that the input matrix, self.bary_eq_vec, (2D array)
        gets transposed inside C, but only when passed in as-is.
         - A hard-pasted copy of the array gets passed in just fine, without transposition
        Checking things in python before shows *NO* difference between the self.bary_eq_vec
        and the hardpasted copy, so the difference is presumably something to do with
        the internal memory array adopted in numpy arrays in different circumstances, and then
        how that is translated to C
        SPECULATION : Is it something to do with bary_eq_vec being reshaped / reduced in dimension
        at some point in the preceeding code? Does this affect memory layout?

        hardpasted = np.array(
        [[ 2.4260322844717380e-01, -3.0090867545871860e+00, -1.1363474792630044e+00,
           9.1908324088463677e-03,  2.0472935970702337e-04, -1.5374992052716444e-03],
         [-2.7888126437563829e+00,  1.7058528476705044e+00,  3.9449682717556311e-01,
          -5.2133449765416028e-03, -7.5295451342754713e-03, -7.8760899838941603e-04],
         [-5.1094716747791424e-01, -3.1254423057577734e+00, -2.0995426392478924e+00,
           7.8191160668487376e-03, -3.3055503238910031e-04,  9.7078064274683739e-04],
         [-2.1582274449687930e+00,  2.0100970703192091e+00,  3.1933652106784710e-01,
          -5.5592970943310026e-03, -7.8839265248180895e-03, -3.4116741888535933e-03]])

        The Good Inputs
        INPUT SET TO HARDPASTED NUMBERS ...

        *** f= *** 0.7916875000004463

         ---ccc--- 0.2426032284471738 -3.0090867545871860 -1.1363474792630044 0.0091908324088464 0.0002047293597070 -0.0015374992052716

         ---ccc--- -2.7888126437563829 1.7058528476705044 0.3944968271755631 -0.0052133449765416 -0.0075295451342755 -0.0007876089983894

         ---ccc--- -0.5109471674779142 -3.1254423057577734 -2.0995426392478924 0.0078191160668487 -0.0003305550323891 0.0009707806427468

         ---ccc--- -2.1582274449687930 2.0100970703192091 0.3193365210678471 -0.0055592970943310 -0.0078839265248181 -0.0034116741888536



        The Bad Inputs
        INPUTS AS ORIGINAL ...

        *** f= *** 0.7916875000004463

         ---ccc--- 0.2426032284471738 -2.7888126437563829 -0.5109471674779142 -2.1582274449687930 -3.0090867545871860 1.7058528476705044

         ---ccc--- -3.1254423057577734 2.0100970703192091 -1.1363474792630044 0.3944968271755631 -2.0995426392478924 0.3193365210678471

         ---ccc--- 0.0091908324088464 -0.0052133449765416 0.0078191160668487 -0.0055592970943310 0.0002047293597070 -0.0075295451342755

         ---ccc--- -0.0003305550323891 -0.0078839265248181 -0.0015374992052716 -0.0007876089983894 0.0009707806427468 -0.0034116741888536

    
    '''

    # Instantiate PreCalc
    P = precalc.PreCalc()

    # Specify input orbfit file
    mpc_orb_json_filepaths = mpc_orb_json_files

    # Call the *end_to_end_precalc()* func that we want to test
    # ---------------------------------------------------
    P.end_to_end_precalc( mpc_orb_json_filepaths )
    # ---------------------------------------------------

    # Check that data has been inserted into the database
    # - NB at this point we are NOT testing for accuracy
    # ---------------------------------------------------
    # Use MPCORB to get the expected designations
    desig_list= [ MPCORB(filepath).designation_data["unpacked_primary_provisional_designation"] for filepath in mpc_orb_json_filepaths ]

    
    # For each object that we expect to be in the database...
    # query the object_coefficients table in the db and see what the returned data looks like
    C   = sql.SQLChecker()
    cur = C.conn.cursor()
    for desig in desig_list:
    

        # Query the object_coefficients table in the db and see what the returned data looks like
        # ---------------------------------------------------
        sqlstr = "SELECT * FROM object_coefficients WHERE unpacked_primary_provisional_designation=?"
        cur.execute(sqlstr , ( desig, ))
        result_raw              = cur.fetchall()[0]
        result_object_coeff_id  = result_raw[0]
        result_desig            = result_raw[1]
        result_sectors          = [ pickle.loads( _ ) for _ in result_raw if _ != None and not isinstance(_,(str,int))]
    
        # Check the designation is as expected
        assert desig == result_desig
    
        # Check that the expected number of sectors is returned
        array_of_days = cheby_checker.Base().JDlist
        count_per_sector = {}
        for s in cheby_checker.Base().map_JD_to_sector_number( array_of_days , cheby_checker.Base().standard_MJDmin):
            n = 1 if s not in count_per_sector else count_per_sector[s] + 1
            count_per_sector[s] = n
        expected_supported_sectors = { k:v for k,v in count_per_sector.items() if v > 4 }
        assert len(result_sectors) == len(expected_supported_sectors)


        # Now query the HEALPIX table
        # ---------------------------------------------------
        # Query the db and examine the results
        sqlstr = "SELECT * FROM objects_by_jdhp WHERE object_coeff_id=?"
        cur = C.conn.cursor()
        cur.execute(sqlstr , ( result_object_coeff_id, ))  #<<-- result_object_coeff_id from coefficient query above
        result_raw     = cur.fetchall()

        # Check the inputs are the same as the outputs
        # NB: Not testing the healpix here...
        for n,r in enumerate(result_raw):
            id,jd,hp,obj_co_id = r
            assert jd==array_of_days[n]
            assert obj_co_id==result_object_coeff_id
"""
    

test_end_to_end_precalc_A()
