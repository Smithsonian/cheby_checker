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


# Import neighbouring packages
# -----------------------------------------------------------------------------
from cheby_checker     import precalc
from cheby_checker     import sql

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
    P.end_to_end_precalc( [mpc_orb_json_filepath] )
    
    # Check that data has been inserted into the database
    # - NB at this point we are NOT testing for accuracy
