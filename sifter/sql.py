# -*- coding: utf-8 -*-
# sifter/sifter/sql

'''
    --------------------------------------------------------------
    sifter's sqlite module.
    
    Jan 2020
    Matt Payne & Mike Alexandersen
    
    This module provides functionalities to
    ...
    
    *WRITE MORE STUFF*
    
    --------------------------------------------------------------
    '''


# Import third-party packages
# --------------------------------------------------------------
import sys, os
import numpy as np
import sqlite3
from sqlite3 import Error
import pickle

# Import neighboring packages
# --------------------------------------------------------------
print(os.path.dirname(os.path.realpath(__file__)))
sys.path.append( os.path.dirname(os.path.realpath(__file__)) )
import precalc




# Data classes/methods
#
# N.B. Many of these sqlite functions are copied straight from
# https://www.sqlitetutorial.net/sqlite-python/
# E.g.
# https://www.sqlitetutorial.net/sqlite-python/creating-database/
# https://www.sqlitetutorial.net/sqlite-python/create-tables/
# ...
#
# -------------------------------------------------------------

def fetch_db_filepath():
    '''
    '''
    B = precalc.Base()
    db_dir = B._fetch_data_directory()
    return os.path.join(db_dir , B.db_filename)

def create_connection(db_file):
    """ Create a database connection to the SQLite database
        specified by db_file
        
        inputs:
        -------
        db_file: database file
        
        return: 
        -------
        Connection object or None
    """
    conn = None
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except Error as e:
        print(e)
    
    return conn


def create_table(conn, create_table_sql):
    """ Create a table from the create_table_sql statement
        
        inputs:
        -------
        conn: Connection object
        
        create_table_sql: a CREATE TABLE statement
        
        return:
        -------
    """
    try:
        c = conn.cursor()
        c.execute(create_table_sql)
    except Error as e:
        print(e)



def create_specific_table(conn):
    """ Create the specific table(s) that we need for *sifter*
        Currently creates:
        (i) tracklets

        inputs:
        -------
        conn: Connection object

        return:
        -------

    """
    
    # Note that I am deliberately setting up the stored *tracklet* data as a "blob"
    # - Because not yet sure what will be in it!
    sql_create_tracklets_table = """ CREATE TABLE IF NOT EXISTS tracklets (
        id integer PRIMARY KEY,
        jd integer NOT NULL,
        hp integer NOT NULL,
        tracklet_name text NOT NULL,
        tracklet blob
        ); """
    
    # create table(s)
    if conn is not None:
        
        # create "jdhp_to_tracklet" table
        create_table(conn, sql_create_tracklets_table)

        # Create compound/combined index on the jd & hp columns
        createSecondaryIndex =  "CREATE INDEX index_jdhp ON tracklets (jd, hp);"
        conn.cursor().execute(createSecondaryIndex)




def upsert_tracklet(conn, jd, hp, tracklet_name, tracklet_dict):
    """
        insert/update tracklet data
        
        N.B ...
        https://stackoverflow.com/questions/198692/can-i-pickle-a-python-dictionary-into-a-sqlite3-text-field
        pdata = cPickle.dumps(data, cPickle.HIGHEST_PROTOCOL)
        curr.execute("insert into table (data) values (:data)", sqlite3.Binary(pdata))

        inputs:
        -------
        conn: Connection object
        
        jd : 
        hp : 
        tracklet_name : 
        tracklet_dict :
        
        return:
        -------
        
        
    """
    
    pdata = pickle.dumps(tracklet_dict, pickle.HIGHEST_PROTOCOL)

    sql =  ''' INSERT OR REPLACE INTO tracklets(jd,hp,tracklet_name,tracklet)
        VALUES(?,?,?,?)
        '''
    cur = conn.cursor()
    cur.execute(sql, (jd, hp, tracklet_name, sqlite3.Binary(pdata),))
    conn.commit()


def delete_tracklet(conn, tracklet_name):
    """
        delete tracklet data
        
        inputs:
        -------
        tracklet_name: string
        
        return:
        -------
        
        
    """
    sql = 'DELETE FROM tracklets WHERE tracklet_name=?'
    cur = conn.cursor()
    cur.execute(sql, (tracklet_name,))
    conn.commit()




def query_tracklets_jdhp(conn, JD, HP):
    """
       Define standard query used to find all tracklets for which jdhp == input
       
       inputs:
       -------
       JD: integer
       HP: integer
       
       return:
       -------
       list of tracklet_names

    """
    cur = conn.cursor()
    cur.execute("SELECT tracklet_name, tracklet FROM tracklets WHERE jd=? AND hp=?", ( int(JD), int(HP) , ))
    
    # return a dictionary-of-dictionaries
    return { row[0]: pickle.loads( row[1] ) for row in cur.fetchall() }

