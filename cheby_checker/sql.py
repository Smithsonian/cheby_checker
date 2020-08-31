# -*- coding: utf-8 -*-
# cheby_checker/cheby_checker/sql.py

'''
    --------------------------------------------------------------
    cheby_checker's sqlite module.
    
    Aug 2020
    Matt Payne
    
    This module provides functionalities to interact with
    (i) The tables related to storing ephemeris representations
        - I.e. chebyshev representations
    (ii) The table(s) related to storing the ITF data for
        - SIFTER
    
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

try:
    from orbit_cheby import orbit_cheby
except ImportError:
    from . import orbit_cheby
assert orbit_cheby.Base(), \
    'Seeing this text at evaluation indicates FAILURE to import orbit_cheby (as Base() should be available)'



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

class DB():

    def __init__(self,):
        self.db_file = self.fetch_db_filepath()
        self.conn    = self.create_connection()


    def fetch_db_filepath(self,):
        '''
        '''
        B = orbit_cheby.Base()
        db_dir = B._fetch_data_directory()
        return os.path.join(db_dir , B.db_filename)

    def create_connection(self,):
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
            conn = sqlite3.connect(self.db_file)
            return conn
        except Error as e:
            print(e)
        
        return conn


    def create_table(self, create_table_sql ):
        """ Create a table from the create_table_sql statement
            
            inputs:
            -------
            conn: Connection object
            
            create_table_sql: a CREATE TABLE statement
            
            return:
            -------
        """
        try:
            c = self.conn.cursor()
            c.execute(create_table_sql)
        except Error as e:
            print(e)


class SQLChecker(DB):

    def __init__(self,):
        super().__init__()


    # ---------------------------------------------
    # Generic functionalities
    # ---------------------------------------------

    def generate_sector_field_names(self,  sector_dict = orbit_cheby.Base().get_required_sector_dict() ):
        '''  Dynamically generate the field-specs that will be required for the coeffs-by-sector  '''
        return [ 'sector_%d_%d' % (i, jd) for i, jd in sector_dict.items() ]

    # ---------------------------------------------
    # Create the *CHECKER* tables
    # ---------------------------------------------
    def create_all_checker_tables():
        self.create_object_desig_table()
        self.create_object_coefficients_table()
        self.create_objects_by_jdhp_table()
        
    def create_object_desig_table(self,):
        """ Create the object_desig table that we need for name-to-integer mapping
            
            inputs:
            -------
            
            return:
            -------
            
            """
        
        
        # Create table ...
        
        sql_statement = """
            CREATE TABLE IF NOT EXISTS object_desig (
            object_id integer PRIMARY KEY,
            primary_unpacked_provisional_designation TEXT UNIQUE);
            """
        # create table
        if self.conn is not None:
            self.create_table( sql_statementc)
        
        # Create indicees
        createSecondaryIndex =  "CREATE INDEX index_desig ON object_desig (primary_unpacked_provisional_designation);"
        conn.cursor().execute(createSecondaryIndex)



    def create_object_coefficients_table(self,):
        """ Create the object_coefficients table(s) that we need for ephemeris calcultions
            
            inputs:
            -------
            
            return:
            -------
            
            """
        
        
        # Create table ...
        # Needs many fields, one per coefficient-set
        #  - Dynamically generate the field-specs that will be required for the coeffs-by-sector ...
        #  - This will look like ... sector_0_2440000 blob , sector_1_2440032 blob, ...
        sector_names = self.generate_sector_field_names()
        sector_spec  = " blob, ".join( sector_names )
        sector_spec  = sector_spec + " blob"
        
        sql_statement = """
            CREATE TABLE IF NOT EXISTS object_coefficients (
            coeff_id integer PRIMARY KEY,
            object_id integer UNIQUE, """ + \
                sector_spec + "); "
            
        # create table
        if self.conn is not None:
            self.create_table( sql_statement )


    def create_objects_by_jdhp_table(self,):
        """ Create the specific objects_by_jdhp table that we need for *mpchecker2*
            
            inputs:
            -------
            
            return:
            -------
            
            """

        # Create table ...
        sql_statement = """ CREATE TABLE IF NOT EXISTS objects_by_jdhp (
            jdhp_id integer PRIMARY KEY,
            jd integer NOT NULL,
            hp integer NOT NULL,
            object_id integer NOT NULL
            ); """
        
        
        # create table(s)
        if self.conn is not None:
            
            # create tables
            self.create_table( sql_statement )
            
            # Create indicees
            createSecondaryIndex =  "CREATE INDEX index_jdhp ON objects_by_jdhp (jd, hp);"
            conn.cursor().execute(createSecondaryIndex)
            createSecondaryIndex =  "CREATE INDEX index_pupd ON objects_by_jdhp (object_id);"
            conn.cursor().execute(createSecondaryIndex)


    # --------------------------------------------------------
    # --- Funcs to write to / update CHECKER db-tables
    # --------------------------------------------------------

    def upsert_MSC(self, M, object_id):
        """
            insert/update multi_sector_cheby object
            
            N.B ...
            https://stackoverflow.com/questions/198692/can-i-pickle-a-python-dictionary-into-a-sqlite3-text-field
            pdata = cPickle.dumps(data, cPickle.HIGHEST_PROTOCOL)
            curr.execute("insert into table (data) values (:data)", sqlite3.Binary(pdata))



            inputs:
            -------
            conn: Connection object
            
            M : MSC-object
             - see orbit_cheby module for detailed specification
             
            object_id : integer
             - The assumption is that this has been generated via the function, insert_desig()
             - I could explicitly query for the object_id via the function, query_number_by_desig(), but it seems unnecessary
            
            return:
            -------
            


        """
        # Sanity checks ?
        # assert ...
        
        # I guess that it will be quicker to do a single insert across all the required fields for a single object
        
        # (i) Get the sector field names required for this specific MSC
        # - Current method seems unnecessarily verbose
        sector_field_names =   self.generate_sector_field_names( sector_dict = \
                                                         { sector_num: orbit_cheby.Base().map_sector_number_to_sector_start_JD(sector_num , orbit_cheby.Base().standard_MJDmin) for sector_num in M.sector_coeffs.keys()}
                                                         )
        sector_field_names.append( 'object_id' )
        
        # (ii) Get the coefficients for each sector
        sector_field_values = [pickle.dumps(coeffs, pickle.HIGHEST_PROTOCOL) for coeffs in M.sector_coeffs.values()]
        sector_field_values.append(object_id)
        
        # (iii) Construct (in a horribly ungraceful manner) an sql insert statement
        sql =  " INSERT OR REPLACE INTO object_coefficients (" + ",".join(sector_field_names) + ") VALUES (" + ",".join(["?" for _ in sector_field_values]) + ");"

        # (iv) Execute the upsert ...
        cur = self.conn.cursor()
        cur.execute(sql, sector_field_values)
        conn.commit()


    def insert_HP(self, JDlist, HPlist, object_id ):
        """
            objects_by_jdhp is structured like ...
            id integer PRIMARY KEY,
            jd integer NOT NULL,
            hp integer NOT NULL,
            object_id integer NOT NULL
     
            inputs:
            -------
            conn: Connection object
            
            object_id : integer
                - The assumption is that this has been generated via the function, insert_desig()
                - I could explicitly query for the object_id via the function, query_number_by_desig(), but it seems unnecessary

        """
        cur = self.conn.cursor()

        # Sense-check
        assert len(JDlist)==len(HPlist), 'len(JDlist)!=len(HPlist) [%d != %d] in insert_HP' % (len(JDlist),len(HPlist))

        # *** Delete old entries ***
        # - My assumption here is that if we are inserting entries for a known object, we have likely calculated a new orbit
        # - Hence for sanitary reasons we should delete the old ones and replace with the new.
        self.delete_JDHP_by_object_id( object_id )
        
        # Insert new ...
        # (a) construct "records" variable which is apparently ammenable to single insert statement ...
        # https://pythonexamples.org/python-sqlite3-insert-multiple-records-into-table/
        # Note the use of "int" to force the numpy variables to play nicely with sql
        records = [ (int(jd), int(hp), object_id) for jd,hp in zip(JDlist,HPlist) ]
        
        # (b) construct sql string
        sqlstr = '''INSERT INTO objects_by_jdhp(jd,hp,object_id) VALUES(?,?,?);'''
        
        # (c) Insert
        cur.executemany(sqlstr, records)
        
        # (d) remember to commit ...
        self.conn.commit()


    def insert_desig(self, primary_unpacked_provisional_designation):
        """
            object_desig is structured like ...
            object_id integer PRIMARY KEY,
            primary_unpacked_provisional_designation TEXT UNIQUE);
            
        """
        # Do the insert
        cur = self.conn.cursor()
        sqlstr = "INSERT INTO object_desig (primary_unpacked_provisional_designation) VALUES ('%s') ON CONFLICT(primary_unpacked_provisional_designation) DO NOTHING" % primary_unpacked_provisional_designation
        cur.execute( sqlstr )
        self.conn.commit()

        # Get the object_id & return it
        # Might be able to replace this query with ... cursor.lastrowid, but not sure what happens in the do-nothing case ...
        cur.execute('''SELECT  object_id FROM object_desig WHERE primary_unpacked_provisional_designation=?;''', (primary_unpacked_provisional_designation,) )
        object_id = cur.fetchall()[0][0]
        return object_id



    # --------------------------------------------------------
    # --- Funcs to delete data
    # --------------------------------------------------------

    def delete_JDHP_by_object_id(self, object_id):
        """
            Delete all rows from "objects_by_jdhp" that match the supplied "primary_unpacked_provisional_designation"
        """
        cur = self.conn.cursor()

        # Construct & execute the sql query
        # - This is matching/joining on object-id# and then deleting only from objects_by_jdhp
        #   (and leaving the entry in object_desig)
        cur.execute( " DELETE FROM objects_by_jdhp WHERE object_id=?;", (int(object_id),) )
        self.conn.commit()



    # --------------------------------------------------------
    # --- Funcs to query db-tables
    # --------------------------------------------------------

    def query_object_coefficients(self,
                                  object_id,
                                  sector_numbers = None):
        """
           Define standard query used to get cheby-coeff data for a named object
           Can optionally select only a subset of sectors
           
           inputs:
           -------
           object_id: integer
            - The assumption is that this has been generated via the function, insert_desig()
           sector_number: integer
            -
           
           return:
           -------
           list of numpy-arrays
            - Each numpy-array item is a set of cheby-coeffs for a specific sector

        """
        cur = self.conn.cursor()
        
        # What sector numbers are we searching for ?
        # - Default is to get data for all of them
        if sector_numbers is None :
            sector_field_names = self.generate_sector_field_names()
        else:
            sector_field_names = self.generate_sector_field_names( sector_dict = {
                                                                sector_num: sector_JD for \
                                                                sector_num, sector_JD in zip(sector_numbers ,
                                                                                          orbit_cheby.Base().map_sector_number_to_sector_start_JD(np.atleast_1d(sector_numbers) ,\
                                                                                                                                                  orbit_cheby.Base().standard_MJDmin))})
        # Construct & execute the sql query
        sqlstr = "SELECT " + ", ".join( sector_field_names ) + " FROM object_coefficients WHERE object_id=?"
        cur.execute(sqlstr , ( object_id, ))

        # Parse the result ...
        result = cur.fetchall()[0]
        return { sfn:pickle.loads( coeff )  for sfn, coeff in zip(sector_field_names, result) if coeff != None }



    def query_number_by_desig(self, primary_unpacked_provisional_designation):
        """
            Given a desig, return a number
            
            May be of little use in practice, but helpful for development
        """
        # Get the object_id & return it
        cur = self.conn.cursor()
        cur.execute('''SELECT object_id FROM object_desig WHERE primary_unpacked_provisional_designation=?;''', (primary_unpacked_provisional_designation,) )
        object_id = cur.fetchall()[0][0]
        return object_id

    def query_desig_by_number(self, object_ids):
        """
        Given a list of object_ids, return the associated primary_unpacked_provisional_designations

        inputs:
        -------
        conn: Connection object

        object_ids: list of integer object_ids

        return:
        -------
        dictionary
         - maps object_id (key) to designation (value)

        """
        # Get the object_id & return it
        cur = self.conn.cursor()
        cur.execute(f'''SELECT object_ID,primary_unpacked_provisional_designation FROM object_desig WHERE object_ID in ({','.join(['?']*len(object_ids))});''', (*(int(_) for _ in object_ids),) )
        # key is object _id, value is desig
        return {_[0]:_[1] for _ in cur.fetchall()}



    def query_coefficients_by_jd_hp(self, JD, HPlist , sector_numbers = None):
        '''
            For a given (single) JD and list of Healpix,
            the query returns the relevant coefficients from the object_coefficients table
            
        '''
        
        # What sector numbers are we searching for ?
        # - Default is to get data for all of them
        if sector_numbers is None :
            sector_field_names = self.generate_sector_field_names()
        else:
            sector_field_names = self.generate_sector_field_names( sector_dict = {
                                                             sector_num: sector_JD for \
                                                             sector_num, sector_JD in zip(sector_numbers ,
                                                                                          orbit_cheby.Base().map_sector_number_to_sector_start_JD(np.atleast_1d(sector_numbers) ,\
                                                                                                                                                  orbit_cheby.Base().standard_MJDmin))})
        # Construct & execute the sql query
        sqlstr =    "SELECT object_coefficients.object_id," + \
                    ", ".join( sector_field_names ) + \
                    f" FROM object_coefficients INNER JOIN objects_by_jdhp ON objects_by_jdhp.object_id = object_coefficients.object_id WHERE objects_by_jdhp.jd=? and objects_by_jdhp.hp in ({','.join(['?']*len(HPlist))})"
        cur = conn.cursor()
        cur.execute(sqlstr , (int(JD), *(int(_) for _ in HPlist), ))

        # Parse the result into a dict-of-dicts
        results      = cur.fetchall()
        return { r[0] : { sfn:pickle.loads( coeff )  for sfn, coeff in zip(sector_field_names, r[1:]) if coeff != None } for r in results}



        # INNER JOIN
        # DISTINCT
        # WHERE
        cur = conn.cursor()

        # This will get the unique object-IDs
        "SELECT DISTINCT object_id FROM objects_by_jdhp WHERE jd=? and hp in ?"


    def query_jd_hp(self, JD, HPlist):
        '''
            For a given (single) JD and list of Healpix,
            the query returns the relevant object_ids
            
            May be of little use in practice, but helpful for development
            
        '''
        # Connection cursor
        cur = self.conn.cursor()

        # This will get the unique object-IDs
        # I hate this fucking query for many reasons ...
        # (i) Why the fuck do I need to force any input numpy-integers, to be integers ???
        #(ii) Why the fuck do I need to explicitly expand the number of "?" in the "in ()" statement ???
        cur.execute(f"SELECT object_id FROM objects_by_jdhp WHERE jd=? and hp in ({','.join(['?']*len(HPlist))});", (int(JD), *(int(_) for _ in HPlist), ) )
        #object_ids = [_[0] for _ in cur.fetchall()]
        return [_[0] for _ in cur.fetchall()]





class SQLSifter(DB):
    
    def __init__():
        super().__init__()


    # ---------------------------------------------
    # Generic functionalities
    # ---------------------------------------------

    # ---------------------------------------------
    # Create the *CHECKER* tables
    # ---------------------------------------------
    def create_sifter_tables(conn):
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
        
        # remember to commit ...
        conn.commit()


    # --------------------------------------------------------
    # --- Funcs to write to / update SIFTER db-tables
    # --------------------------------------------------------
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
        
        # Make cursor ...
        cur = conn.cursor()
        
        # Construct sql for tracklet insert ...
        sql =  ''' INSERT OR REPLACE INTO tracklets(jd,hp,tracklet_name,tracklet)
            VALUES(?,?,?,?)
            '''
        
        # Insert
        cur.execute(sql, (jd, hp, tracklet_name, sqlite3.Binary( pickle.dumps(tracklet_dict, pickle.HIGHEST_PROTOCOL) ),))

        # remember to commit ...
        conn.commit()


    def upsert_tracklets(conn, jd_list, hp_list, tracklet_name_list, tracklet_dict_list):
        """
            insert/update lists of tracklet data
            
            N.B ...
            https://stackoverflow.com/questions/198692/can-i-pickle-a-python-dictionary-into-a-sqlite3-text-field
            pdata = cPickle.dumps(data, cPickle.HIGHEST_PROTOCOL)
            curr.execute("insert into table (data) values (:data)", sqlite3.Binary(pdata))
            
            inputs:
            -------
            conn: Connection object
            
            jd_list :
            hp_list :
            tracklet_name_list :
            tracklet_dict_list :
            
            return:
            -------
            
            """
        # Make cursor ...
        cur = conn.cursor()

        # construct "records" variable which is apparently ammenable to single insert statement ...
        # https://pythonexamples.org/python-sqlite3-insert-multiple-records-into-table/
        records = [ (int(jd), int(hp), tracklet_name, sqlite3.Binary(pickle.dumps(tracklet_dict, pickle.HIGHEST_PROTOCOL))) \
                   for jd, hp, tracklet_name, tracklet_dict \
                   in zip(jd_list, hp_list, tracklet_name_list, tracklet_dict_list) ]

        sql = '''INSERT OR REPLACE INTO tracklets(jd,hp,tracklet_name,tracklet) VALUES(?,?,?,?);'''

        # Insert
        cur.executemany(sql, records)

        # remember to commit ...
        conn.commit()


    # --------------------------------------------------------
    # --- Funcs to delete from SIFTER db-tables
    # --------------------------------------------------------
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


    def delete_tracklets(conn, tracklet_name_list):
        """
            delete list of tracklet data
            
            inputs:
            -------
            tracklet_name: list-of-strings
            
            return:
            -------
            
            
            """
        sql = '''DELETE FROM tracklets WHERE tracklet_name IN (?);'''
        records = [ (tracklet_name,) for tracklet_name in tracklet_name_list ]
        cur = conn.cursor()
        cur.executemany(sql, records)
        conn.commit()


    # --------------------------------------------------------
    # --- Funcs to query SIFTER db-tables
    # --------------------------------------------------------
    def query_tracklets_jdhp(conn, JD, HP):
        """
           Standard query used to find all tracklets for which jd,hp matches input
           
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
        
        # return a list-of-tuples: (tracklet_name, tracklet_dictionary)
        return [ (row[0] , pickle.loads( row[1] ) ) for row in cur.fetchall() ]


    def query_tracklets_jd_hplist(conn, JD, HP_list):
        """
            Standard query used to find all tracklets for which jd,hp_list matches input
            
            inputs:
            -------
            JD: integer
            HP_list: list-of-integers
            
            return:
            -------
            list of tracklet_names
            
            """
        assert isinstance(JD, int) and isinstance(HP_list, list), 'Cannot parse input types in query_tracklets_jd_hplist ... '

        # Get cursor ...
        cur = conn.cursor()

        # Set up query:
        #N.B. https://stackoverflow.com/questions/18363276/how-do-you-do-an-in-query-that-has-multiple-columns-in-sqlite
        for sql in ['''CREATE TEMPORARY TABLE lookup(jd, hp);''',
                    '''INSERT OR REPLACE INTO lookup(jd,hp) VALUES(?,?);''',
                    '''CREATE INDEX lookup_jd_hp ON lookup(jd, hp);''',
                    '''SELECT tracklet_name, tracklet FROM tracklets JOIN lookup ON tracklets.jd = lookup.jd AND tracklets.hp = lookup.hp;''']:
            
            # Execute queries
            print( ' ... = ', sql )
            if 'INSERT' in sql:
                records = [ (JD, hp) for hp in HP_list ]
                cur.executemany(sql, records )
            else:
                cur.execute(sql)

        # return a list-of-tuples: (tracklet_name, tracklet_dictionary)
        return [ (row[0] , pickle.loads( row[1] ) ) for row in cur.fetchall() ]





