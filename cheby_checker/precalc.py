# -*- coding: utf-8 -*-
# mpchecker2/mpchecker2/precalc.py

'''
    --------------------------------------------------------------
    Precalculation module for checker
    
    May 2020
    Matt Payne
    
    This module provides functionalities to
    (a) Save precalculated orbit-chebyshev coefficients for an object
    (b) Save precalculated nightly-healpix locations for an object
    (c) Load data for (a) and (b) for an object
    
    Note that at present *no* attempt is made to ***ALWAYS RETURN NEOs***
     - This needs to be added-to / improved
     - We could (perhaps) return everything closer than 0.X au (on a given JD)
     - Or return everything faster than 0.Y au/day (on a given JD)
     
    --------------------------------------------------------------
    '''


# Import third-party packages
# --------------------------------------------------------------
import sys, os
import numpy as np
import operator
from collections import OrderedDict, defaultdict
from astropy_healpix import HEALPix
from functools import lru_cache
import json


# Import neighboring packages
# --------------------------------------------------------------
try:
    from orbit_cheby import orbit_cheby
except ImportError:
    from . import orbit_cheby
assert orbit_cheby.Base(), \
    'Seeing this text at evaluation indicates FAILURE to import orbit_cheby (as Base() should be available)'

try:
    from orbit_cheby import sql
except ImportError:
    from . import sql


try:
    from orbit_cheby import obs_pos
except ImportError:
    from . import obs_pos




# --------------------------------------------------------
# Primary External Class for Access to ChebyChecker's
#  Methods for Pre-calculating useful quantities
# --------------------------------------------------------

class PreCalc(orbit_cheby.Base , obs_pos.ObsPos):
    '''
        Primary External Class for accessing ChebyChecker's
        pre-calculated data

        Methods
        -------
        
        upsert()
        query_HPlist()
        
    '''
    def __init__(self):
        
        # Give access to "Base" & "ObsPos" methods & attributes
        orbit_cheby.Base.__init__(self)
        obs_pos.ObsPos.__init__(self)
        
        # connect to db
        self.conn = sql.create_connection(sql.fetch_db_filepath())

    def upsert(self, MSCs , geocenterXYZ ):#, integerJDs=None  ):
        '''
            Main method used to insert/update coefficients for an object
            Also handles the healpix calculations used for efficiency-of-read
            
            NB This function could be a lot cleaner if we calculated the geocenterXYZ internally
            
            Inputs:
            -------
            MSCs : List of Multi-Sector-Cheby class objects
            
            geocenterXYZ : np.array
             - Position of the geocenter (in Heliocentric Equatorial coords)
             - Needs to cover the span of the integer days defined in MSC.JDlist
            
            Returns:
            --------
            ????:
            -
            
        '''
        # ensure that the supplied variable is formatted correctly
        MSC_list = self._rectify_inputs(MSCs)
        
        # ensure that the supplied observatory coords have the right shape
        # N.B. Default list of Julian Dates to use, self.JDlist = (2440000 ==> 1968, 2460000.0 ==> 2023)
        assert geocenterXYZ.shape == (3 , len(self.JDlist) ), f'the shape of geocenterXYZ needs to be {(3, len(self.JDlist))} : instead it is {geocenterXYZ.shape}'
        
        # iterate over each MSC in list ...
        for M in MSC_list:
            
            # Insert the name & get back an object_id in return
            object_id = sql.insert_desig(self.conn , M.primary_unpacked_provisional_designation )
        
            # update list of coefficients
            sql.upsert_MSC(self.conn, M , object_id)
            
            # Use the coefficient-dictionary(ies) to get the HP for each integer-JD in JDlist
            # NB: need to restrict the queried dates to the those supported by the MSC
            indicees = np.where( self.JDlist < MSCs[0].get_valid_range_of_dates()[1] )[0]
            HPlist   = M.generate_HP(self.JDlist[indicees],  geocenterXYZ[:,indicees] , APPROX = True)

            # update HP data
            sql.insert_HP(self.conn, self.JDlist[indicees], HPlist, object_id)

    def get_nightly_precalcs(self,JD, HPlist):
        '''
            Main method used to query the precalculated healpix
            For a given JD, finds the object-integers for each HP in a list of HPs
            
            Returns the object-names and coefficient-dictionaries for
            the objects present in each healpix on the specified julian date
            
            Note that at present *no* attempt is made to ***ALWAYS RETURN NEOs***
            - This needs to be added-to / improved
            - We could (perhaps) return everything closer than 0.X au (on a given JD)

            inputs
            ------
            JD: float or int
             - julian date of the night. If <float> will be silently converted to <int>
             
            HPlist: list-of-integers 
             - healpix to be queried 
             - if integer (single healpix) supplied, is silently converted to list
             
            returns
            -------
            dictionary
            - key   = HP
            - value = list of coeff-dictionaries
        '''
        # Establish a connection
        conn = sql.create_connection( sql.fetch_db_filepath() )
        
        # Get the coefficients for the objects
        # - Outer dict is key-ed on object_id
        # - Inner dicts are key-ed on sector
        dict_of_dicts = sql.query_coefficients_by_jd_hp(conn, JD, HPlist , sector_numbers = None)

        # Get the designations for each of the object_ids returned
        desig_dict = sql.query_desig_by_number(conn, list(dict_of_dicts.keys()))
        
        # Swap the object_ids for the designations
        return orbit_cheby.MSC_Loader( FROM_DATABASE = True , desig_dict=desig_dict , dict_of_dicts_coeffs=dict_of_dicts ).MSCs
    
    
    def get_specific_object(self , primary_unpacked_provisional_designation , sector_numbers = None ):
        '''
            ...
        '''
        # Establish a connection
        conn = sql.create_connection( sql.fetch_db_filepath() )
    
        # Get the designation for the object_id (is in the form of a dictionary)
        desig_dict = sql.query_desig_by_number(conn, list(dict_of_dicts.keys()))
    
        # Query for the coefficients
        dict_of_coeffs = query_object_coefficients(conn,
                                  desig_dict[primary_unpacked_provisional_designation],
                                  sector_numbers = sector_numbers):
    
        # Return an MSC
        M = orbit_cheby.MSC()
        M.from_sector_dict(self, primary_unpacked_provisional_designation, dict_of_coeffs )
        return M

    def _rectify_inputs(self,  MSCs ):
        '''
            Private method called to rectify the input to upsert()
            
            inputs:
            -------
            MSCs: MSC or list-of-MSCs
            
            returns: 
            --------
            name_list: list-of-strings 
             - names of each object being "upcerted" 
             
            MSC_list: list-of-MSCs
             - 
            
        '''
        # If singular quantity, make into lists
        if  isinstance(MSCs, orbit_cheby.MSC ):
            MSC_list = [MSCs]

        # If non-singular, check plausibly formatted
        # (N.B. Not doing detailed content checks at this point)
        elif    isinstance(MSCs, (list, np.ndarray)) and \
                np.all([ isinstance(_, orbit_cheby.MSC) for _ in MSCs]):
            MSC_list = MSCs
        else:
            sys.exit('Cannot process inputs of type ... %r ' % type(MSCs) )

        # return everything in list form
        return MSC_list





