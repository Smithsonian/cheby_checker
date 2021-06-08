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
    from orbit_cheby import sql
    from orbit_cheby import obs_pos
    from orbit_cheby import mpc_nbody
except ImportError:
    from . import orbit_cheby
    from . import sql
    from . import obs_pos
    from . import mpc_nbody

assert orbit_cheby.Base(), \
    'Seeing this text at evaluation indicates FAILURE to import orbit_cheby (as Base() should be available)'



# --------------------------------------------------------
# Primary external class for creating & accessing
# various pre-calculated quantities
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
        self.db = sql.SQLChecker()
        self.conn = sql.create_connection(sql.fetch_db_filepath())
        
        
    def end_to_end_precalc(self,filenames, observatoryXYZ=None):
        '''
        Consider adding a high level function to handle ...
        (i) calling mpc_nbody on 1-or-many ORBFIT files
        (ii) calling MSCLoader on the results of (i)
        (iii) calling PreCalc.upsert() on the results of (ii)
        '''
        # Initiate NbodySim class with input files:
        # Run the integrator, by calling the object.
        Sim = mpc_nbody.NbodySim(   input       = filenames,
                                    filetype    = 'eq',
                                    save_parsed =   False,
                                    CHECK_EPOCHS=   True)
                                    )
        Sim(tstart=self.standard_MJDmin , tstep=20, trange=standard_MJDmax) # <<-- No justification for 20 days ...
        
        # Use the MSC_Loader to do all of the work to decalre and populate a list of MSC objects
        # NEEDS TO BE UPDATED TO USE Sim.output_times & Sim.output_vectors
        MSCs = orbit_cheby.MSC_Loader(FROM_ARRAY = True ,
                                        primary_unpacked_provisional_designations = name,
                                        times_TDB = times,
                                        statearray = states).MSCs
                                        
        # Do the precalculations and upsert
        # NB, In general we will *not* be passing observatory coords
        self.upsert( MSCs , observatoryXYZ=observatoryXYZ )


    def upsert(self, MSCs , observatoryXYZ=None ):
        '''
            Main method used to insert/update coefficients for an object
            Also handles the healpix calculations used for efficiency-of-read
                        
            Inputs:
            -------
            MSCs : List of Multi-Sector-Cheby class objects
            
            observatoryXYZ : np.array (optional)
             - Position of the observatory from which you want to calculate nightly healpix (in Heliocentric Equatorial coords)
             - Needs to cover the span of the integer days defined in MSC.JDlist
             - If not supplied, defaults to the geocenter
             - In general, we want to default to the geocenter
             - But for some spacecraft we *might* want to deviate from this, hence allowing the optional override
            
            Returns:
            --------
            ????:
            -
            
        '''
        # ensure that the supplied variable is formatted correctly
        MSC_list = self._rectify_inputs(MSCs)
        
        # If observatory coords are not supplied, set to be the location of the geocenter
        # N.B. Default list of Julian Dates to use, self.JDlist = (2440000 ==> 1968, 2460000.0 ==> 2023)
        if observatoryXYZ is not None:
            observatoryXYZ = self.get_heliocentric_equatorial_xyz(self.JDlist,
                                                                    obsCode="500",
                                                                    verbose=False)
                                        
        # Ensure that the observatory coords have the right shape
        assert observatoryXYZ.shape == (3 , len(self.JDlist) ), \
            f'the shape of observatoryXYZ needs to be {(3, len(self.JDlist))} : instead it is {observatoryXYZ.shape}'
        
        # iterate over each MSC in list ...
        for M in MSC_list:
            
            # Insert the name & get back an object_id in return
            object_id = self.db.insert_desig(M.primary_unpacked_provisional_designation)
        
            # update list of coefficients
            self.db.upsert_MSC(self.db.conn, M , object_id)
            
            # Use the coefficient-dictionary(ies) to get the HP for each integer-JD in JDlist
            # NB: need to restrict the queried dates to the those supported by the MSC
            indicees = np.where( self.JDlist < MSCs[0].get_valid_range_of_dates()[1] )[0]
            HPlist   = M.generate_HP(self.JDlist[indicees],  observatoryXYZ[:,indicees] , APPROX = True)

            # update HP data
            self.db.insert_HP(self.db.conn, self.JDlist[indicees], HPlist, object_id)

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
        
        # THE TWO QUERIES BELOW SHOULD PROBABLY BE COMBINED INTO A SINGLE SQL QUERY ...
    
        # Get the object_id for the supplied desig (is in the form of a dictionary)
        object_id = sql.query_number_by_desig(conn, list(dict_of_dicts.keys()))
    
        # Query for the coefficients
        dict_of_coeffs = sql.query_object_coefficients(conn,
                                                       object_id,
                                                       sector_numbers = sector_numbers)
        # Return dictionary
        return dict_of_coeffs

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





