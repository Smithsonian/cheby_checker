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
        Base.__init__(self)
        ObsPos.__init__(self)
        
        # connect to db
        self.conn = sql.create_connection(sql.fetch_db_filepath())

    def upsert(self, MSCs ):
        '''
            Main method used to insert/update coefficients for an object
            Also handles the healpix calculations used for efficiency-of-read
            
            Inputs:
            -------
            arg: dictionary or list-of-dictionaries
            - dictionaries must be a valid multi_sector_cheby_dict
            - see orbit_cheby module for detailed specification of multi_sector_cheby_dict
            
            Returns:
            --------
            ????:
            -
            
        '''
        # ensure that the supplied variable is formatted correctly
        MSC_list = self._rectify_inputs(MSCs)
        
        # iterate over each MSC in list ...
        for M in MSC_list:
        
            # update list of coefficients
            sql.upsert_MSC(self.conn, M )
            
            # Use the coefficient-dictionary(ies) to get the HP for each integer-JD in JDlist
            observatoryXYZ = ObsPos.get_heliocentric_equatorial_xyz(self.JDlist , obsCode='500')
            HPlist = M.generate_HP(self.JDlist,     observatoryXYZ)

            # update HP data
            self.upsert_HP(M.primary_unpacked_provisional_designation, self.JDlist, HPlist)

    def get_nightly_precalcs(JD, HPlist):
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
        HPlist = [HPlist] if isinstance(HPlist, int) else HPlist
        return {HP: [ self.get_coeff(number) for number in self._query_HP(JD, HP)] for HP in HPlist}


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





