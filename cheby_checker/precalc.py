# -*- coding: utf-8 -*-
# mpchecker2/mpchecker2/precalc.py

"""
    --------------------------------------------------------------
    Precalculation module for checker

    2020-2022
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
    """


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
from . import orbit_cheby
from . import sql
from . import obs_pos
from . import nbody


# --------------------------------------------------------
# Primary external class for creating & accessing
# various pre-calculated quantities
# --------------------------------------------------------

class PreCalc(orbit_cheby.Base , obs_pos.ObsPos, sql.SQLChecker):
    """
        Primary External Class for accessing ChebyChecker's
        pre-calculated data

        Methods
        -------

        upsert()
        query_HPlist()

    """
    def __init__(self):
        
        # Give access to "Base" & "ObsPos" methods & attributes
        orbit_cheby.Base.__init__(self)
        obs_pos.ObsPos.__init__(self)
        sql.SQLChecker.__init__(self)

        # connect to db
        self.db = sql.SQLChecker()
        

    # ------------------------------------------------------------------
    # High level pre-calculation function(s)
    # ------------------------------------------------------------------

    def end_to_end_precalc(self,filenames, observatoryXYZ=None):
        """
        A high level function to handle ...
        (i) calling nbody on 1-or-many ORBFIT files
        (ii) calling MSCLoader on the results of (i)
        (iii) calling PreCalc.upsert() on the results of (ii)

        # NB, In general we will *NOT* be passing observatory coords

        """
        # Initiate NbodySim class
        N = nbody.NbodySim()
                                    
        # Run the integrator, by calling the object.
        N.run_mpcorb(tstart=self.standard_MJDmin , tstop=self.standard_MJDmax, mpcorb_list = filenames )
        
        # Use the MSC_Loader to do all of the work to declare and populate a list of MSC objects
        MSCs = orbit_cheby.MSC_Loader(NbodySim = N).MSCs
                                        
        # Do the precalculations and upsert
        # NB, In general we will *NOT* be passing observatory coords
        self.upsert( MSCs , observatoryXYZ=observatoryXYZ )


    # ------------------------------------------------------------------
    # Data insert/update function(s)
    # ------------------------------------------------------------------

    def upsert(self, MSCs , observatoryXYZ=None ):
        """
            Main method used to insert/update coefficients for an object
            Also handles the healpix calculations used for efficiency-of-read

            Inputs:
            -------
            MSCs : List of Multi-Sector-Cheby class objects
             - MSC class is defined in orbit_cheby.py

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
        """
        # ensure that the supplied variable is formatted correctly
        MSC_list = self._rectify_inputs(MSCs)
        
        # If observatory coords are not supplied, set to be the location of the geocenter
        # N.B. (1) Default list of Julian Dates to use, self.JDlist = (2440000 ==> 1968, 2464000.0 ==> 2034)
        # N.B. (2) This obs-posn code likely takes a long time to execute in its draft/old version
        if observatoryXYZ is not None:
            observatoryXYZ = self.get_heliocentric_equatorial_xyz(self.JDlist,
                                                                    obsCode="500",
                                                                    verbose=False)
                                        
        # Ensure that the observatory coords have the right shape
        assert observatoryXYZ.shape == (3 , len(self.JDlist) ), \
            f'the shape of observatoryXYZ needs to be {(3, len(self.JDlist))} : instead it is {observatoryXYZ.shape}'
        
        # iterate over each MSC in list ...
        for M in MSC_list:
            
            # update list of coefficients for individual MSC
            # NB This is required for "phase-1" (Ephemeris Service)
            self.upsert_MSC_coefficients(M)
            
            # update healpix locations for individual MSC for valid integer-JD in JDlist
            # NB This is required for "phase-2" (Pointer Service)
            self.upsert_MSC_HP(M, observatoryXYZ)
            

    def _rectify_inputs(self,  MSCs ):
        """
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

        """
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

    def upsert_MSC_coefficients(self, M):
        """
            insert/update multi_sector_cheby object
            
            *** TESTING NEEDS TO BE REDONE FOR THIS FUNCTION: MJP 2021-06-10 ***
            *** Consider whether to relocate to MSC object in orbit_cheby    ***
            
            N.B. (1) ...
            https://stackoverflow.com/questions/198692/can-i-pickle-a-python-dictionary-into-a-sqlite3-text-field
            pdata = cPickle.dumps(data, cPickle.HIGHEST_PROTOCOL)
            curr.execute("insert into table (data) values (:data)", sqlite3.Binary(pdata))

            N.B. (2) ...
            The insert statement wouldn't have to look so terrible if the table was constructed
            differently/more-simply in create_object_coefficients_table

            inputs:
            -------
            
            M : MSC-object
             - see orbit_cheby module for detailed specification
             - here we need M to possess a dictionary-attribute named "sector_coeffs"
                         
            return:
            -------
            None


        """
                
        # (i) Get the sector field names required for this specific MSC
        # - Current method seems unnecessarily verbose
        sector_field_names =   self.generate_sector_field_names( sector_dict = \
                                                         { sector_num: Base().map_sector_number_to_sector_start_JD(sector_num , Base().standard_MJDmin) for sector_num in M.sector_coeffs.keys()}
                                                         )
        
        # (ii) Get the coefficients for each sector
        sector_field_values = [pickle.dumps(coeffs, pickle.HIGHEST_PROTOCOL) for coeffs in M.sector_coeffs.values()]
        
        # (iii) Do the insert
        self.upsert_coefficients(M.primary_unpacked_provisional_designation , sector_field_names, sector_field_values)
     
    def upsert_MSC_HP(self, M, observatoryXYZ):
        """
            Upload nightly HP locations for a single MSC to the database

            *** TESTING NEEDS TO BE CREATED FOR THIS FUNCTION: MJP 2021-06-10 ***
            *** Consider whether to relocate to MSC object in orbit_cheby    ***

            inputs:
            -------

            M : MSC-object
             - see orbit_cheby module for detailed specification
             - here we need M to possess a dictionary-attribute named "sector_coeffs"

            observatoryXYZ : np.array (optional)
             - Position of the observatory from which you want to calculate the
               nightly healpix (in Heliocentric Equatorial coords)

            return:
            -------
            None

        """
        # Use the coefficient-dictionary(ies) to get the HP for each integer-JD in JDlist
        # NB: need to restrict the queried dates to the those supported by the MSC
        indicees = np.where( Base().JDlist < M.get_valid_range_of_dates()[1] )[0]
        HPlist   = M.generate_HP(Base().JDlist[indicees],  observatoryXYZ[:,indicees] , APPROX = True)

        # update HP data
        self.db.insert_HP(self.db.conn, self.JDlist[indicees], HPlist, object_id)


    # ------------------------------------------------------------------
    # Data query function(s)
    # ------------------------------------------------------------------

    def get_nightly_precalcs(self,JD, HPlist):
        """
            Main convenience method to-be-used to query the precalculated healpix
            For a given JD & HP in a list of HPs, returns a list of MSC objects
            that are in those HP on that night

            Note that at present *no* attempt is made to ***ALWAYS RETURN NEOs***
            - This needs to be added-to / improved
            - We could (perhaps) return everything closer than 0.X au (on a given JD)?

            inputs
            ------
            JD: float or int
             - julian date of the night. If <float> will be silently converted to <int>

            HPlist: list-of-integers
             - healpix to be queried
             - if integer (single healpix) supplied, is silently converted to list

            returns
            -------
            # *** THIS IS THE DESIRED SIGNATURE (2021-06-10) :   ***
            # *** MSC_Loader NEEDS TO BE REWRITTEN TO ALLOW THIS ***
            MSCs : List of Multi-Sector-Cheby class objects
        """

        # Get the coefficients for the objects
        # - Outer dict is key-ed on primary_unpacked_provisional_designation
        # - Inner dicts are key-ed on sector
        dict_of_dicts = query_coefficients_by_jd_hp(conn, JD, HPlist , sector_numbers = None)
        
        # Swap the object_ids for the designations
        # *** THIS IS THE DESIRED SIGNATURE (2021-06-10) :   ***
        # *** MSC_Loader NEEDS TO BE REWRITTEN TO ALLOW THIS ***
        return orbit_cheby.MSC_Loader(  FROM_DATABASE = True ,
                                        dict_of_dicts = dict_of_dicts ).MSCs
    
    





