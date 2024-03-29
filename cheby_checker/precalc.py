
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
from astropy.time import Time
import pickle


# Import neighboring packages
# --------------------------------------------------------------
from . import orbit_cheby
from . import sql
from . import nbody
from . import coco
from . import cheby_checker
from . import obs_pos


def end_to_end_precalc_wrapper(mpc_orb_jsonb):
    """

    :param mpc_orb_jsonb:
    :return:
    """
    precalc_object = PreCalc()

    # Calculate and insert the precalc results into the db.
    precalc_object.end_to_end_precalc([mpc_orb_jsonb])


# --------------------------------------------------------
# Primary external class for creating & accessing
# various pre-calculated quantities
# --------------------------------------------------------

class PreCalc(cheby_checker.Base, sql.SQLChecker, obs_pos.ObsPos):
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
        cheby_checker.Base.__init__(self)
        sql.SQLChecker.__init__(self)
        obs_pos.ObsPos.__init__(self)

        # connect to db
        self.db = sql.SQLChecker()


    # ------------------------------------------------------------------
    # High level pre-calculation function(s)
    # ------------------------------------------------------------------

    def end_to_end_precalc(self, mpc_orb_list, observatoryXYZ=None):
        """
        A high level function to handle ...
        (i) calling nbody on 1-or-many ORBFIT files
        (ii) calling MSCLoader on the results of (i)
        (iii) calling PreCalc.upsert() on the results of (ii)

        # NB, In general we will *NOT* be passing observatory coords

        :param mpc_orb_list: List of filenames or JSON blobs to precalc.
        """
        # Initiate NbodySim class
        N = nbody.NbodySim()

        # Run the integrator, by calling the object.
        N.run_mpcorb(tstart=self.standard_MJDmin, tstop=self.standard_MJDmax, mpcorb_list=mpc_orb_list)

        # Multi-Sector Cheby - Like an array of coeffs per sector.
        # Use the MSC_Loader to do all of the work to declare and populate a list of MSC objects
        # This is the step that creates the chebyshev coefficients!
        MSCs = orbit_cheby.MSC_Loader(NbodySim=N).MSCs

        # Do the precalculations and upsert
        # NB, In general we will *NOT* be passing observatory coords
        self.upsert(MSCs, helio_eq_observatoryXYZ=observatoryXYZ)


    # ------------------------------------------------------------------
    # Data insert/update function(s)
    # ------------------------------------------------------------------

    def upsert(self, MSCs, helio_eq_observatoryXYZ=None):
        """
            Main method used to insert/update coefficients for an object
            Also handles the healpix calculations used for efficiency-of-read

            Inputs:
            -------
            MSCs : List of Multi-Sector-Cheby class objects
             - MSC class is defined in orbit_cheby.py

            helio_eq_observatoryXYZ : np.array (optional)
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
        # This gets the position of the center of the earth at a particular time.
        t = Time(self.JDlist, format='jd', scale='tdb')
        t_tdb = t.tdb  # <<-- Convert to TDB
        if helio_eq_observatoryXYZ is None:
            t_utc = t.utc  # <<-- This converts to utc (from tdb, above)
            # TODO: This line can be improved by caching or another way. Likely using fixed JDs from a fixed list.
            helio_eq_observatoryXYZ = np.array([
                self.get_heliocentric_equatorial_xyz(jd, obsCode="500", verbose=False) for jd in t_utc.jd
            ])

            # Just a matrix rotation.
            bary_eq_observatoryXYZ  = coco.equatorial_helio2bary(helio_eq_observatoryXYZ, t_tdb.jd)
        else:
            bary_eq_observatoryXYZ = coco.equatorial_helio2bary(helio_eq_observatoryXYZ, t_tdb.jd)

        # iterate over each MSC in list ...
        for M in MSC_list:
            # update list of coefficients for individual MSC
            # NB This is required for "phase-1" (Ephemeris Service)
            object_coeff_id = self.upsert_MSC_coefficients(M)

            # update healpix locations for individual MSC for valid integer-JD in JDlist
            # NB This is required for "phase-2" (Pointer Service)
            self.upsert_MSC_HP(M, bary_eq_observatoryXYZ, object_coeff_id)


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
        sector_field_names = self.generate_sector_field_names(
            sector_dict=
            {sector_num: self.map_sector_number_to_sector_start_JD(sector_num, self.standard_MJDmin) for sector_num in M.sector_coeffs.keys()}
        )

        # (ii) Get the coefficients for each sector
        sector_field_values = [np.array2string(item, separator=",").replace('\n', '') for item in M.sector_coeffs.values()]

        # (iii) Do the insert
        object_coeff_id = self.upsert_coefficients(M.unpacked_primary_provisional_designation, sector_field_names, sector_field_values)

        return object_coeff_id


    def upsert_MSC_HP(self, M, observatoryXYZ, object_coeff_id):
        """
            Upload nightly HealPix (HP) locations for a single MSC to the database

            *** TESTING NEEDS TO BE CREATED FOR THIS FUNCTION: MJP 2021-06-10 ***
            *** Consider whether to relocate to MSC object in orbit_cheby    ***

            inputs:
            -------

            M : MSC-object
             - see orbit_cheby module for detailed specification
             - here we need M to possess a dictionary-attribute named "sector_coeffs"

            helio_eq_observatoryXYZ : np.array (optional)
             - Position of the observatory from which you want to calculate the
               nightly healpix (in Heliocentric Equatorial coords)
             - Everything in orbit_cheby.MSC needs BARYCENTRIC EQUATORIAL coords
            - shape = ( len(times_tdb), 3)

            return:
            -------
            None

        """
        # Use the coefficient-dictionary(ies) to get the HP for each integer-JD in JDlist
        indicees = np.where( self.JDlist < M.get_valid_range_of_dates()[1] )[0]
        HPlist   = M.generate_HP(self.JDlist[indicees],  observatoryXYZ[indicees, :], APPROX=True)

        # update HP data in db
        self.db.insert_HP(self.JDlist[indicees], HPlist, object_coeff_id)


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
        # - Outer dict is key-ed on unpacked_primary_provisional_designation
        # - Inner dicts are key-ed on sector
        dict_of_dicts = query_coefficients_by_jd_hp(conn, JD, HPlist , sector_numbers = None)

        # Swap the object_ids for the designations
        # *** THIS IS THE DESIRED SIGNATURE (2021-06-10) :   ***
        # *** MSC_Loader NEEDS TO BE REWRITTEN TO ALLOW THIS ***
        return orbit_cheby.MSC_Loader(  FROM_DATABASE = True ,
                                        dict_of_dicts = dict_of_dicts ).MSCs
