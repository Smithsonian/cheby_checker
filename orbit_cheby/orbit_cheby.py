# -*- coding: utf-8 -*-
# sifter/sifter/query

'''
    --------------------------------------------------------------
    sifter's / mpchecker's orbit_cheby module.
    
    Jan 2020
    Matt Payne & Margaret Pan & Mike Alexandersen
    
    This module provides functionalities to evaluate
    dictionaries of chebyshev-coefficients
    
    We are developing a standardized approach regarding
    orbit integration and subsequent interpolation using
    chebyshev-coefficients applied to 32-day sectors
    
    To contain all core functionality required to predict the
    position (e.g. Heliocentric Cartesian) and apparent position
    (e.g. RA, Dec) of minor planets, comets, satellites, etc
    - See https://drive.google.com/open?id=1F86lnaHL01xPAACX2SYVdOgvvLLnxfYD
    
    This does *NOT* do the underlying nbody integrations.
    This DOES fast interpolation using supplied chebyshev dictionaries
    
    
    N.B. There is a strong assumption that coordinates are BARYCENTRIC EQUATORIAL
    - I.e. we are assuming that the nbody integration has supplied cartesian
    coordinates in a BARYCENTRIC EQUATORIAL frame.

    Currently envisaged as a core component that will be drawn-upon by
    - Position_Check
    - MP_Checker
    - ID_Check
    - ...
    
    --------------------------------------------------------------
    '''


# Import third-party packages
# --------------------------------------------------------------
import sys, os
import numpy as np
import operator
from collections import OrderedDict, defaultdict
import astropy
from astropy_healpix import HEALPix as hp
from astropy_healpix import healpy
from functools import lru_cache
import json
import itertools

# Import neighboring packages
# --------------------------------------------------------------
import nbody_reader

# Define top-line parameters
# --------------------------------------------------------------

secsPerDay = 86400

# We will assume that all sectors are of length 32-days
sector_length_days = 32

# We will assume that the earliest standard epoch to be used will be JD:2440000
# - This will set the enumeration of the sectors
standard_MJDmin = 40000 # 1968
standard_MJDmax = 64000 # 2034

# Healpix settings
HP_nside    = 16
HP_order    ='nested'
HPix        = hp(nside=HP_nside, order=HP_order)
HP_npix     = HPix.npix

class MSC():
    ''' 
            Multi-Sector Cheby Class
            
            Will hold chebyshev coeeficients for a given object
             
            *** TO DO *** 
            on init:
             - from filepath
             - from array ( & name/label)
             - from database
            on init (2) 
             - populate standard variables (global_standard_JDmin, supported_standard_JDmin, ... )
            on eval / look-up
             - map JD-to-sector
    '''

    # allow different init depending on source ...
    def __init__(self, **kwargs):
        
        # Initialization of default standard PARAMETERS / OPTIONS we use for chebyshevs, etc
        # - These may be overwritten by kwargs
        self.minorder       = 17                                    # : Fitting Chebys
        self.maxorder       = 25                                    # : Fitting Chebys
        self.maxerr         = 1e-8                                  # : Fitting Chebys
        self.FORCE_DATES    = False                                 # : Sector Ranges
        self.TDB_init       = None                                  # : Cheby Span
        self.TDB_final      = None                                  # : Cheby Span
        self.FROM_FILE      = False                                 # : Ingest method
        self.filepath       = False                                 # : Ingest method
        self.FROM_ARRAY     = False                                 # : Ingest method
        self.unpacked_provisional_designation            = None     # : Ingest method
        self.times_TDB      = None                                  # : Ingest method
        self.statearray     = None                                  # : Ingest method
        self.FROM_DATABASE  = None                                  # : Ingest method
        
        
        # Do a loop over the parameters and see if any attributes need updating (because they were supplied)
        for attribute, value in kwargs.items():
            if hasattr(self, attribute):
                setattr(self, attribute, value)
        
        # Initialization of the default fundamental MSC sectors
        # - These will be populated by the *_populate_from...()* functions below
        self.sectors                            = []
       
        # Allow initialization of ARRAYS / COEFFICIENTS from different sources
        if  self.FROM_FILE and self.filepath != None :                                                          ##<<-- Mainly for development
            self._populate_from_nbody_text(self.filepath)
        
        elif self.FROM_ARRAY and self.unpacked_provisional_designation != None and self.statearray != None :    ##<<-- From NBody
            self._populate_from_nbody_array(self.unpacked_provisional_designation, self.times_TDB, self.statearray)
        
        elif self.FROM_DATABASE:                                                                                ##<<-- From database
            self._populate_from_database(**kwargs)
        
        else:                                                                                                   ##<<-- An empty instance ...
            self._generate_empty()
                


    def _generate_empty(self,  ):
        ''' only making this a function on the off-chance that I need to dump other init. stuff in here '''
        print('Defaulting to the production of an empty MSCD')
    

    # Functions to handle different instantiaion methods
    #  - Take input and generate the chebyshev coefficients
    # --------------------------------------------------------------

    def _populate_from_nbody_text( self, text_filepath ):
        '''
            Read data (assumed from an nbody integration) from a text file 
            Use it to generate an array of positions, etc
            Then continue on to use func *_populate_from_nbody_array()*
            
            Mainly of use during development
            
            inputs:
            -------
            text_filepath: string
            - filepath to data from nbody simulation 
            
            returns:
            --------
            boolean
            - whatever output is from *_populate_from_nbody_array()*
            
        '''
        # Read the nbody-json file
        name, times, states  = nbody_reader.parse_nbody_txt( text_filepath )

        # Use the create-from-array function to do the rest
        return self._populate_from_nbody_array(name, times, states)



    def _populate_from_nbody_array(self,
                                   unpacked_provisional_designation ,
                                   times_TDB,
                                   states,
                                   ):
        '''
            Slightly/Highly rewritten version of Margaret's *datediv* routine
            
            inputs:
            -------
            name: string
             - unpacked_provisional_identification of the object
            times_TDB: np.array
             - TDB/TT
            states: np.array
             -
            
            returns:
            --------
            
            '''
        
        
        # I am going to adopt a fairly strict standardized approach:
        # *** Only sectors within standardized ranges will be generated ***
        # E.g. sector-zero ==>> standard_MJD0 -to- standard_MJD0 + sector_length_days, ...
        #
        # (0) Catch problem of incomplete info
        if   self.FORCE_DATES and (self.TDB_init == None or self.TDB_final == None):
            sys.exit('Cannot force when TDB_init and TDB_final not supplied')
        # (1) Allow the user to force specific start/end dates
        elif self.FORCE_DATES and self.TDB_init != None and self.TDB_final != None:
            print('Forcing, but will still use standard %r-day blocks...' % sector_length_days )
            if  self.TDB_init < times_TDB[0] or self.TDB_final > times_TDB[-1]:
                sys.exit(' nbody data does not support the requested dates ')
        # (2) Generate the start/end dates programatically (if not supplied/forced by user)
        elif not self.FORCE_DATES:
            # Here we compare the supplied JDs to the standard end-points
            self.TDB_init  = int(max(standard_MJDmin , times_TDB[0] ))
            self.TDB_final = int(min(standard_MJDmax , times_TDB[-1] ))

        else:
            sys.exit('Do not know how to calculate limits in *_generate_multi_sector_cheby_dict_from_nbody_array()* ')


        # Here we compare the endpoints to sector-endpoints to ensure sectors are fully supported by supplied data
        sector_numbers, sector_JD0s = self.map_JD_to_sector_number_and_sector_start_JD( [TDB_init,TDB_final] )
        sector_numbers[0]           = sector_numbers[0]  if TDB_init - sector_JD0s[0]                        < 1 else sector_numbers[0]  + 1
        sector_numbers[-1]          = sector_numbers[-1] if sector_JD0s[-1] + sector_length_days - TDB_final  < 1 else sector_numbers[-1] - 1


        # Sanity-checks
        assert TDB_init < TDB_final,   \
            ' Problem with TDB_init [%r] > TDB_init [%r] : likely due to supplied JDs falling outside standard ranges' % \
                (TDB_init , TDB_init)
        assert sector_numbers[0] <= sector_numbers[-1] , \
            ' Problem with sector numbers [%r, %r] : likely due to supplied JDs falling outside standard ranges' % \
                (sector_numbers[0] , sector_numbers[-1])


        # Go through sectors
        for ind in range(sector_numbers[0] , sector_numbers[-1] ):
            
            # Identify the indicees of the nbody times for this sector (i.e. those with min < t < max)
            # N.B. a[:,0] == times
            sector_TDB_init    = standard_MJDmin + ind     * sector_length_days
            sector_TDB_final   = standard_MJDmin + (ind+1) * sector_length_days
            indicees           = np.where((times_TDB >=sector_TDB_ini t)   & \
                                          (times_TDB <=sector_TDB_final)    )[0]

            # For each component array, generate chebys
            # Do all coordinates & covariances simultaneously
            maxorder       = min(maxorder,len(indicees))
            self.sectors.append(generate_single_sector_cheb_multi_coord( times_TDB[indicees],
                                                                            states[indicees],
                                                                            minorder,
                                                                            maxorder,
                                                                            maxerr))

        return True



    def _populate_from_database(self, args):
        '''
            Will need method to construct MSC from 1-or-many sectors stored in the sqlite db as coeffs
        '''
        pass


    # Assorted Utility Functions ...
    # --------------------------------------------------------------

    def map_JD_to_sector_number_and_sector_start_JD(self, JD_TDB, JD0=standard_MJDmin):
        '''
            For a given JD_TDB, calculate the sector number it will be in
            Assumes a standard starting julian date
            
            inputs:
            -------
            
            return:
            -------
            '''
        sector_number   = ( np.asarray( JD_TDB ) - JD0).astype(int) // sector_length_days
        sector_JD0      = JD0 + sector_number * sector_length_days
        return sector_number , sector_JD0





    # Functions to evaluate supplied multi-sector-cheby-dictionary
    # --------------------------------------------------------------

    def get_valid_range_of_dates( self,  ):
        '''
            Extract the minimum and maximum dates for
            which the supplied dictionary has valid
            chebyshev-coefficients
            
            Will work on single- or multi-sector cheby-dict
            
            inputs:
            -------
            cheby_dict: dictionary
            - see ... for details
            
            return:
            -------
            t_init : float
            - earliest valid date (time system may be MJD, but is implicit : depends on MPan;'s choice of zero-points ...)
            
            t_final : float
            - latest valid date (time system may be MJD, but is implicit : depends on MPan;'s choice of zero-points ...)
            '''
        return self.TDB_init , self.TDB_final


    def map_times_to_sectors( self, times_tdb ):
        '''
            Given query-times, it is likely to be useful to
            map each time to the relevant single-sector-dictionary
            
            inputs:
            -------
            times_tdb : np.array
            - JD TDB of times at which positions are to be calculated
            
            return:
            -------
            np.array of integers
            - sector # (zero-based) starting from the dictionary's "t_init"
            - length of returned array = len(times)
            '''
        return ( (times_tdb - self.TDB_init ) // sector_length_days ).astype(int)



    def generate_HP( self, times_tdb , observatoryXYZ , APPROX = False, CHECK = False ):
        '''
            Calculate apparent HP-locn from specified observatory-posn(s) at given time(s)
            N.B. observatory-posn(s) must be externally calculated/supplied
            
            inputs:
            -------
            times_tdb : np.array
            - JD TDB of times at which positions are to be calculated
            
            observatoryXYZ: np.array
            - Observatory positions at times.mjd [utc]
            - Dimension =3*len(times)
            
            APPROX: boolean
            - Allow approximate calc ( *no* LTT-correction) of unit-vector
            
            CHECK: boolean
            - Allow validity-checking to be turned on/off
            
            return:
            -------
            np.array of integer healpix
            - length of returned array = len(times)
            
            '''
        
        
        # Get the unit vector from observatory to object
        UV = self.generate_UnitVector(   times_tdb ,
                                    observatoryXYZ,
                                    APPROX = APPROX )
                                            
        # Calc the HP from the UV and return
        return healpy.vec2pix(HP_nside, UV[:,0], UV[:,1], UV[:,2], nest=True if HP_order=='nested' else False )


    def generate_UnitVector( self, times_tdb , observatoryXYZ, APPROX = False ):
        '''
            Calculate apparent UnitVector from specified observatory-posn(s) at given time(s)
            
            N.B. observatory-posn(s) must be :
            (i) externally calculated/supplied
            (ii) equatorial frame (to match assumed frame of orbit)
            
            inputs:
            -------
            times_tdb : np.array
             - JD TDB of times at which positions are to be calculated
            
            observatoryXYZ: np.array
            - Observatory positions at times_tdb
            - shape = (3,len(times_tdb))
            
            return:
            -------
            unit-vectors: np.array 
             - apparent UnitVector from specified observatory-posn(s) at given time(s)
            
            '''
        
        # Get the LTT-corrected position
        # - We allow for the possibility of *NOT* iterating (i.e. doing an approx calc.)
        n_iterations    = 1 if APPROX else 3
        lightDelay      = np.zeros(len(times_tdb))
        for i in range(n_iterations):
            
            # Calculate delayed time (of emission)
            delayedTimes    = times_tdb - lightDelay
            
            # Extract posn of objects at each delayed-time
            objectXYZ       = self.generate_XYZ( delayedTimes )
            
            # Calculate relative sepn-vector from observatory-to-object
            sepn_vectors    = objectXYZ - observatoryXYZ
            
            # Calculate distance to object at each time
            d               = np.linalg.norm(sepn_vectors, ord=2, axis=1, keepdims=True)
            
            # Calculate light-travel-time
            lightDelay      = d.flatten() / (astropy.constants.c * secsPerDay / astropy.constants.au ).value
        
        # Return unit-vector
        return sepn_vectors / d


    def generate_RaDec(self,  times_tdb  , observatoryXYZ, APPROX = False):
        '''
            Calculate apparent RA,DEC (Radians???) from specified observatory-posn(s) at given time(s)
            
            N.B. observatory-posn(s) must be :
            (i) externally calculated/supplied
            (ii) equatorial frame (to match assumed frame of orbit)
            
            inputs:
            -------
            times_tdb : np.array
            - JD TDB of times at which positions are to be calculated
            
            observatoryXYZ: np.array
            - Observatory positions at times_tdb
            - shape = (3,len(times_tdb))
            
            return:
            -------
            unit-vectors: np.array
            - apparent UnitVector from specified observatory-posn(s) at given time(s)
            
        '''
        # Get the unit vector from observatory to object
        UV = self.generate_UnitVector_from_cheby(times_tdb ,
                                            observatoryXYZ,
                                            APPROX = APPROX )
            
        # Convert from unit-vector to RA, Dec and then return
        theta_, RA_   = healpy.vec2ang(UV)
        return RA_ , 0.5*np.pi - theta_



    def generate_XYZ( self, times_tdb ):
        '''
            *** WE ASSUME THE cheby_dict IS VALID FOR ALL SUPPLIED times *************
            *** THIS IMPLIES PRE-CHECKING BY A HIGHER/PRECEEDING FUNCTION ************
            
            
            inputs:
            -------
            times_tdb : np.array
            - JD TDB of times at which positions are to be calculated
            
            return:
            -------
            XYZ posn
            - It is assumed that the chebys & hence the calculated XYZ are in barycentric, equatorial coords
              (this is to facilitate calc of RA,Dec)
            - shape = (3,len(times) )
            - in whatever frame cheby_dict
            '''
        
        # Find which single-sector dictionary to use for a given time
        sectors = self.map_times_to_sectors( times_tdb , multi_sector_cheby_dict )
        
        # Evaluate the chebyshev polynomial using the appropriate single-sector-dict
        #
        # Approach for dictionary / json ...
        #
        # - Seems likely that this current implementation is non-optimal ...
        #XYZs = [ np.array([
        #                   np.polynomial.chebyshev.chebval( t, cheby_dict["sectors"][s]["x"] ),
        #                   np.polynomial.chebyshev.chebval( t, cheby_dict["sectors"][s]["y"] ),
        #                   np.polynomial.chebyshev.chebval( t, cheby_dict["sectors"][s]["z"] )
        #                   ]) for t, s in zip(times.tdb, sectors) ]
        #
        #
        # Approach for arrays ...
        XYZs = [ np.polynomial.chebyshev.chebval( times.mjd[n] , multi_sector_cheby_dict["sectors"][s]['coeffs'][:,0:3] ) for n, s in enumerate( sectors ) ]
        return np.array(XYZs)





# End Of File
# --------------------------------------------------------------


"""
    def generate_multi_sector_cheby_dict_from_nbody_json( json_filepath ,
    TDB_init,
    TDB_final,
    minorder=17,
    maxorder=25,
    maxerr=1e-8,
    FORCE_DATES = False ,
    CHECK = False):
    '''
    Slightly rewritten version of Margaret's *datediv* routine
    
    inputs:
    -------
    json_filepath   -- filepath of json file to be read
    TDB_init,TDB_final -- earliest and latest MJD dates for which you want chebyshev approximations
    
    returns:
    --------
    multi-sector cheby-dict
    '''
    
    # Read the nbody-json file
    nbody_dict = nbody_reader.parse_nbody_json( json_filepath )
    
    # Check whether the supplied data can support the requested date-range
    if FORCE_DATES :
    TDB_init = max(TDB_init , nbody_dict[ nbody_reader.time_fieldname ][0] )
    TDB_final = min(TDB_final , nbody_dict[ nbody_reader.time_fieldname ][-1] )
    else:
    if  TDB_init < nbody_dict[ nbody_reader.time_fieldname ][0] or \
    TDB_final > nbody_dict[ nbody_reader.time_fieldname ][-1]:
    sys.exit(' nbody data does not support the requested dates ')
    
    # Set up a (mostly-empty) multi-sector cheby-dict
    mscd = {
    nbody_reader.object_name    : nbody_dict[nbody_reader.object_name],
    'MJP_TDB_init'      : TDB_init,
    'MJP_TDB_final'     : TDB_final,
    'sectors'           : []
    }
    
    # Split into sectors ...
    numdivs = int(np.ceil((TDB_final-TDB_init)/sector_length_days))
    
    # Go through sectors
    for ind in range(numdivs):
    
    # Set up a (mostly-empty) single-sector cheby-dict
    sscd = { nbody_reader.object_name : nbody_dict[nbody_reader.object_name] }
    
    # Identify the indicees of the nbody times for this sector (i.e. those with min < t < max)
    sscd['MJP_TDB_init']    = TDB_init + ind*sector_length_days
    sscd['MJP_TDB_final']   = min(TDB_final,TDB_init + (ind+1)*sector_length_days)
    indicees                = np.where((nbody_dict[ nbody_reader.time_fieldname ]>=sscd['MJP_TDB_init']) & \
    (nbody_dict[ nbody_reader.time_fieldname ]<=sscd['MJP_TDB_final']) )[0]
    
    
    # Loop over all coordinates & covariances
    lists = [nbody_reader.coord_names, nbody_reader.covar_names]
    for item in itertools.chain(*lists):
    
    # For each component, generate chebys & save into single-sector cheby-dict
    maxorder    = min(maxorder,len(indicees))
    sscd[item]  = generate_single_sector_cheb(  nbody_dict[ nbody_reader.time_fieldname ][indicees],
    nbody_dict[item][indicees],
    minorder,
    maxorder,
    maxerr)
    # append the single-sector cheby-dict into the multi-sector cheby-dict
    mscd['sectors'].append(sscd)
    
    # If being thorough, check that the produced object is valid
    if CHECK:
    check_multi_sector_validity( mscd )
    
    return mscd
    
    
    


    # Functions to check validity ...
    # --------------------------------------------------------------
    def check_validity( cheby_dict ):
        '''
            Check whether the input dictionary has the expected
            structure / variables of a MULTI-SECTOR dictionary

        '''
        # Expected keys & data-types
        expected_keys_and_types = [
        ("name", str),
        ("MJP_TDB_init", (int, float, np.int64, np.float64)),
        ("MJP_TDB_final", (int, float, np.int64, np.float64)),
        ("sectors", (list, np.ndarray) )
        ]

        # Check data is as expected
        # Needs to
        # (i) be a dict
        # (ii) have all the necessary keys
        # (iii) have all the correct data types
        # (iv) individual-sector dictionaries are all valid
        return True if isinstance(cheby_dict , dict) and \
                        np.all([key in cheby_dict and
                        isinstance(cheby_dict[key], typ) for key, typ in
                        expected_keys_and_types]) and \
                        np.all([check_single_sector_validity(sector_dict) for
                           sector_dict in cheby_dict[sectors]]) \
                    else False




    """
