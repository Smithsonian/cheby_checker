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

class MSCD(dict):
    ''' 
            Multi-Sector Cheby Dict
            
             - subclass of standard dict
             - extra functionality implemented to hold other standard info
             
             
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
    def __init__(self, *args):
        
        # Do the complicated initialization
        if   FROM_FILE in args and filepath in args and FROM_FILE:
            self._generate_multi_sector_cheby_dict_from_nbody_text(args)
        elif FROM_ARRAY in args and FROM_ARRAY and name in args and nparray in args :
            self._generate_multi_sector_cheby_dict_from_array(args)
        elif FROM_DATABASE in args and FROM_DATABASE :
            self._generate_multi_sector_cheby_dict_from_database(args)
        else:
            self._generate_empty():
                
        


    def _generate_empty():
        print('Defaulting to the production of an empty MSCD')

# Functions to create / initialize  ...
# --------------------------------------------------------------




def error(a, b):
    ''' get cheb polynomial approximation error at fitted points'''
    
    err = 0
    for i in range(len(a)):
        err = max(err, abs(a[i]-b[i]))
    return err



def _generate_multi_sector_cheby_dict_from_nbody_text( text_filepath ,
                                                     minorder=17,
                                                     maxorder=25,
                                                     maxerr=1e-8,
                                                     FORCE_DATES = False ,
                                                     mindate = None,
                                                     maxdate = None,
                                                     CHECK = False):
    '''
      ...
    '''
    # Read the nbody-json file
    name, a = nbody_reader.parse_nbody_txt( text_filepath )

    # Use the create-from-array function to do the rest
    return generate_multi_sector_cheby_dict_from_nbody_array(name, a)



def _generate_multi_sector_cheby_dict_from_nbody_array( name ,
                                                      array,
                                                     minorder=17,
                                                     maxorder=25,
                                                     maxerr=1e-8,
                                                     FORCE_DATES = False ,
                                                     mindate = None,
                                                     maxdate = None,
                                                     CHECK = False):
    '''
        Slightly rewritten version of Margaret's *datediv* routine
        
        inputs:
        -------
        json_filepath   -- filepath of json file to be read
        mindate,maxdate -- earliest and latest MJD dates for which you want chebyshev approximations
        
        returns:
        --------
        multi-sector cheby-dict
        '''
    
    
    # I am going to adopt a fairly strict standardized approach:
    # *** Only sectors within standardized ranges will be generated ***
    # E.g. sector-zero ==>> standard_MJD0 -to- standard_MJD0 + sector_length_days, ...
    #
    #
    if   FORCE_DATES and (mindate == None or maxdate == None):
        sys.exit('Cannot force when mindate and maxdate not supplied')
    
    elif FORCE_DATES and mindate != None and maxdate != None:
        print('Forcing, but will still use standard %r-day blocks...' % sector_length_days )
        if  mindate < a[0 ,0] or maxdate > a[-1,0]:
            sys.exit(' nbody data does not support the requested dates ')

    elif not FORCE_DATES:
        # Here we compare the supplied JDs to the standard end-points
        mindate = int(max(standard_MJDmin , a[0 ,0] ))
        
        maxdate = int(min(standard_MJDmax , a[-1,0] ))

    else:
        sys.exit('Do not know how to calculate limits in *generate_multi_sector_cheby_dict_from_nbody_text()* ')


    # Here we compare the endpoints to sector-endpoints to ensure sectors are fully supported by supplied data
    sector_numbers, sector_JD0s = map_JD_to_sector_number_and_sector_start_JD( [mindate,maxdate] )
    sector_numbers[0]           = sector_numbers[0]  if mindate - sector_JD0s[0]                        < 1 else sector_numbers[0]  + 1
    sector_numbers[-1]          = sector_numbers[-1] if sector_JD0s[-1] + sector_length_days - maxdate  < 1 else sector_numbers[-1] - 1

    # Sanity-checks
    assert mindate < maxdate,   ' Problem with mindate [%r] > mindate [%r] : likely due to supplied JDs falling outside standard ranges' % (mindate , mindate)
    assert sector_numbers[0] <= sector_numbers[-1] , ' Problem with sector numbers [%r, %r] : likely due to supplied JDs falling outside standard ranges' % (sector_numbers[0] , sector_numbers[-1])

    # Set up a (mostly-empty) multi-sector cheby-dict
    mscd = {
        nbody_reader.object_name : name,
            'MJP_TDB_init'      : mindate,
            'MJP_TDB_final'     : maxdate,
            'sectors'           : []
        }
    print('mscd', mscd)

    # Go through sectors
    for ind in range(sector_numbers[0] , sector_numbers[-1] ):
        
        # Set up a (mostly-empty) single-sector cheby-dict
        sscd = { nbody_reader.object_name : name }
        
        # Identify the indicees of the nbody times for this sector (i.e. those with min < t < max)
        # N.B. a[:,0] == times
        sscd['MJP_TDB_init']    = standard_MJDmin + ind     * sector_length_days
        sscd['MJP_TDB_final']   = standard_MJDmin + (ind+1) * sector_length_days
        indicees                = np.where((a[:,0]>=sscd['MJP_TDB_init']) & \
                                           (a[:,0]<=sscd['MJP_TDB_final']) )[0]
            
            

        # For each component array, generate chebys
        # Do all coordinates & covariances simultaneously
        # N.B. time            == a[indicees,0]
        #      coords at times == a[indicees,1:]
        maxorder       = min(maxorder,len(indicees))
        sscd['coeffs'] = generate_single_sector_cheb_multi_coord(   a[indicees,0],
                                                                    a[indicees,1:],
                                                                    minorder,
                                                                    maxorder,
                                                                    maxerr)
        # save into multi-sector cheby-dict
        mscd['sectors'].append(sscd)

    return mscd


    """
    def generate_multi_sector_cheby_dict_from_nbody_json( json_filepath ,
                                                         mindate,
                                                         maxdate,
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
            mindate,maxdate -- earliest and latest MJD dates for which you want chebyshev approximations
            
            returns:
            --------
            multi-sector cheby-dict
            '''
        
        # Read the nbody-json file
        nbody_dict = nbody_reader.parse_nbody_json( json_filepath )
        
        # Check whether the supplied data can support the requested date-range
        if FORCE_DATES :
            mindate = max(mindate , nbody_dict[ nbody_reader.time_fieldname ][0] )
            maxdate = min(maxdate , nbody_dict[ nbody_reader.time_fieldname ][-1] )
        else:
            if  mindate < nbody_dict[ nbody_reader.time_fieldname ][0] or \
                maxdate > nbody_dict[ nbody_reader.time_fieldname ][-1]:
                sys.exit(' nbody data does not support the requested dates ')

        # Set up a (mostly-empty) multi-sector cheby-dict
        mscd = {
            nbody_reader.object_name    : nbody_dict[nbody_reader.object_name],
                    'MJP_TDB_init'      : mindate,
                    'MJP_TDB_final'     : maxdate,
                    'sectors'           : []
            }

        # Split into sectors ...
        numdivs = int(np.ceil((maxdate-mindate)/sector_length_days))
            
        # Go through sectors
        for ind in range(numdivs):
            
            # Set up a (mostly-empty) single-sector cheby-dict
            sscd = { nbody_reader.object_name : nbody_dict[nbody_reader.object_name] }
            
            # Identify the indicees of the nbody times for this sector (i.e. those with min < t < max)
            sscd['MJP_TDB_init']    = mindate + ind*sector_length_days
            sscd['MJP_TDB_final']   = min(maxdate,mindate + (ind+1)*sector_length_days)
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
    """

    # Assorted Utility Functions ...
    # --------------------------------------------------------------


    def map_JD_to_sector_number_and_sector_start_JD(JD, JD0=standard_MJDmin):
        '''
            For a given JD, calculate the sector number it will be in
            Assumes a standard starting julian date
            
            inputs:
            -------
            
            return:
            -------
            '''
        sector_number   = ( np.asarray( JD ) - JD0).astype(int) // sector_length_days
        sector_JD0      = JD0 + sector_number * sector_length_days
        return sector_number , sector_JD0




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






# Functions to evaluate supplied cheby-dictionaries
# --------------------------------------------------------------

def get_valid_range_of_dates_from_cheby( cheby_dict ):
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
    return cheby_dict["MJP_TDB_init"], cheby_dict["MJP_TDB_final"]


def map_times_to_sectors( times , cheby_dict ):
    '''
        Given query-times, it is likely to be useful to
        map each time to the relevant single-sector-dictionary
        
        inputs:
        -------
        times : astropy.Time object
        - Using astropy.Time object to stop making mistakes w.r.t utc/tt/tdb/...
        - Can be singular or multiple
        
        cheby_dict: dictionary
        - see ... for details
        
        return:
        -------
        np.array of integers
        - sector # (zero-based) starting from the dictionary's "t_init"
        - length of returned array = len(times)
        '''
    return ( (times.mjd - cheby_dict["MJP_TDB_init"] ) // sector_length_days ).astype(int)



def generate_HP_from_cheby( times , multi_sector_cheby_dict , observatoryXYZ , APPROX = False, CHECK = False ):
    '''
        Calculate apparent HP-locn from specified observatory-posn(s) at given time(s)
        N.B. observatory-posn(s) must be externally calculated/supplied
        
        inputs:
        -------
        times : astropy.Time object
        - Using these to stop making mistakes w.r.t utc/tt/tdb/...
        - Can be singular or multiple
        
        cheby_dict: dictionary
        - multi-sector dictionary to represent orbit using checbshev coeffs
        - see ... for details
        
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
    UV = generate_UnitVector_from_cheby(times ,
                                        multi_sector_cheby_dict ,
                                        observatoryXYZ,
                                        APPROX = APPROX )
                                        
                                        # Calc the HP from the UV and return
    return healpy.vec2pix(HP_nside, UV[:,0], UV[:,1], UV[:,2], nest=True if HP_order=='nested' else False )

def generate_UnitVector_from_cheby( times , multi_sector_cheby_dict , observatoryXYZ, APPROX = False ):
    '''
        Calculate apparent UnitVector from specified observatory-posn(s) at given time(s)
        N.B. observatory-posn(s) must be externally calculated/supplied
        
        inputs:
        -------
        times : astropy.Time object
        - Using these to stop making mistakes w.r.t utc/tt/tdb/...
        - Can be singular or multiple
        
        cheby_dict: dictionary
        - multi-sector dictionary to represent orbit using checbshev coeffs
        - see ... for details
        
        observatoryXYZ: np.array
        - Observatory positions at times.jd [utc]
        - Dimension =3*len(times)
        
        return:
        -------
        
        
        '''
    
    # Get the LTT-corrected position
    # - We allow for the possibility of *NOT* iterating (i.e. doing an approx calc.)
    n_iterations    = 1 if APPROX else 3
    lightDelay      = np.zeros(len(times))
    for i in range(n_iterations):
        
        # Calculate delayed time (of emission)
        delayedTimes    = times - lightDelay
        
        # Extract posn of objects at each delayed-time
        objectXYZ       = generate_XYZ_from_cheby( delayedTimes, multi_sector_cheby_dict )
        
        # Calculate relative sepn-vector from observatory-to-object
        sepn_vectors    = objectXYZ - observatoryXYZ
        
        # Calculate distance to object at each time
        d               = np.linalg.norm(sepn_vectors, ord=2, axis=1, keepdims=True)
        
        # Calculate light-travel-time
        lightDelay      = d.flatten() / (astropy.constants.c * 86400 / astropy.constants.au ).value
    
    # Return unit-vector
    return sepn_vectors / d


def generate_RaDec_from_cheby( JD , multi_sector_cheby_dict , observatoryXYZ, APPROX = False):
    '''
        
        '''
    # Get the unit vector from observatory to object
    UV = generate_UnitVector_from_cheby(times ,
                                        multi_sector_cheby_dict ,
                                        observatoryXYZ,
                                        APPROX = APPROX )
        
    # Convert from unit-vector to RA, Dec and then return
    theta_, RA_   = healpy.vec2ang(np.dot( UV , np.transpose(PHYS.rot_mat))) # in equatorial
    return RA_ , 0.5*np.pi - theta_



def generate_XYZ_from_cheby( times , multi_sector_cheby_dict ):
    '''
        *** WE ASSUME THE cheby_dict IS VALID FOR ALL SUPPLIED times *************
        *** THIS IMPLIES PRE-CHECKING BY A HIGHER/PRECEEDING FUNCTION ************
        
        
        inputs:
        -------
        times : astropy.Time object
        - Using these to stop making mistakes w.r.t utc/tt/tdb/...
        - Can be singular or multiple
        
        cheby_dict: dictionary
        - multi-sector (dict-of-dicts) to represent orbit using checbshev coeffs
        - see ... for details
        
        return:
        -------
        XYZ posn
        - shape = (3,len(times) )
        - in whatever frame cheby_dict
        '''
    
    # Find which single-sector dictionary to use for a given time
    sectors = map_times_to_sectors( times , multi_sector_cheby_dict )
    
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






