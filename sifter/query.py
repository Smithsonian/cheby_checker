# -*- coding: utf-8 -*-
# sifter/sifter/query

'''
    --------------------------------------------------------------
    sifter's query module.
    
    Jan 2020
    Matt Payne & Mike Alexandersen
    
    This module provides functionalities to
    (a) ...
    
    *WRITE MORE STUFF*
    
    --------------------------------------------------------------
'''


# Import third-party packages
# --------------------------------------------------------------
import sys, os
import numpy as np
import operator
from collections import OrderedDict, defaultdict
from astropy_healpix import HEALPix as hp
from astropy.time import Time
from functools import lru_cache
import json

import random # only necessary while developing


# Import neighboring packages
# --------------------------------------------------------------
from . import sql
import orbit_cheby as cheby




    
def query( cheby_dict_for_orbit , param_dict , EXPLICIT_CHECK = False ):
    '''
        Overall function to query ITF (& orbits) against supplied orbit(s)

        inputs:
        -------
        cheby_dict_for_orbit : dictionary
        - multi-sector (dict-of-dicts) to represent orbit using checbshev coeffs
        - see ... for details
        
        param_dict: dictionary
        - specify the search parameters/tolerances with which to match (tracklets to input orbit)
        
        return:
        -------
        ??? : tracklet_names ? tracklet_dictionaries ? other ?
    '''
    
    # Check whether the inputs are of the allowed format
    if EXPLICIT_CHECK:
        assert _check_query_inputs(  cheby_dict_for_orbit , param_dict )
    
    # Get Nightly Healpix (evaluates orbital positions on integer days)
    JD_list, HP_list = _get_nightly_healpix( cheby_dict_for_orbit )
    
    # Query pre-calcs for approximate matches ( sql look-up)
    list_of_tracklet_tuples = self._query_precalc(JD_list, HP_list)
    
    # Later could have a refinement step using RoM & AoM
    # _refine_using_RoMAoM()
    
    # Refine approx matches to get "PRECISE" matches (within specified tolerance)
    list_of_tracklet_tuples_precise = self._get_precise_matches( list_of_tracklet_tuples , cheby_dict_for_orbit , param_dict )

    # Return data in desired format (can include screen-dump)
    # - Command-line-options / input-variables may dictate format
    return get_results( list_of_tracklet_tuples_precise , param_dict )


def _check_query_inputs( cheby_dict_for_orbit , param_dict ):
    '''
        Check whether the inputs are of the allowed format
    '''
    # Are they both dictionaries ?
    assert isinstance(cheby_dict_for_orbit , dict) and isinstance(param_dict , dict), \
        ' ... not dictionaries ... types = [%r, %r]' % ( type(cheby_dict_for_orbit), type(param_dict) )
    
    # Does the cheby_dict have the correct content/format ?
    # cheby.check_validity( cheby_dict_for_orbit )
    
    # Does the param_dict have the correct keys ?
    
    return True


def _get_nightly_healpix( cheby_dict_for_orbit ):
    '''
        Generate a list of JD,HP pairs that will be searched 
         - I.e. which healpix to search on which nights 
         
        Do this by evaluating the orbital position at integer JDs
         
        N.B.
        (i)
        We could search every night POPULATED in the ITF
        However, the cheby_dict is only valid over a limited time range 
        So let's just default to searching that valid JD-range 
        
        (ii)
        Not sure what to do about orbital uncertainties / variant orbits
        
        Probably return a list of healpix for each JD
    '''
    
    
    # Establish range of valid integer days from cheby-dictionary
    print( '***DUMMY VALUES: NEED TO POPULATE CORRECTLY !!! ***' )
    JDmin, JDmax = cheby.get_valid_range_of_dates_from_cheby( cheby_dict_for_orbit )
    
    # Use orbit_cheby functionality, *generate_HP_from_cheby()*, to get the HP at nightly intervals
    times = Time( 42000. + np.arange( nSample ) , format='mjd', scale='tdb')
    HP_list = cheby.generate_HP_from_cheby( times , cheby_dict_for_orbit )
    
    return times, HP_list

def _query_precalc( JD_list, HP_list):
    '''
        Query pre-calculated data for *approximate* matches
         - Looks for tracklets in the same HP on the same night
         - Pre-calcs are accessed via sql query
         
        inputs:
        -------
        JD_list: list of integers
         - dates to search
        HP_list: list of integers
         - healpix to search
        
        return:
        -------
        list_of_tracklet_tuples: (tracklet_name, tracklet_dictionary)
         - same structure as the return from sql.query_tracklets
         
    '''
    list_of_tracklet_tuples = []
    for MJD, HP in zip( times.mjd, HP_list):
        list_of_tracklet_tuples.extend( sql.query_tracklets_jdhp(MJD, HP) )
    return list_of_tracklet_tuples

#def._refine_using_RoMAoM()
#    pass

def _get_precise_matches( list_of_tracklet_tuples_approx, cheby_dict_for_orbit , param_dict):
    '''
        Refine the approx matches (from _query_precalc) 
        down to exact matches within the specified tollerances
        
        inputs:
        -------
        list_of_tracklet_tuples_approx: list of tuples
        - each tuple is (tracklet_name, tracklet_dictionary)
        - same structure as list_of_tracklet_tuples from _query_precalc
        
        cheby_dict_for_orbit : dictionary
        - multi-sector (dict-of-dicts) to represent orbit using checbshev coeffs
        - see ... for details
        
        param_dict: dictionary
        - specify the search parameters/tolerances with which to match (tracklets to input orbit)
        
        
        return:
        -------
        list_of_tracklet_tuples_precise: (tracklet_name, tracklet_dictionary)
        - each tuple is (tracklet_name, tracklet_dictionary)
        - same structure as list_of_tracklet_tuples from _query_precalc
        - these are the tracklets which are "close" to the input orbit
        
    '''
    
    """
    # empty lists ...
    list_of_tracklet_tuples_precise = []
    startTimes , endTimes, start_observatoryXYZs, end_observatoryXYZs = [],[],[]
    
    for tracklet_tup in list_of_tracklet_tuples_approx:
        # split tuple into key () and value (tracklet_dict made in *precalc.parse_observation_list()* )
        key, tracklet_dict = tracklet_tup
    
        # Parse obs80 lines
        # *** WE SHOULD PUSH AS MUCH AS POSSIBLE TO PREPROCESSING ***
        parsed = [obs for obs in parse80(tracklet_dict['observations'])]
    
        # Get times
        startTimes.append( parsed.jdutc )
        endTimes.append( parsed.jdutc )

        # - Extract observatory positions at times
        start_observatoryXYZs.append( tracklet_dict['obsXYZ'][0] )
        end_observatoryXYZs.append( tracklet_dict['obsXYZ'][-1] )

    startTimes  = Time( startTimes , format='mjd', scale='tdb')
    endTimes    = Time( endTimes , format='mjd', scale='tdb')
    start_observatoryXYZs = np.array(start_observatoryXYZs)
    end_observatoryXYZs   = np.array(end_observatoryXYZs)

    # - Get precise unit vector to TRACKLET OBSERVATIONS at each time
    # *** WE SHOULD PUSH AS MUCH AS POSSIBLE TO PREPROCESSING ***
    start_trackletUV = ...
    end_trackletUV = ...

    # - Get precise unit vector to ORBIT at each time
    start_orbitUV = generate_UnitVector_from_cheby( startTimes ,
                                                    multi_sector_cheby_dict ,
                                                    observatoryXYZ,
                                                    APPROX = APPROX )

    # - Get the dot-product, and hence angular-sepn, at each time
    start_angles = np.dot(start_orbitUV , start_trackletUV)
    
    # - Select the subset of tracklets which match the criteria defined in the param_dict
    
    """
    
    print( '***DUMMY VALUES: NEED TO POPULATE CORRECTLY !!! ***' )
    if random.randrange(2):
        list_of_tracklet_tuples_precise.append( ( tracklet_name, tracklet_dict ) )
    
    return precise_dictionary


def get_results( list_of_tracklet_tuples_precise , param_dict ):
    '''
        Return nicely formatted results (to be defined)
        - Command-line-options / input-variables may dictate format
        
        inputs:
        -------
        list_of_tracklet_tuples_precise: list of tuples
        - each tuple is (tracklet_name, tracklet_dictionary)
        - this is what is returned from from _get_precise_matches
        
        param_dict: dictionary
        - specify the search parameters/tolerances with which to match (tracklets to input orbit)
        
        
        return:
        -------
        Command-line-options / input-variables may dictate format
        
    '''
    
    # *** DEBUGGING : JUST PRINT EVERYTHING TO SCREEN FOR NOW ***
    for tracklet_name, tracklet_dict in list_of_tracklet_tuples_precise:
        print( tracklet_name )
    
    # *** DEBUGGING : JUST RETURN THE PRECISE LIST FOR NOW ***
    return list_of_tracklet_tuples_precise

