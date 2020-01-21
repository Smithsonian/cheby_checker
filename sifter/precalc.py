# -*- coding: utf-8 -*-
# sifter/sifter/precalc.py

'''
    --------------------------------------------------------------
    sifter's precalculation module.
    
    Jan 2020
    Matt Payne & Mike Alexandersen
    
    This module provides functionalities to
    (a) create & save new precalculations
    (b) load extant precalculations
    
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
from functools import lru_cache
import json


# Import neighboring packages
# --------------------------------------------------------------



# Default for caching stuff using lru_cache
# -------------------------------------------------------------
cache_size_default = 16

# Data classes/methods
# -------------------------------------------------------------

class Base():
    '''
        Parent class to hold some file/directory definitions & methods
    '''

    def __init__(self):

        # Filing ...
        self.HP     = '_healpix.txt'

        # Healpix ...
        self.HP_nside=16
        self.HP_order='nested'
        self.npix = hp(nside=self.HP_nside, order=self.HP_order)

        # Default list of Julian Dates to use (2440000 ==> 1968, 2460000.0 ==> 2023)
        # - Should cover the required span of observations: TO BE CHECKED
        self.JDlist = np.arange(2440000., 2460000.0, 1.0)

    def _fetch_data_directory(self):
        '''
        Returns the default path to the directory where data will be downloaded.
        
        By default, this method will return ~/.mpchecker2/data
        and create this directory if it does not exist.
        
        If the directory cannot be accessed or created, then it returns the local directory (".")
        
        Returns
        -------
        data_dir : str
        Path to location of `data_dir` where data (FITs files) will be downloaded
        '''
            
        data_dir = os.path.join(os.path.expanduser('~'), '.mpchecker_data')

        # If it doesn't exist, make a new data directory
        if not os.path.isdir(data_dir):

            try:
                os.mkdir(data_dir)

            # downloads locally if OS error occurs
            except OSError:
                warnings.warn('Warning: unable to create {}. '
              'Download directory set to be the current working directory instead.'.format(data_dir))
            data_dir = '.'

        return data_dir


# Need other class / func to get unique_integer_ID from unique_string / obs80 / tracklet_ID / whatever ...




class Tracklet(Base):
    '''
        What we need to do "precalculations" on an individual tracklet
    '''

    def __init__(self , observations , unique_integer_ID ):
        # Give access to "Base" methods
        super().__init__()
        # Read the observations
        self.parse_observations(observations)

    def parse_observations(self, observations ):
        '''
            read observational input (probably be in obs80-string formats)
            
            return 
            integer dates & HP for each observation (or just extremal observations?)
            also do RoM & angle
            
        '''
        pass

    def save_results():
        '''
            This should use the r4esults from parse_observations and store them appropriately in a nice file/database structure
        '''
        pass

    def remove_from_stored_results()
        '''
            We need some method to remove 
        '''
        pass


# THIS IS NOT A PRE-CALC
class Query(Base):
    '''
        How to perform a query for an input orbit-specification
    '''
    
    def __init__(self , cheby_dict ):

        # Get Nightly Healpix
        self._get_nightly_healpix()
    
        # Query pre-calcs for approximate matches
        self._query_precalc()
        
        # Perhaps a refinement step using RoM & AoM
        self._refine_using_RoMAoM()
        
        # Refine approx matches to get "PRECISE" matches (within specified tollerance)
        self._get_precise_matches()


    

    def _get_nightly_healpix(self):
        '''
            Perhaps not *EVERY* night, but rather every POPULATED night
            Evaluate the orbital position at integer JDs 
            Not sure what to do about orbital uncertainties / variant orbits
            
            Probably return a list of healpix for each JD
        '''
        return

    def _query_precalc()
        pass

    def._refine_using_RoMAoM()
        pass
    
    def._get_precise_matches()
        pass

    def get_results(self):
        ''' 
            return nicely formatted results
        '''

