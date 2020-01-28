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
from functools import lru_cache
import json


# Import neighboring packages
# --------------------------------------------------------------




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

