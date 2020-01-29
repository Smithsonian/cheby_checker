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
print(os.path.dirname(os.path.realpath(__file__)))
sys.path.append( os.path.dirname(os.path.realpath(__file__)) )
import sql



class Query(Base):
    '''
        How to perform a query for an input orbit-specification
    '''
    
    def __init__(self , cheby_dict ):
        
        # Get Nightly Healpix
        JD_list, HP_list = self._get_nightly_healpix()
        
        # Query pre-calcs for approximate matches
        self.approx_dictionary = self._query_precalc(JD_list, HP_list)
        
        # Later could have a refinement step using RoM & AoM
        # self._refine_using_RoMAoM()
        
        # Refine approx matches to get "PRECISE" matches (within specified tollerance)
        self.precise_dictionary = self._get_precise_matches(self.approx_dictionary)
    
    
    
    def _get_nightly_healpix(self):
        '''
            Perhaps not *EVERY* night, but rather every POPULATED night
            Evaluate the orbital position at integer JDs
            Not sure what to do about orbital uncertainties / variant orbits
            
            Probably return a list of healpix for each JD
            '''
        
        JD_HP_dict = {}
        
        return JD_HP_dict
    
    def _query_precalc(self, JD_list, HP_list):
        '''
            inputs:
            -------
            JD_list: list of integers
             - dates to search
            HP_list: list of integers
             - healpix to search
            
            return:
            -------
            approx_dictionary
             - keys == tracklet_name
             - values == tracklet_dictionaries 
             - same structure as the return from sql.query_tracklets
             
        '''
        approx_dictionary = {}
        for JD, HP in zip(D_list, HP_list):
            approx_dictionary.update( query_tracklets_jdhp(JD, HP) )
        return approx_dictionary
    
    #def._refine_using_RoMAoM()
    #    pass
    
    def._get_precise_matches(self, approx_dictionary, tolerance_dict):
        '''
            Refine the approx matches (from _query_precalc) 
            down to exact matches within the specified tollerances
            
            inputs:
            -------
            approx_dictionary
            - keys == tracklet_name
            - values == tracklet_dictionaries
            - same structure as approx_dictionary from _query_precalc
            
            return:
            -------
            precise_dictionary
            - keys == tracklet_name
            - values == tracklet_dictionaries
            - same structure as approx_dictionary
            
        '''
        precise_dictionary = {}
        
        for tracklet_name, tracklet_dict in combined_dictionary.items():
            # do something involving the contents of the tracklet_dict
            pass
        
        return precise_dictionary

    def get_results(self, precise_dictionary):
        '''
            return nicely formatted results (to be defined)
        '''
        return {} 

