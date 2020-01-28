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
import astropy_healpix.healpy as hp
from functools import lru_cache
import json


# Import neighboring packages
# --------------------------------------------------------------
sys.path.append( os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'sifter') )
import sql



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
        #self.HP     = '_healpix.txt'

        # Healpix ...
        self.HP_nside=16
        self.HP_order='nested'
        self.npix = hp.nside2npix(self.HP_nside)

        # Default list of Julian Dates to use (2440000 ==> 1968, 2460000.0 ==> 2023)
        # - Should cover the required span of observations: TO BE CHECKED
        # *** PROBABLY DON'T EVEN NEED THIS ...
        #self.JDlist = np.arange(2440000., 2460000.0, 1.0)
        
        # sqlite database specs ...
        self.db_filename = 'sifter.db'


    
    def _fetch_data_directory(self, ):
        '''
            Returns the default path to the directory where data will be downloaded.
            
            By default, this method will return ~/.sifter_data/data
            and create this directory if it does not exist.
            
            If the directory cannot be accessed or created, then it returns the local directory (".")
            
            Returns
            -------
            data_dir : str
            Path to location of `data_dir` where data (FITs files) will be downloaded
            '''
        
        data_dir = os.path.join(os.path.expanduser('~'), '.sifter_data')
        
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






class Tracklet(Base):
    '''
        What we need to do "precalculations" on an individual tracklet
    '''
    
    def __init__(self , observations=None ):
        
        # Give access to "Base" methods
        super().__init__()
        
        # connect to db
        self.conn = sql.create_connection( sql.fetch_db_filepath() )
    
        # if observations supplied, parse them ...
        if observations != None:
            parse_observations(observations)
    
    def parse_observations(self, observations ):
        '''
            read observational input (probably be in obs80-string formats)

            Inputs:
            -------
            observations : list ???
             - in obs80 format ???

            Returns
            -------
            JD: integer
             - date
            HP: integer 
             - healpix
            tracklet_name: string 
             - unique name for tracklet
            tracklet_dict: dictionary
             - container for all other data 
             - observations, RoM, angles, ..., ...,  ... 
             - should be everything required for subsequent detailed calculations
            '''
        
        # ***RANDOM DATA ***
        JD = np.random.randint(1000)
        HP = np.random.randint(self.npix)
        
        alphanumeric = np.array( list('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789') )
        tracklet_name = ''.join(list(np.random.choice( alphanumeric , 10 )))
        tracklet_dictionary = {'qwe' : np.random.randint(333) ,
                            'zxcv': ''.join(list(np.random.choice( alphanumeric , 7 )))}
        
        return JD, HP, tracklet_name, tracklet_dictionary
    
    def save_tracklet(self, JD, HP, tracklet_name, tracklet_dictionary):
        '''
            This should use the results from parse_observations and store them appropriately in a nice file/database structure
            
            Inputs:
            -------
            JD : integer
            - day
            
            HP : integer
            - healpix
            
            tracklet_name: string ?
             - unique identifier for tracklet
            
            tracklet_dictionary:
             - all data that we want to save
            
            Returns
            -------
            
        '''
    
        # upload data
        return sql.upsert_tracklet(self.conn, JD,  HP, tracklet_name, tracklet_dict)
    
    def delete_tracklet(self, tracklet_name):
        '''
            We need some method to remove
        '''
        return sql.delete_tracklet(self.conn, tracklet_name)

