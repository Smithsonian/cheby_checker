# -*- coding: utf-8 -*-
# sifter/sifter/query

'''
    --------------------------------------------------------------
    cheby checker.
    
    Jun 2021
    Matt Payne
    
    This module provides access to some parameters / funcs / classes
    that don't sit well anywhere else, and are pretty general and / or
    module-wide
    
    Currently this contains
    (1) "Base" Class
     - Mainly used to specify the hard-wired date ranges and sectors
       that this package will work within
    () ...
    
    --------------------------------------------------------------
    '''


# Import third-party packages
# --------------------------------------------------------------
import sys, os
import numpy as np
from astropy_healpix import HEALPix as hp


# Import neighboring packages
# --------------------------------------------------------------

# Define top-line parameters
# --------------------------------------------------------------

class Base():
    '''
        Parent class to hold fundamental definitions / settings used across ...
        ... various ChebyChecker components.
        
        [[ much of this originally in precalc.py ]]
        
    '''
        
    # --- Notes on the expected mapping of coord-components to position ...
    # --- ... required in order to allow functions ** and ** to work properly ...
    # --- These are fundamental to the proper operation of much of MSC.
    # ------------------------------------------
    '''
        >>> coord_names     = ['x','y','z','vx','vy','vz']
        0   1   2   3    4    5
        >>> covar_names     = [ "_".join([coord_names[i], coord_names[j]]) for i in range(len(coord_names)) for j in range(i,len(coord_names))  ]
        >>> covar_names
        ['x_x', 'x_y', 'x_z', 'x_vx', 'x_vy', 'x_vz', 'y_y', 'y_z', 'y_vx', 'y_vy', 'y_vz', 'z_z', 'z_vx', 'z_vy', 'z_vz', 'vx_vx', 'vx_vy', 'vx_vz', 'vy_vy', 'vy_vz', 'vz_vz']
        6      7      8      9      10      11      12     13     14      15      16      17     18      19      20      21       22       23       24       25       26
        *      *      *                             *      *                              *
    '''

    
    # --- Params / Constants ... ---------------
    # ------------------------------------------
    secsPerDay = 86400
    epsilon    = 1e-6
    
    # --- Healpix settings ---------------------
    # ------------------------------------------
    HP_nside    = 16
    HP_order    ='nested'
    HPix        = hp(nside=HP_nside, order=HP_order)
    HP_npix     = HPix.npix
    '''
    self.HP_nside=16
    self.HP_order='nested'
    self.hp = HEALPix(nside=self.HP_nside, order=self.HP_order)
    self.npix = self.hp.npix
    '''
    
    # --- Date settings ------------------------
    # ------------------------------------------
    
    # We will assume that all sectors are of length 32-days
    sector_length_days = 32
    
    # We will assume a sector is "covered" if there is data ...
    # covering (at least)
    #        +sector_gap < t < sector_length_days-sector_gap days
    # out of the
    #        0 < t < sector_length_days [==32]
    # total length of the sector.
    #  - I have no good justification for this number at present:
    #  - It should be viewed as a free parameter whose value should be experimentally determined.
    sector_gap = 4.0

    # We will assume that the earliest standard epoch to be used will be JD:2440000
    # - This will set the enumeration of the sectors
    # - Could/Should change this so that the end date is checked to be X-days in the future, or something similar
    standard_MJDmin = 2440000 # 1968
    standard_MJDmax = 2464000 # 2034

    # Default list of Julian Dates to use (2440000 ==> 1968, 2464000.0 ==> 2034)
    JDlist = np.arange(standard_MJDmin, standard_MJDmax+1, 1)
    
    
    # --- Data / File / DB settings ------------
    # ------------------------------------------
    # sqlite database specs ...
    db_dir         = '.cheby_checker_data'
    db_filename    = 'cheby_checker.db'
    

    # Data Function(s) ...
    # --------------------------------------------------------------

    def _fetch_data_directory(self):
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
        
        data_dir = os.path.join(os.path.expanduser('~'), self.db_dir)
        print(" *** WARNING TO MPC STAFF *** \n THIS DEVELOPMENTAL CODE IS SAVING TO THE USERS-DIRECTORY \n THIS SHOULD BE CHANGED TO A SINGLE LOCN ON MARSDEN \n *** END OF WARNING ***" )
        
        
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


    # Date Functions ...
    # --------------------------------------------------------------
    @staticmethod
    def map_JD_to_sector_number( JD_TDB , JD0):
        ''' This maps input integer days onto a "sector"
        
            I.e. it maps times into fixed-length sectors (starting from JD0),
            and returns the index of the sector that any input times belong to
        '''
        return ( np.asarray( JD_TDB ) - JD0).astype(int) // Base.sector_length_days
    
    @staticmethod
    def map_sector_number_to_sector_start_JD( sector_number, JD0):
        ''' Function to map a sector number to the starting julian date of that sector
        '''
        return JD0 + sector_number * Base.sector_length_days

    @staticmethod
    def map_JD_to_sector_number_and_sector_start_JD(JD_TDB, JD0):
        '''
            For a given JD_TDB, calculate the sector number it will be in
            Assumes a standard starting julian date
            
            inputs:
            -------
            JD_TDB : number or iterable (i.e. something that numpy can work with!)
            
            returns:
            -------
            sector_number : numpy array
            start_JDs : numpy array
                    
        '''
        sector_number   = Base.map_JD_to_sector_number(JD_TDB , JD0)
        return sector_number , Base.map_sector_number_to_sector_start_JD( sector_number, JD0)

    def map_JDtimes_to_relative_times(self, times_TDB):
        '''
            For a given JD, map the JD to the sector # and 
            the (relative) time within the sector (starting from 0 for each sector)
        '''
        sector_numbers, sector_JD0s = self.map_JD_to_sector_number_and_sector_start_JD( times_TDB, self.standard_MJDmin )
        sector_relative_times = times_TDB - sector_JD0s
        return sector_numbers , sector_relative_times
    
    def get_required_sector_dict(self, ):
        '''
            Understand the basic sector-size used in constructing cheby-coefficients
            - Use that to define a complete list of sectors spanning the JDs in JDlist
            Will look like {0: 2440000, 1: 2440032, ...}
            - N.B. The dictionary automatically handles making a "set" of the keys
        '''
        return {a:b for a,b in zip(*self.map_JD_to_sector_number_and_sector_start_JD( self.JDlist , self.standard_MJDmin ))}



# End Of File
# --------------------------------------------------------------
