# -*- coding: utf-8 -*-
# cheby_checker/cheby_checker/ephem

'''
    --------------------------------------------------------------
    cheby_checker's / ephemeris module.
    
    Jun 2020
    Matt Payne
    
    This module/class is intended to be a fairly light wrapper
    around the PreCalc & MSC classes.
    
    Given a designation (list?) and a set of desired 
    observation-times and observatory-locations, the Ephem class
    gets the data from the database (using PreCalc) and evaluates
    the object positions, etc
    
    Various convenience functions should be developed to
    provide the desired data in the required formats for (e.g.) 
    the online ephemeris service, as well as other downstream 
    functionality (pCheck, MPChecker, etc)
    
    ** As of Aug 2020, there doesn't seem to be a super-strong  **
    ** justification for having this as a stand-alone class.    **
    ** -->> Could just have as a function in some other class?  **
    
    --------------------------------------------------------------
    '''


# Import third-party packages
# --------------------------------------------------------------
import sys, os
import numpy as np

# Import neighboring packages
# --------------------------------------------------------------
try:
    import orbit_cheby
except ImportError:
    from . import orbit_cheby

from data_classes import detn_var_names, Detections

# Ephemeris Object
# --------------------------------------------------------------

class Ephem():
    '''
        ...
    '''
    

    def __init__(self, designations, times, observatoryXYZ = None , obsCode = None):
        '''
           Initializing ephemeris call
           Will get object data from database & instantiate MSC(s)
        '''
    
        # Rectify designations to ensure in list
        self.designations = list(np.atleast_1d(designations))
        
        # Get the observatory positions
        assert observatoryXYZ is not None or obsCode is not None , 'Both observatoryXYZ == None and obsCode == None'
        self.observatoryXYZ   = np.asarray(observatoryXYZ) if obsCode is None else ObsPos().get_heliocentric_equatorial_xyz(times,
                                                                                                                            obsCode=obsCode,
                                                                                                                            verbose=False)

        # Ensure that the times & positions have appropriate dimensions
        self.times              = np.asarray(times)
        assert self.observatoryXYZ.shape == (self.times.shape[0],3)


        # Establish the minimal required set of sectors to be fetched from the database :
        # Turn dates into minimal set of sectors
        sector_numbers = np.asarray(sorted(set( orbit_cheby.Base.map_JD_to_sector_number( self.obsJDs , orbit_cheby.Base.standard_MJDmin) ) ) )
        # Add-in the adjacent sectors (just in case)
        sector_numbers = np.unique(np.concatenate([sector_numbers-1,sector_numbers,sector_numbers+1]))
        # Ensure that all requested sectors are supported
        minAllowed,maxallowed = orbit_cheby.Base.map_JD_to_sector_number( [orbit_cheby.Base.standard_MJDmin,orbit_cheby.Base.standard_MJDmax] , orbit_cheby.Base.standard_MJDmin)
        sector_numbers = sector_numbers[ (sector_numbers>=minAllowed) & (sector_numbers<=maxallowed) ]
        
        # Use MSC_Loader to return the required list of MSCs
        self.MSCs = orbit_cheby.MSC_Loader(
                                           FROM_DB = True ,
                                           primary_unpacked_provisional_designations = self.designations,
                                           sectors = sector_numbers)


    # Convenience functions to generate various ephemeris quantities
    # --------------------------------------------------------------
    def generate_sky_predictions( ):
        '''
           Get the on-sky RA,Dec & associated uncertainties
           
           At some point should probably extend to Rates-of-Motion
           
        '''
        self.prediction_dict = {}
        for M in self.MSCs:
            
            # Get the RA & Dec
            RADEC    = self.M.generate_RaDec( self.obsJDs  , observatoryXYZ=self.obsXYZ)
            #, APPROX = False , DELTASWITCH = False, delta=np.array([0,0,0]))
    
            # Get the covariance matrix in RA & Dec
            covRADEC = self.M.covRaDec( self.obsJDs , self.obsXYZ )

            # Create pseudo-detections from the calculated positions ...
            # Set-up an empty Detections instance of the correct length ...
            Ds = Detections(len(self.times))
            
            # (ii) copy across the times & observatory-positions
            Ds.D[:, detn_var_names['obstime']] = self.times
            Ds.D[:, np.array([detn_var_names['pos1'], detn_var_names['pos2'],detn_var_names['pos3']]) ] = self.observatoryXYZ
            
            # (iii) populate the ra, dec & associated uncert
            Ds.D[:,detn_var_names['ra']]     = RADEC[:,0]
            Ds.D[:,detn_var_names['dec']]    = RADEC[:,1]
            Ds.D[:,detn_var_names['rmsra']]  = covRADEC[:,11] ## <<-- These positions need to be updated !!!!!
            Ds.D[:,detn_var_names['rmsdec']] = covRADEC[:,22]

            # Store the Detections object in a dictionary
            self.prediction_dict[ Ds ]

        return self.prediction_dict




