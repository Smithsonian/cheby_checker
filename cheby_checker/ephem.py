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
    
    --------------------------------------------------------------
    '''


# Import third-party packages
# --------------------------------------------------------------
import sys, os

# Import neighboring packages
# --------------------------------------------------------------


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

        # We'll want to fetch the minimal sub-set of data from the database : establish the required sectors
        # Turn dates into minimal set of sectors
        sector_numbers = np.asarray(sorted(set( orbit_cheby.Base.map_JD_to_sector_number( self.obsJDs , orbit_cheby.Base.standard_MJDmin) ) ) )
        # Add-in the adjacent sectors (just in case)
        sector_numbers = np.unique(np.concatenate([sector_numbers-1,sector_numbers,sector_numbers+1]))
        # Ensure that all requested sectors are supported
        minAllowed,maxallowed = *orbit_cheby.Base.map_JD_to_sector_number( [orbit_cheby.Base.standard_MJDmin,orbit_cheby.Base.standard_MJDmax] , orbit_cheby.Base.standard_MJDmin)
        sector_numbers = sector_numbers[ (sector_numbers>=minAllowed) & (sector_numbers<=maxallowed) ]
        
        # Use the *get_specific_object()* function in PreCalc to return the required list of MSCs
        self.MSCs = [precalc.Precalc.get_specific_object(desig , sector_numbers = sector_numbers) for desig in self.designations]


    # Convenience functions to generate various ephemeris quantities
    # --------------------------------------------------------------
    def generate_sky_predictions( ):
        '''
           Get the on-sky RA,Dec & associated uncertainties
           At some point should probably extend to Rates-of-Motion
        '''

        for M in self.MSCs:
            
            # Get the RA & Dec
            RADEC    = self.M.generate_RaDec( self.obsJDs  , observatoryXYZ=self.obsXYZ)#, APPROX = False , DELTASWITCH = False, delta=np.array([0,0,0]))
    
            # Get the covariance matrix in RA & Dec
            covRADEC = self.M.covRaDec( self.obsJDs , self.obsXYZ )

            # Look at posn_check _generate_predictions) to populate the rest of this function 



