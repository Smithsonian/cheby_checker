# -*- coding: utf-8 -*-
# mpchecker2/mpchecker2/mpchecker2.py

'''
--------------------------------------------------------------
Main MPChecker2 module.

Feb 2020
Matt Payne

To contain all (high-level) functions required to predict
minor planets within the field-of-view of a pointing

As currently sketched-out 
(i)     Pre-Calculations are handled in mpchecker2/mpchecker2/precalc.py
(ii)    MPCHECKER class takes POINTINGs as input and then
        calculates the SUSPECTS within the FoV
(iii)   IDENTIFY class takes SUSPECTS and compares them to data
        (either EXPOSURE or TRACKLETS) and returns IDENTIFICATIONS

--------------------------------------------------------------
'''


# Import third-party packages
# --------------------------------------------------------------
import sys
import numpy as np


# Import neighboring packages
# --------------------------------------------------------------
import precalc
#import orbit_cheby as eph


# Data classes for MPChecker2
# --------------------------------------------------------------

class Pointing(namedtuple('Pointing', 'obstime ra dec radius pos1 pos2 pos3')):
    '''
        To contain the pointing (exposure) which is the input to MPChecker
        
        Will have to contain information like ...
        date,
        obscode/observatory-posn,
        (RA,Dec)/Unit-Vector,
        radius,
        ...
        
        Adding methods to namedtuple (also don't use dict)
        http://zetcode.com/python/namedtuple/
        
    '''
    
    # Magic to prevent use of dict by class ...
    __slots__ = ()
    
    #@property
    #def UnitVec(self):
    #    return sqrt((self.x ** 2 + self.y ** 2))
    
    # Define a string
    def __str__(self):
        return f'Pointing: obstime={self.obstime}  ra={self.ra}  dec={self.dec} pos={self.pos1,self.pos2,self.pos3}'








# Main functional classes for MPChecker2
# --------------------------------------------------------------

class MPCHECKER:
    '''
        Provide all methods associated with MPChecker
         - Get "shortlists" for night
         - Do detailed orbital advance
         - Get "refined" list of actual overlap with FoV
         - Associated refined list with detections/tracklets
        
    '''
    
    
    def __init__(self, ):
        '''
            ...
        '''
        self.suspects = None
    

    def get_objects_which_overlap_FoV( pointing_variable ):
        '''
            We want to calculate which objects are in the FoV 
            (or have some chance of being in the FoV) of a single pointing. 
            
            We load precalculated quantities to help us understand which objects have 
            any chance at all of being close to the FoV
            We
            
            Parameters
            ----------
            
            Returns
            -------
            suspects    :   SUSPECTS-object
            
            Examples
            --------
            >>> ...
        
        '''

        # Make sure pointing_variable is a list-of-pointings
        list_of_pointings = self._rectify_pointings(pointing_variable)

        # Get precalculated nightly datasets
        # - result_dict: key=HP ;  value=list-of-coeff-dictionaries
        for p in list_of_pointings:
            
            # For a given JD, finds the object-integers for each HP in a list of HPs
            # - These are the objects in close HP. Subsequently to be refined.
            MSCs = precalc.PreCalc().get_nightly_precalcs( p.JDutc , p.HPlist )
        
            # Advance all objects of interest to specific time of pointing
            # Calculate topocentric unit-vector using coefficient-dict
            # Calculate sepn between objects & center of pointing
            #  - Should make more specific : have some part of the their sky-plane-uncert in the FoV
            # Check whether this is within the FoV / within the search-radius
            MSCs = np.array( [M for M in MSCs if np.arccos(np.clip( np.dot( M.generate_UnitVector( p.JDutc , np.array([p.pos1, p.pos2, p.pos3] ) ), p.UnitVec() ) , -1, 1)) < p.radius ] )


        return SUSPECTS



    def _rectify_pointings(self,  pointing_variable ):
        '''
            Private method called to rectify the input to "get_objects_which_overlap_FoV()"
            
            inputs:
            -------
            pointing_variable: pointing or list-of-pointings
            
            returns:
            --------
            list_of_pointings : list-of-pointings
            -
        
        '''
        # If singular quantity, make into lists
        if isinstance(pointing_variable, POINTING):
            list_of_pointings = [POINTING]
            
        # If non-singular, check plausibly formatted
        # (N.B. Not doing detailed content checks at this point)
        elif    isinstance(pointing_variable, (list, np.ndarray)) and \
            np.all([ isinstance(_, POINTING) for _ in pointing_variable]):
            list_of_pointings = pointing_variable
        else:
            sys.exit('Cannot process input of type %r in MPCHECKER ' % ( type(pointing_variable)) )
                        
        # return everything in list form
        return list_of_pointings



class IDENTIFY:   ## <<-- What is the difference between identification and attribution ?
    '''
    To match SUSPECTS from POINTINGs to
    a set of detections and/or tracklets
    '''
    
    # Add some content ...
    def __init__(self , suspects, data):
        '''
            
        Parameters
        ----------
        suspects    :   single SUSPECTS-object or list of SUSPECTS-objects
            These are the results of running MPChecker
        data        :   single EXPOSURE-object or list of TRACKLET-objects
            These are the data we wish to identify (attribute?) 
            with previously known object-orbits
            
        Returns
        -------
        ids         :   IDENTIFICATIONS-object
        
        
        Examples
        --------
        >>> ...
        
        '''
        # Single detections ...
        if isinstance(suspects, SUSPECTS) and isinstance(data, exposure):
            ids = identify_SUSPECTS_from_multiple_pointings_with_tracklets_from_multiple_exposures(suspects, data)
        # Multiple detections in tracklets ...
        elif isinstance(suspects, list) and isinstance(data, trackletList):
            ids = identify_SUSPECTS_from_multiple_pointings_with_tracklets_from_multiple_exposures(suspects, data)
        # Erroroneous input ...
        else:
            sys.exit("Incompatible data types input to IDENTIFY class")

        # Does any additional reformatting or proocessing need to be done ?
        
        return ids


    def identify_SUSPECTS_from_a_single_pointing_with_detections_from_a_single_exposure(suspects, exposure):
        '''
        
        Parameters
        ----------
        suspects    :   SUSPECTS-object
        exposure    :   EXPOSURE-object
        
        Returns
        -------
        
        Examples
        --------
        >>> ...
        
        '''
            
        return IDENTIFICATIONS()

    def identify_SUSPECTS_from_multiple_pointings_with_tracklets_from_multiple_exposures(suspectsList, trackletList):
        '''
        Identify SUSPECTS from multiple POINTINGS with a list of TRACKLETS
        Dates of the pointings need to correspond to the dates of the tracklet detections

        Parameters
        ----------
        suspectsList    :   list of SUSPECTS-objects
        trackletList    :   list of TRACKLET-objects

        Returns
        -------

        Examples
        --------
        >>> ...

        '''

        return IDENTIFICATIONS()





