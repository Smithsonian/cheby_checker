# -*- coding: utf-8 -*-
# cheby_checker/cheby_checker/mpchecker2.py

"""
--------------------------------------------------------------
The main "checking" functionalities

Jul 2020
Matt Payne

As currently sketched-out
(i)     Pre-Calculations are handled in precalc.py
(ii)    The Check class provides the methods to undertake
        (a) pCheck == posn_check == Are the supplied
            observations consistent with the predicted posns
            for the named object
        (b) mpchecker == what objects are within the FoV of
            a given pointing
        (c) checkid == what known object would have a posn
            consistent with the supplied detections
        (d) checkidX == was per checkid, but allowing the
            mean anomaly to vary

*** AS OF 2020_07_31, MUCH OF THE CODE IN HERE SHOULD BE    ***
*** VIEWED AS PSEUDO-CODE, JUST SKETCHING OUT APPROXIMATELY ***
*** HOW IT COULD/SHOULD WORK. NO TESTING HAS BEEN DONE      ***

--------------------------------------------------------------
"""

# Import third-party packages
# --------------------------------------------------------------
import numpy as np


# Import neighboring packages
# --------------------------------------------------------------
from . import precalc
from .orbit_cheby import Base
from . import data_classes


# Main functional classes for ...
# --------------------------------------------------------------
class Check(Base):
    """
    Provide all methods associated with "Checking"
     - Get "shortlists" for night
     - Do detailed orbital advance
     - Get "refined" list of actual overlap with FoV
     - Associated refined list with detections/tracklets

     *** AS OF 2020_07_31, MUCH OF THE CODE IN HERE SHOULD BE    ***
     *** VIEWED AS PSEUDO-CODE, JUST SKETCHING OUT APPROXIMATELY ***
     *** HOW IT COULD/SHOULD WORK. NO TESTING HAS BEEN DONE      ***
    """

    def __init__(self,):
        """ Empty initialization """
        self.MSCs = None
    

    def posn_check(self, primary_unpacked_provisional_designation, detections, param_dict=None):
        """
        pCheck: Do ephemeris look-up & calculate residuals w.r.t. detections

        inputs:
        -------
        primary_unpacked_provisional_designation : string
         -
        detections : A single *Detections* object
         -

        returns:
        --------
        None
        """
        # Check parameters are as required

        assert isinstance(primary_unpacked_provisional_designation, str)
        assert isinstance(detections, data_classes.Detections)
        assert param_dict is None or isinstance(param_dict, dict)

        # Call Ephem (ephemeris-querying object)
        Eph = Ephem(primary_unpacked_provisional_designation, detections.obstime, observatoryXYZ=detections.pos)

        # Get the predicted sky-positions
        # - As written PCheck only for single desig, so can just expand out the returned dict which will be of length-1
        # - These 'predictions' will be a 'Detections' object (with some cols as zeroes)
        predictions = Eph.generate_sky_predictions()[unpacked_primary_provisional_designation]

        # Initialize & return Residuals object
        return Residuals(predictions, detections, param_dict=param_dict)

    
    def mpchecker(self, pointings_object):
        """
        We want to calculate which objects are in the FoV
        (or have some chance of being in the FoV) of a single pointing.

        We load precalculated quantities to help us understand which subset of
        known objects have any chance at all of being close to the FoV

        We advance the subset (using chebys) to find whether or not they
        actually are in the FoV

        Inputs
        ----------
        pointing : a Pointing object

        Returns
        -------
        MSCs: a list of MSC objects
        """
        # Check is *Vectorial* object
        # N.B. Pointings & Detections will pass as they inherit from Vectorial
        assert isinstance(pointings_object, data_classes.Vectorial)

        # Get precalculated nightly datasets
        #  - For a given JD, finds the objects that are in HPlist at the night-center
        MSCs = precalc.PreCalc().get_nightly_precalcs( pointings_object.obstime , pointings_object.HPlist )
        
        # Advance all objects of interest to specific time of pointing
        # Calculate topocentric unit-vector using coefficient-dict
        # Calculate sepn between objects & center of pointing
        #
        #  - *** Should make more specific : have some part of the their sky-plane-uncert in the FoV ***
        #
        # Check whether this is within the FoV / within the search-radius
        #    *** Seems like this calc could/should be a method within MSC ***
        MSCs = [M for M in MSCs if np.arccos(np.clip( np.dot(
                                                            M.generate_UnitVector( pointing.obstime , pointing.pos ),
                                                            pointing.unit_vec() ) , -1, 1))
                                    < pointing.radius ]
        return MSCs


    def checkid(self, detections ):
        """
        Find objects close to presented detections

        Inputs
        ----------
        detections : a Detections object

        Returns
        -------
        """
        # Check is *Detections* object
        assert isinstance(detections, data_classes.Detections)

        # Turn each detection into a "pointing" so that we can use mpchecker
        #
        # *** THIS SEEMS SLIGHTLY UNNECESSARY ***
        # *** MOST/ALL OF THE NECESSARY VARIABLES ARE SHARED ACROSS Pointing & Detections
        # *** SHOULD JUST COMBINED / INHERIT FROM COMMON PARENT
        # *** THEN ALLOW MPCHECKER TO WORK ON EITHER
        #
        # NB This is setting pointing_radius == np.max(detection.rmsra, detection.rmsdec)
        #    => mpchecker will be searching for objects within this radius
        pseudo_pointings_list = detections.to_pointings()
        
        # Call mpchecker on each pointing
        # *** We may want to make mpchecker work on lists/arrays of Pointing objects ***
        for p in pseudo_pointings_list:
            MSCs = self.mpchecker( p )
            
        # What's next ?!?!?!?
        # (1) So far we have no concept of a tracklet, just individual detections:
        #     Probably need to have them linked as tracklets...
        # (2) Then need to introduce criteria / parameters to say whether *all*
        #     observations in the tracklet are close-enough


    def checkidX(self, detections ):
        """ A version of checkid that allows varying mean anomaly / time """
        pass
