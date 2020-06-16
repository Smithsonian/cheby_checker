# -*- coding: utf-8 -*-
# sifter/sifter/query

'''
--------------------------------------------------------------
Jun 2020
Matt Payne & Margaret Pan & Mike Alexandersen

Cheby_Checker's position checking module
Intended to be a replacement for the MPC's old "pcheck" routine

Given a primary-provisional-desig, and some observations, the 
task is to 
(a) calculate the residuals between the expected position and the observations, and
(b) decide whether these residuals are "statistically acceptable"

Will do this using the orbit_cheby module & the sql-look-up functionalities


~~~~~~~~~~~~
To Do
 - Import whatever is required
 - Define class 
 - Specify actions:

~~~~~~~~~~~~


--------------------------------------------------------------
'''


'''
    Want there to be a detection object / namedtuple 
    Detection must contain obs-posn @ time of obs
    This will probably go into its own file ...
'''
Detection = namedtuple('Detection', ['obstime',
                                     'rmstime',
                                     'pos1',
                                     'pos2',
                                     'pos3',
                                     'ra',
                                     'dec',
                                     'rmsra',
                                     'rmsdec',
                                     'mag',
                                     'rmsmag'
                                     ])



class PCheck():



    def __init__(self, primary_unpacked_provisional_designation, detections , param_dict = None ):
        '''
        inputs:
        -------
        primary_unpacked_provisional_designation : string 
         - 
        detections : list of detection objects
        
        '''
        # Default evaluation parameters
        self.param_dict = {
            'nSigOrb': 3, # Some parameter that I am likely to interpret as a multiple of the diagonal of the orbital covariance matrix
            'nSigDet': 3, # Some parameter that I am likely to interpret as a multiple of the detection uncertainty
                            }
        if param_dict is not None and isinstance(param_dict, dict) :
            self.param_dict.update( param_dict )

        # Required inputs
        self.primary_unpacked_provisional_designation = primary_unpacked_provisional_designation
        self.detections                               = np.asarray(detections)

        # Probably of general use to extract arrays of quantities ...
        self.obsJDs = A[:,0]
        self.obsXYZ = A[:,2:5]

        # We'll want the MSC object, so just go ahead and populate it using the database
        # - Probably just want to get a sub-set of the data from the database to save time
        sector_numbers = orbit_cheby_Base.map_JD_to_sector_number( self.obsJDs , orbit_cheby_Base.standard_MJDmin)
        self.M = precalc.Precalc.get_specific_object(primary_unpacked_provisional_designation , sector_numbers = sector_numbers)
    
        # We'll probably want the predicted-detections as well, so just make them too
        self._generate_predictions()
    
        # And then we'll want to evaluate the detections against the orbit predictions
        self.evaluate_residuals()
    


    def _generate_predictions():
        '''
            From the orbit, generate the locn of the object at the time of each input detection
        '''
        # Get the RA & Dec
        RADEC    = self.M.generate_RaDec( self.obsJDs  , observatoryXYZ=self.obsXYZ)#, APPROX = False , DELTASWITCH = False, delta=np.array([0,0,0]))
        
        # Get the covariance matrix in the RA & Dec
        covRADEC = self.M.covRaDec( self.obsJDs , self.obsXYZ )

        # Create pseudo-detections from the calculated positions ...
        # (i) create arrray of zeroes of the correct size
        self.predictions = np.zeros( self.detections.shape )
        # (ii) copy across the times & observatory-positions
        indicees = np.array([0,2,3,4])
        self.predictions[:,indicees] = self.detections[:,indicees]
        # (iii) populate the ra, dec & associated uncert
        self.predictions[:,5] = RADEC[:,0]
        self.predictions[:,6] = RADEC[:,1]
        self.predictions[:,7] = covRADEC[:,11]
        self.predictions[:,8] = covRADEC[:,22]

        return self.predictions


    def evaluate_residuals():
        '''
            Compare each detection to the expected value from the orbit 
            Is the residual (difference between them) acceptable?
            
            # I think that the proper solution to this involves understanding whether two ellipses are separated
            # This is doable but complicated:
            # https://www.geometrictools.com/Documentation/IntersectionOfEllipses.pdf
            # http://www.iri.upc.edu/files/scidoc/1852-New-algebraic-conditions-for-the-identification-of-the-relative-position-of-two-coplanar-ellipses.pdf
            # For the sake of rapid development (and execution!) I'll just stick to (what I believe is) an effectively circular approximation

            
            *** SOME VERSION OF THIS LIKELY TO BE USED IN MPCHECKER / CHECKID ***
            *** MIGHT WANT TO LOCATE THIS IN SOME MORE OBVIOUS GENERAL LOCN *****
            
            Could also try and connect this to the logic Federica uses to identify bad-tracklets 
            after doing a full refit
            
        '''
        # Could use astropy ... https://github.com/astropy/astropy/issues/4209
        nRA = 5
        nDec = 6
        nURA  = 7
        nUDec = 8
        offsetRA  = (self.predictions[:,nRA]  - self.detections[:,nRA]  ) * np.cos( self.detections[:,nDec] )
        offsetDec = (self.predictions[:,nDec] - self.detections[:,nDec] )
        offsetTot = np.sqrt( offsetRA**2 + offsetDec**2 )

        # Comparing offsets/residuals to the allowed tolerancess : I have put littel thought into this.
        # - Seems something like a circular approximation
        # - Should be tested & investigated in far more detail to ensure it is better than the old pcheck
        #
        # Here I get the max-unc in the orbit (either RA or Dec) to use as an equivalent circular radius for each detection
        maxUncOrb = np.min( self.predictions[:,[nURA,nUDec]], axis=1 )
        # Here I get the max-unc in the detections (either RA or Dec) to use as an equivalent circular radius for each detection
        maxUncDet = np.min( self.predictions[:,[nURA,nUDec]], axis=1 )
        # Here I construct a max-allowed sepn (factoring in allowed scaling from self.param_dict)
        maxUnc = np.sqrt( (maxUncOrb * self.param_dict['nSigOrb'])**2 + (maxUncDet * self.param_dict['nSigDet'])**2)
        # Check whether the sepn is larger than allowed
        badBool = offsetTot > maxUnc

        # Put everything into a standard results dict
        # - Could / should reformat this as necessary
        self.results_dict = {
            'detections' : self.detections,
            'predictions': self.predictions,
            'offsetRA'   : offsetRA,
            'offsetDec'  : offsetDec,
            'offsetTot'  : offsetTot,
            'badBool'    : badBool,
        }
        return self.results_dict



