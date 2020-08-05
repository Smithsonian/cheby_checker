# -*- coding: utf-8 -*-
# cheby_checker/cheby_checker/data_classes
'''
Data classes for Cheby Checker

NB (1)
    There may be a better place for some of these, but for
    now I'm just collecting them together here.
NB (2)
    There may be some value to combining the Pointing &
    Detections into a single class, or deriving them from a common
    Direction/Vector object
NB (3)
    As of python 3.7 (?) data-classes exist as a nice, separate,
    simple class. However marsden only has 3.6, so I don't want
    to rely on that at the moment.
NB (4)
    Related to NB (3), some of the class definitions below are a
    bit "hacky" as I experiment with what will work best, and as
    I discover new (to me) ways to make classes work efficiently
    
'''

# Import third-party packages
# --------------------------------------------------------------
from collections import namedtuple 
import numpy as np
import sys
from astropy_healpix import healpy
#from functools import lru_cache

# Import neighboring packages
# --------------------------------------------------------------
try:
    from orbit_cheby import orbit_cheby
except ImportError:
    from . import orbit_cheby


# Simple dictionary-definitions of sets of variable names
# --------------------------------------------------------------

# Set up a dictionary to hold some variable names.
# - Add in some convenience-variable names
var_names = { 'Convenience' : ['pos','unit_vec','HP','HPlist'] }

# Common variables for Pointing & Detections: Place in a "Vectorial" (parent) object
var_names['Vectorial'] = ['obstime','ra','dec','pos1','pos2','pos3']

# Variables for "Pointing"
additional_pntng_var_names = ['radius']
var_names['Pointings'] = var_names['Vectorial'] + additional_pntng_var_names

# Variables for "Detections"
additional_detn_var_names = [ 'rmstime','rmsra','rmsdec', 'mag','rmsmag']
var_names['Detections'] = var_names['Vectorial'] + additional_detn_var_names
                                           
# Variables for "Residuals"
var_names['Residuals'] =    ['offsetRA',
                            'offsetDec',
                            'offsetTot',
                            'maxUncOrb',
                            'maxUncDet',
                            'maxUnc',
                            'badBool']

# Set up a name:number mapping between variables and posn in array
def var_map(vars):
    return {v:n for n,v in enumerate(vars)}
var_maps = {k:var_map(vars) for k,vars in var_names.items()}

# Data Class Definitions
# "Vectorial", "Pointing", "Detections", "Residuals"
# --------------------------------------------------------------

class Vectorial(orbit_cheby.Base):
    '''
        May use as a (hidden) parent class from which
        Pointing & Detections will inherit
        
        NB (1): Using __slots__ to get rid of (slowish) __dict__
        NB (2): Overwriting __getattribute__ to provide
                custom slices & combos (see below)
        NB (3): Want an intrinsically array-based method,
                but with the ability to take namedtuple-like
                slices of the array
    '''
    __slots__ =  'iama', 'arr', *var_names['Convenience'], *var_names['Vectorial']

    def __init__(self,):
        self.iama = 'Vectorial'

    # Provide means to initialize 2D array of quantities
    def _custom_array_init(self, arg ):

        # Allow set-up of empty array of length-N
        if isinstance(arg, int ):
            self.arr = np.zeros( (arg, len(var_names[self.iama])) )
         
        # This will interpret the following
        # - Single Iterable of length == len(var_names[self.iama])
        # - Iterable of iterables (inner iterables of length len(var_names[self.iama]))
        # NB
        # - I could have used a structured-array or record-array, but I find the initialization
        # - to be annoying as hell, so gave up. Insted I am handling attributable access to arr
        # -  via __getattribute__
        else:
            self.arr = np.atleast_2d( arg ).astype(float)
            assert self.arr.ndim == 2 ,                               f'Problem:self.arr.ndim={self.arr.ndim}'
            assert self.arr.shape[1] == len(var_names[self.iama]),    f'Problem:self.arr.shape={self.arr.shape}'
            
        return self.arr

    # Overwrite string method
    def __str__(self):
        return f'{self.iama}: obstime={self.obstime} ra={self.ra} dec={self.dec} pos={self.pos1,self.pos2,self.pos3}'

    # Provide convenient access to additional calculated variables
    def __getattribute__(self, name):
        """
            Overwriting default because
            (a) want to allow array slicing
            (b) want to provide convenience quantities
                ('pos', 'unit_vec', HP, HPlist, ...)
            """
        # Default behaviour for these two ...
        if name in ['iama', 'arr']:
            return object.__getattribute__(self, name)
            
        # Allow array slices via name, using the name:number mapping
        elif name in var_maps[self.iama]:
            return self.arr[:,var_maps[self.iama][name]]

        # The next ~four are dependent on the above array slicing
        elif name == 'pos':
            return np.array( [self.pos1, self.pos2, self.pos3] )
        elif name == 'unit_vec' :
            cd = np.cos(self.dec)
            return np.array( [ cd*np.cos(self.ra), cd*np.sin(self.ra), np.sin(self.dec)] )
        elif name == 'HP':
            return healpy.vec2pix(self.HP_nside, self.unit_vec[0], self.unit_vec[1], self.unit_vec[2], nest=True if self.HP_order=='nested' else False )
        elif name == 'HPlist':
            return healpy.query_disc(nside=self.HP_nside , vec=self.unit_vec, radius=self.radius, nest=True if self.HP_order=='nested' else False , inclusive=True)

        # For anything else, again try default
        else:
            # Default behaviour
            return object.__getattribute__(self, name)


 
class Pointings(Vectorial):
    ''' The pointing (exposure) which is the input to MPChecker
        
        Store in class centered around numpy array.
        Inherits from Vectorial
        Uses Vectorial methods to provide namedtuple-like slicing
        
        Need to initialise using
         - Single Iterable of length == len(pntng_var_names)
         - Iterable of iterables (inner iterables of length len(pntng_var_names))
        where pntng_var_names.keys = ( obstime, ra, dec, pos1, pos2, pos3, radius)
        
        ### ,namedtuple('Pointing', pntng_var_names.keys())
    '''
    
    # Add in pntng_var_names in addn to those from Vectorial
    __slots__ = *var_names['Pointings'],

    # Initialize via method provided in Vectorial
    def __init__(self, arg):
        self.iama = 'Pointings'
        self._custom_array_init( arg )





class Detections(Vectorial):
    ''' Store detections in class centered around numpy array.
        Inherits from Vectorial
        Uses Vectorial methods to provide namedtuple-like slicing
        
        Need to initialise using
         - Single Iterable of length == len(detn_var_names)
         - Iterable of iterables (inner iterables of length len(detn_var_names))
        where detn_var_names.keys = ( obstime, ra, dec, pos1, pos2, pos3, 'rmstime','rmsra','rmsdec', 'mag','rmsmag')

    '''
    
    # Add in detn_var_names ( & pntng_var_names) in addn to those from Vectorial
    __slots__ = *var_names['Detections'], *var_names['Pointings'],
    
    # Initialize via method provided in Vectorial
    def __init__(self, arg):
        self.iama = 'Detections'
        self._custom_array_init( arg )
        
    def __getattribute__(self, name):
        """ Attribute override to provide an equivalent 'radius' for a detection
        """
        if name == 'radius':
            # If a detection is used as a pointing, we need a 'radius' quantity
            return np.max(detections.rmsra[n], detections.rmsdec[n])
        else:
            # Default behaviour
            return object.__getattribute__(self, name)




class Residuals:
   ''' For the analysis of residuals: differences between Detections & Predictions
        
       Stores as a lean, mean, np.array-machine ...
       Allows slices using names (like a namedtuple)
       
       *** I REMAIN UNCERTAIN AS TO WHETHER THERE'S ANY POINT HAVING THIS AS ***
       *** A SEPARATE CLASS, OR WHETHER IT SHOULD BE A METHOD ON DETECTION   ***
   '''
   
   #__slots__ = ( *resid_var_names )
   
   def __init__(self, detections, predictions):
       """
           
           inputs:
           -------
           
       """
       
       # Check we have *Detections* objects of same shape
       assert isinstance(detections, Detections)
       assert isinstance(predictions, Detections)
       assert detections.D.shape == predictions.D.shape

       # Default evaluation parameters
       self.param_dict = {
           'nSigOrb': 2, # Some parameter that I am likely to interpret as a multiple of the diagonal of the orbital covariance matrix
           'nSigDet': 2, # Some parameter that I am likely to interpret as a multiple of the detection uncertainty
       }
       # Import any supplied parameters
       if param_dict is not None and isinstance(param_dict, dict) :
           self.param_dict.update( param_dict )

       # Evaluate the residuals & return array object
       self.R = self._evaluate_residuals(detections.D, predictions.D)
           
   def __getattribute__(self, name):
       """
           Overwriting default because
           (a) always want to take fresh slice (in case R is ever updated)
           (b) want to allow named access to the quantities in resid_var_names
           """
       if name in resid_var_names:
           return self.R[:,resid_var_names[name]]
       else:
           # Default behaviour
           return object.__getattribute__(self, name)


   def _evaluate_residuals(self, detections, predictions):
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
       # Calculate offsets in RA & Dec between predictions & detections
       # Could use astropy ... https://github.com/astropy/astropy/issues/4209
       nRA   = detn_var_names['ra']
       nDec  = detn_var_names['dec']
       nURA  = detn_var_names['rmsra']
       nUDec = detn_var_names['rmsdec']
       offsetRA  = (predictions.D[:,nRA]  - detections.D[:,nRA]  ) * np.cos( detections.D[:,nDec] )
       offsetDec = (predictions.D[:,nDec] - detections.D[:,nDec] )
       offsetTot = np.sqrt( offsetRA**2 + offsetDec**2 )
       
       '''
           # Could try to follow the approach suggested by M.Pan ...
           
           "But for other objects I had thought maybe one could do the whole comparison in a coordinate grid with axes parallel to and perpendicular to the major/minor axes of the error ellipses on the positions calculated from the orbit.
           That is, project the RA and Dec RMS's of the new observations onto axes parallel and perpendicular to the error ellipse major/minor axes, and check that the differences between the observed and predicted positions in this error-ellipse centered coordinate system are smaller than the lengths of the major/minor axes.
           That would be more like error rectangles rather than error ellipses, but would be better at preserving any very long aspect ratios in the error ellipses."


            # But then there is the separate discussion with Federica concerning the paper from Milani,
            in which he doesn't even bother with the RA,DEC covariance matrix, but instead constructs
            something like a convex hull based on the positions of a selection of orbits drawn from the
            covariance ellipsoid in cartesian/keplerian coords.
             - in some file / notbook there are some experiments as I attempt to make something like
            a uniform or random sampling from the surface of a 6D sphere.
           
       '''
       
       # Comparing offsets/residuals to the allowed tolerances : I have put very little thought into this.
       # - Seems something like a circular approximation
       # - Should be tested & investigated in far more detail to ensure it is better than the old pcheck
       #
       # Here I get the max-unc in the orbit-predictions (either RA or Dec) to use as an equivalent circular radius for each detection
       maxUncOrb = np.min( self.predictions[:,[nURA,nUDec]], axis=1 )
       # Here I get the max-unc in the detections (either RA or Dec) to use as an equivalent circular radius for each detection
       maxUncDet = np.min( self.predictions[:,[nURA,nUDec]], axis=1 )
       # Here I construct a max-allowed sepn (factoring in allowed scaling from self.param_dict)
       maxUnc = np.sqrt( (maxUncOrb * self.param_dict['nSigOrb'])**2 + (maxUncDet * self.param_dict['nSigDet'])**2 )
       # Check whether the sepn is larger than allowed
       badBool = offsetTot > maxUnc
       
       # Put everything into a *Residuals* object
       # - Could / should reformat this as necessary
       self.residuals = Residuals( np.array(offsetRA, offsetDec, offsetTot, maxUncOrb, maxUncDet, maxUnc, badBool).T )
       return self.results_dict
   
   
   







