    # -*- coding: utf-8 -*-
# sifter/sifter/query

'''
    --------------------------------------------------------------
    sifter's / mpchecker's orbit_cheby module.
    
    Jan 2020
    Matt Payne & Margaret Pan & Mike Alexandersen
    
    This module provides functionalities to evaluate
    dictionaries of chebyshev-coefficients
    
    We are developing a standardized approach regarding
    orbit integration and subsequent interpolation using
    chebyshev-coefficients applied to 32-day sectors
    
    To contain all core functionality required to predict the
    position (e.g. Heliocentric Cartesian) and apparent position
    (e.g. RA, Dec) of minor planets, comets, satellites, etc
    - See https://drive.google.com/open?id=1F86lnaHL01xPAACX2SYVdOgvvLLnxfYD
    
    This does *NOT* do the underlying nbody integrations.
    This DOES fast interpolation using supplied chebyshev dictionaries
    
    
    N.B. There is a strong assumption that coordinates are BARYCENTRIC EQUATORIAL
    - I.e. we are assuming that the nbody integration has supplied cartesian
    coordinates in a BARYCENTRIC EQUATORIAL frame.

    Currently envisaged as a core component that will be drawn-upon by
    - Position_Check
    - MP_Checker
    - ID_Check
    - ...
    
    --------------------------------------------------------------
    '''


# Import third-party packages
# --------------------------------------------------------------
import sys, os
import numpy as np
import operator
from collections import OrderedDict, defaultdict
import astropy
from astropy_healpix import HEALPix as hp
from astropy_healpix import healpy
from functools import lru_cache
import json
import itertools
import numdifftools as nd

# Import neighboring packages
# --------------------------------------------------------------
import nbody_reader

# Define top-line parameters
# --------------------------------------------------------------

secsPerDay = 86400

# We will assume that all sectors are of length 32-days
sector_length_days = 32

# We will assume that the earliest standard epoch to be used will be JD:2440000
# - This will set the enumeration of the sectors
standard_MJDmin = 2440000 # 1968
standard_MJDmax = 2460000 # 2034

# Healpix settings
HP_nside    = 16
HP_order    ='nested'
HPix        = hp(nside=HP_nside, order=HP_order)
HP_npix     = HPix.npix

'''
>>> coord_names     = ['x','y','z','vx','vy','vz']
                        0   1   2   3    4    5
>>> covar_names     = [ "_".join([coord_names[i], coord_names[j]]) for i in range(len(coord_names)) for j in range(i,len(coord_names))  ]
>>> covar_names
['x_x', 'x_y', 'x_z', 'x_vx', 'x_vy', 'x_vz', 'y_y', 'y_z', 'y_vx', 'y_vy', 'y_vz', 'z_z', 'z_vx', 'z_vy', 'z_vz', 'vx_vx', 'vx_vy', 'vx_vz', 'vy_vy', 'vy_vz', 'vz_vz']
  6      7      8      9      10      11      12     13     14      15      16      17     18      19      20      21       22       23       24       25       26
  *      *      *                             *      *                              *
'''





class MSC_Loader():
    '''
        Multi-Sector Cheby -Loader Function
        
        Will load/create/return instances of MSC : Multi-Sector Cheby objects
        
        Will handle multi-particle instantiation: will return **LIST** of *MSC* objects
        
    '''

    # allow different init depending on source ...
    def __init__(self, **kwargs):
        print('INIT MSC_Loader...')
        
        # Initialization of default standard PARAMETERS / OPTIONS we use for chebyshevs, etc
        # - These may be overwritten by kwargs
        self.FORCE_DATES    = False                                 # : Sector Ranges
        self.TDB_init       = None                                  # : Cheby Span
        self.TDB_final      = None                                  # : Cheby Span
        self.FROM_FILE      = False                                 # : Ingest method
        self.filepath       = False                                 # : Ingest method
        self.FROM_ARRAY     = False                                 # : Ingest method
        self.unpacked_provisional_designations           = None     # : Ingest method
        self.times_TDB      = None                                  # : Ingest method
        self.statearray     = None                                  # : Ingest method
        self.FROM_DATABASE  = None                                  # : Ingest method
        
        # The list of MSCs that will be instantiated & returned
        self.MSCs           = []
        
        
        # Do a loop over the parameters and see if any attributes need updating (because they were supplied)
        for attribute, value in kwargs.items():
            if hasattr(self, attribute):
                setattr(self, attribute, value)

        # Allow initialization of ARRAYS / COEFFICIENTS from different sources
        # (i) From text file ( Mainly for development )
        if  self.FROM_FILE and self.filepath != None :
            self._populate_from_nbody_text(self.filepath)
        
        # (ii) From numpy array (from nbody)
        elif self.FROM_ARRAY and self.unpacked_provisional_designations != None and isinstance(self.statearray , (list,tuple,np.ndarray)) :
            self._populate_from_nbody_array(self.unpacked_provisional_designations, self.times_TDB, self.statearray)
        
        # (iii) From database (of pre-calculated chebyshevs)
        elif self.FROM_DATABASE:
            self._populate_from_database(**kwargs)
        
        # (iv) An empty instance ...
        else:
            self._generate_empty()
                


    def _generate_empty(self,  ):
        '''  '''
        print('Defaulting to the production of a list of empty MSCs')
        self.MSCs.append(MSC())



    # Functions to handle different instantiaion methods
    #  - Take input and generate the chebyshev coefficients
    # --------------------------------------------------------------

    def _populate_from_nbody_text( self, text_filepath ):
        '''
            Read data (assumed from an nbody integration) from a text file 
            Use it to generate an array of positions, etc
            Then continue on to use func *_populate_from_nbody_array()*
            
            Mainly of use during development
            
            inputs:
            -------
            text_filepath: string
            - filepath to data from nbody simulation 
            
            returns:
            --------
            boolean
            - whatever output is from *_populate_from_nbody_array()*
            
        '''
        # Read the nbody-json file
        name, times, states  = nbody_reader.parse_nbody_txt( text_filepath )

        # Use the create-from-array function to do the rest
        return self._populate_from_nbody_array(name, times, states)



    def _populate_from_nbody_array(self,
                                   unpacked_provisional_designations ,
                                   times_TDB,
                                   states,
                                   ):
        '''
            Will initialize MSC(s) from supplied arrays
            This handles the partitioning of data into multiple MSCs as required
            
            inputs:
            -------
            name: string
             - unpacked_provisional_identification of the object
            times_TDB: np.array
             - TDB/TT
            states: np.array
             -
            
            returns:
            --------
            
        '''
        
        # *** Strict standardized approach: Only sectors within standardized ranges will be generated ***
        # E.g. sector-zero ==>> standard_MJD0   -to-   standard_MJD0 + sector_length_days, ...
        #
        # (0) Catch problem of incomplete info
        if   self.FORCE_DATES and (self.TDB_init == None or self.TDB_final == None):
            sys.exit('Cannot force when TDB_init and TDB_final not supplied')
        
        # (1) Allow the user to force specific start/end dates
        elif self.FORCE_DATES and self.TDB_init != None and self.TDB_final != None:
            print('Forcing, but will still use standard %r-day blocks...' % sector_length_days )
            if  self.TDB_init < times_TDB[0] or self.TDB_final > times_TDB[-1]:
                sys.exit(' nbody data does not support the requested dates ')
        
        # (2) Generate the start/end dates programatically (if not supplied/forced by user)
        elif not self.FORCE_DATES:
            # Here we compare the supplied JDs to the standard end-points
            self.TDB_init  = int(max(standard_MJDmin , times_TDB[0] ))
            self.TDB_final = int(min(standard_MJDmax , times_TDB[-1] ))

        else:
            sys.exit('Do not know how to calculate limits in *_generate_multi_sector_cheby_dict_from_nbody_array()* ')




        # Check that the dimensionality of the coordinates is consistent with the number of names
        # - N.B. We expect states.shape = (Nt, Nc, Np) or (Nt,Nc), where
        #        Nt = Number of times
        #        Nc = Number of coords/components being fitted
        #        Np = Number of particles
        self.unpacked_provisional_designations = np.atleast_1d(unpacked_provisional_designations)
        # (i) Single Object
        if len(self.unpacked_provisional_designations) == 1 and states.ndim == 2 :
            pass
        # (ii) Multiple Objects
        elif len(self.unpacked_provisional_designations) > 1 and states.ndim == 3 and (len(self.unpacked_provisional_designations) == states.shape[-1]) :
            pass
        # (iii) Problem
        else:
            sys.exit('Inconsistent dimensionality : len(self.unpacked_provisional_designations) = %r and states.ndim = %r and states.shape = %r' % \
                     (len(self.unpacked_provisional_designations) , states.ndim ,  states.shape) )






        # Loop over each of the objects and create a MSC-object for each ...
        for i, unpacked in enumerate(self.unpacked_provisional_designations):
            
            # Get the slice of states corresponding to the particular named object
            state_slice = states if states.ndim == 2 else states[:,:,i]
            
            # Create the MSC (using the appropriate *from_coord_arrays()* function )
            M = MSC()
            M.from_coord_arrays(unpacked , self.TDB_init , self.TDB_final , times_TDB , state_slice )
            self.MSCs.append( M )

        return self.MSCs



    def _populate_from_database(self, args):
        '''
            Will need method to construct MSC from 1-or-many sectors stored in the sqlite db as coeffs
        '''
        pass










class MSC():
    '''
            Multi-Sector Cheby Class
            
            Will hold chebyshev coefficients for a **SINGLE** object
            
             
    '''
    def __init__(self, **kwargs):
    

        # Initialization of default standard PARAMETERS / OPTIONS we use for chebyshevs, etc
        self.minorder       = 17                                    # : Fitting Chebys
        self.maxorder       = 25                                    # : Fitting Chebys
        self.maxerr         = 1e-8                                  # : Fitting Chebys
        
        
        # Time-related quantities
        
        # Fundamental identifiying data for MSC
        self.unpacked_provisional_designation   = None
        self.sector_coeffs                      = {}        # the all-important cheby-coeffs

    def from_coord_arrays(self, unpacked_provisional_designation, TDB_init , TDB_final , times_TDB , states ):
        '''
           Populate the MSC starting from supplied numpy-arrays
            
        '''
        
        # Store the important quantities
        self.unpacked_provisional_designation, self.TDB_init , self.TDB_final = unpacked_provisional_designation, TDB_init , TDB_final

        # Sanity-checks on the supplied start & end times (generally the Loader should have got things right) ...
        assert self.TDB_init >= standard_MJDmin and self.TDB_final <= standard_MJDmax, ' TDB_init & TDB_final not in allowed range'
        assert self.TDB_init < self.TDB_final,   \
            ' Problem with TDB_init [%r] > TDB_final [%r] : likely due to supplied JDs falling outside standard ranges' % \
                (self.TDB_init , self.TDB_final)

        # Sanity-check on dimensionality of input states
        assert states.ndim == 2, 'states.ndim = %d' % states.ndim


        # Here we compare the endpoints to sector-endpoints to ensure sectors are fully supported by supplied data
        sector_numbers, sector_JD0s = self.map_JD_to_sector_number_and_sector_start_JD( [self.TDB_init,self.TDB_final] )
        init           = (sector_numbers[0], sector_JD0s[0] ) if self.TDB_init - sector_JD0s[0]   < 1                       else (sector_numbers[0] + 1 , sector_JD0s[0] + sector_length_days)
        final          = (sector_numbers[-1],sector_JD0s[-1]) if sector_JD0s[-1] + sector_length_days - self.TDB_final  < 1 else (sector_numbers[-1] - 1, sector_JD0s[-1]- sector_length_days)

        # Sanity check
        assert init[0] <= final[0] , \
            ' Problem with sector numbers [%r, %r] : likely due to supplied JDs falling outside standard ranges' % \
                (init[0] , final[0])


        # Go through sector-numbers
        for sector_num in range(init[0] , final[0] ):
            
            # Identify the indicees of the nbody times for this sector (i.e. those with min < t < max)
            # N.B. a[:,0] == times
            sector_TDB_init    = standard_MJDmin + sector_num     * sector_length_days
            sector_TDB_final   = standard_MJDmin + (sector_num+1) * sector_length_days
            indicees           = np.where((times_TDB >=sector_TDB_init )   & \
                                          (times_TDB <=sector_TDB_final)    )[0]

            # Order used for cheby fitting
            self.maxorder   = min(self.maxorder,len(indicees))
            
            # To try and improve accuracy, make times relative to standard_MJDmin
            relative_times = times_TDB - standard_MJDmin
            
            # Calc the coeffs: Do all coordinates & covariances simultaneously
            # N.B. (1) For a single particle, states.shape = (33, 27)
            # - Where 33 = Number times, and 27=Number of components (x,y,z,u,v,w,...)
            # N.B. (2) cheb_coeffs.shape ~ (18, 27)
            # - Where 18 = Number of coeffs and 27=Number of components (x,y,z,u,v,w,...)
            cheb_coeffs = self.generate_cheb_for_sector( relative_times[indicees], states[indicees] )
        
            # Save the fitted coefficients into the sector_coeff dict
            self.sector_coeffs[sector_num] =  cheb_coeffs




    # Function(s) to fit supplied chebys to coords/data
    # --------------------------------------------------------------
    def generate_cheb_for_sector(self, t, y):
        
        '''
            Get lowest order sufficiently accurate Chebyshev polynomial fit
            for a single sector of data
            
            Note recursion
            
            Inputs:
            t: np.array
             - indep variable, expected to be times for this application
             - to try and improve accuracy, assume times relative to standard_MJDmin
             
            y: np.array
             - dep. var. for which we fit cheby-polynomials
            
            returns:
            --------
            np.array
            
        '''
        order           = self.minorder
        chebCandidate   = np.polynomial.chebyshev.chebfit(t, y, int(np.ceil(order)) )
        quickEval       = np.polynomial.chebyshev.chebval(t, chebCandidate).T
        if np.max( np.abs(quickEval - y) ) <= self.maxerr or int(np.ceil(order)) == self.maxorder :
            return chebCandidate
        else:
            return self.generate_cheb_for_sector(t,y)







    # Assorted Utility Functions ...
    # --------------------------------------------------------------

    def map_JD_to_sector_number(self, JD_TDB , JD0=standard_MJDmin):
        return ( np.asarray( JD_TDB ) - JD0).astype(int) // sector_length_days
    
    def map_sector_number_to_sector_start_JD(self, sector_number, JD0=standard_MJDmin):
        return JD0 + sector_number * sector_length_days
    
    def map_JD_to_sector_number_and_sector_start_JD(self, JD_TDB, JD0=standard_MJDmin):
        '''
            For a given JD_TDB, calculate the sector number it will be in
            Assumes a standard starting julian date
            
            inputs:
            -------
            
            return:
            -------
            '''
        sector_number   = self.map_JD_to_sector_number(JD_TDB , JD0=JD0)
        return sector_number , self.map_sector_number_to_sector_start_JD( sector_number, JD0=JD0)



    def get_valid_range_of_dates( self,  ):
        '''
            Return the minimum and maximum dates for which this MSC is valid
        '''
        return self.TDB_init , self.TDB_final

    """
    def map_times_to_sectors( self, times_tdb ):
        '''
            Given query-times, it is likely to be useful to
            map each time to the relevant single-sector-dictionary
            
            inputs:
            -------
            times_tdb : np.array
            - JD TDB of times at which positions are to be calculated
            
            return:
            -------
            np.array of integers
            - sector # (zero-based) starting from the dictionary's "t_init"
            - length of returned array = len(times)
            '''
        if RELATIVE:
            return ( (times_tdb ) // sector_length_days ).astype(int)
        else:
            return ( (times_tdb - self.TDB_init ) // sector_length_days ).astype(int)
    """
    


    # Functions to evaluate supplied multi-sector-cheby-dictionary
    # --------------------------------------------------------------

    def generate_HP( self, times_tdb , observatoryXYZ , APPROX = False, CHECK = False ):
        '''
            Calculate apparent HP-locn from specified observatory-posn(s) at given time(s)
            N.B. observatory-posn(s) must be externally calculated/supplied
            
            inputs:
            -------
            times_tdb : np.array
            - JD TDB of times at which positions are to be calculated
            
            observatoryXYZ: np.array
            - Observatory positions at times.mjd [utc]
            - Dimension =3*len(times)
            
            APPROX: boolean
            - Allow approximate calc ( *no* LTT-correction) of unit-vector
            
            CHECK: boolean
            - Allow validity-checking to be turned on/off
            
            return:
            -------
            np.array of integer healpix
            - length of returned array = len(times)
            
            '''
        
        
        # Get the unit vector from observatory to object
        # NB UV.shape =  (3, len(times_tdb) )
        UV = self.generate_UnitVector(  times_tdb ,
                                        observatoryXYZ,
                                        APPROX = APPROX )

        # Calc the HP from the UV and return
        return healpy.vec2pix(HP_nside, UV[0], UV[1], UV[2], nest=True if HP_order=='nested' else False )


    def generate_UnitVector( self, times_tdb , observatoryXYZ, APPROX = False , DELTASWITCH = False, delta=np.array([0,0,0]) ):
        '''
            Calculate apparent UnitVector from specified observatory-posn(s) at given time(s)
            
            N.B. observatory-posn(s) must be :
            (i) externally calculated/supplied
            (ii) equatorial frame (to match assumed frame of orbit)
            
            inputs:
            -------
            times_tdb : np.array
             - JD TDB of times at which positions are to be calculated
            
            observatoryXYZ: np.array
            - Observatory positions at times_tdb
            - shape = (3,len(times_tdb))
            
            return:
            -------
            unit-vectors: np.array 
             - apparent UnitVector from specified observatory-posn(s) at given time(s)
            
        '''
        
        # Get the LTT-corrected position
        # - We allow for the possibility of *NOT* iterating (i.e. doing an approx calc.)
        n_iterations    = 1 if APPROX else 3
        
        # Init light-delay to 0
        lightDelay      = np.zeros( len(times_tdb) )
        
        for i in range(n_iterations):

            # Calculate delayed time (of emission)
            # N.B. delayedTimes.shape = (len(times_tdb) )
            delayedTimes    = times_tdb - lightDelay
            
            # Extract posn of objects at each delayed-time
            # N.B. objectXYZ.shape    = (3, len(times_tdb) )
            objectXYZ       = self.generate_XYZ( delayedTimes )
            if DELTASWITCH :
                objectXYZ += np.stack([delta for i in range(objectXYZ.shape[-1])], axis=1)
            
            # Calculate relative sepn-vector from observatory-to-object
            sepn_vectors    = objectXYZ - observatoryXYZ
            
            # Calculate distance to object at each time.
            # N.B. sepn_vectors.shape == (3, len(times_tdb)), d.shape=(Np,len(times_tdb))
            d               = np.linalg.norm(sepn_vectors, axis=0)
            
            # Calculate light-travel-time
            lightDelay      = d / (astropy.constants.c * secsPerDay / astropy.constants.au ).value
    
        # Return unit-vector
        return sepn_vectors / d

    def dUVdXYZ( self, times_tdb , observatoryXYZ ):
        '''
            Gradient of the UnitVector w.r.t. the Cartesian X,Y,Z positions
            
            inputs:
            -------
            
            returns:
            --------
        '''
        d = 1e-6
        _dX = ( self.generate_UnitVector( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH = True, delta=np.array([d, 0, 0]) ) \
               -self.generate_UnitVector( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH = True, delta=np.array([-d,0, 0]) ) )
        _dY = ( self.generate_UnitVector( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH = True, delta=np.array([0, d, 0]) ) \
               -self.generate_UnitVector( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH = True, delta=np.array([0,-d, 0]) ) )
        _dZ = ( self.generate_UnitVector( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH = True, delta=np.array([0, 0, d]) ) \
               -self.generate_UnitVector( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH = True, delta=np.array([0, 0,-0]) ) )
    
        return np.stack(np.array( (_dX, _dY, _dZ) ), axis=1).T / (2*d)
    
    def covUV(self, times_tdb , observatoryXYZ ):
        '''
            Evaluate the covariance in unit vectors
            This is calculated using the covariance in XYZ & the gradient of the UV components w.r.t. XYZ
            '''
        dUV     = self.dUVdXYZ( times_tdb , observatoryXYZ )
        cov_XYZ = self.covXYZ( times_tdb  )
        return np.array( [ np.linalg.multi_dot([dUV[i].T , cov_XYZ[i], dUV[i]]) for i in range(len(dUV)) ] )

    def generate_RaDec(self,  times_tdb  , observatoryXYZ=None, APPROX = False , DELTASWITCH = False, delta=np.array([0,0,0])):
        '''
            Calculate apparent RA,DEC (Radians???) from specified observatory-posn(s) at given time(s)
            
            N.B. observatory-posn(s) must be :
            (i) externally calculated/supplied
            (ii) equatorial frame (to match assumed frame of orbit)
            
            inputs:
            -------
            times_tdb : np.array
            - JD TDB of times at which positions are to be calculated
            
            observatoryXYZ: np.array
            - Observatory positions at times_tdb
            - shape = (3,len(times_tdb))
            
            return:
            -------
            RA_Dec: np.array
            - *** degrees ***
            
        '''
        # Get the unit vector from observatory to object
        UV = self.generate_UnitVector(times_tdb ,
                                      observatoryXYZ,
                                      APPROX = APPROX,
                                      DELTASWITCH=DELTASWITCH,
                                      delta = delta)
            
        # Convert from unit-vector to RA, Dec and then return
        return np.array(healpy.vec2ang( UV.T , lonlat = True ))

    def dRaDecdXYZ( self, times_tdb , observatoryXYZ ,         d = 1e-6):
        '''
            Gradient of Ra & Dec w.r.t. the Cartesian X,Y,Z positions
            
            Assumed d = 1e-8 => 1e-8 AU => 1e3m
            
            inputs:
            -------
            
            returns:
            --------
            '''
        _dX = ( self.generate_RaDec( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH = True, delta=np.array([d, 0, 0]) ) \
               -self.generate_RaDec( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH = True, delta=np.array([-d,0, 0]) ) )
        _dY = ( self.generate_RaDec( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH = True, delta=np.array([0, d, 0]) ) \
               -self.generate_RaDec( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH = True, delta=np.array([0,-d, 0]) ) )
        _dZ = ( self.generate_RaDec( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH = True, delta=np.array([0, 0, d]) ) \
               -self.generate_RaDec( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH = True, delta=np.array([0, 0,-0]) ) )
    
        return np.stack(np.array( (_dX, _dY, _dZ) ), axis=1).T / (2*d)
    
    def covRaDec(self, times_tdb , observatoryXYZ ):
        '''
            Evaluate the covariance in RA, Dec 
            This is calculated using the covariance in XYZ & the gradient of RA,Dec w.r.t. XYZ
        '''
        dRD     = self.dRaDecdXYZ( times_tdb , observatoryXYZ )
        cov_XYZ = self.covXYZ( times_tdb  )
        return np.array( [ np.linalg.multi_dot([dRD[i].T , cov_XYZ[i], dRD[i]]) for i in range(len(dRD)) ] )

    

    def generate_XYZ( self, times_tdb  ):
        '''
            Evaluate the XYZ positions at the supplied times
            Convenience wrapper around *evaluate_components()* func
            Ensures we only evaluate XYZ components of the coefficients 
            
            inputs:
            -------
            times_tdb : np.array
            - JD TDB of times at which positions are to be calculated
            
            return:
            -------
            XYZ_posns : np.ndarray
             - N.B. XYZ_posns.shape = (3, Np, len(times_tdb) )
               where Np = Number of particles (==1 for single-particle case)
        '''
        
        
        # Just select/evaluate the first 3-sets of components (X,Y,Z)
        # - See link for specifying slices & ellipsis  : https://docs.scipy.org/doc/numpy/user/basics.indexing.html
        if self.sector_coeffs[ list(self.sector_coeffs.keys())[0] ].ndim == 2  :
            slice_spec = (slice(0,len(self.sector_coeffs[ 0 ])), slice(0,3))
        else:
            sys.exit('self.sector_coeffs[ 0 ].ndim = %d : unable  to proceed if ndim != 2 ' % self.sector_coeffs[ 0 ].ndim  )
        
        # Evaluate only the XYZ coefficients
        # NB: to try to increase accuracy, use times relative to standard_MJDmin
        return self.evaluate_components( times_tdb - standard_MJDmin , slice=slice_spec )

    def covXYZ( self, times_tdb  ):
        '''
            Evaluate the covariance in XYZ positions at the supplied times
            Convenience wrapper around *evaluate_components()* func
            Follows the approach in *generate_XYZ()* used to select required coefficients
            
            returns 
            -------
            np.ndarray
             - shape = ( len(times_tdb) ,3, 3 )
        '''
        # Select/evaluate the appropriate covariance components (...)
        # ['x_x', 'x_y', 'x_z', 'x_vx', 'x_vy', 'x_vz', 'y_y', 'y_z', 'y_vx', 'y_vy', 'y_vz', 'z_z', 'z_vx', 'z_vy', 'z_vz', 'vx_vx', 'vx_vy', 'vx_vz', 'vy_vy', 'vy_vz', 'vz_vz']
        #    6      7      8      9      10      11      12     13     14      15      16      17     18      19      20      21       22       23       24       25       26
        #    *      *      *                             *      *                              *
        cov = self.evaluate_components( times_tdb - standard_MJDmin , slice=(slice(0,len(self.sector_coeffs[ 0 ])), np.array([6,7,8,12,13,17]) ) )
        
        # Reproduce & reshape to get 3x3 matrix for each time
        # 'x_x', 'x_y', 'x_z',  'y_y', 'y_z',  'z_z' -->> 'x_x', 'x_y', 'x_z',  'x_y', 'y_y', 'y_z',  'x_z', 'y_z','z_z'
        #  0      1       2      3      4       5          0      1      2       1      3      4       2      4     5
        #                                                  0      1      2       3      4      5       6      7     8
        #
        # Transpose to make it come out as shape = (Nt,3,3), so that its comparable to dRaDecdXYZ & dUVdXYZ
        return np.array( [cov[0:3,:], np.vstack( (cov[1], cov[3:5,:]) ), np.vstack( (cov[2], cov[4:6,:]) ) ] ).T


    def dXYZdt(self,  times_tdb  , dt=1e-5):
        ''' 
            Use calculate the gradient in XYZ at supplied times
            
            *** DO I EVER NEED dX/dt ??? ***
            
            returns:
            --------

            
        '''
        # numdifftools ~10x slower than direct method below ...
        # return nd.Gradient( self.generate_XYZ )(times_tdb)
        #
        # dt~1e-5 => ~1sec
        return ( self.generate_XYZ( times_tdb + dt ) - self.generate_XYZ( times_tdb - dt ) ) / (2*dt)
    

    def evaluate_components( self, relative_times , slice=() ):
        '''
            
            inputs:
            -------
            times_tdb : np.array
            - to try to increase accuracy, use times relative to standard_MJDmin
            
            return:
            -------
            components calculated from coefficient evaluation at times_tdb
        '''
        # Make sure that the supplied times are *RELATIVE*
        # - Going to do this by demanding that the times be larger than ~0
        # - But because of the LTT calculations reqd in *generate_UnitVector()*, I'm going to allow extra negative time ...
        #   1000AU => ~6days LTT, which should be more than plenty.
        assert relative_times[0] >= -6 and relative_times[-1] < (standard_MJDmax - standard_MJDmin), ' times do not seem to be relative: relative_times[0] = %r' % (relative_times[0] )
        
        # Find which single-sector dictionary to use for each given time
        # Then make a dictionary with key=sector-number, and value=list-of-times
        # N.B. we use JD0 = 0 because we have already been supplied with *relative* times
        times_for_each_sector= defaultdict(list)
        for t , s in zip( relative_times , self.map_JD_to_sector_number(relative_times, JD0=0) ):
            times_for_each_sector[s].append(t)
    
        # Evaluate the chebyshev polynomial
        # - Note we do all of the evaluations for a single sector in one go
        # - So we only need to loop over the sectors
        for n, s in enumerate(sorted(list(set(times_for_each_sector)))):
            XYZs = np.polynomial.chebyshev.chebval( times_for_each_sector[ s ] , self.sector_coeffs[ s ][slice] )
            if n == 0 :
                combinedXYZs = XYZs
            else :
                combinedXYZs = np.append( combinedXYZs, XYZs, axis = -1 )

        return combinedXYZs





# End Of File
# --------------------------------------------------------------


"""
    def generate_multi_sector_cheby_dict_from_nbody_json( json_filepath ,
    TDB_init,
    TDB_final,
    minorder=17,
    maxorder=25,
    maxerr=1e-8,
    FORCE_DATES = False ,
    CHECK = False):
    '''
    Slightly rewritten version of Margaret's *datediv* routine
    
    inputs:
    -------
    json_filepath   -- filepath of json file to be read
    TDB_init,TDB_final -- earliest and latest MJD dates for which you want chebyshev approximations
    
    returns:
    --------
    multi-sector cheby-dict
    '''
    
    # Read the nbody-json file
    nbody_dict = nbody_reader.parse_nbody_json( json_filepath )
    
    # Check whether the supplied data can support the requested date-range
    if FORCE_DATES :
    TDB_init = max(TDB_init , nbody_dict[ nbody_reader.time_fieldname ][0] )
    TDB_final = min(TDB_final , nbody_dict[ nbody_reader.time_fieldname ][-1] )
    else:
    if  TDB_init < nbody_dict[ nbody_reader.time_fieldname ][0] or \
    TDB_final > nbody_dict[ nbody_reader.time_fieldname ][-1]:
    sys.exit(' nbody data does not support the requested dates ')
    
    # Set up a (mostly-empty) multi-sector cheby-dict
    mscd = {
    nbody_reader.object_name    : nbody_dict[nbody_reader.object_name],
    'MJP_TDB_init'      : TDB_init,
    'MJP_TDB_final'     : TDB_final,
    'sectors'           : []
    }
    
    # Split into sectors ...
    numdivs = int(np.ceil((TDB_final-TDB_init)/sector_length_days))
    
    # Go through sectors
    for ind in range(numdivs):
    
    # Set up a (mostly-empty) single-sector cheby-dict
    sscd = { nbody_reader.object_name : nbody_dict[nbody_reader.object_name] }
    
    # Identify the indicees of the nbody times for this sector (i.e. those with min < t < max)
    sscd['MJP_TDB_init']    = TDB_init + ind*sector_length_days
    sscd['MJP_TDB_final']   = min(TDB_final,TDB_init + (ind+1)*sector_length_days)
    indicees                = np.where((nbody_dict[ nbody_reader.time_fieldname ]>=sscd['MJP_TDB_init']) & \
    (nbody_dict[ nbody_reader.time_fieldname ]<=sscd['MJP_TDB_final']) )[0]
    
    
    # Loop over all coordinates & covariances
    lists = [nbody_reader.coord_names, nbody_reader.covar_names]
    for item in itertools.chain(*lists):
    
    # For each component, generate chebys & save into single-sector cheby-dict
    maxorder    = min(maxorder,len(indicees))
    sscd[item]  = generate_single_sector_cheb(  nbody_dict[ nbody_reader.time_fieldname ][indicees],
    nbody_dict[item][indicees],
    minorder,
    maxorder,
    maxerr)
    # append the single-sector cheby-dict into the multi-sector cheby-dict
    mscd['sectors'].append(sscd)
    
    # If being thorough, check that the produced object is valid
    if CHECK:
    check_multi_sector_validity( mscd )
    
    return mscd
    
    
    


    # Functions to check validity ...
    # --------------------------------------------------------------
    def check_validity( cheby_dict ):
        '''
            Check whether the input dictionary has the expected
            structure / variables of a MULTI-SECTOR dictionary

        '''
        # Expected keys & data-types
        expected_keys_and_types = [
        ("name", str),
        ("MJP_TDB_init", (int, float, np.int64, np.float64)),
        ("MJP_TDB_final", (int, float, np.int64, np.float64)),
        ("sectors", (list, np.ndarray) )
        ]

        # Check data is as expected
        # Needs to
        # (i) be a dict
        # (ii) have all the necessary keys
        # (iii) have all the correct data types
        # (iv) individual-sector dictionaries are all valid
        return True if isinstance(cheby_dict , dict) and \
                        np.all([key in cheby_dict and
                        isinstance(cheby_dict[key], typ) for key, typ in
                        expected_keys_and_types]) and \
                        np.all([check_single_sector_validity(sector_dict) for
                           sector_dict in cheby_dict[sectors]]) \
                    else False




    """
