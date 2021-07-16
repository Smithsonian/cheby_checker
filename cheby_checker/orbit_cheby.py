# -*- coding: utf-8 -*-
# cheby_checker/cheby_checker/orbit_cheby.py

'''
    --------------------------------------------------------------
    cheby_checker's orbit_cheby module.
    
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


# Import neighboring packages
# --------------------------------------------------------------
try:
    import nbody_reader
    import nbody
except ImportError:
    from . import nbody_reader
    from . import nbody




class MSC_Loader(Base):
    '''
        Multi-Sector Cheby -Loader Function
        
        Will load/create/return instances of MSC : Multi-Sector Cheby objects
        
        Will handle multi-particle instantiation: will return **LIST** of *MSC* objects
        
    '''

    # allow different init depending on source ...
    def __init__(self, **kwargs):
        print('INIT MSC_Loader...')

        # Give access to "Base" methods & attributes
        Base.__init__(self)

        # Initialization of default standard PARAMETERS / OPTIONS we use for chebyshevs, etc
        # - These may be overwritten by kwargs
        
        self.filepath       = None                                  # : Ingest method (From textfile)

        self.primary_unpacked_provisional_designations   = None     # : Ingest method (From np.array)
        self.times_TDB      = None                                  # : Ingest method (From np.array)
        self.statearray     = None                                  # : Ingest method (From np.array)
        
        self.NbodySim       = None                                  # : Ingest method (From nbody.NbodySim)

        self.FROM_DB        = None                                  # : Ingest method (From sqlite db)
        self.dict_of_dicts  = None                                  # : Ingest method (From sqlite db)
        
        # The list of MSCs that will be instantiated & returned
        self.MSCs           = []
        
        # Do a loop over the parameters and see if any attributes need updating (because they were supplied)
        # I.e. This looks through kwargs for (e.g.) 'filepath', 'NbodySim', 'FROM_DB', etc etc 
        for attribute, value in kwargs.items():
            if hasattr(self, attribute):
                setattr(self, attribute, value)

        # Allow initialization from different sources ...
        #
        # (i) From database (of pre-calculated chebyshevs)
        if self.FROM_DB and self.dict_of_dicts:
            self._populate_from_database( self.dict_of_dicts )

        # (ii) From nbody.NbodySim
        elif self.NbodySim is not None :
            if "primary_unpacked_provisional_designations" not in self.NbodySim.__dict__:
                self.NbodySim.primary_unpacked_provisional_designations = [ str(_) for _ in range(self.NbodySim.output_n_particles)]
                print(f'Populating from NbodySim : ***NO*** designation information : Using {self.NbodySim.primary_unpacked_provisional_designations}')
            self._populate_from_nbody_array(self.NbodySim.primary_unpacked_provisional_designations , self.NbodySim.output_times, self.NbodySim.output_vectors)

        # (iii) From numpy array (from nbody: Mainly for development?)
        elif    self.primary_unpacked_provisional_designations is not None and \
                isinstance(self.statearray , (list,tuple,np.ndarray)) and \
                isinstance(self.times_TDB , (list,tuple,np.ndarray)) and \
                self.statearray.shape[0] == self.times_TDB.shape[0]:
            self._populate_from_nbody_array(self.primary_unpacked_provisional_designations, self.times_TDB, self.statearray)

        # (iv) From text file ( Mainly for development!)
        elif self.filepath != None :
            self._populate_from_nbody_text(self.filepath)

        # (v) An empty instance ...
        else:
            self._generate_empty()
                


    def _generate_empty(self,  ):
        '''  '''
        print('\n','*'*6,'Defaulting to the production of a list of empty MSCs','*'*6,'\n')
        print('\t','This occurs (a) on erroneous input, and (b) when no input supplied\n')
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
                                   primary_unpacked_provisional_designations ,
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
        



        # Check that the dimensionality of the coordinates is consistent with the number of names
        # - N.B. We expect states.shape = (Nt, Nc, Np) or (Nt,Nc), where
        #        Nt = Number of times
        #        Nc = Number of coords/components being fitted
        #        Np = Number of particles
        #
        self.primary_unpacked_provisional_designations = np.atleast_1d(primary_unpacked_provisional_designations)
        #   (i) Single Object
        if len(self.primary_unpacked_provisional_designations) == 1 and states.ndim == 2 :
            pass
        #  (ii) Multiple Objects (Nt, Nc, Np) [2 flavors to try and cope with the inevitable fuck-ups in coordinate organization]:
        elif len(self.primary_unpacked_provisional_designations) >= 1 and states.ndim == 3 and (len(self.primary_unpacked_provisional_designations) == states.shape[2]) :
            pass
        # (iii) Multiple Objects (Nt, Np, Nc) [2 flavors to try and cope with the inevitable fuck-ups in coordinate organization]
        #  **** NEED TO GET RID OF (ii) or (iii) : IT WILL GO WRONG EVENTUALLY IF/WHEN Nc==Np ...
        elif len(self.primary_unpacked_provisional_designations) >= 1 and states.ndim == 3 and (len(self.primary_unpacked_provisional_designations) == states.shape[1]) :
            pass
        #  (iv) Problem
        else:
            sys.exit('Inconsistent dimensionality : len(self.unpacked_provisional_designations) = %r and states.ndim = %r and states.shape = %r' % \
                     (len(self.primary_unpacked_provisional_designations) , states.ndim ,  states.shape) )






        # Loop over each of the objects and create an MSC-object for each ...
        for i, unpacked in enumerate(self.primary_unpacked_provisional_designations):
            
            # Get the slice of states corresponding to the particular named object
            #  **** NEED TO GET RID OF (ii) or (iii) : IT WILL GO WRONG EVENTUALLY IF/WHEN Nc==Np ...
            state_slice = states if states.ndim == 2 else states[:,:,i] if len(self.primary_unpacked_provisional_designations) == states.shape[2] else states[:,i,:]
            
            # Create the MSC (using the appropriate *from_coord_arrays()* function )
            M = MSC()
            M.from_coord_arrays(unpacked , times_TDB , state_slice )
            self.MSCs.append( M )

        return self.MSCs



    def _populate_from_database(self, dict_of_dicts ):
        '''
            Method to construct MSCs from 1-or-many sectors stored in the sqlite db as coeffs

            inputs:
            -------
            dict_of_dicts : dictionary-of-dictionaries
             - As returned by query_coefficients_by_jd_hp
             - Outer dict is key-ed on primary_unpacked_provisional_designation
             - Inner dicts are key-ed on sector

            returns:
            --------
            list of MSC objects

        '''
        
        # Loop over objects and make list of MSCs
        MSCs = []
        for primary_unpacked_provisional_designation in np.atleast_1d(primary_unpacked_provisional_designations) :
        
            # Create the MSC (using the appropriate *from_database()* function )
            M = MSC()
            M.from_database( primary_unpacked_provisional_designation )
            self.MSCs.append( M )

        return MSCs



### ****
### ****
### ****
### ****
### ****
### ****
### ****
### ****
### ****
### ****
### ****
### ****
### ****
### ****
### ****
### ****




class MSC(Base):
    '''
            Multi-Sector Cheby Class
            
            Will hold chebyshev coefficients for a **SINGLE** object
            
             
    '''
    def __init__(self, **kwargs):
    
        # Give access to "Base" methods & attributes
        Base.__init__(self)

        # Initialization of default standard PARAMETERS / OPTIONS we use for chebyshevs, etc
        self.minorder       = 7                                     # : Fitting Chebys
        self.maxorder       = 25                                    # : Fitting Chebys
        self.maxerr         = 1e-8                                  # : Fitting Chebys
        
        # Fundamental identifiying data for MSC
        self.primary_unpacked_provisional_designation   = None
        self.sector_coeffs                      = {}        # the all-important cheby-coeffs
    
    def from_database(self, primary_unpacked_provisional_designation , sector_numbers = None):
        '''
            Used to initialize (create) an MSC object from data in the sqlite database 
            
            (as extracted using the *get_nightly_precalcs()* function in precalc.py )
            
            inputs:
            -------
            primary_unpacked_provisional_designation : string
            
            returns:
            --------
            True
             - Doesn't directly return, just populates the MSC object

        '''
        # unpacked primary provID of the object (e.g. 2020 AA)
        self.primary_unpacked_provisional_designation = primary_unpacked_provisional_designation
        
        # query the database for the required sector coefficients
        coeff_dict = precalc.get_specific_object(primary_unpacked_provisional_designation ,
                                                 sector_numbers = sector_numbers)
    
        # create the entries in the sector_coeffs dict (N.B. have to transform the keys)
        self.sector_coeffs = { int(k.split("_")[1]) : v for k,v in coeff_dict.items() }
    
        # set the min & max supported times / sectors
        self.sector_init, self.sector_final = np.min(list(self.sector_coeffs.keys())) , np.max(list(self.sector_coeffs.keys()))
        self.TDB_init  = Base.map_sector_number_to_sector_start_JD( self.sector_init  , self.standard_MJDmin)
        self.TDB_final = Base.map_sector_number_to_sector_start_JD( self.sector_final , self.standard_MJDmin) + self.sector_length_days - self.epsilon
    
        # check things are as expected
        assert self.TDB_init >= self.standard_MJDmin and self.TDB_final <= self.standard_MJDmax, \
            'Problem with limits in from_database: self.TDB_init = [%f] , self.standard_MJDmin = [%f] self.TDB_final = [%f], self.standard_MJDmax = [%f]' % (self.TDB_init , self.standard_MJDmin , self.TDB_final , self.standard_MJDmax)
    

    def from_coord_arrays(self, primary_unpacked_provisional_designation, times_TDB , states ):
        '''
           Populate the MSC starting from supplied numpy-arrays
           Expected to be used when passing in the data from an NBody integration (REBOUND)
            
            inputs:
            -------
            primary_unpacked_provisional_designation : 
            - 
            
            TDB_init :
            -
            
            TDB_final :
            -
            
            times_TDB :
            -
            
            states :
            -

            returns:
            --------
            True
            - Doesn't directly return, just populates the MSC object
            
        '''
        # Store name
        self.primary_unpacked_provisional_designation = primary_unpacked_provisional_designation
        
        # Sanity-check on dimensionality of input states
        assert states.ndim == 2, 'states.ndim = %d' % states.ndim

        # Generate relative times within sectors (each sector starts from t=0)
        sector_numbers , sector_relative_times = self.map_JDtimes_to_relative_times(times_TDB)
        
        # Go through sector-numbers & populate "sector_coeffs" dict
        sector_numbers = sorted(sector_numbers)
        for sector_number in set(sector_numbers) :
            
            # Identify the indicees of the nbody times for this sector
            indicees = np.where( sector_numbers == sector_number )[0]
            
            # Order used for cheby fitting
            self.maxorder   = min(self.maxorder,len(indicees))
            
            # Calc the coeffs: Do all coordinates & covariances simultaneously
            # N.B. (1) For a single particle, states.shape = (33, 27)
            # - Where 33 = Number times, and 27=Number of components (x,y,z,u,v,w,...)
            # N.B. (2) cheb_coeffs.shape ~ (18, 27)
            # - Where 18 = Number of coeffs and 27=Number of components (x,y,z,u,v,w,...)
            cheb_coeffs, maxFitErr = self.generate_cheb_for_sector( sector_relative_times[indicees], states[indicees] )

            # Check that the maxFitErr is within acceptable bounds
            if maxFitErr > self.maxerr:
                sys.exit(f'MSC.from_coord_arrays: maxFitErr = {maxFitErr} (i.e. > self.maxerr = {self.maxerr})')
            
            # Check that the sector is "well covered"
            # - I.e. N_points >= N_cheb_coeffs  &   t[0] < sector_gap  &  t[-1] > sector_length_days - sector_gap
            # N.B. I have not yet investigated/justified these values
            #      But I have verified that things can go terribly wrong unless *some* constraint is put in place
            # N.B. This check creates the possibility of patchy / non-contiguous sectors
            elif    states[indicees].shape[0] >= cheb_coeffs.shape[0] \
                    and sector_relative_times[indicees][0]  < self.sector_gap \
                    and sector_relative_times[indicees][-1] > self.sector_length_days - self.sector_gap :
            
                # Save the fitted coefficients into the sector_coeff dict
                self.sector_coeffs[sector_number] =  cheb_coeffs
                    
            else:
                if sector_number not in [sector_numbers[0], sector_numbers[-1]]:
                    print(f'Warning : sector_number {sector_number} is not well covered')
                    print(f'This is *not* an  "end" sector so this *will* cause coverage gaps')
                    print(f'states[indicees].shape[0] , cheb_coeffs.shape[0], sector_relative_times[indicees][0], sector_relative_times[indicees][-1]')
                    print(  states[indicees].shape[0] , cheb_coeffs.shape[0], sector_relative_times[indicees][0], sector_relative_times[indicees][-1] )
                    print()
    
        # May be useful to extract start & end sectors / dates
        supported_sector_numbers = sorted(self.sector_coeffs.keys())
        self.sector_init , self.sector_final = supported_sector_numbers[0], supported_sector_numbers[-1]
        self.TDB_init   = self.map_sector_number_to_sector_start_JD( self.sector_init,  self.standard_MJDmin )
        self.TDB_final  = self.map_sector_number_to_sector_start_JD( self.sector_final, self.standard_MJDmin ) + self.sector_length_days - self.epsilon
            

    

    # Function(s) to fit supplied chebys to coords/data
    # --------------------------------------------------------------
    def generate_cheb_for_sector(self, t, y, order=None):
        
        """
            Get lowest order accurate Chebyshev polynomial fit for a single sector
            
            Note recursion
            
            Inputs:
            --------
            t: np.array
             - indep variable, expected to be times for this application
             - t.shape = (Nt,), where Nt is the number of times being evaluated
             - to try and improve accuracy, assume times relative to SECTOR-ZERO TIME
                  ***** I.E. EXPECT THAT 0 < t < Base.sector_length_days *****
             
            y: np.array
             - dep. var. for which we fit cheby-polynomials
             - y.shape = (Nt, Nc), where Nt is defined above and Nc is the number of components being fitted
            
            Returns:
            --------
            np.array
             - shape = (No, Nc), 
             - where Nc is defined above and No is 1+order (and note that "order" can iterate upwards)
            
        """
        # EXPECT THAT 0 < t < Base.sector_length_days (keep times <~32 days)
        assert np.all( t < self.sector_length_days )
        
        order           = self.minorder if order is None else order
        chebCandidate   = np.polynomial.chebyshev.chebfit(t, y, int(np.ceil(order)) )
        quickEval       = np.polynomial.chebyshev.chebval(t, chebCandidate).T
        maxFitErr       = np.max( np.abs(quickEval - y) )
        
        # If fit errors are small, we can finish
        if maxFitErr <= self.maxerr:
            return chebCandidate, maxFitErr
    
        # If we have maxed out the fitting-order, finish but issue warning
        elif int(np.ceil(order)) == self.maxorder :
            print(f'\n\t ***WARNING [orbit_cheby.generate_cheb_for_sector]: int(np.ceil(order)) == self.maxorder == {self.maxorder}')
            print(f'\t t.shape, y.shape : {t.shape, y.shape}')
            print(f'\t t[:3]...t[-3:] : {t[:3]}...{t[-3:]} ')
            return chebCandidate, maxFitErr
        
        # If fit errors are large, try increasing the order of the fit
        else:
            return self.generate_cheb_for_sector(t, y, order + int(np.ceil((self.maxorder - order) / 2)))








    # Date Function(s) ...
    # --------------------------------------------------------------

    def get_valid_range_of_dates( self,  ):
        '''
            Return the minimum and maximum dates for which this MSC is valid
        '''
        return self.TDB_init , self.TDB_final



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
        return healpy.vec2pix(self.HP_nside, UV[0], UV[1], UV[2], nest=True if self.HP_order=='nested' else False )


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
            lightDelay      = d / (astropy.constants.c * self.secsPerDay / astropy.constants.au ).value
    
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
               -self.generate_UnitVector( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH = True, delta=np.array([0, 0,-d]) ) )
    
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
        
        # Is there a way to make this more efficient?
        # - Similar to the nbody code ...
        # - CoV_t = np.linalg.inv( GammaStack_t )
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
        if self.sector_coeffs[ self.sector_init ].ndim == 2  :
            slice_spec = slice(0,3)
        else:
            sys.exit('self.sector_coeffs[ 0 ].ndim = %d : unable  to proceed if ndim != 2 ' % self.sector_coeffs[ 0 ].ndim  )

        # Evaluate only the XYZ coefficients
        return self.evaluate_components( times_tdb  , component_slice_spec=slice_spec )

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
        cov = self.evaluate_components( times_tdb , component_slice_spec=np.array([6,7,8,12,13,17])  )
        
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
    

    def evaluate_components( self, times_TDB , component_slice_spec=slice(None) ):
        """
            
            inputs:
            -------
            times_TDB : np.array
            
            component_slice_spec:
            
            return:
            -------
            components calculated from coefficient evaluation at times_tdb
        """
        
        # Generate relative times within sectors (each sector starts from t=0)
        sector_numbers , sector_relative_times = self.map_JDtimes_to_relative_times(times_TDB)
        
        # Find which single-sector dictionary to use for each given time
        # Then make a dictionary with key=sector-number, and value=list-of-times
        # N.B. we use JD0 = 0 because we have already been supplied with *relative* times
        times_for_each_sector= defaultdict(list)
        for t , s, srt  in zip( times_TDB,  sector_numbers , sector_relative_times ):
            times_for_each_sector[s].append(srt)
        
            # Put in a check for validity here, by asking whether the sectors-calculated here ...
            # ... are actually in the main sector_coeffs-dict
            assert s >= self.sector_init and s<= self.sector_final, \
                'Sector #%d (for time %f) not supported in sector_coeffs (min,max supported times = %f,%f)' % \
                (s, t   , self.TDB_init, self.TDB_final)
        
        # Evaluate the chebyshev polynomial
        # - Note we do all of the evaluations for a single sector in one go
        # - So we only need to loop over the sectors
        for n, s in enumerate(sorted(list(set(times_for_each_sector)))):
            evaluatedComponents = np.polynomial.chebyshev.chebval( times_for_each_sector[ s ] , self.sector_coeffs[ s ][ : , component_slice_spec] )
            if n == 0 :
                combinedComponents = evaluatedComponents
            else :
                combinedComponents = np.append( combinedComponents, evaluatedComponents, axis = -1 )

        return combinedComponents





# End Of File
# --------------------------------------------------------------
