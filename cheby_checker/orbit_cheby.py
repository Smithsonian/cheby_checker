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
from . import nbody
from .cheby_checker import Base


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
                self.NbodySim.primary_unpacked_provisional_designations = [ str(_) for _ in range(self.NbodySim.input_n_particles)]
                print(f'Populating from NbodySim : ***NO*** designation information : Using {self.NbodySim.primary_unpacked_provisional_designations}')
            self._populate_from_nbody_array(self.NbodySim.primary_unpacked_provisional_designations ,
                                            self.NbodySim.output_times,
                                            self.NbodySim.output_states)

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

    """
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
    """

    def _populate_from_nbody_array( self,
                                    primary_unpacked_provisional_designations ,
                                    times_TDB,
                                    states,
                                    covariances = None ):
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
             - This should be ONLY the main fitted vaiables (e.g. 3-posn, 3-vel, + N-non-grav-coeffs)
             - Expect shape == ( len(times_TDB) , len(primary_unpacked_provisional_designations) , 6/7/8/9)
             - In early stages of development, will assume grav-only, hence Nc == 6
            covariances: np.array
             - This should be the covariances corresponding to the states, so dimension should be (Nt, Np, Nc, Nc)
            
            returns:
            --------
            
        '''
        

       
        # Ensure the passed desig variable is array-like
        self.primary_unpacked_provisional_designations = np.atleast_1d(primary_unpacked_provisional_designations)
        


        # Check that the dimensionality of the inputs is consistent
        # - N.B. We expect states.shape = (Nt, Np, Nc), e.g. = (41, 1, 6)
        #        Nt = Number of times
        #        Nc = Number of particles
        #        Np = Number of coords/components being fitted: should be in [6,7,8,9]
        # - N.B. We expect covariances.shape = (Nt, Np, Nc, Nc), e.g. = (41, 1, 6, 6)
        #        Nt = Number of times
        #        Nc = Number of particles
        #        Np = Number of coords/components being fitted: should be in [6,7,8,9]
 
        # Check for consistent times & array dimensions
        assert len(times_TDB) == states.shape[0]
        assert covariances is None or len(self.primary_unpacked_provisional_designations) == covariances.shape[0]

        # Check for consistent designations & array dimensions
        assert len(self.primary_unpacked_provisional_designations) == states.shape[1]
        assert covariances is None or len(self.primary_unpacked_provisional_designations) == covariances.shape[1]

        # Check for consistent number of state components (and covariance if supplied)
        assert states.shape[2] in [6] # [6,7,8,9]  # ASSUME GRAV-ONLY AT THIS STAGE OF DEVELOPMENT ...
        assert covariances is None or covariances.shape[2] == covariances.shape[3] == states.shape[2]



        # Loop over each of the objects and create an MSC-object for each ...
        for i, unpacked in enumerate(self.primary_unpacked_provisional_designations):
                        
            # Create the MSC (using the appropriate *from_coord_arrays()* function )
            # NB: We need to extract the appropriate slice of states corresponding to the particular named object
            M = MSC()
            M.from_coord_arrays(unpacked , times_TDB , states[:,i,:] , covariances = None if covariances is None else covariances[:,i,:,:])
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
        self.minorder       = 5                                     # : Fitting Chebys
        self.maxorder       = 25                                    # : Fitting Chebys
        self.maxerr         = 1e-8                                  # : Fitting Chebys
        
        # It's probably going to be useful to track the number of components
        # (the various covariance calculations & reshapings can get complicated)
        self.n_coordinates = None
        self.covar_bool    = False
        
        # Fundamental identifiying data for MSC
        self.primary_unpacked_provisional_designation   = None
        self.sector_coeffs                      = {}        # the all-important cheby-coeffs



    # Function(s) to populate MSC  from various sources ...
    # --------------------------------------------------------------

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
             
            populates:
            ----------
            self.primary_unpacked_provisional_designation : string
            
            self.sector_coeffs : dict
            
            self.sector_init : int
            
            self.sector_final : int
            
            self.TDB_init : float
            
            self.TDB_final : float
            

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
    

    def from_coord_arrays(self, primary_unpacked_provisional_designation, times_TDB , states , covariances = None ):
        '''
           Populate the MSC starting from supplied numpy-arrays
           Expected to be used when passing in the data from an NBody integration (REBOUNDX)
            
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
            
            states: np.array
             - This should be ONLY the main fitted variables (e.g. 3-posn, 3-vel, + N-non-grav-coeffs)
             - Expect shape == ( Nt ,  Nc ), where Nc in [6,7,8,9]
             - IN early development, assume Nc == 6 (gravity-only)
            covariances: np.array
             - This should be the covariances corresponding to the states, so dimension should be (Nt, Nc, Nc)

            returns:
            --------
            True
            - Doesn't directly return, just populates the MSC object
            
        '''
        # Store name
        self.primary_unpacked_provisional_designation = primary_unpacked_provisional_designation
        
        # Sanity-check on dimensionality of input states
        # NBody.output_states.shape == ( N_times , N_particles , N_coords )
        assert states.ndim == 2, 'states.ndim = %d' % states.ndim
        assert states.shape[0] == len(times_TDB)
        assert states.shape[1] in [6] # [6,7,8,9] #<<-- GRAVITY ONLY TO START WITH
        assert covariances is None or covariances.shape[0] == len(times_TDB)
        assert covariances is None or covariances.shape[1] == covariances.shape[2] == states.shape[1]
        
        # Useful quantities & mappings ...
        self.n_coordinates = states.shape[1]
        self._define_locations()
        
        # Reduce the input square covariance matrix down to triangular form
        triangular = None if covariances is None else self._take_triangular(covariances)
        
        # Generate *RELATIVE* times within sectors (each sector starts from t=0)
        sector_numbers , sector_relative_times = self.map_JDtimes_to_relative_times(times_TDB)
        
        # Go through sector-numbers & populate "sector_coeffs" dict
        sector_numbers = sorted(sector_numbers)
        for sector_number in set(sector_numbers) :
            
            # Identify the indicees of the nbody times for this sector
            indicees = np.where( sector_numbers == sector_number )[0]

            # Order used for cheby fitting
            self.maxorder   = min(self.maxorder,len(indicees))
            
            # Calc the coeffs: We fit all coordinates & covariances simultaneously
            # N.B. (1) For a single particle, with gravity only...
            #                               states.shape      = (Nt, 6)
            #                               covariances.shape = (Nt, 36)
            #                               triangular.shape  = (Nt, 21)
            #                               hstack.shape      = (Nt, 6+21=27)
            # N.B. (2) When we are able to handle non-gravs,
            #          6 -> 6/7/8/9,     36 -> 36/49/64/81 ,     21->21/36/45
            array_to_be_fit = states[indicees] if covariances is None else np.hstack( (states[indicees],triangular[indicees]) )
            cheb_coeffs, maxFitErr = self.generate_cheb_for_sector( sector_relative_times[indicees], array_to_be_fit )
            
            # Check that the maxFitErr is within acceptable bounds
            if maxFitErr > self.maxerr:
                sys.exit(f'MSC.from_coord_arrays: maxFitErr = {maxFitErr} (i.e. > self.maxerr = {self.maxerr})')
            
            # Check that the sector is "well covered"
            # - I.e. N_points >= N_cheb_coeffs  &   t[0] < sector_gap  &  t[-1] > sector_length_days - sector_gap
            # N.B. I have not yet investigated/justified these values
            #      But I have verified that things can go terribly wrong unless *some* constraint is put in place
            # N.B. This check creates the possibility of patchy / non-contiguous sectors
            elif    states[indicees].shape[0]               >= cheb_coeffs.shape[0] \
                    and sector_relative_times[indicees][0]  < self.sector_gap \
                    and sector_relative_times[indicees][-1] > self.sector_length_days - self.sector_gap :
            
                # Save the fitted coefficients into the sector_coeff dict
                self.sector_coeffs[sector_number] =  cheb_coeffs
                    
            else:
                if sector_number not in [sector_numbers[0], sector_numbers[-1]]:
                    print(f'*** Warning : sector_number {sector_number} is not well covered')
                    print(f'*** This is *not* an  "end" sector so this *will* cause coverage gaps')
                    print(f'*** I.e. NO COEFFICIENTS have been produced for sector_number {sector_number}')
                    print(f'*** states[indicees].shape[0]          = {states[indicees].shape[0]}')
                    print(f'*** cheb_coeffs.shape[0]               = {cheb_coeffs.shape[0]}')
                    print(f'*** sector_relative_times[indicees] = {sector_relative_times[indicees]}')
                    print(f'*** times_TDB[indicees][0] = {times_TDB[indicees][0]}')
                    print(f'*** times_TDB[indicees][-1]= {times_TDB[indicees][-1]} ' )

                    print()
    
        # May be useful to extract start & end sectors / dates
        supported_sector_numbers = sorted(self.sector_coeffs.keys())
        self.sector_init , self.sector_final = supported_sector_numbers[0], supported_sector_numbers[-1]
        self.TDB_init   = self.map_sector_number_to_sector_start_JD( self.sector_init,  self.standard_MJDmin )
        self.TDB_final  = self.map_sector_number_to_sector_start_JD( self.sector_final, self.standard_MJDmin ) + self.sector_length_days - self.epsilon
        
        #print( 'supported_sector_numbers = ', supported_sector_numbers)
        #print( 'self.sector_init, self.sector_final  = ', self.sector_init, self.sector_final )
        #print( 'self.TDB_init, self.TDB_final  = ', self.TDB_init, self.TDB_final )

        
            

    # Function(s) related to CoVariance Matrix Structure ...
    # --------------------------------------------------------------
    #  - Convert back-and-forth betweeen "square" and "triangular" representations of covariance matrices
    #  - Keep track of where the important components are within the representations
    '''
        >>> coord_names     = ['x','y','z','vx','vy','vz']
        0   1   2   3    4    5
        >>> covar_names     = [ "_".join([coord_names[i], coord_names[j]]) for i in range(len(coord_names)) for j in range(i,len(coord_names))  ]
        >>> covar_names
        ['x_x', 'x_y', 'x_z', 'x_vx', 'x_vy', 'x_vz', 'y_y', 'y_z', 'y_vx', 'y_vy', 'y_vz', 'z_z', 'z_vx', 'z_vy', 'z_vz', 'vx_vx', 'vx_vy', 'vx_vz', 'vy_vy', 'vy_vz', 'vz_vz']
           6      7      8      9      10      11      12     13     14      15      16      17     18      19      20      21       22       23       24       25       26
        *      *      *                             *      *                              *
    '''
    def _define_locations(self,):
        '''
        '''
        assert self.n_coordinates >=6 and self.n_coordinates <=9
        
        self.triangular_mapping = {21:6, 28:7, 36:8, 45:9}
        
        # Create name:location map for states ...
        self.coord_map = { 'x':0,'y':1,'z':2,'vx':3,'vy':4,'vz':5}
        for n in range(6, self.n_coordinates): self.coord_map[ ['ng1','ng2','ng3'][n-6] ] = n
        assert self.n_coordinates == len(self.coord_map)
        
        self.inv_map = { v:k for k,v in self.coord_map.items()}
                
        # name:location map covariance reduced to triangular form
        self.tri_map = { "_".join([ self.inv_map[i], self.inv_map[j]]):True for i in range(len(self.coord_map))
                                                                             for j in range(i,len(self.coord_map))}
        self.tri_map = {k:n+self.n_coordinates for n,k in enumerate(self.tri_map)}
        
        # name:location map when we have combined states with triangular covariance
        self.combi_map = {}
        self.combi_map.update(self.coord_map)
        self.combi_map.update(self.tri_map)

        # Convenient to know location of components for covXYZ
        # E.g. When n_coordinates=6 , component_slice_spec=np.array([6,7,8,12,13,17])
        # - See link for specifying slices & ellipsis  : https://docs.scipy.org/doc/numpy/user/basics.indexing.html
        self.XYZ_slice_spec       = np.arange( 0,3 )   # 0,1,2
        self.XYZUVW_slice_spec   = np.arange( 0,6 )    # 0 1 2 3 4 5
        self.covXYZ_slice_spec    = np.array( [self.combi_map[k] for k in ['x_x', 'x_y', 'x_z','y_y', 'y_z', 'z_z' ] ] )
        self.covXYZUVW_slice_spec = np.array( [self.combi_map[k] for k in ['x_x', 'x_y', 'x_z', 'x_vx', 'x_vy', 'x_vz', 'y_y', 'y_z', 'y_vx', 'y_vy', 'y_vz', 'z_z', 'z_vx', 'z_vy', 'z_vz', 'vx_vx', 'vx_vy', 'vx_vz', 'vy_vy', 'vy_vz', 'vz_vz'] ] )


    def _take_triangular(self, covariances ):
        '''
        Take a stack of square matricees and reduce to just the triangular components
        Does the opposite of _make_square

        input:
        -------
            covariances: np.ndarray
            - assume shape == (Nt, Nd, Nd)
            - where Nt = Number of times evaluated & Nd = Number of covariance dimensions
            - We expect Nd == 6 at first, and later 6/7/8/9 when we allow non-gravs
            
        returns:
        --------
            triangular: np.ndarray
            -Returned object is of shape (Nt, 21) [assuming input shape == (Nt,6,6)
                6->21 , 7->28 , 8->36 , 9->45
        '''
        
        # check the shape is correct
        assert covariances.ndim == 3
        assert covariances.shape[1] == covariances.shape[2] == self.n_coordinates
        
        # get upper-triangular elemenet at each time slice
        # NB1: There may be a better (more numpy-like) way of doing this iteration over the firt dimension
        # NB2: Returned object is of shape (Nt, 21) [assuming input shape == (Nt,6,6)
        #      6->21 , 7->28 , 8->36 , 9->45
        i = np.triu_indices(self.n_coordinates)
        return np.array( [_[i] for _ in covariances ] )
        
        
    def _make_square( self, covariances_tri ):
        '''
        Take a stack of triangular components and reconstruct a stack of square matricees
        Does the opposite of _take_triangular
        '''
        
        # check the shape is correct
        assert covariances_tri.ndim == 2
        assert covariances_tri.shape[1] in self.triangular_mapping
        
        # Set up indicees
        sq_dim = self.triangular_mapping[covariances_tri.shape[1]]
        i, j = np.triu_indices( sq_dim )
        print(sq_dim)
        print(i)
        print(j)
        # Create empty array
        M = np.empty((covariances_tri.shape[0] , sq_dim, sq_dim))
        
        # Populate array
        for n,_ in enumerate(covariances_tri):
            print('n=',n)
            M[n,i,j] = _
            M[n,j,i] = _
            print(M[n,:,:])
        # Return array of dimension (covariances_tri.shape[0] , sfq_dim, sq_dim)
        return M




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
        
        
        # Get the unit vector from observatory to object: Is of shape = (len(times_tdb) , 3)
        UV = self.generate_UnitVector(  times_tdb ,
                                        observatoryXYZ,
                                        APPROX = APPROX )

        # Calc the HP from the UV and return
        return healpy.vec2pix(self.HP_nside, UV[0], UV[1], UV[2], nest=True if self.HP_order=='nested' else False )


    def generate_UnitVector( self, times_tdb , observatoryXYZ,  APPROX = False ,
                                                                DELTASWITCH_XYZ = False,
                                                                DELTASWITCH_UVW = False,
                                                                delta=np.array([0,0,0]) ):
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
            - shape = (len(times_tdb) , )
            
            delta: np.array
             - Assume
            
            return:
            -------
            unit-vectors: np.array 
             - apparent UnitVector from specified observatory-posn(s) at given time(s)
            
        '''
        # Ensure that the input shape is the same as will be returned by generate_XYZ ...
        assert observatoryXYZ.shape == ( len(times_tdb) , 3)

        # Get the LTT-corrected position
        # - We allow for the possibility of *NOT* iterating (i.e. doing an approx calc.)
        n_iterations    = 1 if APPROX else 3
        
        # Init light-delay to 0
        lightDelay      = np.zeros( len(times_tdb) )

        for i in range(n_iterations):

            # Calculate delayed time (of emission)
            # N.B. delayedTimes.shape = (len(times_tdb) )
            delayedTimes    = times_tdb if i == 0 else times_tdb - lightDelay
            
            # Extract posn of objects at each delayed-time
            # N.B. objectXYZ.shape    = (len(times_tdb) , 3)
            objectXYZ       = self.generate_XYZ( delayedTimes )

            #if DELTASWITCH_XYZ :
            #    '''XYZ: Assume simple translation: shift@delayedTimes == shift@times_tdb'''
            #    objectXYZ += np.stack([delta for i in range(objectXYZ.shape[0])], axis=1)
            
            # Calculate relative sepn-vector from observatory-to-object
            sepn_vectors    = objectXYZ - observatoryXYZ

            # Calculate distance to object at each time.
            # NB1 sepn_vectors.shape == (len(times_tdb) , 3)
            # NB2 d.shape            == (len(times_tdb),)
            d               = np.linalg.norm(sepn_vectors, axis=1)

            # Calculate light-travel-time:
            # lightDelay.shape == (len(times_tdb),)
            lightDelay      = d / (astropy.constants.c * self.secsPerDay / astropy.constants.au ).value


            #if DELTASWITCH_UVW :
            #    '''UVW: Seems more complex than XYZ: need lightDelay?? => Extra Iteration??'''
            #    deltaDelay = (lightDelay*delta[:,None]).T
            #    #objectXYZ += deltaDelay
            #    #sepn_vectors = objectXYZ - observatoryXYZ
            #    sepn_vectors += deltaDelay
            #    d             = np.linalg.norm(sepn_vectors, axis=1)
            #    lightDelay    = d / (astropy.constants.c * self.secsPerDay / astropy.constants.au ).value

        # Return unit-vectors: of shape = (N_times , 3)
        return sepn_vectors / d[:,None]

    def dUVdXYZUVW( self, times_tdb , observatoryXYZ , d = 1e-8):
        '''
            Gradient of the UnitVector w.r.t. the Cartesian X,Y,Z positions
            
            inputs:
            -------
            
            returns:
            --------
        '''
        # Generating displacement vectors ... Are of shape = (len(times_tdb) , 3)
        _dX = ( self.generate_UnitVector( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH_XYZ = True, delta=np.array([d, 0, 0]) ) \
               -self.generate_UnitVector( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH_XYZ = True, delta=np.array([-d,0, 0]) ) )
        _dY = ( self.generate_UnitVector( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH_XYZ = True, delta=np.array([0, d, 0]) ) \
               -self.generate_UnitVector( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH_XYZ = True, delta=np.array([0,-d, 0]) ) )
        _dZ = ( self.generate_UnitVector( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH_XYZ = True, delta=np.array([0, 0, d]) ) \
               -self.generate_UnitVector( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH_XYZ = True, delta=np.array([0, 0,-d]) ) )
        _dU = ( self.generate_UnitVector( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH_UVW = True, delta=np.array([d, 0, 0]) ) \
               -self.generate_UnitVector( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH_UVW = True, delta=np.array([-d,0, 0]) ) )
        _dV = ( self.generate_UnitVector( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH_UVW = True, delta=np.array([0, d, 0]) ) \
               -self.generate_UnitVector( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH_UVW = True, delta=np.array([0,-d, 0]) ) )
        _dW = ( self.generate_UnitVector( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH_UVW = True, delta=np.array([0, 0, d]) ) \
               -self.generate_UnitVector( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH_UVW = True, delta=np.array([0, 0,-d]) ) )

        return np.stack(np.array( (_dX, _dY, _dZ, _dU, _dV, _dW) ), axis=1).T / (2*d)
    
    def covUV(self, times_tdb , observatoryXYZ ):
        '''
            Evaluate the covariance in unit vectors
            This is calculated using the covariance in XYZ & the gradient of the UV components w.r.t. XYZ
            '''
        dUV     = self.dUVdXYZ( times_tdb , observatoryXYZ )
        cov_XYZ = self.covXYZ( times_tdb  )
        return np.array( [ np.linalg.multi_dot([dUV[i].T , cov_XYZ[i], dUV[i]]) for i in range(len(dUV)) ] )

    def generate_RaDec(self,  times_tdb  , observatoryXYZ,  APPROX = False ,
                                                            DELTASWITCH_XYZ = False,
                                                            DELTASWITCH_UVW = False,
                                                            delta=np.array([0,0,0])):
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
            - shape = (len(times_tdb), 3)
            
            return:
            -------
            RA_Dec: np.array
            - *** degrees ***
            
        '''
        # Get the unit vector from observatory to object: Is of shape = ( len(times_tdb) , 3)
        UV = self.generate_UnitVector(times_tdb ,
                                      observatoryXYZ,
                                      APPROX = APPROX,
                                      DELTASWITCH_XYZ=DELTASWITCH_XYZ,
                                      DELTASWITCH_UVW=DELTASWITCH_UVW,
                                      delta = delta)
        
        # Convert from unit-vector to RA, Dec and then return
        # NB: Taking transpose, T, makes retuurned shape ==  ( len(times_tdb) , 2)
        return np.array(healpy.vec2ang( UV , lonlat = True )).T

    def dRaDecdXYZUVW( self, times_tdb , observatoryXYZ ,         d = 1e-8):
        '''
            Gradient of Ra & Dec w.r.t. the Cartesian X,Y,Z positions
            
            Assume:
                d = 1e-8    => 1e-8 AU      => 1.5e3m   == 1.5km for XYZ
                            => 1e-8 AU/day  => 0.017m/s          for UVW
            
            inputs:
            -------
            
            returns:
            --------
            '''
        _dX = ( self.generate_RaDec( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH_XYZ = True, delta=np.array([d, 0, 0]) ) \
               -self.generate_RaDec( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH_XYZ = True, delta=np.array([-d,0, 0]) ) )
        _dY = ( self.generate_RaDec( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH_XYZ = True, delta=np.array([0, d, 0]) ) \
               -self.generate_RaDec( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH_XYZ = True, delta=np.array([0,-d, 0]) ) )
        _dZ = ( self.generate_RaDec( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH_XYZ = True, delta=np.array([0, 0, d]) ) \
               -self.generate_RaDec( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH_XYZ = True, delta=np.array([0, 0,-d]) ) )
        _dU = ( self.generate_RaDec( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH_UVW = True, delta=np.array([d, 0, 0]) ) \
               -self.generate_RaDec( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH_UVW = True, delta=np.array([-d,0, 0]) ) )
        _dV = ( self.generate_RaDec( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH_UVW = True, delta=np.array([0, d, 0]) ) \
               -self.generate_RaDec( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH_UVW = True, delta=np.array([0,-d, 0]) ) )
        _dW = ( self.generate_RaDec( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH_UVW = True, delta=np.array([0, 0, d]) ) \
               -self.generate_RaDec( times_tdb , observatoryXYZ, APPROX = True , DELTASWITCH_UVW = True, delta=np.array([0, 0,-d]) ) )

        return np.stack(np.array( (_dX, _dY, _dZ, _dU, _dV, _dW) ), axis=1).T / (2*d)
    
    
    def covRaDec(self, times_tdb , observatoryXYZ ):
        '''
            Evaluate the covariance in RA, Dec 
            This is calculated using the covariance in XYZ & the gradient of RA,Dec w.r.t. XYZ
        '''
        dRaDecdXYZUVW   = self.dRaDecdXYZUVW( times_tdb , observatoryXYZ )
        covXYZUVW       = self.covXYZUVW( times_tdb  )
        
        # Is there a way to make this more efficient?
        # - Similar to the nbody code ...
        # - CoV_t = np.linalg.inv( GammaStack_t )
        return np.array( [ np.linalg.multi_dot([dRaDecdXYZUVW[i].T , covXYZUVW[i], dRaDecdXYZUVW[i]]) for i in range(len(dRaDecdXYZUVW)) ] )

    

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
             - NB1: XYZ_posns.shape = (len(times_tdb) , 3)
                    This ensures consistency with NbodySim & coco.equatorial_helio2bary
             - NB2: no need for N_particle dimension, as MSC is only for single object
        '''
        return self.evaluate_components( times_tdb  , component_slice_spec=self.XYZ_slice_spec )
        
        
    def generate_XYZUVW( self, times_tdb  ):
        '''
            Evaluate the XYZUVW components at the supplied times

            Not used
        '''
        return self.evaluate_components( times_tdb  , component_slice_spec=self.XYZUVW_slice_spec )



    def covXYZUVW( self, times_tdb  ):
        '''
            Evaluate the covariance in XYZUVW coord-components at the supplied times
            
            Convenience wrapper around *evaluate_components()* func
            Uses "self.covXYZUVW_slice_spec" to select required coefficients
            
            returns 
            -------
            np.ndarray
             - shape = ( len(times_tdb) , 6, 6 )
        '''
        # Select/evaluate the appropriate covariance components*                              *
        cov = self.evaluate_components( times_tdb , component_slice_spec = self.covXYZUVW_slice_spec  )
        
        # Reshape to get square matrix for each time
        return self._make_square(cov)


    def dXYZdt(self,  times_tdb  , dt=1e-5):
        ''' 
            Calculate the gradient in XYZ at supplied times
            
            Not Used
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
             - Making the design decision to return the array in the shape (N_times, N_components)
             - This is to make it as similar as possible to the shape that comes out of reboundx/NbodySim
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

        
        # Taking the transpose to return the array in the shape (N_times, N_components)
        # - This is to make it as similar as possible to the shape that comes out of reboundx/NbodySim
        return combinedComponents.T





# End Of File
# --------------------------------------------------------------
