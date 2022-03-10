# -*- coding: utf-8 -*-
# /cheby_checker/cheby_checker/nbody.py

'''
----------------------------------------------------------------------------
nbody is a wrapper around a REBOUNDX integrator:
reboundx/examples/ephem_forces

Nov 2020
Mike Alexandersen & Matthew Payne & Matthew Holman

This module provides functionalities to
(a) parse an input orbit (e.g. from ORBFIT),
(b) do the nbody integration
(c) read the output from the nbody simulations

----------------------------------------------------------------------------
'''


# Import third-party packages
# -----------------------------------------------------------------------------
import sys
import os
import numpy as np
from astropy.time import Time
import getpass
import json
import itertools
import copy

# Import MPC packages
# -----------------------------------------------------------------------------
from mpc_orb.parse import MPCORB


# Import neighbouring packages
# -----------------------------------------------------------------------------
from examples.ephem_forces.ephem_forces import production_integration_function_wrapper

# old conversion library
from . import MPC_library as mpc

from . import coco
from .decorators import timer


# Constants and stuff
# -----------------------------------------------------------------------------
DATA_PATH = os.path.realpath(os.path.dirname(__file__))
DATA_DIR = os.path.join(os.path.dirname(DATA_PATH), 'dev_data')
au_km = 149597870.700  # This is now a definition


# Top-line parameters for parse_nbody_txt/create_nbody_txt
# MJP 2020-11-10: These parameters are *only* used by parse_nbody_txt/create_nbody_txt
# MJP 2020-11-10: (see bottom of this file)
# MJP 2020-11-10: ** It is unclear to me whether parse_nbody_txt/create_nbody_txt need to exist **
# MJP 2020-11-10: ** Could fold into NbodySim as method(s) **
# MJP 2020-11-10: ** Could remove entirely if not really used anymore **
# --------------------------------------------------------------
# Field names
object_name     = 'unpacked_designation'
time_fieldname  = 'MJD_TDB' # To be clarified if this is correct
coord_names     = ['x','y','z','vx','vy','vz']

# This generates a 'triangular' list of 21-combinations of coord-names
# ['x_x', 'x_y', 'x_z', 'x_vx', 'x_vy', 'x_vz',
#         'y_y', 'y_z', 'y_vx', 'y_vy', 'y_vz',
#                'z_z', 'z_vx', 'z_vy', 'z_vz',
#                      'vx_vx', 'vx_vy', 'vx_vz',
#                               'vy_vy', 'vy_vz',
#                                        'vz_vz']
covar_names     = [ "_".join([coord_names[i], coord_names[j]]) for i in range(len(coord_names)) for j in range(i,len(coord_names))  ]

# Number of fields
nCovars     = len( covar_names )
nComponents = nCovars + len(coord_names)
nFields     = nComponents + 1
# --------------------------------------------------------------



class NbodySim():
    '''
    Provides methods to call reboundx-ephemeris NBody integrator
    '''

    def __init__( self,):

        # class variables: declared @ run-time and/or parsed from input
        self.tstart             = None
        self.tstop              = None

        self.geocentric         = False

        self.integration_epoch  = None
        self.mpcorb_list        = []

        self.input_n_particles  = None

        self.save_output        = False
        self.save_output_file   = None
        self.verbose            = False

        self.CHECK_EPOCHS       = True

        self.unpacked_primary_provisional_designation_list = []

        # class variables: variables that will be used to hold the elements
        # - They get populated by *parse_...* & *make_bary_equatorial*
        self.helio_ecl_vec_EXISTS   = False
        self.helio_ecl_vec          = None
        self.helio_ecl_cov_EXISTS   = False
        self.helio_ecl_cov          = None
        self.bary_eq_vec_EXISTS     = False
        self.bary_eq_vec            = None
        self.bary_eq_cov_EXISTS     = False
        self.bary_eq_cov            = None
        self.non_grav_EXISTS        = False
        self.non_grav_dict_list     = []

        #  class variables: will be used to hold output
        self.output_times       = None
        self.output_states      = None
        self.output_covar     = None


    # --------------------------------------------
    # Top-Level Function(s): (expected to be called by user)
    # --------------------------------------------
    #@timer
    def run_mpcorb(self,**kwargs):
        '''
        Run the nbody integrator with the parsed input.
         - Performs integration by calling *integration_function*
         
        Inputs/Requires:
        -------

        Populates:
        --------
        
        Example:
        --------
        NbodySim().run_mpcorb(  tstart = ... ,
                                tstop  = ... ,
                                mpcorb_list = [],
                                )
        '''
        
        # Parse the provided input and populate the appropriate internal variables.
        # Also checks that all necessary variables have been provided
        self._parse_inputs_run_mpcorb(**kwargs)

        # Populate barycentric quantities from the input heliocentric quantities
        self.make_bary_equatorial()

        # Run the Integration
        self._run_integration()
        
        # Save output if requested
        #if save_output:
        #    self.save_output()
            
        return True


    # --------------------------------------------
    # Functions to parse user-inputs
    # --------------------------------------------
    #@timer
    def _parse_inputs_run_mpcorb(self, **kwargs):
        """ Parse inputs based on mpc_orb.json format
        
        Inputs/Requires:
        -------

        Populates:
        --------
        
        Example:
        --------

        """
        # Check the critical_inputs have been specified ...
        critical_inputs = ['tstart', 'tstop', 'mpcorb_list']
        for key in critical_inputs:
            assert key in kwargs and  kwargs[key] is not None
            self.__dict__[key] = kwargs[key]
            
        assert isinstance(self.tstart,(float,int))
        assert isinstance(self.tstop, (float,int))
        assert isinstance(self.mpcorb_list, (list,tuple,np.ndarray) )

        # Get the number of particles
        self.input_n_particles = len(self.mpcorb_list)
        assert self.input_n_particles, f'self.input_n_particles={self.input_n_particles}'
        
        # Parse the mpc_orb.json input ...
        for _ in self.mpcorb_list:
            self._parse_orbfit_json( _ )

        # Basic error-checking ...
        if (self.input_n_particles is None) or (self.input_n_particles == 0) or (self.helio_ecl_vec is None) :
            raise TypeError("Error parsing from input ")
            
        # Don't need to return anything ...
        return True

  
            
    #@timer
    def _parse_orbfit_json(self, mpcorb_file_or_dict ):
        '''
        Parse a file containing OrbFit elements for a single object & epoch.
        Assumes input is in standard mpcorb-json format


        Input:
        -------
        mpcorb_file_or_dict : Filepath or Dictionary
         - Valid mpcorb.json format dictionary or json-file

        Populates:
        --------
        self.integration_epoch
        
        self.helio_ecl_vec
        self.helio_ecl_vec_EXISTS
        
        self.helio_ecl_cov
        self.helio_ecl_cov_EXISTS
        
        self.non_grav_dict_list
        self.non_grav_EXISTS
        
        '''
        
        # Call -----------  MPCORB -----------------
        # NB(1): This parses & validates ...
        # NB(2): MPCORB Handles either FILE or DICTIONARY input
        M = MPCORB(mpcorb_file_or_dict)
        
        # ----------- DESIGNATION -------------------------
        self.unpacked_primary_provisional_designation_list.append( M.designation_data["unpacked_primary_provisional_designation"] )

        # ----------- TIME -------------------------
        # Using Astropy.time for time conversion,
        # because life's too short for timezones and time scales.
        epoch,timesystem    = M.epoch_data["epoch"],M.epoch_data["timesystem"]
        
        allowed_timesystems = {'TDT':'tt'}
        assert timesystem in allowed_timesystems
        this_mpcorb_epoch = Time(float(epoch), format='mjd', scale=allowed_timesystems[timesystem] )

        # In production, we likely want to be using "standard-epochs" ending in 00.
        # - This will allow multiple objects to be simultaneously integrated
        if self.CHECK_EPOCHS==True and "00." not in str(epoch):
            raise TypeError("Supplied file may not contain standard epochs:"
                            f"epoch= {epoch}")
        
        # If we are processing multiple files, then
        # check that all supplied epochs are the same
        if self.integration_epoch is None :
            self.integration_epoch = this_mpcorb_epoch
        else:
            if self.integration_epoch != this_mpcorb_epoch:
                raise TypeError("The epochs vary from file-to-file: "
                                f"jsonfile_dictionary= {jsonfile_dictionary}")


        # ----------- ELEMENTS/VARIABLES -----------
        # Get Cartesian dict out of the overall jsonfile_dictionary
        cart_dict = M.CAR

        # Parse the elements and stack with any others extracted from previous files
        # - The elements are always of length 6
        numparams                   = cart_dict['numparams']
        element_array               = np.atleast_2d( cart_dict['element_array'] ).reshape(1,6)
        self.helio_ecl_vec_EXISTS   = True #  MPCORB *REQUIRES* params ...
        if self.helio_ecl_vec is None:
            self.helio_ecl_vec = element_array
        else :
            self.helio_ecl_vec = np.append(self.helio_ecl_vec , element_array, axis=0)

        # ----------- CHECK NUMBER OF PARAMS -----------
        '''
        I think that the stacking may well screw-up if there are different numbers
        of parameters for different objects in the same integration
        I should check for that and emit an error/warning ...
        '''

        
        # ----------- COVARIANCE  -------------------
        # Parse jsonfile_dictionary to get the cartesian covariance matrix
        # The cov matrix will be of dimension N^2, with N>=6, as the non-gravs are incluuded
        local_helio_ecl_cov         = np.atleast_3d(cart_dict['covariance_array']).reshape(1,numparams,numparams )
        self.helio_ecl_cov_EXISTS   = True #  MPCORB *REQUIRES* cov ...
        if self.helio_ecl_cov is None:
            self.helio_ecl_cov = local_helio_ecl_cov
        else :
            self.helio_ecl_cov = np.append( self.helio_ecl_cov, local_helio_ecl_cov, axis=0)

                
        # ----------- NONGRAVS -----------------------
        self.non_grav_dict_list.append( M.nongrav_data )
        self.non_grav_EXISTS = np.any( _["non_gravs"] for _ in self.non_grav_dict_list )
        
        
        # ----------- return -----------------------
        return True
        
                
    #@timer
    def make_bary_equatorial(self):
        '''
        Transform heliocentric-ecliptic coordinates into
        barycentric equatorial coordinates
        
        requires:
        ----------
        self.helio_ecl_vec_EXISTS   : Boolean
        self.helio_ecl_vec          : 1D np.ndarray
        self.helio_ecl_cov_EXISTS   : Boolean
        self.helio_ecl_cov          : 2D np.ndarray

        populates:
        ----------
        self.bary_eq_vec_EXISTS     = Boolean
        self.bary_eq_vec            = 1D np.ndarray
        self.bary_eq_cov_EXISTS     = Boolean
        self.bary_eq_cov            = 2D np.ndarray
        '''
        if self.helio_ecl_vec_EXISTS :
            # Transform the helio-ecl-coords to bary-eq-coords
            # NB 2-step transformation for the vector (posn,vel)
            self.bary_eq_vec   = coco.equatorial_helio2bary(
                                    coco.ecliptic_to_equatorial(self.helio_ecl_vec),
                                    float(self.integration_epoch.tdb.to_value('jd'))
                                )
            # Set boolean as well (not sure if we'll really use these ...)
            self.bary_eq_vec_EXISTS = True

        if self.helio_ecl_cov_EXISTS:
            # Only need to do a rotation for the CoV
            # MJP 2022-02-05 : Need to "slice" to ensure only grav-components get rotated
            self.bary_eq_cov = coco.ecliptic_to_equatorial(self.helio_ecl_cov)
        
            # Set booleans as well (not sure if we'll really use these ...)
            self.bary_eq_cov_EXISTS = True

        if not self.helio_ecl_vec_EXISTS and not self.helio_ecl_cov_EXISTS:
            raise TypeError("There does not seem to be any valid helio_ecl to transform into bary_eq")
            
        return True
        
    

    """
    def _parse_input_type(self,input):
        '''
        Because we allow file(s) or file-content(s), need to parse
        
        input :
        -------
        string or list-of-strings
         - each string must be either:
           --- valid file-path
           --- valid file-contents (e.g. orbfit output)
        
        returns:
        --------
        list :
         - each list-item consists of file-contents
         
        '''

        # Make a list
        if isinstance(input, str):
            input = [input]
        
        # Check we have a list
        assert isinstance(input, list)
        
        # If we have filepaths, get contents
        # otherwise, assume these already are data similar to file-contents
        contents = []
        for item in input:
            if os.path.isfile(item):
                # Try to open as JSON  ...
                try:
                    with open(item) as json_file:
                        contents.append( json.load(json_file) )

                # If json-fails, then open as standard text file
                except:
                    with open(item,'r') as f:
                        contents.append(f.readlines())
            else:
                contents.append(input)
                
        return contents
    """


                      
    """
    def parse_orbfit_felfile_txt(self, felfile_contents, CHECK_EPOCHS=True):
        '''
        Parse a file containing OrbFit elements for a single object & epoch.
        Assumes input is in the "old" ~text format (filetype ~ .eq0/.eq1)

        Inputs:
        -------
        felfile_contents : string, *contents* of fel/eq formatted OrbFit output

        Populates:
        --------
        self.helio_ecl_vec_EXISTS   : Boolean
        self.helio_ecl_vec          : 1D np.ndarray
        self.helio_ecl_cov_EXISTS   : Boolean
        self.helio_ecl_cov          : 1D np.ndarray
        self.time                   : astropy Time object
        '''
        
        # Parse the contents of the orbfit file
        try:

            # Get Cartesian Elements out of the file contents
            # NB - There can be multiple blocks, so we get the LAST block
            cart_head = '! Cartesian position and velocity vectors\n'
            carLoc   = len(felfile_contents) - 1 - list(reversed(felfile_contents)).index(cart_head)
            carEls   = felfile_contents[carLoc:carLoc + 25]
            (_, car_x, car_y, car_z, car_dx, car_dy, car_dz
                       ) = carEls[1].split()

            # Form the heliocentric ecliptic cartesian coefficients in to an array
            # and stack with any others extracted from previous files
            # NB ... self.helio_ecl_vec.ndim ==2, self.helio_ecl_vec.shape = (N,6)
            local_helio_ecl_vec = np.atleast_2d([float(car_x),   float(car_y),  float(car_z), \
                                                float(car_dx),  float(car_dy), float(car_dz)] )
            if self.helio_ecl_vec is None:
                self.helio_ecl_vec = local_helio_ecl_vec
            else :
                self.helio_ecl_vec = np.append(self.helio_ecl_vec , local_helio_ecl_vec, axis=0)
                                                    
                                            
            self.helio_ecl_vec_EXISTS = True
                                                      
            # Using Astropy.time for time conversion,
            # because life's too short for timezones and time scales.
            _, mjd_tdt, _   = carEls[2].split()
            local_Time      = Time(float(mjd_tdt), format='mjd', scale='tt')
            
            # In production, we likely want to be using "standard-epochs" ending in 00.
            # - This will allow multiple objects to be simultaneously integrated
            if CHECK_EPOCHS==True and "00." not in mjd_tdt:
                raise TypeError("Supplied file may not contain standard epochs:"
                                f"mjd_tdt= {mjd_tdt}")
                
            if self.time is None :
                self.time = local_Time
            else:
                # If we are processing multiple files, then
                # check that all supplied epochs are the same
                if self.time != local_Time:
                    raise TypeError("The epochs vary from file-to-file: "
                                    f"felfile_contents= {felfile_contents}")
                    
            # Parse carEls (the contents of the orbfit file) to get
            # the cartesian covariance matrix
            # NB local_helio_ecl_cov.ndim ==3 , local_helio_ecl_cov.shape == (1,6,6)
            self.helio_ecl_cov_EXISTS, local_helio_ecl_cov = _parse_Covariance_List(carEls)
            if self.helio_ecl_cov is None:
                self.helio_ecl_cov = local_helio_ecl_cov
            else :
                self.helio_ecl_cov = np.append( self.helio_ecl_cov, local_helio_ecl_cov  )

        except Exception as e:
            print( "Supplied felfile_contents could not be processed.\n"
                    f"felfile_contents={felfile_contents}\n"
                    f"e = {e}\n")
    """
    

    # --------------------------------------------
    # Functions to run NBody integration & parse results
    # --------------------------------------------

    #@timer
    def _run_integration( self,):
        '''
        Run the nbody integrator with the parsed input.
        
        Designed as INTERNAL function: NOT expected to be called directly by the user
        NB: Separating *_run_integration* from (e.g.) run_mpcorb to allow for multiple top level run_* types (e.g. run_array)

        Input:
        ------
        Assumes all required class variables have been populated at preceeding steps
        Requires:
        
        self.tstart : float
         -
        self.tstop  : float
         -
        self.integration_epoch : Astropy.Time object
         -
        self.geocentric = boolean,
         - use geo- (True) or heliocentric (False)
         
        self.n_particles : int
         -
        self.bary_eq_vec : np.ndarray
         - 2D, shape = (N_particles, 6)
         - barycentric equatorial elements (xyzuvw) for N_particles to be integrated
         
        self.non_grav_dict_list = list
         -
        self.epsilon = 1e-8
         -
        self.tstep_max = 32
         =

        Output:
        -------
        reparsed_input      = numpy array
         - input elements, reparsed into array
        n_particles         = integer,
         - the input number of particles
        times               = numpy array,
         - all the output times (including sub-steps)
        states              = numpy array,
         - output elements of dimensions (n_times, n_particles_out, 6)
        output_covariance   = numpy array,
         - output elements of dimensions (n_times, n_particles_out*6, 6)
        n_times             = integer,
         - number of time outputs
        n_particles_out     = integer,
         - number of output particles (different why?)
         
         
         
        *** NOTES ON PROBLEMS WITH MULTI-PARTICLE INPUT ***
        When testing on 4-particle input, I found that the integration ran wild,
        with the solutions zooming off to 10,000au (individually they stayed @ 3au)
        The problem was found to be that the input matrix, self.bary_eq_vec, (2D array)
        gets transposed inside C, but only when passed in as-is.
         - A hard-pasted copy of the array gets passed in just fine, without transposition
        Checking things in python before shows *NO* difference between the self.bary_eq_vec
        and the hardpasted copy, so the difference is presumably something to do with
        the internal memory array adopted in numpy arrays in different circumstances, and then
        how that is translated to C
        SPECULATION : Is it something to do with bary_eq_vec being reshaped / reduced in dimension
        at some point in the preceeding code? Does this affect memory layout?
        
        hardpasted = np.array(
        [[ 2.4260322844717380e-01, -3.0090867545871860e+00, -1.1363474792630044e+00,
           9.1908324088463677e-03,  2.0472935970702337e-04, -1.5374992052716444e-03],
         [-2.7888126437563829e+00,  1.7058528476705044e+00,  3.9449682717556311e-01,
          -5.2133449765416028e-03, -7.5295451342754713e-03, -7.8760899838941603e-04],
         [-5.1094716747791424e-01, -3.1254423057577734e+00, -2.0995426392478924e+00,
           7.8191160668487376e-03, -3.3055503238910031e-04,  9.7078064274683739e-04],
         [-2.1582274449687930e+00,  2.0100970703192091e+00,  3.1933652106784710e-01,
          -5.5592970943310026e-03, -7.8839265248180895e-03, -3.4116741888535933e-03]])

        The Good Inputs
        INPUT SET TO HARDPASTED NUMBERS ...

        *** f= *** 0.7916875000004463

         ---ccc--- 0.2426032284471738 -3.0090867545871860 -1.1363474792630044 0.0091908324088464 0.0002047293597070 -0.0015374992052716

         ---ccc--- -2.7888126437563829 1.7058528476705044 0.3944968271755631 -0.0052133449765416 -0.0075295451342755 -0.0007876089983894

         ---ccc--- -0.5109471674779142 -3.1254423057577734 -2.0995426392478924 0.0078191160668487 -0.0003305550323891 0.0009707806427468

         ---ccc--- -2.1582274449687930 2.0100970703192091 0.3193365210678471 -0.0055592970943310 -0.0078839265248181 -0.0034116741888536



        The Bad Inputs
        INPUTS AS ORIGINAL ...

        *** f= *** 0.7916875000004463

         ---ccc--- 0.2426032284471738 -2.7888126437563829 -0.5109471674779142 -2.1582274449687930 -3.0090867545871860 1.7058528476705044

         ---ccc--- -3.1254423057577734 2.0100970703192091 -1.1363474792630044 0.3944968271755631 -2.0995426392478924 0.3193365210678471

         ---ccc--- 0.0091908324088464 -0.0052133449765416 0.0078191160668487 -0.0055592970943310 0.0002047293597070 -0.0075295451342755

         ---ccc--- -0.0003305550323891 -0.0078839265248181 -0.0015374992052716 -0.0007876089983894 0.0009707806427468 -0.0034116741888536
        '''
        
        
        # Now run the nbody integrator:
        if self.verbose:
            print(f"self.tstart={self.tstart}, self.tstop={self.tstop}, epoch=float(self.integration_epoch.tdb.to_value('jd'))={float(self.integration_epoch.tdb.to_value('jd'))},")

            
        self.output_times, self.output_states, partial_derivatives_wrt_state, partial_derivatives_wrt_NG, return_value = \
            production_integration_function_wrapper(    self.tstart,
                                                        self.tstop,
                                                        float(self.integration_epoch.tdb.to_value('jd')) , # Converting Astropy.Time ...
                                                        self.bary_eq_vec,
                                                        non_grav_dict_list = [],#self.non_grav_dict_list,
                                                        tstep = 20,
                                                        geocentric = 0,#self.geocentric,
                                                        epsilon = 1e-7,  #1e-8
                                                        tstep_min = 0.1, #0.02
                                                        tstep_max = 32.)

        # Reshape the partial derivative arrays
        partial_derivatives_wrt_state, partial_derivatives_wrt_NG = self.reshape_partial_deriv_arrays(  partial_derivatives_wrt_state,partial_derivatives_wrt_NG)
        self.partial_derivatives_wrt_state = partial_derivatives_wrt_state
        
        
        # Calculate the covariance matrix (at each timestep) from the \partial X / \partial X_0 data
        if self.bary_eq_cov_EXISTS is not None:
            self.output_covar = self._get_covariance_from_tangent_vectors(  self.bary_eq_cov,
                                                                                partial_derivatives_wrt_state,
                                                                                partial_derivatives_wrt_NG = partial_derivatives_wrt_NG )
        else:
            self.output_covar = None
        
        return True
                

    #@timer
    def reshape_partial_deriv_arrays( self,  partial_derivatives_wrt_state,  partial_derivatives_wrt_NG):
        '''
        (1) partial_derivatives_wrt_state
        From: partial_derivatives_wrt_state.shape = (N_times, 6*n_particles, 6)
        To  : partial_derivatives_wrt_state.shape = (N_times, n_particles, 6, 6)
        
        (2) partial_derivatives_wrt_NG
        *** NOT YET IMPLEMENTED partial_derivatives_wrt_NG ***
        '''
        # (1) (N_times, 6*n_particlea, 6) -> (N_times, n_particlea, 6, 6)
        #print('partial_derivatives_wrt_state: input shape = ', partial_derivatives_wrt_state.shape)
        partial_derivatives_wrt_state = partial_derivatives_wrt_state.reshape(partial_derivatives_wrt_state.shape[0],-1,6,6)
        #print('partial_derivatives_wrt_state: final shape = ', partial_derivatives_wrt_state.shape)
        #print(partial_derivatives_wrt_state[-1])
        
        # (2)
        if partial_derivatives_wrt_NG is not None:
            raise Error('Have not coded partial_derivatives_wrt_NG into _get_covariance_from_tangent_vectors ')

        return partial_derivatives_wrt_state , partial_derivatives_wrt_NG
        
        
    #@timer
    def _get_covariance_from_tangent_vectors(self, init_covariances, partial_derivatives_wrt_state ,  partial_derivatives_wrt_NG=None ):
        '''
        Follow Milani et al 1999
        Gamma_t = [partial X / partial X_0] Gamma_0 [partial X / partial X_0]^T
        
        NB
        init_covariances.shape              =  (#particles, 6, 6)
        partial_derivatives_wrt_state.shape =  (#times, #particles, 6, 6)

        '''
        
        if partial_derivatives_wrt_NG is not None :
            raise Error('Have not coded partial_derivatives_wrt_NG into _get_covariance_from_tangent_vectors ')

        # init_covariances.shape              =  (#particles, 6, 6)
        # partial_derivatives_wrt_state.shape == (#times, #particles, 6, 6)
        assert init_covariances.ndim == 3
        assert partial_derivatives_wrt_state.ndim == 4
        assert partial_derivatives_wrt_NG is None or partial_derivatives_wrt_NG.ndim  == 4
        assert init_covariances.shape[-2:] == (6,6) # GRAVITY-ONLY ...

        # Take the inverse of the covariance matrix to get the normal matrix
        # NB1: We make a stack of identical matricees to use in the matrix multiplication below
        # NB2: We want the stack height == Number of Times == partial_derivatives_wrt_state.shape[0]
        #     => GammaStack0.shape == partial_derivatives_wrt_state.shape == (#times, #particles, 6, 6)
        Gamma0      = np.linalg.inv(init_covariances)
        GammaStack0 = np.tile( Gamma0, (partial_derivatives_wrt_state.shape[0],1,1,1) )

        # We need each of the individual pd arrays to be individually transposed
        # NB, The "(0,1,3,2)" tuple fixes dimensions 0 & 1 , but indicates that dimensions 2 & 3 will be swapped/transposed
        #     Because this just swaps the 6x6 dimensions, the shape remains == (#times, #particles, 6, 6)
        pds_transposed  = partial_derivatives_wrt_state.transpose( (0,1,3,2) )
        
        # Do matrix multiplication: using the pd's to get the CoV as a func of time
        # NB matmul/@ automagically knows how to work on a stack of matricees
        # - https://stackoverflow.com/questions/34142485/difference-between-numpy-dot-and-python-3-5-matrix-multiplication
        # GammaStack_t.shape = (#times, #particles, 6, 6)
        GammaStack_t    = pds_transposed @ GammaStack0 @ partial_derivatives_wrt_state

        # Magically, np.linalg.inv also knows how to deal with a stack of arrays/matrices
        # - https://stackoverflow.com/questions/11972102/is-there-a-way-to-efficiently-invert-an-array-of-matrices-with-numpy
        # CoV_t.shape = (#times, #particles, 6, 6)
        # The try-except loop is required to catch occasional problems with ...LinAlgError("Singular matrix")... due to singular input matrix
        try:
            CoV_t           = np.linalg.inv( GammaStack_t )
        except:
            CoV_t           = np.linalg.pinv( GammaStack_t )

        return CoV_t
            

    # -----------------------------------
    # Save Funcs ...
    # -----------------------------------
    #@timer
    def save_output(self, output_file='simulation_states.dat'):
        """
        Save all the outputs to file.

        Inputs:
        -------
        output_file : string, filename to write elements to.

        The file is overwritten if it already exists.
        """
        outfile = open(output_file, 'w')

        outfile.write(f'\b]')
        outfile.write(f'\n#Input N_particles: {self.input_n_particles:}')
        outfile.write('\n#Start time, timestep, time range: '
                      f'{self.time_parameters}')
        outfile.write(f'\n#Output N_times: {self.output_n_times:}')
        outfile.write(f'\n#Output N_particles: {self.output_n_particles:}')
        outfile.write('\n#')
        outfile.write('\n#Time               ')
        for j in np.arange(self.output_n_particles):
            outfile.write('x                  y                  '
                          'z                   dx                  '
                          '  dy                    dz                  ')
        for i, timei in enumerate(self.output_times):
            outfile.write(f'\n{timei:} ')
            for j in np.arange(self.output_n_particles):
                for k in np.arange(6):
                    outfile.write(f'{self.output_vectors[i, j, k]:} ')
        outfile.write('\n#End')


    #@timer
    def save_elements(self, save_file='save_file.tmp'):
        """
        Save the barycentric equatorial cartesian elements to file.

        Inputs:
        -------
        output_file : string, filename to write elements to.

        The file is overwritten if it already exists.
        """
        self.tstart = self.time.tdb.jd
        outfile = open(save_file, 'w')
        outfile.write(f"tstart {self.tstart:}\n")
        outfile.write("tstep +20.0\n")
        outfile.write("trange 600.\n")
        outfile.write("geocentric 0\n")
        outfile.write("state\n")
        
        # For whatever reason, we are writing this over two lines
        # - perhaps to compare against JPL?
        for vec in self.bary_eq_vec:
            for n,coeff in enumerate(vec):
                suffix = '\n' if n in [2,5] else ''
                outfile.write(f"{coeff: 18.15e} " + suffix)



    # -----------------------------------
    # Misc Funcs ...
    # -----------------------------------
    #@timer
    def _parse_Covariance_List(Els):
        '''
        Convenience function for reading and splitting the covariance
        lines of an OrbFit file.
        Not intended for user usage.
        '''
        # Set-up array of zeroes
        # NB: CoV.ndim ==3 , CoV.shape == (1,6,6)
        CoV        = np.zeros( (1,6,6) )
        CoV_EXISTS = False

        # Populate triangular part of matrix directly
        ElCov=[]
        for El in Els:
            if El[:4] == ' COV':
                ElCov.append(El)

        if len(ElCov) == 7:
            _, CoV[0,0,0],CoV[0,0,1],CoV[0,0,2] = ElCov[0].split() # c11, c12, c13
            _, CoV[0,0,3],CoV[0,0,4],CoV[0,0,5] = ElCov[1].split() # c14, c15, c16
            _, CoV[0,1,1],CoV[0,1,2],CoV[0,1,3] = ElCov[2].split() # c22, c23, c24
            _, CoV[0,1,4],CoV[0,1,5],CoV[0,2,2] = ElCov[3].split() # c25, c26, c33
            _, CoV[0,2,3],CoV[0,2,4],CoV[0,2,5] = ElCov[4].split() # c34, c35, c36
            _, CoV[0,3,3],CoV[0,3,4],CoV[0,3,5] = ElCov[5].split() # c44, c45, c46
            _, CoV[0,4,4],CoV[0,4,5],CoV[0,5,5] = ElCov[6].split() # c55, c56, c66

            # Populate the symmetric part
            for i in range(1,6):
                for j in range(i):
                    CoV[0,i,j]=CoV[0,j,i]

            # Set boolean
            CoV_EXISTS = True
        return CoV_EXISTS, CoV






    # -----------------------------------
    # Functions to read nbody results
    # -----------------------------------
    #@timer
    def parse_nbody_txt( text_filepath ):
        '''
            Read a text file
            
            returns:
            --------
            name: str
             - name of object being processed: hopefully is an UNPACKED DESIGNATION
            array:
             - time, state, triangular-cov
        '''
        if os.path.isfile( text_filepath ):
            
            # Simple read of file -> array
            array =  np.loadtxt(text_filepath)

            # Get name
            with open(text_filepath) as fh:
                # remove '\n' from the end, and '# ' from the start
                name = fh.readline().strip()[2:]

        #      name  times,   states
        return name, array[:,0], array[:,1:]





    # -----------------------------------
    # Functions to create *FAKE DATA*
    # - Seems like this belongs in TEST DIRECTORY-FUNCTIONALITIES ...
    # -----------------------------------
    #@timer
    def create_nbody_txt( text_filepath ):
        '''
            Convenience function to create a text-file of coordinates
            ****This is only for the purposes of data exploration****
            
            #  I am constructing an array that is ~20,000 rows long
            # Each row stores 1-time, 6-posn&vel, and 21-CoVarCoeffs
            #
            #    ------ 28 ------>>>
            #   |         [[ t, x, y, z, vx, vy, vz, x_x, x_y, ... ]
            #   |         ...
            # 20,000      ...
            #   |         ...
            #   |
            #   V         ...
            #   V         ]
            #

        '''
        
        # times
        times = np.arange(40000, 60000, 1)
        
        # Set up empty array and put times into zeroth element
        a = np.zeros( ( len(times) , nFields)  )
        a[:,0] = times
        
        # Fake some slowly varying coordinate data
        r = np.random.uniform(low=2.5, high=3.5)
        p_days = r**(3./2.) * 365.25
        a[:,1] = r*np.cos( 2.*np.pi * times/p_days )
        a[:,2] = r*np.sin( 2.*np.pi * times/p_days )
        
        # Leave z == 0 to check whether zeros cause problems for cheby
        #a[3]
        
        # Populate vel & CoVar components
        for n, v in enumerate(['vx','vy','vz']):
            a[:,4+n]  = 0.1*n*a[:,1]
        for n, c in enumerate(covar_names):
            a[:,7+n]  = 1.e-9*n*a[:,1]

        # Save to file
        np.savetxt(text_filepath , a , header='2022 AA' )





# End
