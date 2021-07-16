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

# Import neighbouring packages
# -----------------------------------------------------------------------------
try:  # Import ephem_forces from whereever REBX_DIR is set to live
    sys.path.append(os.environ['REBX_DIR'])
    from examples.ephem_forces.ephem_forces import production_integration_function_wrapper
except (KeyError, ModuleNotFoundError):
    from reboundx.examples.ephem_forces.ephem_forces import production_integration_function_wrapper

import MPC_library as mpc

# Payne's dev laptop set up differently ...:
if getpass.getuser() in ['matthewjohnpayne']:
    sys.path.append('/Users/matthewjohnpayne/Envs/mpc_orb/mpc_orb/')
from parse import MPCORB



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





# Data classes/methods
# -----------------------------------------------------------------------------


class ParseElements():
    '''
    Class for parsing input elements and returning them in the format required by
    NbodySim class
    (see below for definition of NBodySim class)
    
    Can be instantiated empty
    
    Can be instantiated with "input"
    
    input :
    -------
    string or list-of-strings
     - each string must be either:
       --- valid file-path
       --- valid file-contents (e.g. orbfit output)
     
    '''

    def __init__(   self,
                    input=None,
                    filetype='json',
                    save_parsed=False ,
                    save_file='save_file.tmp',
                    CHECK_EPOCHS=True):
    
        # The variables that will be used to hold the elements
        # - They get populated by *parse_orbfit_felfile_txt* & *make_bary_equatorial*
        self.helio_ecl_vec_EXISTS   = False
        self.helio_ecl_vec          = None
        self.helio_ecl_cov_EXISTS   = False
        self.helio_ecl_cov          = None
        self.bary_eq_vec_EXISTS     = False
        self.bary_eq_vec            = None
        self.bary_eq_cov_EXISTS     = False
        self.bary_eq_cov            = None
        self.time                   = None
        self.non_grav_EXISTS        = False
        self.non_grav_array         = []
        
    
        print(f"input={input}")

        # If input provided, parse it:
        try :
            assert input is not None

            # Interpret input & read any files that need to be read
            list_of_file_contents = self._parse_input_type(input)
            
            # Extract the required info from each of the file-contents
            print("list_of_file_contents:", type(list_of_file_contents) , [type(_) for _ in list_of_file_contents] )
            for file_contents in list_of_file_contents :
            
                # Option to process orbfit-json files is not fully written
                if filetype == 'json':
                    self.parse_orbfit_json(file_contents, CHECK_EPOCHS=CHECK_EPOCHS)
                    
                # Option to process orbfit-eq0/eq1 fel-files
                if (filetype == 'fel') | (filetype == 'eq'):
                    self.parse_orbfit_felfile_txt(file_contents, CHECK_EPOCHS=CHECK_EPOCHS)
                
            # Convert all input coords to Barycentric Equatorial
            self.make_bary_equatorial()
            
            # Optionally save the parsed data to file
            if save_parsed:
                self.save_elements(save_file=save_file)
                
        # If no input provided / it couldn't be parsed, instantiate empty object
        except Exception as e:
            print("ParseElements.__init__: instantiating empty object.")
            print(f"Exception:{e.__str__()}")

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
            
    def parse_orbfit_json(self, jsonfile_dictionary, CHECK_EPOCHS=True ):
        '''
        Parse a file containing OrbFit elements for a single object & epoch.
        Assumes input is in standard mpcorb-json format
        '''
        
        # Instantiate MPCORB: NB This parses & validates ...
        M = MPCORB(jsonfile_dictionary)
        
        # ----------- TIME -------------------------
        # Using Astropy.time for time conversion,
        # because life's too short for timezones and time scales.
        epoch,timesystem    = M.epoch_data["epoch"],M.epoch_data["timesystem"]
        
        allowed_timesystems = {'TDT':'tt'}
        assert timesystem in allowed_timesystems
        local_Time      = Time(float(epoch), format='mjd', scale=allowed_timesystems[timesystem] )

        # In production, we likely want to be using "standard-epochs" ending in 00.
        # - This will allow multiple objects to be simultaneously integrated
        if CHECK_EPOCHS==True and "00." not in epoch:
            raise TypeError("Supplied file may not contain standard epochs:"
                            f"epoch= {epoch}")
        if self.time is None :
            self.time = local_Time
        else:
            # If we are processing multiple files, then
            # check that all supplied epochs are the same
            if self.time != local_Time:
                raise TypeError("The epochs vary from file-to-file: "
                                f"jsonfile_dictionary= {jsonfile_dictionary}")

        # ----------- ELEMENTS/VARIABLES -----------
        # Get Cartesian dict out of the overall jsonfile_dictionary
        cart_dict = M.CAR

        # Parse the elements and stack with any others extracted from previous files
        # - The elements are always of length 6
        numparams                   = cart_dict['numparams']
        element_array               = np.atleast_2d( cart_dict['element_array'] ).reshape(1,6)
        self.helio_ecl_vec_EXISTS   = True if numparams else False
        if self.helio_ecl_vec is None:
            self.helio_ecl_vec = element_array
        else :
            self.helio_ecl_vec = np.append(self.helio_ecl_vec , element_array, axis=0)
                                                

        
        # ----------- COVARIANCE  -------------------
        # Parse jsonfile_dictionary to get the cartesian covariance matrix
        # The cov matrix will be of dimension N^2, with N>=6, as the non-gravs are incluuded
        local_helio_ecl_cov         = np.atleast_3d(cart_dict['covariance_array']).reshape(1,numparams,numparams )
        self.helio_ecl_cov_EXISTS   = True #if local_helio_ecl_cov else False
        if self.helio_ecl_cov is None:
            self.helio_ecl_cov = local_helio_ecl_cov
        else :
            self.helio_ecl_cov = np.append( self.helio_ecl_cov, local_helio_ecl_cov  )

                
        # ----------- NONGRAVS -----------------------
        self.non_grav_array.append( M.nongrav_data )
        self.non_grav_EXISTS = np.any( _["non_gravs"] for _ in self.non_grav_array )

                

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
            self.bary_eq_vec   = equatorial_helio2bary(
                                    ecliptic_to_equatorial(self.helio_ecl_vec),
                                    self.time.tdb.jd
                                )
            # Set boolean as well (not sure if we'll really use these ...)
            self.bary_eq_vec_EXISTS = True

        if self.helio_ecl_cov_EXISTS:
            # Only need to do a rotation for the CoV
            self.bary_eq_cov = ecliptic_to_equatorial(self.helio_ecl_cov)
        
            # Set booleans as well (not sure if we'll really use these ...)
            self.bary_eq_cov_EXISTS = True

        if not self.helio_ecl_vec_EXISTS and not self.helio_ecl_cov_EXISTS:
            raise TypeError("There does not seem to be any valid helio_ecl to transform into bary_eq")
            
        return True
        
        



# Functions
# -----------------------------------------------------------------------------
    
def ecliptic_to_equatorial(input, backwards=False):
    '''
    Rotates a cartesian vector or Cov-Matrix from mean ecliptic to mean equatorial.
    
    Backwards=True converts backwards, from equatorial to ecliptic.
    
    inputs:
    -------
    input : 2-D or 3-D arrays
     - If 2-D, then input.shape[1] must be in [3,6]
     - If 3-D, then input.shape[1,2] must be  (6,6)
     
    NB: The inputs are 2D & 3D (rather than 1 or 2) so
    that we can "stack" lots of 1D vectors, or 2D Covariance-matricees
     
    output:
    -------
    output : np.ndarray
     - same shape as input
    '''

    # Ensure we have an array
    input = np.atleast_1d(input)
    
    # The rotation matricees we may use
    direction = -1 if backwards else +1
    R3 = mpc.rotate_matrix(mpc.Constants.ecl * direction)
    R6 = np.block( [ [R3, np.zeros((3,3))],[np.zeros((3,3)),R3] ])
    
    # Vector input => Single rotation operation
    # NB ... self.helio_ecl_vec.ndim ==2, self.helio_ecl_vec.shape = (N,6)
    if   input.ndim == 2 and input.shape[1] in [3,6]:
        R      = R6 if input.shape[1] == 6 else R3
        output = R.dot(input.T).T
        
    # Matrix (CoV) input => R & R.T
    # NB: CoV.ndim ==3 , CoV.shape == (N,6,6)
    elif input.ndim == 3 and input.shape == (1,6,6):
        R = R6
        output = R @ input @ R.T

    # Unknown input
    else:
        sys.exit(f'Does not compute: input.ndim=={input.ndim} , input.shape={input.shape}')

    assert output.shape == input.shape
    return output


def equatorial_helio2bary(input_xyz, jd_tdb, backwards=False):
    '''
    Convert from heliocentric to barycentic cartesian coordinates.
    backwards=True converts backwards, from bary to helio.
    input:
        input_xyz - np.ndarray of shape (X,3) or (X,6)
        backwards - boolean
    output:
        output_xyz  - np.ndarray
                    - same shape as input_xyz

    input_xyz MUST BE EQUATORIAL!!!
    '''
    direction = -1 if backwards else +1
    
    # Ensure we have an array of the correct shape to work with
    assert input_xyz.ndim == 2, f" input_xyz={input_xyz}\n input_xyz.shape={input_xyz.shape}\n input_xyz.ndim={input_xyz.ndim}"
    assert input_xyz.shape[1] in [3,6]
    
    # Position & Motion of the barycenter w.r.t. the heliocenter (and vice-versa)
    delta, delta_vel = mpc.jpl_kernel[0, 10].compute_and_differentiate(jd_tdb)
    
    # Work out whether we need xyz or xyzuvw
    delta = delta.reshape((1,3)) if input_xyz.shape[1] == 3 else np.block([delta,delta_vel]).reshape((1,6))

    # Shift vectors & return
    return input_xyz + delta * direction / au_km





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





'''
#
#
#
#
#
#
#
#
#
#
#
#
#
#
'''







class NbodySim():
    '''
    Class for containing all of the N-body related stuff.
    '''

    def __init__(   self,
                    input       =   None,
                    filetype    =   'eq',
                    save_parsed =   False,
                    CHECK_EPOCHS=   True):

        # Expected class variables
        self.PE                 = None # This will be ParseElements object ...
        self.geocentric         = False
        self.input_vectors      = None
        self.input_n_particles  = None
        self.output_times       = None
        self.output_vectors     = None
        self.output_n_times     = None
        self.output_n_particles = None
        self.time_parameters    = None

        #If input filename provided, process it using ParseElements :
        if input:
            self.PE= ParseElements( input       = input,
                                    filetype    = filetype,
                                    save_parsed = save_parsed,
                                    CHECK_EPOCHS= CHECK_EPOCHS)


    def __call__(self,
                tstart      =0,
                tstep       =20,
                trange      =600,
                epoch       =None,
                vectors     =None,
                covariances =None,
                save_output =None,
                verbose     =False):
        '''
        Allows the NbodySim object to be called as a function
         - Performs integration by calling *integration_function* within *run_nbody*
        '''
        self.time_parameters = tstart, tstep, trange
        
        # Get the initial conditions for the particle integration
        if vectors is None or epoch is None :
            try:
                epoch       = self.PE.time.tdb.jd
                vectors     = self.PE.bary_eq_vec
                covariances = self.PE.bary_eq_cov
                nongrav_arr = self.PE.non_grav_array
            except AttributeError as e :
                print(f"Invalid ParseElements ? : {e}")

        if vectors is None or epoch is None :
            raise TypeError("If you didn't parse a particle from an input "
                                "file, you must supply 'vectors' and 'epoch'.")

        # Call the routine that runs the nbody integration
        (self.epoch,
        self.input_vectors,
        self.input_covariances,
        self.input_n_particles,
        self.output_times,
        self.output_vectors,
        self.output_covariances,
        self.output_n_times,
        self.output_n_particles) = self.run_nbody(  epoch,
                                                    vectors,
                                                    covariances,
                                                    nongrav_arr,
                                                    tstart,
                                                    tstep,
                                                    trange,
                                                    geocentric  = self.geocentric,
                                                    verbose     = verbose)
        print(f'###!!!{type(self.output_times):}!!!###' if verbose else '')
        
        # Save output if requested
        if save_output is not None:
            if isinstance(save_output, str):
                self.save_output(output_file=save_output)
            else:
                self.save_output()


    # --------------------------------------------
    # Functions to run integration & parse results
    # --------------------------------------------

    def run_nbody(  self,
                    epoch,
                    input_states,
                    input_covariances,
                    tstart,
                    tend,
                    geocentric=False,
                    verbose=False):
        '''
        Run the nbody integrator with the parsed input.

        Input:
        ------
        epoch :          float
         - Common epoch at which all input vectors are defined
        input_states :   np.ndarray
         - 2D, shape = (N_particles, 6)
         - elements (xyzuvw) for N_particles to be integrated
        input_covariances :   np.ndarray or None
         - 3D, shape = (N_particles, 6,6)
        tstart = float,
         - Julian Date at start of integration.
        tstep = float or integer,
         - major time step of integrator.
        trange = float or integer,
         - rough total time of integration.
        geocentric = boolean,
         - use geo- (True) or heliocentric (False)

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
        '''
        
        # Get the number of particles
        assert input_states.ndim == 2 and input_states.shape[1] == 6, \
            f"input_states.ndim={input_states.ndim},input_states.shape={input_states.shape}"
        n_particles = input_states.shape[0]
        
        # Now run the nbody integrator:
        outtime, states, partial_derivatives_wrt_state, partial_derivatives_wrt_NG, return_value = \
            production_integration_function_wrapper(    tstart,
                                                        tend,
                                                        epoch,
                                                        geocentric,
                                                        n_particles,
                                                        input_states,
                                                        non_grav_dict_list = None,
                                                        epsilon = 1e-8,
                                                        tstep_max = 32.)
                                                        
        # Reshape the partial derivative arrays
        partial_derivatives_wrt_state, partial_derivatives_wrt_NG = self.reshape_partial_deriv_arrays(  partial_derivatives_wrt_state,
                                                                                                        partial_derivatives_wrt_NG)
        
        
        
        # Calculate the covariance matrix (at each timestep) from the \partial X / \partial X_0 data
        if input_covariances is not None:
            final_covariance_arrays = self._get_covariance_from_tangent_vectors(input_covariances,
                                                                                partial_derivatives_wrt_state,
                                                                                partial_derivatives_wrt_NG = partial_derivatives_wrt_NG )
        else:
            final_covariance_arrays = None
            
        return( epoch,
                input_states,
                input_covariances,
                n_particles,
                outtime,
                states,
                final_covariance_arrays)

    def reshape_partial_deriv_arrays( self,  partial_derivatives_wrt_state,  partial_derivatives_wrt_NG):
        '''
        (1) partial_derivatives_wrt_state
        From: partial_derivatives_wrt_state.shape = (N_times, 6*n_particlea, 6)
        To  : partial_derivatives_wrt_state.shape = (N_times, n_particlea, 6, 6)
        
        (2) partial_derivatives_wrt_NG
        *** NOT YET IMPLEMENTED partial_derivatives_wrt_NG ***
        '''
        
        # (1) (N_times, 6*n_particlea, 6) -> (N_times, n_particlea, 6, 6)
        partial_derivatives_wrt_state = partial_derivatives_wrt_state.reshape(partial_derivatives_wrt_state.shape[0],-1,6,6)
        
        # (2)
        if partial_derivatives_wrt_NG is not None:
            raise Error('Have not coded partial_derivatives_wrt_NG into _get_covariance_from_tangent_vectors ')

        return partial_derivatives_wrt_state , partial_derivatives_wrt_NG
        
        
    def _get_covariance_from_tangent_vectors(self, init_covariances, partial_derivatives_wrt_state ,  partial_derivatives_wrt_NG=None ):
        '''
        Follow Milani et al 1999
        Gamma_t = [partial X / partial X_0] Gamma_0 [partial X / partial X_0]^T
        '''
        
        if partial_derivatives_wrt_NG is not None:
            raise Error('Have not coded partial_derivatives_wrt_NG into _get_covariance_from_tangent_vectors ')

        assert init_covariances.ndim == 3
        assert partial_derivatives_wrt_state.ndim == 4
        assert partial_derivatives_wrt_NG is None or partial_derivatives_wrt_NG.ndim  == 4
        #print(f'partial_derivatives_wrt_state.shape = {partial_derivatives_wrt_state.shape}')

        # Take the inverse of the covariance matrix to get the normal matrix
        # NB: We make a stack of identical matricees to use in the matrix multiplication below
        Gamma0      = np.linalg.inv(init_covariances)
        print(f'Gamma0.shape = {Gamma0.shape}')
        GammaStack0 = np.tile( Gamma0, (partial_derivatives_wrt_state.shape[0],1,1,1) )
        #print(f'GammaStack0.shape = {GammaStack0.shape}')

        # We need each of the individual pd arrays to be individually transposed
        # NB, tuple fixes dimensions 0 & 1 , but indicates that dimensions 2 & 3 will be swapped/transposed
        pds_transposed  = partial_derivatives_wrt_state.transpose( (0,1,3,2) )
        #print(f'pds_transposed.shape = {pds_transposed.shape}')

        # Do matrix multiplication: using the pd's to get the CoV as a func of time
        # NB matmul/@ automagically knows how to work on a stack of matricees
        # - https://stackoverflow.com/questions/34142485/difference-between-numpy-dot-and-python-3-5-matrix-multiplication
        GammaStack_t    = pds_transposed @ GammaStack0 @ partial_derivatives_wrt_state
        #print(f'GammaStack_t.shape = {GammaStack_t.shape}')

        # Magically, np.linalg.inv also knows how to deal with a stack of arrays/matrices
        # - https://stackoverflow.com/questions/11972102/is-there-a-way-to-efficiently-invert-an-array-of-matrices-with-numpy
        CoV_t           = np.linalg.inv( GammaStack_t )
        #print(f'CoV_t.shape = {CoV_t.shape}')


        return CoV_t
            

    # -----------------------------------
    # Misc Funcs ...
    # -----------------------------------

    def save_output(self, output_file='simulation_states.dat'):
        """
        Save all the outputs to file.

        Inputs:
        -------
        output_file : string, filename to write elements to.

        The file is overwritten if it already exists.
        """
        outfile = open(output_file, 'w')
        outfile.write(f'#Input vectors: [')
        for coo in self.input_vectors:
            outfile.write(f'{coo} ')
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





 


'''
#
#
#
#
#
#
#
#
#
#
#
#
#
#
'''



# Functions to read nbody results
# --------------------------------------------------------------

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





# Functions to create *FAKE DATA*
# --------------------------------------------------------------


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
