# -*- coding: utf-8 -*-
# mpc_nbody/mpc_nbody/parse_input.py

'''
----------------------------------------------------------------------------
mpc_nbody's module for parsing OrbFit + ele220 elements

Mar 2020
Mike Alexandersen & Matthew Payne & Matthew Holman

This module provides functionalities to
(a) read an OrbFit .fel/.eq file with heliocentric ecliptic cartesian els
(b) read ele220 element strings
(c) convert the above to barycentric equatorial cartesian elements

This is meant to prepare the elements for input into the n-body integrator
----------------------------------------------------------------------------
'''

# Import third-party packages
# -----------------------------------------------------------------------------
import os, sys
import numpy as np
from astropy.time import Time
import getpass

if getpass.getuser() in ['matthewjohnpayne']:  # Payne's dev laptop set up differently ...:
    sys.path.append('/Users/matthewjohnpayne/Envs/mpcvenv/')
import mpcpp.MPC_library as mpc

# Import neighbouring packages
# -----------------------------------------------------------------------------

# Default for caching stuff using lru_cache
# -----------------------------------------------------------------------------

# Constants and stuff
# -----------------------------------------------------------------------------
DATA_PATH = os.path.realpath(os.path.dirname(__file__))
au_km = 149597870.700  # This is now a definition

# Data classes/methods
# -----------------------------------------------------------------------------


class ParseElements():
    '''
    Class for parsing elements and returning them in the format required by
    NbodySim in mpc_nbody
    
    Can be instantiated empty
    
    Can be instantiated with "input"
    
    input :
    -------
    string or list-of-strings
     - each string must be either:
       --- valid file-path
       --- valid file-contents (e.g. orbfit output)
     
    '''

    def __init__(self, input=None, filetype='eq', save_parsed=False , save_file='save_file.tmp', CHECK_EPOCHS=True):
    
        # The variables that will be used to hold the elements
        # - They get populated by *parse_orbfit* & *make_bary_equatorial*
        self.helio_ecl_vec_EXISTS   = False
        self.helio_ecl_vec          = None
        self.helio_ecl_cov_EXISTS   = False
        self.helio_ecl_cov          = None
        self.bary_eq_vec_EXISTS     = False
        self.bary_eq_vec            = None
        self.bary_eq_cov_EXISTS     = False
        self.bary_eq_cov            = None
        self.time                   = None
        
    
        
        # If input provided, parse it:
        try :
            assert input is not None
            
            # Interpret input & read any files that need to be read
            list_of_file_contents = self._parse_input_type(input)
            
            # Extract the required info from each of the file-contents
            for file_contents in list_of_file_contents :
            
                # Option to process ele220 files is not fully written
                if filetype == 'ele220':
                    self.parse_ele220(file_contents, CHECK_EPOCHS=CHECK_EPOCHS)
                    
                # The assumed standard input will be orbfit files ...
                if (filetype == 'fel') | (filetype == 'eq'):
                    self.parse_orbfit(file_contents, CHECK_EPOCHS=CHECK_EPOCHS)
                
            # Convert all input coords to Barycentric Equatorial
            self.make_bary_equatorial()
            
            # Optionally save the parsed data to file
            if save_parsed:
                self.save_elements(save_file=save_file)
                
        # If no input provided / it couldn't be parsed, instantiate empty object
        except:
            print("instantiating empty object.")
                  
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
        if isinstance(input, str): input = list(str)
        
        # Check we have a list
        assert isinstance(input, list)
        
        # If we have filepaths, get contents
        # otherwise, assume these already are data similar to file-contents
        contents = []
        for item in input:
            if os.path.isfile(item):
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
        for n,coeff in enumerate(self.bary_eq_vec):
            suffix = '\n' if n in [2,5] else ''
            outfile.write(f"{coeff: 18.15e} " + suffix)

    def parse_ele220(self, ele220file_contents, CHECK_EPOCHS=True ):
        '''
        Parse a file containing a single ele220 line.
        Currently returns junk data.
        NOT ACTUALLY IMPLEMENTED YET!!!
        '''
        # make fake data & set appropriate variables
        self._get_and_set_junk_data()

    def parse_orbfit(self, felfile_contents, CHECK_EPOCHS=True):
        '''
        Parse a file containing OrbFit elements for a single object & epoch.

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
            local_helio_ecl_vec = np.array([float(car_x),   float(car_y),  float(car_z), \
                                            float(car_dx),  float(car_dy), float(car_dz)] )
            if self.helio_ecl_vec is None:
                self.helio_ecl_vec = np.array([local_helio_ecl_vec])
            else :
                self.helio_ecl_vec = np.stack((self.helio_ecl_vec , local_helio_ecl_vec))
                                                    
                                            
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
                self.Time = local_Time
            else:
                # If we are processing multiple files, then
                # check that all supplied epochs are the same
                if self.Time != local_Time:
                    raise TypeError("The epochs vary from file-to-file: "
                                    f"felfile_contents= {felfile_contents}")
                    
            # Parse carEls (the contents of the orbfit file) to get
            # the cartesian covariance matrix
            self.helio_ecl_cov_EXISTS, local_helio_ecl_cov = _parse_Covariance_List(carEls)
            if self.helio_ecl_cov is None:
                self.helio_ecl_cov = np.array([local_helio_ecl_cov])
            else :
                self.helio_ecl_cov = np.stack( (self.helio_ecl_cov, local_helio_ecl_cov ) )
            print('jgkkjgkljgljgljhl!!!')

        except Exception as e:
            print( "Supplied felfile_contents could not be processed. "
                    f"felfile_contents={felfile_contents}"
                    f"e = {e}")

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
        
        
    def _get_and_set_junk_data(self, BaryEqDirect=False ):
        """Just make some junk data for saving."""
        self.time                           = Time(2458849.5, format='jd', scale='tdb')
        v   = np.array( [3., 2., 1., 0.3, 0.2, 0.1] )
        CoV = 0.01 * np.ones((6,6))
        
        # Default is to make helio-ecl, then calc bary-eq from that
        if not BaryEqDirect:
            self.helio_ecl_vec              = v
            self.helio_ecl_vec_EXISTS       = True
            
            self.helio_ecl_cov              = CoV
            self.helio_ecl_cov_EXISTS       = True
        
            self.make_bary_equatorial()
            
        # Alternative is to directly set bary-eq
        else:
            self.bary_eq_vec                = v
            self.bary_eq_vec_EXISTS         = True
            
            self.bary_eq_cov                = CoV
            self.bary_eq_cov_EXISTS         = True



# Functions
# -----------------------------------------------------------------------------
    
def ecliptic_to_equatorial(input, backwards=False):
    '''
    Rotates a cartesian vector or Cov-Matrix from mean ecliptic to mean equatorial.
    
    Backwards=True converts backwards, from equatorial to ecliptic.
    
    inputs:
    -------
    input : 1-D or 2-D arrays
     - If 1-D, then len(input) must be 3 or 6
     - If 2-D, then input.shape must be (6,6)
     
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
    if   input.ndim == 1 and input.shape[0] in [3,6]:
        R      = R6 if input.shape[0] == 6 else R3
        output = R @ input
        
    # Matrix (CoV) input => R & R.T
    elif input.ndim == 2 and input.shape == (6,6):
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
        input_xyz - np.ndarray length 3 or 6
        backwards - boolean
    output:
        output_xyz  - np.ndarray
                    - same shape as input_xyz

    input_xyz MUST BE EQUATORIAL!!!
    '''
    direction = -1 if backwards else +1

    # Ensure we have an array of the correct shape to work with
    input_xyz = np.atleast_1d(input_xyz)
    assert input_xyz.ndim == 1
    assert input_xyz.shape[0] in [3,6]
    
    # Position & Motion of the barycenter w.r.t. the heliocenter (and vice-versa)
    delta, delta_vel = mpc.jpl_kernel[0, 10].compute_and_differentiate(jd_tdb)
    
    # Work out whether we need xyz or xyzuvw
    delta = delta if input_xyz.shape[0] == 3 else np.block([delta,delta_vel])
    
    # Shift vectors & return
    return input_xyz + delta * direction / au_km





def _old_parse_Covariance_List(Els):
    '''
    Convenience function for reading and splitting the covariance
    lines of an OrbFit file.
    Not intended for user usage.
    '''
    ElCov  = []
    covErr = ""
    for El in Els:
        if El[:4] == ' COV':
            ElCov.append(El)
    if len(ElCov) == 7:
        _, c11, c12, c13 = ElCov[0].split()
        _, c14, c15, c16 = ElCov[1].split()
        _, c22, c23, c24 = ElCov[2].split()
        _, c25, c26, c33 = ElCov[3].split()
        _, c34, c35, c36 = ElCov[4].split()
        _, c44, c45, c46 = ElCov[5].split()
        _, c55, c56, c66 = ElCov[6].split()
    if len(ElCov) != 7:
        c11, c12, c13, c14, c15, c16, c22 = "", "", "", "", "", "", ""
        c23, c24, c25, c26, c33, c34, c35 = "", "", "", "", "", "", ""
        c36, c44, c45, c46, c55, c56, c66 = "", "", "", "", "", "", ""
        covErr = ' Empty covariance Matrix for '
    return (covErr, c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26,
            c33, c34, c35, c36, c44, c45, c46, c55, c56, c66)
    
def _parse_Covariance_List(Els):
    '''
    Convenience function for reading and splitting the covariance
    lines of an OrbFit file.
    Not intended for user usage.
    # MJP : 20200901 : Suggest to just make & return the required matrix
    '''
    # Set-up array of zeroes
    CoV        = np.zeros( (6,6) )
    CoV_EXISTS = False
    
    # Populate triangular part of matrix directly
    ElCov=[]
    for El in Els:
        if El[:4] == ' COV':
            ElCov.append(El)
    if len(ElCov) == 7:
        _, CoV[0,0],CoV[0,1],CoV[0,2] = ElCov[0].split() # c11, c12, c13
        _, CoV[0,3],CoV[0,4],CoV[0,5] = ElCov[1].split() # c14, c15, c16
        _, CoV[1,1],CoV[1,2],CoV[1,3] = ElCov[2].split() # c22, c23, c24
        _, CoV[1,4],CoV[1,5],CoV[2,2] = ElCov[3].split() # c25, c26, c33
        _, CoV[2,3],CoV[2,4],CoV[2,5] = ElCov[4].split() # c34, c35, c36
        _, CoV[3,3],CoV[3,4],CoV[3,5] = ElCov[5].split() # c44, c45, c46
        _, CoV[4,4],CoV[4,5],CoV[5,5] = ElCov[6].split() # c55, c56, c66
        
        # Populate the symmetric part
        for i in range(1,6):
            for j in range(i):
                print(f'Setting Cov[{i,j}] = CoV{[j,i]}')
                CoV[i,j]=CoV[j,i]
                
        # Set boolean
        CoV_EXISTS = True
    return CoV_EXISTS, CoV
 
