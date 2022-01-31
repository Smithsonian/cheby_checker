'''
----------------------------------------------------------------------------
coco.py
 - Coordinate Conversions
 
2020/21/22
Mike Alexandersen & Matthew Payne & Matthew Holman

This module provides functionalities to perform coordinate conversions / rotations / shifts of frame
----------------------------------------------------------------------------
'''

# Import third-party packages
# -----------------------------------------------------------------------------
import sys
import os
import numpy as np

# Import neighbouring packages
# -----------------------------------------------------------------------------
# cheby_checker/                 # <<-- repo
# cheby_checker/cheby_checker    # <<-- python
# cheby_checker/tests            # <<-- tests
this_dir = os.path.abspath(os.path.dirname( __file__ ))
repo_dir = os.path.abspath(os.path.dirname( this_dir ))
test_dir = os.path.join(repo_dir, 'tests')
code_dir = os.path.join(repo_dir, 'cheby_checker')

# old conversion library
import MPC_library as mpc

from cheby_checker import Base


# Coordinate Conversion Functions
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
    return input_xyz + delta * direction / Base().au_km

