# -*- coding: utf-8 -*-
# cheby_checker/cheby_checker/obs_pos.py

'''
    --------------------------------------------------------------
    Parse an obscode and figure out where the observer is at a given time.

    Jan 2020
    Matt Payne & Mike Alexandersen

    *WRITE MORE STUFF*

    --------------------------------------------------------------
    '''

# Import third-party packages
# --------------------------------------------------------------
import sys
import getpass
import numpy as np
from astropy import units as u

# Import neighboring packages
# --------------------------------------------------------------
# Different machines set up differently ...
# ... adding paths to force stuff to work while developing
if getpass.getuser() in ['matthewjohnpayne']:
    sys.path.append('/Users/matthewjohnpayne/Envs/mpcvenv/')
else:
    pass
from mpcpp import MPC_library as mpc


class ObsPos():
    '''Class containing all the functionality for figuring out where
       the observer is. '''

    def __init__(self):
        # MPC Observatory list:
        self.obsCodes = mpc.Observatory()

    def get_heliocentric_equatorial_xyz(self, jd_utc, obsCode=None,
                                        verbose=False):
        '''
            Get the heliocentric EQUATORIAL vector coordinates of the
            observatory at the time jd_utc.
        '''
        obsCode = self.check_obsCode(obsCode, verbose)
        return self.obsCodes.getObservatoryPosition(obsCode, jd_utc,
                                                    old=False) * u.au

    def get_heliocentric_ecliptic_xyz(self, jd_utc, obsCode=None,
                                      verbose=False):
        '''
            Get the heliocentric ECLIPTIC vector coordinates of the
            observatory at the time jd_utc.
        '''
        obsCode = self.check_obsCode(obsCode, verbose)
        helio_OBS_equ = self.get_heliocentric_equatorial_xyz(jd_utc, obsCode)
        helio_OBS_ecl = self.equatorial_to_ecliptic(helio_OBS_equ)
        return helio_OBS_ecl

    def equatorial_to_ecliptic(self, input_xyz, backwards=False):
        '''
            Convert an cartesian vector from mean equatorial to mean ecliptic.
            backwards=True converts backwards, from ecliptic to equatorial.
            input:
                input_xyz - np.array length 3
                backwards - boolean
            output:
                output_xyz - np.array length 3
        '''
        direction = -1 if backwards else +1
        rotation_matrix = mpc.rotate_matrix(-mpc.Constants.ecl * direction)
        output_xyz = np.dot(rotation_matrix, input_xyz.reshape(-1, 1)).flatten()
        return output_xyz

    def check_obsCode(self, obsCode=None, verbose=False):
        '''
            Check whether a valid Observatory Code has been supplied.
            If None, use 500 (Geocentre).
        '''
        if obsCode is None:
            return '500'
        if obsCode in ['XXX', '']:  # Observations with no ObsCode
            print('Bad ObsCode. Will use geocenter.\n' if verbose else '', end='')
            return '500'				# pretend Geocentre.
        if isinstance(obsCode, int):
            obsCode = str(obsCode)
        if len(obsCode) != 3:
            raise NotImplementedError("Bad Observatory Code!\n"
                                      "Observatory Code "
                                      "must be a three character string!\nFor "
                                      "four character Observatory Codes, please "
                                      "bug M. Alexandersen or M. Payne.")
        return obsCode

#End
