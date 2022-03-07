# -*- coding: utf-8 -*-
# cheby_checker/cheby_checker/obs_pos.py

"""
--------------------------------------------------------------
Parse an obscode and figure out where the observer is at a given time.

Jan 2020
Matt Payne & Mike Alexandersen

*WRITE MORE STUFF*

--------------------------------------------------------------
"""

# Import third-party packages
# --------------------------------------------------------------
import numpy as np

# Import local modules
# --------------------------------------------------------------
from . import MPC_library as mpc


class ObsPos:
    """Class containing all the functionality for figuring out where
       the observer is. """

    def __init__(self):
        # MPC Observatory list:
        self.obsCodes = mpc.Observatory()

    def get_heliocentric_equatorial_xyz(self, jd_utc, obsCode=None,
                                        verbose=False):
        """
        Get the heliocentric EQUATORIAL vector coordinates of the
        observatory at the time jd_utc.
        """
        obsCode = self.check_obsCode(obsCode, verbose)
        return self.obsCodes.getObservatoryPosition(obsCode, jd_utc,
                                                    old=False) #* u.au

    def get_heliocentric_ecliptic_xyz(self, jd_utc, obsCode=None, verbose=False):
        """
        Get the heliocentric ECLIPTIC vector coordinates of the
        observatory at the time jd_utc.
        """
        obsCode = self.check_obsCode(obsCode, verbose)
        helio_OBS_equ = self.get_heliocentric_equatorial_xyz(jd_utc, obsCode)
        return self.equatorial_to_ecliptic(helio_OBS_equ)

    def equatorial_to_ecliptic(self, input_xyz, backwards=False):
        """
        Convert a cartesian vector from mean equatorial to mean ecliptic.
        backwards=True converts backwards, from ecliptic to equatorial.
        input:
            input_xyz - np.array length 3
            backwards - boolean
        output:
            output_xyz - np.array length 3
        """
        direction = -1 if backwards else +1
        rotation_matrix = mpc.rotate_matrix(-mpc.Constants.ecl * direction)
        return np.dot(rotation_matrix, input_xyz.reshape(-1, 1)).flatten()

    def check_obsCode(self, obsCode=None, verbose=False):
        """
        Check whether a valid Observatory Code has been supplied.
        If None, use 500 (Geocentre).
        """
        if obsCode is None:
            return '500'
        if obsCode in ['XXX', '']:  # Observations with no ObsCode
            print('Bad ObsCode. Will use geocenter.\n' if verbose else '', end='')
            return '500'

        if isinstance(obsCode, int):
            obsCode = str(obsCode)

        if len(obsCode) != 3:
            raise NotImplementedError(
                f"Bad Observatory Code!\n Observatory Code given: {obsCode}, must be a three digit number!" \
                "\nFor four digit Observatory Codes, please bug M. Alexandersen or M. Payne.")

        return obsCode
