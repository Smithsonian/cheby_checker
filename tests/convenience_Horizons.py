'''
Some convenience classes/functions for querying Horizons
Helps with tests of accuracy in various functions
'''

# import standard packages
# -----------------------------------------------------------------------------
import numpy as np

# import third-party packages
# -----------------------------------------------------------------------------
from astroquery.jplhorizons import Horizons


def nice_Horizons(target, centre, epochs, id_type, refplane='earth'):
    '''
    Mike Alexandersen
    Convenience function to reformat data returned by Horizons
    Only require the inputs I actually want to vary.
    Return in the format I actually want, not an astropy table.
    '''
    horizons_table  = Horizons(target, centre, epochs=epochs, id_type=id_type)
    horizons_vector = horizons_table.vectors(refplane=refplane)
    horizons_xyzv   = horizons_vector['x', 'y', 'z', 'vx', 'vy', 'vz']
    return np.array(list(horizons_xyzv.as_array()[0]))

def nice_Horizons_radec(target, centre, epochs, id_type, refplane='earth'):
    '''
    Convenience function to reformat data returned by Horizons
    Only require the inputs I actually want to vary.
    Return in the format I actually want, not an astropy table.
     - Returned shape == (N_epochs, 2)
    '''
    horizons_table  = Horizons(target, centre, epochs=epochs, id_type=id_type)
    horizons_eph    = horizons_table.ephemerides()
    return np.array( [  list(horizons_eph['RA'].data),
                        list(horizons_eph['DEC'].data) ] ).T
    
def nice_Horizons_LTT(target, centre, epochs, id_type, refplane='earth'):
    '''
    Mike Alexandersen
    Convenience function to reformat data returned by Horizons
    Only require the inputs I actually want to vary.
    Return in the format I actually want, not an astropy table.
    '''
    horizons_table  = Horizons(target, centre, epochs=epochs, id_type=id_type)
    horizons_vector = horizons_table.vectors(refplane=refplane)
    horizons_xyzv   = horizons_vector['lt']
    return np.array(list(horizons_xyzv.as_array()[0]))

def read_Horizons_state_from_text( two_lines):
    '''
    Extract these ...
     X =-2.590350154796811E+00 Y =-7.949342693459856E-02 Z = 1.245107691757731E-01
    VX=-1.454708370733871E-03 VY=-9.503445860627428E-03 VZ=-3.846514535533382E-03

    '''
    xyz_line = two_lines[0]
    uvw_line = two_lines[1]
    
    x = float(xyz_line.split('=')[1].strip('Y '))
    y = float(xyz_line.split('=')[2].strip('Z '))
    z = float(xyz_line.split('=')[3].strip())

    u = float(uvw_line.split('=')[1].strip('VY '))
    v = float(uvw_line.split('=')[2].strip('VZ '))
    w = float(uvw_line.split('=')[3].strip())
    
    return np.array( [x,y,z,u,v,w] )

def extract_first_state_from_text( text_block ):
    '''
    Extract lines that look like...
X =-2.590350154796811E+00 Y =-7.949342693459856E-02 Z = 1.245107691757731E-01
VX=-1.454708370733871E-03 VY=-9.503445860627428E-03 VZ=-3.846514535533382E-03

    from a big block that looks like ...
    
        *******************************************************************************
JPL/HORIZONS                  12345 (1993 FT8)             2022-Jan-28 14:39:42
Rec #:   12345 (+COV) Soln.date: 2021-Nov-10_08:38:58   # obs: 1959 (1993-2021)
 
IAU76/J2000 helio. ecliptic osc. elements (au, days, deg., period=Julian yrs):
 
  EPOCH=  2457108.5 ! 2015-Mar-27.00 (TDB)         Residual RMS= .2812
   EC= .1603033905689926   QR= 2.056207695854036   TP= 2457050.1973502915
   OM= 106.4549280993016   W=  314.1929318541605   IN= 3.350816780296945
   A= 2.448750742541829    MA= 14.99600220651154   ADIST= 2.841293789229623
   PER= 3.832              N= .25720961            ANGMOM= .02657056
   DAN= 2.14602            DDN= 2.68596            L= 60.6968709
   B= -2.401858            MOID= 1.06974006        TP= 2015-Jan-27.6973502915
 
Asteroid physical parameters (km, seconds, rotational period in hours):
   GM= n.a.                RAD= 1.506              ROTPER= n.a.
   H= 14.52                G= .150                 B-V= n.a.
                           ALBEDO= .407            STYP= n.a.
 
ASTEROID comments:
1: soln ref.= JPL#32, OCC=0
2: source=ORB
*******************************************************************************


*******************************************************************************
Ephemeris / WWW_USER Fri Jan 28 14:39:42 2022 Pasadena, USA      / Horizons
*******************************************************************************
Target body name: 12345 (1993 FT8)                {source: JPL#32}
Center body name: Sun (10)                        {source: DE441}
Center-site name: BODY CENTER
*******************************************************************************
Start time      : A.D. 2020-Jan-01 12:00:00.0000 TDB
Stop  time      : A.D. 2020-Jan-01 12:00:00.0000 TDB
Step-size       : DISCRETE TIME-LIST
*******************************************************************************
Center geodetic : 0.00000000,0.00000000,0.0000000 {E-lon(deg),Lat(deg),Alt(km)}
Center cylindric: 0.00000000,0.00000000,0.0000000 {E-lon(deg),Dxy(km),Dz(km)}
Center radii    : 696000.0 x 696000.0 x 696000.0 k{Equator, meridian, pole}
Small perturbers: Yes                             {source: SB441-N16}
Output units    : AU-D
Output type     : GEOMETRIC cartesian states
Output format   : 3 (position, velocity, LT, range, range-rate)
Reference frame : ICRF
*******************************************************************************
Initial IAU76/J2000 heliocentric ecliptic osculating elements (au, days, deg.):
  EPOCH=  2457108.5 ! 2015-Mar-27.00 (TDB)         Residual RMS= .2812
   EC= .1603033905689926   QR= 2.056207695854036   TP= 2457050.1973502915
   OM= 106.4549280993016   W=  314.1929318541605   IN= 3.350816780296945
  Equivalent ICRF heliocentric cartesian coordinates (au, au/d):
   X= 3.047919278950221E-01  Y= 1.902892265722551E+00  Z= 7.692605770652556E-01
  VX=-1.255238959074424E-02 VY= 2.052146789677108E-03 VZ= 1.612315394505861E-03
Asteroid physical parameters (km, seconds, rotational period in hours):
   GM= n.a.                RAD= 1.506              ROTPER= n.a.
   H= 14.52                G= .150                 B-V= n.a.
                           ALBEDO= .407            STYP= n.a.
*******************************************************************************
JDTDB
   X     Y     Z
   VX    VY    VZ
   LT    RG    RR
*******************************************************************************
$$SOE
2458850.000000000 = A.D. 2020-Jan-01 12:00:00.0000 TDB [del_T=     69.183915 s]
 X =-2.590350154796811E+00 Y =-7.949342693459856E-02 Z = 1.245107691757731E-01
 VX=-1.454708370733871E-03 VY=-9.503445860627428E-03 VZ=-3.846514535533382E-03
 LT= 1.498492268422344E-02 RG= 2.594558933811760E+00 RR= 1.558928955626413E-03
$$EOE
*******************************************************************************
 
TIME

  Barycentric Dynamical Time ("TDB" or T_eph) output was requested. This
continuous relativistic coordinate time is equivalent to the relativistic
proper time of a clock at rest in a reference frame comoving with the
solar system barycenter but outside the system's gravity well. It is the
independent variable in the solar system relativistic equations of motion.

  TDB runs at a uniform rate of one SI second per second and is independent
of irregularities in Earth's rotation.

  Calendar dates prior to 1582-Oct-15 are in the Julian calendar system.
Later calendar dates are in the Gregorian system.

REFERENCE FRAME AND COORDINATES

  International Celestial Reference Frame (ICRF)

    The ICRF is an adopted reference frame whose axes are defined relative to
    fixed extragalactic radio sources distributed across the sky.

    The ICRF was aligned with the prior FK5/J2000 dynamical system at the ~0.02
    arcsecond level but is not identical and has no associated standard epoch.

  Symbol meaning [1 au= 149597870.700 km, 1 day= 86400.0 s]:

    JDTDB    Julian Day Number, Barycentric Dynamical Time
    del_T    Time-scale conversion difference TDB - UT (s)
      X      X-component of position vector (au)
      Y      Y-component of position vector (au)
      Z      Z-component of position vector (au)
      VX     X-component of velocity vector (au/day)
      VY     Y-component of velocity vector (au/day)
      VZ     Z-component of velocity vector (au/day)
      LT     One-way down-leg Newtonian light-time (day)
      RG     Range; distance from coordinate center (au)
      RR     Range-rate; radial velocity wrt coord. center (au/day)

ABERRATIONS AND CORRECTIONS

 Geometric state vectors have NO corrections or aberrations applied.

Computations by ...

    Solar System Dynamics Group, Horizons On-Line Ephemeris System
    4800 Oak Grove Drive, Jet Propulsion Laboratory
    Pasadena, CA  91109   USA

    General site: https://ssd.jpl.nasa.gov/
    Mailing list: https://ssd.jpl.nasa.gov/email_list.html
    System news : https://ssd.jpl.nasa.gov/horizons/news.html
    User Guide  : https://ssd.jpl.nasa.gov/horizons/manual.html
    Connect     : browser        https://ssd.jpl.nasa.gov/horizons/app.html#/x
                  API            https://ssd-api.jpl.nasa.gov/doc/horizons.html
                  command-line   telnet ssd.jpl.nasa.gov 6775
                  e-mail/batch   https://ssd.jpl.nasa.gov/ftp/ssd/hrzn_batch.txt
                  scripts        https://ssd.jpl.nasa.gov/ftp/ssd/SCRIPTS
    Author      : Jon.D.Giorgini@jpl.nasa.gov
*******************************************************************************

    '''
    
    # find "$$SOE" in text_block
    SOE_index = [ n for n, line in enumerate(text_block) if "$$SOE" in line ][0]
    
    # no results ?
    if not SOE_index:
        return False

    # extract the two lines that contain the state information ...
    for XYZ_index,line in enumerate(text_block[SOE_index:]):
        if line.strip()[0]=='X':
            break

    state_lines = text_block[SOE_index:][XYZ_index:XYZ_index+2]
    
    # extract the coordinates from the lines
    xyzuvw = read_Horizons_state_from_text( state_lines )

    return xyzuvw
