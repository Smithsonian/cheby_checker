'''
Tests of the "convenience_Horizons" routines that are used for testing.

Some of these tests really function as demos/documentation to
remind myself/ourselves of how these Horizons functions are
intended to work
'''

# Import standard packages
# --------------------------------------------------------------
import numpy as np

# Import the module to be tested
# --------------------------------------------------------------
import convenience_Horizons as Horizons



def test_read_Horizons_state_from_text():
    '''
    This is NOT testing a built-in Horizons function
    This is just testing a little convenience routine created by MJP
    This convenience routine is ONLY used as part of the testing code for Cheby Checker
    '''
    
    # input text
    lines = """
X =-2.590350154796811E+00 Y =-7.949342693459856E-02 Z = 1.245107691757731E-01
VX=-1.454708370733871E-03 VY=-9.503445860627428E-03 VZ=-3.846514535533382E-03
""".split('\n')[1:-1]
    
    # use the target function to extract the coordinaets
    result = Horizons.read_Horizons_state_from_text( lines )

    # check that the results are as expected
    expected_array = np.array([ float('-2.590350154796811E+00'), float('-7.949342693459856E-02'), float('1.245107691757731E-01'),
                                float('-1.454708370733871E-03'), float('-9.503445860627428E-03'), float('-3.846514535533382E-03') ] )

    assert np.allclose(expected_array, result, rtol=1e-08, atol=1e-08)


def test_extract_first_state_from_text():
    '''
    This is NOT testing a built-in Horizons function
    This is just testing a little convenience routine created by MJP
    This convenience routine is ONLY used as part of the testing code for Cheby Checker
    '''
    # input text
    lines = """
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
""".split('\n')[1:-1]

    #print(lines)
    # use the target function to extract the coordinaets
    result = Horizons.extract_first_state_from_text( lines )

    # check that the results are as expected
    expected_array = np.array([ float('-2.590350154796811E+00'), float('-7.949342693459856E-02'), float('1.245107691757731E-01'),
                                float('-1.454708370733871E-03'), float('-9.503445860627428E-03'), float('-3.846514535533382E-03') ] )

    assert np.allclose(expected_array, result, rtol=1e-08, atol=1e-08)





def test_nice_Horizons_A():
    '''
    Testing Mike A's convenience wrapper around Horizon query functionality
     - Much of this test is being done to provide some reminder
       to myself/ourselves as to how to use the Horizons tool
    
    Deliberately *not* using all of the functionalities of pytest here.
    Just want to keep it simple and keep it obvious what everything is supposed to be doing.
    
    Here we extract the
        HELIOCENTRIC state
    for
        Asteroid number 12345 (== 1993 FT8)
    in an
        EQUATORIAL FRAME (refplane='earth')
    
    '''
    
    # Define the variables that will be used in the query
    target  = '12345'   # <<-- Asteroid number 12345 == 1993    FT8
    centre  = '500@10'
    epochs  = '2458850.0'
    id_type = 'smallbody'
    refplane= 'earth'

    # Hardpaste the expected results from a by-hand query of horizons
    hardpasted_results = """
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
""".split('\n')[1:-1]

    # Extract the hardpasted results into an array
    # (using convenience func, "extract_first_state_from_text", tested above)
    expected_array = Horizons.extract_first_state_from_text( hardpasted_results )

    # Call the nice_Horizons function (i.e. the focus of the test)
    result = Horizons.nice_Horizons(target, centre, epochs, id_type, refplane=refplane )
    
    # Check that the results are as expected
    assert np.allclose(expected_array, result, rtol=1e-11, atol=1e-11)






def test_nice_Horizons_B():
    '''
    Testing Mike A's convenience wrapper around Horizon query functionality
     - Much of this test is being done to provide some reminder
       to myself/ourselves as to how to use the Horizons tool
    
    Deliberately *not* using all of the functionalities of pytest here.
    Just want to keep it simple and keep it obvious what everything is supposed to be doing.
    
    Here we extract the
        HELIOCENTRIC state
    for
        GEOCENTER
    in an
        EQUATORIAL FRAME (refplane='earth')
    
    '''
    
    # Define the variables that will be used in the query
    target  = '399'       # Earth
    centre  = '500@10'
    epochs  = '2458850.0'
    id_type = 'majorbody'
    refplane= 'earth'

    # Hardpaste the expected results from a by-hand query of horizons
    hardpasted_results = """
*******************************************************************************
 Revised: April 12, 2021                 Earth                              399
 
 GEOPHYSICAL PROPERTIES (revised Aug 15, 2018):
  Vol. Mean Radius (km)    = 6371.01+-0.02   Mass x10^24 (kg)= 5.97219+-0.0006
  Equ. radius, km          = 6378.137        Mass layers:
  Polar axis, km           = 6356.752          Atmos         = 5.1   x 10^18 kg
  Flattening               = 1/298.257223563   oceans        = 1.4   x 10^21 kg
  Density, g/cm^3          = 5.51              crust         = 2.6   x 10^22 kg
  J2 (IERS 2010)           = 0.00108262545     mantle        = 4.043 x 10^24 kg
  g_p, m/s^2  (polar)      = 9.8321863685      outer core    = 1.835 x 10^24 kg
  g_e, m/s^2  (equatorial) = 9.7803267715      inner core    = 9.675 x 10^22 kg
  g_o, m/s^2               = 9.82022         Fluid core rad  = 3480 km
  GM, km^3/s^2             = 398600.435436   Inner core rad  = 1215 km
  GM 1-sigma, km^3/s^2     =      0.0014     Escape velocity = 11.186 km/s
  Rot. Rate (rad/s)        = 0.00007292115   Surface area:
  Mean sidereal day, hr    = 23.9344695944     land          = 1.48 x 10^8 km
  Mean solar day 2000.0, s = 86400.002         sea           = 3.62 x 10^8 km
  Mean solar day 1820.0, s = 86400.0         Love no., k2    = 0.299
  Moment of inertia        = 0.3308          Atm. pressure   = 1.0 bar
  Mean temperature, K      = 270             Volume, km^3    = 1.08321 x 10^12
  Mean effect. IR temp, K  = 255             Magnetic moment = 0.61 gauss Rp^3
  Geometric albedo         = 0.367           Vis. mag. V(1,0)= -3.86
  Solar Constant (W/m^2)   = 1367.6 (mean), 1414 (perihelion), 1322 (aphelion)
 HELIOCENTRIC ORBIT CHARACTERISTICS:
  Obliquity to orbit, deg  = 23.4392911  Sidereal orb period  = 1.0000174 y
  Orbital speed, km/s      = 29.79       Sidereal orb period  = 365.25636 d
  Mean daily motion, deg/d = 0.9856474   Hill's sphere radius = 234.9
*******************************************************************************


*******************************************************************************
Ephemeris / WWW_USER Fri Jan 28 16:02:02 2022 Pasadena, USA      / Horizons
*******************************************************************************
Target body name: Earth (399)                     {source: DE441}
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
Output units    : AU-D
Output type     : GEOMETRIC cartesian states
Output format   : 3 (position, velocity, LT, range, range-rate)
Reference frame : ICRF
*******************************************************************************
JDTDB
   X     Y     Z
   VX    VY    VZ
   LT    RG    RR
*******************************************************************************
$$SOE
2458850.000000000 = A.D. 2020-Jan-01 12:00:00.0000 TDB [del_T=     69.183915 s]
 X =-1.749585912701602E-01 Y = 8.877645495087018E-01 Z = 3.848482875671789E-01
 VX=-1.721190438300784E-02 VY=-2.874039035670773E-03 VZ=-1.245648654352060E-03
 LT= 5.678966496273616E-03 RG= 9.832825679666131E-01 RR=-1.981645766688001E-05
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
""".split('\n')[1:-1]

    # Extract the hardpasted results into an array
    # (using convenience func, "extract_first_state_from_text", tested above)
    expected_array = Horizons.extract_first_state_from_text( hardpasted_results )

    # Call the nice_Horizons function (i.e. the focus of the test)
    result = Horizons.nice_Horizons(target, centre, epochs, id_type, refplane=refplane)
    print('result=\n' , result)
    # Check that the results are as expected
    assert np.allclose(expected_array, result, rtol=1e-11, atol=1e-11)






def test_nice_Horizons_C():
    '''
    Testing Mike A's convenience wrapper around Horizon query functionality
     - Much of this test is being done to provide some reminder
       to myself/ourselves as to how to use the Horizons tool
    
    Deliberately *not* using all of the functionalities of pytest here.
    Just want to keep it simple and keep it obvious what everything is supposed to be doing.
    
    Here we extract the
        TOPOCENTRIC (F51) state
    for
        Asteroid number 54321 (== 2000 JA81)
    in an
        EQUATORIAL FRAME (refplane='earth')
    

    '''
    
    # Define the variables that will be used in the query
    target  = '54321'       # <<-- Asteroid number 54321 == 2000 JA81
    centre  = 'F51'
    epochs  = '2458850.0'
    id_type = 'smallbody'
    refplane= 'earth'

    # Hardpaste the expected results from a by-hand query of horizons
    hardpasted_results = """
*******************************************************************************
JPL/HORIZONS                  54321 (2000 JA81)            2022-Jan-28 16:08:57
Rec #:   54321 (+COV) Soln.date: 2021-Oct-08_04:39:24   # obs: 1315 (1979-2021)
 
IAU76/J2000 helio. ecliptic osc. elements (au, days, deg., period=Julian yrs):
 
  EPOCH=  2456698.5 ! 2014-Feb-10.00 (TDB)         Residual RMS= .27282
   EC= .2508846058943067   QR= 1.938411174247326   TP= 2456093.1011463138
   OM= 91.32740861093403   W=  91.37096816741918   IN= 6.76912753748867
   A= 2.58760024089671     MA= 143.3507250314229   ADIST= 3.236789307546094
   PER= 4.1625             N= .236787236           ANGMOM= .026786318
   DAN= 2.43937            DDN= 2.41026            L= 182.707997
   B= 6.7671807            MOID= .95572901         TP= 2012-Jun-14.6011463138
 
Asteroid physical parameters (km, seconds, rotational period in hours):
   GM= n.a.                RAD= n.a.               ROTPER= n.a.
   H= 14.45                G= .150                 B-V= n.a.
                           ALBEDO= n.a.            STYP= n.a.
 
ASTEROID comments:
1: soln ref.= JPL#33, OCC=0
2: source=ORB
*******************************************************************************


*******************************************************************************
Ephemeris / WWW_USER Fri Jan 28 16:08:57 2022 Pasadena, USA      / Horizons
*******************************************************************************
Target body name: 54321 (2000 JA81)               {source: JPL#33}
Center body name: Earth (399)                     {source: DE441}
Center-site name: Pan-STARRS 1, Haleakala
*******************************************************************************
Start time      : A.D. 2020-Jan-01 12:00:00.0000 TDB
Stop  time      : A.D. 2020-Jan-01 12:00:00.0000 TDB
Step-size       : DISCRETE TIME-LIST
*******************************************************************************
Center geodetic : 203.744100,20.7071888,3.0763821 {E-lon(deg),Lat(deg),Alt(km)}
Center cylindric: 203.744100,5971.48324,2242.1878 {E-lon(deg),Dxy(km),Dz(km)}
Center pole/equ : ITRF93                          {East-longitude positive}
Center radii    : 6378.1 x 6378.1 x 6356.8 km     {Equator, meridian, pole}
Small perturbers: Yes                             {source: SB441-N16}
Output units    : AU-D
Output type     : GEOMETRIC cartesian states
Output format   : 3 (position, velocity, LT, range, range-rate)
EOP file        : eop.220127.p220422
EOP coverage    : DATA-BASED 1962-JAN-20 TO 2022-JAN-27. PREDICTS-> 2022-APR-21
Reference frame : ICRF
*******************************************************************************
Initial IAU76/J2000 heliocentric ecliptic osculating elements (au, days, deg.):
  EPOCH=  2456698.5 ! 2014-Feb-10.00 (TDB)         Residual RMS= .27282
   EC= .2508846058943067   QR= 1.938411174247326   TP= 2456093.1011463138
   OM= 91.32740861093403   W=  91.37096816741918   IN= 6.76912753748867
  Equivalent ICRF heliocentric cartesian coordinates (au, au/d):
   X= 2.934573285149345E+00  Y=-8.702901499041770E-01  Z=-7.535748078855007E-01
  VX= 3.948600090813408E-03 VY= 7.155151609877323E-03 VZ= 2.568700850735469E-03
Asteroid physical parameters (km, seconds, rotational period in hours):
   GM= n.a.                RAD= n.a.               ROTPER= n.a.
   H= 14.45                G= .150                 B-V= n.a.
                           ALBEDO= n.a.            STYP= n.a.
*******************************************************************************
JDTDB
   X     Y     Z
   VX    VY    VZ
   LT    RG    RR
*******************************************************************************
$$SOE
2458850.000000000 = A.D. 2020-Jan-01 12:00:00.0000 TDB [del_T=     69.183916 s]
 X = 3.140272938432556E-01 Y = 1.401450872643150E+00 Z = 5.824305212783573E-01
 VX= 6.595793009138184E-03 VY= 5.366428257622971E-04 VZ= 1.577496239642071E-03
 LT= 8.950941094980980E-03 RG= 1.549807407979245E+00 RR= 2.414570690380698E-03
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
""".split('\n')[1:-1]

    # Extract the hardpasted results into an array
    # (using convenience func, "extract_first_state_from_text", tested above)
    expected_array = Horizons.extract_first_state_from_text( hardpasted_results )

    # Call the nice_Horizons function (i.e. the focus of the test)
    result = Horizons.nice_Horizons(target, centre, epochs, id_type, refplane=refplane )
    print('result=\n' , result)
    # Check that the results are as expected
    assert np.allclose(expected_array, result, rtol=1e-11, atol=1e-11)






def test_nice_Horizons_D():
    '''
    Similar to test_nice_Horizons_C, but ECLIPTIC insted of equatorial
    
    Here we extract the
        TOPOCENTRIC (F51) state
    for
        Asteroid number 54321 (== 2000 JA81)
    in an
        ECLIPTIC FRAME (refplane='ecliptic')
    

    '''
    
    # Define the variables that will be used in the query
    target  = '54321'       # <<-- Asteroid number 54321 == 2000 JA81
    centre  = 'F51'
    epochs  = '2458850.0'
    id_type = 'smallbody'
    refplane= 'ecliptic'

    # Hardpaste the expected results from a by-hand query of horizons
    hardpasted_results = """
*******************************************************************************
JPL/HORIZONS                  54321 (2000 JA81)            2022-Jan-28 16:19:21
Rec #:   54321 (+COV) Soln.date: 2021-Oct-08_04:39:24   # obs: 1315 (1979-2021)
 
IAU76/J2000 helio. ecliptic osc. elements (au, days, deg., period=Julian yrs):
 
  EPOCH=  2456698.5 ! 2014-Feb-10.00 (TDB)         Residual RMS= .27282
   EC= .2508846058943067   QR= 1.938411174247326   TP= 2456093.1011463138
   OM= 91.32740861093403   W=  91.37096816741918   IN= 6.76912753748867
   A= 2.58760024089671     MA= 143.3507250314229   ADIST= 3.236789307546094
   PER= 4.1625             N= .236787236           ANGMOM= .026786318
   DAN= 2.43937            DDN= 2.41026            L= 182.707997
   B= 6.7671807            MOID= .95572901         TP= 2012-Jun-14.6011463138
 
Asteroid physical parameters (km, seconds, rotational period in hours):
   GM= n.a.                RAD= n.a.               ROTPER= n.a.
   H= 14.45                G= .150                 B-V= n.a.
                           ALBEDO= n.a.            STYP= n.a.
 
ASTEROID comments:
1: soln ref.= JPL#33, OCC=0
2: source=ORB
*******************************************************************************


*******************************************************************************
Ephemeris / WWW_USER Fri Jan 28 16:19:22 2022 Pasadena, USA      / Horizons
*******************************************************************************
Target body name: 54321 (2000 JA81)               {source: JPL#33}
Center body name: Earth (399)                     {source: DE441}
Center-site name: Pan-STARRS 1, Haleakala
*******************************************************************************
Start time      : A.D. 2020-Jan-01 12:00:00.0000 TDB
Stop  time      : A.D. 2020-Jan-01 12:00:00.0000 TDB
Step-size       : DISCRETE TIME-LIST
*******************************************************************************
Center geodetic : 203.744100,20.7071888,3.0763821 {E-lon(deg),Lat(deg),Alt(km)}
Center cylindric: 203.744100,5971.48324,2242.1878 {E-lon(deg),Dxy(km),Dz(km)}
Center pole/equ : ITRF93                          {East-longitude positive}
Center radii    : 6378.1 x 6378.1 x 6356.8 km     {Equator, meridian, pole}
Small perturbers: Yes                             {source: SB441-N16}
Output units    : AU-D
Output type     : GEOMETRIC cartesian states
Output format   : 3 (position, velocity, LT, range, range-rate)
EOP file        : eop.220127.p220422
EOP coverage    : DATA-BASED 1962-JAN-20 TO 2022-JAN-27. PREDICTS-> 2022-APR-21
Reference frame : Ecliptic of J2000.0
*******************************************************************************
Initial IAU76/J2000 heliocentric ecliptic osculating elements (au, days, deg.):
  EPOCH=  2456698.5 ! 2014-Feb-10.00 (TDB)         Residual RMS= .27282
   EC= .2508846058943067   QR= 1.938411174247326   TP= 2456093.1011463138
   OM= 91.32740861093403   W=  91.37096816741918   IN= 6.76912753748867
  Equivalent ICRF heliocentric cartesian coordinates (au, au/d):
   X= 2.934573285149345E+00  Y=-8.702901499041770E-01  Z=-7.535748078855007E-01
  VX= 3.948600090813408E-03 VY= 7.155151609877323E-03 VZ= 2.568700850735469E-03
Asteroid physical parameters (km, seconds, rotational period in hours):
   GM= n.a.                RAD= n.a.               ROTPER= n.a.
   H= 14.45                G= .150                 B-V= n.a.
                           ALBEDO= n.a.            STYP= n.a.
*******************************************************************************
JDTDB
   X     Y     Z
   VX    VY    VZ
   LT    RG    RR
*******************************************************************************
$$SOE
2458850.000000000 = A.D. 2020-Jan-01 12:00:00.0000 TDB [del_T=     69.183916 s]
 X = 3.140272938432556E-01 Y = 1.517483592803339E+00 Z =-2.309558662379520E-02
 VX= 6.595793009138184E-03 VY= 1.119852134073137E-03 VZ= 1.233860245870196E-03
 LT= 8.950941094980980E-03 RG= 1.549807407979244E+00 RR= 2.414570690380698E-03
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

  Ecliptic at the standard reference epoch

    Reference epoch: J2000.0
    X-Y plane: adopted Earth orbital plane at the reference epoch
               Note: IAU76 obliquity of 84381.448 arcseconds wrt ICRF X-Y plane
    X-axis   : ICRF
    Z-axis   : perpendicular to the X-Y plane in the directional (+ or -) sense
               of Earth's north pole at the reference epoch.

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
""".split('\n')[1:-1]

    # Extract the hardpasted results into an array
    # (using convenience func, "extract_first_state_from_text", tested above)
    expected_array = Horizons.extract_first_state_from_text( hardpasted_results )

    # Call the nice_Horizons function (i.e. the focus of the test)
    result = Horizons.nice_Horizons(target, centre, epochs, id_type, refplane=refplane)
    print('result=\n' , result)
    # Check that the results are as expected
    assert np.allclose(expected_array, result, rtol=1e-11, atol=1e-11)









def test_nice_Horizons_E():
    '''
    Here we use Horizons to get the Heliocentric EQUATORIAL position of the Observatory
    (NB we use a hack, setting the target as the Sun, and the center as the observatory)
    
    Here we extract the
        TOPOCENTRIC (F51) state
    for
        The Sun
    in an
        EQUATORIAL FRAME (refplane='earth')
    

    '''
    
    # Define the variables that will be used in the query
    target  = '10'
    centre  = 'F51'
    epochs  = '2458850.0'
    id_type = 'majorbody'
    refplane='earth'

    # Hardpaste the expected results from a by-hand query of horizons
    hardpasted_results = """
*******************************************************************************
 Revised: July 31, 2013                  Sun                                 10

 PHYSICAL PROPERTIES (updated 2018-Aug-15):
  GM, km^3/s^2          = 132712440041.93938  Mass, 10^24 kg        = ~1988500
  Vol. mean radius, km  = 695700              Volume, 10^12 km^3    = 1412000
  Solar radius (IAU)    = 696000 km           Mean density, g/cm^3  = 1.408
  Radius (photosphere)  = 696500 km           Angular diam at 1 AU  = 1919.3"
  Photosphere temp., K  = 6600 (bottom)       Photosphere temp., K  = 4400(top)
  Photospheric depth    = ~500 km             Chromospheric depth   = ~2500 km
  Flatness, f           = 0.00005             Adopted sid. rot. per.= 25.38 d
  Surface gravity       =  274.0 m/s^2        Escape speed, km/s    =  617.7
  Pole (RA,DEC), deg.   = (286.13, 63.87)     Obliquity to ecliptic = 7.25 deg.
  Solar constant (1 AU) = 1367.6 W/m^2        Luminosity, 10^24 J/s = 382.8
  Mass-energy conv rate = 4.260 x 10^9 kg/s   Effective temp, K     = 5772
  Sunspot cycle         = 11.4 yr             Cycle 24 sunspot min. = 2008 A.D.

  Motion relative to nearby stars = apex : R.A.= 271 deg.; DEC.= +30 deg.
                                    speed: 19.4 km/s (0.0112 au/day)
  Motion relative to 2.73K BB/CBR = apex : l= 264.7 +- 0.8; b= 48.2 +- 0.5 deg.
                                    speed: 369 +-11 km/s
*******************************************************************************


*******************************************************************************
Ephemeris / WWW_USER Fri Jan 28 16:31:17 2022 Pasadena, USA      / Horizons
*******************************************************************************
Target body name: Sun (10)                        {source: DE441}
Center body name: Earth (399)                     {source: DE441}
Center-site name: Pan-STARRS 1, Haleakala
*******************************************************************************
Start time      : A.D. 2020-Jan-01 12:00:00.0000 TDB
Stop  time      : A.D. 2020-Jan-01 12:00:00.0000 TDB
Step-size       : DISCRETE TIME-LIST
*******************************************************************************
Center geodetic : 203.744100,20.7071888,3.0763821 {E-lon(deg),Lat(deg),Alt(km)}
Center cylindric: 203.744100,5971.48324,2242.1878 {E-lon(deg),Dxy(km),Dz(km)}
Center pole/equ : ITRF93                          {East-longitude positive}
Center radii    : 6378.1 x 6378.1 x 6356.8 km     {Equator, meridian, pole}
Output units    : AU-D
Output type     : GEOMETRIC cartesian states
Output format   : 3 (position, velocity, LT, range, range-rate)
EOP file        : eop.220127.p220422
EOP coverage    : DATA-BASED 1962-JAN-20 TO 2022-JAN-27. PREDICTS-> 2022-APR-21
Reference frame : ICRF
*******************************************************************************
JDTDB
   X     Y     Z
   VX    VY    VZ
   LT    RG    RR
*******************************************************************************
$$SOE
2458850.000000000 = A.D. 2020-Jan-01 12:00:00.0000 TDB [del_T=     69.183916 s]
 X = 1.749807755042866E-01 Y =-8.877977147373203E-01 Z =-3.848633185157228E-01
 VX= 1.742085896638510E-02 VY= 3.013989047237847E-03 VZ= 1.245251026128347E-03
 LT= 5.679196211202660E-03 RG= 9.833223418736220E-01 RR=-1.085591308641175E-04
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
""".split('\n')[1:-1]

    # Extract the hardpasted results into an array
    # (using convenience func, "extract_first_state_from_text", tested above)
    expected_array = Horizons.extract_first_state_from_text( hardpasted_results )

    # Call the nice_Horizons function (i.e. the focus of the test)
    result = Horizons.nice_Horizons(target, centre, epochs, id_type, refplane=refplane)
    print('result=\n' , result)
    # Check that the results are as expected
    assert np.allclose(expected_array, result, rtol=1e-11, atol=1e-11)







def test_nice_Horizons_F():
    '''
    Similar to test_nice_Horizons_E, but ECLIPTIC instead of equatorial
    
    Here we use Horizons to get the Heliocentric ECLIPTIC position of the Observatory
        (NB we use a hack, setting the target as the Sun, and the center as the observatory)
    
    Here we extract the
        TOPOCENTRIC (F51) state
    for
        The Sun
    in an
        ECLIPTIC FRAME (refplane='ecliptic')
    

    '''
    
    # Define the variables that will be used in the query
    target  = '10'
    centre  = 'F51'
    epochs  = '2458850.0'
    id_type = 'majorbody'
    refplane='ecliptic'
    
    # Hardpaste the expected results from a by-hand query of horizons
    hardpasted_results = """
*******************************************************************************
 Revised: July 31, 2013                  Sun                                 10

 PHYSICAL PROPERTIES (updated 2018-Aug-15):
  GM, km^3/s^2          = 132712440041.93938  Mass, 10^24 kg        = ~1988500
  Vol. mean radius, km  = 695700              Volume, 10^12 km^3    = 1412000
  Solar radius (IAU)    = 696000 km           Mean density, g/cm^3  = 1.408
  Radius (photosphere)  = 696500 km           Angular diam at 1 AU  = 1919.3"
  Photosphere temp., K  = 6600 (bottom)       Photosphere temp., K  = 4400(top)
  Photospheric depth    = ~500 km             Chromospheric depth   = ~2500 km
  Flatness, f           = 0.00005             Adopted sid. rot. per.= 25.38 d
  Surface gravity       =  274.0 m/s^2        Escape speed, km/s    =  617.7
  Pole (RA,DEC), deg.   = (286.13, 63.87)     Obliquity to ecliptic = 7.25 deg.
  Solar constant (1 AU) = 1367.6 W/m^2        Luminosity, 10^24 J/s = 382.8
  Mass-energy conv rate = 4.260 x 10^9 kg/s   Effective temp, K     = 5772
  Sunspot cycle         = 11.4 yr             Cycle 24 sunspot min. = 2008 A.D.

  Motion relative to nearby stars = apex : R.A.= 271 deg.; DEC.= +30 deg.
                                    speed: 19.4 km/s (0.0112 au/day)
  Motion relative to 2.73K BB/CBR = apex : l= 264.7 +- 0.8; b= 48.2 +- 0.5 deg.
                                    speed: 369 +-11 km/s
*******************************************************************************


*******************************************************************************
Ephemeris / WWW_USER Fri Jan 28 16:43:25 2022 Pasadena, USA      / Horizons
*******************************************************************************
Target body name: Sun (10)                        {source: DE441}
Center body name: Earth (399)                     {source: DE441}
Center-site name: Pan-STARRS 1, Haleakala
*******************************************************************************
Start time      : A.D. 2020-Jan-01 12:00:00.0000 TDB
Stop  time      : A.D. 2020-Jan-01 12:00:00.0000 TDB
Step-size       : DISCRETE TIME-LIST
*******************************************************************************
Center geodetic : 203.744100,20.7071888,3.0763821 {E-lon(deg),Lat(deg),Alt(km)}
Center cylindric: 203.744100,5971.48324,2242.1878 {E-lon(deg),Dxy(km),Dz(km)}
Center pole/equ : ITRF93                          {East-longitude positive}
Center radii    : 6378.1 x 6378.1 x 6356.8 km     {Equator, meridian, pole}
Output units    : AU-D
Output type     : GEOMETRIC cartesian states
Output format   : 3 (position, velocity, LT, range, range-rate)
EOP file        : eop.220127.p220422
EOP coverage    : DATA-BASED 1962-JAN-20 TO 2022-JAN-27. PREDICTS-> 2022-APR-21
Reference frame : Ecliptic of J2000.0
*******************************************************************************
JDTDB
   X     Y     Z
   VX    VY    VZ
   LT    RG    RR
*******************************************************************************
$$SOE
2458850.000000000 = A.D. 2020-Jan-01 12:00:00.0000 TDB [del_T=     69.183916 s]
 X = 1.749807755042866E-01 Y =-9.676283142792063E-01 Z = 4.045892447000447E-05
 VX= 1.742085896638510E-02 VY= 3.260613297708340E-03 VZ=-5.640051197420877E-05
 LT= 5.679196211202660E-03 RG= 9.833223418736221E-01 RR=-1.085591308641179E-04
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

  Ecliptic at the standard reference epoch

    Reference epoch: J2000.0
    X-Y plane: adopted Earth orbital plane at the reference epoch
               Note: IAU76 obliquity of 84381.448 arcseconds wrt ICRF X-Y plane
    X-axis   : ICRF
    Z-axis   : perpendicular to the X-Y plane in the directional (+ or -) sense
               of Earth's north pole at the reference epoch.

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
""".split('\n')[1:-1]

    # Extract the hardpasted results into an array
    # (using convenience func, "extract_first_state_from_text", tested above)
    expected_array = Horizons.extract_first_state_from_text( hardpasted_results )

    # Call the nice_Horizons function (i.e. the focus of the test)
    result = Horizons.nice_Horizons(target, centre, epochs, id_type, refplane=refplane)
    print('result=\n' , result)
    # Check that the results are as expected
    assert np.allclose(expected_array, result, rtol=1e-11, atol=1e-11)






