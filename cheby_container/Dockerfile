#
# Container to ensure that all of the rebound/reboundx stuff works as cheby_checker.nbody is critically dependent on this
# There is an earlier / simpler version of this container in the mpc_dev branch of the reboundx/examples/ephem_forces repo
# MJP: 2021-07-16
#
#
FROM ubuntu:18.04
MAINTAINER Matt Payne @ MPC <matthewjohnpayne@gmail.com>


# --- SYSTEM INSTALLS -----------
RUN apt-get -y update && apt-get install -y 	gcc  		\
    						make 		\ 
    						git		\ 	
						python3.6 	\
						python3-pip	\
						perl		\
						wget 		\
						bzip2		\
						vim		\
						sudo		

 
# --- REBOUND INSTALLS ----------
#RUN git clone https://github.com/hannorein/rebound.git
RUN git clone https://github.com/Smithsonian/rebound.git
WORKDIR /rebound
RUN make
ENV REB_DIR=/rebound
ENV LD_LIBRARY_PATH=/rebound/src

WORKDIR /
RUN git clone -b mpc_dev https://github.com/Smithsonian/reboundx.git
WORKDIR /reboundx/examples/ephem_forces/
RUN make
ENV REBX_DIR=/reboundx
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/reboundx/src
ENV JPL_PLANET_EPHEM=/reboundx/examples/ephem_forces/linux_p1550p2650.440
ENV JPL_SB_EPHEM=/reboundx/examples/ephem_forces/sb441-n16.bsp

# Get the ephemeris file in to the container (goes into the ephem_forces sub-dir of reboundx). 
# *** OBVIOUSLY THIS SHOULD BE PART OF REBOUUNDX REPO, BUT NO LFS ON PUBLIC FORKS!?!?!? ***
#RUN wget ftp://ssd.jpl.nasa.gov/pub/eph/small_bodies/asteroids_de441/sb441-n16.bsp
COPY sb441-n16.bsp /reboundx/examples/ephem_forces/sb441-n16.bsp


# --- MPC_ORB INSTALL (from GitHub) --------
# This is used to facilitate the reading of mpc_orb.json format orbits data
WORKDIR /
RUN git clone https://github.com/matthewjohnpayne/mpc_orb.git
WORKDIR mpc_orb
RUN python3.6 -m pip install -e .


# --- PYTHON INSTALLS -----------
# Anaconda installing
# TODO: Use requirements.txt
RUN python3.6 -m pip install numpy pytest astropy astropy_healpix astroquery scipy novas novas_de405 jplephem
# TODO: Install mpcpp? https://bitbucket.org/mpcdev/mpcpp/src/master/

# --- CHEBY_CHECKER INSTALLS ----
# While in development, will bind a local version (see build_container.py)
#WORKDIR /
#RUN git clone https://github.com/matthewjohnpayne/cheby_checker.git
WORKDIR /cheby_checker
#RUN python3.6 -m pip install -e .


# --- Run deployment code on entry :
COPY pythonServer.py /pythonServer.py
ENTRYPOINT ["python3.6", "/pythonServer.py"]




