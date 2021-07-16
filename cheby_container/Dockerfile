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
						python3.8 	\
						python3-pip	\
						perl		\
						wget 		\
						bzip2		\
						vim		\
						sudo		

# --- PYTHON INSTALLS -----------
# Anaconda installing 
RUN python3.6 -m pip install numpy 
RUN python3.8 -m pip install numpy pytest astropy astropy_healpix astroquery scipy novas novas_de405 jplephem 	
 
# --- REBOUND INSTALLS ----------
RUN git clone https://github.com/hannorein/rebound.git
WORKDIR /rebound
RUN make
ENV REB_DIR=/rebound
ENV LD_LIBRARY_PATH=/rebound/src

WORKDIR /
RUN git clone -b mpc_dev https://github.com/matthewjohnpayne/reboundx.git
WORKDIR /reboundx/examples/ephem_forces/
RUN perl -pi -e 's/OPENGL=1/OPENGL=0/g' Makefile
RUN make
ENV REBX_DIR=/reboundx
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/reboundx/src

RUN wget ftp://ssd.jpl.nasa.gov/pub/eph/small_bodies/asteroids_de441/sb441-n16.bsp


# --- CHEBY_CHECKER INSTALLS ----
WORKDIR /
RUN git clone https://github.com/matthewjohnpayne/cheby_checker.git



# --- Run deployment code on entry :
COPY pythonServer.py /pythonServer.py
ENTRYPOINT ["python3.8", "/pythonServer.py"]




