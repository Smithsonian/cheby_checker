# sifter code-base goes in this directory

# Some notes ...
#
# Misc.
# - Do we need some way to track "names"/"unique-IDs" of each tracklet ?
#   --- In MPChecker I map names-to-integers
#   --- Might need trid+obscode+date to make unique
# - What (if anything) needs to be different between "back-fill" of 3e6 tracklets and the on-going processing of any tracklet added to the ITF ? 
#
# Base object ...
# - File / directory locations
#   --- While developing, let's just save locally to files on disk
#   --- Need to make future decision about optimal I/O method 
# - obs80 parsing methods?
# - date ranges of interest
# - healpix parameters 
# - need definitions of "precise match" 
#  
# Tracklet object ... 
# - all observations or just extremal ? 
# - healpix of extremal obs
# - integer days of extremal obs 
# - assume ~ 3e6 tracklets
# - assume ~1e4 days with data
# - implied only 3e2 per night 
# 
# Input orbit object ...
# - cheby representation of orbit: assume supplied as input ?
# - calculate healpix location on integer days
# - do we want "neighbors": probably 
# - what do we do about really fast things 
# - do we want to use velocity / RoM information ?
# - do we want to use uncertainty information (CoV-Matrix / Var-Orbits) ? 
#   --- can we inform healpix size/scale using uncertainties in some way ?
# - need to be able to query precalculations for (each) integer-day and integer-healpix
# - do we want to calculate RoM / angles / etc, to help with rough-cuts? 
# 
# Precalculations 
# - Use tracklet objects 
# - Calculate the "heapix" location of the first and last observation of each tracklet on a regular (daily?) interval, and then store the information by healpix (i.e. having a list of the objects present/close-to each healpix on each night).
# - Perhaps store a tracklet "against" all nearby/neighbouring tracklets ?
#  
#  
