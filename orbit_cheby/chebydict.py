import os
import subprocess
import json
import math
import numpy as np


def error(a, b):

    # get cheb polynomial approximation error at fitted points
    
    err = 0
    for i in range(len(a)):
        err = max(err, abs(a[i]-b[i]))
    return err


def getcheb(x,y,minorder,maxorder,maxerr):

    '''
    Get lowest order sufficiently accurate Chebyshev polynomial fit. 

    Inputs:
    x,y -- lists of values to which to fit Chebyshev polynomial
    minorder,maxorder -- lowest and highest orders of polynomials to consider fitting
    maxerr -- maximum amount by which an acceptable fitted polynomial can be off 
    '''
    
    maxorder = min([maxorder,len(x)])
    order = (maxorder-minorder)/2 + minorder
    chebCandidate = np.polynomial.chebyshev.chebfit(x, y, math.ceil(order))
    quickEval = np.polynomial.chebyshev.chebval(x, chebCandidate)

    if error(quickEval, y) <= maxerr:
        if math.floor(order) == minorder:
            return chebCandidate
        else:
            return getcheb(x,y,minorder,round(order),maxerr)
    else:
        if math.ceil(order) == maxorder:
            return chebCandidate
        else:
            return getcheb(x,y,round(order),maxorder,maxerr)


def datediv(desig_up,coords,covs,mindate,maxdate,divsize,coordorder=['x','y','z','vx','vy','vz'],coordnames=['x','y','z','vx','vy','vz'],coventries = ['00','01','02','11','12','22'],minorder=5,maxorder=25,maxerr=.000001,coverrfrac=.01):

    '''
    Generates dictionaries of Chebyshev coefficients approximating the desired coordnames and covariance values over intervals of size divsize.

    Inputs:
    desig_up -- unpacked designation of object being fit
    coords -- list of lists of cartesian coordinates; each sublist has as entries MJD time and list of position/velocity coords in order given by coordorder
    covs -- list of lists of cartesian covariances; each sublist has as entries MJD time and 6 by 6 covariance matrix
    mindate,maxdate -- earliest and latest MJD dates for which you want chebyshev approximations
    divsize -- length in days of time interval over which to fit a single chebyshev polynomial
    coordnames -- list of coords for which you want chebyshev approximations
    coventries -- list of covariance matrix entries for which you want chebyshev approximations
    minorder,maxorder -- lowest and highest orders of polynomials to consider fitting
    maxerr -- maximum amount (in AU) by which an acceptable fitted polynomial can be off for cartesian coordinates
    coverrfrac -- fractional error by which an acceptable fitted polynomial can be off for covariance values

    Outputs:
    Subdirectory [desig_up]jsondir in current working directory containing JSON-format dictionaries of chebyshev coefficients, one file for each divsize-sized time interval. Files are named [desig_up]_[MJD start date]_[MJD end date].json .
    '''
    
    numdivs = math.ceil((maxdate-mindate)/divsize)

    # create results directory if it doesn't exist yet
    if not os.path.exists(desig_up+'jsondir'):
        p = subprocess.Popen('mkdir '+desig_up+'jsondir',shell=True)
        p.communicate()

    for ind in range(numdivs):

        # collect coord/cov data for the date interval currently being fitted
        mindatenow = mindate + ind*divsize
        maxdatenow = min(maxdate,(ind+1)*divsize + mindate)
        coordsnow = [i for i in coords if (i[0]>=mindatenow and i[0]<=maxdatenow)]
        covsnow = [i for i in covs if (i[0]>=mindatenow and i[0]<=maxdatenow)]
        times = [i[0] for i in coordsnow]

        maxorder = min(maxorder,len(coordsnow))

        # initialize results dictionary, then generate desired chebyshev fits and store in dictionary
        dictnow = {}
        dictnow['name'] = desig_up
        dictnow['t_init'] = int(mindatenow)
        dictnow['t_final'] = int(maxdatenow)

        for name in coordnames:
            try:
                ind = [i for i,s in enumerate(coordorder) if s==name][0]
            except:
                print('coordinate '+name+' not in coordlist')
                continue
            var = [s[1][ind] for s in coordsnow]
            dictnow[name] = getcheb(times,var,minorder,maxorder,maxerr).tolist()

        for entry in coventries:
            covinds = (int(entry[0]),int(entry[1]))
            var = [s[1][covinds[0]][covinds[1]] for s in covsnow]
            errnow = min([abs(i) for i in var if i != 0]) * coverrfrac 
            dictnow['cov'+coordorder[covinds[0]]+coordorder[covinds[1]]] = getcheb(times,var,minorder,maxorder,errnow).tolist()

        # create json, write to file
        with open(desig_up+'jsondir/'+desig_up+'_'+str(mindatenow)+'_'+str(maxdatenow)+'.json','w') as fh:
            json.dump(dictnow,fh)
