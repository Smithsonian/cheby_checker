'''
MJP: 2021-11-08

Code to perform a position-check / residual calculation 

Given a set of observations of a specified object, 
(i) use ephemeris look-up to calculate the expected observation positions as a function of time
(ii) compare simulations to actual observations and calculate residuals

N.B. We exlicitly want to separate out the position-check / residual calculation from 
     the later mpchecker functionality
'''


