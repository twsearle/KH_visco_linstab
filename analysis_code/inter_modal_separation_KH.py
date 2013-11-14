# -----------------------------------------------------------------------------
#
#   Graph the Inter modal separation of an eigenvalue spectrum
#   Last modified: Thu 14 Nov 15:18:44 2013
#
# -----------------------------------------------------------------------------
"""
Chapter 7 of Boyd's book Chebyshev and Fourier spectral methods suggests that a
good way of checking eigenvalue convergence is by comparing the ordinal
difference between adjacent eigenvalues. Should be able to use this to filter
out the physical ones.

So far I have been unable to get this to work properly. I have struggled to
correlate the eigenvalues with peaks in either difference.

"""

# MODULES
from scipy import *
from matplotlib import pylab as plt
from matplotlib import rc
import ConfigParser

# PARAMETERS ------------------------------------------------------------------


k = 7.0

Re = 10000
Weiss = 0.001
beta = 0.1

delta = 0.1

M1Half = 100
M2Half = 200

args = {'k': k, 'M': M1Half, 'Re': Re, 'b': beta, 'Wi': Weiss, 'delta': delta}
filename1 = 'Alex-ev-kx{k}-M{M}-Re{Re}-beta{b}-Wi{Wi}-delta{delta}.dat'.format(**args)
args['M'] = M2Half
filename2 = 'Alex-ev-kx{k}-M{M}-Re{Re}-beta{b}-Wi{Wi}-delta{delta}.dat'.format(**args)

M1 = 2*M1Half
M2 = 2*M2Half

# -----------------------------------------------------------------------------

# FUNCTIONS

def calc_sigma(evs):
    """
    Calculates the scaling factor to calculate the ordinal and nearest
    differences.

    """

    sigma = zeros(vecLen, dtype='d')
    sigma[0] = absolute(evs[0] - evs[1])

    print len(evs), vecLen
    for j in range(1, vecLen-1):
        term1 = absolute(evs[j] - evs[j-1])
        term2 = absolute(evs[j+1] - evs[j])
        sigma[j] = 0.5 * (term1 + term2)
        if sigma[j] is 0.0:
            print "Degenerate mode at: ", j
    del j
    return sigma

def ordinal_difference(evs, evs2):
    """
    Finds an array of ordinal differences given two spectra. 

    do_j = |lambda_j(N1) - lambda_j(N2)| / sigma_j

    Only works if eigenvalues are evenly spaced. I think this means that it
    only works if you can guarantee that the eigenvalues in the different
    spectra correspond to each other. We can't guarantee this.

    """
    
    sigma = calc_sigma(evs)
 
    # print sigma
    # ordinal difference
    ordDiff = zeros(vecLen, dtype='d')
    for j in range(vecLen):
        ordDiff[j] = absolute(evs[j] - evs2[j]) / sigma[j]
    del j
    return ordDiff

def nearest_difference(evs1, evs2):
    """
    Find the distance between the two nearest eigenvalues.

    dn_j = min |lambda_j(N1) - lamda_i(N2)| / sigma_j
          -----
          k in [1, N2]

    """

    sigma = calc_sigma(evs1)
    nearestDiff = zeros((vecLen), dtype='d')
    for j in range(vecLen):
        minimum = infty
        for i in range(vecLen2):
            diff = absolute(evs1[j] - evs2[i]) / sigma[j]
            if diff < minimum:
                minimum = diff
        del i
        nearestDiff[j] = minimum
    del j

    return nearestDiff


# MAIN 


evs1 = genfromtxt(filename1)
evs1 = evs1[:,0] + 1.j*evs1[:,1]
evs2 = genfromtxt(filename2)
evs2 = evs2[:,0] + 1.j*evs2[:,1]

# vecLen set by the lower resolution eigenvalues
vecLen = len(evs1)
vecLen2 = len(evs2)
#print evs1

ordDiff = ordinal_difference(evs1, evs2)
nearestDiff = nearest_difference(evs1, evs2)

#make plots prettier:
inches_per_Lx = 1.4
inches_per_Ly = 2.2
fig_width =  8
fig_height = 4*inches_per_Ly      
fig_size =  [fig_width,fig_height]
rc('figure', figsize=fig_size)

savetxt('ordDiff.dat', ordDiff)
savetxt('nearestDiff.dat', nearestDiff)

plt.semilogy(r_[0:vecLen], ordDiff, 'rx')
plt.semilogy(r_[0:vecLen], nearestDiff, 'b.')
plt.show()
