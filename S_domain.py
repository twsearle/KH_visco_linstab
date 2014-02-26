#----------------------------------------------------------------------------#
#   Linear Stability Analysis of velocity profile in 2D 
#   Last modified: Mon 22 Jul 2013 12:10:19 BST
#----------------------------------------------------------------------------#

"""Using the fully spectral method, find the stability of a 2D velocity 
   profile"""

#MODULES
import sys
import time
from scipy import *
from scipy import linalg
import cPickle as pickle
import ConfigParser

#FUNCTIONS

def mkDiffY():
    """Makes a matrix to differentiate array wrt to y"""

    y = YPOINTS 
    # The C function:
    C = ones(M, dtype='d')
    C[0] = 2.0
    C[M-1] = 2.0

    # Set up the differentiation matrix
    Dy = zeros((M, M) , dtype= 'd')
    Dy[0, 0] = (2.0*(M-1)**2 + 1) / 6.0
    Dy[M-1, M-1] = -(2.0*(M-1)**2 + 1) /  6.0
    for n in range(M):
        for m in range(M):
            if n != m:
                Dy[n, m] =  C[n] / C[m]
                Dy[n, m] = Dy[n, m] * pow(-1, n + m) / (y[n] - y[m]) 		 
            if n == m and n != 0 and n != (M-1):
                Dy[n, m] = -y[m] / (2 * (1 - pow(y[m],2)) )
    return Dy

def mkEqnMat():
    """make the matrix holding the linear stability equations."""
    equation_mat = zeros((6*M,6*M), dtype='complex')

    ####################### x direction       ############################
    #du
    equation_mat[0:M,0:M] = -Re*1.j*kx*MMU + beta*(-(kx**2)*eye(M,M) + MDYY)
    #dv
    equation_mat[0:M,M:2*M] = -Re*MMDYU
    #dp
    equation_mat[0:M,2*M:3*M] = -1.j*kx*eye(M,M)
    #dtxx
    equation_mat[0:M,3*M:4*M] = (1.-beta)*1.j*kx*eye(M,M)
    #dtyy
    equation_mat[0:M,4*M:5*M] = 0 
    #dtxy
    equation_mat[0:M,5*M:6*M] = (1.-beta)*MDY
    ####################### y direction       ############################
    #du
    equation_mat[M:2*M,0:M] = 0
    #dv
    equation_mat[M:2*M,M:2*M] = -Re*1.j*kx*MMU + beta*(-(kx**2)*eye(M,M) + MDYY)
    #dp
    equation_mat[M:2*M,2*M:3*M] = -MDY
    #dtxx
    equation_mat[M:2*M,3*M:4*M] = 0 
    #dtyy
    equation_mat[M:2*M,4*M:5*M] = (1.-beta)*MDY
    #dtxy
    equation_mat[M:2*M,5*M:6*M] = (1.-beta)*1.j*kx*eye(M,M)
    ####################### incompressibility ############################
    #du
    equation_mat[2*M:3*M,0:M] = 1.j*kx*eye(M,M)
    #dv
    equation_mat[2*M:3*M,M:2*M] = MDY
    #dp
    equation_mat[2*M:3*M,2*M:3*M] = 0
    #dtxx
    equation_mat[2*M:3*M,3*M:4*M] = 0
    #dtyy
    equation_mat[2*M:3*M,4*M:5*M] = 0
    #dtxy
    equation_mat[2*M:3*M,5*M:6*M] = 0
    ####################### xx stress         ############################
    #du
    equation_mat[3*M:4*M,0:M]     = - Weiss*2.j*kx*MMTXX - Weiss*2*dot(MMTXY,MDY)\
                                    - 2.j*kx*eye(M,M)
    #dv
    equation_mat[3*M:4*M,M:2*M]   = Weiss*diagflat(dot(MDY,tauxx))
    #dp
    equation_mat[3*M:4*M,2*M:3*M] = 0
    #dtxx
    equation_mat[3*M:4*M,3*M:4*M] = Weiss*1.j*kx*MMU + eye(M,M)
    #dtyy
    equation_mat[3*M:4*M,4*M:5*M] = 0
    #dtxy
    equation_mat[3*M:4*M,5*M:6*M] = - Weiss*2*MMDYU
    ####################### yy stress         ############################
    #du
    equation_mat[4*M:5*M,0:M]     = 0 
    #dv
    equation_mat[4*M:5*M,M:2*M]   = - Weiss*2.j*kx*MMTXY - 2*MDY
    #dp
    equation_mat[4*M:5*M,2*M:3*M] = 0
    #dtxx
    equation_mat[4*M:5*M,3*M:4*M] = 0
    #dtyy
    equation_mat[4*M:5*M,4*M:5*M] = Weiss*1.j*kx*MMU + eye(M,M)
    #dtxy
    equation_mat[4*M:5*M,5*M:6*M] = 0
    ####################### xy stress         ############################
    #du
    equation_mat[5*M:6*M,0:M]     = - MDY
    #dv
    equation_mat[5*M:6*M,M:2*M]   = Weiss*diagflat(dot(MDY,tauxy)) \
                                   - Weiss*1.j*kx*MMTXX - 1.j*kx*eye(M,M)
    #dp
    equation_mat[5*M:6*M,2*M:3*M] = 0 
    #dtxx
    equation_mat[5*M:6*M,3*M:4*M] = 0
    #dtyy
    equation_mat[5*M:6*M,4*M:5*M] = - Weiss*MMDYU
    #dtxy
    equation_mat[5*M:6*M,5*M:6*M] = Weiss*1.j*kx*MMU + eye(M,M)

    #Apply BC's to equation matrix:
    equation_mat[0,:] = zeros(6*M)
    equation_mat[0,0] = 1
    equation_mat[M-1,:] = zeros(6*M) 
    equation_mat[M-1,M-1] = 1

    #Apply BC for v in equation matrix
    equation_mat[M,:] = zeros(6*M)
    equation_mat[M,M] = 1
    equation_mat[2*M-1,:] = zeros(6*M)
    equation_mat[2*M-1,2*M-1] = 1
    
    return equation_mat

#MAIN
config = ConfigParser.RawConfigParser()
fp = open('2D_flow_settings.cfg')
config.readfp(fp)
M = config.getint('settings', 'M')
Re = config.getfloat('settings', 'Re')
beta = config.getfloat('settings','beta')
Weiss = config.getfloat('settings','Weiss')
#kx = config.getfloat('settings', 'kx')
DELTA = config.getfloat('settings', 'DELTA')

fp.close()

ksettings = r_[0.22]

base_filename = '-M{M}-Re{Re}-beta{beta}-Wi{Weiss}-delta{delta}.pickle'.format(\
                  M=M,Re=Re,beta=beta,Weiss=Weiss, delta=DELTA)

element_number = r_[0:M]
YPOINTS = cos(pi*element_number/(M-1))

################make the velocity and stress profile###################
U0 = 1

#The length of the system is 2 so L = 1 in normalisation
vel_function = lambda y: U0*tanh(y/DELTA) / tanh(1/DELTA)

U = vel_function(YPOINTS)
#enforce bc's
U[0] = U0
U[M-1] = -U0

MDY = mkDiffY()

diffU = dot(MDY, U)

#Calculate stresses:

tauxx = Weiss*2*dot(MDY,U)*dot(MDY,U)
tauyy = zeros(M, dtype='d')
tauxy = dot(MDY,U)
#######################################################################

leading_eigs = zeros((len(ksettings),3))

for kindx, kx in enumerate(ksettings):
    print 'kx: {kx}'.format(kx=kx)
    print 'kindx: {kindx}'.format(kindx=kindx)
    #diag flat constructs a diagonal matrix with the vector on the diagonal
    MDY = mkDiffY()
    MDYY = dot(MDY,MDY)
    MMDYU = diagflat(dot(MDY,U))
    MMU   = diagflat(U)
    MMTXX = diagflat(tauxx)
    MMTXY = diagflat(tauxy)

    eqn_mat = mkEqnMat()

    #Make RHS matrix
    RHS = zeros((6*M,6*M), dtype= 'd')

    RHS[:2*M, :2*M] = Re*eye(2*M,2*M)
    RHS[2*M:3*M, 2*M:3*M] = 0
    RHS[3*M:, 3*M:] = -Weiss*eye(3*M,3*M)

    # apply bc's for small perturbation equation on RHS
    RHS[0,0] = 0
    RHS[M-1,M-1] = 0

    RHS[M,M] = 0
    RHS[2*M-1,2*M-1] = 0

    # Solve for the eigenvalues
    eigenvals = linalg.eig(eqn_mat, RHS, right=False, overwrite_a=True, overwrite_b=True)

    # Save output
    eigarray = vstack((real(eigenvals), imag(eigenvals))).T

    #pickle.dump((eigenvals,eigvecs), open('evecs-kx{kx}{fn}'.format(kx=kx, fn=base_filename), 'w'))

    #remove nans and infs from eigenvalues
    eigarray = eigarray[~isnan(eigarray).any(1), :]
    eigarray = eigarray[~isinf(eigarray).any(1), :]

    savetxt('ev-kx{kx}{fn}.dat'.format(kx=kx, fn=base_filename[:-7]), eigarray)

    large_eigs = zeros((len(eigarray[:,0]),2))
    for i in range(len(eigarray[:,0])):
        if (eigarray[i,0] >= 0) and (eigarray[i,0]<50):
            large_eigs[i,:] = eigarray[i,:]

    #print large_eigs
    selected_row = argmax(large_eigs[:,0])
    leading_eigs[kindx,0] =  kx
    leading_eigs[kindx,1:] = large_eigs[selected_row,:]
    
#savetxt('lead-ev{fn}.dat'.format(fn=base_filename[:-7]), leading_eigs)
