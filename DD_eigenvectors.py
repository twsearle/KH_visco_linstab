##------------------------------------

from scipy import *
from scipy.linalg import inv, solve, eig, det
import cPickle as pickle
import sys


##------------------------------------

M = 150

Re = 0.0
Wi = 10.
bbeta = 0.1

_delta = 0.1
alpha = 0.28
##------------------------------------
parameter_string =  \
"""Parameters:       M = {M}
                     Re = {Re}
                     Wi = {Wi}
                     beta = {bbeta}
                     _delta ={_delta}""".format(M=M,Re=Re,Wi=Wi,bbeta=bbeta,
                                                _delta=_delta)
print parameter_string

II = identity(M,dtype='d')

cbar = ones(M,dtype='d')
cbar[0] = 2.0
cbar[M-1] = 2.0


ygl = zeros(M,dtype='d')
for m in range(M):
    ygl[m] = cos(pi*m/(M-1))

zzz_up = zeros(M,dtype='d')
for m in range(M):
    zzz_up[m] = 0.5*(ygl[m]+1.0)

zzz_down = zeros(M,dtype='d')
for m in range(M):
    zzz_down[m] = 0.5*(ygl[m]-1.0)


D1 = zeros((M,M),dtype='d')
for l in range(M):
    for j in range(M):
        if l != j:
            D1[l,j] = cbar[l]*((-1)**(l+j))/(cbar[j]*(ygl[l]-ygl[j]))

for j in range(1,M-1):
    D1[j,j] = -0.5*ygl[j]/(1.0-ygl[j]*ygl[j])

D1[0,0] = (2.0*(M-1)*(M-1)+1.0)/6.0
D1[M-1,M-1] = -D1[0,0]

D1 = 2*D1
D2 = dot(D1,D1)

## Laminar profile

U0up = zeros(M,dtype='d')
for m in range(M):
    U0up[m] = tanh(zzz_up[m]/_delta)/tanh(1.0/_delta)
Txyup = dot(D1,U0up)
Txxup = 2*Wi*Txyup*Txyup
U0up_p = dot(D1,U0up)

U0down = zeros(M,dtype='d')
for m in range(M):
    U0down[m] = tanh(zzz_down[m]/_delta)/tanh(1.0/_delta)
Txydown = dot(D1,U0down)
Txxdown = 2*Wi*Txydown*Txydown
U0down_p = dot(D1,U0down)


fileName = 'b{bbeta}-M{M}-Re{Re}-Wi{Wi}-del{delta}.dat'.format(bbeta=bbeta, M=M,
                                                         Re=Re,  Wi=Wi,
                                                           delta=_delta)

LPL = D2 - alpha*alpha*II 

## LHS

## ux[0:M] uy[M:2*M] p[2*M:3*M]
## sxx[3*M:4*M] sxy[4*M:5*M] syy[5*M:6*M]

LHS = zeros((12*M,12*M),dtype='D')

##### Upper half

# NS x

LHS[0:M,0:M]     = -Re*1j*alpha*U0up*II + bbeta*LPL
LHS[0:M,M:2*M]   = -Re*U0up_p*II
LHS[0:M,2*M:3*M] = -1j*alpha*II
LHS[0:M,3*M:4*M] = (1.0-bbeta)*1j*alpha*II
LHS[0:M,4*M:5*M] = (1.0-bbeta)*D1

# NS y

LHS[M:2*M,M:2*M]   = -Re*1j*alpha*U0up*II + bbeta*LPL
LHS[M:2*M,2*M:3*M] = -D1
LHS[M:2*M,4*M:5*M] = (1.0-bbeta)*1j*alpha*II
LHS[M:2*M,5*M:6*M] = (1.0-bbeta)*D1

# Incompressibility

LHS[2*M:3*M,0:M] = 1j*alpha*II
LHS[2*M:3*M,M:2*M] = D1

## Oldroyd-B xx

LHS[3*M:4*M,0:M]     = -2*Wi*1j*alpha*Txxup*II - 2*Wi*dot(Txyup*II,D1) - 2*1j*alpha*II
LHS[3*M:4*M,M:2*M]   = Wi*dot(D1,Txxup)*II
LHS[3*M:4*M,3*M:4*M] = II + 1j*alpha*Wi*U0up*II
LHS[3*M:4*M,4*M:5*M] = -2*Wi*U0up_p*II

## Oldroyd-B xy

LHS[4*M:5*M,0:M]     = -D1
LHS[4*M:5*M,M:2*M]   = Wi*dot(D1,Txyup)*II - Wi*1j*alpha*Txxup*II - 1j*alpha*II
LHS[4*M:5*M,4*M:5*M] = II + 1j*alpha*Wi*U0up*II
LHS[4*M:5*M,5*M:6*M] = -Wi*U0up_p*II

## Oldroyd-B yy

LHS[5*M:6*M,M:2*M]  = -2*Wi*Txyup*1j*alpha*II - 2*D1
LHS[5*M:6*M,5*M:6*M] = II + 1j*alpha*Wi*U0up*II


##### Lower half

# NS x

LHS[6*M:7*M,6*M:7*M]     = -Re*1j*alpha*U0down*II + bbeta*LPL
LHS[6*M:7*M,7*M:8*M]   = -Re*U0down_p*II
LHS[6*M:7*M,8*M:9*M] = -1j*alpha*II
LHS[6*M:7*M,9*M:10*M] = (1.0-bbeta)*1j*alpha*II
LHS[6*M:7*M,10*M:11*M] = (1.0-bbeta)*D1

# NS y

LHS[7*M:8*M,7*M:8*M]   = -Re*1j*alpha*U0down*II + bbeta*LPL
LHS[7*M:8*M,8*M:9*M] = -D1
LHS[7*M:8*M,10*M:11*M] = (1.0-bbeta)*1j*alpha*II
LHS[7*M:8*M,11*M:12*M] = (1.0-bbeta)*D1

# Incompressibility

LHS[8*M:9*M,6*M:7*M] = 1j*alpha*II
LHS[8*M:9*M,7*M:8*M] = D1

## Oldroyd-B xx

LHS[9*M:10*M,6*M:7*M]     = -2*Wi*1j*alpha*Txxdown*II - 2*Wi*dot(Txydown*II,D1) - 2*1j*alpha*II
LHS[9*M:10*M,7*M:8*M]   = Wi*dot(D1,Txxdown)*II
LHS[9*M:10*M,9*M:10*M] = II + 1j*alpha*Wi*U0down*II
LHS[9*M:10*M,10*M:11*M] = -2*Wi*U0down_p*II

## Oldroyd-B xy

LHS[10*M:11*M,6*M:7*M]     = -D1
LHS[10*M:11*M,7*M:8*M]   = Wi*dot(D1,Txydown)*II - Wi*1j*alpha*Txxdown*II - 1j*alpha*II
LHS[10*M:11*M,10*M:11*M] = II + 1j*alpha*Wi*U0down*II
LHS[10*M:11*M,11*M:12*M] = -Wi*U0down_p*II

## Oldroyd-B yy

LHS[11*M:12*M,7*M:8*M]  = -2*Wi*Txydown*1j*alpha*II - 2*D1
LHS[11*M:12*M,11*M:12*M] = II + 1j*alpha*Wi*U0down*II

## RHS

RHS = zeros((12*M,12*M),dtype='D')

RHS[0:M,0:M]         = Re*II
RHS[M:2*M,M:2*M]     = Re*II
RHS[3*M:4*M,3*M:4*M] = -Wi*II
RHS[4*M:5*M,4*M:5*M] = -Wi*II
RHS[5*M:6*M,5*M:6*M] = -Wi*II

RHS[6*M:7*M,6*M:7*M]     = Re*II
RHS[7*M:8*M,7*M:8*M]     = Re*II
RHS[9*M:10*M,9*M:10*M]   = -Wi*II
RHS[10*M:11*M,10*M:11*M] = -Wi*II
RHS[11*M:12*M,11*M:12*M] = -Wi*II

## Boundary conditions

LHS[0]     = zeros(12*M,dtype='D')
LHS[M-1]   = zeros(12*M,dtype='D')
LHS[M]     = zeros(12*M,dtype='D')
LHS[2*M-1] = zeros(12*M,dtype='D')

LHS[6*M]     = zeros(12*M,dtype='D')
LHS[7*M-1]   = zeros(12*M,dtype='D')
LHS[7*M]     = zeros(12*M,dtype='D')
LHS[8*M-1] = zeros(12*M,dtype='D')


LHS[0,0]         = 1.0
LHS[M,M]         = 1.0
LHS[7*M-1,7*M-1] = 1.0
LHS[8*M-1,8*M-1] = 1.0

LHS[M-1,M-1] = 1.0
LHS[M-1,6*M] = -1.0
LHS[2*M-1,2*M-1] = 1.0
LHS[2*M-1,7*M] = -1.0

#LHS[6*M,0:M] = D1[M-1]
#LHS[6*M,6*M:7*M] = -D1[0]

#LHS[7*M,M:2*M] = D1[M-1]
#LHS[7*M,7*M:8*M] = -D1[0]

LHS[6*M,0:M]   = bbeta*D1[M-1]
LHS[6*M,2*M-1] = 1j*alpha*bbeta
LHS[6*M,5*M-1] = (1.0-bbeta)
LHS[6*M,6*M:7*M]   = -bbeta*D1[0]
LHS[6*M,7*M] = -1j*alpha*bbeta
LHS[6*M,10*M] = -(1.0-bbeta)

LHS[7*M,M:2*M] = 2*bbeta*D1[M-1]
LHS[7*M,3*M-1] = -1.0
LHS[7*M,6*M-1] = (1.0-bbeta)
LHS[7*M,7*M:8*M] = -2*bbeta*D1[0]
LHS[7*M,8*M] = 1.0
LHS[7*M,11*M] = -(1.0-bbeta)


RHS[0]     = zeros(12*M,dtype='D')
RHS[M-1]   = zeros(12*M,dtype='D')
RHS[M]     = zeros(12*M,dtype='D')
RHS[2*M-1] = zeros(12*M,dtype='D')

RHS[6*M]   = zeros(12*M,dtype='D')
RHS[7*M-1] = zeros(12*M,dtype='D')
RHS[7*M]   = zeros(12*M,dtype='D')
RHS[8*M-1] = zeros(12*M,dtype='D')


(eigenvals,eigenvecs) = eig(LHS,RHS,left=False,right=True)

large_evs = zeros(len(eigenvals))
for i in xrange(12*M):
    if (real(eigenvals[i]) > 0) and (real(eigenvals[i]) < 50):
        large_evs[i] = real(eigenvals[i])
del i

lead_index = argmax(large_evs)

eigarray = vstack((real(eigenvals), imag(eigenvals))).T

print 'for the chosen eigenvector the leading index is: {e}'.format(e=lead_index)
print 'with value: ', eigarray[lead_index,:]
savetxt('ev-k{k}-{file}'.format(k=alpha, file=fileName), eigarray)

du =   concatenate((eigenvecs[   :M,   lead_index], eigenvecs[6*M:7*M,
                                                              lead_index]))
dv =   concatenate((eigenvecs[  M:2*M, lead_index], eigenvecs[7*M:8*M,
                                                              lead_index]))
dp =   concatenate((eigenvecs[2*M:3*M, lead_index], eigenvecs[8*M:9*M,
                                                              lead_index]))
dtxx = concatenate((eigenvecs[3*M:4*M, lead_index], eigenvecs[9*M:10*M,
                                                              lead_index]))
dtxy = concatenate((eigenvecs[4*M:5*M, lead_index], eigenvecs[10*M:11*M,
                                                              lead_index]))
dtyy = concatenate((eigenvecs[5*M:6*M, lead_index], eigenvecs[11*M:12*M,
                                                              lead_index]))

psi = -dv / (1.j*alpha)

pickle.dump((psi,du,dv,dp,dtxx,dtyy,dtxy), open('evec-k{k}-{file}'.format(k=alpha, file=fileName), 'w'))


