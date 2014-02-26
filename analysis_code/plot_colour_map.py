#
#   Plot eigenvector
#
"""Convert grid of GL points/Fourier points back to real space, then plot in a
   colour map"""

from matplotlib import pyplot as plt
import cPickle as pickle
from scipy import *
from scipy import interpolate
# PARAMETERS

M = 150

Re = 0.0
Wi = 10.
bbeta = 0.1

_delta = 0.1
alpha = 0.28

XPOINTS = 100

#FUNCTIONS
def my_interp(vec):

    f = interpolate.interp1d(y_points[::-1],vec[::-1], bounds_error=False,
                         kind='linear')
    return f(y_new)

def fourier_transform(vec):
    mat = zeros((2*M, XPOINTS), dtype='complex')
    for m in xrange (2*M):
        for n in xrange(XPOINTS):
            mat[m,n] = vec[m] * exp(1.j*alpha*x_points[n])

    del m,n
    return mat

# MAIN

ygl = zeros(M,dtype='d')
for m in range(M):
    ygl[m] = cos(pi*m/(M-1))

zzz_up = zeros(M,dtype='d')
for m in range(M):
    zzz_up[m] = 0.5*(ygl[m]+1.0)

zzz_down = zeros(M,dtype='d')
for m in range(M):
    zzz_down[m] = 0.5*(ygl[m]-1.0)

y_points = concatenate((zzz_up, zzz_down))

x_points = zeros(XPOINTS,dtype='d') 
for n in range(XPOINTS):
    #               2lambda     * fractional position
    x_points[n] = (2.*pi/alpha) * ((1.*n)/XPOINTS)
del n

y_new = zeros(2*M)
for m in range(2*M):
    y_new[m] = -1 + ((1.*m)/M)
del m

fileName = \
'evec-k{k}-b{b}-M{M}-Re{Re}-Wi{Wi}-del{delta}.dat'.format(k=alpha,
                                                              b=bbeta, M=M,
                                                              Re=Re,  Wi=Wi,
                                                              delta=_delta)
titleString = \
' for k = {k}, b = {b}, Re = {Re}, Wi = {Wi}, delta = {delt}'.format(k=alpha,
                                                                b=bbeta,
                                                                Re=Re,  Wi=Wi,
                                                                delt=_delta)

(psi,du,dv,dp,dtxx,dtyy,dtxy) = pickle.load(open(fileName,'r'))

#f = interpolate.interp1d(y_points[::-1],psi[::-1], bounds_error=False,
#                         kind='linear')

#plt.plot(y_points,real(psi), 'o', y_points, imag(psi))
#plt.show()

psi = my_interp(psi)
du = my_interp(du)
dv = my_interp(dv)
dp = my_interp(dp)
dtxx = my_interp(dtxx)
dtyy = my_interp(dtyy)
dtxy = my_interp(dtxy)

dtxx2D = fourier_transform(dtxx)
normalisation = dtxx2D[M,XPOINTS/2]
dtxx2D = dtxx2D/normalisation
psi2D = fourier_transform(psi)/normalisation
du2D   = fourier_transform(du)/normalisation
dv2D = fourier_transform(dv)/normalisation
dp2D = fourier_transform(dp)/normalisation
dtyy2D = fourier_transform(dtyy)/normalisation
dtxy2D = fourier_transform(dtxy)/normalisation

grid_x, grid_y = meshgrid(x_points,y_new)

#Proof that the old coordinate system is upside down
#savetxt('rgh.txt', vstack((y_points,y_new)).T)

plt.figure()
plt.imshow(real(psi2D), origin='lower', extent=[0,((2.*pi)/alpha),-1,1], aspect=4)
plt.contour(grid_x, grid_y, real(psi2D))
plt.colorbar(orientation='horizontal')
plt.axhline(y=_delta, ls='--', linewidth=0.5,  color='black')
plt.axhline(y=-_delta, ls='--', linewidth=0.5,  color='black')
plt.title('psi'+titleString)
plt.savefig('psi.pdf')

plt.figure()
plt.imshow(real(dtxx2D), origin='lower', extent=[0,((2.*pi)/alpha),-1,1], aspect=4)
plt.contour(grid_x, grid_y, real(dtxx2D))
plt.colorbar(orientation='horizontal')
plt.axhline(y=_delta, ls='--', linewidth=0.5,  color='black')
plt.axhline(y=-_delta, ls='--', linewidth=0.5,  color='black')
plt.title('txx'+titleString)
plt.savefig('txx.pdf')

plt.figure()
plt.imshow(real(dtxy2D), origin='lower', extent=[0,((2.*pi)/alpha),-1,1], aspect=4)
plt.contour(grid_x, grid_y, real(dtxy2D))
plt.colorbar(orientation='horizontal')
plt.axhline(y=_delta, ls='--', linewidth=0.5,  color='black')
plt.axhline(y=-_delta, ls='--', linewidth=0.5,  color='black')
plt.title('txy'+titleString)
plt.savefig('txy.pdf')
