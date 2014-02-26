#------------------------------------------------------------------------------
#   Make a 2D plot using matplotlib
#   Last modified: Wed 23 Oct 15:13:03 2013
#
#------------------------------------------------------------------------------

from scipy import *
import matplotlib.pyplot as plt
from matplotlib import rc
# This next import is supposed to make sure fonts are good, regardless of os
from matplotlib.font_manager import fontManager, FontProperties

# MAIN

fig_width_pt = 452.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27                # Convert pt to inch
golden_ratio = (sqrt(5)+1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt   # width in inches
fig_height = fig_width/golden_ratio      # height in inches
fig_size =  [fig_width,fig_height]

#Use the Latex distribution for text.
rc('text', usetex=True)
#Make sure font is the same as that in the template/style file!
rc('font', family='serif', size='12')
rc('figure', figsize=fig_size)

lowBeta = genfromtxt('stability-beta0.1-Re10000.0.dat')
x1,y1 = lowBeta[:,0]/0.1, lowBeta[:,2]*0.1
midBeta = genfromtxt('stability-beta0.5-Re10000.0.dat')
x2,y2 = midBeta[:,0]/0.1, midBeta[:,2]*0.1
highBeta = genfromtxt('stability-beta0.9-Re10000.0.dat')
x3,y3 = highBeta[:,0]/0.1, highBeta[:,2]*0.1
outFilename = 'KH_high_Re_vary_Wi.pdf'

fig = plt.figure()

ax = plt.subplot(111)

ax.plot(x1,y1,'-',linewidth=1)
ax.plot(x2,y2,'-',linewidth=1)
ax.plot(x3,y3,'-',linewidth=1)

# Insert a horizontal line to mark the solution for the Newtonian instability
plt.axhline(y=0.1880, ls='--', linewidth=0.6,  color='black')

# Insert point at yaxis intercept
theIntercept = lowBeta[0,2]*0.1
ax.plot([0], [theIntercept], 'o', color='red')

#Change the yaxis ticks to include the value for the intercept
#ticks = ax.yaxis.get_majorticklocs()
#print ticks
#ticks = concatenate((ticks[:-1],[theIntercept],ticks[-1:]))
#print ticks
#ax.yaxis.set_ticks(ticks)

plt.legend( (r'$\beta = 0.1$',r'$\beta = 0.5$',r'$\beta = 0.9$'), loc='best' )
ax.set_xlim(0.0,200)
ax.set_ylim(0.18,0.195)
ax.set_xlabel(r'Wi')
ax.set_ylabel(r'$\lambda_{r}$')

plt.savefig(outFilename)

plt.show()
