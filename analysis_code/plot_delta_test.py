#------------------------------------------------------------------------------
#   Make a 2D plot using matplotlib
#   Last modified: Wed 23 Oct 15:15:13 2013
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

FileList =[ r"../high_R_delta_test/b0.1-M150-Re20000.0-Wi0.0005-del0.05.dat",
            r"../high_R_delta_test/b0.1-M150-Re10000.0-Wi0.001-del0.1.dat",
            r"../high_R_delta_test/b0.1-M150-Re5000.0-Wi0.002-del0.2.dat" ]

low = genfromtxt(FileList[0])
x1,y1 = low[:,0]*0.1, low[:,1]*0.1
mid = genfromtxt(FileList[1])
x2,y2 = mid[:,0]*0.1, mid[:,1]*0.1
high = genfromtxt(FileList[2])
x3,y3 = high[:,0]*0.1, high[:,1]*0.1
outFilename = r'high_Re_vary_delta.pdf'

fig = plt.figure()

ax = plt.subplot(111)

ax.plot(x1*0.05,y1*0.05,'-',linewidth=1)
ax.plot(x2*0.1,y2*0.1,'-',linewidth=1)
ax.plot(x3*0.2,y3*0.2,'-',linewidth=1)

plt.legend( (r'$\Delta = 0.05$',r'$\Delta = 0.1$',r'$\Delta = 0.2$'), loc='best' )
#ax.set_xlim(0.0,100)
#ax.set_ylim(1,2.)
ax.set_xlabel(r'$k$')
ax.set_ylabel(r'$\lambda_{r}$')

plt.savefig(outFilename)

plt.show()
