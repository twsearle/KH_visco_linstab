#------------------------------------------------------------------------------
#   Make a 2D plot using matplotlib
#   Last modified: Wed 23 Oct 15:05:51 2013
#
#------------------------------------------------------------------------------

from scipy import *
import matplotlib.pyplot as plt
from matplotlib import rc

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


lowBeta = genfromtxt('beta_vs_Wi_stability.dat')
outFilename = 'onset_beta_Wi.pdf'

plt.plot(lowBeta[:,0], lowBeta[:,1]/0.1, 'ro') 

#plt.axis([0, 100, 0, 1])

plt.ylabel(r'$Wi$ at the onset of instability')
plt.xlabel(r'$\beta$')

fig = plt.gcf()

plt.savefig(outFilename)

plt.show()
