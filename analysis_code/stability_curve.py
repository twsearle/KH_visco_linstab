#------------------------------------------------------------------------------
#   Process data to make a stability curve.
#   Last modified: Sun 19 May 16:45:45 2013
#
#------------------------------------------------------------------------------

from scipy import *
import glob
import re

#-----------------------------------------

M = 150

Re = 10000.0
bbeta = 0.1

_delta = 0.1

#-----------------------------------------

# MAIN

fileList = glob.glob('*.dat')

outstream = open('stability-beta{bbeta}-Re{Re}.dat'.format(bbeta=bbeta, Re=Re),
               'w')

for filename in fileList:
    splitString = re.split('-',filename)
    MString = splitString[0]
    varString = splitString[1]

    m = re.match('M',MString)
    if m is  None:
        continue

    files_M = int(MString[m.end():])

    if files_M is M:
        print files_M
        print varString 
        varMatch = re.match('Wi',varString)
        if varMatch is None:
            continue
        filesVar = float(varString[varMatch.end():-4])
        print filesVar

        dispersionData = genfromtxt(filename)
        maxIndex = argmax(dispersionData[:,1])

        outstream.write('{var}, {k}, {eig}\n'.format(
                        var=filesVar,
                        k=dispersionData[maxIndex,0],
                        eig=dispersionData[maxIndex,1]))

outstream.close()



