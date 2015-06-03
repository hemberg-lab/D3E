'''

D3E-Cmd
Discrete Distributional Differential Expression Command Line Tool

Author: Mihails Delmans (md656@cam.ac.uk)
Advisor: Martin Hemberg (mh26@sanger.ac.uk)
Version: 1.0

Tested with:
scipy 0.15.1
numpy 1.8.0rc1
sympy.mpmath 0.18

'''

import sys

tFileName = sys.argv[1]
cFileName = sys.argv[2]
oFileName = sys.argv[3]
k = float(sys.argv[4])


cFile = open(cFileName)
cFile.readline()
cFile.readline()

pMin = float('Inf')

for line in cFile:
	p = float(line.split()[21])
	if p < pMin:
		pMin = p

tFile = open(tFileName)
oFile = open(oFileName,'w')

h = tFile.readline()
tFile.readline()

oFile.write(h + '\n')

for line in tFile:
	tabs = line.split()
	p = float(tabs[21])
	if p < pMin * k:
		oFile.write(line)
