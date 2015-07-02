'''

D3E-Merge-Output
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
from os import listdir

inDir = sys.argv[1]
outFileName = sys.argv[2]
label1 = sys.argv[3]
label2 = sys.argv[4]

outputFile = open(outFileName, 'w')
hr = False

for fileName in listdir(inDir):
	if fileName.startswith(label1+'_'+label2):
		inFile = open(inDir + fileName,'r')
		h = inFile.readline()
		
		if not hr:
			outputFile.write(h)
			hr = True

		inFile.readline()

		for line in inFile:
			if line.strip():
				outputFile.write(line)
