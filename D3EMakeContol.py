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
from D3EUtil import readData

inFileName = sys.argv[1]
outFileName = sys.argv[2]
label = sys.argv[3]

inFile = open(inFileName)
outFileName = open(outFileName,'w')

p1, p2, ids, lineStatus = readData(inFile, label, label, normalise=False, removeZeros=False, useSpikeIns = False, verbose = True)

n = len(p1[0]) / 2

outFileName.write('GeneID\t' + '\t'.join( [ label+'_1' ] * n + [label+'_2'] * (len(p1[0])-n) ) + '\n')

for idx, p in zip(ids, p1):
	outFileName.write(idx + '\t' + '\t'.join( str(int(x)) for x in p ) + '\n')
