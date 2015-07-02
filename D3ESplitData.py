'''
D3E-Split-Data
Discrete Distributional Differential Expression Command Line Tool

Author: Mihails Delmans (md656@cam.ac.uk)
Advisor: Martin Hemberg (mh26@sanger.ac.uk)
Version: 1.0

Tested with:
scipy 0.13.0b1
numpy 1.8.0rc1
sympy.mpmath 0.18

Copyright 2015 Mihails Delmans, Martin Hemberg

This file is part of D3E.

D3E is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

D3E is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with D3E.  If not, see <http://www.gnu.org/licenses/>.

'''

from D3EUtil import readData
import sys

inFileName = sys.argv[1]
label1 = sys.argv[2]
label2 =sys.argv[3]
outDir = sys.argv[4]
n = int(sys.argv[5])

inFile = open(inFileName)

data1, data2, ids, lineStatus = readData(inFile, label1, label2, normalise = True)

i = 0

for p1, p2, idx in zip(data1, data2, ids):
	if i % n == 0:
		ouputFile = open( outDir + label1 + '_' + label2 + '_' + str(i/n) + '.txt', 'w' )
		ouputFile.write('GeneID\t' + '\t'.join( [label1] * len(data1[0]) + [label2] * len(data2[0]) ) + '\n' )
	ouputFile.write('\t'.join ( [idx] + [str(x) for x in p1+p2] ) + '\n')
	i = i + 1