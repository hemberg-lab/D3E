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
