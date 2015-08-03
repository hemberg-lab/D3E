'''
D3E-Cmd
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
from __future__ import division

from D3EUtil import readData, getParamsBayesian, getParamsMoments, cramerVonMises, logStatus, goodnessOfFit, distributionTest
from D3EUtil import Params, BioParams, Status 

from argparse import ArgumentParser, FileType
from numpy import mean, log2
from scipy.stats import variation

def checkCramerVonMises(pValue, comment, idx):
	if pValue == -1:
		logStatus( Status(1, idx, "Could not estimate Cramer-von Mises statistic(" + comment + ")") )

parser = ArgumentParser(description='3DExpress')
parser.add_argument(action = 'store', type=FileType('r'), dest = 'inputFile', help='Input file...')
parser.add_argument(action = 'store', type=FileType('w'), dest = 'outputFile', help='Output file...')
parser.add_argument(action = 'store', type = str, nargs = 1, dest = 'label1', help = 'Label for the first cell type / condition')
parser.add_argument(action = 'store', type = str, nargs = 1, dest = 'label2', help = 'Label for the second cell type / condition')

parser.add_argument('-m', '--mode', action = 'store', type = int, dest='mode', default = 1, choices =[0,1], help = 'Mode of operation\n '\
													'0: Apply Method of moments\n'
													'1: Apply Bayesian approach (default)\n')

parser.add_argument('-t', '--test', action = 'store', type = int, dest = 'test', default = 0, choices = [0,1,2], help = 'Statistical test\n '\
													'0: Cramer-von Mises test (default)\n'
													'1: Kolmogorov-Smirnov test\n'
													'2: Anderson-Darling test\n')
parser.add_argument('-n','--normalise', action = 'store', type = int, dest='normalise', default = 1, choices = [0,1], help='Normalise the data')
parser.add_argument('-z', '--removeZeros', action = 'store', type = int, dest='removeZeros', default = 0, choices = [0,1], help='Remove zeros')
parser.add_argument('-s', '--useSpikeIns', action = 'store', type = int, dest='useSpikeIns', default = 0, choices = [0,1], help='Use spike-ins for normalisation. Requires at least one row with id starting with "spike" ')
parser.add_argument('--spikeInStart', action = 'store', type = str, nargs = 1, dest = 'spikeInStart', default = 'spike', help = 'A string to identify spike-ins.')
parser.add_argument('-v', action = 'store_const', const = True, dest = 'verbose', default = True, help = 'verbose')
args = parser.parse_args()

data1, data2, ids, lineStatus = readData(args.inputFile, args.label1[0], args.label2[0], args.normalise, args.removeZeros, args.useSpikeIns, args.spikeInStart)

if args.verbose:
	for status in lineStatus:
		logStatus(status)


if args.mode == 0:
	args.outputFile.write('#GeneID\ta1\tb1\tg1\tGOF1\ta2\tb2\tg2\tGOF2\ts1\tf1\td1\ts2\tf2\td2\tRs\tRf\tRd\tp-value\tmu1\tcv1\tmu2\tcv2\n\n')
elif args.mode == 1:
	args.outputFile.write('#GeneID\ta1\tb1\tg1\tGOF1\ta2\tb2\tg2\tGOF2\ts1\tf1\td1\ts2\tf2\td2\tRs\tRf\tRd\tpSize\tpFreq\tpDuty\tp-value\tmu1\tcv1\tmu2\tcv2\n\n')

for p1,p2,idx in zip(data1, data2, ids):

	difference = distributionTest(p1,p2, args.test)
	
	if difference == -1:
		logStatus( Status(1, idx, "Could not estimate p-value. Further analysis aborted.") )
		continue

	if args.mode == 0:
		params1 = getParamsMoments(p1)
		params2 = getParamsMoments(p2)

		if (params1.alpha <=0  or params1.beta <= 0 or params1.gamma <=0 or
			params2.alpha <=0  or params2.beta <= 0 or params2.gamma <=0):
				logStatus( Status(1, idx, "Could not estimate parameters.") )
				
				params1 = Params(alpha = float('nan'), beta = float('nan'), gamma = float('nan'), c = float('nan') )
				params2 = Params(alpha = float('nan'), beta = float('nan'), gamma = float('nan'), c = float('nan') )
				bioParams1 = BioParams( size = float('nan'), freq = float('nan'), duty = float('nan') )
				bioParams2 = BioParams( size = float('nan'), freq = float('nan'), duty = float('nan') )
		else:
			bioParams1 = BioParams(size = params1.gamma / params1.beta, freq = 1 / params1.alpha, duty = params1.alpha / (params1.alpha + params1.beta) )
			bioParams2 = BioParams(size = params2.gamma / params2.beta, freq = 1 / params2.alpha, duty = params2.alpha / (params2.alpha + params2.beta) )
	
	elif args.mode == 1:

		params1, bioParams1 = getParamsBayesian(p1)
		params2, bioParams2 = getParamsBayesian(p2)

		sizeP = cramerVonMises(bioParams1.size.sample, bioParams2.size.sample)
		checkCramerVonMises(sizeP, 'for size', idx)

		freqP = cramerVonMises(bioParams1.freq.sample, bioParams2.freq.sample)
		checkCramerVonMises(freqP, 'for f', idx)

		dutyP = cramerVonMises(bioParams1.duty.sample, bioParams2.duty.sample)
		checkCramerVonMises(dutyP, 'for t', idx)

	gof1 = goodnessOfFit(p1, params1)
	checkCramerVonMises(gof1, 'Goodnes of fit for a first cell type', idx)

	gof2 = goodnessOfFit(p2, params2)
	checkCramerVonMises(gof1, 'Goodnes of fit for a second cell type', idx)

	if args.mode == 0:
		args.outputFile.write('{0:s}\t{1:4.4f}\t{2:4.4f}\t{3:4.4f}\t{4:4.4e}\t'
								'{5:4.4f}\t{6:4.4f}\t{7:4.4f}\t{8:4.4e}\t'
								'{9:4.4f}\t{10:4.4f}\t{11:4.4f}\t{12:4.4f}\t{13:4.4f}\t{14:4.4f}\t'
								'{15:4.4f}\t{16:4.4f}\t{17:4.4f}\t'
								'{18:4.4e}\t{19:4.4f}\t{20:4.4f}\t{21:4.4f}\t{22:4.4f}\n'.format(idx,
								params1.alpha, params1.beta, params1.gamma, gof1,
								params2.alpha, params2.beta, params2.gamma, gof2,
								bioParams1.size, bioParams1.freq, bioParams1.duty,
								bioParams2.size, bioParams2.freq, bioParams2.duty,
								log2( bioParams1.size / bioParams2.size ),
								log2( bioParams1.freq / bioParams2.freq ),
								log2( bioParams1.duty / bioParams2.duty ),
								difference, mean(p1), variation(p1), mean(p2), variation(p2)) )
	elif args.mode == 1:
		args.outputFile.write('{0:s}\t{1:4.4f}\t{2:4.4f}\t{3:4.4f}\t{4:4.4e}\t'
								'{5:4.4f}\t{6:4.4f}\t{7:4.4f}\t{8:4.4e}\t'
								'{9:4.4f}\t{10:4.4f}\t{11:4.4f}\t{12:4.4f}\t{13:4.4f}\t{14:4.4f}\t'
								'{15:4.4f}\t{16:4.4f}\t{17:4.4f}\t'
								'{18:4.4e}\t{19:4.4e}\t{20:4.4e}\t{21:4.4e}\t{22:4.4f}\t{23:4.4f}\t{24:4.4f}\t{25:4.4f}\n'.format(idx,
								params1.alpha.mean(), params1.beta.mean(), params1.gamma.mean(), gof1,
								params2.alpha.mean(), params2.beta.mean(), params2.gamma.mean(), gof2,
								bioParams1.size.mean(), bioParams1.freq.mean(), bioParams1.duty.mean(),
								bioParams2.size.mean(), bioParams2.freq.mean(), bioParams2.duty.mean(),
								log2( bioParams2.size.mean() / bioParams1.size.mean()),
								log2( bioParams2.freq.mean() / bioParams1.freq.mean()),
								log2( bioParams2.duty.mean() / bioParams1.duty.mean()),
								sizeP,freqP, dutyP, difference, mean(p1), variation(p1), mean(p2), variation(p2) ) )

	logStatus( Status(0, idx, "Analysis complete.") )

args.outputFile.close()
