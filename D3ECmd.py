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

'''
from __future__ import division

from D3EUtil import readData, getParamsBayesian, getParamsMoments, cramerVonMises, logStatus, goodnessOfFit
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

parser.add_argument('-m', '--mode', action = 'store', type = int, dest='mode', default = 2, choices =[1,2], help = 'Mode of operation\n ' \
													'1: Apply Method of moments\n'
													'2: Apply Bayesian approach\n')
parser.add_argument('-n','--normalise', action = 'store', type = int, dest='normalise', default = 1, choices = [0,1], help='Normalise the data')
parser.add_argument('-z', '--removeZeros', action = 'store', type = int, dest='removeZeros', default = 0, choices = [0,1], help='Remove zeros')
parser.add_argument('-s', '--useSpikeIns', action = 'store', type = int, dest='useSpikeIns', default = 0, choices = [0,1], help='Use spike-ins for normalisation. Requires at least one row with id starting with "spike" ')
parser.add_argument('-v', action = 'store_const', const = True, dest = 'verbose', default = False, help = 'verbose')
args = parser.parse_args()

data1, data2, ids, lineStatus = readData(args.inputFile, args.label1[0], args.label2[0], args.normalise, args.removeZeros, args.useSpikeIns)

if args.verbose:
	for status in lineStatus:
		logStatus(status)


if args.mode == 1:
	args.outputFile.write('#GeneID\ta1\tb1\tg1\tGOF1\t\ta2\tb2\tg2\tGOF2\t\ts1\tf1\td1\t\ts2\tf2\td2\t\tRs\tRf\tRd\t\tp-value\t\tmu1\tcv1\t\tmu2\tcv2\n\n')
elif args.mode == 2:
	args.outputFile.write('#GeneID\ta1\tb1\tg1\tGOF1\t\ta2\tb2\tg2\tGOF2\t\ts1\tf1\td1\t\ts2\tf2\td2\t\tRs\tRf\tRd\t\tpSize\tpFreq\tpDuty\t\tp-value\t\tmu1\tcv1\t\tmu2\tcv2\n\n')

for p1,p2,idx in zip(data1, data2, ids):

	difference = cramerVonMises(p1,p2)
	
	if difference == -1:
		logStatus( Status(1, idx, "Could not estimate Cramer-von Mises statistic. Further analysis aborted.") )
		continue

	if args.mode == 1:
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
	
	elif args.mode == 2:

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

	if args.mode == 1:
		args.outputFile.write('{0:s}\t\t{1:4.4f}\t{2:4.4f}\t{3:4.4f}\t{4:4.4e}\t\t'
								'{5:4.4f}\t{6:4.4f}\t{7:4.4f}\t{8:4.4e}\t\t'
								'{9:4.4f}\t{10:4.4f}\t{11:4.4f}\t\t{12:4.4f}\t{13:4.4f}\t{14:4.4f}\t\t'
								'{15:4.4f}\t{16:4.4f}\t{17:4.4f}\t\t'
								'{18:4.4e}\t\t{19:4.4f}\t{20:4.4f}\t\t{21:4.4f}\t{22:4.4f}\n'.format(idx,
								params1.alpha, params1.beta, params1.gamma, gof1,
								params2.alpha, params2.beta, params2.gamma, gof2,
								bioParams1.size, bioParams1.freq, bioParams1.duty,
								bioParams2.size, bioParams2.freq, bioParams2.duty,
								log2( bioParams1.size / bioParams2.size ),
								log2( bioParams1.freq / bioParams2.freq ),
								log2( bioParams1.duty / bioParams2.duty ),
								difference, mean(p1), variation(p1), mean(p2), variation(p2)) )
	elif args.mode == 2:
		args.outputFile.write('{0:s}\t\t{1:4.4f}\t{2:4.4f}\t{3:4.4f}\t{4:4.4e}\t\t'
								'{5:4.4f}\t{6:4.4f}\t{7:4.4f}\t{8:4.4e}\t\t'
								'{9:4.4f}\t{10:4.4f}\t{11:4.4f}\t\t{12:4.4f}\t{13:4.4f}\t{14:4.4f}\t\t'
								'{15:4.4f}\t{16:4.4f}\t{17:4.4f}\t'
								'{18:4.4e}\t{19:4.4e}\t{20:4.4e}\t\t{21:4.4e}\t\t{22:4.4f}\t{23:4.4f}\t\t{24:4.4f}\t{25:4.4f}\n'.format(idx,
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
