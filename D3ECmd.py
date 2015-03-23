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
from ThreeDExpressUtil import readData, getParamsBayesian, cramerVonMises, logStatus, goodnessOfFit
from ThreeDExpressUtil import Params, BioParams, Status 

from argparse import ArgumentParser, FileType
from numpy import mean, log2

def checkCramerVonMises(pValue, comment, idx):
	if pValue == -1:
		logStatus( Status(1, idx, "Could not estimate Cramer-von Mises statistic(" + comment + ")") )

parser = ArgumentParser(description='3DExpress')
parser.add_argument(action = 'store', type=FileType('r'), dest = 'inputFile', help='Input file...')
parser.add_argument(action = 'store', type=FileType('w'), dest = 'outputFile', help='Output file...')
parser.add_argument(action = 'store', type = str, nargs = 1, dest = 'label1', help = 'Label for the first cell type / condition')
parser.add_argument(action = 'store', type = str, nargs = 1, dest = 'label2', help = 'Label for the second cell type / condition')

parser.add_argument('-m', '--mode', action = 'store', type = int, dest='mode', default = 1, choices =[0,1,2,3], help = 'Mode of operation\n ' \
													'0: Output pValue only\n'
													'1: Output parameters and pValue (default)\n'\
													'2: Additionally outputs GOF\n'
													'3: Additionally outputs GOF and means\n')
parser.add_argument('-n','--normalise', action = 'store', type = int, dest='normalise', default = 1, choices = [0,1], help='Normalise the data')
parser.add_argument('-z', '--removeZeros', action = 'store', type = int, dest='removeZeros', default = 0, choices = [0,1], help='Remove zeros')
parser.add_argument('-s', '--useSpikeIns', action = 'store', type = int, dest='useSpikeIns', default = 0, choices = [0,1], help='Use spike-ins for normalisation. Requires at least one row with id starting with "spike" ')
parser.add_argument('-v', action = 'store_const', const = True, dest = 'verbose', default = False, help = 'verbose')
args = parser.parse_args()

data1, data2, ids, lineStatus = readData(args.inputFile, args.label1[0], args.label2[0], args.normalise, args.removeZeros, args.useSpikeIns)

if args.verbose:
	for status in lineStatus:
		logStatus(status)

for p1,p2,idx in zip(data1, data2, ids):

	difference = cramerVonMises(p1,p2)
	
	if difference == -1:
		logStatus( Status(1, idx, "Could not estimate Cramer-von Mises statistic. Further analysis aborted.") )
		continue

	if args.mode > 0:
		params1, bioParams1 = getParamsBayesian(p1)
		params2, bioParams2 = getParamsBayesian(p2)

		tauP = cramerVonMises(bioParams1.tau.sample, bioParams2.tau.sample)
		checkCramerVonMises(tauP, 'for tau', idx)

		fP = cramerVonMises(bioParams1.f.sample, bioParams2.f.sample)
		checkCramerVonMises(fP, 'for f', idx)

		tP = cramerVonMises(bioParams1.t.sample, bioParams2.t.sample)
		checkCramerVonMises(tP, 'for t', idx)

	if args.mode > 1:
		gof1 = goodnessOfFit(p1, params1)
		checkCramerVonMises(gof1, 'Goodnes of fit for a first cell type', idx)

		gof2 = goodnessOfFit(p2, params2)
		checkCramerVonMises(gof1, 'Goodnes of fit for a second cell type', idx)

	if args.mode == 0:
		args.outputFile.write('{0:s}\t{1:2.2e}\t{2:2.2f}\t{3:2.2f}\n'.format(idx, difference, mean(p1), mean(p2)))
	elif args.mode == 1:
		args.outputFile.write('{0:s}\t\t{1:2.2f}\t{2:2.2f}\t{3:2.2f}\t\t{4:2.2f}\t{5:2.2f}\t{6:2.2f}\t\t'
								'{7:2.2f}\t{8:2.2f}\t{9:2.2f}\t\t{10:2.2f}\t{11:2.2f}\t{12:2.2f}\t\t'
								'{13:2.2f}\t{14:2.2f}\t{15:2.2f}\t'
								'{16:2.2e}\t{17:2.2e}\t{18:2.2e}\t\t{19:2.2e}\n'.format(idx,
								params1.alpha.mean(), params1.beta.mean(), params1.gamma.mean(),
								params2.alpha.mean(), params2.beta.mean(), params2.gamma.mean(),
								bioParams1.tau.mean(), bioParams1.f.mean(), bioParams1.t.mean(),
								bioParams2.tau.mean(), bioParams2.f.mean(), bioParams2.t.mean(),
								log2( bioParams2.tau.mean() / bioParams1.tau.mean()),
								log2( bioParams2.f.mean() / bioParams1.f.mean()),
								log2( bioParams2.t.mean() / bioParams1.t.mean()),
								tauP,fP, tP, difference) )
	elif args.mode == 2:
		args.outputFile.write('{0:s}\t\t{1:2.2f}\t{2:2.2f}\t{3:2.2f}\t{4:2.2e}\t\t'
								'{5:2.2f}\t{6:2.2f}\t{7:2.2f}\t{8:2.2e}\t\t'
								'{9:2.2f}\t{10:2.2f}\t{11:2.2f}\t\t{12:2.2f}\t{13:2.2f}\t{14:2.2f}\t\t'
								'{15:2.2f}\t{16:2.2f}\t{17:2.2f}\t'
								'{18:2.2e}\t{19:2.2e}\t{20:2.2e}\t\t{21:2.2e}\n'.format(idx,
								params1.alpha.mean(), params1.beta.mean(), params1.gamma.mean(), gof1,
								params2.alpha.mean(), params2.beta.mean(), params2.gamma.mean(), gof2,
								bioParams1.tau.mean(), bioParams1.f.mean(), bioParams1.t.mean(),
								bioParams2.tau.mean(), bioParams2.f.mean(), bioParams2.t.mean(),
								log2( bioParams2.tau.mean() / bioParams1.tau.mean()),
								log2( bioParams2.f.mean() / bioParams1.f.mean()),
								log2( bioParams2.t.mean() / bioParams1.t.mean()),
								tauP,fP, tP, difference) )
	elif args.mode == 3:
		args.outputFile.write('{0:s}\t\t{1:2.2f}\t{2:2.2f}\t{3:2.2f}\t{4:2.2e}\t\t'
								'{5:2.2f}\t{6:2.2f}\t{7:2.2f}\t{8:2.2e}\t\t'
								'{9:2.2f}\t{10:2.2f}\t{11:2.2f}\t\t{12:2.2f}\t{13:2.2f}\t{14:2.2f}\t\t'
								'{15:2.2f}\t{16:2.2f}\t{17:2.2f}\t'
								'{18:2.2e}\t{19:2.2e}\t{20:2.2e}\t\t{21:2.2e}\t\t{22:2.2f}\t{23:2.2f}\n'.format(idx,
								params1.alpha.mean(), params1.beta.mean(), params1.gamma.mean(), gof1,
								params2.alpha.mean(), params2.beta.mean(), params2.gamma.mean(), gof2,
								bioParams1.tau.mean(), bioParams1.f.mean(), bioParams1.t.mean(),
								bioParams2.tau.mean(), bioParams2.f.mean(), bioParams2.t.mean(),
								log2( bioParams2.tau.mean() / bioParams1.tau.mean()),
								log2( bioParams2.f.mean() / bioParams1.f.mean()),
								log2( bioParams2.t.mean() / bioParams1.t.mean()),
								tauP,fP, tP, difference, mean(p1), mean(p2)) )

	logStatus( Status(0, idx, "Analysis complete.") )

args.outputFile.close()
