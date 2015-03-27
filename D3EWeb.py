'''
D3E-Web
Discrete Distributional Differential Expression Web Tool

Author: Mihails Delmans (md656@cam.ac.uk)
Advisor: Martin Hemberg (mh26@sanger.ac.uk)
Version: 1.0

Tested with:
scipy 0.13.0b1
numpy 1.8.0rc1
sympy.mpmath 0.18
'''

from __future__ import division
from D3EUtil import readData, getParamsMoments, cramerVonMises, logStatus, goodnessOfFit
from D3EUtil import Params, BioParams, Status 

from numpy import log2, mean
from scipy.stats import variation

import StringIO
import json
import sys
import fileinput

def getBioParams(params):
	alpha = params.alpha
	beta = params.beta
	gamma = params.gamma

	return BioParams( gamma / beta, alpha, alpha / (alpha + beta) )

inputString = fileinput.input()[0].rstrip('\n')

inputJSON = json.loads(inputString)
inputFile = StringIO.StringIO(inputJSON['input']);
label1 = inputJSON['label1']
label2 = inputJSON['label2']
normalise = int(inputJSON['normalise'])
removeZeros = inputJSON['removeZeros']
useSpikeIns = inputJSON['useSpikeIns']

data1, data2, ids, lineStatus = readData(inputFile, label1, label2, normalise, removeZeros, useSpikeIns)

errors = []
rows = []
cols = []

for status in lineStatus:
	if status.code == 1:
		errors.append( {'geneId': status.idx, 'message': 'Warning: ' + status.message}  )
	elif status.code == 2:
		errors.append( {'geneId': status.idx, 'message': 'Error: ' + status.message}  )


rowKeys = ["geneId"]
rowKeys.extend( ["col" + str(i) for i in range(1,21)] )

labels = ["Gene ID", "a1", "b1", "g1", "a2", "b2", "g2",
					 "s1", "f1", "d1", "s2", "f2", "d2",
					 "Rs", "Rf", "Rd", "mu1", "mu2", "cv1", "cv2",  "p"]

for key, label in zip(rowKeys, labels):
	cols.append( {"key" : key, "label" : label} )



for p1,p2,idx in zip(data1, data2, ids):

	difference = cramerVonMises(p1,p2)
	
	if difference == -1:
		errors.append({"geneId" : idx, "type" : "Warning",
		"message" : "Could not estimate Cramer-von Mises statistic. Further analysis aborted." })
		continue

	params1 = getParamsMoments(p1)
	params2 = getParamsMoments(p2)

	if (params1.alpha <= 0 or params2.alpha <= 0 or params1.beta <= 0 or params2.beta <= 0 or
		params1.gamma <= 0 or params2.gamma <= 0):
			errors.append({"geneId" : idx, "type" : "Warning",
			"message" : "Could not estimate parameters for atleast one cell type." })
	else:
		bioParams1 = getBioParams(params1)
		bioParams2 = getBioParams(params2)

		rowString = ('{0:s}\t{1:2.2f}\t{2:2.2f}\t{3:2.2f}\t{4:2.2f}\t{5:2.2f}\t{6:2.2f}\t'
		    '{7:2.2f}\t{8:2.2f}\t{9:2.2f}\t{10:2.2f}\t{11:2.2f}\t{12:2.2f}\t'
								'{13:2.2f}\t{14:2.2f}\t{15:2.2f}\t{16:2.2f}\t{17:2.2f}\t{18:2.2f}\t{19:2.2f}\t'
								'{20:2.2e}\n').format(idx,
								params1.alpha, params1.beta, params1.gamma,
								params2.alpha, params2.beta, params2.gamma,
								bioParams1.size, bioParams1.freq, bioParams1.duty,
								bioParams2.size, bioParams2.freq, bioParams2.duty,
								log2( bioParams1.size / bioParams2.size),
								log2( bioParams1.freq / bioParams2.freq),
								log2( bioParams1.duty / bioParams2.duty),
								mean(p1), mean(p2), variation(p1), variation(p2),
								difference)
		rowVals = rowString.split('\t')

		rows.append( dict( zip( rowKeys, rowVals ) ) )


jsonOutput = { "rows" : rows, "columns" : cols, "errors" : errors }
print json.dumps(jsonOutput)

