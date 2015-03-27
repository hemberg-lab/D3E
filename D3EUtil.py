'''
D3E
Discrete Distributional Differential Expression

A set of functions, which are used by D3E Tool

Author: Mihails Delmans (md656@cam.ac.uk)
Advisor: Martin Hemberg (mh26@sanger.ac.uk)
Version: 1.0

Tested with:
scipy 0.13.0b1
numpy 1.8.0rc1
sympy.mpmath 0.18

'''

from __future__ import division
from scipy.special import kv, gammaln
from scipy.stats import gmean
from decimal import Decimal, getcontext
from collections import namedtuple
from rvar import RVar
from numpy import log, array, zeros, median, rint, power, hstack, hsplit, seterr, mean
from numpy.random import beta, poisson
import sympy.mpmath as mp

seterr(all='ignore')

mp.mp.dps = 30 
mp.mp.pretty = True

Params = namedtuple("Params", ["alpha", "beta", "gamma", "c"])
BioParams = namedtuple("BioParams", ["size", "freq", "duty"])
Status = namedtuple("LineStatus", ["code", "idx", "message"])

def logStatus(status):
	statusType = ['Log','Warning','Error']
	print status.idx + ' - ' + statusType[status.code] + ': ' + status.message

# Read a header of an input file, get indeces of colums that match specified labels
def _readHeader(header, label1, label2):

	if header.lower.startswith('geneid'):
		tabs = header.split('\t')

		colIdx1 = [i for i,x in enumerate(tabs) if x == label1 ]
		colIdx2 = [i for i,x in enumerate(tabs) if x == label2 ]
	else:
		return [],[], 1

	return colIdx1, colIdx2, 0


def _normalisationWeights(data):
	nGenes = data.shape[0]
	nCells = data.shape[1]
	geneGM = gmean(data+1, 1)
	weights = array( [median( (data[:,i]+1) / geneGM) for i in range(nCells)] )
	return weights

# Read input data file, extract read counts, normalise, report errors
def readData(inputFile, label1, label2, normalise=True, removeZeros=False, useSpikeIns = False, verbose = False, ):

	data1 = []
	data2 = []
	spikeIns = []
	lineStatus = []
	ids = []
	empty = []

	header = inputFile.readline()

	colIdx1, colIdx2, status = _readHeader(header, label1, label2)

	if status == 1:
		lineStatus.append( Status(2, 'Header', 'Invalid header format') )
		return [],[],[],lineStatus
	else:
		if not colIdx1:
			lineStatus.append( Status(2, 'Header', "No colums with label '" + label1 + "' found") )
			return [],[],[],lineStatus
		elif not colIdx2:
			lineStatus.append( Status(2, 'Header', "No colums with label '" + label2 + "' found") )
			return [],[],[],lineStatus
		else:
			lineStatus.append( Status(0, 'Header', "Read OK") )

	for line in inputFile:
		if line.strip():
			cols = line.split()

			idx = cols[0]

			p1 = [ float(cols[x]) for x in colIdx1 ]
			p2 = [ float(cols[x]) for x in colIdx2 ]


			if max(p1) == 0 or max(p2) == 0:
				lineStatus.append( Status(1, idx,  "Null expression detected") )
				continue
			else:
				lineStatus.append( Status(0, idx,  "Line read OK") )


			if idx.lower().startswith('spike'):
				spikeIns.append( p1 + p2 )
				lineStatus.append( Status(0, idx,  "Spike-in detected") )
				continue

			data1.append(p1)
			data2.append(p2)
			ids.append(idx)

	if normalise:

		if useSpikeIns:
			if len(spikeIns) == 0:
				lineStatus.append(Status(2, 'Spike-ins','No spike-ins data detected') )
				return [],[],[],lineStatus

			weights = _normalisationWeights(array(spikeIns))
		
		else:
			weights = _normalisationWeights( hstack( ( array(data1), array(data2) ) ) )

		splitColumn = array(data1).shape[1]

		data1 = [ (array(readLine) / weights[:splitColumn]).tolist() for readLine in data1 ]
		data2 = [ (array(readLine) / weights[splitColumn:]).tolist() for readLine in data2 ]

	if removeZeros:

		dataFiltered1 = []
		dataFiltered2 = []
		idsFiltered = []

		for p1, p2, idx in zip(data1,data2,ids):
			p1 = filter(lambda x: x!=0, p1)
			p2 = filter(lambda x: x!=0, p2)

			if len(p1) != 0 and len(p2) != 0:
				dataFiltered1.append(p1)
				dataFiltered2.append(p2)
				idsFiltered.append(idx)
			else:
				lineStatus.append( Status(1, idx,  "Empty expression after zero removal") )

		data1 = dataFiltered1
		data2 = dataFiltered2
		ids = idsFiltered
	
	inputFile.close()

	return data1, data2, ids, lineStatus


# Sorts set x and returns a rank of elements in a sorted set
def _sortRank(x): 

	xs = sorted(x)
	r = []

	i = 1 
	while i <= len(xs):
		rStart = i
		rEnd = i

		if i < len(xs):
			while xs[i-1] == xs[i]:
				i += 1
				rEnd = i
				if i >= len(xs):
					break

		i += 1
		n = (rEnd-rStart + 1)
		r.extend([ (rStart + rEnd) / 2 for j in range(n)])

	return xs, r

# Perform a Cramer - von Mises test of two samples x and y. H0: samples x and y are drawn from the same distribution. Returns a p-value.
# Anderson, Theodore W., On the Distribution of the Two-Sample Cramer-von Mises Criterion, 1962, The Annals of Mathematical Statistics 33, 1148-1159
# Anderson, Theodore W., Donald A. Darling, 1952, Asymptotic Theory of Certain Goodness of Fit Criteria Based on Stochastic Processes, The Annals of Mathematical Statistics 23, 193-212.
def cramerVonMises(x, y):

	try:
		x = sorted(x);
		y = sorted(y);

		pool = x + y;

		ps, pr = _sortRank(pool)

		rx = array ( [ pr[ind] for ind in [ ps.index(element) for element in x ] ] )
		ry = array ( [ pr[ind] for ind in [ ps.index(element) for element in y ] ] )

		n = len(x)
		m = len(y)

		i = array(range(1, n+1))
		j = array(range(1, m+1))

		u = n * sum ( power( (rx - i), 2 ) ) + m * sum ( power((ry - j), 2) )

		t = u / (n*m*(n+m)) - (4*n*m-1) / (6 * (n+m))
		Tmu = 1/6 + 1 / (6*(n+m))
		Tvar = 1/45 * ( (m+n+1) / power((m+n),2) ) * (4*m*n*(m+n) - 3*(power(m,2) + power(n,2)) - 2*m*n) / (4*m*n)
		t = (t - Tmu) / power(45*Tvar, 0.5) + 1/6

		if t < 0:
			return -1
		elif t <= 12:
			a = 1-mp.nsum(lambda x : ( mp.gamma(x+0.5) / ( mp.gamma(0.5) * mp.fac(x) ) ) *
					 mp.power( 4*x + 1, 0.5 ) * mp.exp ( - mp.power(4*x + 1, 2) / (16*t) ) *
					  mp.besselk(0.25 , mp.power(4*x + 1, 2) / (16*t) ), [0,100] ) / (mp.pi*mp.sqrt(t))
			return float(mp.nstr(a,3))
		else:
			return 0
	except Exception as e:
		print e
		return -1



# Generation of a sample from Poisson-beta distribution with given parameters parms and size n
def randPoissonBeta(params,n):

	x = beta( params.alpha, params.beta, n )
	p = poisson( x * params.gamma )

	return p

# Estimation of the parameters of a sample drawn from a Poisson-beta distribution using momet matching technique from 
# Peccoud, Jean, Bernard Ycart, 1995, Markovian modelling of gene product synthesis, Theoretical Population Biology 48, 222-234.
def getParamsMoments(p):

	try:
		rm1 = sum(p) / len(p)
		rm2 = sum( [pow(x,2) for x in p] ) / len(p)
		rm3 = sum( [pow(x,3) for x in p] ) / len(p) 

		fm1 = rm1
		fm2 = rm2 - rm1
		fm3 = rm3 - 3*rm2 + 2*rm1

		r1 = fm1
		r2 = fm2 / fm1
		r3 = fm3 / fm2

		alpha = 2*r1 * (r3 - r2) / (r1*r2 - 2*r1*r3 + r2*r3)
		beta = 2 * (r2 - r1) * (r1 - r3) * (r3 - r2) / ((r1*r2 - 2*r1*r3 + r2*r3) *(r1 - 2*r2 + r3))
		gamma = (-r1*r2 + 2*r1*r3 - r2*r3) / (r1 - 2*r2 +r3)
	except:
		return Params(-1,-1,-1,-1)

	return Params(alpha, beta, gamma, 0)

# Esimation of the parameters of a sample drawn from a Poisson-beta distribution using bayesian inference method
# Kim, Jong Kyoung, John C. Marioni, 2013, Inferring the kinetics of stochastic gene expression from single-cell RNA-sequencing data, Genome Biology 14, R7.
def getParamsBayesian(p, iterN=1000):

	HyperParams = namedtuple("HyperParams", ["k_alpha", "theta_alpha", "k_beta", "theta_beta", "k_gamma", "theta_gamma"])

	parFit = getParamsMoments(p)

	hyperParams = HyperParams(k_alpha = 1, theta_alpha = 100, k_beta = 1, theta_beta = 100, k_gamma = 1, theta_gamma = max(p) )
	
	if parFit.alpha > 0 and parFit.beta > 0 and parFit.gamma > 0:

		params = Params(alpha = RVar(parFit.alpha), beta = RVar(parFit.beta), gamma = RVar(parFit.gamma), c = [ RVar(0.5) for i in range(len(p)) ] )	
	else:
		params = Params(alpha = RVar(0.5), beta = RVar(0.5), gamma = RVar(mean(p)+1e6), c = [ RVar(0.5) for i in range(len(p)) ] )

	bioParams = BioParams(size = RVar(), freq = RVar(), duty = RVar() )

	save = False

	for i in range(iterN):

		if i > iterN / 2:
			save = True

			alpha = params.alpha.value
			beta = params.beta.value
			gamma = params.gamma.value

			bioParams.size.sample.append( gamma / beta )
			bioParams.freq.sample.append( alpha )
			bioParams.duty.sample.append( alpha / (alpha + beta)  )

		for c,pi in zip(params.c,p):
			c.setSampleFunction(lambda x: (params.alpha.value - 1 ) * log(x) + (params.beta.value - 1) * log(1-x) + pi * log(x) - params.gamma.value * x)
			c.draw(saveToSample = save)

		params.gamma.setSampleFunction( lambda x: (hyperParams.k_gamma-1)*log(x) - x / hyperParams.theta_gamma + log(x) * sum(p) - x * sum( [c.value for c in params.c] ) )
		params.gamma.draw(saveToSample = save)
		
		params.alpha.setSampleFunction( lambda x: (hyperParams.k_alpha-1)*log(x) - x / hyperParams.theta_alpha + len(p)*(gammaln(x + params.beta.value) - gammaln(x)) + (x-1) * sum([ log(c.value) for c in params.c ]) )
		params.alpha.draw(saveToSample = save)

		params.beta.setSampleFunction( lambda x:  (hyperParams.k_beta-1)*log(x) - x / hyperParams.theta_beta + len(p)*(gammaln(x + params.alpha.value) - gammaln(x)) + (x-1) * sum([ log(1-c.value) for c in params.c ]) )
		params.beta.draw(saveToSample = save)

	return params, bioParams

def goodnessOfFit(p, params):

	try:
		alpha = params.alpha.mean()
		beta = params.beta.mean()
		gamma = params.gamma.mean()
	except:
		alpha = params.alpha
		beta = params.beta
		gamma = params.gamma

	pr = randPoissonBeta( Params(alpha,beta,gamma,0), len(p))

	return cramerVonMises(pr,p)




