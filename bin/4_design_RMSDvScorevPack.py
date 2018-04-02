# Marco Mravic DeGrado Lab Nov 2016
# Use this to align a directory full of rosetta membrane ab initio (low-res centroid) structures (folding trial) to some input model
# Default is C-alphas and only the TM regions from span file, but can change to backbone or full protein

# May also lead into structural clustering with mini-batch k-means

# python ../../../bin/4_folding_RMSDvScore_clustering.py ./ ../../looped_EEEH-1.pdb
from prody import *
import sys, os, numpy as np

inDir 	= sys.argv[1]
ref_pdb = parsePDB( sys.argv[2], subset='ca' )

#selStrDecoys = 'resnum 3 to 43'

# Look through the files, to parse score and align
trials, scores, rmsds, packstat = {}, [], [], []
print 'Aligning pdbs from', inDir, 'to', sys.argv[1]

for f in os.listdir( inDir ):
	path = os.path.join( inDir, f )

	if path[-3:] != 'pdb': continue

	# find score
	with open( path ) as fin:
		for i in fin:
			if i[:4] == 'pose':
				score = round( float( i.rsplit()[-1] ) , 2)
			if i[:4] == 'pack':
				pack = round( float( i.rsplit()[-1] ) , 2)
	try:
		score
	except NameError:
		print 'skipping scoreless file', f
		continue

	#print path, 
	# calculate best alignment of TM domains
	mobile 	= parsePDB( path, subset='ca' )
	seq = mobile.getSequence()[1:-1]
#	mobile	= pdb.select( selStrDecoys )
	superpose( mobile, ref_pdb )
	rmsd = round( calcRMSD( mobile, ref_pdb ), 2)
#	print rmsd, score, pack


	trials[f] = (rmsd,score,pack, seq)
	scores.append( score )
	rmsds.append( rmsd )

scores = np.array( scores )
rmsds = np.array(  rmsds )
print 
import matplotlib.pyplot as plt


seqs = []
for k,v in sorted( trials.items(), reverse=True, key=lambda x: x[1][2] ):
	if v[2] > 0.69 and v[1] <-105:

		print k, v[0], v[1], v[2], v[3]
		if v[3] not in seqs:
			seqs.append( v[3] )
for k in seqs:
	print k


plt.plot(rmsds,scores, 'ro')
plt.show()

