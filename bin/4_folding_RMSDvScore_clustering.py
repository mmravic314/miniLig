# Marco Mravic DeGrado Lab Nov 2016
# Use this to align a directory full of rosetta membrane ab initio (low-res centroid) structures (folding trial) to some input model
# Default is C-alphas and only the TM regions from span file, but can change to backbone or full protein

# May also lead into structural clustering with mini-batch k-means

# python ../../../bin/4_folding_RMSDvScore_clustering.py ./ ../../looped_EEEH-1.pdb
from prody import *
import sys, os, numpy as np

inDir 	= sys.argv[1]
ref_pdb = parsePDB( sys.argv[2], subset='ca' )

selStrDecoys = 'resnum 3 to 43'

# Look through the files, to parse score and align
trials, scores, rmsds = {}, [], []
print 'Aligning pdbs from', inDir, 'to', sys.argv[1]

for f in os.listdir( inDir ):
	path = os.path.join( inDir, f )

	# find score
	with open( path ) as fin:
		for i in fin:
			if i[:4] == 'pose':
				score = round( float( i.rsplit()[-1] ) , 2)
	try:
		score
	except NameError:
		print 'skipping scoreless file', f
		continue

	# calculate best alignment of TM domains
	pdb 	= parsePDB( path, subset='ca' )
	mobile	= pdb.select( selStrDecoys )
	superpose( mobile, ref_pdb )
	rmsd = round( calcRMSD( mobile, ref_pdb ), 2)
	print path, rmsd, score


	trials[f] = (rmsd,score)
	scores.append( score )
	rmsds.append( rmsd )

scores = np.array( scores )
rmsds = np.array(  rmsds )
print 
import matplotlib.pyplot as plt



for k,v in sorted( trials.items(), key=lambda x: x[1][1] ):
	print k, v[0], v[1]


plt.plot(rmsds,scores, 'ro')
plt.show()

