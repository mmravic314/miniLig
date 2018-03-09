#For match file from master, look up helices and find if they are connected... find if loop is connecting segment

import sys, os, subprocess as sp, numpy as np, shutil 
from prody import *
from operator import itemgetter
from itertools import groupby
from PDButil import UnNatAA
from collections import defaultdict, Counter

inPDB = parsePDB( sys.argv[1], subset = 'bb' )

# python ~/bin/loopFinder.py 56_helices.pdb 56_helices.m ~/termanal/support.default/151218_masterDB_parsedPDB/ ~/tertBuilding/CMP_bobo/56_helices/
# python ~/bin/loopFinder.py BA-Loop.pdb BA-Loop.m ~/termanal/support.default/v2_162901_bc_30-scPDB_oPDB/ BA-Loop_Matches/


seqs = []
loops = defaultdict(list)

# look into file
lineNum = 0
with open( sys.argv[2] ) as file:
	for m in file:
		
		rmsd = float(  m[:7].strip() )
		if rmsd > 0.9: 
			break

		fragMatch = [ tuple( x.strip(',').strip('()').split(',') ) for x in m.split('[')[-1].split(']')[0].split() ]

		gapL = int ( fragMatch[1][0] ) - int( fragMatch[0][1] ) -1

		if gapL > 10 or gapL < 0:
			lineNum +=1
			continue

		#Convert match residue index ranges to array of all ( for prody selection, e.g. [38,39,40,41,59,60,61,62] )
		all_match = []
		for k in fragMatch:
			rng = np.arange( int( k[0] ), 1 + int( k[1] )  )
			all_match.extend( rng )

		resis = ' '.join( [ str( x ) for x in sorted( all_match[:] )] )


		selStr 		= 'bb resnum ' + resis
		selRng 		= 'resnum %d to %d' % ( all_match[0], all_match[-1] )
		look_upPath = os.path.join( sys.argv[3], os.path.basename( m.split()[1] ).split('.')[0] + '.pdb' )	
		mPdb 		= parsePDB( look_upPath )

		r = 0
		for res in mPdb.iterResidues():
			res.setResnum( r )
			r += 1


		wholeset 	= mPdb.select( selRng ).copy()
		ends		= mPdb.select( selStr ).copy()

		seq =  ''.join( [ UnNatAA[ x.getResname() ]  for x in wholeset.iterResidues() ] )

		if seq in seqs: continue
		seqs.append( seq )

		print 'loop Length:', gapL, 'RMSD:', rmsd, 'match%d' % ( lineNum ),# 'residue indices in match:', fragMatch
		print	seq
		print

		wholeset.setTitle( str(rmsd) )
		ends.setTitle( str(rmsd) )
		loopPath = os.path.join( sys.argv[4], 'loopMatch' + str(lineNum) + '.pdb' )
		endsPath = os.path.join( sys.argv[4], 'endsMatch' + str(lineNum) + '.pdb' )

		# find transformation matrix from original PDB (aligned to target) to native PDB
		betterMat	 = superpose( ends, inPDB )[1]

		wholeset 	 = applyTransformation( betterMat, wholeset.copy())
		wholeset.setTitle(loopPath)

		writePDB( loopPath, wholeset )
		writePDB( endsPath, ends )
		## Align loop file to input PDB coords

		loops[gapL].append(wholeset.copy())

		lineNum +=1

# find biggest cluster
maxV, bestLen = 0, 0
for loop_len, pdbs in sorted( loops.items() ):
	print loop_len, len(pdbs)
	if len(pdbs) > maxV:
		maxV, bestLen = len(pdbs), loop_len

# HACK to also do a specific loop number:
bestLen = 6

## now cluster biggest cluster...
# finding a good cut off may be difficult
from scipy.cluster.hierarchy import linkage, fcluster
from itertools import combinations

print 'most common loop:', bestLen, '....now clustering...'



ssd_array 	= []
for n, pair in enumerate( combinations( loops[bestLen], 2 ) ):
		pair = pair[0].select('bb').copy(), pair[1].select('bb').copy()
		ssd_array.append( calcRMSD( *pair ) )

ssd_array 	= np.array( ssd_array )
cutoff		= 1.8
linkMat 	= linkage( ssd_array , method='complete', metric='euclidean')
h_clust 	= fcluster( linkMat, cutoff, criterion='distance')
numClust 	= len( set(h_clust) )
print 'at cutoff %.2f Angstrom, Unique clusters found:' % cutoff, numClust, '\n'

clust_cnts = tuple( Counter( h_clust).most_common()[:2] )
c1, c2, c1_Percent, c2_percent  = clust_cnts[0][0], clust_cnts[1][0], 100*float(clust_cnts[0][1] )/maxV, 100*float( clust_cnts[1][1] )/maxV
print 'top 2 clusters (percent): %d (%0.1f); %d (%0.1f)' % (c1, c1_Percent, c2, c2_percent)
print 'loops per cluster:',  Counter( h_clust).most_common()



# cluster 1 & 2 fetch and put into folder
cluster1, cluster2 = [], []
for n in np.arange(len(h_clust)):

	index = h_clust[n]
	if index == c1:
		cluster1.append( loops[bestLen][n] )

	if index == c2:
		cluster2.append( loops[bestLen][n] )

print 'moving cluster members to new dir for c1 & c2, +fastas'
## now move these loops to a newDir
c1_dir 		= os.path.join(sys.argv[4], 'cluster1-l%d_%d' % (bestLen,clust_cnts[0][1] ) )

if not os.path.exists( c1_dir ):
	os.mkdir(c1_dir)
c1_fasta 	= os.path.join(c1_dir , 'cluster1_fasta.txt' )
c1_fasta_txt = ''
for i in cluster1:
	oldPath = i.getTitle()
	newPath = os.path.join( c1_dir, os.path.basename( i.getTitle() ) )
	seq = i.select('ca').getSequence()
	c1_fasta_txt += '>>%s\n%s\n' % (oldPath, seq)
	shutil.copy( oldPath, newPath )
oFile = open( c1_fasta, 'w')
oFile.write(c1_fasta_txt)

## now move these loops to a newDir
c2_dir 		= os.path.join(sys.argv[4], 'cluster2-l%d_%d' % (bestLen,clust_cnts[0][1] ) )

if not os.path.exists( c2_dir ):
	os.mkdir(c2_dir)
c2_fasta 	= os.path.join(c2_dir, 'cluster2_fasta.txt' )
c2_fasta_txt = ''
for i in cluster2:
	oldPath = i.getTitle()
	newPath = os.path.join( c2_dir, os.path.basename( i.getTitle() ) )
	seq = i.select('ca').getSequence()
	c2_fasta_txt += '>>%s\n%s\n' % (oldPath, seq)
	shutil.copy( oldPath, newPath )
oFile = open( c2_fasta, 'w')
oFile.write(c2_fasta_txt)






