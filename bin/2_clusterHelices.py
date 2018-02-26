import sys, os, numpy as np
from prody import *
from scipy.cluster.hierarchy import linkage, fcluster
from itertools import combinations
from collections import Counter, defaultdict
from Bio.Cluster.cluster import kmedoids

helixDir 		= sys.argv[1]
rmsd_mat_Store 	= sys.argv[2]

 
helixDir = [os.path.join( helixDir, x) for x in os.listdir(helixDir) if 'pdb' in x]

helix_array = []

for f in helixDir:
	pdb = parsePDB( f, subset='bb' )
	helix_array.append( pdb )

ssd_array= []
for n, pair in enumerate( combinations( helix_array, 2 ) ):
	print calcRMSD( *pair )
