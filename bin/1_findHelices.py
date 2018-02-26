import sys, os, numpy as np
from prody import *
from collections import defaultdict

#

# input 1: match file
# input 2: path to database pdb dir
# input 3: query PDB
# input 4: path the dir with match structures e.g. 'match00001.pdb'

matchFile 	= sys.argv[1]
db_dir 		= sys.argv[2]
queryFile 	= sys.argv[3]
matchDir 	= sys.argv[4]

matchHelix  = {}

query 		= parsePDB( queryFile )
query_seq 	= query.select('ca').getSequence()

probeAtom 	= query.select('name SG resnum 114').copy()

# query positions

#### functions & stuff #####

def list2selStr( array ):
	return 'resnum %s' % ' '.join( [str(int(x)) for x in array] )

def findCorrectRanges( seg1_rng, seg2_rng, opdb, match_seq ):
	
	match_1, match_2 	= match_seq[:5], match_seq[5:]
	seg_1seq, seg2_seq 	= opdb.select( 'ca resnum %s' % seg1_rng).getSequence(), opdb.select( 'ca resnum %s' % seg2_rng).getSequence()

	#print match_1, match_2
	# sample and fix resi rng 1
	if seg_1seq == match_1:
		seg1_flg = 1
	else:
		seg1_flg = 0

		for i in np.arange( 1, 31 ):
			baseRng = tuple([ int(x) for x in seg1_rng.split(' to ') ] )
			fwdRng 	= ' '.join( [ str(x) for x in np.arange( baseRng[0], baseRng[1] + 1) + i ] )
			fwdObj	= opdb.select( 'ca resnum %s' % fwdRng)

			if baseRng[0] <= i:
				rvRng = ''
				rvObj = 0
			else:
				rvRng 	= ' '.join( [ str(x) for x in np.arange( baseRng[0], baseRng[1] + 1) - i ] )
				rvObj 	= opdb.select('ca resnum %s' % rvRng)

			try:
				fwdSeq = fwdObj.getSequence()
			except ValueError:
				fwdSeq = ''
			try:
				rvSeq  = rvObj.getSequence()
			except (ValueError, AttributeError):
				rvSeq = ''

			if fwdSeq == match_1:
				seg1_rng, seg_1seq = fwdRng, fwdSeq
				seg1_flg += 1
				break
			if rvSeq == match_1:
				seg1_rng, seg_1seq = rvRng, fwdSeq
				seg1_flg += 1
				break

	if seg2_seq == match_2:
		seg2_flg = 1
	else:
		seg2_flg = 0

		for i in np.arange( 1, 31 ):
			baseRng = tuple([ int(x) for x in seg2_rng.split(' to ') ] )
			fwdRng, rvRng = ' '.join( [ str(x) for x in np.arange( baseRng[0], baseRng[1] + 1) + i ] ), ' '.join( [ str(x) for x in np.arange( baseRng[0], baseRng[1] + 1) - i ] )
			fwdRng 	= ' '.join( [ str(x) for x in np.arange( baseRng[0], baseRng[1] + 1) + i ] )
			fwdObj	= opdb.select( 'ca resnum %s' % fwdRng)

			if baseRng[0] <= i:
				rvRng = ''
				rvObj = 0
			else:
				rvRng 	= ' '.join( [ str(x) for x in np.arange( baseRng[0], baseRng[1] + 1) - i ] )
				rvObj 	= opdb.select('ca resnum %s' % rvRng)

			try:
				fwdSeq = fwdObj.getSequence()
			except ValueError:
				fwdSeq = ''
			try:
				rvSeq  = rvObj.getSequence()
			except (ValueError, AttributeError):
				rvSeq = ''

	#		print 2, fwdSeq, rvSeq, fwdRng, rvRng

			if fwdSeq == match_2:
				seg2_rng, seg2_seq = fwdRng, fwdSeq
				seg2_flg += 1
				break
			if rvSeq == match_2:
				seg2_rng, seg2_seq = rvRng, rvSeq
				seg2_flg += 1
				break

#	print seg1_seq, seg2_seq, seg1_flg , seg2_flg, 'fixed', '\n'
	if seg1_flg and seg2_flg:
		return 1, seg1_rng, seg2_rng
	else:
		return 0, 0, 0



############ end #############


################## MAIN ##################

# skip repeat matches... of sequences already seen...
seqs_visited 	= defaultdict(list)
helix_DistData, helix_LenData 	= {}, {} 
helixDir 		= 'located_helices_8res'
helixDir_full	= 'located_helices_FULL'

if not os.path.exists( helixDir ):
	os.mkdir( helixDir )

if not os.path.exists( helixDir_full ):
	os.mkdir( helixDir_full )

for n, line in enumerate( open(matchFile), 1 ):

	# find match pdb and skip if sequence is has already been observed, hash if not
	index 		= "%05d"% (n)
	match_pdb 	= parsePDB( os.path.join( matchDir, 'match%s.pdb' % index) )
	match_seq 	= match_pdb.select('ca').getSequence()

	if match_seq in seqs_visited[match_seq[:3]]:
			continue
	else:
			seqs_visited[match_seq[:3]].append(match_seq)

	# parse data in lines of match file
	line 		= line.split()
	rmsd, opdb_path, seg1_rng, seg2_rng = float(line[0]), line[1][-10:-1]+'b', line[2][2:-2].split(','), line[3][1:-2].split(',')

	# extract matching fragment from original PDB... but renumber it first
	opdb_path 	= os.path.join( db_dir, opdb_path )
	seg1_rng, seg2_rng = ' to '.join(seg1_rng), ' to '.join(seg2_rng)

	opdb 		= parsePDB( opdb_path )
	stResi 		= opdb.getResnums()[0]
	cnt			= 0
	for resi in opdb.iterResidues():
		#resi.setResnum( resi.getResnum() - stResi )
		resi.setResnum(cnt)
		cnt+=1

	
	seg1_seg2 = opdb.select( 'resnum %s' % (seg1_rng) ) + opdb.select( 'resnum %s' % (seg2_rng) )
	seg1_seq, seg2_seq = opdb.select( 'ca resnum %s' % (seg1_rng) ).getSequence(),  opdb.select( 'ca resnum %s' % (seg2_rng) ).getSequence()
	#print index, opdb_path, seg1_rng, seg2_rng

	# check if the match is correctly indexed to original PDB, by sequnce; correct if not...
	if match_seq != seg1_seg2.select('ca').getSequence():	
		#print index, opdb_path, seg1_rng, seg2_rng, '\n'
		#print seg1_seq, seg2_seq
		fix_flg, seg1_rng, seg2_rng = findCorrectRanges( seg1_rng, seg2_rng, opdb, match_seq )
		if not fix_flg:
			print 'ERROR! Sequence could not be found in this PDB\n\n\n'
			continue
		seg1_seg2 = opdb.select( 'resnum %s' % (seg1_rng) ) + opdb.select( 'resnum %s' % (seg2_rng) )
		seg1_seq, seg2_seq = opdb.select( 'ca resnum %s' % (seg1_rng) ).getSequence(),  opdb.select( 'ca resnum %s' % (seg2_rng) ).getSequence()
		


	# transform the opdb to match with the query... then search for residues nearby Cys S
	superpose( seg1_seg2.select('bb'), query.select('bb') )
	nrmsd = calcRMSD( seg1_seg2.select('bb'), query.select('bb') )
	if nrmsd - float(rmsd) > 0.8 : 
		continue


	contacts 			= opdb.select('(same residue as within 7.5 of target) and (not resnum %s %s )' % (seg1_rng,seg2_rng), target=probeAtom)
	if contacts:
		resnums_contacting  = contacts.select('ca').getResnums()
	else:
		continue

	helical_resnums 	= []
	not_helical 		= []

	# find which of these are helical 
	for i in resnums_contacting:
		resi = opdb[i,]
		try:
			phi, psi = calcPhi( resi ), calcPsi( resi )
			if -75 < phi < -45 and -60 < psi < -30:
				helical_resnums.append( i )
			else:
				not_helical.append( i )
		except (ValueError, TypeError ):			# if residue doesn't exist, is at termini, or disconnected, or anything else wierd... skip
				not_helical.append( i )

	if not len(helical_resnums):
		continue

	#print 'helix found nearby', match_seq ,  helical_resnums
	### add 5 residues on each side
	new_segments = []
	for i in helical_resnums:
		resiRng = np.append( np.append(  np.arange( i-5, i ), np.array([i]) ),  np.arange( i+1, i + 6 ) )
		new_segments = np.append( new_segments, resiRng )

	#print 'check these also', list( set(new_segments) )

	helical_resnums = []
	for i in list( set(sorted(new_segments)) ):
		if i in not_helical:
			continue
		resi = opdb[int(i),]
		try:
			phi, psi = calcPhi( resi ), calcPsi( resi )
			if -80 < phi < -40 and -80 < psi < -40:
				helical_resnums.append( i )
		except (ValueError, TypeError):			# if residue doesn't exist, is at termini, or disconnected, or anything else wierd... skip
			pass

	#print 'nearby helical residues', helical_resnums
	if len(helical_resnums) < 8:
		continue


	helical_resnums = sorted(helical_resnums)
	# find continuous segments 8-residues or more
	segments 	= []
	dyn_segList = [ helical_resnums[0] ]
	for step in np.arange( 1 ,len( helical_resnums )   ) :
	#	print step, helical_resnums[step] ,  helical_resnums[step-1], helical_resnums[step] - helical_resnums[step-1]

		if  helical_resnums[step] - helical_resnums[step - 1]  == 1:
			dyn_segList.append( helical_resnums[step] )
#			if step == len( helical_resnums ) - 2:
#				dyn_segList.append( helical_resnums[step + 1] )

		else:
			if len( dyn_segList ) >= 8:			# if long enough segment found, save it and reset queue
				segments.append( dyn_segList )
				dyn_segList = [ helical_resnums[step] ]

			else: 
				dyn_segList = [ helical_resnums[step] ]


	# if queue at the end of loop is long enough, save it
	if len(dyn_segList) >= 8:			
		segments.append( dyn_segList )

	if len (segments) < 1:
		continue

	print index, opdb_path, seg1_rng, seg2_rng
	print '>> bonafide continuous helices (n>7)', segments

	
	## from each segment, extract the 8 residue span 
	#  with the closest mean CA-CA distance to the sheet CA's 
	query_CA = query.select('ca')
	segN = 0
	for seg in segments:
		maxDist, res_set, seg_8 = 0, 0, seg[:8]
		print seg, 'og', len(seg)
		for i in np.arange( len(seg ) - 7):
			
			segi = opdb.select( 'ca %s' % list2selStr( seg_8 + i ) )
			# calc mean inverse distance between backbones 
			meanInvDist =  round( np.mean( 1/buildDistMatrix( segi, query_CA )), 5 )
#			print i + seg, list2selStr(seg + i), meanInvDist
			if meanInvDist	> maxDist:
				maxDist = meanInvDist
				res_set = i

		# set a distance cut-off so we don't look at poorly interacting helices. 
		if maxDist < 0.08: continue

		res_set 	= res_set + seg_8
		segi 		= opdb.select( '%s' % list2selStr(res_set) )
		all_helix 	= opdb.select( '%s' % list2selStr(seg) )

		helix_name = 'm%s_%d.pdb' % (index, segN)

		helix_DistData[ helix_name ] = maxDist
		helix_LenData[ helix_name ]	 = len(seg)
		writePDB( os.path.join(helixDir, helix_name), segi )
		writePDB( os.path.join(helixDir_full, 'full_'+helix_name), all_helix )

		print helix_name, segi.select( 'ca' ).getSequence(), maxDist, '\n'
		segN += 1




print np.mean( helix_DistData.values() ), np.std( helix_DistData.values() ), 'Helices Saved: ', len(helix_DistData.values()), '\n'

print '# helix_ID | mean_inverse_distance_2_target | length_helix_in_parent'
for k,v in sorted( helix_DistData.items() , key = lambda x: int( x[0][1:6]) ):
	print k, v, helix_LenData[k]




