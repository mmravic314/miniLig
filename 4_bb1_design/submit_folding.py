import os, sys, subprocess as sp 

qsub_header = '''#!/usr/bin/python

#$ -S /usr/bin/python
#$ -l mem_free=1G
#$ -l arch=linux-x64
#$ -l netapp=1G
#$ -l h_rt=00:20:00
#$ -cwd
#$ -j y
#$ -o /netapp/home/mmravic/
#$ -t 1-2
'''

rosiDB 		= os.path.join( sys.argv[2], 'database/' )
#rosiFOLD	= os.path.join( sys.argv[2], 'source/bin/AbinitioRelax.macosclangrelease' )
rosiFOLD       = os.path.join( '$ROSETTA/', 'source/bin/AbinitioRelax.linuxgccrelease' ) 
#rosiSCORE      = os.path.join( sys.argv[1], 'source/bin/score.macosclangrelease' )
#rosiSCORE	= os.path.join( sys.argv[1], 'source/bin/score.linuxgccrelease' ) 

for wrkDir in [ x for x in os.listdir( sys.argv[1] ) if 'refold' in x  ]:
	qsub_header = '''#!/usr/bin/python
#$ -S /bin/bash
#$ -l mem_free=1G
#$ -l arch=linux-x64
#$ -l netapp=1G
#$ -l h_rt=00:20:00
#$ -cwd
#$ -j y
#$ -o /netapp/home/mmravic/miniLig/4_bb1_design/%s/logs
#$ -t 1-75000
pwd
''' % wrkDir



	inFasta   =  '%s_seq.txt' % os.path.basename( wrkDir )[7:] 
	disulfide = 'disulfide.txt' 
	outputsSS = 'outputs_SS/'
	outputs   = 'outputs/' 
	qsub_cmdTx= qsub_header

	add_nulls = lambda number, zero_count : "{0:0{1}d}".format(number, zero_count)
	try:
		output_suffix 			= add_nulls( os.environ["SGE_TASK_ID"], 7 ) 
	except KeyError:
		output_suffix                      = add_nulls( 1, 7 )

	sil_outSS = os.path.join( outputsSS, 'model_%s.out' % output_suffix )
	pdb_outSS = os.path.join( outputsSS, 'model_%s.out' % output_suffix )
	sil_out   = os.path.join( outputs, 'model_%s.out' % output_suffix )
	sil_outSS = os.path.join( outputs, 'modelSS_%s.out' % output_suffix )
	pdb_out   = os.path.join( outputs, 'model_%s.out' % output_suffix )

	print '\n***', wrkDir, '***\n'
#	inFasta = 'fasta'
	os.chdir( wrkDir )
#	inFasta = os.path.join( wrkDir, inFasta )
	#print open(os.path.join( wrkDir, inFasta ) ).read()
#	print os.listdir( '.' )

	##### Now check for a manually generated flag file
	#		 to see if this directory has already been processed
	if 'flag.txt' in os.listdir( '.' ):
		continue
	#  check this dir is ready to process, skip if files missing 
	if 'frag_3' not in os.listdir( '.' ):
		os.chdir( sys.argv[1] )
		continue


	noSSrefold_cmd = [ rosiFOLD, 
'-database', 				rosiDB, 
'-in:file:fasta', 			inFasta,
'-in:file:frag3', 			'frag_3',
'-in:file:frag9', 			'frag_9',
'-abinitio:relax', 			'True',
'-abinitio::increase_cycles', '10',
'-abinitio::rg_reweight', 	'0.5',
'-abinitio::rsd_wt_helix', 	'0.5',
'-abinitio::rsd_wt_loop', 	'0.5',
'-use_filters', 		'true', 
'-nstruct',			'1',
'-relax::fast',
'-out:overwrite', 
'-out:file:silent', 		sil_out,
'-out:path',                    outputs,
#'-out:pdb',                    pdb_out,
#'-psipred_ss2',					ss2F,
'-psipred_ss2',					'psipred-ss2',
'-rebuild_disulf',              'true', 
'-detect_disulf',               'true', 
#'-in:fix_disulf',               disulfide, 
]

	print  '\n', noSSrefold_cmd, '\n'
	qsub_cmdTx += ' '.join( noSSrefold_cmd ) + '\n'
	# Run folding
#	sp.call( noSSrefold_cmd )


	SSrefold_cmd = [ rosiFOLD, 
'-database', 				rosiDB, 
'-in:file:fasta', 			inFasta,
'-in:file:frag3', 			'frag_3',
'-in:file:frag9', 			'frag_9',
'-abinitio:relax', 			'True',
'-abinitio::increase_cycles', '10',
'-abinitio::rg_reweight', 	'0.5',
'-abinitio::rsd_wt_helix', 	'0.5',
'-abinitio::rsd_wt_loop', 	'0.5',
'-use_filters', 		'true', 
'-nstruct',			'1',
'-relax::fast',
'-out:overwrite', 
'-out:file:silent', 		sil_outSS,
'-out:path',                    outputsSS,
#'-out:pdb',                    pdb_out,
'-psipred_ss2',					'psipred-ss2',
'-rebuild_disulf',              'true', 
'-detect_disulf',               'true', 
'-in:fix_disulf',               disulfide, 
]
	print SSrefold_cmd
#	sp.call( SSrefold_cmd )
	qsub_cmdTx += ' '.join( SSrefold_cmd ) + '\n'
	qsub_cmdTx +='for i in ./logs/*; do gzip $i; done\n'
	print qsub_cmdTx
	qsub_file = open( 'fold_designs.sh', 'w' )
	qsub_file.write(qsub_cmdTx)
	qsub_file.close()
	


	os.chdir( sys.argv[1] )
		
