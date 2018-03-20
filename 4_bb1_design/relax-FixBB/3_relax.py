### python 3_relax.py ~/rosetta/ pdL1-Mini1_EEEH.pdb  3_relaxTHIN.xml Resfile.txt  3_relax-Mini1_EEEH.cst relax_outputs/

# 1) path 2 rosetta main
# 2) path to input structure file
# 3) path to XML script for protocol
# 4) resfile
# 5) constraint file for backbone atom minimization
import sys, os, subprocess as sp, re


################## MAIN #######################
# Non-variable args
rosetta_database_path   = os.path.join( sys.argv[1] , 'database/' )
rosetta_scriptsEXE_path = os.path.join( sys.argv[1], 'source/bin/rosetta_scripts.macosclangrelease' )
design_script_path      = sys.argv[3]
struc_path 				= sys.argv[2]
resfile_path			= sys.argv[4]
cst_path				= sys.argv[5]

# Variable args

output_prefix			= sys.argv[6]
try:
	output_suffix 			= '_out%s' % (str(  os.environ["SGE_TASK_ID"]) )
except KeyError:
	output_suffix                      = '_out%s' % ( 'Local' )

cmd = [
		rosetta_scriptsEXE_path,
		'-database', rosetta_database_path,
		'-parser:protocol', design_script_path,
		'-in:file:s', struc_path,
		'-nstruct', '2',
		'-out:prefix', output_prefix,   
		'-out:suffix', output_suffix,                               
		'-out:no_nstruct_label',
		'-out:overwrite',
		'-use_input_sc',
        '-packing:resfile', resfile_path,
		'-parser:script_vars', 'cst_file=%s' % ( cst_path )
]

print
print cmd
print

sp.call( cmd )
