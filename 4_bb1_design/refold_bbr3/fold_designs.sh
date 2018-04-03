#!/usr/bin/python
#$ -S /bin/bash
#$ -l mem_free=1G
#$ -l arch=linux-x64
#$ -l netapp=1G
#$ -l h_rt=00:20:00
#$ -cwd
#$ -j y
#$ -o /netapp/home/mmravic/miniLig/4_bb1_design/refold_bbr3/logs
#$ -t 1-2
pwd
$ROSETTA/source/bin/AbinitioRelax.linuxgccrelease -database /netapp/home/mmravic/bin/Rosetta/database/ -in:file:fasta bbr3_seq.txt -in:file:frag3 frag_3 -in:file:frag9 frag_9 -abinitio:relax True -abinitio::increase_cycles 10 -abinitio::rg_reweight 0.5 -abinitio::rsd_wt_helix 0.5 -abinitio::rsd_wt_loop 0.5 -use_filters true -nstruct 1 -relax::fast -out:overwrite -out:file:silent outputs/model_0000001.out -out:path outputs/ -psipred_ss2 psipred-ss2 -rebuild_disulf true -detect_disulf true
$ROSETTA/source/bin/AbinitioRelax.linuxgccrelease -database /netapp/home/mmravic/bin/Rosetta/database/ -in:file:fasta bbr3_seq.txt -in:file:frag3 frag_3 -in:file:frag9 frag_9 -abinitio:relax True -abinitio::increase_cycles 10 -abinitio::rg_reweight 0.5 -abinitio::rsd_wt_helix 0.5 -abinitio::rsd_wt_loop 0.5 -use_filters true -nstruct 1 -relax::fast -out:overwrite -out:file:silent outputs/modelSS_0000001.out -out:path outputs_SS/ -psipred_ss2 psipred-ss2 -rebuild_disulf true -detect_disulf true -in:fix_disulf disulfide.txt
for i in ./logs/*; do gzip $i; done
