#!/bin/csh
# wrapper to run MNE-Python script in SLURM/sbatch
# Version 0.18 in Python 3 on NEW cluster

# for example call:
# sbatch -o wrapper_test.out 
# --export=pycmd="FPVS_filter_raw.py",subj_idx=13

# echo "PBS_ARRAYID " $PBS_ARRAYID
# echo $myvar

# activate conda environment for corresponding version
# of MNE-Python
echo "The fun begins here."

echo "Node number name host:"
echo $SLURM_NODEID $SLURMD_NODENAME $SLURM_SUBMIT_HOST

# show conda version
echo "Conda version:"
conda --version

conda activate mne0.19

# show path to conda environment
echo "Conda environment:"
echo $CONDA_PREFIX

# show where python executable is located
echo "Python location:"
which python

# pycmd: python script to execute
# subj_idx: subject index to be processed
# var: variable string (e.g. for filenames)

# /imaging/local/software/miniconda/envs/mne0.18/bin/python $pycmd $subj_idx
python $pycmd $subj_idx