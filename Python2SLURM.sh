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

conda activate mne0.20

# show path to conda environment
echo "Conda environment:"
echo $CONDA_PREFIX

# show where python executable is located
echo "Python location:"
which python

setenv MESA_GL_VERSION_OVERRIDE 3.3

# pycmd: python script to execute
# subj_idx: subject index to be processed
# var: variable string (e.g. for filenames)

# python $pycmd $subj_idx

# from: https://docs.enthought.com/mayavi/mayavi/tips.html
# also: https://patricksnape.github.io/2014/offscreen_rendering/
echo "Using xvfb-run"
xvfb-run -d --server-args="-screen 0 1024x768x24" python $pycmd $subj_idx