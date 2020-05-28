#!/imaging/local/software/miniconda/envs/mne0.20/bin/python
"""
Prepare Source Spaces for FPVS Frequency Sweep.

==========================================

OH, January 2020
"""

# TO DO: turn into function, use info from config file for FPVS

import sys

import os
from os import path as op
import numpy as np

from copy import deepcopy

# import matplotlib
# matplotlib.use('Agg') #  for running graphics on cluster ### EDIT

from matplotlib import pyplot as plt

from importlib import reload

import mne

import config_sweep as config
reload(config)

print(mne.__version__)


filename = "/imaging/local/software/mne_python/set_MNE_2.7.3_FS_6.0.0_environ.py"
# for Python 3 instead of execfile
exec(compile(open(filename, "rb").read(), filename, 'exec'))

subjects_dir = config.subjects_dir


def make_source_space(sbj_id):
    """ Make Source Space. """

    subject = config.mri_subjects[sbj_id]

    if subject == '':

        print('No subject name for MRI specified - doing nothing now.')

        return

    ### set up source space
    src = mne.setup_source_space(subject, spacing=config.src_spacing,
                                 subjects_dir=subjects_dir,
                                 add_dist=True)

    # check if bem directory exists, and create it if not
    bem_path = op.join(subjects_dir, subject, 'bem')
    if not os.path.isdir(bem_path):
        print('Creating directory %s.' % bem_path)
        os.mkdir(bem_path)

    src_fname = op.join(subjects_dir, subject, 'bem', subject + '_' + str(config.src_spacing) + '-src.fif')

    print("###\nWriting source spaces to " + src_fname + "\n###")
    mne.write_source_spaces(src_fname, src, overwrite=True)


# get all input arguments except first
if len(sys.argv) == 1:

    sbj_ids = np.arange(0, len(config.map_subjects)) + 1

else:

    # get list of subjects IDs to process
    sbj_ids = [int(aa) for aa in sys.argv[1:]]


for ss in sbj_ids:

    make_source_space(ss)

print('Done.')
