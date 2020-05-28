#!/imaging/local/software/miniconda/envs/mne0.20/bin/python
"""
=========================================================
Make and plot sensitivity maps for FPVS.

Saving figures doesn't work on cluster yet.
run FPVS_SensitivityMaps.py <config.do_subjs>
=========================================================

"""
from time import sleep

import os
# needed to run on SLURM
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

sleep(1)
print('00\n')

import sys
from os import path as op

from copy import deepcopy

import numpy as np

print('0a')
from xvfbwrapper import Xvfb
vdisplay = Xvfb(width=1920, height=1080)
vdisplay.start()

print('0')

from mayavi import mlab

mlab.options.offscreen = True

print('0a')
import matplotlib
matplotlib.use('Agg')  # possibly for running on cluster

from importlib import reload

from time import sleep

import mne

import config_sweep as config
reload(config)

# whether to morph the STCs or not
do_morph = 1

subjects_dir = config.subjects_dir

# where figures will be written to
bem_dir = '/group/erp/data/olaf.hauk/MEG/FPVS/data/MRI/BEM_figs'

ch_types = ['grad', 'mag', 'eeg']

print('1')

def run_make_sensitivity_maps(sbj_id):
    """Compute sensitivity maps for one subject.

    Plot to figure.
    Return dictionary with STCs for ch_types.
    """

    subject = config.mri_subjects[sbj_id]

    if subject == '':

        print('No subject name for MRI specified - doing nothing now.')

        return

    print('Making Forward Solution for %s.' % subject)

    sbj_path = op.join(config.data_path, config.map_subjects[sbj_id][0])

    fwd_fname = op.join(sbj_path, subject + '_EEGMEG-fwd.fif')
    print('Reading forward solution: %s.' % fwd_fname)

    fwd_eegmeg = mne.read_forward_solution(fwd_fname)

    maps = {}  # will contains STCs
    for ch_type in ch_types:

        # channel types will be picked
        fwd = deepcopy(fwd_eegmeg)

        print('2')

        maps[ch_type] = mne.sensitivity_map(fwd=fwd, ch_type=ch_type, mode='free')

        if do_morph:

            print('3')

            # morph STCs for group averaging
            print('Computing morphing matrix.')
            morph_mat = mne.compute_source_morph(src=maps[ch_type], subject_from=subject,
                                                 subject_to='fsaverage', subjects_dir=subjects_dir)

            map_morph = morph_mat.apply(maps[ch_type])

            map_fname = op.join(sbj_path, subject + '_SM_mph_%s.stc' % ch_type)

            print('Saving morphed sensitivity map to %s.' % map_fname)
            map_morph.save(map_fname)

    return maps


# get all input arguments except first
if len(sys.argv) == 1:

    sbj_ids = np.arange(0, len(config.map_subjects)) + 1

else:

    # get list of subjects IDs to process
    sbj_ids = [int(aa) for aa in sys.argv[1:]]

for ss in sbj_ids:

    # dictionary for different channel types
    maps = run_make_sensitivity_maps(ss)

    subject = config.mri_subjects[ss]

    # channel types for which plotted sensitivity maps
    for ch_type in ch_types:

        time_label = '%s sensitivity' % ch_type

        print('4')

        fig = maps[ch_type].plot(time_label=time_label, subjects_dir=subjects_dir,
                                 clim=dict(kind='percent', lims=[0, 50, 100]))

        fig.show_view('lateral')

        fig_fname = op.join(bem_dir, subject + '_SM_%s.jpg' % ch_type)
        print('Saving figure to %s' % fig_fname)

        mlab.savefig(fig_fname)

print('Done.')
