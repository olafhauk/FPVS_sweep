#!/imaging/local/software/miniconda/envs/mne0.20/bin/python
"""
=========================================================
BEM Model and Source Space for FPVS.

Needs Freesurfer. Run freesurfer_6.0.0 before starting
mne_python_v0.19.
=========================================================

"""
# OH, March 2020

print(__doc__)

import sys

import os
from os import path as op

import numpy as np

import importlib
from importlib import reload

import mne

import config_sweep as config
reload(config)

# Whether to create or update symbolic links or not
# Note: This would overwrite any links manually created,
# e.g. after shrinking inner skull
do_links = False

# set Freesurfer environment variables
filename = "/imaging/local/software/mne_python/set_MNE_2.7.3_FS_6.0.0_environ.py"
# for Python 3 instead of execfile
exec(compile(open(filename, "rb").read(), filename, 'exec'))

# where MRIs are
subjects_dir = config.subjects_dir

# where figures will be written to
bem_fig_dir = '/group/erp/data/olaf.hauk/MEG/FPVS/data/MRI/BEM_figs'

# default conductivities
conductivity_1 = (0.3,)  # for single layer
conductivity_3 = (0.3, 0.006, 0.3)  # for three layers


def run_make_bem(sbj_id):

    subject = config.mri_subjects[sbj_id]

    if subject == '':

        print('No subject name for MRI specified - doing nothing now.')

        return

    if do_links:

        print('Updating links to surfaces.')
        cmd_list = []
        sd = config.subjects_dir
        sb = subject

        cmd_list.append('ln -sf %s/%s/bem/watershed/%s_inner_skull_surface \
                        %s/%s/bem/inner_skull.surf' % (sd, sb, sb, sd, sb))
        cmd_list.append('ln -sf %s/%s/bem/watershed/%s_outer_skull_surface \
                        %s/%s/bem/outer_skull.surf' % (sd, sb, sb, sd, sb))
        cmd_list.append('ln -sf %s/%s/bem/watershed/%s_outer_skin_surface \
                        %s/%s/bem/outer_skin.surf' % (sd, sb, sb, sd, sb))
        cmd_list.append('ln -sf %s/%s/bem/watershed/%s_brain_surface \
                        %s/%s/bem/brain_surface.surf' % (sd, sb, sb, sd, sb))

        [os.system(cmd) for cmd in cmd_list]

    else:

        print('###\nNOT updating symbolic links to surfaces.\n###')

    print('Making BEM model for %s.' % subject)

    print('Creating BEM surfaces.')

    # one-shell BEM for MEG
    print('###\nMEG\n###')
    model = mne.make_bem_model(subject=subject, ico=4,
                               conductivity=conductivity_1,
                               subjects_dir=subjects_dir)

    bem = mne.make_bem_solution(model)

    bem_fname = op.join(subjects_dir, subject, 'bem', subject + '_MEG-bem.fif')

    print('Writing BEM solution for MEG to file %s.' % bem_fname)
    mne.bem.write_bem_solution(bem_fname, bem)

    # three-shell BEM for EEG+MEG
    print('EEG+MEG')
    model = mne.make_bem_model(subject=subject, ico=4,
                               conductivity=conductivity_3,
                               subjects_dir=subjects_dir)

    bem = mne.make_bem_solution(model)

    bem_fname = op.join(subjects_dir, subject, 'bem', subject + '_EEGMEG-bem.fif')

    print('Writing BEM solution for EEG+MEG to file %s.' % bem_fname)

    mne.bem.write_bem_solution(bem_fname, bem)

# get all input arguments except first
if len(sys.argv) == 1:

    sbj_ids = np.arange(0, len(config.map_subjects)) + 1

else:

    # get list of subjects IDs to process
    sbj_ids = [int(aa) for aa in sys.argv[1:]]

for ss in sbj_ids:

    run_make_bem(ss)

print('Done.')
