#!/imaging/local/software/miniconda/envs/mne0.20/bin/python
"""
=========================================================
Plot sensors, scalp and source space to check alignment.
=========================================================

Currently (30.3.20) only plots correctly in 0.19
(in 0.20 only plots fiducials).
"""

print(__doc__)

import sys

from os import path as op
# sys.path.insert(1, '/imaging/local/software/mne_python/v0.11/')

import numpy as np

# import mayavi
from mayavi import mlab

mlab.options.offscreen = False

# import matplotlib
# matplotlib.use('Agg')  # possibly for running on cluster

# End

from importlib import reload

import mne

import config_sweep as config
reload(config)

print(mne.__version__)

subjects_dir = config.subjects_dir

# where figures will be written to
bem_dir = '/group/erp/data/olaf.hauk/MEG/FPVS/data/MRI/BEM_figs'

# to enable saving of mlab figures
# mlab.options.offscreen = True

def run_plot_coregistration(sbj_id):
    """Plot head and sensors to check coregistration."""
    subject = config.mri_subjects[sbj_id]

    if subject == '':

        print('No subject name for MRI specified - doing nothing now.')

        return

    print('Plotting coordinate alignment for subject %d, %s.' % (sbj_id, subject))

    # path to subject's EEG/MEG data
    sbj_path = op.join(config.data_path, config.map_subjects[sbj_id][0])

    # Source space
    src_fname = op.join(subjects_dir, subject, 'bem', subject + '_' + str(config.src_spacing) + '-src.fif')

    src = mne.read_source_spaces(src_fname)

    # coordinate transformation
    trans_fname = op.join(sbj_path, subject + '_cor-trans.fif')

    # raw file for info
    # doesn't matter which raw file, as long as transed
    raw_stem = config.sss_map_fnames[sbj_id][1][0]

    raw_fname = op.join(sbj_path, raw_stem[:-4] + '_f_' + config.raw_ICA_suff + '.fif')

    info = mne.io.read_info(raw_fname)

    # mayavi.engine.current_scene.scene.off_screen_rendering = True

    # fig_trans = mne.viz.plot_alignment(info, trans_fname, subjects_dir=subjects_dir, subject=subject,
    #                                    dig=True, meg='sensors', eeg='original', mri_fiducials=True, src=src)

    fig_bem = mne.viz.plot_bem(subjects_dir=subjects_dir, subject=subject, orientation='sagittal', brain_surfaces='white')

    mlab.show()

    # print('View:')
    # my_view = mlab.view()
    # print(my_view)

    # # change camera angle
    # # mlab.view(0., 90., 90.)

    # fig_fname = op.join(bem_dir, subject + '_alignment.jpg')

    # print('Saving figure %s.' % fig_fname)

    # mlab.savefig(fig_fname)

    # mlab.close(all=True)

    return fig_trans


# get all input arguments except first
if len(sys.argv) == 1:

    sbj_ids = np.arange(0, len(config.map_subjects)) + 1

else:

    # get list of subjects IDs to process
    sbj_ids = [int(aa) for aa in sys.argv[1:]]

for ss in sbj_ids:

    run_plot_coregistration(ss)

print('Done.')
