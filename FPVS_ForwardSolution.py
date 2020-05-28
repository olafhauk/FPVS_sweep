#!/imaging/local/software/miniconda/envs/mne0.20/bin/python
"""
=========================================================
Make forward solution for FPVS.
=========================================================

"""

import sys
from os import path as op

import numpy as np

from importlib import reload

import mne

import config_sweep as config
reload(config)

subjects_dir = config.subjects_dir


def run_make_forward_solution(sbj_id):

    subject = config.mri_subjects[sbj_id]

    if subject == '':

        print('No subject name for MRI specified - doing nothing now.')

        return

    print('Making Forward Solution for %s.' % subject)

    sbj_path = op.join(config.data_path, config.map_subjects[sbj_id][0])

    # doesn't matter which raw file, as long as transed
    raw_stem = config.sss_map_fnames[sbj_id][1][0]
    raw_fname = op.join(sbj_path, raw_stem[:-4] + '_f_' +
                        config.raw_ICA_suff + '.fif')

    src_fname = op.join(subjects_dir, subject, 'bem', subject + '_' + str(config.src_spacing) + '-src.fif')

    print('Source space from: %s' % src_fname)
    src = mne.read_source_spaces(src_fname)

    # coordinate transformation
    trans_fname = op.join(sbj_path, subject + '_cor-trans.fif')

    # one-shell BEM for MEG
    bem_fname = op.join(subjects_dir, subject, 'bem', subject + '_MEG-bem.fif')

    print('BEM: %s' % bem_fname)
    bem = mne.bem.read_bem_solution(bem_fname)

    fwd_fname = op.join(sbj_path, subject + '_MEG-fwd.fif')
    print('Making forward solution: %s.' % fwd_fname)

    fwd_meg = mne.make_forward_solution(raw_fname, trans=trans_fname, src=src, bem=bem,
                                        meg=True, eeg=False, mindist=5.0, verbose=True)

    mne.write_forward_solution(fname=fwd_fname, fwd=fwd_meg, overwrite=True)

    ### three-shell BEM for MEG
    bem_fname = op.join(subjects_dir, subject, 'bem', subject + '_EEGMEG-bem.fif')
    print('BEM: %s' % bem_fname)
    bem = mne.bem.read_bem_solution(bem_fname)

    fwd_fname = op.join(sbj_path, subject + '_EEGMEG-fwd.fif')
    print('Making forward solution: %s.' % fwd_fname)

    fwd_eegmeg = mne.make_forward_solution(info=raw_fname, trans=trans_fname, src=src, bem=bem,
                                           meg=True, eeg=True, mindist=5.0, verbose=True)

    mne.write_forward_solution(fname=fwd_fname, fwd=fwd_eegmeg, overwrite=True)

# get all input arguments except first
if len(sys.argv) == 1:

    sbj_ids = np.arange(0, len(config.map_subjects)) + 1

else:

    # get list of subjects IDs to process
    sbj_ids = [int(aa) for aa in sys.argv[1:]]

for ss in sbj_ids:

    run_make_forward_solution(ss)

print('Done.')
