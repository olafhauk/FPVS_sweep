#!/imaging/local/software/miniconda/envs/mne0.20/bin/python
"""
=========================================================
Make inverse operator for FPVS.
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


def run_make_inverse_operator(sbj_id):

    subject = config.mri_subjects[sbj_id]

    if subject == '':

        print('No subject name for MRI specified - doing nothing now.')

        return

    print('Making Inverse Operator for %s.' % subject)

    sbj_path = op.join(config.data_path, config.map_subjects[sbj_id][0])

    # doesn't matter which raw file, as long as transed
    raw_stem = config.sss_map_fnames[sbj_id][1][0]
    raw_fname = op.join(sbj_path, raw_stem[:-4] + '_f_' +
                        config.raw_ICA_suff + '.fif')

    info = mne.io.read_info(raw_fname)

    fwd_fname = op.join(sbj_path, subject + '_MEG-fwd.fif')
    print('Reading MEG forward solution: %s.' % fwd_fname)

    fwd_meg = mne.read_forward_solution(fwd_fname)

    fwd_fname = op.join(sbj_path, subject + '_EEGMEG-fwd.fif')
    print('Reading EEG/MEG forward solution: %s.' % fwd_fname)

    fwd_eegmeg = mne.read_forward_solution(fwd_fname)

    fname_cov = op.join(sbj_path, 'rest_sss_f_raw-cov.fif')

    print('Reading covariance matrix: %s.' % fname_cov)
    noise_cov = mne.cov.read_cov(fname=fname_cov)

    # make inverse operator
    loose = 0.2
    depth = None

    invop_meg = mne.minimum_norm.make_inverse_operator(info, fwd_meg, noise_cov,
                                                       fixed='auto', loose=loose, depth=depth,
                                                       rank='info')

    invop_eegmeg = mne.minimum_norm.make_inverse_operator(info, fwd_eegmeg, noise_cov,
                                                          fixed='auto', loose=loose, depth=depth,
                                                          rank='info')

    inv_fname = op.join(sbj_path, subject + '_MEG-inv.fif')
    print('Writing MEG inverse operator: %s.' % inv_fname)
    mne.minimum_norm.write_inverse_operator(fname=inv_fname, inv=invop_meg)

    inv_fname = op.join(sbj_path, subject + '_EEGMEG-inv.fif')
    print('Writing EEG/MEG inverse operator: %s.' % inv_fname)
    mne.minimum_norm.write_inverse_operator(fname=inv_fname, inv=invop_eegmeg)

# get all input arguments except first
if len(sys.argv) == 1:

    sbj_ids = np.arange(0, len(config.map_subjects)) + 1

else:

    # get list of subjects IDs to process
    sbj_ids = [int(aa) for aa in sys.argv[1:]]

for ss in sbj_ids:

    run_make_inverse_operator(ss)

print('Done.')
