#!/imaging/local/software/miniconda/envs/mne0.20/bin/python
"""
Compute noise covariance matrix for FPVS from rest runs.

Still all to do!
==========================================

OH, March 2020
"""

import sys

from os import path as op
import numpy as np

from importlib import reload

import mne

import config_sweep as config
reload(config)


print(mne.__version__)


def run_make_covmat(sbj_id):
    """Compute spectra for one subject."""
    # path to subject's data
    sbj_path = op.join(config.data_path, config.map_subjects[sbj_id][0])

    # raw-filename mappings for this subject
    sss_map_fname = config.sss_map_fnames[sbj_id]

    # Concatenate raws for two rest runs.
    # Assumes that rest files are the first two files in sss_map_fname.
    raw_all = []  # will contain list of all raw files
    for raw_stem_in in sss_map_fname[1][:2]:

        fname_raw = op.join(sbj_path, raw_stem_in[:-7] + 'sss_f_raw')

        # Read raw data info
        raw1 = mne.io.read_raw_fif(fname_raw + '.fif', preload=True)

        raw_all.append(raw1)

    # concatenate raws
    print('Concatenating %d raw files.' % len(raw_all))
    raw = mne.concatenate_raws(raw_all)

    cov = mne.cov.compute_raw_covariance(raw, reject=config.reject,
                                         method='auto', rank='info')

    fname_cov = op.join(sbj_path, 'rest_sss_f_raw-cov.fif')

    print('Writing covariance matrix to: %s.' % fname_cov)
    mne.cov.write_cov(fname=fname_cov, cov=cov)

    return


# get all input arguments except first
if len(sys.argv) == 1:

    sbj_ids = np.arange(0, len(config.map_subjects)) + 1

else:

    # get list of subjects IDs to process
    sbj_ids = [int(aa) for aa in sys.argv[1:]]


for ss in sbj_ids:

    # raw, psds, psds_as_evo, freqs = run_PSD_raw(ss)
    data_runs = run_make_covmat(ss)
