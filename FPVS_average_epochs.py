#!/imaging/local/software/miniconda/envs/mne0.20/bin/python
"""
Average epochs from FPVS sweeps for ERP analysis.

Read epochs created by FPVS_epochs_sweeps.py and average.
==========================================

OH, May 2020
"""

import sys

from os import path as op
import numpy as np

from copy import deepcopy

from importlib import reload

import mne

import config_sweep as config
reload(config)


print(mne.__version__)

# conditions
conds = config.do_conds

def run_average_epochs(sbj_id):
    """Average epochs for one subject."""
    # path to subject's data
    sbj_path = op.join(config.data_path, config.map_subjects[sbj_id][0])

    # base frequencies as strings
    freqs_all = [str(ff) for ff in config.fpvs_freqs]

    # for evoked created with and without Notch filter for base frequency
    for do_notch in [0, 1]:

        if do_notch:  # if Notch filter at base frequency requested

            # add to epoch file name
            str_notch = '_nch'

        else:

            str_notch = ''

        for cond in conds:  # conditions

            if cond == 'face':  # hack, no frequency sweep for faces

                freqs = ['6.0']

            else:  # for all word condition, use all sweep frequencies

                freqs = freqs_all

            for freq in freqs:  # frequencies

                epo_fname = op.join(
                    sbj_path, 'EPO', '%s_f_%s_%s%s-epo.fif' %
                    (cond, config.raw_ICA_suff, ''.join(freq.split('.')),
                     str_notch))

                print('Reading epochs from %s.' % epo_fname)

                epochs = mne.read_epochs(epo_fname)

                # averaging epochs
                evoked = epochs.average()

                # projection necessary for source estimation
                evoked.set_eeg_reference(ref_channels='average',
                                         projection=True)

                evo_fname = op.join(
                    sbj_path, 'AVE', '%s_f_%s_%s%s-ave.fif' %
                    (cond, config.raw_ICA_suff, ''.join(freq.split('.')),
                     str_notch))

                print('Writing evoked data to %s.' % evo_fname)
                mne.write_evokeds(evo_fname, evoked)

    return


# get all input arguments except first
if len(sys.argv) == 1:

    sbj_ids = np.arange(0, len(config.map_subjects)) + 1

else:

    # get list of subjects IDs to process
    sbj_ids = [int(aa) for aa in sys.argv[1:]]


for ss in sbj_ids:

    # raw, psds, psds_as_evo, freqs = run_PSD_raw(ss)
    data_runs = run_average_epochs(ss)
