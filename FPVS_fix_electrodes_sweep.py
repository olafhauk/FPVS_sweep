#!/imaging/local/software/miniconda/envs/mne0.19/bin/python
"""
mne_check_eeg_locations for FPVS Frequency Sweep.

OH, modified October 2019
"""

import sys
import os
from os import path as op
import numpy as np

from importlib import reload

import config_sweep as config
reload(config)


def run_fix_electrodes(sbj_id):
    """Apply mne_check_eeg_locations to one subject."""
    # path to raw data for maxfilter
    map_subject = config.map_subjects[sbj_id]

    # raw-filename mappings for this subject
    sss_map_fname = config.sss_map_fnames[sbj_id]

    for raw_fname_out in sss_map_fname[1]:

        fname_sss = op.join(config.data_path, map_subject[0],
                            raw_fname_out + '.fif')

        print('Fixing electrode locations for %s.' % fname_sss)
        os.system(config.check_cmd % fname_sss)

# get all input arguments except first
if len(sys.argv) == 1:

    sbj_ids = np.arange(0, len(config.map_subjects)) + 1

else:

    # get list of subjects IDs to process
    sbj_ids = [int(aa) for aa in sys.argv[1:]]


for ss in sbj_ids:

    run_fix_electrodes(ss)
