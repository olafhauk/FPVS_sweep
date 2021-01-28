#!/imaging/local/software/miniconda/envs/mne0.20/bin/python
"""
Apply ICA for FPVS Frequency Sweep.

Decompostion computed in FPVS_Compute_ICA.py
Based on Fiff_Apply_ICA.py.
==========================================
OH, modified October 2019
"""

import sys

from os import path as op
import numpy as np

from importlib import reload

import mne

print('MNE Version: %s\n\n' % mne.__version__)  # just in case
print(mne)

import config_sweep as config
reload(config)

# conditions
conds = config.do_conds

###############################################
### Parameters
###############################################


# "emulate" the args from ArgParser in Fiff_Apply_ICA.py
# filenames depend on subject, the rest are variables
class CreateArgs:
    """Parser for input arguments."""

    def __init__(self, FileRawIn, FileICA, FileRawOut):
        self.FileRawIn = FileRawIn
        self.FileICA = FileICA
        self.FileRawOut = FileRawOut


def run_Apply_ICA(sbj_id):
    """Apply previously computed ICA to raw data."""

    # path to subject's data
    sbj_path = op.join(config.data_path, config.map_subjects[sbj_id][0])

    # raw-filename mappings for this subject
    tmp_fnames = config.sss_map_fnames[sbj_id][1]

    # only use files for correct conditions
    sss_map_fnames = []
    for cond in conds:
        for [fi, ff] in enumerate(tmp_fnames):
            if cond in ff:
                sss_map_fnames.append(ff)

    for raw_stem_in in sss_map_fnames:

        # -ica.fif will be appended
        FileICA = op.join(sbj_path, raw_stem_in[:-7] + 'sss_f_raw')

        FileRawIn = op.join(sbj_path, raw_stem_in[:-7] + 'sss_f_raw')

        # -ica.fif will be appended
        # one file per subject
        FileICA = op.join(sbj_path, config.map_subjects[sbj_id][0] + '_sss_f_raw')

        # _ica_raw.fif will be appended
        FileRawOut = op.join(sbj_path, raw_stem_in[:-7] + 'sss_f_' +
                             config.raw_ICA_suff)

        # define variables for the following ICA pipeline
        # this would be from command line arguments of Fiff_Apply_ICA.py
        args = CreateArgs(FileRawIn, FileICA, FileRawOut)

        # Now turn the "fake command line parameters" into variables for the
        # analysis pipeline
        ###
        # ANALAYSIS PARAMETERS
        ###

        # raw-filenames to be subjected to ICA for this subject
        raw_fname_in = args.FileRawIn + '.fif'

        # save raw with ICA applied and artefacts removed
        if args.FileRawOut == '':

            raw_fname_out = args.FileRawIn + config.raw_ICA_suff + '.fif'

        else:

            raw_fname_out = args.FileRawOut + '.fif'

        # file with ICA decomposition
        if args.FileICA == '':

            ica_fname_in = args.FileRawIn + '-ica.fif'

        else:

            ica_fname_in = args.FileICA + '-ica.fif'

        ###
        # APPLY ICA
        ###

        print('Reading raw file %s' % raw_fname_in)
        raw = mne.io.read_raw_fif(raw_fname_in, preload=True)

        print('Reading ICA file %s' % ica_fname_in)
        ica = mne.preprocessing.read_ica(ica_fname_in)

        print('Applying ICA to raw file')
        ica.apply(raw)

        print('Saving raw file with ICA applied to %s' % raw_fname_out)
        raw.save(raw_fname_out, overwrite=True)

        # check if EEG in raw data
        if not raw.__contains__('eeg'):

            args.ChanTypes = ['meg']

            print('No EEG found in raw data, continuing with MEG only.')


# get all input arguments except first
if len(sys.argv) == 1:

    sbj_ids = np.arange(0, len(config.map_subjects)) + 1

else:

    # get list of subjects IDs to process
    sbj_ids = [int(aa) for aa in sys.argv[1:]]


for ss in sbj_ids:

    run_Apply_ICA(ss)
