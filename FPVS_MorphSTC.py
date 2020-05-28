#!/imaging/local/software/miniconda/envs/mne0.20/bin/python
"""
Morph individual STCs for FPVS.

==========================================

OH, April 2020
"""

import sys

from os import path as op

import numpy as np

from importlib import reload

import mne

import FPVS_functions
reload(FPVS_functions)

import config_sweep as config
reload(config)

print(mne.__version__)

# for some plots of SNRs
unit_scalings = dict(eeg=1., mag=1., grad=1.)

# Base frequencies for frequency sweep for words (not faces)
freqs_all = [str(ff) for ff in config.fpvs_freqs]

# separate filename prefixes for ICAed and non-ICAed data
prefix = ''
if 'ica' in config.raw_ICA_suff:
    prefix = 'ICA'

# which data to morph
types = ['psd', 'psd_z', 'psd_sum_odd', 'psd_sum_base', 'psd_harm_odd',
         'psd_harm_base', 'psd_harm_topos_odd', 'psd_harm_topos_base']

subjects_dir = config.subjects_dir

def morph_stcs(sbj_ids):
    """Morph STCs for sbj_ids."""
    print('Morphing STCs for subjects:')
    print(*sbj_ids)

    # get condition names and frequency names
    # assumed to be consistent across subjects
    sss_map_fname = config.sss_map_fnames[sbj_ids[0]]
    conds = []  # names of conditions
    for raw_stem_in in sss_map_fname[1][2:]:

        conds.append(raw_stem_in[:4])

    conds = np.unique(conds)

    # for Evoked data are in one file for all frequencies
    # for STC data are in separate files per condition and freq
    for sbj_id in sbj_ids:  # across all subjects, EDIT ###

        # path to subject's data
        sbj_path = op.join(config.data_path,
                           config.map_subjects[sbj_id][0])

        subject = config.mri_subjects[sbj_id]

        if subject == '':

            print('No subject name for MRI specified - doing nothing now.')

            return

        # getting an STC for source space for this subject
        fname_stc = op.join(
            sbj_path, 'STC', '%sPSDTopo_%s_%s-lh.stc' % (prefix, 'face',
                                                         '6.0')
        )
        stc_from = mne.read_source_estimate(fname_stc)

        print('Computing morphing matrix for %s.' % subject)

        morph_mat = mne.compute_source_morph(
            src=stc_from, subject_from=subject, subject_to=config.stc_morph,
            subjects_dir=subjects_dir
        )

        for cond in conds:  # conditions

            print('Condition: %s.' % cond)

            if cond == 'face':  # no frequency sweep for faces

                freqs = ['6.0']  # base frequency for this condition (Hz, str)

            else:  # for all word condition, use all sweep frequencies

                # base frequencies for this condition (Hz as string)
                freqs = freqs_all

            for (fi, freq) in enumerate(freqs):

                    print('Reading PSD results from STC files:')

                    fname_stc = op.join(
                        sbj_path, 'STC', '%sPSDTopo_%s_%s-lh.stc' %
                        (prefix, cond, freq)
                    )
                    print(fname_stc)
                    stc = mne.read_source_estimate(fname_stc)
                    stc_mph = morph_mat.apply(stc)
                    fname_mph = op.join(
                        sbj_path, 'STC', '%sPSDTopo_%s_%s_mph-lh.stc' %
                        (prefix, cond, freq)
                    )
                    stc_mph.save(fname_mph)
                    #
                    fname_stc = op.join(
                        sbj_path, 'STC', '%sPSDTopoZ_%s_%s-lh.stc' %
                        (prefix, cond, freq)
                    )
                    print(fname_stc)
                    stc = mne.read_source_estimate(fname_stc)
                    stc_mph = morph_mat.apply(stc)
                    fname_mph = op.join(
                        sbj_path, 'STC', '%sPSDTopoZ_%s_%s_mph-lh.stc' %
                        (prefix, cond, freq)
                    )
                    stc_mph.save(fname_mph)
                    #
                    fname_stc = op.join(
                        sbj_path, 'STC', '%sPSDHarm_%s_%s-lh.stc' %
                        (prefix, cond, freq)
                    )
                    print(fname_stc)
                    stc = mne.read_source_estimate(fname_stc)
                    stc_mph = morph_mat.apply(stc)
                    fname_mph = op.join(
                        sbj_path, 'STC', '%sPSDHarm_%s_%s_mph-lh.stc' %
                        (prefix, cond, freq)
                    )
                    stc_mph.save(fname_mph)
                    #
                    fname_stc = op.join(
                        sbj_path, 'STC', '%sPSDHarmBase_%s_%s-lh.stc' %
                        (prefix, cond, freq)
                    )
                    print(fname_stc)
                    stc = mne.read_source_estimate(fname_stc)
                    stc_mph = morph_mat.apply(stc)
                    fname_mph = op.join(
                        sbj_path, 'STC', '%sPSDHarmBase_%s_%s_mph-lh.stc' %
                        (prefix, cond, freq)
                    )
                    stc_mph.save(fname_mph)
                    #
                    fname_stc = op.join(
                        sbj_path, 'STC', '%sPSDSumTopoOdd_%s_%s-lh.stc' %
                        (prefix, cond, freq)
                    )
                    print(fname_stc)
                    stc = mne.read_source_estimate(fname_stc)
                    stc_mph = morph_mat.apply(stc)
                    fname_mph = op.join(
                        sbj_path, 'STC', '%sPSDSumTopoOdd_%s_%s_mph-lh.stc' %
                        (prefix, cond, freq)
                    )
                    stc_mph.save(fname_mph)
                    #
                    fname_stc = op.join(
                        sbj_path, 'STC', '%sPSDSumTopoBase_%s_%s-lh.stc' %
                        (prefix, cond, freq)
                    )
                    print(fname_stc)
                    stc = mne.read_source_estimate(fname_stc)
                    stc_mph = morph_mat.apply(stc)
                    fname_mph = op.join(
                        sbj_path, 'STC', '%sPSDSumTopoBase_%s_%s_mph-lh.stc' %
                        (prefix, cond, freq)
                    )
                    stc_mph.save(fname_mph)
                    #
                    fname_stc = op.join(
                        sbj_path, 'STC', '%sPSDSumToposOdd_%s_%s-lh.stc' %
                        (prefix, cond, freq)
                    )
                    print(fname_stc)
                    stc = mne.read_source_estimate(fname_stc)
                    stc_mph = morph_mat.apply(stc)
                    fname_mph = op.join(
                        sbj_path, 'STC', '%sPSDSumToposOdd_%s_%s_mph-lh.stc' %
                        (prefix, cond, freq)
                    )
                    stc_mph.save(fname_mph)
                    #
                    fname_stc = op.join(
                        sbj_path, 'STC', '%sPSDSumToposBase_%s_%s-lh.stc' %
                        (prefix, cond, freq)
                    )
                    print(fname_stc)
                    stc = mne.read_source_estimate(fname_stc)
                    stc_mph = morph_mat.apply(stc)
                    fname_mph = op.join(
                        sbj_path, 'STC', '%sPSDSumToposBase_%s_%s_mph-lh.stc' %
                        (prefix, cond, freq)
                    )
                    stc_mph.save(fname_mph)

    return


# get all input arguments except first
if ((len(sys.argv) == 1) or
   (int(sys.argv[1]) > np.max(list(config.map_subjects.keys())))):

    # IDs don't start at 0
    sbj_ids = config.do_subjs

else:

    # get list of subjects IDs to process
    sbj_ids = [int(aa) for aa in sys.argv[1:]]


# requires all subjects to average across
morph_stcs(sbj_ids)
