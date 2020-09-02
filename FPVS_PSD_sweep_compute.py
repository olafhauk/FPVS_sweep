#!/imaging/local/software/miniconda/envs/mne0.20/bin/python
"""
Compute PSD for average raw sensor and source data for FPVS.

Reads average raw data from FPVS_get_sweeps.py.
Plot figures.
Compute z-scores.
Compute TFR if specified.
==========================================

OH, October 2019
removed plotting parts Feb 2020
added source space Apr 2020
"""

import sys

from os import path as op
import numpy as np

import os
# needed to run on SLURM
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

from copy import deepcopy

from mayavi import mlab
mlab.options.offscreen = True

import matplotlib
# for running graphics on cluster ### EDIT
# required even if not plotting?
matplotlib.use('Agg')

from importlib import reload

import mne

import config_sweep as config
reload(config)

import FPVS_functions as Ff
reload(Ff)


print(mne.__version__)

# perform TFR of raw data or not
# do_tfr = config.do_tfr

# separate filename prefixes for ICAed and non-ICAed data
prefix = ''
if 'ica' in config.raw_ICA_suff:
    prefix = 'ICA'

freqs_all = [str(ff) for ff in config.fpvs_freqs]

# conditions
# conds = ['face', 'pwhf', 'pwlf', 'lfhf']
conds = config.do_conds

def run_PSD_raw(sbj_id):
    """Compute spectra for one subject."""
    subject = config.mri_subjects[sbj_id]

    # path to subject's data
    sbj_path = op.join(config.data_path, config.map_subjects[sbj_id][0])

    inv_fname = op.join(sbj_path, subject + '_EEGMEG-inv.fif')

    print('Reading EEG/MEG inverse operator: %s.' % inv_fname)
    invop = mne.minimum_norm.read_inverse_operator(inv_fname)

    # raw-filename mappings for this subject
    sss_map_fname = config.sss_map_fnames[sbj_id]

    # # get condition names and frequency names
    # conds = []  # names of conditions
    # for raw_stem_in in sss_map_fname[1][2:]:

    #     conds.append(raw_stem_in[:4])

    # conds = np.unique(conds)

    print('Frequencies used: ')
    print(freqs_all)

    # initialise sum across harmonics for conditions
    sum_harms_odd = {}  # for oddball frequency
    sum_harms_base = {}  # for base frequencies
    # topographies for harmonics
    topos_harms_odd = {}  # for oddball frequency
    topos_harms_base = {}  # for base frequencies
    for cond in conds:

        sum_harms_odd[cond] = {}
        sum_harms_base[cond] = {}
        topos_harms_odd[cond] = {}
        topos_harms_base[cond] = {}

    # Go through conditions and frequencies
    # EDIT
    for cond in conds:  # conditions

        print('###\nCondition: %s.\n###' % cond)

        # create list of Evoked objects for all frequencies per condition
        psd_as_evo, psd_z_as_evo, sum_odd_as_evo, sum_base_as_evo,\
            topos_odd_as_evo, topos_base_as_evo, psd_harm_as_evo,\
            psd_harm_base_as_evo = [], [], [], [], [], [], [], []

        if cond == 'face':  # hack, no frequency sweep for faces

            freqs = ['6.0']

            # round to make sure combine_harmonics finds the right
            # frequencies
            oddfreq = round(config.fpvs_odd_freq['faces'], 2)

            # number of bins for z-scores
            snr_bins = config.psd_snr_bins['faces']

            nave = 5  # for inverse

        else:  # for all word condition, use all sweep frequencies

            freqs = freqs_all

            oddfreq = round(config.fpvs_odd_freq['words'], 2)

            # number of bins for z-scores
            snr_bins = config.psd_snr_bins['words']

            nave = 3  # for inverse

        for freq in freqs:  # frequencies

            basefreq = float(freq)  # hack, float-to-string-to-float-again

            # initialise for this base frequency
            sum_harms_odd[cond][freq] = []
            sum_harms_base[cond][freq] = []
            topos_harms_odd[cond][freq] = []
            topos_harms_base[cond][freq] = []

            # input average raw data; remove dot from frequency string
            fname = 'rawavg_%s_%s_%s.fif' % (cond, ''.join(freq.split('.')),
                                             config.raw_ICA_suff)

            fname_raw_in = op.join(sbj_path, fname)

            print('Reading average raw data from %s:' % fname_raw_in)

            raw = mne.io.read_raw_fif(fname_raw_in, preload=True)

            print('Resample to %s Hz.' % config.psd_resample)

            if config.psd_resample is not None:

                raw.resample(sfreq=config.psd_resample)

            # reduce raw data to relevant channels
            raw.pick_types(meg=True, eeg=True, eog=False, ecg=False,
                           stim=False, misc=False, chpi=False)

            # info = raw.info

            # Compute PSD for raw data

            # # find smallest power of 2 larger than number of samples
            # n_fft = 2**(len(raw.times) - 1).bit_length()

            n_fft = config.psd_nfft

            print('n_fft: %d' % n_fft)

            fmin = config.psd_fmin
            fmax = config.psd_fmax
            print('###\nComputing psd_welch() from %f to %f Hz.' %
                  (fmin, fmax))

            # print('Computing psd_welch() in sensor space.')
            # psds, psd_freqs = mne.time_frequency.psd_welch(raw, fmin=fmin,
            #                                                fmax=fmax,
            #                                                n_fft=n_fft)

            print('Computing psd_welch() in sensor and source space.')
            stc_psd, evo_psd = mne.minimum_norm.compute_source_psd(
                raw=raw, inverse_operator=invop, lambda2=1 / 9., method='MNE',
                fmin=fmin, fmax=fmax, n_fft=n_fft, overlap=.5,
                nave=nave, bandwidth='hann', low_bias=True, return_sensor=True)

            # turn power to amplitudes
            evo_psd.data = np.sqrt(evo_psd.data)

            stc_psd.data = np.sqrt(stc_psd.data)

            stc_psd.subject = subject

            fname_stc = op.join(sbj_path, 'STC', '%sPSDTopo_%s_%s' %
                                (prefix, cond, freq))

            stc_psd.save(fname_stc)

            # frequencies in PSD
            psd_freqs = evo_psd.times

            print('Frequencies from %f to %f.' % (psd_freqs[0], psd_freqs[-1]))

            freq_resol = psd_freqs[1] - psd_freqs[0]
            print('Frequency resolution:\n%f.\n###' % freq_resol)

            # # To plot PSDs like Evoked, pretend sample frequency
            # info['sfreq'] = 1. / (psd_freqs[1] - psd_freqs[0])

            # # Baseline-corrected PSDs as Evoked object
            # as_evo = mne.EvokedArray(psds, info,
            #                          tmin=psd_freqs[0],
            #                          comment=('PSDTopo_' + cond + '_' + freq))

            evo_psd.comment = 'PSDTopo_' + cond + '_' + freq

            psd_as_evo.append(evo_psd)

            # Z-score PSDs with neighbouring frequency bins

            print(type(evo_psd))

            print('Computing Z-scores for Evoked.')
            # inputing Evoked object
            psd_z = Ff.psd_z_score(evo_psd, snr_bins, mode='z',
                                   n_gap=config.psd_n_gap)

            # The following not needed anywhere?
            # print('Computing Z-scores for STC.')

            # psd_z_stc = Ff.psd_z_score(stc_psd, snr_bins, mode='z',
            #                            n_gap=config.psd_n_gap)

            # psd_z_stc.subject = subject

            # fname_stc = op.join(sbj_path, 'STC', '%sPSDTopoZ_%s_%s' %
            #                     (prefix, cond, freq))

            # psd_z_stc.save(fname_stc)

            # # z-scored PSDs as Evoked object
            # as_evo = mne.EvokedArray(psds_z, info,
            #                          tmin=psd_freqs[0],
            #                          comment=('PSDTopoZ_' + cond + '_' + freq))

            psd_z.comment = 'PSDTopoZ_' + cond + '_' + freq

            psd_z_as_evo.append(psd_z)

            # Compute the sum across harmonics of oddball frequency for this
            # condition and base frequency

            # TO DO: first sum amplitudes, then compute z-score
            # TO DO: summed topographies as z-scores should be taken from
            # summed epochs, frequency 0

            # # summed topography across harmonics of oddball frequency
            # print('Summing topographies for %d harmonics for oddball'
            #       ' frequency.' % config.fpvs_n_harms_odd)

            # # sum amplitudes
            # sum_harms_odd[cond][freq], topos_evoked, freqs_harm =\
            #     Ff.combine_harmonics_topos(
            #         psd=evo_psd, freqs=psd_freqs, basefreq=basefreq,
            #         oddfreq=oddfreq, n_harms=config.fpvs_n_harms_odd,
            #         method='sum')

            # topos_evoked.comment = ' '.join(str(freqs_harm))

            # topos_harms_odd[cond][freq] = topos_evoked

            # # sum amplitudes for STCs
            # sum_harms_odd_stc, topos_harms_stc, freqs_harm =\
            #     Ff.combine_harmonics_topos(
            #         psd=stc_psd, freqs=psd_freqs, basefreq=basefreq,
            #         oddfreq=oddfreq, n_harms=config.fpvs_n_harms_odd,
            #         method='sum')

            # sum_harms_odd_stc.subject = subject
            # # hack to keep frequencies of harmonics
            # topos_harms_stc.subject = ' '.join(str(freqs_harm))

            # fname_stc = op.join(sbj_path, 'STC', '%sPSDSumTopoOdd_%s_%s' %
            #                     (prefix, cond, freq))

            # sum_harms_odd_stc.save(fname_stc)

            # fname_stc = op.join(sbj_path, 'STC', '%sPSDSumToposOdd_%s_%s' %
            #                     (prefix, cond, freq))

            # topos_harms_stc.save(fname_stc)

            # # summed topography across harmonics of base frequency
            # print('Summing topographies for %d harmonics for base'
            #       ' frequency.' % config.fpvs_n_harms_base)

            # sum_harms_base[cond][freq], topos_evoked, freqs_harms =\
            #     Ff.combine_harmonics_topos(
            #         psd=evo_psd, freqs=psd_freqs, basefreq=None,
            #         oddfreq=basefreq, n_harms=config.fpvs_n_harms_base,
            #         method='sum')

            # # hack to keep frequencies of harmonics
            # topos_evoked.comment = ' '.join(str(freqs_harm))

            # topos_harms_base[cond][freq] = topos_evoked

            # sum_harms_base_stc, topos_harms_stc, freqs_harms =\
            #     Ff.combine_harmonics_topos(
            #         psd=stc_psd, freqs=psd_freqs, basefreq=None,
            #         oddfreq=basefreq, n_harms=config.fpvs_n_harms_base,
            #         method='sum')

            # sum_harms_base_stc.subject = subject
            # # hack to keep frequencies of harmonics
            # topos_harms_stc.subject = ' '.join(str(freqs_harms))

            # fname_stc = op.join(sbj_path, 'STC', '%sPSDSumTopoBase_%s_%s' %
            #                     (prefix, cond, freq))

            # sum_harms_base_stc.save(fname_stc)

            # fname_stc = op.join(sbj_path, 'STC', '%sPSDSumToposBase_%s_%s' %
            #                     (prefix, cond, freq))

            # topos_harms_stc.save(fname_stc)

            # # needs another dimension for Evoked object
            # to_evo_odd = np.expand_dims(sum_harms_odd[cond][freq], 1)
            # to_evo_base = np.expand_dims(sum_harms_base[cond][freq], 1)

            # # Combined harmonics as Evoked object
            # # Note: also for oddball frequency the "latency" is the base
            # # frequency, because that's our experimental manipulation
            # as_evo_odd = mne.EvokedArray(to_evo_odd, info,
            #                              tmin=basefreq,
            #                              comment=('PSDSumTopoOdd_' + cond +
            #                                       '_' + freq))

            # as_evo_base = mne.EvokedArray(to_evo_base, info,
            #                               tmin=basefreq,
            #                               comment=('PSDSumTopoBase_' + cond +
            #                                        '_' + freq))

            # sum_odd_as_evo.append(as_evo_odd)
            # sum_base_as_evo.append(as_evo_base)

            # sum_odd_as_evo.append(sum_harms_odd[cond][freq])
            # sum_base_as_evo.append(sum_harms_base[cond][freq])
            # topos_odd_as_evo.append(topos_harms_odd[cond][freq])
            # topos_base_as_evo.append(topos_harms_base[cond][freq])

            print('Summing PSDs across %d harmonics for oddball frequency' %
                  config.fpvs_n_harms_odd)

            # get PSDs around harmonics
            psd_harm, topo, topos, freqs_harm = Ff.psds_across_harmonics(
                psd=evo_psd, freqs=psd_freqs, basefreq=basefreq,
                oddfreq=oddfreq, n_harms=config.fpvs_n_harms_odd,
                n_bins=snr_bins, n_gap=config.psd_n_gap, method='sum')

            # get PSDs around harmonics for z-scores
            # needed to get z-scored topographies for harmonics
            psd_harm_z, topo_z, topos_z, freqs_harm_z = Ff.psds_across_harmonics(
                psd=psd_z, freqs=psd_freqs, basefreq=basefreq,
                oddfreq=oddfreq, n_harms=config.fpvs_n_harms_odd,
                n_bins=snr_bins, n_gap=config.psd_n_gap, method='sum')

            # compute z-score after summing
            psd_harm = Ff.psd_z_score(
                psd_harm, snr_bins, mode='z', n_gap=config.psd_n_gap)

            # Topography of z-scored summed harmonics at centre frequency
            topo_evo = deepcopy(psd_harm)
            topo_evo.crop(tmin=0., tmax=0.)
            sum_harms_odd[cond][freq] = topo_evo

            # z-scored topographies for individual harmonics
            topos.comment = ' '.join(str(freqs_harm_z))

            topos_harms_odd[cond][freq] = topos_z

            # STCs odd

            psd_harm_stc, topo, topos, freqs_harm = Ff.psds_across_harmonics(
                psd=stc_psd, freqs=psd_freqs, basefreq=basefreq,
                oddfreq=oddfreq, n_harms=config.fpvs_n_harms_odd,
                n_bins=snr_bins, n_gap=config.psd_n_gap, method='sum')

            # compute z-score after summing
            psd_harm_stc = Ff.psd_z_score(
                psd_harm_stc, snr_bins, mode='z', n_gap=config.psd_n_gap)

            psd_harm_stc.subject = subject

            fname_stc = op.join(sbj_path, 'STC', '%sPSDHarm_%s_%s' %
                                (prefix, cond, freq))

            psd_harm_stc.save(fname_stc)

            # MNE of z-scored summed harmonics at centre frequency
            topo_stc = deepcopy(psd_harm_stc)
            topo_stc.crop(tmin=0., tmax=0.)
            sum_harms_odd_stc = topo_stc

            sum_harms_odd_stc.subject = subject

            # topographies for individual harmonics
            topos_harms_stc = topos
            # hack to keep frequencies of harmonics
            topos_harms_stc.subject = ' '.join(str(freqs_harm))

            fname_stc = op.join(sbj_path, 'STC', '%sPSDSumTopoOdd_%s_%s' %
                                (prefix, cond, freq))

            sum_harms_odd_stc.save(fname_stc)

            fname_stc = op.join(sbj_path, 'STC', '%sPSDSumToposOdd_%s_%s' %
                                (prefix, cond, freq))

            topos_harms_stc.save(fname_stc)

            # BASE FREQUENCY

            psd_harm.comment = 'PSDHarm_' + cond + '_' + freq

            print('Summing PSDs across %d harmonics for base frequency' %
                  config.fpvs_n_harms_base)

            # Sanity check - do it for base frequency
            # i.e. basefreq as oddfreq here, for all its harmonics
            psd_harm_base, topo, topos, freqs_harm = Ff.psds_across_harmonics(
                psd=evo_psd, freqs=psd_freqs, basefreq=None, oddfreq=basefreq,
                n_harms=config.fpvs_n_harms_base, n_bins=snr_bins,
                n_gap=config.psd_n_gap, method='sum')

            # sum across harmonics for z-scores
            # needed to get z-scored topographies for harmonics
            psd_harm_base_z, topo_z, topos_z, freqs_harm_z = Ff.psds_across_harmonics(
                psd=psd_z, freqs=psd_freqs, basefreq=None, oddfreq=basefreq,
                n_harms=config.fpvs_n_harms_base, n_bins=snr_bins,
                n_gap=config.psd_n_gap, method='sum')

            # compute z-score after summing
            psd_harm_base = Ff.psd_z_score(
                psd_harm_base, snr_bins, mode='z', n_gap=config.psd_n_gap)

            # Topography of z-scored summed harmonics at centre frequency
            topo_evo = deepcopy(psd_harm_base)
            topo_evo.crop(tmin=0., tmax=0.)
            sum_harms_base[cond][freq] = topo_evo

            # z-scored topographies for individual harmonics
            topos.comment = ' '.join(str(freqs_harm))

            topos_harms_base[cond][freq] = topos_z

            # STCs base

            psd_harm_base_stc, topo, topos, freqs_harm = Ff.psds_across_harmonics(
                psd=stc_psd, freqs=psd_freqs, basefreq=None,
                oddfreq=basefreq, n_harms=config.fpvs_n_harms_base,
                n_bins=snr_bins, n_gap=config.psd_n_gap, method='sum')

            # compute z-score after summing
            psd_harm_base_stc = Ff.psd_z_score(
                psd_harm_base_stc, snr_bins, mode='z', n_gap=config.psd_n_gap)

            psd_harm_base_stc.subject = subject

            fname_stc = op.join(sbj_path, 'STC', '%sPSDHarmBase_%s_%s' %
                                (prefix, cond, freq))

            psd_harm_base_stc.save(fname_stc)

            psd_harm_base.comment = 'PSDHarmBase_' + cond + '_' + freq

            # MNE of z-scored summed harmonics at centre frequency
            topo_stc = deepcopy(psd_harm_base_stc)
            topo_stc.crop(tmin=0., tmax=0.)
            sum_harms_base_stc = topo_stc

            sum_harms_base_stc.subject = subject

            # topographies for individual harmonics
            topos_harms_stc = topos
            # hack to keep frequencies of harmonics
            topos_harms_stc.subject = ' '.join(str(freqs_harm))

            fname_stc = op.join(sbj_path, 'STC', '%sPSDSumTopoBase_%s_%s' %
                                (prefix, cond, freq))

            sum_harms_base_stc.save(fname_stc)

            fname_stc = op.join(sbj_path, 'STC', '%sPSDSumToposBase_%s_%s' %
                                (prefix, cond, freq))

            topos_harms_stc.save(fname_stc)

            # tmin = -(snr_bins + 1) * freq_resol  # include baseline
            # info['sfreq'] = 1. / freq_resol  # to display samples as time pts

            # as_evo = mne.EvokedArray(psd_harm, info, tmin=tmin,
            #                          comment=('PSDHarm_' + cond + '_' + freq))
            # psd_harm_as_evo.append(as_evo)

            # as_evo = mne.EvokedArray(psd_harm_base, info, tmin=tmin,
            #                          comment=('PSDHarmBase_' + cond + '_' +
            #                                   freq))
            # psd_harm_base_as_evo.append(as_evo)

            psd_harm_as_evo.append(psd_harm)

            psd_harm_base_as_evo.append(psd_harm_base)

            sum_odd_as_evo.append(sum_harms_odd[cond][freq])
            sum_base_as_evo.append(sum_harms_base[cond][freq])
            topos_odd_as_evo.append(topos_harms_odd[cond][freq])
            topos_base_as_evo.append(topos_harms_base[cond][freq])


        # Save Evoked objects for later group stats:

        print('Saving PSD results as evoked files:')

        fname_evo = op.join(sbj_path, 'AVE', '%sPSDTopo_%s%s' %
                            (prefix, cond, '-ave.fif'))
        print(fname_evo)
        mne.write_evokeds(fname_evo, psd_as_evo)

        fname_evo = op.join(sbj_path, 'AVE', '%sPSDTopoZ_%s%s' %
                            (prefix, cond, '-ave.fif'))
        print(fname_evo)
        mne.write_evokeds(fname_evo, psd_z_as_evo)

        fname_evo = op.join(sbj_path, 'AVE', '%sPSDHarm_%s%s' %
                            (prefix, cond, '-ave.fif'))
        print(fname_evo)
        mne.write_evokeds(fname_evo, psd_harm_as_evo)

        fname_evo = op.join(sbj_path, 'AVE', '%sPSDHarmBase_%s%s' %
                            (prefix, cond, '-ave.fif'))
        print(fname_evo)
        mne.write_evokeds(fname_evo, psd_harm_base_as_evo)

        fname_evo = op.join(sbj_path, 'AVE', '%sPSDSumTopoOdd_%s%s' %
                            (prefix, cond, '-ave.fif'))
        print(fname_evo)
        mne.write_evokeds(fname_evo, sum_odd_as_evo)

        fname_evo = op.join(sbj_path, 'AVE', '%sPSDSumTopoBase_%s%s' %
                            (prefix, cond, '-ave.fif'))
        print(fname_evo)
        mne.write_evokeds(fname_evo, sum_base_as_evo)

        fname_evo = op.join(sbj_path, 'AVE', '%sPSDSumToposOdd_%s%s' %
                            (prefix, cond, '-ave.fif'))
        print(fname_evo)
        mne.write_evokeds(fname_evo, topos_odd_as_evo)

        fname_evo = op.join(sbj_path, 'AVE', '%sPSDSumToposBase_%s%s' %
                            (prefix, cond, '-ave.fif'))
        print(fname_evo)
        mne.write_evokeds(fname_evo, topos_base_as_evo)

    return


# get all input arguments except first
if len(sys.argv) == 1:

    sbj_ids = np.arange(0, len(config.map_subjects)) + 1

else:

    # get list of subjects IDs to process
    sbj_ids = [int(aa) for aa in sys.argv[1:]]


for ss in sbj_ids:

    # raw, psds, psds_as_evo, freqs = run_PSD_raw(ss)
    run_PSD_raw(ss)
