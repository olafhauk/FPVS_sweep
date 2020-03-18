#!/imaging/local/software/miniconda/envs/mne0.19/bin/python
"""
Compute PSD for average raw data for FPVS Frequency Sweep.

Average raw data from FPVS_get_sweeps.py.
Plot figures.
Compute z-scores.
Compute TFR if specified.
==========================================

OH, October 2019
removed plotting parts Feb 2020
"""

import sys

from os import path as op
import numpy as np

# import matplotlib
# matplotlib.use('Agg') #  for running graphics on cluster ### EDIT

from importlib import reload

import mne

import config_sweep as config
reload(config)

import FPVS_functions as Ff
reload(Ff)


print(mne.__version__)

# perform TFR of raw data or not
# do_tfr = config.do_tfr

print('Sunshine.')


def run_PSD_raw(sbj_id):
    """Compute spectra for one subject."""
    # path to subject's data
    sbj_path = op.join(config.data_path, config.map_subjects[sbj_id][0])

    # raw-filename mappings for this subject
    sss_map_fname = config.sss_map_fnames[sbj_id]

    # get condition names and frequency names
    conds = []  # names of conditions
    for raw_stem_in in sss_map_fname[1][2:]:

        conds.append(raw_stem_in[:4])

    conds = np.unique(conds)

    freqs_all = [str(ff) for ff in config.fpvs_freqs]

    print('Frequencies used: ')
    print(freqs_all)

    # initialise sum across harmonics for conditions
    sum_harms_odd = {}  # for oddball frequency
    sum_harms_base = {}  # for base frequencies
    for cond in conds:

        sum_harms_odd[cond] = {}
        sum_harms_base[cond] = {}

    # Go through conditions and frequencies
    for cond in conds:  # conditions

        print('###\nCondition: %s.\n###' % cond)

        # create list of Evoked objects for all frequencies per condition
        psds_as_evo, psds_z_as_evo, sum_odd_as_evo, sum_base_as_evo,\
            psd_harm_as_evo, psd_harm_base_as_evo = [], [], [], [], [], []

        if cond == 'face':  # hack, no frequency sweep for faces

            freqs = ['6.0']

            snr_bins = config.psd_snr_bins['faces']  # number of bins for z-scores

        else:  # for all word condition, use all sweep frequencies

            freqs = freqs_all

            snr_bins = config.psd_snr_bins['words']  # number of bins for z-scores

        for freq in freqs:  # frequencies

            sum_harms_odd[cond][freq] = []  # initialise for this base frequency
            sum_harms_base[cond][freq] = []  # initialise for this base frequency

            # input average raw data; remove dot from frequency string
            fname = 'rawavg_%s_%s_%s.fif' % (cond, ''.join(freq.split('.')),
                                             config.raw_ICA_suff)

            fname_raw_in = op.join(sbj_path, fname)

            print('Reading average raw data from %s:' % fname_raw_in)

            raw = mne.io.read_raw_fif(fname_raw_in, preload=True)

            # reduce raw data to relevant channels
            raw.pick_types(meg=True, eeg=True, eog=False, ecg=False,
                           stim=False, misc=False, chpi=False)

            # Compute PSD for raw data
            n_fft = len(raw.times)

            fmin = config.psd_fmin
            fmax = config.psd_fmax
            print('###\nComputing psd_welch() from %f to %f Hz.' % (fmin, fmax))

            psds, psd_freqs = mne.time_frequency.psd_welch(raw, fmin=fmin,
                                                           fmax=fmax,
                                                           n_fft=n_fft)

            freq_resol = psd_freqs[1] - psd_freqs[0]
            print('Frequency resolution:\n%f.\n###' % freq_resol)

            info = raw.info

            # To plot PSDs like Evoked, pretend sample frequency
            info['sfreq'] = 1. / (psd_freqs[1] - psd_freqs[0])

            # convert PSD to amplitudes (rather than power)
            psds = np.sqrt(psds)

            # Z-score PSDs with neighbouring frequency bins

            print('Computing Z-scores.')
            psds_z = Ff.psd_z_score(psds, snr_bins, mode='z',
                                    n_gap=config.psd_n_gap)

            # Baseline-corrected PSDs as Evoked object
            as_evo = mne.EvokedArray(psds, info,
                                     tmin=psd_freqs[0],
                                     comment=('PSDTopo_' + cond + '_' + freq))
            psds_as_evo.append(as_evo)

            # z-scored PSDs as Evoked object
            as_evo = mne.EvokedArray(psds_z, info,
                                     tmin=psd_freqs[0],
                                     comment=('PSDTopoZ_' + cond + '_' + freq))
            psds_z_as_evo.append(as_evo)

            # Plot PSD as spectrum plus topographies (plot_joint())
            print('Plotting PSDs.')

            # Compute the sum across harmonics of oddball frequency for this
            # condition and base frequency
            if cond == 'face':

                # round to make sure combine_harmonics finds the right
                # frequencies
                oddfreq = round(config.fpvs_odd_freq['faces'], 2)

            else:

                oddfreq = round(config.fpvs_odd_freq['words'], 2)

            basefreq = float(freq)  # hack, float-to-string-to-float-again

            # summed topography across harmonics of oddball frequency
            print('Summing topographies for %d harmonics for oddball'
                  ' frequency.' % config.fpvs_n_harms_odd)

            sum_harms_odd[cond][freq] =\
                Ff.combine_harmonics_topos(psds=psds_z,
                                           freqs=psd_freqs,
                                           basefreq=basefreq,
                                           oddfreq=oddfreq,
                                           n_harms=config.fpvs_n_harms_odd,
                                           method='sum')

            # summed topography across harmonics of base frequency
            print('Summing topographies for %d harmonics for base'
                  ' frequency.' % config.fpvs_n_harms_base)

            sum_harms_base[cond][freq] =\
                Ff.combine_harmonics_topos(psds=psds_z,
                                           freqs=psd_freqs,
                                           basefreq=999.,
                                           oddfreq=basefreq,
                                           n_harms=config.fpvs_n_harms_base,
                                           method='sum')

            # needs another dimension for Evoked object
            to_evo_odd = np.expand_dims(sum_harms_odd[cond][freq], 1)
            to_evo_base = np.expand_dims(sum_harms_base[cond][freq], 1)

            # Combined harmonics as Evoked object
            # Note: also for oddball frequency the "latency" is the base
            # frequency, because that's our experimental manipulation
            as_evo_odd = mne.EvokedArray(to_evo_odd, info,
                                         tmin=basefreq,
                                         comment=('PSDSumTopoOdd_' + cond + '_' + freq))

            as_evo_base = mne.EvokedArray(to_evo_base, info,
                                          tmin=basefreq,
                                          comment=('PSDSumTopoBase_' + cond + '_' + freq))

            sum_odd_as_evo.append(as_evo_odd)
            sum_base_as_evo.append(as_evo_base)

            print('Summing PSDs across %d harmonics for oddball frequency' %
                  config.fpvs_n_harms_odd)

            # get PSDs around harmonics
            psd_harm = Ff.psds_across_harmonics(psds=psds_z, freqs=psd_freqs,
                                                basefreq=basefreq,
                                                oddfreq=oddfreq,
                                                n_harms=config.fpvs_n_harms_odd,
                                                n_bins=snr_bins,
                                                n_gap=config.psd_n_gap,
                                                method='sum')

            print('Summing PSDs across %d harmonics for base frequency' %
                  config.fpvs_n_harms_base)

            # Sanity check - do it for base frequency
            # i.e. basefreq as oddfreq here, for all its harmonics
            psd_harm_base = Ff.psds_across_harmonics(psds=psds_z, freqs=psd_freqs,
                                                     basefreq=999.,
                                                     oddfreq=basefreq,
                                                     n_harms=config.fpvs_n_harms_base,
                                                     n_bins=snr_bins,
                                                     n_gap=config.psd_n_gap,
                                                     method='sum')

            tmin = -(snr_bins + 1) * freq_resol  # include baseline
            info['sfreq'] = 1. / freq_resol  # to display samples as time points

            as_evo = mne.EvokedArray(psd_harm, info, tmin=tmin,
                                     comment=('PSDHarm_' + cond + '_' + freq))
            psd_harm_as_evo.append(as_evo)

            as_evo = mne.EvokedArray(psd_harm_base, info, tmin=tmin,
                                     comment=('PSDHarmBase_' + cond + '_' + freq))
            psd_harm_base_as_evo.append(as_evo)

        # Save Evoked objects for later group stats:

        print('Saving PSD results as evoked files:')

        # separate filename prefixes for ICAed and non-ICAed data
        prefix = ''
        if 'ica' in config.raw_ICA_suff:
            prefix = 'ICA'

        fname_evo = op.join(sbj_path, '%sPSDTopo_%s%s' %
                            (prefix, cond, '-ave.fif'))
        print(fname_evo)
        mne.write_evokeds(fname_evo, psds_as_evo)

        fname_evo = op.join(sbj_path, '%sPSDTopoZ_%s%s' %
                            (prefix, cond, '-ave.fif'))
        print(fname_evo)
        mne.write_evokeds(fname_evo, psds_z_as_evo)

        fname_evo = op.join(sbj_path, '%sPSDHarm_%s%s' %
                            (prefix, cond, '-ave.fif'))
        print(fname_evo)
        mne.write_evokeds(fname_evo, psd_harm_as_evo)

        fname_evo = op.join(sbj_path, '%sPSDHarmBase_%s%s' %
                            (prefix, cond, '-ave.fif'))
        print(fname_evo)
        mne.write_evokeds(fname_evo, psd_harm_base_as_evo)

        fname_evo = op.join(sbj_path, '%sPSDSumTopoOdd_%s%s' %
                            (prefix, cond, '-ave.fif'))
        print(fname_evo)
        mne.write_evokeds(fname_evo, sum_odd_as_evo)

        fname_evo = op.join(sbj_path, '%sPSDSumTopoBase_%s%s' %
                            (prefix, cond, '-ave.fif'))
        print(fname_evo)
        mne.write_evokeds(fname_evo, sum_base_as_evo)

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
