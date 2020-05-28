#!/imaging/local/software/miniconda/envs/mne0.19/bin/python
"""
Compute PSDs in Source Space for average raw data for FPVS Frequency Sweep.

Average raw data from FPVS_get_sweeps.py.
Based on FOVS_PSD_sweep.py.
# Plot figures.
# Compute z-scores.
# Compute TFR if specified.
==========================================

OH, Possible better to add to existing function
"""

import sys

from os import path as op
import numpy as np

from copy import deepcopy

# import matplotlib
# matplotlib.use('Agg') #  for running graphics on cluster ### EDIT

from matplotlib import pyplot as plt

from importlib import reload

import mne
from mne.report import Report

import config_sweep as config
reload(config)

import FPVS_functions as Ff
reload(Ff)


print(mne.__version__)

# perform TFR of raw data or not
# do_tfr = config.do_tfr

# sub-directory for figures per subject
# separate for ICAed and non-ICAed data
if 'ica' in config.raw_ICA_suff:
    figs_dir = 'Figures_ICA'
else:
    figs_dir = 'Figures'

close_fig = 1  # close figures only if close_fig==1

# plt.ion() # interactive plotting

# for some plots of SNRs
unit_scalings = dict(eeg=1., mag=1., grad=1.)

print('Sunshine.')


def run_PSD_raw(sbj_id):
    """Compute spectra for one subject."""
    # initialise html report for one subject
    report = Report(subject=str(sbj_id), title='FPVS PSDs')

    subject = config.mri_subjects[sbj_id]

    # path to subject's data
    sbj_path = op.join(config.data_path, config.map_subjects[sbj_id][0])

    inv_fname = op.join(sbj_path, subject + '_EEGMEG-inv.fif')
    print('Reading EEG/MEG inverse operator: %s.' % inv_fname)

    invop = mne.minimum_norm.read_inverse_operator(inv_fname)

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
    sum_harms = {}
    for cond in conds:

        sum_harms[cond] = {}

    # Go through conditions and frequencies
    for cond in conds:  # conditions

        # create list of Evoked objects for all frequencies per condition
        psds_as_evo, psds_z_as_evo, sum_as_evo, psd_harm_as_evo,\
            psd_harm_base_as_evo = [], [], [], [], []

        if cond == 'face':  # hack, no frequency sweep for faces

            freqs = ['6.0']

        else:  # for all word condition, use all sweep frequencies

            freqs = freqs_all

        for freq in freqs:  # frequencies

            sum_harms[cond][freq] = []  # initialise for this base frequency

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

            psds

            # print('Computing psd_welch().')
            # psds, psd_freqs = mne.time_frequency.psd_welch(raw, fmin=0.,
            #                                                fmax=40.,
            #                                                n_fft=n_fft)

        #     print('Frequency resolution: %f.' % (psd_freqs[1] - psd_freqs[0]))

        #     info = raw.info

        #     # To plot PSDs like Evoked, pretend sample frequency
        #     info['sfreq'] = 1000. / (psd_freqs[1] - psd_freqs[0])

        #     # convert PSD to amplitudes (rather than power)
        #     psds = np.sqrt(psds)

        #     # "BASELINE-CORRECT" PSDs with neighbouring frequency bins
        #     print("\nCorrecting baseline.\n")
        #     if config.psd_base_bins > 0:

        #         psds = Ff.psd_correct_baseline(psds, config.psd_base_bins,
        #                                     config.psd_n_gap)

        #     print('Computing SNRs.')
        #     psds_z = Ff.psd_convert_to_snr(psds, config.psd_snr_bins,
        #                                 config.psd_n_gap)

        #     # Baseline-corrected PSDs as Evoked object
        #     as_evo = mne.EvokedArray(psds, info,
        #                              tmin=psd_freqs[0] / 1000.,
        #                              comment=('PSD ' + cond + ' ' + freq))
        #     psds_as_evo.append(as_evo)

        #     # z-scored PSDs as Evoked object
        #     as_evo = mne.EvokedArray(psds_z, info,
        #                              tmin=psd_freqs[0] / 1000.,
        #                              comment=('PSD_Z ' + cond + ' ' + freq))
        #     psds_z_as_evo.append(as_evo)

        #     # Plot PSD as spectrum plus topographies (plot_joint())
        #     print('Plotting PSDs.')

        #     # label for condition and base frequency
        #     label_str = '%s_%s' % (cond, ''.join(freq.split('.')))

        #     file_label = 'PSDtopo_%s' % label_str

        #     figs = Ff.plot_psd_as_evo(psds_as_evo[-1], sbj_path, file_label,
        #                            close_fig)

        #     chtypes = ['mag', 'grad']  # sequence of plots in plot_joint
        #     if len(figs) == 3:
        #         chtypes = chtypes + ['eeg']

        #     for [fig, chtype] in zip(figs, chtypes):

        #         sec_label = '%s_%s' % (label_str, chtype)

        #         report.add_figs_to_section(fig, sec_label, section=sec_label,
        #                                    scale=1)

        #     file_label = 'PSDtopoZ_%s' % label_str
        #     figs = Ff.plot_psd_as_evo(psds_z_as_evo[-1], sbj_path, file_label,
        #                            close_fig, scalings=unit_scalings)

        #     for [fig, chtype] in zip(figs, chtypes):

        #         sec_label = '%s_%s_Z' % (label_str, chtype)

        #         report.add_figs_to_section(fig, sec_label, section=sec_label,
        #                                    scale=1)

        #     # Compute the sum across harmonics of oddball frequency for this
        #     # condition and base frequency

        #     print('Combining harmonics.')

        #     if cond == 'face':

        #         # round to make sure combine_harmonics finds the right
        #         # frequencies
        #         oddfreq = round(config.fpvs_odd_freq['faces'], 2)

        #     else:

        #         oddfreq = round(config.fpvs_odd_freq['words'], 2)

        #     basefreq = float(freq)  # hack, float-to-string-to-float-again

        #     sum_harms[cond][freq] = Ff.combine_harmonics_topos(psds=psds_z,
        #                                                     freqs=psd_freqs,
        #                                                     basefreq=basefreq,
        #                                                     oddfreq=oddfreq,
        #                                                     n_harms=config.fpvs_n_harms,
        #                                                     method='sum')

        #     # needs another dimension for Evoked object
        #     to_evo = np.expand_dims(sum_harms[cond][freq], 1)

        #     # Combined harmonics as Evoked object
        #     as_evo = mne.EvokedArray(to_evo, info,
        #                              tmin=basefreq / 1000.,
        #                              comment=('Sum ' + cond + ' ' + freq))
        #     sum_as_evo.append(as_evo)

        #     # get PSDs around harmonics
        #     psd_harm = Ff.psds_across_harmonics(psds=psds_z, freqs=psd_freqs,
        #                                      basefreq=basefreq,
        #                                      oddfreq=oddfreq,
        #                                      n_harms=config.fpvs_n_harms,
        #                                      n_bins=config.psd_snr_bins,
        #                                      method='sum')

        #     # Sanity check - do it for base frequency
        #     # i.e. basefreq as oddfreq here, for all its harmonics
        #     psd_harm_base = Ff.psds_across_harmonics(psds=psds_z, freqs=psd_freqs,
        #                                           basefreq=999.,
        #                                           oddfreq=basefreq,
        #                                           n_harms=3,
        #                                           n_bins=config.psd_snr_bins,
        #                                           method='sum')

        #     tmin = -config.psd_snr_bins * 0.001  # include baseline
        #     info['sfreq'] = 1000.  # to display samples as time points

        #     as_evo = mne.EvokedArray(psd_harm, info, tmin=tmin,
        #                              comment=('PSD Harm ' + cond + ' ' + freq))
        #     psd_harm_as_evo.append(as_evo)

        #     as_evo = mne.EvokedArray(psd_harm_base, info, tmin=tmin,
        #                              comment=('PSD Harm ' + cond + ' ' + freq))
        #     psd_harm_base_as_evo.append(as_evo)

        #     # Plotting PSDs across harmonics
        #     fig = psd_harm_as_evo[-1].plot(spatial_colors=True,
        #                                    scalings=unit_scalings)

        #     # path to sub-directory for figures
        #     figs_path = op.join(sbj_path, figs_dir)

        #     fname_fig = op.join(figs_path, 'HarmPSD_%s_%s.pdf'
        #                         % (cond, freq))

        #     fig.savefig(fname_fig)

        #     sec_label = '%s_%s_HarmPSD' % (label_str, chtype)
        #     report.add_figs_to_section(fig, sec_label, section=sec_label,
        #                                scale=1)

        #     plt.close(fig)

        #     # Plotting PSDs across harmonics
        #     fig = psd_harm_base_as_evo[-1].plot(spatial_colors=True,
        #                                         scalings=unit_scalings)

        #     fname_fig = op.join(figs_path, 'HarmPSDbase_%s_%s.pdf'
        #                         % (cond, freq))

        #     fig.savefig(fname_fig)

        #     sec_label = '%s_%s_HarmPSDbase' % (label_str, chtype)
        #     report.add_figs_to_section(fig, sec_label, section=sec_label,
        #                                scale=1)

        #     plt.close(fig)

        #     times = [basefreq / 1000.]

        #     # Filename stem for figure; channel type to be added later
        #     fname_fig = op.join(figs_path, 'SumTopo_%s_%s' % (cond, freq))

        #     # For html report section label
        #     sec_label = '%s_%s_Sum' % (label_str, chtype)

        #     Ff.plot_evo_topomap(sum_as_evo[-1], times, chtypes, fname_fig, report,
        #                         label_str)

        # # Save HTML report
        # fname_report = op.join(figs_path, str(sbj_id) +
        #                        '_report.html')

        # report.save(fname_report, overwrite=True, open_browser=False)

        # # Save Evoked objects for later group stats:

        # print('Saving PSD results as evoked files:')

        # fname_evo = op.join(sbj_path, 'PSDtopo_%s%s' % (cond, '-ave.fif'))
        # print(fname_evo)
        # mne.write_evokeds(fname_evo, psds_as_evo)

        # fname_evo = op.join(sbj_path, 'PSDtopoZ_%s%s' % (cond, '-ave.fif'))
        # print(fname_evo)
        # mne.write_evokeds(fname_evo, psds_z_as_evo)

        # fname_evo = op.join(sbj_path, 'HarmPSD_%s%s' % (cond, '-ave.fif'))
        # print(fname_evo)
        # mne.write_evokeds(fname_evo, psd_harm_as_evo)

        # fname_evo = op.join(sbj_path, 'HarmPSDbase_%s%s' % (cond, '-ave.fif'))
        # print(fname_evo)
        # mne.write_evokeds(fname_evo, psd_harm_base_as_evo)

        # fname_evo = op.join(sbj_path, 'SumTopo_%s%s' % (cond, '-ave.fif'))
        # print(fname_evo)
        # mne.write_evokeds(fname_evo, sum_as_evo)

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
